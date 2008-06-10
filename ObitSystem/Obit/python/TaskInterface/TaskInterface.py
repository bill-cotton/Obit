# $Id: TaskInterface.py,v 1.1.1.1 2004/07/21 16:02:53 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2004
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
#
#                   Python tasking interface
#
#  This module contains classes useful for a tasking interface to python.
#  An TaskInterface object contains input parameters for a given program
#  (AIPS or other task) and starts the task passes this information to it.
#  Validity checking and saving/restoring objects is supported.
#     The parameters for a given task are defined in a Task Definition File
#  (TDF) which gives the order, names, types, ranges and dimensionalities.
#  A TDF is patterened after AIPS HELP files.
#     The Task Definition File can be derived from the AIPS Help file with the
#  addition of:
#   - A line starting with **TYPE** giving the type (AIPS or OBIT)
#   - A line before the beginning of each parameter definition of the form:
#   **PARAM** [type] [dim]
#       where [type] is float or str (string) and [dim] is the 
#       dimensionality as a blank separated list of integers, e.g.
#       **PARAM** str 12 5       (5 strings of 12 characters)
#   HINT: No matter what POPS thinks, all strings are multiples of 4 characters
#   For non AIPS usage dbl (double), int (integer=long), boo (boolean)
#   are defined.
#   AIPS Parameter names are all upper case and all numerics are floats
#   (i.e. require ".0" on interger values)
#
#    Currently AIPS and Obit (sorta) tasks are supported.
#    All tasks are run synchronously on the local host.
#    A RunTime error is thrown if an error is encountered.
#
#    The typical usage sequence is:
#     1) Create object using TaskInterface.TaskInterface(name, file)
#        where name is the task name and file the name of the TDF
#        The parameters will be initialized with the default values
#        example: >>> xxx=TaskInterface.TaskInterface("NOBAT","NOBAT.TDF")
#     2) Set parameter values using setValue
#        example: >>> xxx.setValue("DETIME",0.0)
#     3) If running AIPS set user ID and POPS number
#        example: >>> xxx.Auserid=100; xxx.ANPops=1
#     4) Inputs can be reviewed:
#        >>> xxx.Input()   ( also .Help() and .Explain())
#        or checked (type, dimension, range):
#        >>> xxx.Check()
#     5) The task can be executed:
#        >>> xxx.Go()
#     6) Values can be retrieved:
#        >>> detime = xxx.getValue("DETIME")
#     7) Current state of the object can be saved and restored using
#        Put() and Get().  The default value is <task_name>.inp
#        (done automatically by Go()) or in a named file.
#     8) Parameter values can be reset to the default using 

# Notes on python AIPS/Obit tasking framework
# 
# All:
#   - task Put/Get automatically (where?)
# 
# for Obit:
#   - Proper bin directory
# 
# for AIPS:
#    - set $DA??, $LOAD to run correctly
#    - better managment for TDF, StartAIPSTasks files

import os, pickle, string, struct

class TaskInterfaceParam:
    """
    Task parameter and interface object

    This class contains the information for a single task parameter
    which may be a single, or array of values all of the same type.
    This is description, the actual data is kept in a separate dictionary
    using self.name as the key.
    """
    def __init__(self, name) :
        """Create new TaskInterfaceParam"""
        self.name    = name   # Parameter name
        self.type    = None   # type string, "str","float" (AIPS)
                              # also "dbl", "bool", "int" (OBIT)
        self.dim     = [1]    # array dimensionality (or length of string)
        self.range   = [-1.0e20, 1.0e20] # range of allowed values for numeric
        self.default = 0.0    # Scalar or list of default values
        self.help    = None   # Short description string
        self.AIPScode=' '     # AIPS parameter code:
        # ' '   GO
        # '*'   GO  TELL
        # '?'       TELL
        # '&'   GO  TELL  GET
        # '%'       TELL  GET
        # '$'   GO        GET
        # '@'             GET
        
        
    # end init
    
    def setDefault(self, value):
        """ Set values to default 
        Inputs:
        value = data value dictionary
        """
        ################################################################
        value[self.name] = []
        for x in self.default:
            value[self.name].append(x);
        # end default
        
    def check(self, value):
        """ Check is data right type and in range
        Inputs:
        value = data value dictionary
        returns True if OK, else throws Runtime exception
        """
        ################################################################
        # Check dimensionality
        total = 1
        start = 0
        if (self.type==str):  # Strings are different
            start = 1
        for x in self.dim[start:]:
            total = total * x
        dlen = 1
        if value[self.name].__class__==list:
            dlen = len(value[self.name])
        if total!=dlen:
            print "Dimension",dlen,"total actual=",total
            msg = "Dimensionality wrong for parameter "+self.name
            raise RuntimeError,msg
        
        # Check type and value
        numeric = ((self.type==float) | (self.type=='dbl')) | (self.type==int)
        myvalue = value[self.name]    # extract value array from dictionary
        if myvalue.__class__==list:
            for x in myvalue:
                if (x.__class__!=self.type):
                    print "actual type",x.__class__,"should be",self.type
                    msg = "Value type wrong for parameter "+self.name
                    raise RuntimeError,msg
                if numeric:
                    if (x<self.range[0]) | (x>self.range[1]):
                        print x,"Not in range",self.range
                        msg = "Value out of range for parameter "+self.name
                        raise RuntimeError,msg
        else:  # Check scalar
            if (myvalue.__class__!=self.type):
                print "actual type",myvalue.__class__,"should be",self.type
                msg = "Value type wrong for parameter "+self.name
                raise RuntimeError,msg
            if numeric:
                if (myvalue<self.range[0]) | (myvalue>self.range[1]):
                    print myvalue,"Not in range",self.range
                    msg = "Value out of range for parameter "+self.name
                    raise RuntimeError,msg
        return True  # OK
        # End check

    def setAIPS(self, fd, value):
        """
        Write parameter value to AIPS TD file opened as fd
        Inputs:
        fd    = open TD file positioned for write
        value = data value dictionary
        """
        ################################################################
        # Is this a "GO" parameter
        code = self.AIPScode
        if (code==' ') | (code=='*') | (code=='&') | (code=='$'):
            myvalue = value[self.name]    # extract value array from dictionary
            # list and scalar separate
            if myvalue.__class__==list: # List
                # AIPS only supports float and strings
                if self.type==float:     # Float
                    for x in myvalue:
                        fd.write(struct.pack('f',x))
                elif self.type==str:     # String
                    for x in myvalue:
                        # Blank pad to proper length
                        xx = string.ljust(string.strip(x),self.dim[0])  
                        fmt = str(self.dim[0])+'s'
                        fd.write(struct.pack(fmt,xx))
            else:  # Scalar
                if self.type==float:     # Float
                    fd.write(struct.pack('f',myvalue))
                elif self.type==str:     # String
                    # Blank pad to proper length
                    xx = string.ljust(string.strip(myvalue),self.dim[0])  
                    fmt = str(self.dim[0])+'s'
                    fd.write(struct.pack(fmt,xx))
        # end setAIPS
        

    def getAIPS(self, fd, value):
        """
        Read parameter value from AIPS TD file opened as fd
        Inputs:
        fd    = open TD file positioned for read
        value = data value dictionary
        """
        ################################################################
        # Is this a "GET" parameter
        code = self.AIPScode
        if (code=='&') | (code=='%') | (code=='$') | (code=='@'):
            myvalue = value[self.name]    # extract value array from dictionary
            # list and scalar separate
            if myvalue.__class__==list: # List
                out = []
                # AIPS only supports float and strings
                if self.type==float:     # Float
                    for x in myvalue:
                        y = struct.unpack('f',fd.read(4))
                        out.append(y[0])
                elif self.type==str:     # String
                    for x in myvalue:
                        y = fd.read(self.dim[0])
                        out.append(y[0])
            else:  # Scalar
                if self.type==float:     # Float
                    out = struct.unpack('f',fd.read(4))[0]
                elif self.type==str:     # String
                    out = string.strip(fd.read(self.dim[0]))[0]
            value[self.name] = out       # Save value from file
        # end getAIPS
        

    def setOBIT(self, fd, value):
        """
        Write parameter value to OBIT Input file opened as fd
        Inputs:
        fd    = open input file positioned for write
        value = data value dictionary
        """
        ################################################################
         # Is this a "GO" parameter
        code = self.AIPScode
        if (code==' ') | (code=='*') | (code=='&') | (code=='$'):
            myvalue = value[self.name]    # extract value array from dictionary
            data = string.strip(str(myvalue),"[]")   # remove braces
            data = string.replace(data,",", " ")     # remove commas
            
            if self.type==float:     # Float
                fd.write("$Key = "+self.name+" Flt  ("+str(self.dim[0])+")\n ")
                print "$Key = "+self.name+" Flt  ("+str(self.dim[0])+") "
            elif self.type==str:     # String
                fd.write("$Key  "+self.name+" Str ("+str(self.dim[0])+") \n")
                print "$Key  "+self.name+" Str ("+str(self.dim[0])+") "
            elif self.type==bool:   # Boolean
                fd.write("$Key  "+self.name+" Boo ("+str(self.dim[0])+") \n")
                print "$Key  "+self.name+" Boo ("+str(self.dim[0])+") "
            elif self.type==int:     # Integer
                fd.write("$Key Int "+self.name+" ("+str(self.dim[0])+") \n ")
                print "$Key Int "+self.name+" ("+str(self.dim[0])+") "
                fd.write(data+" \n")
                print data+" "
        # end setOBIT 
        

class TaskInterface:
    """
    Task inputs and execution from Python class
    
    This class manipulates inputs objects and online documentation.
    Objects can be instantiated from an AIPSlike input description file
    and stored to or from a file.  Tasks can also be started and
    run synchronously.
    Values are saved and retrieved using setValue and getValue.

    Constructor:
    TaskInterface.TaskInterface(name, file)

    Inputs:
       name = name of task
       file = if given and non None, the name of the TDF file
              if not specified a skeleton object is returned.
    """
    def __init__(self, name, file=None) :
        """
        Create new TaskInterface
        For AIPS Usage set members:
            Auserid = AIPS user ID
            ANpops  = Pops number (1 good default)
         
        Inputs:
            name = name (task) of the object
            file = file name of the task definition file
                   if None do not read file
        """
        self.name = name       # Object name
        self.sync = False      # Run task synchronously?
        self.host = 'local'    # Host on which to run, "local"=local
        self.file = file       # Name of file from which defined
        self.type = 'OBIT'     # Task type to start
        self.help = None       # Short task description string
        self.param= []         # List of parameter definition objects
        self.value= {}         # Dictionary of actual parameter values
        self.Auserid=1         # AIPS user ID
        self.ANpops=1          # AIPS Pops number

        # Parse file
        if file==None:
            return
        fd = open(file)
        line = fd.readline();
        while line!="":
            if line[0:3]==";! ":   # Get description
                self.help = string.strip(line[2:])
            elif line[0:8]=='**TYPE**':     # Line with task type
                self.type = string.strip(line[8:])
            elif line[0:9]=='**PARAM**':    # beginning of parameter definition
                # Get type and dimension
                parts = string.split(line)
                type  = parts[1]
                if type=="float":
                    type = float
                elif type=="str":
                    type = str
                elif type=="int":
                    type = int
                elif type=="boo":
                    type = bool
                # Dimensionality
                dim = []
                total = 1
                for x in parts[2:]:
                    total *= int(x)
                    dim.append(int(x));
                # Want number of strings, not no. characters
                if type==str:
                    total /= dim[0];
                # Read first line of parameter definition
                line  = fd.readline();
                # If no range given all are blank
                name  = string.rstrip(line[0:8])
                if line[10:22]!="            ":
                    r1 = float(line[10:22])
                else:
                    r1 = -1.0e30
                if line[22:34]!="            ":
                    r2 = float(line[22:34])
                else:
                    r2 = 1.0e30
                help  = string.rstrip(line[35:])
                acode = line[9:10]
                p = TaskInterfaceParam(name)
                p.type    = type
                p.dim     = dim
                p.AIPScode= acode
                p.range   = [r1,r2]
                p.help    = help
                # Default by type
                p.default = []
                if type==float:
                    for i in range(total):
                        p.default.append(0.0)
                if type==str:
                    for i in range(total):
                        p.default.append("")
                if type==int:
                     for i in range(total):
                        p.default.append(0)
                if type==bool:
                     for i in range(total):
                        p.default.append(False)
                if type=="dbl":
                    for i in range(total):
                        p.default.append(0.0)
                self.param.append(p)   # Add to list
                # Set to default
                self.value[p.name] = p.default
            line = fd.readline();      # next line
        fd.close()
        # end init
        
    def Go(self):
        """
        Execute the task
        Currently all run synchronously
        A RuntimeError is thrown if the inputs are invalid or the task fails.
        The contents are checked and saved before executing the task
        Obit: copy parameters into temporary input (ObitParser) file
              and execute program
        AIPS: write task info in TD file and execute script StartAIPSTask
              passing it the name of a symbolic link
              Needs $LOAD, $DA?? defined
        """
        ################################################################
        # Validity check
        self.Check()
        # Save parameters
        self.Put()
        # By task type
        if self.type=="OBIT":
            pid = os.getpid()
            tmpInput = "/tmp/"+self.name+"Input."+str(pid)
            fd = open(tmpInput,"w")
            # Loop over parameters
            for x in self.param:
                x.setOBIT(fd, self.value)
            fd.close()
            # run it
            os.spawnv(os.P_WAIT,self.name,[""," input = "+tmpInput])
            # Note: using os.P_NOWAIT will cause an asynchronous call
            # returning the PID of the process, then wait for the process
            # with os.waitpid(pid,0)
            #pid=os.spawnv(os.P_NOWAIT,self.name,[""," input = "+tmpInput])
            #os.waitpid(pid,0)
            #print "PID",pid
            os.unlink(tmpInput)                      # remove input file
            #
        elif self.type=="AIPS":
            # Note, this will really only work on the local host
            da00 = os.getenv("DA00")
            AIPSTD = da00+"/TDD000004;"
            npops = self.ANpops         # Pops number
            if (npops<1) | (npops>15):
                raise RuntimeError,"Illegal POPS number "+str(npops)
            lenCB = 256 * 4             # Length of control section
            lenTD = 1024 * 4            # Length of task data section
            fd = open(AIPSTD,"r+b")
            # Write control record
            fd.seek((npops-1)*20,0)
            tname = string.ljust(string.strip(self.name),8)
            fd.write(struct.pack('8s',tname))   # task name
            fd.write(struct.pack('l',-999))     # initial return code
            fd.write(struct.pack('2l',0,0))
            #
            # Task data area
            fd.seek(lenCB+(npops-1)*lenTD,0)
            # Task independent data
            userNo = self.Auserid # User number
            TVnumber = 0          # No AIPS TV allowed
            TKnumber = 0          # No AIPS graphics device
            MsgKill  = 0          # All AIPS messages
            IsBatch  = 32000      # Run in batch mode
            DbgAIPS  = 0          # No debugger
            DoWait   = 1.0        # Run in "wait" mode
            version  = "TST:"     # AIPS version
            # values to file
            fd.write(struct.pack('l',userNo))
            fd.write(struct.pack('l',TVnumber))
            fd.write(struct.pack('l',TKnumber))
            fd.write(struct.pack('l',MsgKill))
            fd.write(struct.pack('l',IsBatch))
            fd.write(struct.pack('l',DbgAIPS))
            fd.write(struct.pack('2l',0,0))    # reserved
            fd.write(struct.pack('f',DoWait))
            fd.write(struct.pack('4s',version))
            # Loop over parameters
            for x in self.param:
                x.setAIPS(fd, self.value)
            fd.close()
            #
            # Now start up task, link with name*pops number in /tmp
            # and run from there using script StartAIPSTask
            # Where is the executable really?
            loadDir = os.getenv("LOAD")
            taskExec = loadDir+"/"+string.upper(self.name)+".EXE"
            # Symbolic link to be executed so task knows its POPS no.
            pid = os.getpid()
            taskLink = "/tmp/"+string.upper(self.name)+str(npops)+"."+str(pid)
            os.symlink(taskExec, taskLink)           # Set execution link
            # Let 'er rip
            os.spawnv(os.P_WAIT,"StartAIPSTask",["",taskLink])
            os.unlink(taskLink)                      # remove link executed
            
            # Get results from TD file
            fd = open(AIPSTD,"rb")
            # Read return code
            fd.seek((npops-1)*20+8,0)
            controlRec = fd.read(4)
            retCode = struct.unpack('l',controlRec)
            print "debug return code",retCode[0]
            # If it's not 0 it failed
            if retCode[0]!=0:
                raise RuntimeError,"Task "+self.name+" returns code "+str(retCode[0])

            # Loop over parameters
            fd.seek(lenCB+(npops-1)*lenTD+40,0)
            for x in self.param:
                x.getAIPS(fd, self.value)
            fd.close()
        # end Go

    def Wait(self):
        """
        Suspend until associated task terminates
        Not yet implemented
        """
        ################################################################
        # end Wait

    def Input(self):
        """
        Print the contents of an input object
        """
        ################################################################
        print 'Input Object',self.name
        print self.help
        print "From file",self.file,", type ",self.type
        for x in self.param:
            print x.name,x.type,x.dim,x.help,"\n     ",self.value[x.name]
        # end Input

    def Help(self):
        """
        Print the contents of the help section of associated task
        definition file.
        """
        ################################################################
        # Parse file
        fd = open(self.file)
        line = fd.readline();
        # Skip inputs section
        while line!="":
            if line[0:5]=='-----':
                break
            line = fd.readline();
        #
        # Print help section
        line = fd.readline();
        while line!="":
            if line[0:5]=='-----':
                break
            if line[0:1]!=';':
                print string.rstrip(line)
            line = fd.readline();
        fd.close()
        # end Help

    def Explain(self):
        """ Print the contents of the explain section of associated task
        definition file.
        """
        ################################################################
         # Parse file
        fd = open(self.file)
        line = fd.readline();
        # Skip inputs section
        while line!="":
            if line[0:5]=='-----':
                break
            line = fd.readline();
        #
        # Print help section
        line = fd.readline();
        while line!="":
            if line[0:1]!=';':
                print string.rstrip(line)
            line = fd.readline();
        fd.close()
       # end Explain

    def setDefault(self):
        """ Set values to default """
        for x in self.param:
        ################################################################
            x.setDefault(self.value)
        # end setDefault
        
    def setValue(self, name, value):
        """
        Set parameter value, consistency with definition is checked
        Throws Runtime exception on error
        
        Inputs:
        name  = parameter name
        value = scalar or list of values 
        """
        ################################################################
        # Get type of value
        type = value.__class__
        length = 1                 # Number of elements in value
        if type==list:
            type = value[0].__class__
            length = len(value)
        found = False    # Haven't found name yet
        # Look up name in parameter list
        for x in self.param:
            if x.name==name:   # Is this it?
                found = True
                # Check that input scalar or list as needed
                if value.__class__!=x.type:   # Wrong list/scalar
                    print "Want value type",x.type,"given",type
                    msg = name+" not scalar or list as needed"
                    raise RuntimeError,msg
                # Check value type
                if type!=x.type:   # Bad data type
                    print "Correct type",x.type,"given",type
                    msg = "Wrong data type for "+name
                    raise RuntimeError,msg
                # Check Dimension
                total = 1
                start = 0
                if (x.type==str):  # Strings are different
                    start = 1;
                for y in x.dim[start:]:
                    total = total * y
                if total!=length:    # Bad dimensionality
                    print "total dimensionality",total,"given",length
                    msg = "Wrong dimensionality for "+name
                    raise RuntimeError,msg
                # OK if it gets here - set value
                self.value[name] = value
                # Check validity (range)
                x.check(self.value)
        # Was it found?
        if found==False:
            msg = name+" not a parameter for "+self.name
            raise RuntimeError,msg
        # end setValue
    
    def getValue(self, name):
        """
        Returns parameter value
        Throws Runtime exception on error (not found)
        Inputs:
        name  = parameter name
        Returns a scalar or list of values
        """
        ################################################################
        # Look up name in parameter list
        for x in self.param:
            if x.name==name:   # Is this it?
                return self.value[name]
        # Oops not in list
        msg = name+" not a parameter for "+self.name
        raise RuntimeError,msg
        # end getValue
        
    def Check(self):
        """ Check all parameters"""
        ################################################################
        for x in self.param:
            x.check(self.value)
        # End check

    def Get (self, file=None):
        """
        Read an input from persistent (disk) form
        Inputs:
            file = disk file name to read
                   default is self.name.inp
        """
        ################################################################
        tfile = file
        if file==None:
            tfile = self.name+".inp"
        fd = open(tfile, "r")
        self = pickle.load(fd)
        fd.close()
        # end Get

    def Put (self, file=None):
        """
        Write an input to persistent (disk) form
        Inputs:
            file = disk file name to write
                   default is self.name.inp
        """
        ################################################################
        tfile = file
        if file==None:
            tfile = self.name+".inp"
        fd = open(tfile, "w")
        pickle.dump(self,fd)
        fd.close()
        # end Put
        
    def PIsA (inTaskInterface):
        """ Tells if input really a Python Obit TaskInterface
        
        return true, false (1,0)
        inTaskInterface   = Python TaskInterface object
        """
    ################################################################
        # Checks
        if inTaskInterface.__class__ != TaskInterface:
            print "Class actually is",inTaskInterface.__class__
            return 0
        return 1
    # end PIsA

