"""
                   Obit tasking interface

  This module contains classes useful for an Obit tasking interface to python.
  An ObitTask object contains input parameters for a given Obit program.
     The parameters for a given task are defined in a Task Definition File
  (TDF) which gives the order, names, types, ranges and dimensionalities.
  A TDF is patterened after AIPS HELP files.
     The Task Definition File can be derived from the AIPS Help file with the
  addition of:
   - A line before the beginning of each parameter definition of the form:
   **PARAM** [type] [dim] **DEF** [default]
       where [type] is float or str (string) and [dim] is the 
       dimensionality as a blank separated list of integers, e.g.
       **PARAM** str 12 5       (5 strings of 12 characters)
       default (optional) is the default value
   HINT: No matter what POPS thinks, all strings are multiples of 4 characters
   For non AIPS usage dbl (double), int (integer=long), boo (boolean
   "T" or "F") are defined.

"""
#-----------------------------------------------------------------------
# Copyright (C) 2005,2007,2010 Associated Universities, Inc. Washington DC, USA.
# Copyright (C) 2005 Joint Institute for VLBI in Europe
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

# Bits from AIPS.
from Proxy.AIPS import ehex

# Bits from the generic Task implementation.
from Proxy.Task import Task

# Generic Python stuff.
import fcntl, glob, os, pickle, select, struct, string, pty


class _ObitTaskParams:
    """ Obit Task parameter class """
    def __parse(self, name):
        """Determine the proper attributes for the Obit task NAME by
        parsing its TDF file."""

        task = None
        strlen = None
        deff = [0]
        gotDesc = False
        min = None
        max = None
        adverb = None
        hlp = None
        short_help = None

        path = name + '.TDF'
        if (not os.access(path,os.F_OK)) and  (os.getenv("OBITEXEC") != None):
            path = os.environ['OBITEXEC'] + '/TDF/' + name + '.TDF'
        if (not os.access(path,os.F_OK)) and  (os.getenv("OBITEXEC") != None):
            path = os.environ['OBITEXEC'] + '/tdf/' + name + '.TDF'
        if (not os.access(path,os.F_OK)) and  (os.getenv("OBIT") != None):
            path = os.environ['OBIT'] + '/TDF/' + name + '.TDF'
        # Check OBITSD if needbe
        if  (not os.access(path,os.F_OK)) and (os.getenv('OBITSD') != None):
            path = os.getenv('OBITSD') +'/TDF/' + '/' + name+ '.TDF'
        # Check standard linux directory if needbe
        if  (not os.access(path,os.F_OK)):
            path = '/usr/lib/obit/tdf/' + '/' + name+ '.TDF'
        
        # Better have found it by here
        if  (not os.access(path,os.F_OK)):
            # Oh bugger
            msg = "Task '%s' task definition file not found" % (name)
            print msg
            raise RuntimeError, msg

        input = open(path)
        line = " "
        while (line):
            #for line in input:
            line = input.readline() # read next line
            if not line:  # EOF?
                break
            # A line of dashes terminates the parameter definitions.
            if line.startswith('--------'):
                break;

            # Comment lines start with ';'.
            if line.startswith(';'):
                continue

            if not task:
                task = line.split()[0]
                min_start = line.find('LLLLLLLLLLLL')
                min_end = line.rfind('L')
                max_start = line.find('UUUUUUUUUUUU')
                max_end = line.rfind('U')
                dir_start = min_start - 2
                dir_end = min_start - 1
                continue

            if line.startswith(task):
                # get short help string
                if  not short_help:
                    short_help = line
                continue

            if line.startswith(' ') or line.startswith('\n'):
                # Continuation of parameter description?
                if adverb and len(line)>35:
                    self.hlp_dict[adverb].append(string.rstrip(line[35:]))
                continue

            # Description of parameter?
            if line.startswith('**PARAM**'):
                gotDesc = True
                # Get type and dimension.
                parts = line.split()
                 # Dimensionality.
                dim = []
                total = 1
                count = 2
                gotDef = False
                dtype = None
                #print "DEBUG parts",parts
                for x in parts[2:]:
                    if x == "**DEF**":
                        gotDef = True
                        break
                    count = count + 1
                    total *= int(x)
                    dim.append(int(x));
                # Default value?
                if gotDef:
                    ddeff = parts[count+1]
                # Want number of strings, not number of characters.
                if parts[1] == 'str':
                    total = total / dim[0]
                # Type.
                type = parts[1]
                if type == 'float':
                    dtype = 'Flt'
                    type = float
                    strlen = None
                    if gotDef:
                        ddeff = float(ddeff)
                        deff = total * [ddeff]
                    else:
                        deff = total * [0.0]
                elif type == 'str':
                    dtype = 'Str'
                    type = str
                    strlen = dim[0]
                    if gotDef:
                        deff = total * [ddeff]
                    else:
                        deff = total * [strlen * ' ']
                elif type == 'int':
                    dtype = 'Int'
                    type = int
                    strlen = None
                    if gotDef:
                        ddeff = int(ddeff)
                        deff = total * [ddeff]
                    else:
                        deff = total * [0]
                elif type == 'boo':
                    type = bool
                    dtype = 'Boo'
                    strlen = None
                    if gotDef:
                        ddeff = ddeff=="T"
                        deff = total * [ddeff]
                    else:
                        deff = total * [False]
                elif type == 'dbl':
                    dtype = 'Dbl'
                    strlen = None
                    if gotDef:
                        ddeff = float(ddeff)
                        deff = total * [ddeff]
                    else:
                        deff = total * [0.0]
                        #print "DEBUG line",line,type,dim,deff
                else:
                    msg = "UNKNOWN TYPE: %s " % (type)
                    raise RuntimeError, msg
                continue

            # If just parsed PARAM line get parameter.
            if gotDesc:
                gotDesc = False
                adverb = line.split()[0]
                code = line[min_start - 1:min_start]
                hlp = [string.rstrip(line[35:])]
                if not code:
                    code = ' '
                try:
                    min = float(line[min_start:min_end].strip())
                    max = float(line[max_start:max_end].strip())
                except:
                    min = None
                    max = None

                # Assume type/dimension is one just read.
                # If only one entry, convert deff to scalar.
                #print "DEBUG adverb", adverb, deff, dim, code
                if len(deff) == 1:
                    deff = deff[0]
                self.default_dict[adverb] = deff # default
                self.dim_dict[adverb] = dim # dimensionality
                if code in ' *&$' or len(adverb) > 9:
                    self.input_list.append(adverb)
                if code in '&%$@':
                    self.output_list.append(adverb)
                if strlen:
                    self.strlen_dict[adverb] = strlen
                if min != None:
                    self.min_dict[adverb] = min
                if max != None:
                    self.max_dict[adverb] = max
                if hlp != None:
                    self.hlp_dict[adverb] = hlp
                if dtype != None:
                    self.type_dict[adverb] = dtype


        # Save short description
        if short_help != None:
            self.short_help = short_help

        # Parse HELP section.
        while(line):
            line = input.readline() # read next line
            if not line:  # EOF?
                break
            # A line of dashes terminates the help message.
            if line.startswith('--------'):
                break;
            self.help_string = self.help_string + line

        # Parse EXPLAIN section. - rest of file
        while(line):
            line = input.readline() # read next line
            if not line:  # EOF?
                break
            #for line in input:
            self.explain_string = self.explain_string + line

        # Add "retCode" to output parmeters, -1 => task not finished
        adverb = "retCode"
        self.default_dict[adverb] = -1 # default
        self.dim_dict[adverb] = [1]
        self.hlp_dict[adverb] = ["Task return code, 0=OK, ",\
                                 "-1=> not finished or failed"]
        self.output_list.append(adverb)

    def __init__(self, name, version):
        """ Create Server Obit task structure
        
        default_dict   = Dictionary with default values of parameters
        input_list     = List of input parameters in order
        output_list    = List of output parameters in order
        min_dict       = Parameter minimum values as a List
        max_dict       = Parameter maximum values as a List
        hlp_dict       = Parameter descriptions (list of strings)
                          as a dictionary
        strlen_dict    = String parameter lengths as dictionary
        help_string    = Task Help documentation as list of strings
        explain_string = Task Explain documentation as list of strings
        short_help     = One line description of task

        """
        self.default_dict = {}
        self.input_list = []
        self.output_list = []
        self.min_dict = {}
        self.max_dict = {}
        self.strlen_dict = {}
        self.hlp_dict = {}
        self.dim_dict = {}
        self.type_dict = {}
        self.short_help = ''
        self.help_string = ''
        self.explain_string = ''
  
        self.name = name
        if version in ['OLD', 'NEW', 'TST']:
            self.version = os.path.basename(os.environ[version])
        else:
            self.version = version

        # The following causes a huge amount of trouble and little gain
        #path = os.environ['HOME'] + '/.ObitTalk/' \
        #       + self.version + '/' + name + '.pickle'
        #try:
        #    unpickler = pickle.Unpickler(open(path))
        #    self.default_dict = unpickler.load()
        #    self.input_list = unpickler.load()
        #    self.output_list = unpickler.load()
        #    self.min_dict = unpickler.load()
        #    self.max_dict = unpickler.load()
        #    self.strlen_dict = unpickler.load()
        #    self.dim_dict = unpickler.load()
        #    self.hlp_dict = unpickler.load()
        #    self.help_string = unpickler.load()
        #    self.explain_string = unpickler.load()
        #    self.short_help = unpickler.load()
        #except (IOError, EOFError):
        #    self.__parse(name)
        #
        #    # Make sure the directory exists.
        #    if not os.path.exists(os.path.dirname(path)):
        #        os.makedirs(os.path.dirname(path))
        #
        #    pickler = pickle.Pickler(open(path, mode='w'))
        #    pickler.dump(self.default_dict)
        #    pickler.dump(self.input_list)
        #    pickler.dump(self.output_list)
        #    pickler.dump(self.min_dict)
        #    pickler.dump(self.max_dict)
        #    pickler.dump(self.strlen_dict)
        #    pickler.dump(self.dim_dict)
        #    pickler.dump(self.hlp_dict)
        #    pickler.dump(self.help_string)
        #    pickler.dump(self.explain_string)
        #    pickler.dump(self.short_help)

        # instead - parse TDF file each time
        self.__parse(name)

    # Provide a dictionary-like interface to deal with the
    # idiosyncrasies of XML-RPC.
    def __getitem__(self, key):
        return self.__dict__[key]


class ObitTask(Task):
    """ Server-side Obit task interface """
    def __init__(self):
        Task.__init__(self)
        self._params = {}
        self._popsno = {}

    def params(self, name, version):
        """Return parameter set for version VERSION of task NAME."""
        return _ObitTaskParams(name, version)

    def __write_adverb(self, params, file, adverb, value, dtype):
        """Write Obit input text file."""

        assert(adverb in params.input_list)

        # Get type, may be scalar or list
        #dtype = type(value)
        #if dtype == list:
        #    dtype = type(value[0])

        # Convert to string for numeric types
        if type(value) == list:
            data = string.join(map(str, value))
        else:
            data = str(value)

        dim = params.dim_dict[adverb]   # Dimensionality array
        dimStr = "(" + str(dim[0]) + ")"
        if (len(dim) > 1):
            if (dim[1] > 1):
                dimStr = "(" + str(dim[0]) + "," + str(dim[1]) + ")"

        if dtype == "Flt":
            file.write("$Key = " + adverb + " Flt " + dimStr + "\n")
            file.write(data + "\n")     # Write data to file
        elif dtype == "Dbl":
            file.write("$Key = " + adverb + " Dbl " + dimStr + "\n")
            file.write(data + "\n")     # Write data to file
        elif dtype == "Str":
            file.write("$Key = " + adverb + " Str " + dimStr + "\n")
            if type(value) == list:
                for x in value:
                    file.write(x + "\n") # Write data to file
            else:
                #print "DEBUG write_adverb", adverb, dtype, dim, value
                file.write(value + "\n") # Write data to file
        elif dtype == "Boo":
            file.write("$Key = " + adverb + " Boo " + dimStr + "\n")
            if type(value) == list:
                #print "DEBUG value", adverb, value
                for x in value:
                    if x:
                        file.write(" T") # Write data to file.
                    else:
                        file.write(" F")
            else:
                if value:
                    file.write(" T")    # Write data to file.
                else:
                    file.write(" F")
            file.write("\n")
        elif dtype == "Int":
            file.write("$Key = " + adverb + " Int " + dimStr + "\n")
            file.write(data + "\n" )    # Write data to file.
        else:
            print "DEBUG ObitTask adverb", adverb, dim, dtype
            raise AssertionError, type(value)

    def __read_adverb(self, params, file, adverb):
        """read value from task output file."""

        assert(adverb in params.output_list)

        gotIt = False    # Not yet found entry
        count = 0        # no values parset yet
        total = 1
        value = []       # to accept
        for line in file:
            # DEBUG
            #print line
            # Look for header for parameter
            if line.startswith("$Key = " + adverb):
                # remove unwanted symbols
                line = line.replace('=',' ')
                line = line.replace('(',' ')
                line = line.replace(',',' ')
                line = line.replace(')',' ')
                gotIt = True
                parts = string.split(line)
                # How many values
                total = 1
                # DEBUG print parts  
                for x in parts[3:]:
                    total *= int(x)
                dtype  = parts[2]
                if type=="str":
                    total = total / parts[3] # Number of strings.

            # Read data one value per line after 'gotIt'.
            elif gotIt:
                # DEBUG print "gotIt", type, line
                if dtype == 'Flt':
                    value.append(float(line))
                elif dtype == 'Dbl':
                    value.append(float(line))
                elif dtype == 'Str':
                    value.append(line)
                elif dtype == 'Int':
                    value.append(int(line))
                elif dtype == 'Boo':
                    if line.startswith('T'):
                        value.append(True)
                    else:
                        value.append(False)
                count = count + 1

            if gotIt and count >= total:  # Finished?
                    break

        # Convert to scalar if only one.
        if len(value)==1:
            value = value[0]
        # Done
        # DEBUG print "fetch adverb", adverb, value
        return value

    def spawn(self, name, version, userno, msgkill, isbatch, input_dict):
        """Start the task.
        
        Writes task input parameters, data directories and other
        information in the task parameter file and starts the task
        asynchronously returning immediately.  Messages must be
        retrieved calling messages.
        If the debug parameter is True then a "Dbg" copy of the task
        input file is created which will not br automatically destroyed.
        name        task name
        version     version of task
        userno      AIPS user number
        msgkill     AIPS msgkill level, not used in Obit tasks
        isbatch     True if this is a batch process , not used in Obit tasks
        input_dict  Input parameters as dictionary
        Returns task id
        """

        params = _ObitTaskParams(name, version)
        popsno = _allocate_popsno()
        index = popsno - 1

        # Set input and output text parameter/result files
        tmpInput = "/tmp/" + params.name + "Input." + str(popsno)
        tmpOutput = "/tmp/" + params.name + "Output." + str(popsno)

        in_file = open(tmpInput, mode="w")

        # Add user id input file
        in_file.write("$Key = AIPSuser Int (1)\n")
        in_file.write(str(userno)+"\n")

        # Add directories to input file - both types for 'AIPS'
        inType = input_dict['DataType']
        outType = input_dict['DataType']
        if "outDType" in input_dict:
            outType = input_dict["outDType"]
        if  (inType=='AIPS') or (outType=='AIPS'):
            AIPSdirs = input_dict["AIPSdirs"]
            nAIPS = len(AIPSdirs)
            in_file.write("$Key = nAIPS Int (1)\n")
            in_file.write(str(nAIPS)+"\n")
            in_file.write("$Key = AIPSdirs Str (128,"+str(nAIPS)+")\n")
            for x in AIPSdirs:
                in_file.write(x+"\n")
            FITSdirs = input_dict["FITSdirs"]
            nFITS = len(FITSdirs)
            in_file.write("$Key = nFITS Int (1)\n")
            in_file.write(str(nFITS)+"\n")
            in_file.write("$Key = FITSdirs Str (128,"+str(nFITS)+")\n")
            for x in FITSdirs:
                in_file.write(x+"\n")

        elif (input_dict['DataType']=='FITS') or (input_dict['outDType']=='FITS'):
            FITSdirs = input_dict["FITSdirs"]
            nFITS = len(FITSdirs)
            in_file.write("$Key = nFITS Int (1)\n")
            in_file.write(str(nFITS)+"\n")
            in_file.write("$Key = FITSdirs Str (128,"+str(nFITS)+")\n")
            for x in FITSdirs:
                in_file.write(x+"\n")

        for adverb in params.input_list:
            self.__write_adverb(params, in_file, adverb, input_dict[adverb], \
                                params.type_dict[adverb])

        in_file.close()

        # If debugging add a link to the input file to preserve it.
        if input_dict['DEBUG']:
            tmpDebug = tmpInput + 'Dbg'
            if os.access(tmpDebug, os.F_OK):
                os.unlink(tmpDebug)     # Remove any old version file.
            os.link(tmpInput, tmpDebug) # Add new link.
            # Tell about it.
            print "Saving copy of Obit task input in " + tmpDebug

        path = name
        if (not os.access(path,os.F_OK)) and  (os.getenv("OBITEXEC") != None):
            path = os.environ['OBITEXEC'] + '/bin/' + name
        if (not os.access(path,os.F_OK)) and (os.getenv("OBIT") != None):
            path = os.environ['OBIT'] +'/bin/' + name
        # Check OBITSD if needbe
        if  (not os.access(path,os.F_OK)) and (os.getenv('OBITSD') != None):
            path = os.getenv('OBITSD') +'/bin/' +  name
        
        # Check standard linux directory if needbe
        if  (not os.access(path,os.F_OK)):
            path = "/usr/lib/obit/bin/"+ name

        # Better have found it by here
        if  (not os.access(path,os.F_OK)):
            # Oh bugger
            msg = "Task '%s' executable not found" % (name)
            print msg
            raise RuntimeError, msg
        
        arglist = [name, "-input", tmpInput, "-output", tmpOutput,
                   "-pgmNumber", str(popsno), "-AIPSuser", str(userno)]
        tid = Task.spawn(self, path, arglist)
        self._params[tid] = params
        self._popsno[tid] = popsno
        return tid

    def messages(self, tid):
        """Return task's messages.
        
        Return a list of messages each as a tuple (1, message)
        tid   = Task id in pid table of process
        """

        # Add a default priority to the messages
        messages = Task.messages(self, tid)
        return [(1, msg) for msg in messages]
    
    def wait(self, tid):
        """Wait for the task to finish.
        
        Waits for Obit task to finishes, reads the task output
        parameter file and deleted tas input and output parameter
        file. 
        Returns output parameters in adictionary.
        tid   = Task id in pid table of process
        """

        assert(self.finished(tid))

        params = self._params[tid]
        popsno = self._popsno[tid]
        index = popsno - 1

        tmpInput = "/tmp/" + params.name + "Input." + str(popsno)
        tmpOutput = "/tmp/" + params.name + "Output." + str(popsno)

        output_dict = {}
        out_file = open(tmpOutput, mode='r')
        for adverb in params.output_list:
            # Need to parse whole file each time as order not specified
            output_dict[adverb] = self.__read_adverb(params, out_file, adverb)
        out_file.close()

        if os.access(tmpInput, os.F_OK):
            os.unlink(tmpInput)         # Remove input file.
        if os.access(tmpOutput, os.F_OK):
            os.unlink(tmpOutput)        # Remove output file.

        _free_popsno(popsno)  # Releast Pops number

        # Check if terminated normally
        retCode = output_dict["retCode"]
        if retCode != 0:
            msg = "Task '%s' returns '%d'" % (params.name, retCode)
            raise RuntimeError, msg

        del self._params[tid]
        del self._popsno[tid]
        Task.wait(self, tid)

        return output_dict



def _allocate_popsno():
    """ allocate pops number 

    In order to prevent multiple Obit instances from using the same POPS
    number, every Obit instance creates a lock file in /tmp.  These lock
    files are named Obitx.yyy, where x is the POPS number (in extended
    hex) and yyy is the process ID of the Obit instance.
    Creates a file /tmp/Obit_pops.pid to indicate this pops number is
    in use.
    """
    for popsno in range(1,16):
        # In order to prevent a race, first create a lock file for
        # POPSNO.
        try:
            path = '/tmp/Obit' + ehex(popsno, 1, 0) + '.' + str(os.getpid())
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0644)
            os.close(fd)
        except:
            continue

        # Get a list of likely lock files and iterate over them.
        # Leave out our own lock file though.
        files = glob.glob('/tmp/Obit' + ehex(popsno, 1, 0) + '.[0-9]*')
        files.remove(path)
        for file in files:
            #print "DEBUG try lock file ",file
            # If I can't write it it probably not mine
            if not os.access(file, os.W_OK):
                #print "DEBUG write access failed for ",file
                break
            # If the part after the dot isn't an integer, it's not a
            # proper lock file.
            try:
                pid = int(file.split('.')[1])
            except:
                continue

            # Check whether the Obit instance is still alive.
            try:
                os.kill(pid, 0)
            except:
                # The POPS number is no longer in use.  Try to clean
                # up the lock file.  This might fail though if we
                # don't own it.
                #print "DEBUG kill 0 failed for ",pid
                try:
                    os.unlink(file)
                except:
                    pass
            else:
                # The POPS number is in use.
                break
            
        else:
            # The POPS number is still free.
            #print "DEBUG assign popsno ",popsno
            return popsno

        # Clean up our own mess.
        os.unlink(path)

    raise RuntimeError, "No free Obit POPS number available on this system"

def _free_popsno(popsno):
    """ Deallocate pops number
    
    Unlinks /tmp/Obit_pops file
    """
    path = '/tmp/Obit' + ehex(popsno, 1, 0) + '.' + str(os.getpid())
    os.unlink(path)
