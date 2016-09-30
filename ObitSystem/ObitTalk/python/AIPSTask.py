"""

This module provides the AIPSTask class.  It adapts the Task class from
the Task module to be able to run classic AIPS tasks:

>>> imean = AIPSTask('imean')

The resulting class instance has all associated adverbs as attributes:

>>> print imean.ind
0.0
>>> imean.ind = 1
>>> print imean.indisk
1.0
>>> imean.indi = 2.0
>>> print imean.ind
2.0

It also knows the range for these attributes:

>>> imean.ind = -1
Traceback (most recent call last):
  ...
ValueError: value '-1.0' is out of range for attribute 'indisk'
>>> imean.ind = 10.0
Traceback (most recent call last):
  ...
ValueError: value '10.0' is out of range for attribute 'indisk'

>>> imean.inc = 'UVDATA'

>>> print imean.inclass
UVDATA

>>> imean.blc[1:] = [128, 128]
>>> print imean.blc
[None, 128.0, 128.0, 0.0, 0.0, 0.0, 0.0, 0.0]

>>> imean.blc = AIPSList([256, 256])
>>> print imean.blc
[None, 256.0, 256.0, 0.0, 0.0, 0.0, 0.0, 0.0]

It doesn't hurt to apply AIPSList to a scalar:
>>> AIPSList(1)
1

And it works on matrices (lists of lists) too:
>>> AIPSList([[1,2],[3,4],[5,6]])
[None, [None, 1, 2], [None, 3, 4], [None, 5, 6]]

It should also work for strings:
>>> AIPSList('foobar')
'foobar'
>>> AIPSList(['foo', 'bar'])
[None, 'foo', 'bar']
"""
# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2007,2011,2016 Associated Universities, Inc. Washington DC, USA.
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

# Generic Python stuff
import pydoc, select, fcntl, signal

# Global AIPS defaults.
from AIPS import AIPS

# Generic Task implementation.
from Task import Task, List

# Generic Python stuff.
import glob, os, pickle, sys

def count_entries(l):
    """Count number of last non blank/zero entries in list l"""
    count = len(l)
    c     = count
    if type(l[0])==int:  # integers
        for j in range(c-1,0,-1):
            if (l[j]==0):
                count -= 1
            else:
                break
        count = max(4,count) # At least 4
    elif type(l[0])==float:  # float
        for j in range(c-1,0,-1):
            if (l[j]==0.0):
                count -= 1
            else:
                break
        count = max(2,count)  # At least 2
                    
    elif type(l[0])==str:  # string
        for j in range(c-1,0,-1):
            if (l[j]==None) or ((len(l[j])>0) and (not l[j].isspace())):
                count -= 1
            else:
                break
        count = max(1,count)  # At least 1
    else:           # Something else
        count = len(l)

    count = min (count, len(l))
    # Trap AIPS lists
    if l[0]==None:
        count -= 1
    return count
#  end count_entries

class AIPSTask(Task):

    """This class implements running AIPS tasks.
    
    The AIPSTask class, handles client-side task related operations.
    Actual task definition and operations are handled by server-side
    proxies. For local operations, the server-side functionality is
    implemented in the same address space but remote operation is
    through an xmlrpc interface.  Tasks are run as separate processes
    in all cases.

    Each defined disk has an associated proxy, either local or remote.
    A proxy is a module with interface functions,
    local proxies are class modules from subdirectory Proxy with the
    same name (i.e. ObitTask) and the server functions are implemented
    there.  Remote proxies are specified by a URL and a proxy from the
    xmlrpclib module is used.

    When an object is created, the task secific parameters and
    documentation are retrieved by parsing the task Help file.and the
    POPSDAT.HLP file for parameter definitions.  This is performed on
    the server-side. 
    """

    # Package.
    _package = 'AIPS'

    # List of adverbs referring to data.
    _data_adverbs = ['indata', 'outdata',
                     'in2data', 'in3data', 'in4data', 'out2data']
    
    # List of adverbs referring to disks.
    _disk_adverbs = ['indisk', 'outdisk',
                     'in2disk', 'in3disk', 'in4disk', 'out2disk']
    
    # List of adverbs referring to file names.
    _file_adverbs = ['infile', 'infile2', 'outfile', 'outprint',
                     'ofmfile', 'boxfile', 'oboxfile']
    
    # Default version.
    version = os.environ.get('VERSION', 'NEW')

    # Default user number.
    userno = 0

    # Default verbosity level.
    msgkill = 0

    # Default to batch mode.
    isbatch = 32000

    # Run synchronous?
    doWait = False

    # Logging file?
    logFile = ""

    def __init__(self, name, **kwds):
        """  Create AIPS task object
        
        Creates task object and calls server function to parse the
        task help and POPSDAT.HLP files to obtain task specific
        parametrs and documentation.
        Following is a list of class members:
        _default_dict   = Dictionary with default values of parameters
        _input_list     = List of input parameters in order
        _output_list    = List of output parameters in order
        _min_dict       = Parameter minimum values as a List
        _max_dict       = Parameter maximum values as a List
        _hlp_dict       = Parameter descriptions (list of strings)
                          as a dictionary
        _strlen_dict    = String parameter lengths as dictionary
        _help_string    = Task Help documentation as list of strings
        _explain_string = Task Explain documentation as list of strings
        _short_help     = One line description of task
        _message_list   = list of execution messages
        Current parameter values are given as class members.
        """
        if not self._task_type:
            self._task_type = 'AIPS'
        Task.__init__(self)
        self._name = name
        self._input_list = []
        self._output_list = []
        self._message_list = []
        self._remainder = ""     # Partial message buffer
        # Optional arguments.
        if 'version' in kwds:
            self.version = kwds['version']
        else:
            if 'AIPS_VERSION' in os.environ:
                self.version = os.environ["AIPS_VERSION"]

        # Update default user number.
        if self.__class__.userno == 0:
            self.__class__.userno = AIPS.userno

        # See if there is a proxy that can hand us the details for
        # this task.
        params = None
        for proxy in AIPS.proxies:
            try:
                inst = getattr(proxy, self.__class__.__name__)
                params = inst.params(name, self.version)
            except Exception, exception:
                print exception
                if AIPS.debuglog:
                    print >>AIPS.debuglog, exception
                continue
            break
        if not params:
            msg = "%s task '%s' is not available" % (self._package, name)
            raise RuntimeError, msg

        # The XML-RPC proxy will return the details as a dictionary,
        # not a class.
        self._default_dict = params['default_dict']
        self._input_list = params['input_list']
        self._output_list = params['output_list']
        self._min_dict = params['min_dict']
        self._max_dict = params['max_dict']
        self._hlp_dict = params['hlp_dict']
        self._strlen_dict = params['strlen_dict']
        self._help_string = params['help_string']
        self._explain_string = params['explain_string']
        self._short_help = params['short_help']
        if self._task_type=='OBIT':
            self._type_dict = params['type_dict']
            self._dim_dict  = params['dim_dict']
        for adverb in self._default_dict:
            if type(self._default_dict[adverb]) == list:
                value = self._default_dict[adverb]
                self._default_dict[adverb] = List(self, adverb, value)

        # Initialize all adverbs to their default values.
        self.__dict__.update(self._default_dict)

        # The maximum value for disk numbers is bogus.
        for name in self._disk_adverbs:
            if name in self._max_dict:
                self._max_dict[name] = float(len(AIPS.disks) - 1)
                
        return                          # __init__

    def __eq__(self, other):
        """ Check if two task objects are for the same task """
        if self.__class__ != other.__class__:
            return False
        if self._name != other._name:
            return False
        if self.userno != other.userno:
            return False
        for adverb in self._input_list:
            if self.__dict__[adverb] != other.__dict__[adverb]:
                return False
            continue
        return True
    
    def copy(self):
        """ Return a copy of a given task object"""
        task = AIPSTask(self._name, version=self.version)
        task.userno = self.userno
        for adverb in self._input_list:
            task.__dict__[adverb] = self.__dict__[adverb]
            continue
        return task

    def defaults(self):
        """Set adverbs to their defaults."""
        self.__dict__.update(self._default_dict)
        
    def __display_adverbs(self, adverbs, file=None):
        """Display task ADVERBS values and descriptions"""
        inpList = self._short_help
        inpList = inpList + "Adverbs     Values                        "+\
                  "           Comments\n"
        inpList = inpList + "------------------------------------------"+\
                  "-----------------------------------------------------\n"

        for adverb in adverbs:
            #if self.__dict__[adverb] == '':
            #    print "'%s': ''" % adverb
            #else:
            #    value = PythonList(self.__dict__[adverb])
            #    print "'%s': %s" % (adverb, value)
            #    pass
            #continue
            i = 0
            j = 0
            hlps = self._hlp_dict[adverb]
            hlp = hlps[j]
            value = PythonList(self.__dict__[adverb])
            s1 = str(adverb)+"             "
            # Python is SO hateful
            if str(value).startswith('['):  # list
                # How many nonzero/nonblank entries
                ecount = count_entries(value)
                s2=""
                while i<ecount and 2+len(s2)+len(str(value[i])[:47])<50:
                    s2 = s2 + str(value[i])[:47]+", "
                    i = i+1
                # remove final comma
                if i==ecount:
                    s2 = s2[:len(s2)-2]
                # Left justify
                s2 = s2 + "                                         "
                inpList = inpList + "%12s%50s%s\n" % (s1[:11], s2[:49], hlp)
                # Loop until done with values
                doh = 0
                while i<ecount:
                    # continuation of description?
                    j = j+1
                    if j<len(hlps):
                        hlp =  hlps[j]
                    else:
                        hlp = " "
                    s2=""
                    while i<ecount and 2+len(s2)+len(str(value[i])[:47])<50:
                        s2 = s2 + str(value[i])[:47]+", "
                        i = i+1
                    # remove final comma
                    if i==ecount:
                        s2 = s2[:len(s2)-2]
                    # Left justify
                    s2 = s2 + "                                         "
                    inpList = inpList + "            %50s%s\n" % (s2[:49], hlp)
                    doh += 1;
                    if doh>ecount:
                        break
            else:  # Scalar
                s2 = str(value)+"                                   "
                inpList = inpList + "%12s%50s%s\n" % (s1[:11], s2[:49], hlp)
                    
            # Any more parameter description lines?
            s1 = "     "
            s2 = "     "
            j = j + 1
            while j<len(hlps):
                hlp =  hlps[j]
                j = j + 1
                inpList = inpList + "%12s%50s%s\n" % (s1, s2, hlp)

        if file:
            fd = open(file, "a")
            fd.write(inpList)
            fd.close()
        else:
            pydoc.ttypager(inpList)
        del inpList

    def inputs(self, file=None):
        """Display all inputs for this task."""
        self.__display_adverbs(self._input_list, file=file)

    def outputs(self, file=None):
        """Display all outputs for this task."""
        self.__display_adverbs(self._output_list, file=file)

    def _retype(self, value):
        """ Recursively transform a 'List' into a 'list' """

        if type(value) == List:
            value = list(value)
            for i in range(1, len(value)):
                value[i] = self._retype(value[i])

        return value

    def spawn(self):
        """Spawn the task.
         
        Writes task input parameters, task parameter file and starts
        the task asynchronously returning immediately.  Messages must be
        retrieved calling messages.
        Returns (proxy, tid)
        """

        if self.userno == 0:
            raise RuntimeError, "AIPS user number is not set"

        input_dict = {}
        for adverb in self._input_list:
            input_dict[adverb] = self._retype(self.__dict__[adverb])

        # Figure out what proxy to use for running the task, and
        # translate the related disk numbers.
        url = None
        proxy = None
        found = False
        for adverb in self._disk_adverbs:
            if adverb in input_dict:
                disk = int(input_dict[adverb])
                if disk == 0:
                    continue
                if not url and not proxy:
                    url = AIPS.disks[disk].url
                    proxy = AIPS.disks[disk].proxy()
                    found = True
                if AIPS.disks[disk].url != url:
                    raise RuntimeError, \
                          "AIPS disks are not on the same machine"
                input_dict[adverb] = float(AIPS.disks[disk].disk)
        if not found:
            raise RuntimeError, \
                  "Unable to determine where to execute task"

        # Adjust disks for proxy
        self.adjust_disk(input_dict, url)

        inst = getattr(proxy, self.__class__.__name__)
        tid = inst.spawn(self._name, self.version, self.userno,
                         self.msgkill, self.isbatch, input_dict)

        self._message_list = []
        return (proxy, tid)

    def finished(self, proxy, tid):
        """Determine if task has finished 
        
        Determine whether the task specified by PROXY and TID has
        finished.
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process
        """

        inst = getattr(proxy, self.__class__.__name__)
        return inst.finished(tid)

    def messages(self, proxy=None, tid=None):
        """Return task messages
        
        Returns list of messages and appends them to the object's
        message list.        
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process
       """

        # Bombs on remote callif not proxy and not tid:
        if not tid:
            return self._message_list

        inst = getattr(proxy, self.__class__.__name__)
        messbuff = inst.messages(tid)
        #print "MessBuff",messbuff
        # Parse messages into complete lines
        messages = self.parseMessage(messbuff)
        if not messages:
            return None
        for message in messages:
            self._message_list.append(message[1])
            if message[0] > abs(self.msgkill):
                #print message[1]
                pass
            continue
        return [message[1] for message in messages]

    def wait(self, proxy, tid):
        """Wait for the task to finish.
        
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process
       """

        while not self.finished(proxy, tid):
            pass
            #self.messages(proxy, tid)
        inst = getattr(proxy, self.__class__.__name__)
        output_dict = inst.wait(tid)
        for adverb in self._output_list:
            self.__dict__[adverb] = output_dict[adverb]
            continue
        return

    def feed(self, proxy, tid, banana):
        """Feed the task a  BANANA.

        Pass a message to a running task's sdtin
        proxy   = Proxy giving access to server
        tid     = Task id in pid table of process
        bananna = text message to pass to task input
        """
        
        inst = getattr(proxy, self.__class__.__name__)
        return inst.feed(tid, banana)

    def abort(self, proxy, tid, sig=signal.SIGTERM):
        """Abort the task specified by PROXY and TID.
        
        Calls abort function for task tid on proxy.
        None return value
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process to be terminated
        sig   = signal to sent to the task
        """

        inst = getattr(proxy, self.__class__.__name__)
        return inst.abort(tid, sig)

    def go(self):
        """Run the task.
        
        Writes task input parameters in the task parameter file and
        starts the task synchronously returning only when the task
        terminates. Messages are displayed as generated by the task,
        saved in an array returned from the call and, if the task
        member logFile is set, written to this file.
        """

        (proxy, tid) = self.spawn()
        #log = []
        count = 0
        rotator = ['|\b', '/\b', '-\b', '\\\b']
        # Logging to file?
        if len(self.logFile)>0:
            AIPS.log = file(self.logFile,'a')
        else:
            AIPS.log = None
        try:
            try:
                while not self.finished(proxy, tid):
                    messages = self.messages(proxy, tid)
                    if messages:
                        for message in messages:
                            print message
                            if AIPS.log:
                                if type(message)==str:
                                    x=AIPS.log.write('%s\n' % message)
                                else:
                                    x=AIPS.log.write('%s\n' % message[1])
                        #print "DEBUG message",messages
                        #log.extend(messages)
                        if AIPS.log:
                            AIPS.log.flush()
                    elif sys.stdout.isatty():
                        #sys.stdout.write(rotator[count % 4])
                        #sys.stdout.flush()
                        pass
                    events = select.select([sys.stdin.fileno()], [], [], 0)
                    if sys.stdin.fileno() in events[0]:
                        flags = fcntl.fcntl(sys.stdin.fileno(), fcntl.F_GETFL)
                        flags |= os.O_NONBLOCK
                        fcntl.fcntl(sys.stdin.fileno(), fcntl.F_SETFL, flags)
                        message = sys.stdin.read(1024)
                        flags &= ~os.O_NONBLOCK
                        fcntl.fcntl(sys.stdin.fileno(), fcntl.F_SETFL, flags)
                        self.feed(proxy, tid, message)
                        rotator = []
                        pass
                    count += 1
                    continue
                pass
            except KeyboardInterrupt, exception:
                self.abort(proxy, tid)
                raise exception

            self.wait(proxy, tid)
        finally:
            pass
        if AIPS.log:
            AIPS.log.close()
        return #log

    def __call__(self):
        return self.go()

    def __getattr__(self, name):
        if name in self._data_adverbs:
            class _AIPSData: pass
            value = _AIPSData()
            prefix = name.replace('data', '')
            value.name = Task.__getattr__(self, prefix + 'name')
            value.klass = Task.__getattr__(self, prefix + 'class')
            value.disk = Task.__getattr__(self, prefix + 'disk')
            value.seq = Task.__getattr__(self, prefix + 'seq')
            return value
        return Task.__getattr__(self, name)
    
    def __setattr__(self, name, value):
        if name in self._data_adverbs:
            prefix = name.replace('data', '')
            Task.__setattr__(self, prefix + 'name', value.name)
            Task.__setattr__(self, prefix + 'class', value.klass)
            Task.__setattr__(self, prefix + 'disk', value.disk)
            Task.__setattr__(self, prefix + 'seq', value.seq)
        else:
            # We treat 'infile', 'outfile' and 'outprint' special.
            # Instead of checking the length of the complete string,
            # we only check the length of the final component of the
            # pathname.  The backend will split of the direcrory
            # component and use that as an "area".
            attr = self._findattr(name)
            #file_adverbs = ['infile', 'outfile', 'outprint']
            if attr in self._file_adverbs and type(value) == str and \
                   os.path.dirname(value):
                if len(os.path.basename(value)) > self._strlen_dict[attr] - 2:
                    msg = "string '%s' is too long for attribute '%s'" \
                          % (value, attr)
                    raise ValueError, msg
                self.__dict__[attr] = value
            else:
                Task.__setattr__(self, name, value)
                pass
            pass
        return

    def adjust_disk(self, dict, url):
        """Adjusts disk numbers and sets list of disks for proxy

        Also converts lists to normal python lists
        dict = task input dictionary
        url  = url of proxy, None = local
        """

        # Fix lists
        for x in dict:
            if type(dict[x])== self.List or type(dict[x]) == list:
                tlist=[]
                for y in dict[x][1:]:
                    # List of AIPS lists?
                    if (type(y)== self.List or type(y) == list) and (y[0]==None):
                        y = y[1:]
                    tlist.append(y)
                dict[x] = tlist
        
        # AIPS data
        AIPSdirs = []
        for x in AIPS.disks:
            if x!=None and x.url==url:
                AIPSdirs.append(x.dirname)
        #
        # Save data directories
        dict["AIPSdirs"] = AIPSdirs

        # Adjust disk numbers, in AIPS everything is a float
        for x in self._disk_adverbs:
            if x in dict:
                diskno = int(dict[x])
                i = 1;
                # Look for matching AIPS directory name
                for y in AIPSdirs:
                    if AIPS.disks[diskno]:
                        if y==AIPS.disks[diskno].dirname:
                            dict[x] = float(i)
                            break
                    i = i+1
        # DEBUG
        #print "DEBUG",dict
        
        return

    pass                                # end class AIPSTask


class AIPSMessageLog:

    # Default user number.
    userno = -1

    def __init__(self):
        # Update default user number.
        if self.userno == -1:
            self.userno = AIPS.userno
        return

    def zap(self):
        """Zap message log."""

        proxy = AIPS.disks[1].proxy()
        inst = getattr(proxy, self.__class__.__name__)
        return inst.zap(self.userno)

    pass                                # class AIPSMessageLog


def AIPSList(list):
    """Transform a Python array into an AIPS array.

    Returns a list suitable for using 1-based indices.
    """

    try:
        # Make sure we don't consider strings to be lists.
        if str(list) == list:
            return list
    except:
        pass

    try:
        # Insert 'None' at index zero, and transform LIST's elements.
        _list = [None]
        for l in list:
            _list.append(AIPSList(l))
            continue
        return _list
    except:
        # Apparently LIST isn't a list; simply return it unchanged.
        return list

def PythonList(list):
    """Transform an AIPS array into a Python array.

    Returns a list suitable for using normal 0-based indices.
    """

    try:
        if list[0] != None:
            return list

        _list = []
        for l in list[1:]:
            _list.append(PythonList(l))
            continue
        return _list
    except:
        # Apparently LIST isn't a list; simply return it unchanged.
        return list


# Tests.
if __name__ == '__main__':
    import doctest, sys
    results = doctest.testmod(sys.modules[__name__])
    sys.exit(results[0])
