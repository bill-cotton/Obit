"""

This module provides the bits and pieces to implement an AIPSTask
proxy object.

"""
# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2006,2007 Associated Universities, Inc. Washington DC, USA.
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


# Global AIPS defaults.
from Proxy.AIPS import AIPS, ehex

# The results from parsing POPSDAT.HLP.
from Proxy.Popsdat import Popsdat

# Bits from the generic Task implementation.
from Proxy.Task import Task

# Generic Python stuff.
import glob, os, pickle, signal, struct, string

class _AIPSTaskParams:
    """ AIPS Task parameter class """
    def __parse(self, name):
        """Determine the proper attributes for the AIPS task NAME by
        parsing its HELP file. and POPSDAT.HLP file"""

        # Pretend we know nothing yet.
        task = None
        adverb = None
        short_help = None

        popsdat = Popsdat(self.version)

        path = self.version + '/HELP/' + name.upper() + '.HLP'
        input = open(path)

        # Parse INPUTS section.
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

            # Empty lines start with '\n'.
            if line.startswith('\n'):
                continue

            # Short task description
            if task and line.startswith(task+" "):
                # get short help string
                if  not short_help:
                    short_help = line
                continue

            # Continuation lines start with ' '.
            if line.startswith(' '):
                # Continuation of parameter description?
                if adverb and len(line)>max_end:
                    self.hlp_dict[adverb].append(string.rstrip(line[max_end:]))
                continue

            if not task:
                min_start = line.find('LLLLLLLLLLLL')
                min_end = line.rfind('L')
                max_start = line.find('UUUUUUUUUUUU')
                max_end = line.rfind('U')
                dir_start = min_start - 2
                dir_end = min_start - 1
                if not min_start == -1 and not max_start == -1:
                    task = line.split()[0]
                continue

            adverb = line.split()[0].lower()
            code = line[min_start - 1:min_start]
            if not code:
                code = ' '
            try:
                min = float(line[min_start:min_end])
                max = float(line[max_start:max_end])
            except:
                min = None
                max = None

            # Get help text from end of line
            hlp = [string.rstrip(line[max_end:])]
 
            match_key = None
            if adverb in popsdat.default_dict:
                 match_key = adverb
            else:
                # Some HELP files contain typos.
                for key in popsdat.default_dict:
                    if key.startswith(adverb):
                        if match_key:
                            msg = "adverb '%s' is ambiguous" % adverb
                            raise AttributeError, msg
                        else:
                            match_key = key
            if not match_key:
                match_key = key
            self.default_dict[adverb] = popsdat.default_dict[match_key]

            if code in ' *&$':
                self.input_list.append(adverb)
            if code in '&%$@':
                self.output_list.append(adverb)
            if adverb in popsdat.strlen_dict:
                self.strlen_dict[adverb] = popsdat.strlen_dict[adverb]
            if min != None:
                self.min_dict[adverb] = min
            if max != None:
                self.max_dict[adverb] = max
            if hlp != None:
                self.hlp_dict[adverb] = hlp

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
            self.explain_string = self.explain_string + line

    def __init__(self, name, version):
        """ Create Server AIPS  task structure
        
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
        self.short_help = ''
        self.explain_string = ''
        self.help_string = ''

        self.name = name
        if version in os.environ:
            self.version = os.environ[version]
        else:
            self.version = os.environ['AIPS_ROOT'] + '/' + version
            pass

        # Bad juju
        #path = os.environ['HOME'] + '/.ObitTalk/' \
        #       + os.path.basename(self.version) + '/' \
        #       + name.lower() + '.pickle'
        #try:
        #    unpickler = pickle.Unpickler(open(path))
        #    self.default_dict = unpickler.load()
        #    self.input_list = unpickler.load()
        #    self.output_list = unpickler.load()
        #    self.min_dict = unpickler.load()
        #    self.max_dict = unpickler.load()
        #    self.strlen_dict = unpickler.load()
        #    self.hlp_dict = unpickler.load()
        #    self.help_string = unpickler.load()
        #    self.explain_string = unpickler.load()
        #    self.short_help = unpickler.load()
        #except (IOError, EOFError):

        # Always parse
        self.__parse(name)

        #    # Make sure the directory exists.
        #    if not os.path.exists(os.path.dirname(path)):
        #        os.makedirs(os.path.dirname(path))
        
        #    pickler = pickle.Pickler(open(path, mode='w'))
        #    pickler.dump(self.default_dict)
        #    pickler.dump(self.input_list)
        #    pickler.dump(self.output_list)
        #    pickler.dump(self.min_dict)
        #    pickler.dump(self.max_dict)
        #    pickler.dump(self.strlen_dict)
        #    pickler.dump(self.hlp_dict)
        #    pickler.dump(self.help_string)
        #    pickler.dump(self.explain_string)
        #    pickler.dump(self.short_help)
        
    # Provide a dictionary-like interface to deal with the
    # idiosyncrasies of XML-RPC.
    def __getitem__(self, key):
        return self.__dict__[key]


class AIPSTask(Task):
    """ Server-side AIPS task interface """
    # List of adverbs referring to file names.
    _file_adverbs = ['infile', 'infile2', 'outfile', 'outprint',
                     'ofmfile', 'boxfile', 'oboxfile']
    
    def __init__(self):
        Task.__init__(self)
        self._params = {}
        self._popsno = {}
        self._userno = {}
        self._msgno = {}
        self._msgkill = {}
        self.POPSpid = os.getpid()

    def params(self, name, version):
        """Return parameter set for version VERSION of task NAME."""
        return _AIPSTaskParams(name, version)

    def __write_adverb(self, params, file, adverb, value):
        """Write (sub)value VALUE of adverb ADVERB into TD file FILE."""

        assert(adverb in params.input_list)

        if type(value) == float:
            file.write(struct.pack('f', value))
        elif type(value) == str:
            #print "DEBUG params.strlen_dict", adverb, params.strlen_dict[adverb]
            strlen = ((params.strlen_dict[adverb] + 3) // 4) * 4
            fmt = "%ds" % strlen
            file.write(struct.pack(fmt, value.ljust(strlen)))
        elif type(value) == list:
            for subvalue in value:
                self.__write_adverb(params, file, adverb, subvalue)
        else:
            print "offending adverb",adverb,", value",value,", type",type(value)
            raise AssertionError, type(value)

    def __read_adverb(self, params, file, adverb, value=None):
        """Read (sub)value for adverb ADVERB from TD file FILE."""

        assert(adverb in params.output_list)

        # We use the default value for type checks.
        if value == None:
            value = params.default_dict[adverb]

        if type(value) == float:
            (value,) = struct.unpack('f', file.read(4))
        elif type(value) == str:
            strlen = ((params.strlen_dict[adverb] + 3) // 4) * 4
            fmt = "%ds" % strlen
            (value,) = struct.unpack(fmt, file.read(strlen))
            value.strip()
        elif type(value) == list:
            newvalue = [None]
            for subvalue in value[1:]:
                subvalue = self.__read_adverb(params, file, adverb, subvalue)
                newvalue.append(subvalue)
                continue
            value = newvalue
        else:
            raise AssertionError, type(value)
        return value

    def spawn(self, name, version, userno, msgkill, isbatch, input_dict):
        """Start the task.
        
        Writes task input parameters in theTD file and starts the task
        asynchronously returning immediately.
        Messages must be  retrieved calling messages.
        Attempts to use singel hardcoded AIPS TV
        name        task name
        version     version of task
        userno      AIPS user number
        msgkill     AIPS msgkill level,
        isbatch     True if this is a batch process
        input_dict  Input parameters as dictionary
        Returns task id
        """

        params = _AIPSTaskParams(name, version)
        popsno = _allocate_popsno()
        index = popsno - 1

        try:
            # A single hardcoded TV will do until support for multiple
            # TVs is implemented.
            ntvdev = 1

            # Construct the environment for the task.  For the 'infile',
            # 'outfile' and 'outprint' adverbs, we split off the directory
            # component of the pathname and use that as the area.
            env = os.environ.copy()
            area = 'a'
            for adverb in self._file_adverbs:
                if adverb in input_dict:
                    assert(ord(area) <= ord('z'))
                    dirname = os.path.dirname(input_dict[adverb])
                    if dirname:
                        if not os.path.isdir(dirname):
                            msg = "Directory '%s' does not exist" % dirname
                            raise RuntimeError, msg
                        env[area] = dirname
                        basename = os.path.basename(input_dict[adverb])
                        input_dict[adverb] = area + ':' + basename
                        area = chr(ord(area) + 1)
                        pass
                    pass
                continue
            # Send output to the TV running on this machine.
            env['TVDEV' + ehex(ntvdev, 2, 0)] = 'sssin:localhost'
            # DEBUG
            # print "aips environment",env
            
            td_name = os.environ['DA00'] + '/TD' + AIPS.revision + '000004;'
            td_file = open(td_name, mode='r+b')
            
            td_file.seek(index * 20)
            td_file.write(struct.pack('8s', name.upper().ljust(8)))
            td_file.write(struct.pack('l', -999))
            td_file.write(struct.pack('2l', 0, 0))
            
            td_file.seek(1024 + index * 4096)
            td_file.write(struct.pack('i', userno))
            td_file.write(struct.pack('i', ntvdev))
            td_file.write(struct.pack('i', 0))
            td_file.write(struct.pack('i', msgkill + 32000 - 1))
            td_file.write(struct.pack('i', isbatch))
            td_file.write(struct.pack('i', 0))
            td_file.write(struct.pack('2i', 0, 0))
            td_file.write(struct.pack('f', 1.0))
            td_file.write(struct.pack('4s', '    '))
            for adverb in params.input_list:
                self.__write_adverb(params, td_file, adverb,
                                    input_dict[adverb])
                continue

            td_file.close()

            # Create the message file if necessary and record the
            # number of messages currently in it.
            user = ehex(userno, 3, 0)
            ms_name = os.environ['DA01'] + '/MS' + AIPS.revision \
                      + user + '000.' + user + ';'
            if not os.path.exists(ms_name):
                ms_file = open(ms_name, mode='w')
                ms_file.truncate(1024)
                ms_file.close()
                os.chmod(ms_name, 0664)
                pass
            ms_file = open(ms_name, mode='r')
            (msgno,) = struct.unpack('i', ms_file.read(4))
            ms_file.close()

            path = params.version + '/' + os.environ['ARCH'] + '/LOAD/' \
                   + name.upper() + ".EXE"
            tid = Task.spawn(self, path, [name.upper() + str(popsno)], env)
            
        except Exception, exception:
            _free_popsno(popsno)
            raise exception
        
        self._params[tid] = params
        self._popsno[tid] = popsno
        self._userno[tid] = userno
        self._msgkill[tid] = msgkill
        self._msgno[tid] = msgno
        return tid

    def __read_message(self, file, msgno):
        """ Read a message from an AIPS message file

        Returns (task, popsno, priority, message)
        file  = AIPS message file name
        msgno = Which messahe (0-rel)
        """
        file.seek((msgno / 10) * 1024 + 8 + (msgno % 10) * 100)
        (tmp, task, message) = struct.unpack('i8x5s3x80s', file.read(100))
        (popsno, priority) = (tmp / 16, tmp % 16)
        task = task.rstrip()
        message = message.rstrip()
        return (task, popsno, priority, message)

    def messages(self, tid):
        """Return task's messages.
        
        Return a list of messages each as a tuple (1, message)
        tid   = Task id in pid table of process
        """

        # Make sure we read the messages, even if we throw them away
        # later to prevent the task from blocking.
        messages = Task.messages(self, tid)
        #print "Raw",messages

        # Check that tid in bounds
        #if tid>len(self._popsno):
        #    print "DEBUG Proxy/AIPSTask", messages, tid, len(self._popsno)
        #    return messages

        # Convert into list of lines
        tmessage = []
        for msg in messages:
            lmess = msg.splitlines(True)
            for mm in lmess:
                tmessage.append(mm)

        # Strip out all formal messages.
        start = '%-5s%d' % (self._params[tid].name.upper(), self._popsno[tid])
        lmessages = [msg for msg in tmessage if msg.startswith(start)]
        #print "Filtered",lmessages
        pmessages = [msg for msg in tmessage if not msg.startswith(start)]
        #print "Print",pmessages

        # These messages will be looked up in the AIPS message log
        #messages = [(1, msg) for msg in lmessages]
        messages = []

        user = ehex(self._userno[tid], 3, 0)
        ms_name = os.environ['DA01'] + '/MS' + AIPS.revision \
                  + user + '000.' + user + ';'
        ms_file = open(ms_name, mode='r')

        (msgno,) = struct.unpack('i', ms_file.read(4))
        while self._msgno[tid] < msgno:
            (task, popsno, priority, msg) = \
                   self.__read_message(ms_file, self._msgno[tid])
            # Filter
            if popsno == self._popsno[tid]:
                messages.append((priority, '%-5s%d: %s\n' % (task, popsno, msg)))
                pass
            self._msgno[tid] += 1
            continue

        ms_file.close()
        # Add "print" messages
        if len(pmessages)>0:
            for msg in pmessages:
                messages.append((1,msg))
        #print "returned messages",messages
        return messages

    def wait(self, tid):
        """Wait for the task to finish.

        When task returns, the output parameters are parser from the
        TD file.
        Returns output parameters in adictionary.
        tid   = Task id in pid table of process
        """

        assert(self.finished(tid))

        params = self._params[tid]
        popsno = self._popsno[tid]
        index = popsno - 1

        td_name = os.environ['DA00'] + '/TDD000004;'
        try:
            td_file = open(td_name, mode='rb')
            
            td_file.seek(index * 20 + 8)
            (result,) = struct.unpack('i', td_file.read(4))
            if result != 0:
                msg = "Task '%s' returns '%d'" % (params.name, result)
                raise RuntimeError, msg
            
            td_file.seek(1024 + index * 4096 + 40)
            output_dict = {}
            for adverb in params.output_list:
                output_dict[adverb] = self.__read_adverb(params, td_file, adverb)
                
            td_file.close()

        finally:
            _free_popsno(popsno)
            pass

        del self._params[tid]
        del self._popsno[tid]
        del self._userno[tid]
        del self._msgno[tid]
        Task.wait(self, tid)

        return output_dict

    def abort(self, tid, sig=signal.SIGTERM):
        """Abort the task specified by PROXY and TID.
        
        Calls abort function for task tid on proxy.
        None return value
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process to be terminated
        sig   = signal to sent to the task
                AIPS seems to ignore SIGINT, so use SIGTERM instead.
        """

        _free_popsno(self._popsno[tid])

        del self._params[tid]
        del self._popsno[tid]
        del self._userno[tid]
        del self._msgno[tid]

        return Task.abort(self, tid, sig)

    pass                                # class AIPSTask

class AIPSMessageLog:
    def __init__(self):
        return

    def _open(self, userno):
        user = ehex(userno, 3, 0)
        ms_name = os.environ['DA01'] + '/MS' + AIPS.revision \
                  + user + '000.' + user + ';'
        return open(ms_name, mode='r+')

    def zap(self, userno):
        """Zap message log."""

        ms_file = self._open(userno)
        ms_file.write(struct.pack('i', 0))
        return True                # Return something other than None.

    pass                                # class AIPSMessageLog



def _allocate_popsno():
    """ allocate pops number 

    In order to prevent multiple AIPS instances from using the same POPS
    number, every AIPS instance creates a lock file in /tmp.  These lock
    files are named AIPSx.yyy, where x is the POPS number (in extended
    hex) and yyy is the process ID of the AIPS instance.
    Create a file in /tmp/AIPS_pops.pid to indicate this pops number is
    in use.
    """
    for popsno in range(1,16):
        # In order to prevent a race, first create a lock file for
        # POPSNO.
        try:
            path = '/tmp/AIPS' + ehex(popsno, 1, 0) + '.' + str(os.getpid())
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0666)
            os.close(fd)
        except:
            continue

        # Get a list of likely lock files and iterate over them.
        # Leave out our own lock file though.
        files = glob.glob('/tmp/AIPS' + ehex(popsno, 1, 0) + '.[0-9]*')
        files.remove(path)
        for file in files:
            # If the part after the dot isn't an integer, it's not a
            # proper lock file.
            try:
                pid = int(file.split('.')[1])
            except:
                continue

            # Check whether the AIPS instance is still alive.
            try:
                os.kill(pid, 0)
            except:
                # The POPS number is no longer in use.  Try to clean
                # up the lock file.  This might fail though if we
                # don't own it.
                try:
                    os.unlink(file)
                except:
                    pass
            else:
                # The POPS number is in use.
                break
        else:
            # The POPS number is still free.
            return popsno

        # Clean up our own mess.
        os.unlink(path)

    raise RuntimeError, "No free AIPS POPS number available on this system"

def _free_popsno(popsno):
    """ Deallocate pops number
    
    Unlinks /tmp/AIPS_pops file
    """
    path = '/tmp/AIPS' + ehex(popsno, 1, 0) + '.' + str(os.getpid())
    os.unlink(path)


