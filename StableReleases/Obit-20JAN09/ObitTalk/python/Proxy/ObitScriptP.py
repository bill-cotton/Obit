"""

This module provides the bits and pieces to implement an ObitScript
proxy object.

"""
# Copyright (C) 2007 Associated Universities, Inc. Washington DC, USA.
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


# Global ObitScript defaults.
from Proxy.AIPS import AIPS, ehex

# Bits from the generic Task implementation.
from Proxy.Task import Task

# Generic Python stuff.
import glob, os, signal, struct, string

class ObitScript(Task):
    """ Server-side ObitScript script interface """
    
    def __init__(self):
        Task.__init__(self)
        self._popsno = {}
        self._userno = {}
        self._msgno = {}
        self._msgkill = {}
        self.POPSpid = os.getpid()

    def __write_script(self, name, userno, popsno, file, in_dict):
        """Write script text into a text file
        
        Writes Obit initialization and shutdown into script
        """
        # Open script file
        sc_file=open(file, mode="w")
        AIPSDirs = in_dict["AIPSDirs"]
        nAIPS = len(AIPSDirs)
        FITSDirs = in_dict["FITSDirs"]
        nFITS = len(FITSDirs)

        # Obit initialization
        sc_file.write("import Obit, OErr, OSystem, UV, Image, Table, InfoList, ODisplay\n")
        sc_file.write("from OTObit import *\n")
        sc_file.write("if OSystem.PIsInit():\n")
        sc_file.write("    OSystem.Shutdown() # Remove any prior\n")
        sc_file.write("# Initialize Obit\n")
        sc_file.write("err=OErr.OErr()\n")
        sc_file.write("ObitSys=OSystem.OSystem ('"+name+"', "+str(popsno)+", "+str(userno)+",\n")
        sc_file.write("    "+str(nAIPS)+", "+str(AIPSDirs)+", \n")
        sc_file.write("    "+str(nFITS)+", "+str(FITSDirs)+", True, False, err)\n")
        sc_file.write("OErr.printErrMsg(err, 'Error with Obit startup')\n")
        sc_file.write("AIPS.AIPS.userno="+str(userno)+"\n\n")
        
        #  Add script text
        for x in in_dict["script"]:
            if not x.endswith("\n"):  # Make sure terminated
                x += "\n"
            sc_file.write(x)
        
        # Obit shutdown
        sc_file.write("\n# Shutdown Obit\n")
        sc_file.write("OErr.printErr(err)\n")
        sc_file.write("OSystem.Shutdown(ObitSys)\n")

        sc_file.close()
        # end __write_script

    def spawn(self, name, version, userno, msgkill, isbatch, input_dict):
        """Start the script.in an externak ObitTalk
        
        Writes script test into temporary file in /tmp and executes ObitTalk
        asynchronously returning immediately.
        Messages must be  retrieved calling messages.
        name        script name
        version     version of any AIPS tasks
        userno      AIPS user number
        msgkill     AIPStask msgkill level,
        isbatch     True if this is a batch process
        input_dict  Input info as dictionary
        Returns script id
        """

        popsno = _allocate_popsno()
        index = popsno - 1
        self.name = name
        self.userno = userno

        try:
            # Construct the environment for the script.  
            # Set AIPS and FITS directory variables
            env = os.environ.copy()
            # AIPS directories
            i = 0
            for dirname in input_dict["AIPSDirs"]:
                i = i+1
                area = 'DA' + ehex(i, 2, '0')
                env[area] = dirname
            # end AIPS dirs
            # FITS directories
            i = 0
            for dirname in input_dict["FITSDirs"]:
                i = i+1
                if i==1:
                    area = "FITS"
                else:
                    area = 'FITS' + ehex(i, 2, '0')
                env[area] = dirname
            # End FITS Dirs
            
            # print "script environment",env
            # print "PYTHONPATH",env["PYTHONPATH"]

            # Script text file
            sc_name = "/tmp/" + name+ "Script." + str(popsno)+".py"
            
            # Write script to file
            self.__write_script (self.name, userno, popsno, sc_name, input_dict)

            # If debugging add a link to the input file to preserve it.
            if input_dict['debug']:
                tmpDebug = "/tmp/" + name+ "Script." + str(popsno)+"Dbg.py"
                if os.access(tmpDebug, os.F_OK):
                    os.unlink(tmpDebug)        # Remove any old version file.
                os.link(sc_name, tmpDebug) # Add new link.
                # Tell about it.
                print "Saving copy of Obit task input in " + tmpDebug
            
            # Start script in separate ObitTalk process.
            path = 'ObitTalk'
            tid = Task.spawn(self, path, ["ObitTalk", sc_name], env)
            
        except Exception, exception:
            _free_popsno(popsno)
            raise exception
        
        self._popsno[tid]  = popsno
        self._userno[tid]  = userno
        self._msgkill[tid] = msgkill
        self._msgno[tid]   = 0

        return tid

    def messages(self, tid):
        """Return script's messages.
        
        Return a list of messages each as a tuple (1, message)
        tid   = Script id in pid table of process
        """

        # Make sure we read the messages, even if we throw them away
        # later to prevent the script from blocking.
        messages = Task.messages(self, tid)

        return [(1, msg) for msg in messages]

    def wait(self, tid):
        """Wait for the script to finish.

        Waits until script is finished
        tid   = Script id in pid table of process
        """

        assert(self.finished(tid))
        
        # Cleanup
        popsno = self._popsno[tid]
        _free_popsno(self._popsno[tid])
        del self._popsno[tid]
        del self._userno[tid]
        del self._msgno[tid]
        # wait
        Task.wait(self, tid)
        
        # Delete Script text file
        sc_name = "/tmp/" + self.name+ "Script." + str(popsno)+".py"
        if os.access(sc_name, os.F_OK):
            os.unlink(sc_name)
        
        return True  # Return something other than None

    def abort(self, tid, sig=signal.SIGTERM):
        """Abort the script specified by PROXY and TID.
        
        Calls abort function for script tid on proxy.
        No return value
        proxy = Proxy giving access to server
        tid   = Script id in pid table of process to be terminated
        sig   = signal to sent to the script
                ObitScript seems to ignore SIGINT, so use SIGTERM instead.
        """

        _free_popsno(self._popsno[tid])

        del self._popsno[tid]
        del self._userno[tid]
        del self._msgno[tid]

        return Task.abort(self, tid, sig)

    pass                                # class ObitScript

def _allocate_popsno():
    """ allocate pops number 

    In order to prevent multiple ObitScript instances from using the same POPS
    number, every ObitScript instance creates a lock file in /tmp.  These lock
    files are named ObitScriptx.yyy, where x is the POPS number (in extended
    hex) and yyy is the process ID of the ObitScript instance.
    Create a file in /tmp/ObitScript_pops.pid to indicate this pops number is
    in use.
    """
    for popsno in range(1,16):
        # In order to prevent a race, first create a lock file for
        # POPSNO.
        try:
            path = '/tmp/ObitScript' + ehex(popsno, 1, 0) + '.' + str(os.getpid())
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0666)
            os.close(fd)
        except:
            continue

        # Get a list of likely lock files and iterate over them.
        # Leave out our own lock file though.
        files = glob.glob('/tmp/ObitScript' + ehex(popsno, 1, 0) + '.[0-9]*')
        files.remove(path)
        for file in files:
            # If the part after the dot isn't an integer, it's not a
            # proper lock file.
            try:
                pid = int(file.split('.')[1])
            except:
                continue

            # Check whether the ObitScript instance is still alive.
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

    raise RuntimeError, "No free ObitScript POPS number available on this system"

def _free_popsno(popsno):
    """ Deallocate pops number
    
    Unlinks /tmp/ObitScript_pops file
    """
    path = '/tmp/ObitScript' + ehex(popsno, 1, 0) + '.' + str(os.getpid())
    os.unlink(path)


