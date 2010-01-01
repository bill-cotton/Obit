"""

This module provides the ObitScript class.
This class allows running Obit/python scripts either
locally or remotely

ObitScripts are derived from Task and share most of execution properties.
In particular, ObitScripts can be executed either locally or remotely.
In this context a script is a character string containing a sequence of
ObitTalk or other python commands and may be included when the script
object is created or attached later.
An example:
script="import OSystem\nprint 'Welcome user',OSystem.PGetAIPSuser()\n"
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

# Generic Python stuff
import pydoc, select, fcntl, signal
import glob, os, pickle, sys

# Available proxies.
import LocalProxy
from xmlrpclib import ServerProxy

# Global AIPS and FITS defaults.
from Task import Task
from AIPS import AIPS
from FITS import FITS

class ObitScript(Task):

    """This class implements running Obit/python Script
    
    The ObitScript class, handles client-side script related operations.
    Actual script operations are handled by server-side proxies.
    For local operations, the server-side functionality is
    implemented in the same address space but remote operation is
    through an xmlrpc interface.  

    An ObitScript has an associated proxy, either local or remote.
    A proxy is a module with interface functions,
    local proxies are class modules from subdirectory Proxy with the
    same name (i.e. ObitScript) and the server functions are implemented
    there.  Remote proxies are specified by a URL and a proxy from the
    xmlrpclib module is used.

    """

    # Package.
    _package = 'ObitScript'

    # Default version.
    version = os.environ.get('VERSION', 'NEW')

    # Default user number.
    userno = AIPS.userno

    # Default to batch mode.
    isbatch = 32000

    # Run synchronous?
    doWait = False

    # Logging file?
    logFile = ""

    # AIPS directories
    AIPSDirs = []

    # FITS directories
    FITSDirs = []

    # Script
    script = []

    # Debugging
    debug = False

    # Default verbosity level.
    msgkill = 0

    def __init__(self, name, **kwds):
        """  Create ObitScript task object
        
        Creates Script Object.
        name  = name of script object
        Optional Keywords:
            script   = Script to execute as string or list of strings
            file     = Name of text file containing script
            URL      = URL on which the script is to be executed
                       Default = None = local execution
            AIPSDirs = List of AIPS directories on URL
                       Default = current AIPS directories on url
            FITSDirs = List of FITS directories on URL
                       Default = current FITS directories on url
            AIPSUser = AIPS user number for AIPS data files
                       Default is current
            version  = AIPS version string, Default = current
        Following is a list of class members:
            url      = URL of execution server, None=Local
            proxy    = Proxy for URL
            script   = Script as text string
            userno   = AIPS user number
            AIPSDirs = List of AIPS directories on URL
            FITSDirs = List of FITS directories on URL
            AIPSUser = AIPS user number for AIPS data files
            version  = AIPS version string
            _message_list = messages from Script execution
        """
        Task.__init__(self)
        
        self._name = name
        self.url = None         # URL - local
        self.proxy = LocalProxy  # proxy
        self._message_list = []
        self.AIPSDirs = []
        self.FITSDirs = []
        self.script   = []
        self._remainder = ""   # Partial message buffer
        
        # Optional arguments.
        if 'URL' in kwds:
            self.url = kwds['URL']
            self.proxy = ServerProxy(self.url)
        elif 'url' in kwds:
            self.url = kwds['url']
            self.proxy = ServerProxy(self.url)
        else:
            self.url = None
            self.proxy = LocalProxy
        
        # AIPS directories on target host
        self.AIPSDirs= []
        if 'AIPSDirs' in kwds:
            for x in  kwds['AIPSDirs']:
                self.AIPSDirs.append(x)
        else:   # Current disks on self.url
            for x in AIPS.disks:
                if x!=None and x.url==self.url:
                    self.AIPSDirs.append(x.dirname)
        
        # FITS directories on target host
        self.FITSDirs= []
        if 'FITSDirs' in kwds:
            for x in  kwds['FITSDirs']:
                self.FITSDirs.append(x)
        else:   # Current disks on self.url
            for x in FITS.disks:
                if x!=None and x.url==self.url:
                    self.FITSDirs.append(x.dirname)
        
        if 'AIPSUser' in kwds:
            self.userno = kwds['AIPSUser']
        else:
            self.userno = AIPS.userno

        # Script given
        if 'script' in kwds:
            scr =  kwds['script']
            if type(scr) == list:  # List of strings
                self.script = scr
            else: # Simple string
                self.script.append(scr)

        # File name given
        if 'file' in kwds:
            input = open(kwds['file'])
            line = " "
            while (line):
                line = input.readline() # read next line
                if not line:  # EOF?
                    break
                self.script.append(line)
            input.close()
            
        # Update default user number.
        if self.__class__.userno == 0:
            self.__class__.userno = AIPS.userno

        return                          # __init__

    def spawn(self):
        """Spawn the script.
         
        Starts script asynchronously returning immediately
        Messages must be retrieved calling messages.
        Returns (proxy, tid)
        """

        if not self.script:
            raise RuntimeError, "ObitScript: no script given"

        if self.userno == 0:
            raise RuntimeError, "ObitScript user number is not set"

        # Convert directory lists to list
        AIPSDirs = []
        for x in self.AIPSDirs:
            AIPSDirs.append(x)
        FITSDirs = []
        for x in self.FITSDirs:
            FITSDirs.append(x)
        
        # Package info as a dict for xmrrpc
        input_dict = {"script":self.script, "debug":self.debug, \
                      "AIPSDirs":AIPSDirs, "FITSDirs":FITSDirs}

        # Use proxy to start
        inst = getattr( self.proxy, self.__class__.__name__)
        tid = inst.spawn(self._name, self.version, self.userno,
                         self.msgkill, self.isbatch, input_dict)

        self._message_list = []
        return (self.proxy, tid)

    def finished(self, proxy, tid):
        """Determine if script has finished 
        
        Determine whether the script specified by PROXY and TID has
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

        # Bombs on remote call if not proxy and not tid:
        if not tid:
            return self._message_list

        inst = getattr(proxy, self.__class__.__name__)
        messbuff = inst.messages(tid)
        # Parse messages into complete lines
        messages = self.parseMessage(messbuff)
        if not messages:
            return None
        for message in messages:
            self._message_list.append(message[1])
             #    if message[0] > abs(self.msgkill):
             #        #print message[1]
             #        pass
             #    continue
        return [message[1] for message in messages]

    def wait(self, proxy, tid):
        """Wait for the script to finish.
        
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process
       """

        while not self.finished(proxy, tid):
            pass
            #self.messages(proxy, tid)
        inst = getattr(proxy, self.__class__.__name__)
        inst.wait(tid)
        return

    def feed(self, proxy, tid, banana):
        """Feed the script a  BANANA.

        Pass a message to a running script's sdtin
        proxy   = Proxy giving access to server
        tid     = Script task id in pid table of process
        bananna = text message to pass to script input
        """
        
        inst = getattr(proxy, self.__class__.__name__)
        return inst.feed(tid, banana)

    def abort(self, proxy, tid, sig=signal.SIGTERM):
        """Abort the script specified by PROXY and TID.
        
        Calls abort function for task tid on proxy.
        None return value
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process to be terminated
        sig   = signal to sent to the task
        """

        inst = getattr(proxy, self.__class__.__name__)
        return inst.abort(tid, sig)

    def go(self):
        """Execute the script.
        
        Writes task input parameters in the task parameter file and
        starts the task synchronously returning only when the task
        terminates. Messages are displayed as generated by the task,
        and, if the task member logFile is set, written to this file.
        """

        (proxy, tid) = self.spawn()
        log = []
        count = 0
        # Logging to file?
        if len(self.logFile)>0:
            ObitScript.log = file(self.logFile,'a')
        else:
            ObitScript.log = None
        try:
            try:
                while not self.finished(proxy, tid):
                    messages = self.messages(proxy, tid)
                    if messages:
                        for message in messages:
                            print message
                            if ObitScript.log:
                                if type(message)==str:
                                    x=ObitScript.log.write('%s\n' % message)
                                else:
                                    x=ObitScript.log.write('%s\n' % message[1])
                        #print "DEBUG message",messages
                        log.extend(messages)
                        if ObitScript.log:
                            ObitScript.log.flush()
                    elif sys.stdout.isatty():
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
        if ObitScript.log:
            ObitScript.log.close()
        return log

    def inputs(self):
        """List script"""
        scrList = "Listing of script "+self._name+"\n"
        if self.url:  # Remote host?
            scrList += "Execution host: "+self.url+"\n"
        for x in self.script:
            if not x.endswith("\n"):  # Make sure terminated
                x += "\n"
            scrList += x
        pydoc.ttypager(scrList)
        del scrList

    def outputs(self):
        """Not defined."""
        print "output not defined for a script"

    def help(self):
        """List script."""
        self.inputs()

    def explain(self):
        """List script"""
        self.inputs()

    def __call__(self):
        return self.go()

    def __getattr__(self, name):
        return Task.__getattr__(self, name)
    
    def __setattr__(self, name, value):
        self.__dict__[name] = value
        return

class ObitScriptMessageLog:

    # Default user number.
    userno = -1

    def __init__(self):
        # Update default user number.
        if self.userno == -1:
            self.userno = ObitScript.userno
        return

    def zap(self):
        """Zap message log."""

        proxy = ObitScript.disks[1].proxy()
        inst = getattr(proxy, self.__class__.__name__)
        return inst.zap(self.userno)

    pass                                # class ObitScriptMessageLog



# Tests.
if __name__ == '__main__':
    import doctest, sys
    results = doctest.testmod(sys.modules[__name__])
    sys.exit(results[0])

