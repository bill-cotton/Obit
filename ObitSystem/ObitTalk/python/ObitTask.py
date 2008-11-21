"""

This module provides the ObitTask class.  It adapts the Task class from
the Task module to be able to run Obit tasks:

>>> imean = ObitTask('Template')

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

"""
# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2005,2007 Associated Universities, Inc. Washington DC, USA.
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


import LocalProxy
# AIPSTask implementation.
from AIPSTask import AIPSTask
from AIPS import AIPS
from FITS import FITS

# Generic Task implementation.
from Task import Task, List

# Generic Python stuff.
import glob, os, pickle, sys, signal

class ObitTask(AIPSTask):

    """This class implements running Obit tasks.

    The ObitTask class, derived from the AIPSTask class, handles
    client-side task related operations.  Actual task definition
    and operations are handled by server-side proxies.
    For local operations, the server-side functionality is implemented
    in the same address space but remote operation is through an
    xmlrpc interface.  Tasks are run as separate processes in all
    cases.

    Each defined disk has an associated proxy, either local or remote.
    A proxy is a module with interface functions,
    local proxies are class modules from subdirectory Proxy with the
    same name (i.e. AIPSTask) and the server functions are implemented
    there.  Remote proxies are specified by a URL and a proxy from the
    xmlrpclib module is used.

    When an object is created, the task secific parameters and
    documentation are retrieved by parsing the task TDF file.
    This is performed on the server-side.
    """

    # Package.
    _package = 'Obit'

    # List of adverbs referring to disks.
    _disk_adverbs = ['inDisk', 'in2Disk', 'in3Disk', 'in4Disk', 'outDisk', 'out2Disk']
    _indisk_adverbs = ['inDisk', 'in2Disk', 'in3Disk', 'in4Disk']
    _outdisk_adverbs = ['outDisk', 'out2Disk']

    # Default version.
    version = 'OBIT'

    # Debugging?
    debug = False

    # Run synchronous?
    doWait = False

    # Logging file?
    logFile = ""

    def __init__(self, name):
        """  Create Obit task object
        
        Creates task object and calls server function to parse TDF
        file to obtain task specific parametrs and documentation.
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
        retCode         = Task return code, 0=Finished OK
        debug           = If true save task parameter file
        logFile         = if given, the name of the file in which to
                          write messages
        doWait          = True if synchronous  operation wished
        Current parameter values are given as class members.
        """
        AIPSTask.__init__(self, name)
        self._remainder = ""   # Partial message buffer
        if self.userno == 0:
            self.userno = 1

        # Set disk directory names
        #print "DEBUG ObitTask._init__ before ",FITS.disks
        dirs = []
        for x in FITS.disks:
            if x!=None:
                #print "DEBUG ObitTask._init__ ",x.dirname
                dirs.append(x.dirname)
        #print "DEBUG ObitTask._init__ before ",FITS.disks
        self.__dict__["FITSdirs"] = dirs
        #print "DEBUG ObitTask._init__ after ",self.__dict__["FITSdirs"]
        dirs = []
        for x in AIPS.disks:
            if x!=None:
                #print "DEBUG ObitTask._init__ ",x.dirname
                dirs.append(x.dirname)
        self.__dict__["AIPSdirs"] = dirs
        #print "DEBUG ObitTask._init__ after ",self.__dict__["AIPSdirs"]

    def abort(self, proxy, tid, sig=signal.SIGKILL):
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
        
        Writes task input parameters, data directories and other
        information in the task parameter file and starts the task
        synchronously returning only when the task terminates.
        Messages are displayed as generated by the task,
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
            while not self.finished(proxy, tid):
                messages = self.messages(proxy, tid)
                if messages:
                    for message in messages:
                        #log.append(message)
                        print message[1]
                        if AIPS.log:
                            if type(message)==str:
                                x=AIPS.log.write('%s\n' % message)
                            else:
                                x=AIPS.log.write('%s\n' % message[1])
                        continue
                    pass
                if AIPS.log:
                    AIPS.log.flush()
                elif sys.stdout.isatty():
                    # Boo-hiss sys.stdout.write(rotator[count % 4])
                    # sys.stdout.flush()
                    pass
                count += 1
                continue
        except KeyboardInterrupt, exception:
            #self.abort(proxy, tid, sig=signal.SIGKILL)
            self.abort(proxy, tid)
            raise exception

        self.wait(proxy, tid)
        if AIPS.log:
            AIPS.log.close()
        return #log

    def spawn(self):
        """Spawn the task.
        
        Writes task input parameters, data directories and other
        information in the task parameter file and starts the task
        asynchronously returning immediately.  Messages must be
        retrieved calling messages.
        Returns (proxy, tid)
        """

        if self.userno == 0:
            raise RuntimeError, "AIPS user number is not set"

        self.__dict__["retCode"] = -1  # init return code
        
        input_dict = {}
        for adverb in self._input_list:
            input_dict[adverb] = self._retype(self.__dict__[adverb])

        # Debugging?
        input_dict["DEBUG"] = self.debug

        # Figure out what proxy to use for running the task, and
        # translate the related disk numbers.
        url = None
        found = False
        allLocal = True # Local FITS files?
        proxy = None
        for adverb in self._disk_adverbs:
            if adverb in input_dict:
                disk = int(input_dict[adverb])
                if disk == 0:
                    continue
                if not url and not proxy:
                    if self.__dict__['DataType'] == 'AIPS':
                        url = AIPS.disks[disk].url
                        proxy = AIPS.disks[disk].proxy()
                        found = True
                        if AIPS.disks[disk].url != url:
                            raise RuntimeError, \
                                  "AIPS disks are not on the same machine"
                        input_dict[adverb] = int(AIPS.disks[disk].disk)
                    elif self.__dict__['DataType'] == 'FITS':
                        url = FITS.disks[disk].url
                        proxy = FITS.disks[disk].proxy()
                        found = True
                        allLocal = allLocal and (url==None) # Local?
                        if FITS.disks[disk].url != url:
                            raise RuntimeError, \
                                  "FITS disks are not on the same machine"
                        input_dict[adverb] = int(FITS.disks[disk].disk)
        # If FITS and all disks=0, run locally
        if (self.__dict__['DataType'] == 'FITS') and allLocal:
            found = True
            url   = None
            proxy = LocalProxy
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

    def messages(self, proxy=None, tid=None):
        """Return messages for the task specified by PROXY and TID.
        
        Returns list of messages and appends them to the object's
        message list.
        proxy = Proxy giving access to server
        tid   = Task id in pid table of process
        """

        # Bombs on remote calling (proxy!=None) and not tid:
        if not tid:
            return self._message_list

        inst = getattr(proxy, self.__class__.__name__)
        messbuff = inst.messages(tid)
        # Parse messages into complete lines
        messages = self.parseMessage(messbuff)
        if not messages:
            return None
        self._message_list.extend(messages)
        return messages

    def adjust_disk(self, dict, url):
        """Adjusts disk numbers and sets list of disks for proxy
        
        dict = task input dictionary
        url  = url of proxy, None = local
        """

        # AIPS data
        AIPSdirs = []
        for x in AIPS.disks:
            if x!=None and x.url==url:
                AIPSdirs.append(x.dirname)
        #
        # FITS data 
        FITSdirs = []
        for x in FITS.disks:
            if x!=None and x.url==url:
                FITSdirs.append(x.dirname)
        #
        
        # Save data directories
        dict["FITSdirs"] = FITSdirs
        dict["AIPSdirs"] = AIPSdirs
        
        # Adjust disk numbers
        DataType = dict["DataType"]
        for x in self._disk_adverbs:
            # Input or output data
            if x in self._indisk_adverbs:
                DataType = dict["DataType"]
            if x in self._outdisk_adverbs:
                if "outDType" in dict:
                    DataType = dict["outDType"]
            if x in dict:
                diskno = dict[x]
                if DataType=="AIPS":
                    i = 1;
                    # Look for matching AIPS directory name
                    for y in AIPSdirs:
                        if AIPS.disks[diskno]:
                            if y==AIPS.disks[diskno].dirname:
                                dict[x] = i
                                break
                        i = i+1
                elif DataType=="FITS":
                    i = 0;
                    # Look for matching FITS directory name
                    if FITS.disks[diskno]:
                        for y in FITSdirs:
                            if y==FITS.disks[diskno].dirname:
                                if dict[x]>0: # Don't translate disk 0 => CWD
                                    dict[x] = i
                                break
                            i = i+1
        # DEBUG
        #print "DEBUG dict",dict
        
        return

# Tests.
if __name__ == '__main__':
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
