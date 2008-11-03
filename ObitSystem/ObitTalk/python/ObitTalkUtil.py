"""
This module provides some utility functions for dealing with ObitTalk.

"""
# Copyright (C) 2008 Associated Universities, Inc.
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
import os, pydoc, string
# Global AIPS defaults
from AIPS import AIPS
from AIPS import AIPSDisk
import AIPSUtil
# Global FITS defaults
from FITS import FITS
from FITS import FITSDisk
import AIPSDir, FITSDir
import OErr
err = OErr.OErr()

def SetEnviron(AIPS_ROOT=None, AIPS_VERSION=None, OBIT_EXEC=None, \
               DA00=None, ARCH="LINUX", aipsdirs=None, fitsdirs=None):
    """ Set environment variables for AIPS, Obit and data access

    AIPS_ROOT    = name of the directory at the root of the AIPS directories
    AIPS_VERSION = name of the directory at root of specific AIPS version
    OBIT_EXEC    = name of directory with Obit bin and TDF directories
    DA00         = name of directory for AIPS DA00
    ARCH         = AIPS architecture string
    aipsdirs     = list of tuples defining AIPS data directories
                   first element is url, None = local,
                       for remote data use form 'http://mymachine.org:8000/RPC2'
                   second element is full path to AIPS data directory
    fitsdirs     = list of tuples defining FITS data directories
                   first element is url, None = local
                   second element is full path to FITS directory
    """
    # Set AIPS_ROOT and AIPS_VERSION for access to AIPS Software
    if AIPS_ROOT:
        os.environ["AIPS_ROOT"] = AIPS_ROOT  # python environment
        cmd = "export AIPS_ROOT="+AIPS_ROOT  # shell environment
        #print cmd
        z=os.system(cmd)
    if AIPS_VERSION:
        os.environ["AIPS_VERSION"] = AIPS_VERSION  # python environment
        cmd = "export AIPS_VERSION="+AIPS_VERSION  # shell environment
        #print cmd
        z=os.system(cmd)

    # Set Obit directory
    if OBIT_EXEC:
        os.environ["OBITEXEC"] = OBIT_EXEC  # python environment
        cmd = "export OBITEXEC="+OBIT_EXEC   # shell environment
        #print cmd
        z=os.system(cmd)

    # Set AIPS DA00
    if DA00:
        os.environ["DA00"] = DA00   # python environment
        cmd = "export DA00="+DA00   # shell environment
        #print cmd
        z=os.system(cmd)

    # Set AIPS DA00
    if ARCH:
        os.environ["ARCH"] = ARCH   # python environment
        cmd = "export ARCH="+ARCH   # shell environment
        #print cmd
        z=os.system(cmd)

    # Set AIPS directories
    if aipsdirs:
        for dirname in aipsdirs:
            url  = dirname[0]
            # Already known?
            known = False
            for adisk in AIPS.disks[1:]:
                if (adisk.url==url) and (adisk.dirname==dirname[1]):
                    known = True
                    break
            if not known:
                # New local or remote AIPS disk
                disk = len(AIPS.disks)
                AIPS.disks.append(AIPSDisk(url, disk, dirname[1]+'/'))
                AIPSDir.PSetDir(disk, dirname[1]+'/', err, URL=url)
                # Define environment variable for local disks
                if url==None:
                    dskno = "DA"+AIPSUtil.ehex(disk,width=2,padding='0')
                    os.environ[dskno] = dirname[1]+'/'        # Python environment
                    cmd = "export "+dskno+"="+dirname[1]+'/'  # shell environment
                    #print cmd
                    z=os.system(cmd)

    # Set NVOL so AIPS Tasks will know how many AIPS disks
    naips = len(AIPS.disks)-1
    os.environ['NVOL'] = str(naips)   # Python evironment

    # Set FITS directories (URL, disk name)
    import FITSDir  # Need Obit/python
    if fitsdirs:
        from FITS import FITSDisk
        # First is $FITS
        if len(fitsdirs)>0:
            cmd = "export FITS="+fitsdirs[0][1]+'/'
            z=os.system(cmd)                         # shell environment
            os.environ["FITS"] = fitsdirs[0][1]+'/'  # Python environment
            # Override previous?
            if len(FITS.disks)>0:
                FITS.disks[1].dirname = fitsdirs[0][1]+'/'
                FITSDir.FITSdisks[0] = fitsdirs[0][1]
                disk = 0
                FITSDir.PSetDir(fitsdirs[0][1]+'/', disk, err, URL=url)
            else:
                url = fitsdirs[0][0]
                FITSDir.PAddDir(fitsdirs[0][1]+'/', err, URL=url)
            # Other FITS directories
            for dirname in fitsdirs[1:]:
                url  = dirname[0]
                # Already known?
                known = False
                for fdisk in FITS.disks[1:]:
                    if (fdisk.url==url) and (fdisk.dirname==dirname[1]):
                        known = True
                        break
                if not known:
                    # New local or remote AIPS disk
                    disk = len(FITS.disks)
                    FITSDir.PAddDir(dirname[1]+'/', err, URL=url)
                    # Define environment variable for local disks
                    if url==None:
                        dskno = "FITS"+AIPSUtil.ehex(disk-1,width=2,padding='0')
                        os.environ[dskno] = dirname[1]+'/'           # Python environment
                        cmd = "export "+dskno+"="+dirname[1]+'/'     # shell environment
                        #print cmd
                        z=os.system(cmd)
    # end SetEnviron

def ListAIPSDirs():
    """ List currently defined AIPS directories
    
    """
    dirlist = "\nAIPS Directories \n"
    diskno = 1
    for adisk in AIPS.disks[1:]:
        if adisk.url:
            url = adisk.url
        else:
            url = "localhost"
        dirlist =  dirlist+string.rjust(str(diskno),3)+"  "+adisk.dirname+ \
                  ", URL="+url+"\n"
        diskno += 1
    # User pager
    pydoc.ttypager(dirlist)
    # end  ListAIPSDirs

def ListFITSDirs():
    """ List currently defined FITS directories
    
    """
    dirlist = "\nFITS Directories \n"
    diskno = 1
    for fdisk in FITS.disks[1:]:
        if fdisk.url:
            url = fdisk.url
        else:
            url = "localhost"
        dirlist =  dirlist+string.rjust(str(diskno),3)+"  "+fdisk.dirname+ \
                  ", URL="+url+"\n"
        diskno += 1
    # User pager
    pydoc.ttypager(dirlist)
    # end  ListFITSDirs
