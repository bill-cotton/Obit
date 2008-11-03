# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2007 Associated Universities, Inc.
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

"""

This module provides the AIPSImage and AIPSUVData classes.  These
classes implement most of the data oriented verb-like functionality
from classic AIPS.

"""

# Global AIPS defaults.
from AIPS import AIPS

# Generic Python stuff.
import sys

# This code is way too clever.  Instead of implementing each and every
# function call provided by a proxy, class _Method implements a
# callable object that invokes a named method of a proxy, passing a
# description of the AIPS data it should operate on as the first
# argument.  To make matters worse, the same implementation is used
# for class AIPSImage and class AIPSData and dispatching is done based
# on the name of the class.  This way adding a function to the proxy
# automatically makes it callable as a method of both AIPSImage and
# AIPSUVData.

def _whoami():
    """Return the name of the function that called us."""
    return sys._getframe(1).f_code.co_name


class _AIPSDataMethod:

    """This class implements dispatching function calls to a proxy."""

    def __init__(self, inst, name):
        self.inst = inst
        self.name = name

    def __call__(self, *args):
        func = self.inst._method(self.name)
        return func(self.inst.desc, *args)


class _AIPSDataDesc:

    """This class implements the description of AIPS data that is used
       when dispatching function calls to a proxy."""

    def __init__(self, name, klass, disk, seq):
        self.userno = AIPS.userno
        self.name   = name
        self.klass  = klass
        self.disk   = disk
        self.seq    = seq

    # Provide a dictionary-like interface to deal with the
    # idiosyncrasies of XML-RPC.
    def __getitem__(self, key):
        return self.__dict__[key]


class _AIPSData:

    """This class describes generic AIPS data."""

    def __init__(self, name, klass, disk, seq):
        self.desc  = _AIPSDataDesc(name, klass, AIPS.disks[disk].disk, seq)
        self.proxy = AIPS.disks[disk].proxy()
        self.disk  = disk
        self.url   = AIPS.disks[disk].url
        self.dirs  = None  # AIPS data directories
        return

    name = property(lambda self: self.desc.name,
                    doc='Name of this data set.')
    klass = property(lambda self: self.desc.klass,
                     doc='Class of this data set.')
    disk = property(lambda self: self.desc.disk,
                    doc='Disk where this data set is stored.')
    seq = property(lambda self: self.desc.seq,
                   doc='Sequence number of this data set.')
    userno = property(lambda self: self.desc.userno,
                      doc='User number used to access this data set.')

    def __repr__(self):
        repr = "%s('%s', '%s', %d, %d)" % \
               (self.__class__.__name__,
                self.name, self.klass, self.disk, self.seq)
        return repr

    def __str__(self):
        return self.__repr__()

    def __getattr__(self, name):
        if name in self.desc.__dict__:
            return self.desc.__dict__[name]
        return _AIPSDataMethod(self, name)

    def table(self, type, version):
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.disk = thedisk
        return _AIPSTable(self, type, version, dirs)

    def _method(self, name):
        return getattr(getattr(self.proxy, self.__class__.__name__), name)

    def exists(self):
        """Check whether this image or data set exists.

        Returns True if the image or data set exists, False otherwise.
        """
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc)

    def verify(self):
        """Verify whether this image or data set can be accessed."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc)

    def header(self):
        """Get the header for this image or data set.

        Returns the header as a dictionary."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc)

    def clearstat(self, code=4):
        """Clear any AIPS status

        Returns the header as a dictionary."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc, code)

    def tables(self):
        """Get the list of extension tables."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc)

    def table_highver(self, type):
        """Get the highest version of an extension table.

        Returns the highest available version number of the extension
        table TYPE."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.dirs = dirs
        self.desc.disk = thedisk
        return self._method(_whoami())(self.desc, type)

    def zap(self):
        """Destroy this image or data set."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc)

    def header_table(self, type, version):
        """Get the header of an extension table.

        Returns the header of version VERSION of the extension table
        TYPE."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc, type, version)

    def getrow_table(self, type, version, rowno):
        """Get a row from an extension table.

        Returns row ROWNO from version VERSION of extension table TYPE
        as a dictionary."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc, type, version, rowno)

    def zap_table(self, type, version):
        """Destroy an extension table.

        Deletes version VERSION of the extension table TYPE.  If
        VERSION is 0, delete the highest version of table TYPE.  If
        VERSION is -1, delete all versions of table TYPE."""
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc, type, version)


class AIPSImage(_AIPSData):

    """This class describes an AIPS image."""
    def __init__(self, name, klass, disk, seq):
        self.desc    = _AIPSDataDesc(name, klass, AIPS.disks[disk].disk, seq)
        self.proxy   = AIPS.disks[disk].proxy()
        self.disk    = disk
        self.url     = AIPS.disks[disk].url
        self.dirs    = None  # AIPS data directories
        self.myClass = "AIPSImage"
        self.FileType= "AIPS"
        self.Aname   = name
        self.Aclass  = klass
        self.Aseq    = seq
        self.Disk    = disk
        self.Atype   = "MA"
        self.Otype   = "Image"
        return

    def display(self, dispURL="http://localhost:8765/RPC2"):
        """Display an image.

        Displays image on ObitView server on dispURL
        dispURL = URL of ObitView server on which to display
        Returns True if successful
        """
        (thedisk,dirs) = adjust_disk(self.disk, self.url)
        self.desc.disk = thedisk
        self.desc.dirs = dirs
        return self._method(_whoami())(self.desc, dispURL)



class AIPSUVData(_AIPSData):

    """This class describes an AIPS UV data set."""
    def __init__(self, name, klass, disk, seq):
        self.desc    = _AIPSDataDesc(name, klass, AIPS.disks[disk].disk, seq)
        self.proxy   = AIPS.disks[disk].proxy()
        self.disk    = disk
        self.url     = AIPS.disks[disk].url
        self.dirs    = None  # AIPS data directories
        self.myClass = "AIPSUVData"
        self.FileType= "AIPS"
        self.Aname   = name
        self.Aclass  = klass
        self.Aseq    = seq
        self.Disk    = disk
        self.Atype   = "UV"
        self.Otype   = "UV"
        return


class _AIPSTableMethod(_AIPSDataMethod):

    """ This class implements dispatching table oriented function
    calls to a proxy."""

    def __init__(self, inst, name):
        _AIPSDataMethod.__init__(self, inst, name)

    def __call__(self, *args):
        func = self.inst.data._method(self.name + '_table')
        return func(self.inst.data.desc,
                    self.inst.name, self.inst.version, *args)


class _AIPSTable:

    """This class describes a generic AIPS extension table."""

    def __init__(self, data, name, version):
        self.data = data
        self.name = name
        self.version = version

    def __getattr__(self, name):
        return _AIPSTableMethod(self, name)


class AIPSCat:
    def __init__(self, disk):
        proxy = AIPS.disks[disk].proxy()
        url   = AIPS.disks[disk].url
        (thedisk,dirs) = adjust_disk(disk, url)
        self.catalog = proxy.AIPSCat.cat(thedisk, AIPS.userno, url, dirs)
        return

    def __repr__(self):
        # Print something useful if the catalog is empty.
        if len(self.catalog) == 0:
            return 'Empty'

        return ''.join (['%d %s\n' % entry for entry in self.catalog]).strip()

    def info(self, slot):
        """ Returns information on catalog slot slot

        Returns None if not found
        slot = catalog slot number in AIPS catalog self
        """
        out = None
        if len(self.catalog) == 0:
            return out
        for entry in self.catalog:
            if slot==entry[0]:
                return entry[1]
        return out   # Failed to find
    # end info


def adjust_disk(disk, url):
    """Adjusts disk numbers and sets list of directories on url
    
    Returns (outdisk, dirs)
    where outdisk is the disk number of disk on url and dirs is
    the list of AIPS directories
    disk = Disk number to be converted
    url  = url of proxy, None = local
    """
    
    # AIPS data directories
    AIPSdirs = []
    for x in AIPS.disks:
        if x!=None and x.url==url:
            AIPSdirs.append(x.dirname)
    
    # Adjust disk number
    i = 1;
    outdisk = 0
    ldisk = int(disk)
    # Look for matching AIPS directory name
    for y in AIPSdirs:
        if AIPS.disks[ldisk] and y==AIPS.disks[ldisk].dirname:
            outdisk = i
            break
        i = i+1
    
    return (outdisk, AIPSdirs)
# end adjust_disk
