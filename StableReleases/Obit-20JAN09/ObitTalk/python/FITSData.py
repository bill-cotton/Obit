# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2005 Associated Universities, Inc. Washington DC, USA.
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

This module provides the FITSImage and FITSUVData classes.  These
classes implement most of the data oriented verb-like functionality
from classic FITS.

"""

# Global FITS defaults.
from FITS import FITS

# Generic Python stuff.
import sys
import LocalProxy

# This code is way too clever.  Instead of implementing each and every
# function call provided by a proxy, class _Method implements a
# callable object that invokes a named method of a proxy, passing a
# description of the FITS data it should operate on as the first
# argument.  To make matters worse, the same implementation is used
# for class FITSImage and class FITSData and dispatching is done based
# on the name of the class.  This way adding a function to the proxy
# automatically makes it callable as a method of both FITSImage and
# FITSUVData.

def _whoami():
    """Return the name of the function that called us."""
    return sys._getframe(1).f_code.co_name


class _FITSDataMethod:

    """This class implements dispatching function calls to a proxy."""

    def __init__(self, inst, name):
        self.inst = inst
        self.name = name

    def __call__(self, *args):
        func = self.inst._method(self.name)
        return func(self.inst.desc, *args)


class _FITSDataDesc:

    """This class implements the description of FITS data that is used
       when dispatching function calls to a proxy."""

    def __init__(self, filename, dirname):
        self.filename = filename
        self.dirname  = dirname

    # Provide a dictionary-like interface to deal with the
    # idiosyncrasies of XML-RPC.
    def __getitem__(self, key):
        return self.__dict__[key]


class _FITSData:

    """This class describes generic FITS data."""

    def __init__(self, name, disk):
        self.desc = _FITSDataDesc(name, FITS.disks[disk].dirname)
        self.proxy = FITS.disks[disk].proxy()
        self.disk = disk
        return

    filename = property(lambda self: self.desc.filename,
                    doc='Filename of this data set.')
    disk = property(lambda self: self.desc.disk,
                    doc='Disk where this data set is stored.')

    def __repr__(self):
        repr = "%s('%s', %d)" % \
               (self.filename, self.disk)
        return repr

    def __str__(self):
        return self.__repr__()

    def __getattr__(self, name):
        if name in self.desc.__dict__:
            return self.desc.__dict__[name]
        return _FITSDataMethod(self, name)

    def table(self, type, version):
        return _FITSTable(self, type, version)

    def _method(self, name):
        return getattr(getattr(self.proxy, self.__class__.__name__), name)

    def exists(self):
        """Check whether this image or data set exists.

        Returns True if the image or data set exists, False otherwise."""
        return self._method(_whoami())(self.desc)

    def verify(self):
        """Verify whether this image or data set can be accessed."""
        return self._method(_whoami())(self.desc)

    def header(self):
        """Get the header for this image or data set.

        Returns the header as a dictionary."""
        return self._method(_whoami())(self.desc)

    def tables(self):
        """Get the list of extension tables."""
        return self._method(_whoami())(self.desc)

    def table_highver(self, type):
        """Get the highest version of an extension table.

        Returns the highest available version number of the extension
        table TYPE."""
        return self._method(_whoami())(self.desc, type)

    def zap(self):
        """Destroy this image or data set."""
        return self._method(_whoami())(self.desc)

    def header_table(self, type, version):
        """Get the header of an extension table.

        Returns the header of version VERSION of the extension table
        TYPE."""
        return self._method(_whoami())(self.desc, type, version)

    def getrow_table(self, type, version, rowno):
        """Get a row from an extension table.

        Returns row ROWNO from version VERSION of extension table TYPE
        as a dictionary."""
        return self._method(_whoami())(self.desc, type, version, rowno)

    def zap_table(self, type, version):
        """Destroy an extension table.

        Deletes version VERSION of the extension table TYPE.  If
        VERSION is 0, delete the highest version of table TYPE.  If
        VERSION is -1, delete all versions of table TYPE."""
        return self._method(_whoami())(self.desc, type, version)


class FITSImage(_FITSData):

    """This class describes an FITS image."""
    def __init__(self, filename, disk):
        if disk==0:
            proxy   = LocalProxy
            dirname = "./"
        else:
            proxy = FITS.disks[disk].proxy()
            dirname = FITS.disks[disk].dirname
        self.desc     = _FITSDataDesc(filename, dirname)
        self.proxy    = proxy
        self.filename = filename
        self.FileName = filename
        self.Fname    = filename
        self.disk     = disk
        self.Disk     = disk
        self.myClass  = "FITSImage"
        self.FileType = "FITS"
        self.Otype    = "Image"
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


class FITSUVData(_FITSData):

    """This class describes an FITS UV data set."""
    def __init__(self, filename, disk):
        if disk==0:
            proxy   = LocalProxy
            dirname = "./"
        else:
            proxy = FITS.disks[disk].proxy()
            dirname = FITS.disks[disk].dirname
        self.desc     = _FITSDataDesc(filename, dirname)
        self.proxy    = FITS.disks[disk].proxy()
        self.filename = filename
        self.FileName = filename
        self.Fname    = filename
        self.disk     = disk
        self.Disk     = disk
        self.myClass  = "FITSUVData"
        self.FileType = "FITS"
        self.Otype    = "UV"
        return


class _FITSTableMethod(_FITSDataMethod):

    """ This class implements dispatching table oriented function
    calls to a proxy."""

    def __init__(self, inst, name):
        _FITSDataMethod.__init__(self, inst, name)

    def __call__(self, *args):
        func = self.inst.data._method(self.name + '_table')
        return func(self.inst.data.desc,
                    self.inst.name, self.inst.version, *args)


class _FITSTable:

    """This class describes a generic FITS extension table."""

    def __init__(self, data, name, version):
        self.data = data
        self.name = name
        self.version = version

    def __getattr__(self, name):
        return _FITSTableMethod(self, name)

class FITSCat:
    def __init__(self, disk, dir=None):
        if disk==0:
            proxy = LocalProxy
            url   = " "
        else:
            proxy = FITS.disks[disk].proxy()
            url   = FITS.disks[disk].url
        userno = 1
        (thedisk,dirs) = adjust_disk(disk, url)
        self.catalog = proxy.FITSCat.cat(thedisk,dir, userno, url, dirs)
        return

    def __repr__(self):
        # Print something useful if the catalog is empty.
        if len(self.catalog) == 0:
            return 'Empty'

        return ''.join (['%d %s\n' % entry for entry in self.catalog]).strip()

    def info(self, slot):
        """ Returns information on catalog slot slot

        Returns None if not found
        slot = catalog slot number in FITS catalog self
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
    the list of FITS directories
    disk = Disk number to be converted
    url  = url of proxy, None = local
    """
    
    # FITS data directories
    FITSdirs = []
    for x in FITS.disks:
        if x!=None and x.url==url:
            FITSdirs.append(x.dirname)
    
    # Adjust disk number
    i = 1;
    outdisk = 0
    ldisk = int(disk)
    # Look for matching FITS directory name
    for y in FITSdirs:
        if FITS.disks[ldisk] and y==FITS.disks[ldisk].dirname:
            outdisk = i
            break
        i = i+1
    
    return (outdisk, FITSdirs)
# end adjust_disk
