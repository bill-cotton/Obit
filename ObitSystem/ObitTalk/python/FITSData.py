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

    def __init__(self, filename, disk):
        self.filename = filename
        self.disk = disk

    # Provide a dictionary-like interface to deal with the
    # idiosyncrasies of XML-RPC.
    def __getitem__(self, key):
        return self.__dict__[key]


class _FITSData:

    """This class describes generic FITS data."""

    def __init__(self, name, disk):
        self.desc = _FITSDataDesc(name, FITS.disks[disk].disk)
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

    def __getattr__(self, filename):
        if filename in self.desc.__dict__:
            return self.desc.__dict__[filename]
        return _FITSDataMethod(self, filename)

    def table(self, type, version):
        return _FITSTable(self, type, version)

    def _method(self, name):
        return getattr(getattr(self.proxy), filename)

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
    pass


class FITSUVData(_FITSData):

    """This class describes an FITS UV data set."""
    pass


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
