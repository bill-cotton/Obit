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

This module provides the FITSDisk and FITS classes.  

"""

# Generic Python stuff.
import os

# Available proxies.
import LocalProxy
from xmlrpclib import ServerProxy


class FITSDisk:
    
    """Class representing a (possibly remote) FITS disk.  An instance
       of this class stores an FITS disk number and the URL of the
       proxy through which it can be accessed.  For local FITS disks
       the URL will be None."""

    def __init__(self, url, disk, dirname):
        self.url = url
        self.disk = disk
        self.dirname = dirname

    def proxy(self):
        """Return the proxy through which this FITS disk can be
           accessed."""
        if self.url:
            return ServerProxy(self.url)
        else:
            return LocalProxy


class FITS:
    
    """Container for several FITS-related default values.
    Uses environment variables: FITS, FITS01, FITS02...
    """

    # List of available proxies.
    proxies = [ LocalProxy ]

    # FITS disk mapping. 
    disks = [ None ]                    # Disk numbers are one-based.

    # Look for FITS
    area = 'FITS'
    if area in os.environ:
        dirname = os.getenv(area)
        disks.append(FITSDisk(None, 1, dirname))

    # FITS01, FITS02...
    # Who will ever need more than 19 FITS disks?    
    for disk in xrange(2, 20):
        area = 'FITS%02d' % (disk-1)
        if not area in os.environ:
            break
        dirname = os.getenv(area)
        disks.append(FITSDisk(None, disk, dirname))
        continue

    # Message log.
    log = None
