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

This module provides the AIPSDisk and AIPS classes.  Together they
provide some basic infrastructure used by the AIPSTask and AIPSData
modules.

"""

# Generic Python stuff.
import os
from AIPSUtil import *

# Available proxies.
import LocalProxy
from xmlrpclib import ServerProxy


class AIPSDisk:
    
    """Class representing a (possibly remote) AIPS disk.  An instance
       of this class stores an AIPS disk number and the URL of the
       proxy through which it can be accessed.  For local AIPS disks
       the URL will be None."""

    def __init__(self, url, disk, dirname):
        self.url = url
        self.disk = disk
        self.dirname = dirname

    def proxy(self):
        """Return the proxy through which this AIPS disk can be
           accessed."""
        if self.url:
            return ServerProxy(self.url)
        else:
            return LocalProxy


class AIPS:
    
    """Container for several AIPS-related default values."""

    # Default AIPS user ID.
    userno = 0

    # List of available proxies.
    proxies = [ LocalProxy ]

    # AIPS disk mapping. 
    disks = [ None ]                    # Disk numbers are one-based.

    # AIPS seems to support a maximum of 35 disks.
    for disk in xrange(1, 35):
        area = 'DA' + ehex(disk, 2, '0')
        dirname = os.getenv(area)
        if not area in os.environ:
            break
        disks.append(AIPSDisk(None, disk, dirname))
        continue

    # Message log.
    log = None

    # Debug log.
    debuglog = None

