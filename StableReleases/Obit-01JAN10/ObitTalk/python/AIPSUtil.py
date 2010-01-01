# Copyright (C) 2005 Joint Institute for VLBI in Europe
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

This module provides some utility functions for dealing with AIPS.

"""

def ehex(n, width=0, padding=None):
    """Convert a number into "extended hex".

    Returns the extended hex presentation for N, optionally padding it
    up to WIDTH with PADDING."""

    ehex_digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    result = ''

    while n > 0:
        result = ehex_digits[n % len(ehex_digits)] + result
        n /= len(ehex_digits)
        width -= 1
        continue

    # Pad if requested to do so.
    if padding != None:
        while width > 0:
            result = str(padding) + result
            width -= 1
            continue
        pass

    return result
