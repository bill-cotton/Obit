"""
This module provides support for starting an interactive ObitTalk
session.

"""
# Copyright (C) 2006-2010 Associated Universities, Inc.
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

import os
# Global AIPS defaults
from AIPS import AIPS
# Global FITS defaults
from FITS import FITS

# The main classes ObitTalk provides.
from AIPSTask import *
from AIPSData import *
from ObitTask import *
from FITSData import *
from ObitScript import *
 
# Use our own, somewhat restricted rlcompleter.
import readline, otcompleter
readline.parse_and_bind("tab: complete")

# This is not a batch job.
AIPSTask.isbatch = 0

# Separate the blurb below from what the Python interpreter spits out.
print ""

# Is there any startup script in $HOME/.obitrc.py or ./.obitrc.py
path = ".obitrc.py"
if (not os.path.exists(path)) and ("HOME" in os.environ):
    path = os.environ["HOME"]+"/.obitrc.py"

print "Welcome to ObitTalk"
# Using AIPS?
if (("AIPS_ROOT" in os.environ) and ("AIPS_VERSION" in os.environ)) \
       or (os.path.exists(path)):
    while True:
        try:
            input = raw_input("Please enter your AIPS user ID number: ")
            AIPS.userno = int(input)
        except KeyboardInterrupt:
            print ""
            print "AIPS user ID number is not set"
            break
        except:
            print "That is not a valid AIPS user ID number"
            continue
        else:
            break
# End AIPS statrup

# Override help() such that it prints something useful for instances
# of AIPSTask.
_help = help
def help(obj):
    if isinstance(obj, AIPSTask):
        obj.help()
    else:
        _help(obj)
        pass
    return

# Import interactive utilities
from OTObit import *
import OTObit

# Set Task Name
OSystem.PSetPgmName("ObitTalk")

# Execute startup script
if os.path.exists(path):
    try:
        execfile (path)
    except Exception, exception:
        print exception
    else:
        pass


# suppress spurious messages about LD_PRELOAD
if os.environ.has_key('LD_PRELOAD'):
    del os.environ['LD_PRELOAD']
    
