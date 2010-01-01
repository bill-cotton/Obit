# Copyright (C) 2007 Associated Universities, Inc
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

This module provides instances to dispatch function calls locally,
without doing any RPC.

"""

# The AIPSTask module should always be available.
import Proxy.AIPSTask
AIPSTask = Proxy.AIPSTask.AIPSTask()

# The same goes for the ObitTask module.
import Proxy.ObitTask
ObitTask = Proxy.ObitTask.ObitTask()

import Proxy.ObitScriptP
ObitScript = Proxy.ObitScriptP.ObitScript()

# Local AIPS data directories - this is imported in AIPS - recursive
import AIPS
    
# The AIPSData module depends on Obit.  Since Obit might not be
# available, leave out the AIPSUVData and AIPSImage instances if we
# fail to load the module.
try:
    import Proxy.AIPSData
except:
    pass
else:
    AIPSImage  = Proxy.AIPSData.AIPSImage(None)
    AIPSUVData = Proxy.AIPSData.AIPSUVData(None)
    AIPSCat    = Proxy.AIPSData.AIPSCat()
# FITS
try:
    import Proxy.FITSData
except:
    pass
else:
    FITSImage  = Proxy.FITSData.FITSImage(None)
    FITSUVData = Proxy.FITSData.FITSUVData(None)
    FITSCat    = Proxy.FITSData.FITSCat()
