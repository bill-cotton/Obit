"""
Init Obit package

This module is intended to init the basic Obit system and
may be imported in scripts run from the command line.
AIPS.AIPS.userno should be set prior to import or call
OSystem.PSetAIPSuser (user) afterwards to set userid
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007,2019,2024
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
from __future__ import absolute_import
global Adisk, Fdisk, FITSdisks, nFITS, AIPSdisks, nAIPS
global userno, popsno
import Obit, FITSDir, AIPSDir, OSystem, OErr
import ODisplay, AIPS
import os
if not OSystem.PIsInit ():
    # Init Obit if not already done
    popsno = 1
    userno =  AIPS.AIPS.userno
    err=OErr.OErr()
    ObitSys=None
    Adisk = 1
    Fdisk = 0
    ObitSys=OSystem.OSystem ("ObitPython", popsno, userno,
                             AIPSDir.nAIPS, AIPSDir.AIPSdisks, \
                             FITSDir.nFITS, FITSDir.FITSdisks,
                             True, False, err)
    OErr.printErrMsg(err, "Error with Obit startup")

# ObitTalk classes if available
try:
    import AIPS
    import AIPSData, FITSData, AIPSTask, ObitTask
    from AIPSTask import AIPSTask
    from ObitTask import ObitTask
    userno =  AIPS.AIPS.userno
except:
    pass
else:
    pass

if len(AIPSDir.AIPSdisks)==0:
    AIPSDir.AIPSdisks.append('./')
if len(FITSDir.FITSdisks)==0:
    FITSDir.FITSdisks.append('./')
        
