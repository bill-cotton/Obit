# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2008
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

# Python ObitOTFGetAtmCor utilities
import Obit, OTF, OErr, Table

def P(inOTF, outOTF, err):
    """ Determine gain offset calibration for an OTF dataset.

    The offset is determined from an atmospheric model.
    The gain calibration is determined from the average noise cal values
    and the atmospheric opacity.
    Results are placed in a newly created OTFSoln table.
    Calibration parameters are on the inOTF info member.
    "solInt"   float (1,1,1) Solution interval in days [def 10 sec].
    "Tau0"     float (1,1,1) Zenith opacity in nepers [def 0].
    "aTemp"    float (*,1,1) Effective atmospheric temperature in data units.
                i.e. the additional offset per airmass due to the atmosphere.
                per detector in units of the cal.
    "minEl"    float (1,1,1) Minimum elevation (deg)
    "tRx"      float (*,1,1) Receiver temperature per detector in units of the cal
    "calJy"    float (*,1,1) Noise cal value in Jy, per detector [def 1.0] .
    "Azoff"    float (*,1,1) Offset in deg to add to (cross El) [def 0] .
    "Eloff"    float (*,1,1) Offset in deg to add to El  [def 0.0] .

    returns OTFSoln table with solutions
    inOTF   = Python Obit OTF from which the solution is to be determined
    outOTF  = Python Obit OTF onto which the solution table is to be appended.
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OTF.PIsA(outOTF):
        raise TypeError,"outOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    out = Table.Table(" ")
    if err.isErr: # existing error?
        return out
    out.me = Obit.ObitOTFGetAtmCor (inOTF.me, outOTF.me, err.me)
    return out
    # end P

def PAtmEm(inOTF, outOTF, err):
    """ Set atmospheric emission corrections for an OTF dataset.

    Set atmospheric emission corrections into the detector offsets for an OTF dataset
    Results are placed in a newly created OTFSoln table.
    Calibration parameters are on the inOTF info member.
    "solInt"   float (1,1,1) Solution interval in days [def 10 sec].
    "AtmEm"    float (*,1,1) Equivalent zenith atmospheric brightness in Jy [0]
    "Tau0"     float (1,1,1) Zenith opacity in nepers [def 0].

    returns OTFSoln table with solutions
    inOTF   = Python Obit OTF from which the solution is to be determined
    outOTF  = Python Obit OTF onto which the solution table is to be appended.
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OTF.PIsA(outOTF):
        raise TypeError,"outOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    out = Table.Table(" ")
    if err.isErr: # existing error?
        return out
    out.me = Obit.OTFGetAtmEm (inOTF.me, outOTF.me, err.me)
    return out
    # end PAtmEm


