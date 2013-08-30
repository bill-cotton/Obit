# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2013
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

#  Python ObitOTFSoln2Cal utilities
import Obit, OTF, OErr, Table

def POTFSoln2Cal(inOTF, outOTF, err):
    """ Apply a Soln table to a Cal table and write a new Cal table

    Calibration parameters are on the inOTF info member.
    If an input Cal table is specified then apply Solutions in this routine,
    if no input Cal table, then copy the Soln table to a new Cal table
    in ObitOTFSolnCopyCal.
    "SOLNUSE"   int scalar Input Solution table version 
    "CALIN"     int scalar Input Cal table version 
                iff <0 then no input Cal table, copy Soln records to output.
    "CALOUT"    int scalar) Output Calibration table version
    returns updated OTFCal table
    inOTF   = Python Obit OTF from which the solution table is appended
    outOTF  = Python Obit OTF on which the calibration tables reside
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
    out.me = Obit.OTFSoln2Cal(inOTF.me, outOTF.me, err.me)
    return out
    # end POTFSoln2Cal
    


