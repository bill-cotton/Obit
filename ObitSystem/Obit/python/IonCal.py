# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007
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

# Python interface to Obit Ionospheric calibration utilities
import Obit, UV, Table, OErr

def PIoN2SolNTableConvert (inUV, outSNVer, NITable, pos, err):
    """ Evaluate Ionospheric model table at pos and convert to SN table

    Returns resultant SN table
    inUV     = UV data for output SN table.
    Control parameters on inUV info member
      "doAntOff" OBIT_bool scalar True if correctionss for antenna
                offset from array center wanted [def False]
    outSNVer = Desired output SN table version, 0=> new
    NITable  = Ionospheric model table to evaluate
    pos      = [RA, Dec] shift (deg) in which NITable to be evaluated.
    err      = Obit Error stack, returns if not empty.
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError, 'PIoN2SolNTableConvert: Bad input UV data'
    if not Table.PIsA(NITable):
        raise TypeError, 'PIoN2SolNTableConvert: Bad NI input table'

    # Create output SN table object
    outSNTable = Table.Table("None")

    # Evaluate
    outSNTable.me = Obit.IoN2SolNTableConvert (inUV.me, outSNVer, \
                                               NITable.me, pos, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error Converting NI to SN table")
    return outSNTable
    # end  PIoN2SolNTableConvert


