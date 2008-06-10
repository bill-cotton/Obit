# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2006
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

# Python interface to ObitSoln2Cal utilities
import Obit, UV, Table, OErr

# SN table smoothing
Soln2CalInput={
    'structure':['PSoln2Cal',[('solnVer',  'Input Solution (SN) table version '),
                              ('calIn',    'Input Cal (CL) table version, 0=high, -1=none'),
                              ('calOut',   'Output Calibration table version, 0=>create new'),
                              ('subA',     'Selected subarray (default 1)'),
                              ('interMode','Interpolation mode 2PT, SELF POLY SIMP AMBG CUBE MWF '),
                              ('interParm','interpolation parameters'),
                              ('interNPoly','number of terms in polynomial'),
                              ('maxInter', 'Max. time (day) over which to interpolate.'),
                              ('allPass', 'If true copy unmodified entries as well.'),
                              ('refAnt',  'Ref ant to use. (default 1)')]],
    # defaults
    'solnVer':0,
    'calIn':0,
    'calOut':0,
    'subA':1,
    'interMode':"    ",
    'interNPoly':2,
    'maxInter':1.0,
    'allPass':False,
    'refAnt':1}
def PSoln2Cal (inUV, outUV, err, input=Soln2CalInput):
    """ Apply a gain solution to a calibration table

    inUV     = UV data with solution and input calibration
    outUV    = UV data for output calibration table
    err      = Python Obit Error/message stack
    input    = input parameter dictionary
    
    Input dictionary entries:
    solnVer = Input Solution (SN) table version 
    calIn   = Input Cal (CL) table version, 0=high, -1=none
    calOut  = Output Calibration table version, 0=>create new
    subA    = Selected subarray (default 1)
    interMode =  Interpolation mode 2PT, SELF POLY SIMP AMBG CUBE MWF '),
        "2PT " = linear vector interpolation with no SN smoothing.
        "SELF" = Use only SN solution from same source which is closest in time.
        "POLY" = Fit a polynomial to the SN rates and delays.
                  Use the integral of the rate polynomial for the phases. (NYI)
        "SIMP" = Simple linear phase connection between SN phase
                 entries, assumes phase difference less than 180 degrees.
        "AMBG" = Linear phase connection using rates to resolve phase ambiguities.
        "CUBE" = As AMBG but fit third order polynomial to phases and rates.
        "MWF " = Median window filter of SN table before 2PT interpolation
        "GAUS" = Gaussian smoothing of SN table before 2PT interpolation,
        "BOX " = Boxcar smoothing of SN table before 2PT interpolation,
    interParm =  interpolation parameters, smoothing time (hr)
                 amplitude, phase, delay/rate
    interNPoly = number of terms in polynomial'),
    allPass = If true copy unmodified entries as well (default False)
    refAnt  = Ref ant to use. (default 1)
    """
    ################################################################
    # Get input parameters
    InData    = input["InData"]
    InTable   = input["InTable"]
    isuba     = input["isuba"]
    #
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError, 'PSoln2Cal: Bad input UV data'
    if not UV.PIsA(outUV):
        raise TypeError, 'PSoln2Cal: Bad output UV data'

    # Set control values on UV 
    dim[0] = 1;
    inInfo = UV.PGetList(InData)  # Add control to UV data
    dim[0] = 4
    InfoList.PAlwaysPutString  (inInfo, "interMode", dim, [input["interMode"]])
    dim[0] = 1;
    InfoList.PAlwaysPutInt   (inInfo, "solnVer",  dim, [input["solnVer"]])
    InfoList.PAlwaysPutInt   (inInfo, "calIn",  dim, [input["calIn"]])
    InfoList.PAlwaysPutInt   (inInfo, "calOut",  dim, [input["calOut"]])
    InfoList.PAlwaysPutInt   (inInfo, "subA",  dim, [input["subA"]])
    InfoList.PAlwaysPutInt   (inInfo, "interNPoly",  dim, [input["interNPoly"]])
    InfoList.PAlwaysPutInt   (inInfo, "refAnt",  dim, [input["refAnt"]])
    InfoList.PAlwaysPutBool  (inInfo, "allPass",dim, [input["allPass"]])
    dim[0] = len(input["interParm"])
    InfoList.PAlwaysPutFloat (inInfo, "allPass",dim, input["interParm"])
    # Calibrate
    Obit.UVSoln2Cal (inUV, outUV, err)
    if err.isErr:
        printErrMsg(err, "Error applying SN table to CL table")
    # end  PSoln2Cal


