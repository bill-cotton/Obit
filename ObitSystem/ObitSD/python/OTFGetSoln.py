""" OTF calibration and editing utilities
"""
# $Id: OTFGetSoln.py,v 1.10 2008/02/27 15:47:10 bcotton Exp $
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

# Python ObitOTFGetSoln utilities
import Obit, OTF, OErr, Table

# Commonly used, dangerous variables
dim=[1,1,1,1,1]

def PCal (inOTF, outOTF, err):
    """ Determine offset calibration for an OTF residual for multibeam systems


    Uses an Atmospheric model across the array of detectors.
    Calibration parameters are on the inOTF info member.
    "solInt"   float scalar Solution interval in days [def 1 sec].
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
    out.me = Obit.OTFGetSolnCal(inOTF.me, outOTF.me, err.me)
    return out
    # end PCal 
    

def PGain (inOTF, outOTF, err):
    """ Determine gain calibration for an OTF from a residual data set.

    Gain Calibration is based on the "Cal" values in the data.
    Average value of the noise Cal is computed and entered in the data. (Gain only)
    Additive values are determined from the median values of the residual data.
    Solution type controlled by calType 
    Calibration parameters are on the inOTF info member.
    "solInt"    float scalar Solution interval in days [def 10 sec].
                This should not exceed 1000 samples.  Solutions will be truncated
                at this limit.
    "minRMS"    float scalar minimum allowable fractional solution RMS [def 0.1].
                bad solutions are replaced with pervious good value. [Gain soln]
    "calJy"     float array Calibrator value in Jy per detector [Gain soln] [def 1.0] .
    "minEl"     float scalar Minimum elevation allowed (deg)
    "calType"   string  Calibration type desired
                "Gain" => Gain (multiplicative, cal) solution only
                "Offset" => Offset (additive) solution only
                "GainOffset" both gain and offset calibration (probably bad idea).
                anything else or absent => Gain only.
    returns OTFSoln table with solutions
    inOTF   = Python Obit OTF (residual) from which the solution is to be determined 
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
    out.me = Obit.OTFGetSolnGain(inOTF.me, outOTF.me, err.me)
    return out
    # end PGain


def PFilter (inOTF, outOTF, err):
    """ Determine offset calibration for an OTF by time filtering a residual data set.

    The time series of each detector is filtered to remove structure on time 
    scales shorter than solInt. 
    Scans in excess of 5000 samples will be broken into several.
    Calibration parameters are on the inOTF info member.
    "solInt"    float scalar Solution interval in days [def 10 sec].
                This should not exceed 1000 samples.  Solutions will be truncated
                at this limit.
    "minEl"     float scalar Minimum elevation allowed (deg)
    returns OTFSoln table with solutions
    inOTF   = Python Obit OTF (residual) from which the solution is to be determined
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
    out.me = Obit.OTFGetSolnFilter(inOTF.me, outOTF.me, err.me)
    return out
    # end PFilter

def PPolyBL (inOTF, outOTF, err):
    """ Fits polynomial additive term on median averages of a residual data set.

    Each solution interval in a scan is median averaged
    (average of 9 points around the median) and then a polynomial fitted.
    Calibration parameters are on the inOTF info member.
    "solInt"    float scalar Solution interval in days [def 10 sec].
                This should not exceed 5000 samples.  Solutions will be truncated
                at this limit.
    "minEl"     float scalar Minimum elevation allowed (deg)
    "Order"     Polynomial order [def 1]
    returns OTFSoln table with solutions
    inOTF   = Python Obit OTF (residual) from which the solution is to be determined
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
    out.me = Obit.OTFGetSolnPolyBL(inOTF.me, outOTF.me, err.me)
    return out
    # end PPolyBL

def PMBBase (inOTF, outOTF, err):
    """ Fits polynomial additive term on median averages of a residual data set.

    Each solution interval in a scan is median averaged
    (average of 9 points around the median) and then a polynomial fitted.
    Calibration parameters are on the inOTF info member.
    "SolInt"    float scalar Solution interval in days [def 10 sec].
                This should not exceed 5000 samples.  Solutions will be truncated
                at this limit.
    "minEl"     float scalar Minimum elevation allowed (deg)
    "Order"     Polynomial order [def 1]
    "clipSig"   data outside of range +/- CLIP sigma are blanked [def large]
    "plotDet"   Detector number to plot per scan [def =-1 = none]
    returns OTFSoln table with solutions
    inOTF   = Python Obit OTF (residual) from which the solution is to be determined
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
    out.me = Obit.OTFGetSolnMBBase(inOTF.me, outOTF.me, err.me)
    return out
    # end PMBBase

def POTFGetInstCal (inOTF, outOTF, err):
    """ Determine instrumental calibration for an OTF multibeam OTF dataset.

    Calculates average gain from cal measurements
    calculates median ofset
    Calibration parameters are on the inOTF info member.
    "solInt"   float scalar Solution interval in days [def 1 sec].
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
    out.me = Obit.OTFGetInstCal(inOTF.me, outOTF.me, err.me)
    return out
    # end POTFGetInstCal


def POTFGetSolnPARGain (inOTF, outOTF, err):
    """ Determine Penn Array type gain calibration from data

    Determine instrumental gain from Penn Array-like cal measurements
    Penn Array-like = slow switching, many samples between state changes
    Average On and Off for each detector and difference.
    Mult factor = calJy/(cal_on-cal_off) per detector, may be negative.
    Data scan averaged and repeated for any subsequent scans without cal On data
    Write Soln entries at beginning and end of each scan.
    Note: Any scans prior to a scan with calibration data will be flagged.
    Calibration parameters are on the inOTF info member.
    "calJy"     OBIT_float (*,1,1) Calibrator value in Jy per detector [def 1.0] .
                Duplicates if only one given.
    returns OTFSoln table with solutions
    inOTF   = Python Obit OTF from which the solution is to be determined
              prior calibration/selection applied if requested
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
    out.me = Obit.OTFGetSolnPARGain(inOTF.me, outOTF.me, err.me)
    return out
    # end POTFGetSolnPARGain

def POTFGetSolnPointTab (inOTF, outOTF, err):
    """ Write pointing corrections into a OTFSoln table

    Given a set of pointing correction, write to OTFSoln table
    Calibration parameters are on the inOTF info member.
        "POffset"   OBIT_float (3,?,1) Table of pointing offsets
                    Triplets (time(day), Azoff(asec), Decoff(asec))
                    MUST be in ascending time order
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
    out.me = Obit.OTFGetSolnPointTab(inOTF.me, outOTF.me, err.me)
    return out
    # end POTFGetSolnPointTab


def POTFGetDummyCal (inOTF, outOTF, inter, ver, ncoef, err):
    """ Create dummy OTFCal table table (applying will not modify data)

    returns OTFCal table with solutions
    inOTF   = input Python Obit OTF 
    outOTF  = Python Obit OTF onto which the cal table is to be appended.
    inter   = time interval (sec) between entries
    ver     = OTFCal table version
    ncoef   = Number of coefficients (across array feeds) in table
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
    out = Table.Table(" ")
    if err.isErr: # existing error?
        return out
    #
    # Set interval on table
    dim[0] = 1; dim[1] = 1
    inInfo = Obit.OTFGetList(inOTF.me)    # return Obit object
    Obit.InfoListPutFloat(inInfo, "solInt",  dim, [inter/86400.0], err.me);
    out.me = Obit.OTFGetDummyCal(inOTF.me, outOTF.me, ver, ncoef, err.me)
    # Cleanup Obit objects
    inInfo    = Obit.InfoListUnref(inInfo)
    #
    return out
    # end OTFGetDummyCal

def PFlag (inOTF, model, outOTF, FGVer, err):
    """ Flag data on basis of comparison of statistics of model, residuals

    Determine flagging from the statistics of time stream data.
    If a model is supplied, this is determined from the residuals to the model.
    The residual data will be derived from inOTF minus the model.
    Detector/intervals in excess of the larger of maxRMS or maxRatio times the 
    equivalent model RMS are flagged in OTFFlag table FGVer on outOTF.
    An evaluation is made independently in each flagInt and detector.
    Control parameters are on the inOTF info member.
    "flagInt"   float Flaffing interval in days [def 10 sec].
                This should not exceed 1000 samples.  Intervals will be truncated
                at this limit.
    "maxRMS"    float  Maximum allowable  RMS in Jy[ def 1.0].
    "maxRatio"  float  Max. allowable ratio to equivalent model RMS [2]
    "minEl"     float  Minimum elevation allowed (deg)
    inOTF   =  Python Obit OTF from which the solution is to be determined 
    model   =  Input CLEAN model, if not provided (None) then the model RMS
               is assumed 0.0
               The pixels values in the image array should be the estimated
               antenna response to the source in that direction.
    outOTF  =  Python Obit OTF onto which the solution table is to be appended.
    FGVer   =  Flag version for output, if 0 create new
    err     =  Python Obit Error/message stack
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
    if err.isErr: # existing error?
        return
    if model==None:
        Obit.OTFGetSolnFlagNoModel(inOTF.me,outOTF.me, FGVer, err.me)
    else:
        Obit.OTFGetSolnFlag(inOTF.me, model.me, outOTF.me, FGVer, err.me)
    # end PFlag



