# $Id$
""" Utility to estimate and subtract quasi-constant RFI from UV data
"""
#-----------------------------------------------------------------------
#  Copyright (C) 2009
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

# Interferometric RFI excision class
import Obit, OErr, InfoList, UV, Table, types

# Python shadow class to ObitUVRFIXize class
 
class UVRFIXizePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.UVRFIXizeUnref(Obit.UVRFIXize_me_get(self.this))
            # In with the new
            Obit.UVRFIXize_me_set(self.this,value)
            return
        if name=="SkyModel":
            PSetSkyModel(self, value)
            return 
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != UVRFIXize:
            return
        if name == "me" : 
            return Obit.UVRFIXize_me_get(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        if name=="RFI":
            return PGetRFI(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != UVRFIXize:
            return
        return "<C UVRFIXize instance> " + Obit.UVRFIXizeGetName(self.me)
#
class UVRFIXize(UVRFIXizePtr):
    """ Python Obit Image class

    This class does RFI estimation and removal

    UVRFIXize Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List
                (readonly)
    skyModel    - SkyModel, use PGetSkyModel
    """
    def __init__(self, name) :
        self.this = Obit.new_UVRFIXize(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_UVRFIXize(self.this)


def newObit(name, err):
    """ Create and initialize an UVRFIXize structure

    Create sky model object
    Returns the Python UVRFIXize object
    name     = name desired for object (labeling purposes)
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = UVRFIXize (name)
    return out      # seems OK
    # end newObit

def PGetList (inUVRFIXize):
    """ Return the member InfoList

    returns InfoList
    inUVRFIXize  = Python UVRFIXize object
    """
    ################################################################
     # Checks
    if not PIsA(inUVRFIXize):
        raise TypeError,"inUVRFIXize MUST be a Python Obit UVRFIXize"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.UVRFIXizeGetList(inUVRFIXize.me)
    return out
    # end PGetList

# Define UVRFIXize - most parameters can be defined here
RFIXizeInput={'structure':['RFIXize',[('solInt','Counter rotated SN table interval [def 60 sec]'),
                                      ('doInvert','If True invert solution [def False]'),
                                      ('timeInt','Data integration time in sec [def 10 sec]'),
                                      ('timeAvg','Time interval (min) over which to average residual'),
                                      ('minRot','Min. fringe rotation in an integration, less flagged'),
                                      ('maxRot','Max. fringe rotation in an integration, more no RFI'),
                                      ('minRFI','Minimum RFI amplitude (Jy) to subtract')]],
              # defaults
              'solInt':60.0,
              'doInvert':False,
              'timeInt':10.0,
              'timeAvg':0.95,
              'minRot':0.25,
              'maxRot':2.00,
              'minRFI':50}

def PCreate (name, inUV, residUV, outUV, input=RFIXizeInput):
    """ Create the parameters and underlying structures of a UVRFIXize.

    Returns UVRFIXize created.
    name      = Name to be given to object
                Most control parameters are in InfoList member
    inUV      = Input Python uv data to be corrected
    residUV   = Input Python residual uv data
    outUV     = Output Python uv data to be written
    err       = Python Obit Error/message stack
    input     = control parameters:
    solInt      = Counter rotated SN table interval [def 1 min]
    doInvert    = If True invert solution [def False]
    timeInt     = Data integration time in sec [def 10 sec]
    timeAvg     = Time interval (min) over which to average residual.
    minRot      = Min. fringe rotation (turns) in an integration,
                  data with smaller values flagged if > minRFI
                  Only reduce this if a very good model has been
                  subtracted to make the residual data.
    maxRot      = Max. fringe rotation (turns) in an integration 
                  to estimate RFI 
    minRFI      = Minimum RFI amplitude to subtract
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not UV.PIsA(residUV):
        raise TypeError,"residUV MUST be a Python Obit UV"
    if not UV.PIsA(outUV):
        raise TypeError,"outUV MUST be a Python Obit UV"
    #
    # Create
    out = UVRFIXize(name);
    out.me = Obit.UVRFIXizeCreate(name, inUV.me, residUV.me, outUV.me)
    # Set RFIXize control values on out
    inInfo = PGetList(out)    # 
    dim = [1,1,1,1,1]
    # Set control values on RFIXize
    InfoList.PAlwaysPutFloat   (inInfo, "solInt",  dim, [input["solInt"]])
    InfoList.PAlwaysPutBoolean (inInfo, "doInvert",dim, [input["doInvert"]])
    InfoList.PAlwaysPutFloat   (inInfo, "timeInt", dim, [input["timeInt"]])
    InfoList.PAlwaysPutFloat   (inInfo, "timeAvg", dim, [input["timeAvg"]])
    InfoList.PAlwaysPutFloat   (inInfo, "minRFI",  dim, [input["minRFI"]])
    InfoList.PAlwaysPutFloat   (inInfo, "minRot",  dim, [input["minRot"]])
    InfoList.PAlwaysPutFloat   (inInfo, "maxRot",  dim, [input["maxRot"]])
    #
    return out;
    # end PCreate

def PRemove (inRFI, err):
    """" Estimate and remove RFI from UV data

    This technique will estimate and remove a quasi-continuous interfering
    signal from a uvdata set.
    Estimates the contributions to UV data from RFI and removes them
    RFI is estimated by counterrotating a set of residual UV data by the
    interferometer field center fringe rate and the averaging.
    This operation is performed on each channel/IF and polarization datastream.
    The outUV data passed to PCreate will be filled with a corrected
    copy of inUV.

    This operation involves three steps:
        1) Counter-rotate/average residual 
        The residual data is counter-rotated by the field center fringe rate and then
        averaged to timeAvg min and written to scratch file RFIUV.
        2) Filter to estimate RFI \\
        The counter-rotated and averaged residuals are filtered to produce a file which 
        contains the estimated of the RFI contribution to the signal.
        3) Correct input data by estimated RFI \\
        The closest RFI visibility estimate to each visibility in myUV is subtracted to
        write outUV.
    inRFI = Input UVRFIXize object
    err   = Obit error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inRFI):
        raise TypeError,"inRFI MUST be a Python Obit UVRFIXize"
    #
    # Call constituent pieces
    PCounterRot (inRFI, err)
    PFilter (inRFI, err)
    PCorrect (inRFI, err)
    # end PRemove

def PCounterRot (inRFI, err):
    """" Counter rotate/average residual data

    See PRemove for details
    inRFI = Input UVRFIXize object
    err   = Obit error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inRFI):
        raise TypeError,"inRFI MUST be a Python Obit UVRFIXize"
    #
    Obit.UVRFIXizeCounterRot (inRFI.me, err.me)
    # end PCounterRot

def PFilter (inRFI, err):
    """" Filter residual dataset

    See PRemove for details
    inRFI = Input UVRFIXize object
    err   = Obit error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inRFI):
        raise TypeError,"inRFI MUST be a Python Obit UVRFIXize"
    #
    Obit.UVRFIXizeFilter (inRFI.me, err.me)
    # end PFilter

def PCorrect (inRFI, err):
    """" Correct inUV using estimated RFI and write outUV

    See PRemove for details
    inRFI = Input UVRFIXize object
    err   = Obit error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inRFI):
        raise TypeError,"inRFI MUST be a Python Obit UVRFIXize"
    #
    Obit.UVRFIXizeCorrect (inRFI.me, err.me)
    # end PCorrect 

def PGetRFI (inRFI):
    """ Return the data set containing estimated RFI

    returns Residual dataset, this is only valid after calling
    PRemove or PCounterRot.
    If PFilter (or PRemove) has beed called the residual will have been
    filtered to give the estimate of the RFI.
    inRFI  = Python UVRFIXize object
    """
    ################################################################
     # Checks
    if not PIsA(inRFI):
        raise TypeError,"inRFI MUST be a Python Obit UVRFIXize"
    #
    RFIUV    = UV("RFI")
    RFIUV.me = Obit.UVRFIXizeGetRFI (inRFI.me);
    return RFIUV
    # end PGetRFI

def PGetName (inRFI):
    """ Tells Image object name (label)

    returns name as character string
    inRFI  = Python UVRFIXize object
    """
    ################################################################
     # Checks
    if not PIsA(inRFI):
        raise TypeError,"inRFI MUST be a Python Obit UVRFIXize"
    #
    return Obit.UVRFIXizeGetName(inRFI.me)
    # end PGetName

def PIsA (inRFI):
    """ Tells if input really a Python Obit UVRFIXize

    return true, false (1,0)
    inRFI   = Python UVRFIXize object
    """
    ################################################################
    # Checks
    if inRFI.__class__ != UVRFIXize:
        return 0
    return Obit.UVRFIXizeIsA(inRFI.me)
    # end PIsA
