""" This class is for performing CLEAN on images.

This implements an OTF record based CLEAN
In a record-based clean, the accumulated model is
occasionally subtracted from the OTF data and a new
residual image formed.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2006,2007
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

# Python shadow class to ObitDConCleanOTFRec class
import Obit, OErr, Image, FArray, Table, InfoList, OWindow

class CleanOTFRecPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    #def __del__(self):
    #    if self.thisown == 1 :
    #        # If Obit has been unloaded don't bother
    #        if Obit.__class__ == Obit:
    #            Obit.delete_CleanOTFRec(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.CleanOTFRec_me_set(self.this,value)
            return
        if name=="Dirty":
            PSetDirty(self, value)
            return 
        if name=="Beam":
            PSetBeam(self, value)
            return 
        if name=="Clean":
            PSetClean(self, value)
            return 
        if name=="Weight":
            PSetWeight(self, value)
            return 
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.CleanOTFRec_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="Dirty":
            return PGetDirty(self)
        if name=="Beam":
            return PGetBeam(self)
        if name=="Clean":
            return PGetClean(self)
        if name=="Weight":
            return PGetWeight(self)
        if name=="Size":
            return PGetCleanSize(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C CleanOTFRec instance>"
class CleanOTFRec(CleanOTFRecPtr):
    """ This class is for performing CLEAN on images.

    This implements an OTF record based CLEAN
    In a record-based clean, the accumulated model is
    occasionally subtracted from the OTF data and a new
    residual image formed.

    Arguments to the constructor:
    name   - Name of the CLEAN object (a label)
    inOTF  - Python Obit OTF to be imaged
    err    - Python Obit Error/message stack
    """
    def __init__(self, name, inOTF, err) :
        self.this = Obit.new_CleanOTFRec(name, inOTF, err)
        self.thisown = 1
    def __del__(self):
        if Obit!=None:
            Obit.delete_CleanOTFRec(self.this)

    def DefWindow(self, err):
        """ Set default window (all image)
        
        self   = Python OTF object
        err       = Python Obit Error/message stack
        """
        PDefWindow(self, err)
        # end DefWindow

    def AddWindow(self, window, err):
        """ Add a  window
        
        self   = Python OTF object
        window = set of 4 integers:
                 if window[0]<0 box is round and
                 window[1]=radius, [2,3] = center
                 else rectangular and
                 blc=(window[0],window[1]), trc= blc=(window[2],window[3])
        err    = Python Obit Error/message stack
        """
        PAddWindow(self, window, err)
        # end AddWindow

def PCreate (name, inOTF, err):
    """ Create CleanOTFRec  Object

    returns CleanOTFRec object
    name    = Name for clean
    inOTF   = Python Obit OTF to image
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks 
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return None
    #
    out = CleanOTFRec(name, inOTF.me, err.me)
    
    return out
    # end PCreate


def PGetWindow (inCleanOTFRec):
    """ Return the member OWindow

    returns OWindow
    inCleanOTFRec  = Python CleanOTFRec object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanOTFRec):
        raise TypeError,"inCleanOTFRec MUST be a Python Obit CleanOTFRec"
    #
    out    = OWindow.OWindow()
    out.me = Obit.CleanOTFRecGetWindow(inCleanOTFRec.me)
    return out
    # end PGetWindow

def PSetWindow (inCleanOTFRec, window):
    """ Replace OWindow in the CleanOTFRec

    inCleanOTFRec  = Python CleanOTFRec object
    window      = Python OWindow to attach
    """
    ################################################################
    # Checks
    if not PIsA(inCleanOTFRec):
        raise TypeError,"inCleanOTFRec MUST be a Python ObitCleanOTFRec"
    if not OWindow.PIsA(window):
        raise TypeError,"array MUST be a Python Obit OWindow"
    #
    Obit.CleanOTFRecSetWindow(inCleanOTFRec.me, window.me)
    # end PSetWindow

def PDefWindow (clean, err):
    """ Set default windows on image mosaic member.

    If mosaic member Radius>0 then make round boxes on Fly's eye field
    with this radius, else use rectangular box including all but outer 5 pixels
    On outlier fields, use rectangular box of width OutlierSize.
    Assumes all images in mosaic have descriptors defined.
    clean     = Clean object containing mosaic
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(clean):
        raise TypeError,"mosaic MUST be a Python Obit CleanOTFRec"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.CleanOTFRecDefWindow(clean.me,  err.me)
    # end PDefWindow

def PAddWindow (inCleanOTFRec, window, err):
    """ Add a window to be CLEANed

    inCleanOTFRec = Python CleanOTFRec object
    window     = set of 4 integers:
                 if window[0]<0 box is round and
                 window[1]=radius, [2,3] = center
                 else rectangular and
                 blc=(window[0],window[1]), trc= blc=(window[2],window[3])
    err        = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCleanOTFRec):
        raise TypeError,"inCleanOTFRec MUST be a Python ObitCleanOTFRec"
    if err.isErr: # existing error?
        return 
    #
    Obit.CleanOTFRecAddWindow(inCleanOTFRec.me, window, err.me)
    # end PAddWindow

# Perform Clean
CleanInput={'structure':['Clean',[('CleanOTFRec','CleanOTFRec Object'),
                                  ('outName','Base of output file names'),
                                  ('outDisk','output disk number'),
                                  ('Targets','List of target names to include'),
                                  ('Feeds','List of selected feed numbers, '),
                                  ('Scans','Range of scan numbers,'),
                                  ('BChan','First spectral channel selected. [def all]'),
                                  ('EChan','Highest spectral channel selected. [def all]'),
                                  ('keepCal','True = keep cal-on data [def True]'),
                                  ('Niter','Maximum number of CLEAN iterations'),
                                  ('Patch','Beam patch in pixels [def 100]'),
                                  ('BeamSize','Restoring beam FWHM (deg)'),
                                  ('Gain','CLEAN loop gain'),
                                  ('minFlux','Minimun flux density (Jy)'),
                                  ('doRestore','If True, restore components to image'),
                                  ('noResid','If True do not include residuals in restored image'),
                                  ('fracPeak','Fraction of residual to CLEAN to per major cycle'),
                                  ('Factor','CLEAN depth factor'),
                                  ('autoWindow','Automatically set WIndows?'),
                                  ('doCalib',  '>0 -> calibrate, 2=> also calibrate Weights'),
                                  ('gainUse',  'SN/CL table version number, 0-> use highest'),
                                  ('flagVer',  'Flag table version, 0-> use highest, <0-> none'),
                                  ('timeRange','Selected timerange in days.  0s -> all'),
                                  ('RA',      'Center RA (deg)'),
                                  ('Dec',     'Center Dec (deg)'),
                                  ('xCells',  'Cell spacing in X (asec)'),
                                  ('yCells',  'Cell spacing in Y (asec)'),
                                  ('minWt',   'Minimum sum of weights'),
                                  ('nx',      'Minimum number of cells in X '),
                                  ('ny',      'Minimum number of cells in Y '),
                                  ('Proj',    'Projection string "-SIN", "-ARC", "-TAN"'),
                                  ('ConvType','Convolving function type'),
                                  ('ConvParm','Convolving function parameters'),
                                  ('dispURL', 'URL of display server'),
                                  ('CCVer',   'CC table version number [0 => highest]')]],
            # defaults
            'CleanOTFRec':None,
            'outName':'noName',
            'outDisk':0,
            'Targets':["Any"],
            'Feeds': [0],
            'Scans': [0,1000000],
            'BChan': 1,
            'EChan': 0,
            'keepCal': True,
            'Niter':100,
            'Patch':100,
            'BeamSize':0.0,
            'Gain':0.1,
            'minFlux':0.0,
            'doRestore':True,
            'noResid':False,
            'fracPeak':0.75,
            'Factor':0.0,
            'autoWindow':False,
            'doCalib':0,
            'gainUse':0,
            'flagVer':-1,
            'RA':0.0,
            'Dec':0.0,
            'xCells': 0.0,
            'yCells': 0.0,
            'minWt': 0.1,
            'nx':10,
            'ny':10,
            'Proj':'-SIN',
            'ConvType':3,
            'ConvParm':[0.0,0.0,0.0,0.0],
            'dispURL':'None',
            'CCVer':0}

def PClean (err, input=CleanInput):
    """ Performs OTF record based CLEAN

    The peak in the image is iteratively found and then the beam
    times a fraction of the peak is subtracted and the process is iterated.
    One frequency plane is processed at a time.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
    Input dictionary entries:
    CleanOTFRec = Input CleanOTFRec,
    outName     = Base name for output files, will be prefixed with tmp and tmpBeam
    outDisk     = FITS disk number for output images.
                  and post fixed with .fits for CLEAN image and dirty beam.
    Targets     = List of target names to include
    Feeds       = List of selected feed numbers
    Scans       = Range of scan numbers
    BChan       = First spectral channel selected. [def all]
    EChan       = Highest spectral channel selected. [def all]
    keepCal     = keep cal-on data [def True]
    Niter       = Maximum number of CLEAN iterations
    Patch       = Beam patch in pixels [def 100]
    maxPixel    = Maximum number of residuals [def 20000]
    BeamSize    = Restoring beam (deg)
    Gain        = CLEAN loop gain
    minFlux     = Minimun flux density (Jy)
    doRestore   = If True, restore Clean components to image [def. True]
    noResid     = If True do not include residuals in restored image
    Factor      = CLEAN depth factor
    fracPeak    = Fraction of residual to CLEAN to per major cycle
    autoWindow  = True if autoWindow feature wanted.
    doCalib     = >0 -> calibrate, 2=> also calibrate Weights
    gainUse     = SN/CL table version number, 0-> use highest
    flagVer     = Flag table version, 0-> use highest, <0-> none
    timeRange   = Selected timerange in days.  0s -> all
    RA          = Center RA (deg)
    Dec         = Center Dec (deg)
    xCells      = Cell spacing in X (asec)
    yCells      = Cell spacing in Y (asec)
    minWt       = Minimum sum of weights
    nx          = Minimum number of cells in X 
    ny          = Minimum number of cells in Y 
    Proj        = Projection string "-SIN", "-ARC", "-TAN"
    ConvType    = Convolving function typee: [def=3]
                  0 = pillbox, 3 = Gaussian, 4 = Exp*Sinc,
                  5 = Spherodial wave
    ConvParm    = Convolving function parameters
    dispURL     = URL of display server
    CCVer       = CC table version number

     Gridding convolution functions:
         0 = pillbox, 
         2 = Sinc, 
            Parm[0] = halfwidth in cells,
            Parm[1] = Expansion factor
         3 = Gaussian,
             Parm[0] = halfwidth in cells,[def 3.0]
             Parm[1] = Gaussian with as fraction or raw beam [def 1.0]
         4 = Exp*Sinc
             Parm[0] = halfwidth in cells, [def 2.0]
             Parm[1] = 1/sinc factor (cells) [def 1.55]
             Parm[2] = 1/exp factor (cells) [def 2.52]
             Parm[3] = exp power [def 2.0]
         5 = Spherodial wave
             Parm[0] = halfwidth in cells [def 3.0]
             Parm[1] = Alpha [def 5.0]
             Parm[2] = Expansion factor [not used]

    """
    ################################################################
    # Get input parameters
    inCleanOTFRec  = input["CleanOTFRec"]
    # Checks
    if not PIsA(inCleanOTFRec):
        print "Really is",inCleanOTFRec.__class__
        raise TypeError,"inCleanOTFRec MUST be a Python Obit CleanOTFRec"
    if err.isErr: # existing error?
        return 
    #
    dim = [1,1,1,1,1]
    #
    # Set CLEAN control values on CleanOTFRec
    dim[0] = 1;
    inInfo = PGetList(inCleanOTFRec)    #
    InfoList.PAlwaysPutInt   (inInfo, "Niter",    dim, [input["Niter"]])
    InfoList.PAlwaysPutInt   (inInfo, "Patch",    dim, [input["Patch"]])
    InfoList.PAlwaysPutInt   (inInfo, "CCVer",    dim, [input["CCVer"]])
    InfoList.PAlwaysPutFloat (inInfo, "BeamSize", dim, [input["BeamSize"]])
    InfoList.PAlwaysPutFloat (inInfo, "Gain",     dim, [input["Gain"]])
    InfoList.PAlwaysPutFloat (inInfo, "minFlux",  dim, [input["minFlux"]])
    InfoList.PAlwaysPutFloat (inInfo, "Factor",   dim, [input["Factor"]])
    InfoList.PAlwaysPutBoolean (inInfo, "doRestore",  dim, [input["doRestore"]])
    InfoList.PAlwaysPutBoolean (inInfo, "noResid",    dim, [input["noResid"]])
    InfoList.PAlwaysPutBoolean (inInfo, "autoWindow", dim, [input["autoWindow"]])
    InfoList.PAlwaysPutFloat   (inInfo, "fracPeak",   dim, [input["fracPeak"]])
    # Set imaging controls on OTF
    OTFInfo = PGetList(PGetOTF(inCleanOTFRec))    # 
    InfoList.PAlwaysPutInt   (OTFInfo, "outDisk",  dim, [input["outDisk"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "BChan",    dim, [input["BChan"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "EChan",    dim, [input["ECHan"]])
    InfoList.PAlwaysPutBoolean (OTFInfo, "keepCal",   dim, [input["keepCal"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "doCalib",  dim, [input["doCalib"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "gainUse",  dim, [input["gainUse"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "flagVer",  dim, [input["flagVer"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "nx",  dim, [input["nx"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "ny",  dim, [input["ny"]])
    InfoList.PAlwaysPutFloat (OTFInfo, "RA",  dim, [input["RA"]])
    InfoList.PAlwaysPutFloat (OTFInfo, "Dec",      dim, [input["Dec"]])
    InfoList.PAlwaysPutInt   (OTFInfo, "ConvType", dim, [input["ConvType"]])
    InfoList.PAlwaysPutFloat (OTFInfo, "minWt",    dim, [input["minWt"]])
    xCells = input["xCells"] / 3600.0;
    if xCells>0.0:
        xCells = -xCells
    InfoList.PAlwaysPutFloat (OTFInfo, "xCells",   dim, [input[xCells]])
    yCells = input["yCells"] / 3600.0;
    InfoList.PAlwaysPutFloat (OTFInfo, "yCells",   dim, [input[yCells]])
    dim[0] = len(input["FEEDS"])
    InfoList.PAlwaysPutInt   (OTFInfo, "Feeds",    dim, input["Feeds"])
    dim[0] = len(input["SCANS"])
    InfoList.PAlwaysPutInt   (OTFInfo, "Scans",    dim, input["Scans"])
    dim[0] = len(input["ConvParm"])
    InfoList.PAlwaysPutFloat (OTFInfo, "ConvParm", dim, input["ConvParm"])
    dim[0] = len(input["timeRange"])
    InfoList.PAlwaysPutFloat (OTFInfo, "timeRange",dim, input["timeRange"])
    dim[0] = len(input["outName"])
    InfoList.PAlwaysPutInt   (OTFInfo, "outName",  dim, [input["outName"]])
    dim[0] = len(input["Proj"])
    InfoList.PAlwaysPutInt   (OTFInfo, "Proj",     dim, [input["Proj"]])
    dim[0] = len(input["dispURL"])
    InfoList.PAlwaysPutFloat (inInfo, "dispURL",   dim, input["dispURL"])
    dim[0] = 16; dim[1] = len(input["TARGETS"])
    InfoList.PAlwaysPutString  (OTFInfo, "Targets",dim, input["Targets"])
    #
    # show any errors 
    OErr.printErrMsg(err, "Clean: Error setting parameters")
    #
    # Do operation
    Obit.CleanOTFRecClean(inCleanOTFRec.me, err.me)
    # end PClean

def PRestore (inCln, err):
    """ Restores components

    This is done automatically unless doRestore=False in PClean
    inCln   = Python Obit input OTFClean
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.CleanOTFRecRestore (inCln.me, err.me)
    # end PRestore


def PGetBeam (inCln):
    """ Get Dirty beam image

    returns Dirty beam as Python Obit Image
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = Image.Image("None")
    out.me = Obit.CleanOTFRecGetBeam (inCln.me)
    return out
    # end PGetBeam


def PSetBeam (inCln, image):
    """ Set Dirty beam image

    inCln   = Python Obit input OTFClean
    image   = Python Obit Image for dirty beam
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not Image.ImagePIsA(image):
        raise TypeError,"Image MUST be a Python Obit Image"
    #
    Obit.CleanOTFRecSetBeam (inCln.me, image.me)
    # end PSetBeam


def PGetClean (inCln):
    """ Get Clean image

    returns Clean image as Python Obit Image
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = Image.Image("None")
    out.me = Obit.CleanOTFRecGetClean (inCln.me)
    return out
    # end PGetClean


def PSetClean (inCln, image):
    """ Set Clean image

    inCln   = Python Obit input OTFClean
    image   = Python Obit Image for clean image
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not Image.PIsA(image):
        raise TypeError,"Image MUST be a Python Obit Image"
    #
    Obit.CleanOTFRecSetClean (inCln.me, image.me)
    # end PSetClean


def PGetWeight (inCln):
    """ Get Weight image

    returns Weight image as Python Obit Image
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = Image.Image("None")
    out.me = Obit.CleanOTFRecGetWeight (inCln.me)
    return out
    # end PGetWeight


def PSetWeight (inCln, image):
    """ Set Weight image

    inCln   = Python Obit input OTFClean
    image   = Python Obit Image for clean image
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not Image.PIsA(image):
        raise TypeError,"Image MUST be a Python Obit Image"
    #
    Obit.CleanOTFRecSetWeight (inCln.me, image.me)
    # end PSetWeight


def PGetMosaic (inCln, err):
    """ Return the member ImageMosaic

    returns ImageMosaic
    inCln     = Python CleanOTFRec object
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python CleanOTFRec"
    #
    out    = ImageMosaic.ImageMosaic("None", 1)
    out.me = Obit.UVImagerGetMosaic(inCln.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error obtaining Mosaic member")
    return out
    # end PGetMosaic

def PGetOTF (inCln):
    """ Return the member uvdata

    returns OTF data (as OTF)
    inCln  = Python ImageMosaic object
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit leanOTFRec"
    #
    out    = OTF.OTF("None")
    out.me = Obit.CleanOTFRecGetOTF(inCln.me)
    return out
    # end PGetOTF


def PGetList (inCln):
    """ Get InfoList

    return InfoList
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.CleanOTFRecGetList (inCln.me)
    return out
    # end PGetList


def PGetNiter (inCln):
    """ Get maximum number of CLEAN iterations

    This is only set after the CLEAN has run.
    returns maximum number of CLEAN iterations (int)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFRecGetNiter (inCln.me)
    # end PGetNiter

def PGetGain (inCln):
    """ Get CLEAN loop gain
 
    This is only set after the CLEAN has run.
    returns CLEAN loop gain (float)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFRecGetGain (inCln.me)
    # end PGetGain

def PGetFlux (inCln):
    """  Get min abs flux density for CLEAN

    This is only set after the CLEAN has run.
    returns min abs flux density for CLEAN (float)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFRecGetFlux (inCln.me)
    # end PGetFlux 

def PGetCleanSize (inCln):
    """ Get CLEAN restoring beam size

    This is only set after the CLEAN has run.
    returns CLEAN restoring beam size in pixels (float)
    inCln  = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFRecGetCleanSize (inCln.me)
    # end PGetCleanSize


def PIsA (inCln):
    """ Tells if object thinks it's a Python ObitOTFClean

    return true, false (1,0)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if inCln.__class__ != CleanOTFRec:
        return 0
    return Obit.CleanOTFRecIsA(inCln.me)
    # end PIsA

