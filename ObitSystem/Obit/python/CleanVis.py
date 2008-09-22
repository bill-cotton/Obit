""" Python Obit Visibility-based CLEAN class

CleanVis Members with python interfaces:
List      - used to pass instructions to processing 
Mosaic    - ImageMosaic
SkyModel  - Sky Model
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2005,2007
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

# Obit CleanVis
import Obit, OErr, ImageMosaic, InfoList, UV, OWindow, SkyModel

# Python shadow class to ObitDConCleanVis class
 
class CleanVisPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.CleanVisUnref(Obit.CleanVis_me_get(self.this))
            # In with the new
            Obit.CleanVis_me_set(self.this,value)
            return
        if name=="Mosaic":
            PSetMosaic(self, value)
            return 
        if name=="SkyModel":
            PSetSkyModel(self, value)
            return 
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != CleanVis:
            return
        if name == "me" : 
            return Obit.CleanVis_me_get(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        if name=="Number":
            return PGetNumber(self)
        if name=="Mosaic":
            return PGetMosaic(self)
        if name=="SkyModel":
            return PGetSkyModel(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != CleanVis:
            return
        return "<C CleanVis instance> " + Obit.CleanVisGetName(self.me)
#
class CleanVis(CleanVisPtr):
    """ Python Obit Image class

    This class does visibility-based (Cotton-Schwab) CLEAN

    CleanVis Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List (readonly)
    mosaic    - ImageMosaic, use PGetMosaic
    """
    def __init__(self, name) :
        self.this = Obit.new_CleanVis(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_CleanVis(self.this)

    def DefWindow(self, err):
        """ Set default window (all image)
        
        self   = Python OTF object
        err       = Python Obit Error/message stack
        """
        PDefWindow(self, err)
        # end DefWindow

    def AddWindow(self, field, window, err):
        """ Add a  window
        
        self   = Python OTF object
        field  = Which field (1-rel) is the window in?
        window = set of 4 integers:
                 if window[0]<0 box is round and
                 window[1]=radius, [2,3] = center
                 else rectangular and
                 blc=(window[0],window[1]), trc= blc=(window[2],window[3])
        err    = Python Obit Error/message stack
        """
        PAddWindow(self, field, window, err)
        # end AddWindow

# Commonly used, dangerous variables
dim=[1,1,1,1,1]
blc=[1,1,1,1,1,1,1]
trc=[0,0,0,0,0,0,0]
err=OErr.OErr()

def input(inputDict):
    """ Print the contents of an input Dictionary

    inputDict = Python Dictionary containing the parameters for a routine
    There should be a member of the dictionary ('structure') with a value
    being a list containing:
    1) The name for which the input is intended (string)
    2) a list of tuples consisting of (parameter name, doc string)
       with an entry for each parameter in the dictionary.
       The display of the the inputs dictionary will be in the order of
       the tuples and display the doc string after the value.
       An example:
       Soln2CalInput={'structure':['Soln2Cal',[('InData','Input OTF'),
                                               ('soln','input soln table version'),
                                               ('oldCal','input cal table version, -1=none'),
                                               ('newCal','output cal table')]],
                      'InData':None, 'soln':0, 'oldCal':-1, 'newCal':0}
    """
    ################################################################
    structure = inputDict['structure']  # Structure information
    print 'Inputs for ',structure[0]
    for k,v in structure[1]:
        print '  ',k,' = ',inputDict[k],' : ',v
        
    # end input

def newObit(name, err):
    """ Create and initialize an CleanVis structure

    Create sky model object
    Returns the Python CleanVis object
    name     = name desired for object (labeling purposes)
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = CleanVis (name)
    return out      # seems OK
    # end newObit

def PCopy (inCleanVis, outCleanVis, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inCleanVis, copies pointers
    inCleanVis  = Python CleanVis object to copy
    outCleanVis = Output Python CleanVis object, must be defined
    err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python Obit CleanVis"
    if not PIsA(outCleanVis):
        raise TypeError,"outCleanVis MUST be a Python Obit CleanVis"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.CleanVisCopy (inCleanVis.me, outCleanVis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying CleanVis")
    # end PCopy

def PGetList (inCleanVis):
    """ Return the member InfoList

    returns InfoList
    inCleanVis  = Python CleanVis object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python Obit CleanVis"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.CleanVisGetList(inCleanVis.me)
    return out
    # end PGetList

def PGetMosaic (inCleanVis):
    """ Return the member mosaic

    returns ImageMosaic
    inCleanVis  = Python CleanVis object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python Obit CleanVis"
    #
    out    = ImageMosaic.ImageMosaic("None", 1)
    out.me = Obit.CleanVisGetImageMosaic(inCleanVis.me)
    return out
    # end PGetMosaic

def PSetMosaic (inCleanVis, mosaic):
    """ Replace an ImageMosaic in the CleanVis

    inCleanVis  = Python CleanVis object
    mosaic      = Python ImageMosaic to attach
    """
    ################################################################
    # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python ObitCleanVis"
    if not ImageMosaic.PIsA(mosaic):
        raise TypeError,"array MUST be a Python Obit ImageMosaic"
    #
    Obit.CleanVisSetImageMosaic(inCleanVis.me, mosaic.me)
    # end PSetMosaic

def PGetSkyModel (inCleanVis):
    """ Return the member skymodel

    returns SkyModel
    inCleanVis  = Python CleanVis object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python Obit CleanVis"
    #
    out    = SkyModel.SkyModel("None")
    out.me = Obit.CleanVisGetSkyModel(inCleanVis.me)
    return out
    # end PGetSkyModel

def PSetSkyModel (inCleanVis, skymodel):
    """ Replace an ImageSkyModel in the CleanVis

    inCleanVis  = Python CleanVis object
    skymodel    = Python SkyModel to attach
    """
    ################################################################
    # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python ObitCleanVis"
    if not SkyModel.PIsA(skymodel):
        raise TypeError,"array MUST be a Python Obit ImageSkyModel"
    #
    Obit.CleanVisSetImageSkyModel(inCleanVis.me, skymodel.me)
    # end PSetSkyModel

def PGetWindow (inCleanVis):
    """ Return the member OWindow

    returns OWindow
    inCleanVis  = Python CleanVis object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python Obit CleanVis"
    #
    out    = OWindow.OWindow()
    out.me = Obit.CleanVisGetWindow(inCleanVis.me)
    return out
    # end PGetWindow

def PSetWindow (inCleanVis, window):
    """ Replace OWindow in the CleanVis

    inCleanVis  = Python CleanVis object
    window      = Python OWindow to attach
    """
    ################################################################
    # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python ObitCleanVis"
    if not OWindow.PIsA(window):
        raise TypeError,"array MUST be a Python Obit OWindow"
    #
    Obit.CleanVisSetWindow(inCleanVis.me, window.me)
    # end PSetWindow

def PAddWindow (inCleanVis, field, window, err):
    """ Add a window to a field to be CLEANed

    inCleanVis  = Python CleanVis object
    field         = Which field (1-rel) is the window in?
    window        = set of 4 integers:
                    if window[0]<0 box is round and
                    window[1]=radius, [2,3] = center
                    else rectangular and
                    blc=(window[0],window[1]), trc= blc=(window[2],window[3])
    """
    ################################################################
    # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python ObitCleanVis"
    #
    Obit.CleanVisAddWindow(inCleanVis.me, field, window, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error adding CLEAN window")
    # end PAddWindow

# Define CleanVis - most parameters must be defined here
CleanInput={'structure':['Clean',[('Niter','Maximum number of CLEAN iterations'),
                                  ('minPatch','Minimum beam patch in pixels [def 100]'),
                                  ('maxPixel','Maximum number of residuals [def 20000]'),
                                  ('BMAJ','Restoring beam major axis (deg)'),
                                  ('BMIN','Restoring beam minor axis (deg)'),
                                  ('BPA','Restoring beam position angle (deg)'),
                                  ('Gain','CLEAN loop gain'),
                                  ('minFlux','Minimun flux density (Jy)'),
                                  ('Factor','CLEAN depth factor'),
                                  ('Plane','Plane being processed, 1-rel indices of axes 3-?'),
                                  ('autoWindow','Automatically set WIndows?'),
                                  ('CCVer','CC table version number [0 => highest]'),
                                  ('Mode','Model mode, 0=fastest, 1=DFT, 2=Grid'),
                                  ('doCalSelect','Select/calibrate/edit data?'),
                                  ('Stokes',   'Stokes parameter, blank-> unchanged from input'),
                                  ('BChan',    'First spectral channel selected. [def all]'),
                                  ('EChan',    'Highest spectral channel selected. [def all]'),
                                  ('BIF',      'First IF selected. [def all]'),
                                  ('EIF',      'Highest IF selected. [def all]'),
                                  ('doPol',    '>0 -> calibrate polarization.'),
                                  ('doCalib',  '>0 -> calibrate, 2=> also calibrate Weights'),
                                  ('gainUse',  'SN/CL table version number, 0-> use highest'),
                                  ('flagVer',  'Flag table version, 0-> use highest, <0-> none'),
                                  ('BLVer',    'BL table version, 0> use highest, <0-> none'),
                                  ('BPVer',    'Band pass (BP) table version, 0-> use highest'),
                                  ('Subarray', 'Selected subarray, <=0->all [default all]'),
                                  ('freqID',   'Selected Frequency ID, <=0->all [default all]'),
                                  ('timeRange','Selected timerange in days.  0s -> all'),
                                  ('UVRange',  'Selected UV range in wavelengths. 0s -> all'),
                                  ('Sources',  'Source names selected unless any starts with'),
                                  ('Antennas', 'A list of selected antenna numbers, if any is negative'),
                                  ('corrType', 'Correlation type, 0=cross corr only, 1=both, 2=auto only.'),
                                  ('doBand',   'Band pass application type <0-> none'),
                                  ('Smooth',   'Specifies the type of spectral smoothing'),
                                  ('DoWeight', 'If True apply uniform weighting corrections to uvdata'),
                                  ('PBCor',    'If True make freq. dependent rel. pri. beam corr.'),
                                  ('Robust',   'Briggs robust parameter. (AIPS definition)'),
                                  ('UVTaper',  'UV plane taper, sigma in klambda, deg as [maj, min, pa]'),
                                  ('WtSize',   'Size of weighting grid in cells [image]'),
                                  ('WtBox',    'Size of weighting box in cells [def 0]'),
                                  ('WtFunc',   'Weighting convolution function [def. 1]'),
                                  ('WtPower',  'Power to raise weights to.  [def = 1.0]'),
                                  ('doBeam',  'True if beams are to be made'),
                                  ('Type',    'Underlying file type, 0=FITS, 1=AIPS'),
                                  ('Name',    'Name of image'),
                                  ('Class',   'Root of class name'),
                                  ('Seq',     'Sequence number'),
                                  ('Disk',    'Disk number for underlying files'),
                                  ('FOV',     'Field of view (deg) for Mosaic'),
                                  ('doFull',  'If true, make full field image'),
                                  ('NField',  'Number of fields defined in input'),
                                  ('xCells',  'Cell spacing in X (asec) for all images'),
                                  ('yCells',  'Cell spacing in Y (asec) for all images'),
                                  ('nx',      'Minimum number of cells in X for NField images'),
                                  ('ny',      'Minimum number of cells in Y for NField image'),
                                  ('RAShift', 'Right ascension shift (AIPS convention) for each field'),
                                  ('DecShift','Declination for each field'),
                                  ('Catalog', 'AIPSVZ format catalog for defining outliers, None=do not use'),
                                  ('OutlierFlux', 'Minimum estimated outlyer flux density (Jy)'),
                                  ('OutlierDist', 'Maximum distance to add outlyers (deg)'),
                                  ('OutlierSI',   'Spectral index to estimate flux density'),
                                  ('OutlierSize', 'Size of outlyer field (pixels)'),
                                  ('doRestore', 'Restore image when done? [def True]'),
                                  ('doFlatten', 'Flatten image when done? [def True]')]],
            # defaults
            'Niter':100,
            'minPatch':100,
            'maxPixel':20000,
            'BMAJ':0.0,
            'BMIN':0.0,
            'BPA':0.0,
            'Gain':0.1,
            'minFlux':0.0,
            'Factor':0.0,
            'Plane':[1,1,1,1,1],
            'autoWindow':True,
            'CCVer':0,
            'Mode':0,
            'doCalSelect':False,
            'Stokes':'FULL',
            'BChan':0,
            'EChan':0,
            'BIF':0,
            'EIF':0,
            'doPol':False,
            'doCalib':0,
            'gainUse':0,
            'flagVer':-1,
            'BLVer':-1,
            'BPVer':-1,
            'Subarray':0,
            'freqID':0,
            'timeRange':[0.0,10.0],
            'UVRange':[0.0,0.0],
            'Sources':[""],
            'Antennas': [0],
            'corrType':0,
            'doBand':-1,
            'Smooth':[0.0, 0.0, 0.0],
            'DoWeight': True,                                   
            'PBCor': True,                                   
            'Robust': 0.0,                                   
            'UVTaper': [0.0,0.0,0.0],
            'WtSize':-1,
            'WtBox':0,
            'WtFunc':1,
            'WtPower':1.0,
            'doBeam':False,
            'Type':1,
            'Name':None,
            'Class':None,
            'Seq':0,
            'Disk':1,
            'FOV':0.0,
            'doFull':False,
            'NField':1,
            'xCells':0.0,
            'yCells':0.0,
            'nx':[128],
            'ny':[128],
            'RAShift':[0.0],
            'DecShift':[0.0],
            'Catalog':'None',
            'OutlierFlux':0.1,
            'OutlierDist':1.0,
            'OutlierSI':-0.75,
            'OutlierSize':50,
            'doRestore':True,
            'doFlatten':True }
def PCreate (name, uvdata, err, input=CleanInput):
    """ Create the parameters and underlying structures of a CleanVis.

    Returns CleanVis created.
    name      = Name to be given to object
                Most control parameters are in InfoList member
    uvdata    = Python uv data from which image is to be made
    err       = Python Obit Error/message stack
    input     = control parameters:
    Niter       = Maximum number of CLEAN iterations
    minPatch    = Minimum beam patch in pixels [def 100]
    maxPixel    = Maximum number of residuals [def 20000]
    BMAJ        = Restoring beam major axis (deg)
    BMIN        = Restoring beam minor axis (deg)
    BPA         = Restoring beam position angle (deg)
    Gain        = CLEAN loop gain
    minFlux     = Minimun flux density (Jy)
    Factor      = CLEAN depth factor
    Plane       = Plane being processed, 1-rel indices of axes 3-?
    autoWindow  = True if autoWindow feature wanted.
    CCVer       = CC table version number
    Mode        = Model mode, 0=fastest, 1=DFT, 2=Grid
    doCalSelect = Select/calibrate/edit data?),
    Stokes   = Stokes parameter, blank-> unchanged from input),
    BChan    = First spectral channel selected. [def all]),
    EChan    = Highest spectral channel selected. [def all]),
    BIF      = First IF selected. [def all]),
    EIF      = Highest IF selected. [def all]),
    doPol    = >0 -> calibrate polarization.),
    doCalib  = >0 -> calibrate, 2=> also calibrate Weights),
    gainUse  = SN/CL table version number, 0-> use highest),
    flagVer  = Flag table version, 0-> use highest, <0-> none),
    BLVer    = BL table version, 0> use highest, <0-> none),
    BPVer    = Band pass (BP) table version, 0-> use highest),
    Subarray = Selected subarray, <=0->all [default all]),
    freqID   = Selected Frequency ID, <=0->all [default all]),
    timeRange= Selected timerange in days. [2 floats] 0s -> all),
    UVRange  = Selected UV range in wavelengths. 0s -> all),
    Sources  = Source names selected unless any starts with),
    Antennas = A list of selected antenna numbers, if any is negative),
    corrType = Correlation type, 0=cross corr only, 1=both, 2=auto only.),
    doBand   = Band pass application type <0-> none),
    Smooth   = Specifies the type of spectral smoothing [three floats]
    DoWeight = True if Weighting to be applied
    PBCor    = If True make freq. dependent rel. pri. beam corr.
    Robust   = Briggs robust parameter. (AIPS definition)
    UVTaper  = UV plane taper, sigma in klambda,deg as [maj, min, pa]
    WtSize   = Size of weighting grid in cells [same as image nx]
    WtBox    = Size of weighting box in cells [def 1]
    WtFunc   = Weighting convolution function [def. 1]
               1=Pill box, 2=linear, 3=exponential, 4=Gaussian
               if positive, function is of radius, negative in u and v.
    WtPower  = Power to raise weights to.  [def = 1.0]
               Note: a power of 0.0 sets all the output weights to 1 as modified
               by uniform/Tapering weighting.
               Applied in determinng weights as well as after.
    DoBeam   = True if beams are to be made
    Type     = Underlying file type, 0=FITS, 1=AIPS
    Name     = Name of image, used as AIPS name or to derive FITS filename
    Class    = Root of class, used as AIPS class or to derive FITS filename
    Seq      = Sequence number
    Disk     = Disk number for underlying files
    FOV      = Field of view (deg) for Mosaic
    doFull   = If True, create full field (FOV) image
    NField   = Number of fields defined in input,
               if unspecified derive from data and FOV
    xCells   = Cell spacing in X (asec) for all images,
               if unspecified derive from data
    yCells   = Cell spacing in Y (asec) for all images,
               if unspecified derive from data
    nx       = Minimum number of cells in X for NField images
               if unspecified derive from data
    ny       = Minimum number of cells in Y for NField images
               if unspecified derive from data
    RAShift  = Right ascension shift (AIPS convention) for each field
               if unspecified derive from FOV and data
    DecShift = Declination for each field
               if unspecified derive from FOV and data
    Catalog  = AIPSVZ format catalog for defining outliers, None=do not use
    OutlierFlux = Minimum estimated outlyer flux density (Jy)
    OutlierDist = Maximum distance to add outlyers (deg)
    OutlierSI   = Spectral index to estimate flux density
    OutlierSize = Size of outlyer field (pixels)
    doRestore   = Restore image when done? [def True]
    doFlatten   = Flatten image when done? [def True]
    """
    ################################################################
    # Checks
    if not UV.PIsA(uvdata):
        raise TypeError,"uvData MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    dim = [1,1,1,1,1]
    # Set imaging control values on uvdata
    dim[0] = 1;
    inInfo = UV.PGetList(uvdata)    # 
    InfoList.PPutBoolean (inInfo, "doCalSelect",  dim, [input["doCalSelect"]], err)
    dim[0] = 4
    InfoList.PAlwaysPutString (inInfo, "Stokes",   dim, [input["Stokes"]])
    dim[0] = 1;
    InfoList.PPutInt  (inInfo, "BChan",      dim, [input["BChan"]],    err)
    InfoList.PPutInt  (inInfo, "EChan",      dim, [input["EChan"]],    err)
    InfoList.PPutInt  (inInfo, "BIF",        dim, [input["BIF"]],      err)
    InfoList.PPutInt  (inInfo, "EIF",        dim, [input["EIF"]],      err)
    itemp = int(input["doPol"])
    InfoList.PPutInt  (inInfo, "doPol",      dim, [itemp],             err)
    InfoList.PPutInt  (inInfo, "doCalib",    dim, [input["doCalib"]],  err)
    InfoList.PPutInt  (inInfo, "doBand",     dim, [input["doBand"]],   err)
    InfoList.PPutInt  (inInfo, "gainUse",    dim, [input["gainUse"]],  err)
    InfoList.PPutInt  (inInfo, "flagVer",    dim, [input["flagVer"]],  err)
    InfoList.PPutInt  (inInfo, "BLVer",      dim, [input["BLVer"]],    err)
    InfoList.PPutInt  (inInfo, "BPVer",      dim, [input["BPVer"]],    err)
    InfoList.PPutInt  (inInfo, "Subarray",   dim, [input["Subarray"]], err)
    InfoList.PPutInt  (inInfo, "freqID",     dim, [input["freqID"]],   err)
    InfoList.PPutInt  (inInfo, "corrType",   dim, [input["corrType"]], err)
    dim[0] = 2
    InfoList.PPutFloat (inInfo, "UVRange",   dim, input["UVRange"],    err)
    dim[0] = 3
    InfoList.PPutFloat (inInfo, "Smooth",    dim, input["Smooth"],     err)
    dim[0] = 2
    InfoList.PPutFloat (inInfo, "timeRange", dim, input["timeRange"],  err)
    dim[0] = len(input["Antennas"])
    InfoList.PPutInt  (inInfo, "Antennas",   dim, input["Antennas"],   err)
    dim[0] = 16; dim[1] = len(input["Sources"])
    InfoList.PAlwaysPutString  (inInfo, "Sources", dim, input["Sources"])
    # Weighting parameters
    dim[0] = 1; dim[1] = 1;
    InfoList.PPutBoolean (inInfo,"DoWeight",dim,[input["DoWeight"]],err)
    InfoList.PPutBoolean (inInfo,"PBCor",  dim, [input["PBCor"]],   err)
    InfoList.PPutFloat  (inInfo, "Robust", dim, [input["Robust"]],  err)
    InfoList.PPutInt    (inInfo, "WtBox",  dim, [input["WtBox"]],   err)
    InfoList.PPutInt    (inInfo, "WtFunc", dim, [input["WtFunc"]],  err)
    InfoList.PPutFloat  (inInfo, "WtPower",dim, [input["WtPower"]], err)
    dim[1] = len(input["UVTaper"])
    InfoList.PPutFloat  (inInfo, "Taper",  dim, input["UVTaper"],   err)
    WtSize   = input["WtSize"]
    if (WtSize>0):
        print "WtSize", WtSize
        # Change name for C routine.
        dim[0] = 1;
        InfoList.PPutInt  (inInfo, "nuGrid",  dim, [WtSize], err)
        InfoList.PPutInt  (inInfo, "nvGrid",  dim, [WtSize], err)
    # Define image
    dim[0] = 1; dim[1] = 1;
    InfoList.PPutInt    (inInfo, "imFileType", dim, [input["Type"]],   err)
    InfoList.PPutInt    (inInfo, "imSeq",      dim, [input["Seq"]],    err)
    InfoList.PPutInt    (inInfo, "imDisk",     dim, [input["Disk"]],   err)
    InfoList.PPutFloat  (inInfo, "FOV",      dim, [input["FOV"]],    err)
    InfoList.PPutBoolean (inInfo, "doFull",  dim, [input["doFull"]], err)
    InfoList.PPutInt    (inInfo, "NField",   dim, [input["NField"]], err)
    InfoList.PPutFloat  (inInfo, "xCells",   dim, [input["xCells"]], err)
    InfoList.PPutFloat  (inInfo, "yCells",   dim, [input["yCells"]], err)
    InfoList.PPutFloat  (inInfo, "OutlierFlux", dim, [input["OutlierFlux"]], err)
    InfoList.PPutFloat  (inInfo, "OutlierDist", dim, [input["OutlierDist"]], err)
    InfoList.PPutFloat  (inInfo, "OutlierSI",   dim, [input["OutlierSI"]],   err)
    InfoList.PPutInt    (inInfo, "OutlierSize", dim, [input["OutlierSize"]], err)
    InfoList.PPutFloat  (inInfo, "BMAJ",     dim, [input["BMAJ"]],     err)
    InfoList.PPutFloat  (inInfo, "BMIN",     dim, [input["BMIN"]],     err)
    InfoList.PPutFloat  (inInfo, "BPA",      dim, [input["BPA"]],      err)
    dim[0] = len(input["Name"])
    InfoList.PAlwaysPutString (inInfo, "imName",     dim, [input["Name"]])
    dim[0] = len(input["Class"])
    InfoList.PAlwaysPutString (inInfo, "imClass",    dim, [input["Class"]])
    dim[0] = len(input["Catalog"])
    InfoList.PAlwaysPutString (inInfo, "Catalog",  dim, [input["Catalog"]])
    dim[0] = len(input["nx"])
    InfoList.PAlwaysPutInt    (inInfo, "nx",       dim, input["nx"])
    dim[0] = len(input["ny"])
    InfoList.PAlwaysPutInt    (inInfo, "ny",       dim, input["ny"])
    dim[0] = len(input["RAShift"])
    InfoList.PAlwaysPutFloat  (inInfo, "RAShift",  dim, input["RAShift"])
    dim[0] = len(input["DecShift"])
    InfoList.PAlwaysPutFloat  (inInfo, "DecShift", dim, input["DecShift"])
    #OErr.printErrMsg(err, "CleanVisCreate: Error setting parameters")
    #
    # Create
    out = CleanVis(name);
    out.me = Obit.CleanVisCreate(name, uvdata.me,  err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating CleanVis")
    # Set Clean control values on out
    dim[0] = 1;
    inInfo = PGetList(out)    # 
    InfoList.PPutInt   (inInfo, "Niter",    dim, [input["Niter"]],    err)
    InfoList.PPutInt   (inInfo, "minPatch", dim, [input["minPatch"]], err)
    InfoList.PPutInt   (inInfo, "maxPixel", dim, [input["maxPixel"]], err)
    InfoList.PPutInt   (inInfo, "CCVer",    dim, [input["CCVer"]],    err)
    InfoList.PPutInt   (inInfo, "Mode",     dim, [input["Mode"]],    err)
    InfoList.PPutFloat (inInfo, "BMAJ",     dim, [input["BMAJ"]],     err)
    InfoList.PPutFloat (inInfo, "BMIN",     dim, [input["BMIN"]],     err)
    InfoList.PPutFloat (inInfo, "BPA",      dim, [input["BPA"]],      err)
    InfoList.PPutFloat (inInfo, "Gain",     dim, [input["Gain"]],     err)
    InfoList.PPutFloat (inInfo, "minFlux",  dim, [input["minFlux"]],  err)
    InfoList.PPutFloat (inInfo, "Factor",   dim, [input["Factor"]],   err)
    InfoList.PPutBoolean (inInfo, "doRestore",  dim, [input["doRestore"]], err)
    InfoList.PPutBoolean (inInfo, "doFlatten",  dim, [input["doFlatten"]], err)
    InfoList.PPutBoolean (inInfo, "autoWindow", dim, [input["autoWindow"]],err)
    dim[0] = len(input["Plane"])
    InfoList.PAlwaysPutInt   (inInfo, "Plane",    dim, input["Plane"])
    # show any errors 
    #OErr.printErrMsg(err, "CleanVisCreate: Error setting parameters")
    #
    return out;
    # end PCreate

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
        raise TypeError,"mosaic MUST be a Python Obit CleanVis"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
   #
    Obit.CleanVisDefWindow(clean.me,  err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error with default windows")
    # end PDefWindow

def PClean (inCleanVis, err):
    """ Performs visibiility-based CLEAN

    The peak in the image is iteratively found and then a limited
    region of the beam times a fraction of the peak is subtracted
    and the process is iterated.  Occasionally a proper residual
    image is formed by Fourier transforming the Clean components,
    subtracting them from the visibility data and then imaging the
    residual data.  When done, the componsnts are restored and
    if needed, the image mosaic flattened to a single plane.
    inCleanVis = CleanVis object
    err        = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python Obit CleanVis"
    #
    # Do operation
    Obit.CleanVisDeconvolve(inCleanVis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error in deconvolution")
    # end PClean

def PReimage (inCleanVis, uvdata, err):
    """ See if an image needs to be remade

    See if an image needs to be remade because a source which exceeds
    the flux  threshold is not centered (as determined by moments)
    on the reference pixel (within toler pixel).
    A new (96x96) field is added centered on the offending source and a negative
    clean window added to the position of the source in its original window.
    Avoid duplicates of the same source and ensure that all occurances of this 
    source in any exant field has a negative clean window added.
    Multiple centering sources per facet are allowed
    A boolean entry "autoCenField" with value True is added to the info member of any
    image members added to the mosaic member.
    
    inCleanVis  = Python CleanVis object
    uvdata      = Python uv data from which image is made
    err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python ObitCleanVis"
    if not UV.PIsA(uvdata):
        raise TypeError,"uvData MUST be a Python Obit UV"
    #
    out = Obit.CleanVisReimage(inCleanVis.me, uvdata.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error in Reimage")
    return out
    # end PReimage

def PGetName (inCleanVis):
    """ Tells Image object name (label)

    returns name as character string
    inCleanVis  = Python CleanVis object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanVis):
        raise TypeError,"inCleanVis MUST be a Python Obit CleanVis"
    #
    return Obit.CleanVisGetName(inCleanVis.me)
    # end PGetName

def PIsA (inCleanVis):
    """ Tells if input really a Python Obit CleanVis

    return true, false (1,0)
    inCleanVis   = Python CleanVis object
    """
    ################################################################
    # Checks
    if inCleanVis.__class__ != CleanVis:
        return 0
    return Obit.CleanVisIsA(inCleanVis.me)
    # end PIsA
