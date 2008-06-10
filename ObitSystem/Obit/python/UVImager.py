# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2007
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

# Interferometric imaging class
import Obit, OErr, ImageMosaic, InfoList, UV, types
import ImageDesc

# Python shadow class to ObitUVImager class
 
class UVImagerPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.UVImagerUnref(Obit.UVImager_me_get(self.this))
            # In with the new
            Obit.UVImager_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != UVImager:
            return
        if name == "me" : 
            return Obit.UVImager_me_get(self.this)
        # Virtual members
        if name=="mosaic":
            return PGetMosaic(self)
        if name=="uvdata":
            return PGetUV(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != UVImager:
            return
        return "<C UVImager instance> " + Obit.UVImagerGetName(self.me)
#
class UVImager(UVImagerPtr):
    """ Python Obit Image class

    This class contains gives access to the functions needed to convert
    an Obit UV into one or more dirty images.
    Both FITS and AIPS cataloged data are supported.

    UVImager Members with python interfaces:
    uvdata    - input UV data, use PGetUV
    mosaic    - output Image mosaic, use PGetMosaic
    """
    def __init__(self, name, uvdata, err) :
        self.this = Obit.new_UVImager(name, uvdata, err)
    def __del__(self):
        if Obit!=None:
            Obit.delete_UVImager(self.this)
    # end class definition

# Commonly used, dangerous variables
dim=[1,1,1,1,1]


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

# Create an UVImager
UVCreateImagerInput={
    'structure':['UVCreateImager',[('InData',  'Input UV data'),
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
                                   ('Robust',   'Briggs robust parameter. (AIPS definition)'),
                                   ('UVTaper',  'UV plane taper, sigma in klambda as [u,v]'),
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
                                   ('OutlierSize', 'Size of outlyer field (pixels)')]],
    # defaults
    'InData':None,
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
    'Robust': 0.0,                                   
    'UVTaper': [0.0,0.0],
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
    'OutlierSize':50 }
def PUVCreateImager (err, name = 'myUVImager', input=UVCreateImagerInput):
    """ Create an imager to generate a set of images needed to cover a region.

    Create and return a Python Obit UVImager based on input parameters and
    uv data to be imaged.  Images should be fully defined when returned.
    err     = Python Obit Error/message stack
    name    = Name for output object
    input   = input parameter dictionary
    
    Input dictionary entries:
    InData   = Input Python UV data to image
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
    Robust   = Briggs robust parameter. (AIPS definition)
    UVTaper  = UV plane taper, sigma in klambda as [u,v]
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

    returns = Output UVImager with defined ImageMosaic
    """
    ################################################################
    # Get input parameters
    InData   = input["InData"]
    #
    # Checks
    if not UV.PIsA(InData):
        raise TypeError, 'UVCreateImage: Bad input UV data'
    # Set control values on UV 
    dim[0] = 1;
    inInfo = UV.PGetList(InData)  # Add control to UV data
    # Calibration/editing/selection
    InfoList.PPutBoolean (inInfo, "doCalSelect",  dim, [input["doCalSelect"]], err)
    dim[0] = 4;
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
    InfoList.PPutFloat  (inInfo, "Robust", dim, [input["Robust"]],  err)
    InfoList.PPutInt    (inInfo, "WtBox",  dim, [input["WtBox"]],   err)
    InfoList.PPutInt    (inInfo, "WtFunc", dim, [input["WtFunc"]],  err)
    InfoList.PPutFloat  (inInfo, "WtPower",dim, [input["WtPower"]], err)
    dim[1] = 2
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
    InfoList.PPutInt    (inInfo, "imFileType",dim, [input["Type"]],   err)
    InfoList.PPutInt    (inInfo, "imSeq",    dim, [input["Seq"]],    err)
    InfoList.PPutInt    (inInfo, "imDisk",   dim, [input["Disk"]],   err)
    InfoList.PPutFloat  (inInfo, "FOV",      dim, [input["FOV"]],    err)
    InfoList.PPutBoolean (inInfo, "doFull",  dim, [input["doFull"]], err)
    InfoList.PPutInt    (inInfo, "NField",   dim, [input["NField"]], err)
    InfoList.PPutFloat  (inInfo, "xCells",   dim, [input["xCells"]], err)
    InfoList.PPutFloat  (inInfo, "yCells",   dim, [input["yCells"]], err)
    InfoList.PPutFloat  (inInfo, "OutlierFlux", dim, [input["OutlierFlux"]], err)
    InfoList.PPutFloat  (inInfo, "OutlierDist", dim, [input["OutlierDist"]], err)
    InfoList.PPutFloat  (inInfo, "OutlierSI",   dim, [input["OutlierSI"]],   err)
    InfoList.PPutInt    (inInfo, "OutlierSize", dim, [input["OutlierSize"]], err)
    dim[0] = len(input["Name"])
    InfoList.PAlwaysPutString (inInfo, "imName",   dim, [input["Name"]])
    dim[0] = len(input["Class"])
    InfoList.PAlwaysPutString (inInfo, "imClass",  dim, [input["Class"]])
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
    #
    # Create
    out = UVImager(name, InData.me, err.me)
    # show any errors 
    #OErr.printErrMsg(err, "UVImage: Error creating UVImager object")
    #
    return out
    # end UVCreateImage

def PCopy (inImager, outImager, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inImager, copies pointers
    (probably not very useful)
    inImager  = Python UVImager object to copy
    outImager = Output Python UVImager object, must be defined
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inImager):
        raise TypeError,"inImager MUST be a Python Obit UVImager"
    if not PIsA(outImager):
        raise TypeError,"outImager MUST be a Python Obit UVImager"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.UVImageCopy (inImager.me, outImager.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying UVImage")
    # end PCopy


def PGetMosaic (InImager, err):
    """ Return the member ImageMosaic

    returns ImageMosaic
    InImager  = Python ImageMosaic object
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(InImager):
        raise TypeError,"InImager MUST be a Python Obit UVImager"
    #
    out    = ImageMosaic.ImageMosaic("None", 1)
    out.me = Obit.UVImagerGetMosaic(InImager.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error obtaining Mosaic member")
    return out
    # end PGetMosaic

def PGetUV (InImager):
    """ Return the member uvdata

    returns uv data (as UV)
    InImager  = Python ImageMosaic object
    """
    ################################################################
    # Checks
    if not PIsA(InImager):
        raise TypeError,"InImager MUST be a Python Obit UVImager"
    #
    out    = UV.UV("None")
    out.me = Obit.UVImagerGetUV(InImager.me)
    return out
    # end PGetMosaic


# Define UVWeight input dictionary
UVWeightInput={'structure':['UVWeight',[('InImager','Input UVImager'),
                                        ('Robust',   'Briggs robust parameter. (AIPS definition)'),
                                        ('UVTaper',  'UV plane taper, sigma in klambda as [u,v]'),
                                        ('WtSize',   'Size of weighting grid in cells [image]'),
                                        ('WtBox',    'Size of weighting box in cells [def 0]'),
                                        ('WtFunc',   'Weighting convolution function [def. 1]'),
                                        ('WtPower',  'Power to raise weights to.  [def = 1.0]'),
                                        ('Channel',  'Channel (1-rel) number to image, 0-> all.')]],
               # defaults
               'InImager':None,
               'DoWeight': True,                                   
               'Robust': 0.0,                                   
               'UVTaper': [0.0,0.0],
               'WtSize':-1,
               'WtBox':0,
               'WtFunc':1,
               'WtPower':1.0,
               'Channel': 0 }                                   
def PUVWeight (InImager, err, input=UVWeightInput):
    """ Apply any calibration/editing/selection and do UV weighting

    UV Data is calibrated and weighted, values in input override those given
    when the UVImage was created.  This operation is optionally done in PUVImage.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
    Input dictionary entries:
    InImager   = Input Python UVImager to image
    Robust   = Briggs robust parameter. (AIPS definition)
    UVTaper  = UV plane taper, sigma in klambda as [u,v]
    WtSize   = Size of weighting grid in cells [same as image nx]
    WtBox    = Size of weighting box in cells [def 1]
    WtFunc   = Weighting convolution function [def. 1]
               1=Pill box, 2=linear, 3=exponential, 4=Gaussian
               if positive, function is of radius, negative in u and v.
    WtPower  = Power to raise weights to.  [def = 1.0]
               Note: a power of 0.0 sets all the output weights to 1 as modified
               by uniform/Tapering weighting.
               Applied in determinng weights as well as after.
    Channel  = Channel (1-rel) number to image, 0-> all.
    """
    ################################################################
    # Get input parameters
    Channel  = input["Channel"]
    WtSize   = input["WtSize"]
    #
    # Checks
    if not PIsA(InImager):
        raise TypeError, 'PUVWeight: Bad input UVImager'

    # Set control values on UV 
    dim[0] = 1;
    inInfo = UV.PGetList(PGetUV(InImager))
    InfoList.PPutFloat  (inInfo, "Robust", dim, [input["Robust"]],  err)
    InfoList.PPutInt    (inInfo, "WtBox",  dim, [input["WtBox"]],   err)
    InfoList.PPutInt    (inInfo, "WtFunc", dim, [input["WtFunc"]],  err)
    InfoList.PPutFloat  (inInfo, "WtPower",dim, [input["WtPower"]], err)
    dim[1] = 2
    InfoList.PPutFloat  (inInfo, "Taper",  dim, input["UVTaper"],   err)
    if (WtSize>0):
        print "WtSize", WtSize
        # Change name for C routine.
        dim[0] = 1;
        InfoList.PPutInt  (inInfo, "nuGrid",  dim, [WtSize], err.me)
        InfoList.PPutInt  (inInfo, "nvGrid",  dim, [WtSize], err.me)
    #
    Obit.UVImagerWeight (InImager.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error weighting UV data")
    # show any messages, bail if errors 
    #OErr.printErrMsg(err, "PUVWeight: Error weighting UV data")
    #
    # end PUVWeight

def PUVImage (InImager, err, field=0, doWeight=False, doBeam=True, doFlatten=False):
    """ Form image optionally doing the calibration/weighting and flatten stages

    UV Data is imaged to form dirty beam. Control parameters are those on UV data
    passed at UVImager creation.
    InImager  = UVImager object
    field     = field number (1-rel) to image, 0=> all
    err       = Python Obit Error/message stack
    doWeight  = Do weighting (Uniform, tapering,...) before imaging
                If True then input data is modified.
    doBeam    = Make Beam before image, must be done once.
    doFlatten = Flatten image mosaic onto single image
    """
    ################################################################
    #
    # Checks
    if not PIsA(InImager):
        print "input class actually ",InImager.__class__ 
        raise TypeError, 'PUVImage: Bad input UVImager'

    #
    Obit.UVImagerImage (InImager.me, field, doWeight, doBeam, doFlatten, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error Imaging data")
    # show any messages, bail if errors 
    #OErr.printErrMsg(err, "PUVImage: Error imaging UV data")
    #
    # end PUVImage

def PUVFlatten (InImager, err):
    """ Flatten ImageMosaic member

    The images of the mosaic member are flattened to a single image.
    This operation is optionally done in PUVImage.          
    InImager  = UVImager object
    err       = Python Obit Error/message stack
    """
    ################################################################
    #
    # Checks
    if not PIsA(InImager):
        raise TypeError, 'PUVFlatten : Bad input UVImager'

    #
    Obit.UVImagerFlatten (InImager.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error flattening mosaic")
    # show any messages, bail if errors 
    #OErr.printErrMsg(err, "PUVFlatten: Error imaging UV data")
    #
    # end PUVFlatten 

def PIsA (InImager):
    """ Tells if input really a Python Obit UVImager

    return true, false (1,0)
    InImager = Python UVImager object
    """
    ################################################################
    # Checks
    if InImager.__class__ != UVImager:
        return 0
    return Obit.UVImagerIsA(InImager.me)
    # end PIsA
