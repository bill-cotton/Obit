""" Python Obit SkyModelVMBeam class

This class contains a sky model using an image of the antenna
response and can Fourier transform it and subtract from or
divide into a UV data.
Both FITS and AIPS cataloged files are supported.

SkyModelVMBeam Members with python interfaces:
List      - used to pass instructions to processing
Mosaic    - ImageMosaic
"""
# $Id: SkyModelVMBeam.py 2 2008-06-10 15:32:27Z bill.cotton $
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

# Obit SkyModelVMBeam
import Obit, OErr, ImageMosaic, InfoList, UV, SkyModel

# Python shadow class to ObitSkyModelVMBeam class
 
# class name in C
myClass = "ObitSkyModelVMBeam"
 
#
class SkyModelVMBeam(SkyModel.SkyModel):
    """ Python Obit SkyModelVMBeam class
    
    This class contains an ionospheric sky model and can Fourier
    transform it and subtract from or divide into a UV data.
    Both FITS and AIPS cataloged files are supported.
    
    SkyModelVMBeam Members with python interfaces:
    List      - used to pass instructions to processing
    Mosaic    - ImageMosaic
    """
    def __init__(self, name) :
        self.this = Obit.new_SkyModelVMBeam(name)
        self.myClass = myClass
    def __del__(self):
        if Obit!=None:
            Obit.delete_SkyModelVMBeam(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.SkyModelVMBeamUnref(Obit.SkyModelVMBeam_me_get(self.this))
            # In with the new
            Obit.SkyModelVMBeam_me_set(self.this,value)
            return
        if name=="Mosaic":
            Obit.SkyModelVMSetImageMosaic(self.me, value.me)
            return
        self.__dict__[name] = value
    def __getattr__(self, name):
        if self.__class__ != SkyModelVMBeam:
            return
        if name == "me" : 
            return Obit.SkyModelVMBeam_me_get(self.this)
        if name=="Mosaic":
            return SkyModel.PGetMosaic(self)
        if name=="List":
            return SkyModel.PGetList(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != SkyModelVMBeam:
            return
        return "<C SkyModelVMBeam instance> " + Obit.SkyModelVMBeamGetName(self.me)
    def cast(self, toClass):
        """ Casts object pointer to specified class
        
        self     = object whose cast pointer is desired
        toClass  = Class string to cast to
        """
        ################################################################
        # Get pointer with type of this class
        out =  self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
    
# end class SkyModelVMBeam

def PCreate (name, mosaic, uvData, IBeam, VBeam, QBeam, UBeam, \
             IBeamPh, VBeamPh, QBeamPh, UBeamPh, err):
    """ Create the parameters and underlying structures of a SkyModelVMBeam.
    
    name      = Name to be given to object
    Most control parameters are in InfoList member
    mosaic    = Python ImageMosaic to attach
    uvData    = Python UV data to be operated on
    IBeam     = Python Image I Beam normalized image
    VBeam     = Python Image V Beam fractional image
    QBeam     = Python Image Q Beam fractional image
    UBeam     = Python Image U Beam fractional image
    IBeamPh   = Python Phase Image I Beam 
    VBeamPh   = Python Phase Image V Beam 
    QBeamPh   = Python Phase Image Q Beam 
    UBeamPh   = Python Phase Image U Beam 
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not ImageMosaic.PIsA(mosaic):
        raise TypeError,"uvData MUST be a Python Obit UV"
    #
    out = SkyModelVMBeam("None");
    out.me = Obit.SkyModelVMBeamCreate(name, mosaic.me, uvData.me, \
                                       IBeam.me, VBeam.me, QBeam.me, UBeam.me, \
                                       IBeamPh.me, VBeamPh.me, QBeamPh.me, UBeamPh.me, \
                                       err.me)
    return out;
    # end PCreate

# Structure for passing inputs, add to SkyModel version
UVSubInput = SkyModelVMBeam.cUVSubInput

def PSubUV (err, input=UVSubInput):
    """ Fourier transform Sky model and subtract from uv data

    A SkyModelVMBeam is Fourier transformed and subtracted 
    from InData and written to outData.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
    Input dictionary entries:
    InData      = Input UV data,
    SkyModelVMBeam    = Input SkyModelVMBeam,
    OutData     = Output uv data,
    doCalSelect = Select/calibrate/edit data?
    REPLACE     = Replace data with model?
    Stokes      = Stokes parameter, blank-> unchanged from input),
    doPhase     = If Phase images provided. [def false]
    CCVer       = CC table versions to use [def all 0 => highest]
    BComp       = Start CC to use per table, 1-rel [def 1 ]
    EComp       = Highest CC to use per table, 1-rel [def to end]
    BChan       = First spectral channel selected. [def all]),
    EChan       = Highest spectral channel selected. [def all]),
    BIF         = First IF selected. [def all]),
    EIF         = Highest IF selected. [def all]),
    doPol       = >0 -> calibrate polarization.),
    doCalib     = >0 -> calibrate, 2=> also calibrate Weights),
    flagVer     = Flag table version, 0-> use highest, <0-> none),
    BLVer       = BL table version, 0> use highest, <0-> none),
    BPVer       = Band pass (BP) table version, 0-> use highest),
    Subarray    = Selected subarray, <=0->all [default all]),
    freqID      = Selected Frequency ID, <=0->all [default all]),
    timeRange   = Selected timerange in days. [8 floats] 0s -> all),
    UVRange     = Selected UV range in wavelengths. 0s -> all),
    Sources     = Source names selected unless any starts with),
    Antennas    = A list of selected antenna numbers, if any is negative),
    corrType    = Correlation type, 0=cross corr only, 1=both, 2=auto only.),
    doBand      = Band pass application type <0-> none),
    Smooth      = Specifies the type of spectral smoothing [three floats]
    do3D        = If 3D imaging wanted. [def false]
    Factor      = Model multiplication factor (-1=>add) [def 1]
    PBCor       = If TRUE make relative primary beam corrections. [def false]
    antSize     = Diameter of antennas for PBCor,.[def 25.0]
    minFlux     = Minimum flux density model or pixel [def -1.0e20]
    Type        = Model type (ObitSkyModelVMBeamType) [def OBIT_SkyModelVMBeam_Comps]
                  0=CC Comps, 1=Image, 2=Model
    Mode        = Model mode (ObitSkyModelVMBeamMode) [def OBIT_SkyModelVMBeam_Fastest]
                  0=fastest, 1=DFT, 2=Grid
    MODPTFLX    = Point model flux in Jy, [def 0.0]')
    MODPTXOF    = Point model x offset in deg  [def 0.0]
    MODPTYOF    = Point model y offset in deg  [def 0.0]
    MODPTYP     = Point other parameters  [def all 0.0]
    """
    ################################################################
    # Get input parameters
    inData            = input["InData"]
    inSkyModelVMBeam  = input["SkyModelVMBeam"]
    outData           = input["OutData"]
    # Checks
    if not PIsA(inSkyModelVMBeam):
        raise TypeError,"inSkyModelVMBeam MUST be a Python Obit SkyModelVMBeam"
    if not UV.PIsA(inData):
        raise TypeError,"inData MUST be a Python Obit UV"
    if not UV.PIsA(outData):
        raise TypeError,"outData MUST be a Python Obit UV"
    #
    dim = [1,1,1,1,1]
    #
    # Set control values on SkyModelVMBeam/inData
    dim[0] = 1;
    inInfo = inSkyModelVMBeam.List
    uvInfo = inData.List   # 
    input["SkyModel"] = inSkyModelVMBeam
    #
    # Do operation
    SkyModel.PSubUV(err, input=input)
    if err.isErr:
        OErr.printErrMsg(err, "Error subtraction SkyModelVMBeam from UV data")
    # end PSubUV

# Divide a SkyModelVMBeam into a UV data set
# Structure for passing inputs, add to SkyModel version
UVDivInput = SkyModelVMBeam.cUVDivInput

def PDivUV (err, input=UVSubInput):
    """ Fourier transform Sky model and divide into uv data

    A SkyModelVMBeam is Fourier transformed and divided into
    InData and written to outData.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
    Input dictionary entries:
    InData      = Input UV data,
    SkyModelVMBeam    = Input SkyModelVMBeam,
    OutData     = Output uv data,
    doCalSelect = Select/calibrate/edit data?),
    doPhase     = If Phase images provided. [def false]
    REPLACE     = Replace data with model?
    Stokes      = Stokes parameter, blank-> unchanged from input),
    CCVer       = CC table versions to use [def all 0 => highest]
    BComp       = Start CC to use per table, 1-rel [def 1 ]
    EComp       = Highest CC to use per table, 1-rel [def to end]
    BChan       = First spectral channel selected. [def all]),
    EChan       = Highest spectral channel selected. [def all]),
    BIF         = First IF selected. [def all]),
    EIF         = Highest IF selected. [def all]),
    doPol       = >0 -> calibrate polarization.),
    doCalib     = >0 -> calibrate, 2=> also calibrate Weights),
    flagVer     = Flag table version, 0-> use highest, <0-> none),
    BLVer       = BL table version, 0> use highest, <0-> none),
    BPVer       = Band pass (BP) table version, 0-> use highest),
    Subarray    = Selected subarray, <=0->all [default all]),
    freqID      = Selected Frequency ID, <=0->all [default all]),
    timeRange   = Selected timerange in days. [8 floats] 0s -> all),
    UVRange     = Selected UV range in wavelengths. 0s -> all),
    Sources     = Source names selected unless any starts with),
    Antennas    = A list of selected antenna numbers, if any is negative),
    corrType    = Correlation type, 0=cross corr only, 1=both, 2=auto only.),
    doBand      = Band pass application type <0-> none),
    Smooth      = Specifies the type of spectral smoothing [three floats]
    do3D        = If 3D imaging wanted. [def false]
    Factor      = Model multiplication factor (-1=>add) [def 1]
    PBCor       = If TRUE make relative primary beam corrections. [def false]
    antSize     = Diameter of antennas for PBCor,.[def 25.0]
    minFlux     = Minimum flux density model or pixel [def -1.0e20]
    Type        = Model type (ObitSkyModelVMBeamType) [def OBIT_SkyModelVMBeam_Comps]
                  0=CC Comps, 1=Image, 2=Model
    Mode        = Model mode (ObitSkyModelVMBeamMode) [def OBIT_SkyModelVMBeam_Fastest]
                  0=fastest, 1=DFT, 2=Grid
    MODPTFLX    = Point model flux in Jy, [def 0.0]')
    MODPTXOF    = Point model x offset in deg  [def 0.0]
    MODPTYOF    = Point model y offset in deg  [def 0.0]
    MODPTYP     = Point other parameters  [def all 0.0]
    """
    ################################################################
    # Get input parameters
    inData           = input["InData"]
    inSkyModelVMBeam  = input["SkyModelVMBeam"]
    outData          = input["OutData"]
    # Checks
    if not PIsA(inSkyModelVMBeam):
        raise TypeError,"inSkyModelVMBeam MUST be a Python Obit SkyModelVMBeam"
    if not UV.PIsA(inData):
        raise TypeError,"inData MUST be a Python Obit UV"
    if not UV.PIsA(outData):
        raise TypeError,"outData MUST be a Python Obit UV"
    #
    #
    dim = [1,1,1,1,1]
    #
    # Set control values on SkyModelVMBeam/inData
    dim[0] = 1;
    inInfo = inSkyModelVMBeam.List
    uvInfo = inData.List   # 
    input["SkyModel"] = inSkyModelVMBeam
    #
    # Do operation
    SkyModel.PDivUV(err, input=input)
    if err.isErr:
        OErr.printErrMsg(err, "Error dividing SkyModelVMBeam into UV data")
    # end PDivUV

def PIsA (inSkyModel):
    """ Tells if input really a Python Obit SkyModelVMBeam

    return true, false (1,0)
    inSkyModel   = Python SkyModelVMBeam object
    """
    ################################################################
    # Checks - allow inheritence
    if not str(inSkyModel.__class__).startswith("SkyModelVMBeam"):
        return 0
    sm = inSkyModel.cast(myClass)  # cast pointer
    return Obit.SkyModelVMBeamIsA(sm)
    # end PIsA
