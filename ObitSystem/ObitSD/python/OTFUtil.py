""" OTF Utility module
"""
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

# Python ObitOTF utilities
import Obit, OTF, OErr, FInterpolate, FArray, Image, ImageDesc, InfoList
import Table

def PSubImage(inOTF, outOTF, image, desc, err):
    """ Subtract a 2D ObitFArray from an OTF

    inOTF   = input Python Obit OTF
    outOTF  = output Python Obit OTF, must be previously defined
    image   = Python Obit Image data to subtract, as FArray
    desc    = Python Obit ImageDesc for image
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OTF.PIsA(outOTF):
        raise TypeError,"outOTF MUST be a Python Obit OTF"
    if not FArray.PIsA(image):
        raise TypeError,"image MUST be a Python Obit FArray"
    if not ImageDesc.PIsA(desc):
        raise TypeError,"desc MUST be a Python Obit ImageDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.OTFUtilSubImage(inOTF.me, outOTF.me, image.me, desc.me, err.me)
    # end PSubImage

def PModelImage(inOTF, outOTF, image, desc, err):
    """ Replace OTF data with the model in a 2D ObitFArray

    inOTF   = input Python Obit OTF
    outOTF  = output Python Obit OTF, must be previously defined
    image   = Python Obit Image model to use, as FArray
    desc    = Python Obit ImageDesc for image
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OTF.PIsA(outOTF):
        raise TypeError,"outOTF MUST be a Python Obit OTF"
    if not FArray.PIsA(image):
        raise TypeError,"image MUST be a Python Obit FArray"
    if not ImageDesc.PIsA(desc):
        raise TypeError,"desc MUST be a Python Obit ImageDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.OTFUtilModelImage(inOTF.me, outOTF.me, image.me, desc.me, err.me)
    # end PModelImage

def PScale(inOTF, outOTF, scale, offset, err):
    """ Scale and offset data in an OTF

    out = in*scale + offset
    inOTF   = input Python Obit OTF
    outOTF  = output Python Obit OTF, must be previously defined
    scale   = multiplicative term
    offset  = additive term
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
    if err.isErr: # existing error?
        return
    #
    Obit.OTFUtilScale(inOTF.me, outOTF.me, scale, offset, err.me)
    # end PScale

def PNoise(inOTF, outOTF, scale, offset, sigma, err):
    """ Scale and offset and add Gaussian noise to data in an OTF

    out = in*scale + offset + noise(sigma)
    inOTF   = input Python Obit OTF
    outOTF  = output Python Obit OTF, must be previously defined
    scale   = multiplicative term
    offset  = additive term
    sigma   = Std. deviation of noise to be added
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
    if err.isErr: # existing error?
        return
    #
    Obit.OTFUtilNoise(inOTF.me, outOTF.me, scale, offset, sigma, err.me)
    # end PNoise

def PCreateImage (inOTF, err):
    """ Create basic ObitImage structure and fill out descriptor.

    Imaging parameters are on the inOTF info member.
    "nx"     int scalar Dimension of image in RA [no default].
    "ny"     int scalar Dimension of image in declination[no default]
    "RA"     float scalar Right Ascension of center of image
             Default is observed position center in inOTF 
    "Dec"    float scalar Declination of center of image
             Default is observed position center in inOTF 
    "xCells" float scalar X (=RA) cell spacing in degrees [no default]
    "yCells" float scalar Y (=dec) cell spacing in degrees [no default]
    "Proj"   string (4,1,1) Projection string "-SIN", "-ARC", "-TAN"
             [Default "-SIN"]
    returns Image
    inOTF   = Python Obit OTF from which the solution is to be determined
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    out     = Image.Image("None")
    if err.isErr: # existing error?
        return out
    out.me = Obit.OTFUtilCreateImage (inOTF.me, err.me)
    return out
    # end PCreateImage

def PMakeImage (inOTF,  outImage, err, doBeam=True, Beam=None, Wt=None):
    """ Makes an image from an OTF

    Convolves data onto a grid, accumulated and normalizes.
    Imaging parameters are on the inOTF info member.
    "minWt"  float scalar Minimum summed gridding convolution weight [def 0.1]
    "beamNx" int scalar "X" size of Beam (pixels)
    "beamNy" int scalar "Y" size of Beam (pixels)
    "doScale" bool scalar If true, convolve/scale beam [def True]
              Only use False if providing a dirty beam which already
              includes the effects of gridding.
    "deMode"  bool scalar Subtract image mode from image? [def False]
    "deBias"  bool scalar Subtract calibration bias from image? [def False]
             Note, this doesn't really work the way you would like
    "ConvType" Convolving function Type 0=pillbox,3=Gaussian,4=exp*sinc,5=Sph wave
              def [3]
    "ConvParm"float[10] = Convolving function parameters
    inOTF   = Python Obit input OTF 
    outImage= Python Obit Image
    err     = Python Obit Error/message stack
    doBeam  = if True make "Beam" (PSF) image and use to normalize
    Beam    = Actual instrumental Beam to use, else Gaussian [def None]
    Wt      = Image to save gridding weight array [def None]'),
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not Image.PIsA(outImage):
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return 
    #
    # Dummy image arguments if needed
    if Beam:
        lbeam = Beam   # Actual beam given
    else:
        lbeam = Image.Image("NoBeam")  # Not given
    if Wt:
        lWt = Wt   # Actual weight image given
    else:
        lWt = Image.Image("NoWt")  # Not given
    Obit.OTFUtilMakeImage(inOTF.me, outImage.me, doBeam, lbeam.me, lWt.me, err.me);
    # end PMakeImage


def PIndex (inOTF, err):
    """ Index an OTF

    inOTF   = Python Obit OTF 
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        print "Really is",inOTF.__class__
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.OTFUtilIndex (inOTF.me, err.me)
    # end PIndex

def PFitCal (inOTF, detect, err):
    """ Fits calibrator scans

    Gets calibrator information from the Target table (position, flux density)
    Input values on inOTF
    "Scans"    OBIT_Int (1,1,1) Scan numbers to process, should be 4 scans
    "Tau0"     OBIT_float (1,1,1) Zenith opacity in nepers [def 0].
    "ATemp"    OBIT_float (*,1,1) Effective atmospheric temperature in data units.
        i.e. the additional offset per airmass due to the atmosphere.
        per detector in units of the cal. [def 300]
    "doPlot"   OBIT_Bool (1,1,1) If present and True, plot data and models [def FALSE]
    Output values on inOTF
    "TRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
    "calJy"    OBIT_float (*,1,1) Noise cal value in Jy, per detector
    "RAoff"    OBIT_float (*,1,1) Offset in deg to add to RA
    "Decoff"   OBIT_float (*,1,1) Offset in deg to add to Dec
    "Timeoff"  OBIT_float (*,1,1) Offset in time (day) (actual-expected peak)
    
    inOTF   = Python Obit OTF
    detect  = detector number (0-rel), -1 => all
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.OTFUtilFitCal (inOTF.me, detect, err.me)
    # end PFitCal

def PFitOnOff (inOTF, detect, err):
    """ Fits calibrator On/Off scan pair

   Gets calibrator information from the Target table (flux density)
    Input values on inOTF
    "Scan"     OBIT_Int (2,1,1) Scan numbers to process, should be 4 scans
    Output values on inOTF
    "TRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
    "calJy"    OBIT_float (*,1,1) Noise cal value in Jy, per detector
    
    inOTF   = Python Obit OTF
    detect  = detector number (0-rel), -1 => all
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.OTFUtilFitOnOff (inOTF.me, detect, err.me)
    # end PFitOnOff

def PFitBPOnOff (inOTF, scans, err, BPVer=1):
    """ Fits bandpass from On/Off scan pair

    Gets calibrator information from the Target table (flux density)
    Writes results in OTFBP table
    Input values on inOTF
    inOTF   = Python Obit OTF
    scans   = pair of scans [off,on]
    err     = Python Obit Error/message stack
    BPVer   = output BP table,  0=>new
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.OTFUtilFitBPOnOff (inOTF.me, scans[0], scans[1], BPVer, err.me)
    # end PFitBPOnOff

def PFitTip (inOTF, err):
    """ Fits tipping scans

    Gets reciever/sky brightness/tau0 from tipping scan
    Input values on inOTF
    "Scan"     OBIT_Int (1,1,1) Scan number to process, should be a tipping scan
    "TCal"     OBIT_float (*,1,1) Noise cal value in K, per detector
    "TSky"     OBIT_float (1,1,1) Physical temperature of the sky (K) [def 300]
    "minEl"    OBIT_float (1,1,1) Min. elevation (deg).
                                  This may be needed to cut off low elevations where the
                                  mountains are in the beam.
    "doPlot"   OBIT_Bool (1,1,1) If present and True, plot data and models [def FALSE]
    Output values on inOTF
    "Tau0"     OBIT_float (1,1,1) Zenith opacity in nepers 
    "ATemp"    OBIT_float (*,1,1) Effective atmospheric temperature in data units.
               i.e. the additional offset per airmass due to the atmosphere.
               per detector in units of the cal.
    "TRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
    "TipRMS"   OBIT_float (*,1,1) RMS residual (K) per detector
    
    inOTF   = Python Obit OTF
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.OTFCalUtilFitTip (inOTF.me, err.me)
    # end PFitTip

def PFitNod (inOTF, err):
    """ Fits nodding scans

    Gets cal values in units of Jy
    Input values on inOTF
    "Scan"     OBIT_Int   (1,1,1) Scan number to process, should be a nodding scan
    "calFlux"  OBIT_float (1,1,1) If given, the calibrator flux density, otherwise
                                  get from OTFTarget table
    Output values on inOTF
    "TRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
    "calJy"    OBIT_float (*,1,1) Equivalent Jy per detector of the cal.
    
    inOTF   = Python Obit OTF
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.OTFCalUtilFitNod (inOTF.me, -1, err.me)
    # end PFitNod

def PDiffNod (inOTF, scan, err):
    """ Differences nodding scans

    Differences ons and offs in an nodding scan
    Output values on inOTF
    "OnOff"   OBIT_float (*,1,1) Differences of On-Off for each detector
              in the same order as defined in the data.  For beamswitched
              data this will be twice the source strength.
     
    inOTF   = Python Obit OTF, target position must be in OTFTarget table
              Calibration specified is applied.
    Scan    = scan number to process
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.OTFUtilDiffNod (inOTF.me, scan, err.me)
    # end PDiffNod

def PFlag (inOTF, err,
           timeRange=[0.0,1.0e20], Target="Any", flagVer=1, Chans=[1,0],
           Stokes="1111", Feed=0, Reason=" "):
    """ Adds entry to flag table

    Adds flagging table entry.
    Note: there is currently no sorting so the flagging entries must be
    added in time order
    inOTF     = Python Obit OTF
    err       = Python Obit Error/message stack
    timeRange = pair of floats giving the beginning and end time in days,
                inclusive, of the data to be flagged
    Target    = target name list, "Any" => all sources.
    flagVer   = flagging table version number
    Chans     = pair of ints giving first and last spectral channel numbers
                (1-rel) to be flagged; 0s => all
    Stokes    = String giving stokes to be flagged, 
                "FFF"  where F is '1' to flag corresponding stokes, '0' not.
                Stokes order 'R', 'L', 'RL' or 'X', 'Y', 'XY'
    Feed      = feed number (1-rel) to flag, 0 => all
    Reason    = reason string for flagging (max 24 char)
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    # Set flagging parameters 
    inInfo = OTF.PGetList(inOTF)
    dim = OTF.dim
    dim[0] = 1; dim[1] = 1; dim[2] = 1
    InfoList.PAlwaysPutInt(inInfo, "flagVer", dim, [flagVer])
    InfoList.PAlwaysPutInt(inInfo, "Feed", dim, [Feed])
    dim[0] = 2
    InfoList.PAlwaysPutInt(inInfo, "Chans", dim, Chans)
    dim[0] = 2
    InfoList.PAlwaysPutFloat(inInfo, "timeRange", dim, timeRange)
    dim[0] = len(Target)
    InfoList.PAlwaysPutString(inInfo, "Target", dim, [Target])
    dim[0] = len(Stokes)
    InfoList.PAlwaysPutString(inInfo, "Stokes", dim, [Stokes])
    dim[0] = len(Reason)
    InfoList.PAlwaysPutString(inInfo, "Reason", dim, [Reason])
    #
    Obit.OTFCalUtilFlag (inOTF.me, err.me)
    # end PFlag

def PConvBeam (CCTab, Beam, Template, err):
    """ Convolve a set of Clean components with a beam image

    CCTab      CC Table, following parameters on infoList member
       "BComp" OBIT_int (1,1,1) Start CC to use, 1-rel [def 1 ]
       "EComp" OBIT_int (1,1,1) Highest CC to use, 1-rel [def to end ]
    Beam       Beam Image to convolve with CCs
    Template   Template FArray for output array
    err        Obit Error stack
    Returns an ObitFArray whose size is that of Template, spacing is that of Beam
        with the CCs in CCTab convolved with Beam and the (0,0) position is
        (nx/2,ny/2) (0-rel)
    """
    ################################################################
    # Checks
    if not Table.PIsA(CCTab):
        raise TypeError,"CCTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    out = FArray.FArray("Convolved CCs")
    if err.isErr: # existing error?
        return out
    out.me = Obit.OTFUtilConvBeam (CCTab.me, Beam.me, Template.me, err.me)
    return out
    # end PConvBeam

def PSubModel (inOTF, outOTF, modelImg, beamImg, err):
    """ Subtract a CLEAN sky model from an OTF

    Uses the AIPS CC table on modelImg convolved with the beam shape in beamImg
    to create a sky model which is subtracted from inOTF.
    If outOTF is defined, the subtracted data is written to it, else a scratch
    OTF is returned.
    If modelImg is None then the input data (with any calibration and selection)
    is copied to the output OTF
    Returns model subtracted or copied OTF.
    inOTF    = input Python Obit OTF, any calibration and editing parameters
        on List are applied
    outOTF   = output Python Obit OTF, None for scratch output
    modelImg = CLEAN image with attached AIPS CC table
        CLEAN components should have been scaled to proper units
        following optional parameters on infoList member
        "clip"  OBIT_float (1,1,1) Lowest value to allow in sky model [def. 0]
        "BComp" OBIT_int (1,1,1) Start CC to use, 1-rel [def 1 ]
        "EComp" OBIT_int (1,1,1) Highest CC to use, 1-rel [def to end ]
    beamImg  =  Beam Image to convolve with CCs, unused if modelImg=None
    err      =  Obit Error stack
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
    # set Output
    if outOTF:
        out = outOTF
    else:
        out = OTF.PScratch (inOTF, err)
        OErr.printErrMsg(err, "Error creating scratch OTF")
    # If no model given just copy
    if modelImg==None:
        inOTF.Copy(out,err)
        OErr.printErrMsg(err, "Error copying data")
        return out
    
    # Model check
    if not Image.PIsA(modelImg):
        raise TypeError,"modelImg MUST be a Python Obit Image"
    if not Image.PIsA(beamImg):
        raise TypeError,"beamImg MUST be a Python Obit Image"
    # Sky model by comvolving CCs with beamImg
    CCTab = modelImg.NewTable (Table.READONLY, "AIPS CC", 1, err)
    OErr.printErrMsg(err, "Error accessing iinput AIPS CC table")
    # Copy range of CCs if given
    comp = modelImg.List.get("BComp")
    if comp[0]==0:
        CCTab.List.set("BComp", comp[4])
    comp = modelImg.List.get("EComp")
    if comp[0]==0:
        CCTab.List.set("EComp", comp[4])

    # Make image
    model = PConvBeam (CCTab, beamImg, modelImg.FArray, err)

    # clip model as requested
    clip = modelImg.List.get("clip")
    if clip[0]==0:
        minFlux = clip[4][0]
    else:
        minFlux = 0.0
    print "Clip model below ",  minFlux,
    FArray.PClip (model, minFlux, 1.0e20, 0.0)

    # Subtract model
    PSubImage(inOTF, out, model, modelImg.Desc, err)
    OErr.printErrMsg(err, "Error subtracting sky model")

    # Cleanup
    del model

    return out
    # end PSubModel

def PEditFD(inOTF, outOTF, err):
    """ Frequency domain editing of an OTF

    Compares data with a running median in frequency and makes flag table
    entries for channels with discrepant values or excessive RMSes.
    Data in channels which are mostly, but not completely flagged may also
    be flagged.
    Generates an OTFFlag flagging table on outOTF.
    Multiple threads may be used, up to one per spectrum.
    Control parameters on inOTF List member:
      flagTab  int   FG table version number [ def. 1]
      timeAvg  float Time interval over which to average (min) [def. 1]
      FDwidMW  int   The width of the median window in channels. 
                     An odd number (5) is recommended, [def nchan-1]
      FDmaxRMS float[2] Flag all channels having RMS values > FDmaxRMS[0] 
                        of the channel median sigma.[def 6.]
                        plus FDmaxRMS[1] (def 0.1) of the channel 
                        average in quadrature.
      FDmaxRes float Max. residual flux in sigma allowed [def 6.]
      minGood  float Minimum fraction of time samples at and below 
               which a channel/interval will be flagged.  
               [def 0.25, -1 ->no flagging by fraction of samples.]
    inOTF   = input Python Obit OTF, prior editing and calibration honored
    outOTF  = output Python Obit OTF for flag table
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
    if err.isErr: # existing error?
        return
    #
    Obit.OTFFlagEditFD(inOTF.me, outOTF.me, err.me)
    # end PEditFD

