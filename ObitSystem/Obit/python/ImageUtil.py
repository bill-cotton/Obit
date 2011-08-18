""" Python Obit Image utility module
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2011
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

# Python interface to ObitImageUtil utilities
import Obit, Image, ImageDesc, FArray, UV, Table, TableUtil, History, OErr
import OSystem

def PICreateImage (inUV, fieldNo, doBeam, err):
    """
    Create an image from information on an ObitUV
    
    returns  Python Image

    * inUV     = Python UV object, following read from InfoList member 

       ======== ================== =============================================
       "nChAvg" OBIT_int (1,1,1)   number of channels to average.
                                   This is for spectral line observations and
                                   is ignored if the IF axis on the uv data has
                                   more than one IF.  Default is continuum =
                                   average all freq/IFs. 0=> all.
       "rotate" OBIT_float (?,1,1) Desired rotation on sky (from N thru E) in
                                   deg. [0]
       "nx"     OBIT_int (?,1,1)   Dimension of image in RA [no default].
                                   This and the following are arrays with one
                                   entry per field.
       "nxBeam" OBIT_int (?,1,1)   Dimension of beam in RA, [def. nx]
       "ny"     OBIT_int (?,1,1)   Dimension of image in declination[no default]
       "nyBeam" OBIT_int (?,1,1)   Dimension of beam in declination, [def. ny]
       "xCells" OBIT_float (?,1,1) X (=RA) cell spacing in degrees [no default]
       "yCells" OBIT_float (?,1,1) Y (=dec) cell spacing in degrees [no default]
       "xShift" OBIT_float (?,1,1) Desired shift in X (=RA) in degrees. [0]
       "yShift" OBIT_float (?,1,1) Desired shift in Y (=dec) in degrees. [0]
       "nuGrid" OBIT_int (1,1,1)   Size in pixels of weighting grid for uniform
                                   weighting
       ======== ================== =============================================

    * fieldNo  = Which field (1-rel) in imaging parameter arrays.
    * doBeam   = if TRUE also create beam as the myBeam member of returned image.
    * err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    out    = Image("None")
    out.me = Obit.ImageUtilCreateImage (inUV.me, fieldNo, doBeam, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error Creating Image")
    return out
    # end PCreateImage

def PMakeImage (inUV, outImage, channel, doBeam, doWeight, err):
    """
    Grids UV, FFTs and makes corrections for the gridding convolution.

    * inUV     = Input Python uv data. Should be in form of Stokes to be imaged
      will all calibration and selection applied.
    * outImage = Python Image to be written.  Must be previously instantiated.
      Beam normalization factor is written to output Beam infoList as SUMWTS
    * channel  = Which frequency channel to image, 0->all.
    * doBeam   = if TRUE also make beam.  Will make the myBeam member of outImage.
      If FALSE, and myGrid->BeamNorm 0.0 then reads SUMWTS value from beam infolist
    * doWeigh  = if TRUE Apply uniform weighting corrections to uvdata before imaging
    * err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageUtilMakeImage(inUV.me, outImage.me, channel, doBeam, doWeight, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating Image from UV data")
    # end PMakeImage

def PInterpolateImage (inImage, outImage, err, 
                       inPlane=[1,1,1,1,1], outPlane=[1,1,1,1,1], hwidth=2):
    """
    Interpolates one image onto another's grid.
    
    Pixels in outImage are interpolated from inImage

    * inImage  = Input Python Image.
    * outImage = Python Image to be written.  Must be previously instantiated.
    * err      = Python Obit Error/message stack
    * inPlane  = 5 element int array with 1, rel. plane number [1,1,1,1,1]
      giving location of plane to be interpolated
    * outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
      giving location of plane to be written
    * hwidth   = half width of interpolation kernal [1-4] default 2
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    if len(inPlane) != 5:
        raise TypeError,"inPlane must have 5 elements"
    if len(outPlane) != 5:
        raise TypeError,"outPlane must have 5 elements"
    #
    Obit.ImageUtilInterpolateImage(inImage.me, outImage.me,
                                   inPlane, outPlane, hwidth, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error interpolating Image")
    # end PInterpolateImage

def PPBApply (inImage, pntImage, outImage, err,
              inPlane=[1,1,1,1,1], outPlane=[1,1,1,1,1], antSize=25.0):
    """
    Multiply an image by the primary beam pattern of another
    
    Pixels in outImage are inImage multiplied by the antenna beam pattern
    from pntImage

    * inImage  = Input Python Image.
    * pntImage = Python Image giving pointing position (ObsRA, ObsDec)
    * outImage = Python Image to be written.  Must be previously instantiated.
    * err      = Python Obit Error/message stack
    * inPlane  = 5 element int array with 1, rel. plane number [1,1,1,1,1]
      giving location of plane to be interpolated
    * outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
      giving location of plane to be written
    * antSize  = Antenna diameter assumed, default 25m
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(pntImage):
        print "Actually ",pntImage.__class__
        raise TypeError,"pntImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if len(inPlane) != 5:
        raise TypeError,"inPlane must have 5 elements"
    if len(outPlane) != 5:
        raise TypeError,"outPlane must have 5 elements"
    #
    Obit.ImageUtilPBApply(inImage.me, pntImage.me, outImage.me,
                          inPlane, outPlane, antSize, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error applying Primary beam correction to Image")
    # end PPBApply

def PPBImage (pntImage, outImage, err,
              minGain=0.1, outPlane=[1,1,1,1,1], antSize=25.0):
    """
    Calculate an image with a primary beam pattern
    
    Make an image of the antenna primary beam pattern based on the pointing
    position in an image.

    * pntImage = Python Image giving pointing position (ObsRA, ObsDec)
    * outImage = Python Image to be written.  Must be previously instantiated.
    * err      = Python Obit Error/message stack
    * minGain  = minimum allowed gain (lower values blanked).
    * outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
      giving location of plane to be written
    * antSize  = Antenna diameter assumed, default 25m
    """
    ################################################################
    # Checks
    if not Image.PIsA(pntImage):
        print "Actually ",pntImage.__class__
        raise TypeError,"pntImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if len(outPlane) != 5:
        raise TypeError,"outPlane must have 5 elements"
    #
    Obit.ImageUtilPBImage(outImage.me, outImage.me,
                          outPlane, antSize, minGain, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error with primary beam image")
    # end PPBImage

def PPBCorr (inImage, pntImage, outImage, err,
             inPlane=[1,1,1,1,1], outPlane=[1,1,1,1,1], antSize=25.0):
    """
    Correct (divide) an image by the primary beam pattern of another
    
    Pixels in outImage are inImage multiplied by the antenna beam pattern
    from pntImage

    * inImage  = Input Python Image.
    * pntImage = Python Image giving pointing position (ObsRA, ObsDec)
    * outImage = Python Image to be written.  Must be previously instantiated.
    * err      = Python Obit Error/message stack
    * inPlane  = 5 element int array with 1, rel. plane number [1,1,1,1,1]
      giving location of plane to be interpolated
    * outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
      giving location of plane to be written
    * antSize  = Antenna diameter assumed, default 25m
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(pntImage):
        print "Actually ",pntImage.__class__
        raise TypeError,"pntImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if len(inPlane) != 5:
        raise TypeError,"inPlane must have 5 elements"
    if len(outPlane) != 5:
        raise TypeError,"outPlane must have 5 elements"
    #
    Obit.ImageUtilPBCorr(inImage.me, pntImage.me, outImage.me,
                         inPlane, outPlane, antSize, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error making primary beam correction")
    # end PPBCorr

def PScaleImage (inImage, scale, err):
    """
    Scale the pixel values in an image
    
    Scale image, optionally by plane, scales any CC tables, writes history

    * inImage = Obit Python Image
    * scale   = Scaling factor, if scalar, multiply all pixels,
      otherwise one value per image plane.
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Open image
    Image.POpen (inImage, Image.READWRITE, err)
    #  Get input descriptor to see how many planes
    inDesc = Image.PGetDesc(inImage)
    inDescDict = ImageDesc.PGetDict(inDesc)
    ndim  = inDescDict["naxis"]
    inNaxis = inDescDict["inaxes"]
    # Work buffer
    inImageArray = Image.PGetFArray(inImage)
    ImageBuffer = FArray.PCopy(inImageArray, err)
    # Reset max/min
    inDescDict["minval"] =  1.0e20
    inDescDict["maxval"] = -1.0e20
    inImage.Desc.Dict = inDescDict               # Update descriptor on image
    Image.PDirty(inImage)                        # Force update
    Image.PClose (inImage, err)

    # list of planes to loop over (0-rel)
    if (ndim>0) and (inNaxis[2]>0):  
        planes = range(inNaxis[2])
    else:
        planes = [0]
    
    # Loop over planes
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        # Get image plane
        Image.PGetPlane (inImage, ImageBuffer, doPlane, err)

        # Scaling factor
        if type(scale)==list:
            scl = scale[iPlane]
        else:
            scl = scale

        # Scale
        FArray.PSMul(ImageBuffer, scl)

        # Write output
        Image.PPutPlane (inImage, ImageBuffer, doPlane, err)

        # Scale any CC table
        highVer = Image.PGetHighVer (inImage, "AIPS CC")
        if (iPlane+1<=highVer):
            CCTab = Image.PImageGetTable (inImage, Image.READWRITE, "AIPS CC", iPlane+1, err)
            PCCScale (CCTab, 1, 0, scl, err)
        # end loop over planes
    # Write history
    inHistory  = History.History("history", inImage.List, err)
    # Add this programs history
    inHistory.Open(History.READWRITE, err)
    inHistory.TimeStamp(" Start Obit ScaleImage",err)
    if type(scale)==list:
        i = -1
        for iPlane in planes:
            i = i + 1
            scl = scale[i]
            inHistory.WriteRec(-1,"ScaleImage / scale["+str(i+1)+"] = "+str(scl),err)
    else:
        inHistory.WriteRec(-1,"ScaleImage / scale = "+str(scale),err)
    inHistory.Close(err)
# end PScaleImage

def PCCScale (inCCTab, startComp, endComp, scale, err):
    """
    Scale flux densities in a CC table
    
    Flux densities of CC entries startComp through endComp are scales by scle

    * inCCTab   = Input Python TableCC
    * startComp = first (1-rel) component
    * endComp   = highest component [1-rel= 0=> all
    * scale     = flux density scaling factor
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Table.PIsA(inCCTab):
        raise TypeError,"inCCTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageUtilCCScale(inCCTab.me, startComp, endComp, scale, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error scaling CC table")
    # end PCCScale

def PUVFilter (inImage, outImage, radius, err):
    """
    Filter an image outside of a radius from the origin in FT space.
    
    Intended to filter out out of band noise in single dish images.
    Filters by a function with 1.0/(nx*ny) inside radius and outside tapers
    by an exponential with scale distance 10 pixels.

    * inImage   = Input Image
    * outImage  = Output image, may be inImage
    * radius    = distance from origin in uv space (m)
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageUtilUVFilter(inImage.me, outImage.me, radius, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error UV Filtering image")
    # end PUVFilter

def PImageAdd (in1Image, in2Image, outImage, err, \
               chkPos=False, factor1=1.0, factor2=1.0):
    """
    Adds Pixels in in2Image from in1Image and write to outImage
    
    Adds scaled pixel values, writes history

    * in1Image = input Obit Python Image 1
    * in2Image = input Obit Python Image 2
    * outImage = output Obit Python Image, must be defined but not instantiated
    * err      = Python Obit Error/message stack
    * chkPos   = If true also check the coordinates on each axis
      Check is if pixels are within 0.01 of a pixel
    * factor1  = Scaling factor for in1Image
    * factor2  = Scaling factor for in2Image
    """
    ################################################################
    # Checks
    if not Image.PIsA(in1Image):
        raise TypeError,"in1Image MUST be a Python Obit Image"
    if not Image.PIsA(in2Image):
        raise TypeError,"in2Image MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Clone output from input 1
    in1Image.Clone (outImage, err)
    # Open images
    Image.POpen (in1Image, Image.READONLY, err)
    Image.POpen (in2Image, Image.READONLY, err)
    Image.POpen (outImage, Image.WRITEONLY, err)
    #  Get input descriptor to see how many planes
    in1Desc = in1Image.Desc
    in2Desc = in2Image.Desc
    # Check compatibility
    ImageDesc.PCheckCompat (in1Desc, in2Desc, chkPos=chkPos)
    inDescDict = in1Desc.Dict
    ndim  = inDescDict["naxis"]
    inNaxis = inDescDict["inaxes"]
    # Work buffer
    inImageArray = Image.PGetFArray(in1Image)
    ImageBuffer1 = FArray.PCopy(inImageArray, err)
    ImageBuffer2 = FArray.PCopy(inImageArray, err)

    # list of planes to loop over (0-rel)
    if (ndim>0) and (inNaxis[2]>0):  
        planes = range(inNaxis[2])
    else:
        planes = [0]
    
    # Loop over planes
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        # Get image planes
        Image.PGetPlane (in1Image, ImageBuffer1, doPlane, err)
        Image.PGetPlane (in2Image, ImageBuffer2, doPlane, err)

        # Scale
        FArray.PSMul(ImageBuffer1, factor1)
        FArray.PSMul(ImageBuffer2, factor2)

        # Add
        FArray.PAdd(ImageBuffer1, ImageBuffer2, ImageBuffer2)

        # Write output
        Image.PPutPlane (outImage, ImageBuffer2, doPlane, err)

        # end loop over planes
    # Close
    in2Image.Close(err)
    in2Image.Close(err)
    outImage.Close(err)
    # Error?
    if err.isErr:
        OErr.printErrMsg(err, "Error subtracting Images")
    # Write history
    in1History  = History.History("history", in1Image.List, err)
    in2History  = History.History("history", in2Image.List, err)
    outHistory  = History.History("history", outImage.List, err)
    # Copy Histories
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit PImageAdd",err)
    outHistory.WriteRec(-1, "/ PImageAdd Input 1 History",err)
    outHistory.Close(err)
    info = in1Image.List.Dict
    # FITS? - copy header
    if ("FileType" in info) and (info["FileType"][2][0]==0):
        History.PCopyHeader(in1History, outHistory, err)
    #Not needed History.PCopy(in1History, outHistory, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.WriteRec(-1, "/      ",err)
    outHistory.WriteRec(-1, "/ ******   PImageAdd Input 2 History",err)
    outHistory.Close(err)
    info = in2Image.List.Dict
    # FITS? - copy header
    if ("FileType" in info) and (info["FileType"][2][0]==0):
        History.PCopyHeader(in2History, outHistory, err)
    History.PCopy(in2History, outHistory, err)
    # Add this programs history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit PImageAdd",err)
    outHistory.WriteRec(-1,OSystem.PGetPgmName()+" factor1 = "+str(factor1),err)
    outHistory.WriteRec(-1,OSystem.PGetPgmName()+" factor2 = "+str(factor2),err)
    outHistory.Close(err)
# end PImageAdd

def FFTHeaderUpdate(inIm, naxis, err):
    """
    Fix Image header for an image being FFTed
    
    Update first two axes for the effect of FFT

    * inID  = image with descriptor to update
    * naxis = dimensionality of array being FFTed (not size in inID)
    * err   = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inIm):
        raise TypeError,"inIm MUST be a Python Obit Image"
    header = inIm.Desc.Dict
    # Image to uv plane
    if header["ctype"][0][0:8]=="RA---SIN":
        header["ctype"][0] = "UU-L"+header["ctype"][0][4:]
        header["ctype"][1] = "VV-L"+header["ctype"][0][4:]
        dx = header["cdelt"][0]/57.296
        header["cdelt"][0] = 1.0 / (naxis[0]*dx)
        dy = header["cdelt"][1]/57.296
        header["cdelt"][1] = 1.0 / (naxis[1]*dy)
        header["crpix"][0] = 1.0 + header["inaxes"][0]/2.0
        header["crpix"][1] = 1.0 + header["inaxes"][1]/2.0
        header["crval"][0] = 0.0
        header["crval"][1] = 0.0
    # end image to uv plane
    inIm.Desc.Dict = header
    inIm.UpdateDesc(err)

    #  end FFTHeaderUpdate

import FFT, CArray, FeatherUtil
def PImageFFT (inImage, outAImage, outPImage, err):
    """
    FFTs an Image
    
    FFT inImage and write as real and imaginary as full plane (hermetian) 

    * inImage   = input Obit Python Image 1
      Any BLC and/or TRC set will be honored
    * outAImage = output Obit Python Amplitude image of FFT
      must be defined but not instantiated
    * outPImage = output Obit Python Phase (deg) image of FFT
      must be defined but not instantiated
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(outAImage):
        raise TypeError,"outAImage MUST be a Python Obit Image"
    if not Image.PIsA(outPImage):
        raise TypeError,"outPImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Clone output images
    inImage.Clone(outAImage,err)
    inImage.Clone(outPImage,err)
    OErr.printErrMsg(err, "Error initializing images")

    # Size of FFT
    inImage.Open(Image.READONLY, err)
    inImage.Read(err)
    OErr.printErrMsg(err, "Error reading input")
    inHead = inImage.Desc.Dict
    FFTdim = [FFT.PSuggestSize(inHead["inaxes"][0]), FFT.PSuggestSize(inHead["inaxes"][1])]

    # Create float arrays for FFT size
    inFArray  = FArray.FArray("inF",  naxis=FFTdim)
    outFArray = FArray.FArray("outF", naxis=FFTdim)

    # Pad input into work FArray
    FArray.PPad(inImage.FArray, inFArray, 1.0)
    # and God said "The center of an FFT will be at the corners"
    FArray.PCenter2D(inFArray)
    # Zero output FArray and use as imaginary part
    FArray.PFill(outFArray, 0.0)
    
    # Create FFT for full complex FFT
    FFTfor = FFT.FFT("FFT", 1, 1, 2, FFTdim)
    
    # Create complex arrays for FFT size
    inCArray  = CArray.CArray("inC", naxis=FFTdim)
    outCArray = CArray.CArray("outC", naxis=FFTdim)
    
    # Copy input to scratch CArray
    CArray.PComplex(inFArray, outFArray, inCArray)
    
    # FFT
    FFT.PC2C(FFTfor, inCArray, outCArray)
    
    # Extract amplitude
    CArray.PAmp(outCArray, outFArray)
    # and God said "The center of an FFT will be at the corners"
    FArray.PCenter2D(outFArray)
    
    # Extract output portion and write
    outAImage.Open(Image.WRITEONLY,err)
    outAImage.FArray = FeatherUtil.PExtract (FFTfor, outFArray, outAImage.FArray, err)
    OErr.printErrMsg(err, "Error extracting output amplitude image")
    outAImage.WriteFA(outAImage.FArray, err)
    # Fix header
    FFTHeaderUpdate(outAImage, FFTdim, err)
    outAImage.Close(err)
    OErr.printErrMsg(err, "Error writing output amplitude image")
    
    # Extract phase
    CArray.PPhase(outCArray, outFArray)
    # To degrees
    FArray.PSMul(outFArray, 57.2956)
    # and God said "The center of an FFT will be at the corners"
    FArray.PCenter2D(outFArray)

    # Extract output portion and write
    outPImage.Open(Image.WRITEONLY,err)
    outPImage.FArray = FeatherUtil.PExtract (FFTfor, outFArray, outPImage.FArray, err)
    OErr.printErrMsg(err, "Error extracting output phase image")
    outPImage.WriteFA(outPImage.FArray, err)
    # Fix header
    FFTHeaderUpdate(outPImage, FFTdim, err)
    outPImage.Close(err)
    # Error?
    OErr.printErrMsg(err, "Error writing output phase image")

    # get any BLC, TRC for history
    info = inImage.List.Dict
    blc = [1,1,1,1,1,1,1]
    if 'BLC' in info:
        blc = info["BLC"][2]
    trc = [0,0,0,0,0,0,0]
    if 'TRC' in info:
        trc = info["TRC"][2]

    # Write history
    i = 0
    imtype = ("Amplitude","Phase")
    for outImage in (outAImage, outPImage):
        inHistory  = History.History("history", inImage.List, err)
        outHistory = History.History("history", outImage.List, err)
        # Copy History
        # FITS? - copy header
        if ("FileType" in info) and (info["FileType"][2][0]==0):
            History.PCopyHeader(inHistory, outHistory, err)
        #Not needed History.PCopy(inHistory, outHistory, err)
        # Add this programs history
        outHistory.Open(History.READWRITE, err)
        outHistory.TimeStamp(" Start Obit PImageFFT",err)
        outHistory.WriteRec(-1,OSystem.PGetPgmName()+" BLC = "+str(blc),err)
        outHistory.WriteRec(-1,OSystem.PGetPgmName()+" TRC = "+str(trc),err)
        outHistory.WriteRec(-1,OSystem.PGetPgmName()+" type = "+imtype[i],err)
        i += 1
        outHistory.Close(err)
# end PImageFFT

def PImageT2Spec (inImage, outImage, nTerm, 
                  inCCVer, outCCVer, err,
                  refFreq=1.0e9, terms=None, startCC=1, endCC=0):
    """
    Convert an ObitImage(MF) (TSpec CCs) to an ObitImageWB (Spec CCs)
    
    Output CC Table will have fitted spectra rather than tabulated spectra.
    If an integrated spectrum is given, the sum of the input CC sprectra
    are forced to this spectrum.
    Copies spectral planes and converts specified CC table
    Tabulated spectrum fitted with spectrum weighting by primary beam

    * inImage  = input Obit Python Image 1
      Must have freq axis type = "SPECLNMF"
    * outImage = output Obit Python image
      must be defined but not instantiated
      On return will be replaced bu image created
    * nTerm    = Number of output Spectral terms, 2=SI, 3=also curve.
    * inCCVer  = Input CCTable to convert, 0=> highest
    * outCCVer = Output CCTable, 0=>1
    * err      = Python Obit Error/message stack
    * refFreq  = Reference frequency (Hz) for total spectrum
    * terms    = if not None, parameters of total spectrum
      [flux density at refFreq, spectral index at refFreq, ...]
    * startCC  = First 1-rel component to convert
    * endCC    = Last 1-rel component to convert, 0=> all
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    # Merge CCs to temp cc table
    tmpCCver = Image.PGetHighVer(inImage, "AIPS CC") + 1;
    inTab    = inImage.NewTable(Image.READONLY, "AIPS CC", inCCVer, err)
    noParms  = inTab.Desc.List.Dict["NO_PARMS"][2][0]
    tmpTab   = inImage.NewTable(Image.WRITEONLY, "AIPS CC", tmpCCver, err, noParms=noParms)
    TableUtil.PCCMerge(inTab, tmpTab, err)
    # Fix spectrum if needed
    if terms:
        nterm = len(terms)
        Obit.TableCCUtilFixTSpec(inImage.me, tmpCCver, \
                                 refFreq, nterm, terms,
                                 startCC, endCC, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Adjusting spectrum of CC Table")
    # Convert
    outImage.me = Obit.ImageUtilT2Spec(inImage.me, outImage.me, nTerm, tmpCCver, \
                                       outCCVer, startCC, endCC, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error Converting image/CC Table")
    # Delete temporary CC table
    inImage.ZapTable("AIPS CC", tmpCCver, err)
    # Do history for spectrum modification
    pgmName = OSystem.PGetPgmName()
    outHistory = History.History("history", outImage.List, err)
    History.POpen(outHistory, History.READWRITE, err)
    History.PTimeStamp(outHistory," Start Obit "+pgmName,err)
    if terms:
        History.PWriteRec(outHistory,-1,pgmName+"  nterm  = "+str(nterm),err)
        History.PWriteRec(outHistory,-1,pgmName+"  refFreq = "+str(refFreq ),err)
        History.PWriteRec(outHistory,-1,pgmName+"  terms   = "+str(terms),err)
    History.PClose(outHistory, err)
# end PImageT2Spec
