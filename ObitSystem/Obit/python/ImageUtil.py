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

# Python interface to ObitImageUtil utilities
import Obit, Image, ImageDesc, FArray, UV, Table, History, OErr

def PICreateImage (inUV, fieldNo, doBeam, err):
    """ Create an image from information on an ObitUV

    returns  Python Image
    inUV     = Python UV object, following read from InfoList member 
       "nChAvg" OBIT_int (1,1,1) number of channels to average.
           This is for spectral line observations and is ignored
           if the IF axis on the uv data has more than one IF.
           Default is continuum = average all freq/IFs. 0=> all.
       "rotate" OBIT_float (?,1,1) Desired rotation on sky (from N thru E) in deg. [0]
       "nx"     OBIT_int (?,1,1) Dimension of image in RA [no default].
           This and the following are arrays with one entry per field.
       "nxBeam" OBIT_int (?,1,1) Dimension of beam in RA, [def. nx]
       "ny"     OBIT_int (?,1,1) Dimension of image in declination[no default]
       "nyBeam" OBIT_int (?,1,1) Dimension of beam in declination, [def. ny]
       "xCells" OBIT_float (?,1,1) X (=RA) cell spacing in degrees [no default]
       "yCells" OBIT_float (?,1,1) Y (=dec) cell spacing in degrees [no default]
       "xShift" OBIT_float (?,1,1) Desired shift in X (=RA) in degrees. [0]
       "yShift" OBIT_float (?,1,1) Desired shift in Y (=dec) in degrees. [0]
       "nuGrid"   OBIT_int (1,1,1) Size in pixels of weighting grid for uniform weighting
    fieldNo  = Which field (1-rel) in imaging parameter arrays.
    doBeam   = if TRUE also create beam as the myBeam member of returned image.
    err      = Python Obit Error/message stack
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
    """ Grids UV, FFTs and makes corrections for the gridding convolution.

    inUV     = Input Python uv data. Should be in form of Stokes to be imaged
               will all calibration and selection applied.
    outImage = Python Image to be written.  Must be previously instantiated.
               Beam normalization factor is written to output Beam
               infoList as SUMWTS
    channel  = Which frequency channel to image, 0->all.
    doBeam   = if TRUE also make beam.  Will make the myBeam member of outImage.
               If FALSE, and myGrid->BeamNorm 0.0 then reads SUMWTS value 
               from beam infolist
    doWeigh  = if TRUE Apply uniform weighting corrections to uvdata before imaging
    err      = Python Obit Error/message stack
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
    """ Interpolates one image onto another's grid.

    Pixels in outImage are interpolated from inImage
    inImage  = Input Python Image.
    outImage = Python Image to be written.  Must be previously instantiated.
    err      = Python Obit Error/message stack
    inPlane  = 5 element int array with 1, rel. plane number [1,1,1,1,1]
               giving location of plane to be interpolated
    outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
               giving location of plane to be written
    hwidth   = half width of interpolation kernal [1-4] default 2
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
    """ Multiply an image by the primary beam pattern of another

    Pixels in outImage are inImage multiplied by the antenna beam pattern
    from pntImage
    inImage  = Input Python Image.
    pntImage = Python Image giving pointing position (ObsRA, ObsDec)
    outImage = Python Image to be written.  Must be previously instantiated.
    err      = Python Obit Error/message stack
    inPlane  = 5 element int array with 1, rel. plane number [1,1,1,1,1]
               giving location of plane to be interpolated
    outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
               giving location of plane to be written
    antSize  = Antenna diameter assumed, default 25m
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
    """ Calculate an image with a primary beam pattern

    Make an image of the antenna primary beam pattern based on the pointing
    position in an image.
    pntImage = Python Image giving pointing position (ObsRA, ObsDec)
    outImage = Python Image to be written.  Must be previously instantiated.
    err      = Python Obit Error/message stack
    minGain  = minimum allowed gain (lower values blanked).
    outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
               giving location of plane to be written
    antSize  = Antenna diameter assumed, default 25m
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
    """ Correct (divide) an image by the primary beam pattern of another

    Pixels in outImage are inImage multiplied by the antenna beam pattern
    from pntImage
    inImage  = Input Python Image.
    pntImage = Python Image giving pointing position (ObsRA, ObsDec)
    outImage = Python Image to be written.  Must be previously instantiated.
    err      = Python Obit Error/message stack
    inPlane  = 5 element int array with 1, rel. plane number [1,1,1,1,1]
               giving location of plane to be interpolated
    outPlane = 5 element int array with 1, rel. plane number [1,1,1,1,1]
               giving location of plane to be written
    antSize  = Antenna diameter assumed, default 25m
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
    """ Scale the pixel values in an image

    Scale image, optionally by plane, scales any CC tables, writes history
    inImage   Obit Python Image
    scale     Scaling factor, if scalar, multiply all pixels,
              otherwise one value per image plane.
    err       Python Obit Error/message stack
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
    """ Scale flux densities in a CC table

    Flux densities of CC entries startComp through endComp are scales by scle
    inCCTab   = Input Python TableCC
    startComp = first (1-rel) component
    endComp   = highest component [1-rel= 0=> all
    scale     = flux densiti scaling factor
    err      = Python Obit Error/message stack
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

