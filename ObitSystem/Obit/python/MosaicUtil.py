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

# Python utility package for Mosaicing images by weighting them together
import Image, ImageDesc, ImageUtil, FArray, InfoList, OErr

def PMakeMaster(template, size, SumWtImage, SumWt2, err):
    """ Create a pair of images to accumulation of partial products

    Create an image to contain the Sum of the input Images times the
    weights, and another for the sum of the weights squared.
    The descriptive material is from image template
    template   = Image with position etc, to be copied
    size       = output image size in pixels, e.g. [200,200]
    SumWtImage = First output image, must be defined (i.e. files named)
                 but not fully created.
    SumWt2     = Second output image, like SumWtImage
    err        = Python Obit Error/message stack

    """
    ################################################################
    # Checks
    if not Image.PIsA(template):
        print "Actually ",template.__class__
        raise TypeError,"template MUST be a Python Obit Image"
    if not Image.PIsA(SumWtImage):
        print "Actually ",SumWtImage.__class__
        raise TypeError,"SumWtImage MUST be a Python Obit Image"
    if not Image.PIsA(SumWt2):
        print "Actually ",SumWt2.__class__
        raise TypeError,"SumWt2 MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Get image info from template
    Image.POpen (template, 1, err)
    Image.PRead (template, err)
    desc = Image.PGetDesc(template)
    descDict = ImageDesc.PGetDict(desc)  # Python dict object
    Image.PClose (template, err)
    #OErr.printErrMsg(err, "Error reading input image "+Image.PGetName(template))
    #
    # Create zero filled array for data
    outArray = FArray.FArray("Initial array", size)
    #
    # Modify the descriptor for output.
    naxis = size[0:3]
    # Update reference pixel, pixel shift an integral number
    dim = descDict["inaxes"]
    pixOff = [naxis[0]/2-dim[0]/2, naxis[1]/2-dim[1]/2]
    crpix = descDict["crpix"]
    crpix[0] = crpix[0] + pixOff[0]
    crpix[1] = crpix[1] + pixOff[1]
    # Update size
    dim[0] = naxis[0];
    dim[1] = naxis[1]
    #print "debug dim",dim
    descDict["inaxes"] = dim   
    descDict["bitpix"] = -32  # output floating
    #
    # Do SumWtImage
    desc = Image.PGetDesc(SumWtImage) 
    ImageDesc.PSetDict(desc, descDict) # set output descriptor
    # Write output image
    Image.POpen(SumWtImage, 2, err)
    Image.PWriteFA(SumWtImage, outArray, err)
    Image.PClose(SumWtImage, err)
    #OErr.printErrMsg(err, "Error writing image for "+Image.PGetName(SumWtImage))
    #
    # Do SumWtImage
    desc = Image.PGetDesc(SumWt2) 
    ImageDesc.PSetDict(desc, descDict) # set output descriptor
    # Write output image
    Image.POpen(SumWt2, 2, err)
    Image.PWriteFA(SumWt2, outArray, err)
    Image.PClose(SumWt2, err)
    #OErr.printErrMsg(err, "Error writing image for "+Image.PGetName(SumWt2))
    # end PMakeMaster

def PWeightImage(inImage, factor, SumWtImage, SumWt2, err, minGain=0.1):
    """ Sum an image onto Weighting accumulators using PB corrections

    Calculate the weights for an image from the primary beam pattern
    And accumulate into the correct locations in the accumulation images.
    inImage    = Image to be accumulated
    factor     = Additional multiplication factor, normally 1.0
    SumWtImage = First output image, must be defined (i.e. files named)
                 but not fully created.
    SumWt2     = Second output image, like SumWtImage
    err        = Python Obit Error/message stack
    minGain    = minimum allowed gain (lower values blanked).
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(SumWtImage):
        print "Actually ",SumWtImage.__class__
        raise TypeError,"SumWtImage MUST be a Python Obit Image"
    if not Image.PIsA(SumWt2):
        print "Actually ",SumWt2.__class__
        raise TypeError,"SumWt2 MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Open files
    #Image.POpen(inImage, 1, err)
    Image.POpen(SumWtImage, 3, err)
    Image.POpen(SumWt2, 3, err)
    #  Get output descriptor to see how many planes
    outDesc = Image.PGetDesc(SumWtImage)
    outDescDict = ImageDesc.PGetDict(outDesc)
    outNaxis = outDescDict["inaxes"]
    print "Accumulation naxis",outNaxis
    #  Get input descriptor to see how many planes
    inDesc = Image.PGetDesc(inImage)
    inDescDict = ImageDesc.PGetDict(inDesc)
    ndim  = inDescDict["naxis"]
    inNaxis = inDescDict["inaxes"]
    #print "debug input naxis is ",inNaxis
    # Test if compatible
    if inNaxis[2] < outNaxis[2]:
        print "input has",inNaxis[2],"planes and output",outNaxis[2]
        raise RuntimeError,"input image has too few planes "
    if (ndim>0) and (inNaxis[2]>0):  # list of planes to loop over (0-rel)
        planes = range(inNaxis[2])
    else:
        planes = [0]
    #
    # Loop over planes
    WtArray = None
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        # Get image 
        Image.PGetPlane (inImage, None, doPlane, err)
        #OErr.printErrMsg(err, "Error reading image for "+Image.PGetName(inImage))
        inImageArray = Image.PGetFArray(inImage)
        #
        # Make weight image, first pass
        #if iPlane ==0:
        # DEBUG
        if WtArray == None:
            WtImage = Image.Image("WeightImage")
            Image.PCloneMem(inImage, WtImage, err)
            ImageUtil.PPBImage(inImage, WtImage, err, minGain)
            OErr.printErrMsg(err, "Error making weight image for "+Image.PGetName(inImage))
            WtArray  = WtImage.FArray   # Get array
            # DEBUG
            #tempDesc  = Image.PGetDesc(WtImage)
            #tempArray = Image.PGetFArray(WtImage);
            #Image.PFArray2FITS(tempArray, "PBImage.fits", err, 0, tempDesc)
        #
        # Make image*Wt and Wt^2 memory resident images
        ImageWt = Image.Image("ImageXwt")
        Image.PCloneMem(inImage, ImageWt, err)
        ImageWtArray = Image.PGetFArray(ImageWt)
        FArray.PMul(inImageArray, WtArray, ImageWtArray);
        #
        WtWt = Image.Image("wtXwt")
        Image.PCloneMem(inImage, WtWt, err)
        WtWtArray = Image.PGetFArray(WtWt)
        FArray.PMul(WtArray, WtArray, WtWtArray);
        #
        # Now the interpolated versions to be summed to the accumulation arrays
        InterpWtImage = Image.Image("InterpWtImage")
        Image.PClone2(inImage, SumWtImage, InterpWtImage, err)
        ImageUtil.PInterpolateImage(ImageWt, InterpWtImage, err)
        #OErr.printErrMsg(err, "Error interpolating image "+Image.PGetName(inImage))
        InterpWtWt = Image.Image("InterpWtWt")
        Image.PClone2(inImage, SumWtImage, InterpWtWt, err)
        ImageUtil.PInterpolateImage(WtWt, InterpWtWt, err)
        #OErr.printErrMsg(err, "Error interpolating wt*wt "+Image.PGetName(inImage))

        # Debug
        #if iPlane == 0:
        #    tempDesc    = Image.PGetDesc(InterpWtImage)
        #    tempArray = Image.PGetFArray(InterpWtImage);
        #    Image.PFArray2FITS(tempArray, "debugWeight.fits", err, 1, tempDesc)
        #    # end debug
        #
        # Read accumulation image plane
        Image.PGetPlane(SumWtImage, None, doPlane, err)
        Image.PGetPlane(SumWt2,  None, doPlane, err)
        #OErr.printErrMsg(err, "Error reading accumulation image ")
        #
        # Determine alignment
        inDesc = Image.PGetDesc(InterpWtImage)       # get descriptors
        inDescDict = ImageDesc.PGetDict(inDesc)
        outDesc = Image.PGetDesc(SumWtImage)
        outDescDict = ImageDesc.PGetDict(outDesc)
        naxis = inDescDict["inaxes"]                # find input center pixel in output
        pos1 = [int(naxis[0]*0.5+0.5), int(naxis[1]*0.5+0.5)]
        xpos1 = [float(pos1[0]),float(pos1[1])]
        xpos2 = ImageDesc.PCvtPixel (inDesc, xpos1, outDesc, err)
        #OErr.printErrMsg(err, "Error converting pixel locations for "+Image.PGetName(inImage))
        pos2 = [int(xpos2[0]+0.5), int(xpos2[1]+0.5)]
        #print "DEBUG CvtPixel", naxis, pos1, pos2
        #
        # Accumulate
        SumWtImageArray = Image.PGetFArray(SumWtImage)
        InterpWtArray = Image.PGetFArray(InterpWtImage)
        #print "DEBUG Mean before ",FArray.PMean(SumWtImageArray)
        FArray.PShiftAdd (SumWtImageArray, pos2, InterpWtArray,  pos1, factor, SumWtImageArray)
        #print "DEBUG Mean after ",FArray.PMean(SumWtImageArray)
        SumWt2Array = Image.PGetFArray(SumWt2)
        InterpWtWtArray = Image.PGetFArray(InterpWtWt)

        # Blank weight whereever image is blank or zero
        FArray.PInClip(InterpWtArray, -1.0e-20, 1.0e-20, FArray.PGetBlank())
        FArray.PBlank (InterpWtWtArray, InterpWtArray, InterpWtWtArray);
        # DEBUG
        #print "DEBUG, RMS",InterpWtArray.RMS
        #tempDesc  = Image.PGetDesc(InterpWtWt)
        #Image.PFArray2FITS(InterpWtWtArray, "PBImage2.fits", err, 0, tempDesc)

        FArray.PShiftAdd (SumWt2Array,     pos2, InterpWtWtArray,pos1, factor, SumWt2Array)
        #
        # Debug
        #if iPlane==0:
        #    tempDesc = Image.PGetDesc(SumWtImage)
        #    tempArray = Image.PGetFArray(SumWtImage);
        #    pos=[0,0]
        #    #print "Max",FArray.PMax(tempArray,pos)
        #    print "Mean",FArray.PMean(tempArray)
        #    #Image.PFArray2FITS(tempArray, "debugFITS.fits", err, 1, tempDesc)
        #    Image.PFArray2FITS(tempArray, "debugFITS.fits", err)
        #    # end debug
            
        # Write output
        Image.PPutPlane(SumWtImage, None, doPlane, err)
        Image.PPutPlane(SumWt2, None, doPlane, err)
        #OErr.printErrMsg(err, "Error writing accumulation image ")
        # end loop over planes
    # close output
    #Image.PClose(inImage, err)
    Image.PClose(SumWtImage, err)
    Image.PClose(SumWt2, err)

    # end PWeightImage
    
def PAccumIxWt(im, wt, factor, accum, accumwt, err):
    """ Accumulate im * wt into accum

    Used to accumulate images which don't need PB corrections
    and have a weight image.
    im      = image to accumulate
    wt      = weight image corresponding to accum
    factor  = Additional multiplication factor, normally 1.0
    accum   = image into which to accumulate im*wt
    accumwt = image into which to accumulate wt
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(im):
        print "Actually ",im.__class__
        raise TypeError,"im MUST be a Python Obit Image"
    if not Image.PIsA(wt):
        print "Actually ",wt.__class__
        raise TypeError,"wt MUST be a Python Obit Image"
    if not Image.PIsA(accum):
        print "Actually ",accum.__class__
        raise TypeError,"accum MUST be a Python Obit Image"
    
    #
    # Open files
    #Image.POpen(im, 1, err)
    Image.POpen(accum, Image.READWRITE, err)
    Image.POpen(accumwt, Image.READWRITE, err)
    #  Get output descriptor to see how many planes
    outDesc     = accum.Desc
    outDescDict = outDesc.Dict
    outNaxis    = outDescDict["inaxes"]
    print "Accumulation naxis",outNaxis
    #  Get input descriptor to see how many planes
    inDesc     = im.Desc
    inDescDict = inDesc.Dict
    ndim       = inDescDict["naxis"]
    inNaxis    = inDescDict["inaxes"]
    #print "debug input naxis is ",inNaxis
    # Test if compatible
    if inNaxis[2] < outNaxis[2]:
        print "input has",inNaxis[2],"planes and output",outNaxis[2]
        raise RuntimeError,"input image has too few planes "
    if (ndim>0) and (inNaxis[2]>0):  # list of planes to loop over (0-rel)
        planes = range(inNaxis[2])
    else:
        planes = [0]
    #
    # Loop over planes
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        # Get image 
        Image.PGetPlane (im, None, doPlane, err)
        #OErr.printErrMsg(err, "Error reading image for "+Image.PGetName(im))
        imArray = im.FArray
        # Get Weight
        Image.PGetPlane (wt, None, doPlane, err)
        #OErr.printErrMsg(err, "Error reading image for "+Image.PGetName(wt))
        WtArray = wt.FArray
        #
        # Make image*Wt memory resident image
        ImageWt = Image.Image("ImageXwt")
        Image.PCloneMem(im, ImageWt, err)
        ImageWtArray = ImageWt.FArray
        FArray.PMul(imArray, WtArray, ImageWtArray);

        #
        # Now the interpolated versions to be summed to the accumulation arrays
        InterpWtImage = Image.Image("InterpWtImage")
        Image.PClone2(im, accum, InterpWtImage, err)
        ImageUtil.PInterpolateImage(ImageWt, InterpWtImage, err)
        #OErr.printErrMsg(err, "Error interpolating image "+Image.PGetName(im))
        InterpWt = Image.Image("InterpWt")
        Image.PClone2(im, accum, InterpWt, err)
        ImageUtil.PInterpolateImage(wt, InterpWt, err)
        #OErr.printErrMsg(err, "Error interpolating wt "+Image.PGetName(im))
        
        #
        # Read accumulation image plane
        Image.PGetPlane(accum, None, doPlane, err)
        Image.PGetPlane(accumwt,  None, doPlane, err)
        #OErr.printErrMsg(err, "Error reading accumulation image ")
        #
        # Determine alignment
        inDesc      = InterpWtImage.Desc
        inDescDict  = inDesc.Dict
        outDesc     = accum.Desc
        outDescDict = outDesc.Dict
        naxis       = inDescDict["inaxes"]    # find input center pixel in output
        pos1        = [int(naxis[0]*0.5+0.5), int(naxis[1]*0.5+0.5)]
        xpos1       = [float(pos1[0]),float(pos1[1])]
        xpos2       = ImageDesc.PCvtPixel (inDesc, xpos1, outDesc, err)
        #OErr.printErrMsg(err, "Error converting pixel locations for "+Image.PGetName(im))
        pos2        = [int(xpos2[0]+0.5), int(xpos2[1]+0.5)]
        #
        # Accumulate
        accumArray = accum.FArray
        InterpWtArray = InterpWtImage.FArray
        FArray.PShiftAdd (accumArray, pos2, InterpWtArray,  pos1, factor, accumArray)
        accumwtArray = accumwt.FArray
        InterpWtWtArray = InterpWt.FArray
        # Blank weight whereever image is blank or zero
        FArray.PInClip(InterpWtArray, -1.0e-20, 1.0e-20, FArray.PGetBlank())
        FArray.PBlank (InterpWtWtArray, InterpWtArray, InterpWtWtArray);
        FArray.PShiftAdd (accumwtArray,     pos2, InterpWtWtArray,pos1, factor, accumwtArray)
        #
            
        # Write output
        Image.PPutPlane(accum, None, doPlane, err)
        Image.PPutPlane(accumwt, None, doPlane, err)
        #OErr.printErrMsg(err, "Error writing accumulation image ")
        # end loop over planes
    # close output
    #Image.PClose(im, err)
    Image.PClose(accum, err)
    Image.PClose(accumwt, err)

# End PAccumIxWt

def PNormalizeImage(SumWtImage, SumWt2, outImage, err, minWt=0.1):
    """ Sum an image onto Weighting accumulators

    Normalize SumWtImage by SumWt2 write to outImage
    Minimum allowed value in SumWt2 is minWt
    SumWtImage = First output image, must be defined (i.e. files named)
                 but not fully created.
    SumWt2     = Second output image, like SumWtImage
    outImage   = Output image, must be defined.
    err        = Python Obit Error/message stack
    minWt      = minimum summed weight (lower values blanked).
    """
    ################################################################
    # Checks
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not Image.PIsA(SumWtImage):
        print "Actually ",SumWtImage.__class__
        raise TypeError,"SumWtImage MUST be a Python Obit Image"
    if not Image.PIsA(SumWt2):
        print "Actually ",SumWt2.__class__
        raise TypeError,"SumWt2 MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Open files
    Image.POpen(outImage,   2, err)
    Image.POpen(SumWtImage, 1, err)
    Image.POpen(SumWt2,     1, err)
    #  Get descriptor to see how many planes
    outDesc = Image.PGetDesc(SumWtImage)
    outDescDict = ImageDesc.PGetDict(outDesc)
    outNaxis = outDescDict["inaxes"]
    print "Accumulation naxis",outNaxis
    #  Get input descriptor to see how many planes
    inDesc     = Image.PGetDesc(outImage)
    inDescDict = ImageDesc.PGetDict(outDesc)
    ndim       = inDescDict["naxis"]
    inNaxis    = inDescDict["inaxes"]
    #print "debug input naxis is ",inNaxis
    # Test if compatible
    if inNaxis[2] < outNaxis[2]:
        print "input has",inNaxis[2],"planes and output",outNaxis[2]
        raise RuntimeError,"input image has too few planes "
    if (ndim>0) and (inNaxis[2]>0):  # list of planes to loop over (0-rel)
        planes = range(inNaxis[2])
    else:
        planes = [0]
    #
    # Loop over planes
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        # Get images
        Image.PGetPlane (SumWtImage, None, doPlane, err)
        Image.PGetPlane (SumWt2,     None, doPlane, err)
        OErr.printErrMsg(err, "Error reading images")
        # Clip
        FArray.PClipBlank (SumWt2.FArray, minWt, 1.0e25)
        # Divide
        FArray.PDiv (SumWtImage.FArray, SumWt2.FArray, outImage.FArray)
        # Write
        Image.PPutPlane(outImage, None, doPlane, err)
        OErr.printErrMsg(err, "Error Writing normalized image ")
        # end loop over planes
        # close output
    Image.PClose(outImage, err)
    Image.PClose(SumWtImage, err)
    Image.PClose(SumWt2, err)

    # end PNormalizeImage
