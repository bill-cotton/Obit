# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2014
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
import Image, ImageDesc, ImageUtil, FArray, InfoList, OErr, History
import os

def PMakeMaster(template, size, SumWtImage, SumWt2, err):
    """
    Create a pair of images to accumulation of partial products
    
    Create an image to contain the Sum of the input Images times the
    weights, and another for the sum of the weights squared.
    The descriptive material is from image template

    * template   = Image with position etc, to be copied
    * size       = output image size in pixels, e.g. [200,200]
    * SumWtImage = First output image, must be defined (i.e. files named)
      but not fully created.
    * SumWt2     = Second output image, like SumWtImage
    * err        = Python Obit Error/message stack
    
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
    desc   = Image.PGetDesc(SumWtImage)
    ImageDesc.PSetDict(desc, descDict) # set output descriptor
    # Write output image
    Image.POpen(SumWtImage, 2, err)
    nplane = template.Desc.Dict["inaxes"][2]
    for iplane in range(1,(nplane+1)):
        plane = [iplane,1,1,1,1]
        SumWtImage.PutPlane(outArray, plane, err)
    Image.PClose(SumWtImage, err)
    #OErr.printErrMsg(err, "Error writing image for "+Image.PGetName(SumWtImage))
    #
    # Do SumWtImage
    desc = Image.PGetDesc(SumWt2) 
    ImageDesc.PSetDict(desc, descDict) # set output descriptor
    # Write output image
    Image.POpen(SumWt2, 2, err)
    for iplane in range(1,(nplane+1)):
        plane = [iplane,1,1,1,1]
        SumWt2.PutPlane(outArray, plane, err)
    Image.PClose(SumWt2, err)
    #OErr.printErrMsg(err, "Error writing image for "+Image.PGetName(SumWt2))
    # Write history - sorta
    inHistory  = History.History("history", template.List, err)
    outHistory = History.History("history", SumWtImage.List, err)
    # Copy History
    History.PCopy(inHistory, outHistory, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit PMakeMaster",err)
    outHistory.Close(err)
    # end PMakeMaster

def PWeightImage(inImage, factor, SumWtImage, SumWt2, err, minGain=0.1,
                 iblc=[1,1], itrc=[0,0], restart=0):
    """
    Sum an image onto Weighting accumulators using PB corrections
    
    Calculate the weights for an image from the primary beam pattern
    And accumulate into the correct locations in the accumulation images.

    * inImage    = Image to be accumulated
    * factor     = Additional multiplication factor, normally 1.0
    * SumWtImage = First output image, must be defined (i.e. files named)
      but not fully created.
    * SumWt2     = Second output image, like SumWtImage
    * err        = Python Obit Error/message stack
    * minGain    = minimum allowed gain (lower values blanked).
    * iblc       = BLC in plane to start selection
    * itrc       = TRC in plane to end selection
    * restart   = restart channel no. 0-rel
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
    # Open accumulation files
    #Image.POpen(inImage, 1, err)
    Image.POpen(SumWtImage, 3, err)
    Image.POpen(SumWt2, 3, err)
    #  Get output descriptor to see how many planes
    outDesc     = Image.PGetDesc(SumWtImage)
    outDescDict = ImageDesc.PGetDict(outDesc)
    outNaxis    = outDescDict["inaxes"]
    print "Accumulation naxis",outNaxis
    #  Get input descriptor to see how many planes
    inDesc     = Image.PGetDesc(inImage)
    inDescDict = ImageDesc.PGetDict(inDesc)
    ndim       = inDescDict["naxis"]
    inNaxis    = inDescDict["inaxes"]
    
    # Test if compatible
    if inNaxis[2] < outNaxis[2]:
        print "input has",inNaxis[2],"planes and output",outNaxis[2]
        raise RuntimeError,"input image has too few planes "
    if (ndim>0) and (inNaxis[2]>0):  # list of planes to loop over (0-rel)
        planes = range(restart,inNaxis[2])
    else:
        planes = [0]
    #
    # Set BLC,TRC 
    inImage.List.set("BLC",[iblc[0], iblc[1],1,1,1,1,1])
    inImage.List.set("TRC",[itrc[0], itrc[1],0,0,0,0,0])
    # Loop over planes
    #print "Before loop",os.times()
    WtImage = None
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        if not (iPlane%20):
            print "At plane", iPlane+1,os.times()
        # Make weight image, first pass
        if WtImage == None:
            # Get image 
            Image.PGetPlane (inImage, None, doPlane, err)
            #OErr.printErrMsg(err, "Error reading image for "+Image.PGetName(inImage))
            #
            WtImage = Image.Image("WeightImage")
            Image.PCloneMem(inImage, WtImage, err)
            pln = [max(1,inNaxis[2]/2),1,1,1,1]
            ImageUtil.PPBImage(inImage, WtImage, err, minGain, outPlane=pln)
            OErr.printErrMsg(err, "Error making weight image for "+Image.PGetName(inImage))
            
            # The interpolated versions
            InterpWtImage = Image.Image("InterpWtImage")
            Image.PClone2(inImage, SumWtImage, InterpWtImage, err)
            
            # Interpolated weight image
            InterpWt = Image.Image("InterpWt")
            Image.PClone2(inImage, SumWtImage, InterpWt, err)
            ImageUtil.PInterpolateImage(WtImage, InterpWt, err)
            #OErr.printErrMsg(err, "Error interpolating wt*wt "+Image.PGetName(inImage))
            # Interpolated weight image Squared
            InterpWtWt = Image.Image("InterpWtWt")
            Image.PClone2(inImage, SumWtImage, InterpWtWt, err)
            # Determine alignment
            inDesc = Image.PGetDesc(InterpWtImage)       # get descriptors
            inDescDict = ImageDesc.PGetDict(inDesc)
            outDesc = Image.PGetDesc(SumWtImage)
            outDescDict = ImageDesc.PGetDict(outDesc)
            naxis = inDescDict["inaxes"]                # find input center pixel in output
            pos1 = [int(naxis[0]*0.5+0.5), int(naxis[1]*0.5+0.5)]
            xpos1 = [float(pos1[0]),float(pos1[1])]
            xpos2 = ImageDesc.PCvtPixel (inDesc, xpos1, outDesc, err)
            pos2 = [int(xpos2[0]+0.5), int(xpos2[1]+0.5)]
            # End init wt image
        # Interpolate image plane
        ImageUtil.PInterpolateImage(inImage, InterpWtImage, err, inPlane=doPlane)
        #print "after interpolate",os.times()
        #pos=[0,0]  # DEBUG
        #print "Max  I",FArray.PMax(inImage.FArray,pos),pos
        #print "Max interpolated I",FArray.PMax(InterpWtImage.FArray,pos),pos
        #inImage.Close(err)
        #OErr.printErrMsg(err, "Error interpolating image "+Image.PGetName(inImage))


        # Interpolated image times beam
        FArray.PMul(InterpWtImage.FArray, InterpWt.FArray, InterpWtImage.FArray)
        #print "after mul 1",os.times()
        #
        # Read accumulation image planes
        Image.PGetPlane(SumWt2,  None, doPlane, err)
        Image.PGetPlane(SumWtImage, None, doPlane, err)
        #print "after read accum",os.times()
        #OErr.printErrMsg(err, "Error reading accumulation image ")
        #
        #
        #pos=[0,0]  # DEBUG
        #print "Max I*Wt",FArray.PMax(InterpWtImage.FArray,pos),pos
        # Accumulate
        FArray.PShiftAdd (SumWtImage.FArray, pos2, InterpWtImage.FArray,  pos1, factor, SumWtImage.FArray)
        #print "after shift Add 1",os.times()

        # Square weight image
        FArray.PMul(InterpWt.FArray, InterpWt.FArray, InterpWtWt.FArray)
        #print "after mul 2",os.times()
        
        # Blank weight whereever image is blank or zero
        FArray.PInClip(InterpWt.FArray, -1.0e-20, 1.0e-20, FArray.PGetBlank())
        # Blank weight squared where image * Wt is blanked
        FArray.PBlank (InterpWtWt.FArray, InterpWt.FArray, InterpWtWt.FArray);
        #print "after blank",os.times()
        # Accumulate Wt*Wt
        FArray.PShiftAdd (SumWt2.FArray,     pos2, InterpWtWt.FArray,pos1, factor, SumWt2.FArray)
        #print "after shift Add 2",os.times()
        #
        # Write output
        Image.PPutPlane(SumWt2, None, doPlane, err)
        Image.PPutPlane(SumWtImage, None, doPlane, err)
        #print "after write",os.times()
        #OErr.printErrMsg(err, "Error writing accumulation image ")
        # end loop over planes
    # close output
    #Image.PClose(inImage, err)
    Image.PClose(SumWtImage, err)
    Image.PClose(SumWt2, err)

    # end PWeightImage
    
def PAccumIxWt(im, wt, factor, accum, accumwt, err):
    """
    Accumulate im * wt into accum
    
    Used to accumulate images which don't need PB corrections
    and have a weight image.

    * im      = image to accumulate
    * wt      = weight image corresponding to accum
    * factor  = Additional multiplication factor, normally 1.0
    * accum   = image into which to accumulate im*wt
    * accumwt = image into which to accumulate wt
    * err     = Python Obit Error/message stack
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
    """
    Sum an image onto Weighting accumulators
    
    Normalize SumWtImage by SumWt2 write to outImage
    Minimum allowed value in SumWt2 is minWt

    * SumWtImage = First output image, must be defined (i.e. files named)
      but not fully created.
    * SumWt2     = Second output image, like SumWtImage
    * outImage   = Output image, must be defined.
    * err        = Python Obit Error/message stack
    * minWt      = minimum summed weight (lower values blanked).
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

    # Write history - sorta
    inHistory  = History.History("history", SumWtImage.List, err)
    outHistory = History.History("history", outImage.List, err)
    # Copy History
    History.PCopy(inHistory, outHistory, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit PNormalize",err)
    outHistory.Close(err)
    # end PNormalizeImage

def PGetOverlap(in1Image, in2Image, err):
    """
    Determine the overlap region in in1Image with in2Image
    
    Returns (BLC, TRC) in in1Image of overlap, only BLC pixel if no overlap
    * in1Image   = first input image
    * in2Image   = second input image, need not be same grid
                   but should not be rotated wrt in1Image
    * err        = Python Obit Error/message stack
    
    """
    ################################################################
    # Checks
    if not Image.PIsA(in1Image):
        print "Actually ",inI1mage.__class__
        raise TypeError,"in1Image MUST be a Python Obit Image"
    if not Image.PIsA(in2Image):
        print "Actually ",in21mage.__class__
        raise TypeError,"in2Image MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Is there overlap?
    if ImageDesc.POverlap(in1Image.Desc, in2Image.Desc, err):
        d1 = in1Image.Desc.Dict
        nx1 = d1['inaxes'][0]
        ny1 = d1['inaxes'][1]
        d2 = in2Image.Desc.Dict
        nx2 = d2['inaxes'][0]
        ny2 = d2['inaxes'][1]
        # DEBUG ALL
        #return ([1,1,1,1,1,1,1], [nx1,ny1,0,0,0,0,0])
        # Determine corners of in1Image in in2Image  1-rel
        corn1 = []
        xpos = [float(0), float(0)]
        ypos = ImageDesc.PCvtPixel (in1Image.Desc, xpos, in2Image.Desc, err)
        corn1.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        xpos = [float(0), float(ny1)]
        ypos = ImageDesc.PCvtPixel (in1Image.Desc, xpos, in2Image.Desc, err)
        corn1.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        xpos = [float(nx1), float(ny1)]
        ypos = ImageDesc.PCvtPixel (in1Image.Desc, xpos, in2Image.Desc, err)
        corn1.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        xpos = [float(nx1), float(0)]
        ypos = ImageDesc.PCvtPixel (in1Image.Desc, xpos, in2Image.Desc, err)
        corn1.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        # Determine corners of in2Image in in1Image
        corn2 = []
        xpos = [float(0), float(0)]
        ypos = ImageDesc.PCvtPixel (in2Image.Desc, xpos, in1Image.Desc, err)
        corn2.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        xpos = [float(0), float(ny2)]
        ypos = ImageDesc.PCvtPixel (in2Image.Desc, xpos, in1Image.Desc, err)
        corn2.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        xpos = [float(nx1), float(ny2)]
        ypos = ImageDesc.PCvtPixel (in2Image.Desc, xpos, in1Image.Desc, err)
        corn2.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        xpos = [float(nx2), float(0)]
        ypos = ImageDesc.PCvtPixel (in2Image.Desc, xpos, in1Image.Desc, err)
        corn2.append([int(ypos[0]+0.5), int(ypos[1]+0.5)])
        # 1 entirely inside 2?
        if ((corn1[0][0]>0) and (corn1[2][0]>0) and (corn1[0][1]<=nx2) and (corn1[2][0]<=ny2)):
            return ([1,1,1,1,1,1,1], [nx1, ny1, 0,0,0,0,0])
        # Corner 0 in in2?
        if ((corn1[0][0]>0) and (corn1[0][0]<=nx2) and (corn1[0][1]>0) and (corn1[0][1]<=ny2)):
            blc = [1, 1, 1, 1, 1, 1, 1]
            trc = [min(corn2[2][0],nx1), min(corn2[2][1],ny1), 0,0,0,0,0]
            return (blc, trc)
        # Corner 1 in in2?
        if ((corn1[1][0]>0) and (corn1[1][0]<=nx2) and (corn1[1][1]>0) and (corn1[1][1]<=ny2)):
            blc = [1, min(corn2[3][1], ny1), 1, 1, 1, 1, 1]
            trc = [min (corn2[3][0], nx1), ny1, 0,0,0,0,0]
            return (blc, trc)
        # Corner 2 in in2?
        if ((corn1[2][0]>0) and (corn1[2][0]<=nx2) and (corn1[2][1]>0) and (corn1[2][1]<=ny2)):
            blc = [max(1, corn2[0][0]), max(1, corn2[0][1]), 1, 1, 1, 1, 1]
            trc = [nx1, ny1,  0,0,0,0,0]
            return (blc, trc)
         # Corner 3 in in2?
        if ((corn1[3][0]>0) and (corn1[3][0]<=nx2) and (corn1[3][1]>0) and (corn1[3][1]<=ny2)):
            blc = [max(1,corn2[1][0]), 1, 1, 1, 1, 1]
            trc = [nx1, min(corn2[1][1],ny1), 0,0,0,0,0]
            return (blc, trc)
        # What the hell, make it all
        return ([1,1,1,1,1,1,1], [nx1,ny1,0,0,0,0,0])
    else:
        # Default is no overlap
        print "no overlap"
        return ([1,1,1,1,1,1,1], [1,1,1,1,1,1,1])
    # end PGetOverlap

def PMaskCube(inImage, Mask, outImage, err):
    """
    Blank inImage where Mask is blanked or 0.0
    
    * inImage    = input Image cube
    * Mask       = 1 plane mask Image
    * outImage   = Output Image, must be defined.
    * err        = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(Mask):
        print "Actually ",Mask.__class__
        raise TypeError,"Mask MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #  Clone output
    inImage.Clone(outImage, err)
    # Open files
    Image.POpen(outImage,Image.WRITEONLY, err)
    Image.POpen(inImage, Image.READONLY,  err)
    Image.POpen(Mask,    Image.READONLY,  err)
    OErr.printErrMsg(err, "Error opening images")
    #  how many planes?
    ndim       = inImage.Desc.Dict["naxis"]
    inNaxis    = inImage.Desc.Dict["inaxes"]
    # list of planes to loop over (0-rel)
    if (ndim>2) and (inNaxis[2]>0):  
        planes = range(inNaxis[2])
    else:
        planes = [0]
    # Read Mask plane
    Image.PGetPlane (Mask, None, [1,1,1,1,1], err)
    OErr.printErrMsg(err, "Error reading mask image")
    # Mask where exactly 0.0
    FArray.PInClip(Mask.FArray, -1.0e-25, 1.0e-25, FArray.fblank)
    # Loop over planes
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        # Get image plane
        Image.PGetPlane (inImage, None, doPlane, err)
        OErr.printErrMsg(err, "Error reading input image")
        # Make sure compatable
        if not FArray.PIsCompatable(inImage.FArray, Mask.FArray):
            raise RuntimeError,"inImage and Mask incompatable"
        # Mask where blanked
        FArray.PBlank (inImage.FArray, Mask.FArray, outImage.FArray)
        # Write
        Image.PPutPlane(outImage, None, doPlane, err)
        OErr.printErrMsg(err, "Error Writing blanked image ")
        # end loop over planes
    # close files
    Image.PClose(outImage, err)
    Image.PClose(inImage, err)
    Image.PClose(Mask, err)
    # Write history
    inHistory  = History.History("history", inImage.List, err)
    outHistory = History.History("history", outImage.List, err)
    # Copy History
    History.PCopy(inHistory, outHistory, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit PMaskCube",err)
    outHistory.Close(err)
    # end PMaskCube

def PMaskCube2(inImage, Mask, outImage, err):
    """
    Blank inImage where Mask is blanked or 0.0
    
    * inImage    = input Image cube
    * Mask       = mask Image cube
    * outImage   = Output Image, must be defined.
    * err        = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(Mask):
        print "Actually ",Mask.__class__
        raise TypeError,"Mask MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #  Clone output
    inImage.Clone(outImage, err)
    # Open files
    Image.POpen(outImage,Image.WRITEONLY, err)
    Image.POpen(inImage, Image.READONLY,  err)
    Image.POpen(Mask,    Image.READONLY,  err)
    OErr.printErrMsg(err, "Error opening images")
    #  how many planes?
    ndim       = inImage.Desc.Dict["naxis"]
    inNaxis    = inImage.Desc.Dict["inaxes"]
    # list of planes to loop over (0-rel)
    if (ndim>2) and (inNaxis[2]>0):  
        planes = range(inNaxis[2])
    else:
        planes = [0]
    # Loop over planes
    for iPlane in planes:
        doPlane = [iPlane+1,1,1,1,1]
        # Read Mask plane
        Image.PGetPlane (Mask, None, doPlane, err)
        OErr.printErrMsg(err, "Error reading mask image")
        # Mask where exactly 0.0
        FArray.PInClip(Mask.FArray, -1.0e-25, 1.0e-25, FArray.fblank)
        # Get image plane
        Image.PGetPlane (inImage, None, doPlane, err)
        OErr.printErrMsg(err, "Error reading input image")
        # Make sure compatable
        if not FArray.PIsCompatable(inImage.FArray, Mask.FArray):
            raise RuntimeError,"inImage and Mask incompatable"
        # Mask where blanked
        FArray.PBlank (inImage.FArray, Mask.FArray, outImage.FArray)
        # Write
        Image.PPutPlane(outImage, None, doPlane, err)
        OErr.printErrMsg(err, "Error Writing blanked image ")
        # end loop over planes
    # close files
    Image.PClose(outImage, err)
    Image.PClose(inImage, err)
    Image.PClose(Mask, err)
    # Write history
    inHistory  = History.History("history", inImage.List, err)
    outHistory = History.History("history", outImage.List, err)
    # Copy History
    History.PCopy(inHistory, outHistory, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit PMaskCube",err)
    outHistory.Close(err)
    # end PMaskCube2


