""" Utility module for Feathering images
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2008
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

# Python utility package for feathering images
import Image, ImageDesc, ImageUtil, FFT, FArray, CArray, OErr

def PCreateFFT(image, dir):
    """ Create an FFT object suitable for FFTing an image

    Create an FFT for an efficient size equal to or larger than image
    One needed for each direction to be FFTed.
    returns  Python FFT object
    image    = Image to be FFTed
    dir      = FFT direction 1 -> R2C, 2 -> C2R
    """
    ################################################################
    # Checks
    if not Image.PIsA(image):
        raise TypeError,"image MUST be a Python Obit Image"
    #
    # Get image info from descriptor
    desc = image.Desc
    descDict = desc.Dict
    rank = 2    # Only 2D FFTs
    dim  = descDict["inaxes"]
    # Compute size to be efficient
    i = 0
    effDim = []
    for x in dim:
        effDim.append(FFT.PSuggestSize(x))
        i = i+1
    effDim
    name = "FFT for " + Image.PGetName(image)
    out    = FFT.FFT(name, dir, 2, rank, effDim)
    return out
    # end PCreateFFT

def PCreateFFTArray(inFFT):
    """ Create a half plane CArray suitable for the output of FFTing an image

    returns  Python CArray object of suitable size (2D)
    inFFT  = FFT to be applied
    """
    ################################################################
    # Checks
    if not FFT.PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    #
    # FFT info
    FFTrank = FFT.PGetRank(inFFT)
    FFTdim  = FFT.PGetDim(inFFT)

    naxis = FFTdim[0:2]
    naxis[0] = 1 + naxis[0]/2
    out = CArray.CArray ("Temp CArray for FFT", naxis)
    return out
    # end PCreateFFTArray

def PPad (inFFT, inImage, outImage, err):
    """ Zero Pads an image as needed for an FFT

    Any blanked values are replaced with zeroes
    inFFT    = Gives size of FFT needed
    inImage  = Python Obit Image to be padded.
    outImage = Python Obit Image for output
               Must previously exist
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FFT.PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Read input, Copy Descriptor
    inImage.Open(Image.READONLY, err)
    inImage.Read(err)
    desc = inImage.Desc
    inImage.Close(err)
    #OErr.printErrMsg(err, "Error reading input image "+Image.PGetName(inImage))
    
    # FFT info
    FFTrank = FFT.PGetRank(inFFT)
    FFTdim  = FFT.PGetDim(inFFT)

    # Get data arrays
    inArray = inImage.FArray
    naxis = FFTdim[0:2]   # Make output big enough for FFT
    outArray = FArray.FArray("Padded array", naxis)

    # Copy/pad/deblank array
    # debugPPadArray (inFFT, inArray, outArray)
    FArray.PPad(inArray, outArray, 1.0)

    # Reset output image size of first two axes
    naxis = outArray.Naxis                    # output size
    descDict = desc.Dict                      # Input Python dict object
    dim   = descDict["inaxes"]                # input size
    crpix = descDict["crpix"]
    # Update reference pixel, pixel shift an integral number
    pixOff = [naxis[0]/2-dim[0]/2, naxis[1]/2-dim[1]/2]
    crpix[0] = crpix[0] + pixOff[0]
    crpix[1] = crpix[1] + pixOff[1]
    # Update size
    dim[0] = naxis[0];
    dim[1] = naxis[1]
    descDict["inaxes"] = dim   
    descDict["bitpix"] = -32  # output floating
    desc = outImage.Desc    # Replace descriptor on output
    desc.Dict = descDict

    # Write output image
    outImage.Open(Image.WRITEONLY, err)
    outImage.WriteFA(outArray, err)
    outImage.Close(err)
    outImage.FArray = outArray  # Now attach array to image to
    # keep it from being zeroed by write
    #OErr.printErrMsg(err, "Error writing padded image for "+Image.PGetName(inImage))

    # end PPad

def PBigger (naxis, inImage, outImage, err):
    """ Increases the size of an image and zero pads

    Any blanked values are replaced with zeroes
    naxis    = dimension array of desired output
    inImage  = Python Obit Image to be padded.
    outImage = Python Obit Image for output
               Must previously exist
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Read input, Copy Descriptor
    inImage.Open(Image.READONLY, err)
    inImage.Read(err)
    desc = inImage.Desc
    inImage.Close(err)
    #OErr.printErrMsg(err, "Error reading input image "+Image.PGetName(inImage))
    
    # Get data arrays
    inArray = inImage.FArray
    outArray = FArray.FArray("Padded array", naxis)

    # Copy/pad/deblank array
    FArray.PPad(inArray, outArray, 1.0)

    # Reset output image size of first two axes
    naxis = outArray.Naxis                    # output size
    descDict = desc.Dict                      # Input Python dict object
    dim   = descDict["inaxes"]                # input size
    crpix = descDict["crpix"]
    # Update reference pixel, pixel shift an integral number
    pixOff = [naxis[0]/2-dim[0]/2, naxis[1]/2-dim[1]/2]
    crpix[0] = crpix[0] + pixOff[0]
    crpix[1] = crpix[1] + pixOff[1]
    # Update size
    dim[0] = naxis[0];
    dim[1] = naxis[1]
    descDict["inaxes"] = dim   
    descDict["bitpix"] = -32  # output floating
    desc = outImage.Desc    # Replace descriptor on output
    desc.Dict = descDict

    # Write output image
    outImage.Open(Image.WRITEONLY, err)
    outImage.WriteFA(outArray, err)
    outImage.Close(err)
    outImage.FArray = outArray  # Now attach array to image to
    # keep it from being zeroed by write
    #OErr.printErrMsg(err, "Error writing padded image for "+Image.PGetName(inImage))

    # end PBigger

def PPadArray (inFFT, inArray, outArray):
    """ Zero Pads an array as needed for an FFT

    Any blanked values are replaced with zeroes
    inFFT    = Gives size of FFT needed
    inArray  = Python FArray to be padded.
    outArray = Python FArray containing inArray but zero filled.
               Must previously exist.
    """
    ################################################################
    # Checks
    if not FFT.PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    if not FArray.PIsA(inArray):
        print "Actually ",inArray.__class__
        raise TypeError,"inArray MUST be a Python Obit FArray"
    if not FArray.PIsA(outArray):
        print "Actually ",outArray.__class__
        raise TypeError,"outArray MUST be a Python Obit FArray"
    #
    # FFT info
    FFTrank = FFT.PGetRank(inFFT)
    FFTdim  = FFT.PGetDim(inFFT)

    # Array info
    ArrayNdim  = inArray.Ndim
    ArrayNaxis = inArray.Naxis

    # Zero fill output
    FArray.PFill(outArray, 0.0)

    # Insert inArray into outArray - center as well as possible
    pos1 = [FFTdim[0]/2, FFTdim[1]/2]
    pos2 = [ArrayNaxis[0]/2, ArrayNaxis[1]/2]
    FArray.PShiftAdd (outArray, pos1, inArray, pos2, 1.0, outArray)

    # Replace any blanks with zeroes
    FArray.PDeblank(outArray, 0.0)
    # end PPadArray

def PExtract (inFFT, inArray, outArray, err):
    """ Extract a Real array from one padded for FFTs

    Any blanked values are replaces with zeroes
    returns outArray
    inFFT    = Gives size of FFT used
    inArray  = Python FArray with FFT results.
    outArray = Python FArray describing results
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FFT.PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    if not FArray.PIsA(inArray):
        print "Actually ",inArray.__class__
        raise TypeError,"inArray MUST be a Python Obit FArray"
    if not FArray.PIsA(outArray):
        print "Actually ",outArray.__class__
        raise TypeError,"outArray MUST be a Python Obit FArray"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # FFT info
    FFTrank = FFT.PGetRank(inFFT)
    FFTdim  = FFT.PGetDim(inFFT)
    
    # Target Array info
    ArrayNdim  = outArray.Ndim
    ArrayNaxis = outArray.Naxis

    # Get window to extract
    cen = [FFTdim[0]/2, FFTdim[1]/2];
    blc = [0,0]; trc=[0,0]
    blc[0] = cen[0] - ArrayNaxis[0] / 2; trc[0] = cen[0] - 1 + ArrayNaxis[0] / 2
    blc[1] = cen[1] - ArrayNaxis[1] / 2; trc[1] = cen[1] - 1 + ArrayNaxis[1] / 2

    # Extract
    out = FArray.PSubArr(inArray, blc, trc, err)
    return out
    # end PExtract

def PDeblank (inImage, value):
    """ Replace blanks in the FArray for image inImage

    Any blanked values are replaced with value
    inImage  = Python Image whose FArray is to be deblanked
    value    = value to replace blanks, e.g. 0.0
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    #
    # Get data array
    inArray  = Image.PGetFArray(inImage)
    # Replace any blanks with value
    FArray.PDeblank(inArray, value)
    # end PDeblank

def PMakeBeamMask (inImage, inFFT, err):
    """ Make uv plane weighting array

    Creates an FArray the size of a plane in inImage, FFT,
    takes real part and normalizes the central value to one
    Resulting array is returned.
    inImage  = Python Image whose FArray is to be converted to a weight mask
    inFFT    = Python Obit fortward FFT object
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not FFT.PIsA(inFFT):
        print "Actually ",inFFT.__class__
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    #
    # Make copy of data array
    inArray  = inImage.FArray
    outArray = FArray.PClone(inArray, err)
    #OErr.printErrMsg(err, "Error duplicating FArray for "+Image.PGetName(inImage))

    # Add model
    PCreateModel(inImage, outArray)
    
    # Pad for FFT
    FFTdim  = FFT.PGetDim(inFFT)
    FFTArray = FArray.PClone(inArray, err)
    naxis = FFTdim[0:2]   # Make output big enough for FFT
    FFTArray = FArray.FArray("FFT array", naxis)
    PPadArray (inFFT, outArray, FFTArray)
    del outArray   # Cleanup
   
    # Swaparoonie
    FArray.PCenter2D (FFTArray)

    # FFT
    uvArray = PCreateFFTArray(inFFT)
    PFFTR2C (inFFT, FFTArray, uvArray)
    del FFTArray  # Cleanup

    # Extract Real part
    naxis = CArray.PGetNaxis(uvArray)[0:2]
    maskArray = FArray.FArray("Mask array for "+Image.PGetName(inImage), naxis)
    CArray.PReal(uvArray, maskArray)
    del uvArray   # Cleanup

    # Normalize
    pos = [0, 1+naxis[1]/2]
    peak = FArray.PMax(maskArray, pos)
    norm = 1.0 / peak
    FArray.PSMul(maskArray, norm)

    return maskArray
    # end PMakeBeamMask

def PFFTR2C (inFFT, inArray, outArray):
    """ Real to half plane complex FFT

    inFFT    = Python Obit FFT object
    inArray  = Python FArray To be FFTed
    outArray = Python CArray to contain the FFT
               Must previously exist
    """
    ################################################################
    if not FFT.PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    if not FArray.PIsA(inArray):
        print "Actually ",inArray.__class__
        raise TypeError,"inArray MUST be a Python Obit FArray"
    if not CArray.PIsA(outArray):
        print "Actually ",outArray.__class__
        raise TypeError,"outArray MUST be a Python Obit CArray"
    #
    FFT.PR2C(inFFT, inArray, outArray)
    # end PFFTR2C

def PFFTC2R (inFFT, inArray, outArray):
    """  Half plane complex to Real FFT

    inFFT    = Python Obit FFT object
    inArray  = Python CArray To be FFTed
    outArray = Python FArray to contain the FFT
               Must previously exist
    """
    ################################################################
    if not FFT.PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    if not CArray.PIsA(inArray):
        print "Actually ",inArray.__class__
        raise TypeError,"inArray MUST be a Python Obit CArray"
    if not FArray.PIsA(outArray):
        print "Actually ",outArray.__class__
        raise TypeError,"outArray MUST be a Python Obit FArray"
    #
    FFT.PC2R(inFFT, inArray, outArray)
    # end PFFTC2R

def PCreateModel(image, outArray):
    """ Fill an FArray with a model the size and shape of the resolution in an image

    A model image is inserted in outArray derived from the restoring beam in
    image.  The size and geometry of OutArray must be those described by image
    image    = Python Obit Image with info
    outArray = Python FArray to receive model
               Must previously exist
    """
    ################################################################
    # Checks
    if not Image.PIsA(image):
        raise TypeError,"image MUST be a Python Obit Image"
    if not FArray.PIsA(outArray):
        print "Actually ",outArray.__class__
        raise TypeError,"outArray MUST be a Python Obit FArray"
    #
    # Check compatability
    array = Image.PGetFArray(image)
    if not FArray.PIsCompatable (array, outArray):
        array = Obit.FArrayUnref(array)    # Cleanup
        raise TypeError,"image and array incompatible"
    del array    # Cleanup

    # Get image info from descriptor
    desc     = image.Desc
    descDict = desc.Dict
    beamMaj = descDict["beamMaj"]
    beamMin = descDict["beamMin"]
    beamPA  = descDict["beamPA"]
    cdelt   = descDict["cdelt"]
    crpix   = descDict["crpix"]
    crota   = descDict["crota"]
    inaxes  = descDict["inaxes"]

    # Check that beam OK
    if beamMaj < 0.0001/3600.0:
        raise TypeError,"No beam provided for "+image
    print "Beam",beamMaj*3600.0,beamMin*3600.0,beamPA

    # Zero array
    FArray.PFill (outArray, 0.0)

    amp = 1.0
    Cen = [crpix[0]-1.0, crpix[1]-1.0] # zero ref
    Cen = [float(inaxes[0]/2), float(inaxes[1]/2)] # zero ref
    GauMod = [beamMaj/abs(cdelt[0]), beamMin/abs(cdelt[0]), beamPA-90.0]
                          
    FArray.PEGauss2D(outArray, amp, Cen, GauMod)
    # end PCreateModel

def PAccumImage(FFTfor, inImage, wtArray, accArray, workArray, err):
    """ Accumulate the weighted FT of an FArray

    inImage is FFTed, multiplied by wtArray and accumulated into accArray
    FFTfor   = FFT object to FT inArray
    inImage  = Image to be accumulated
               must be a size compatable with FFTfor
               returned with contents swapped for FFTs
    wtArray  = FArray containing accumulation weights
               must be a size compatable with FT of inArray
    accArray = CArray in which the results are to be accumulated
               must be a size compatable with FT of inArray
    workArray= CArray for temporary storage of FT of inArray
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FFT.PIsA(FFTfor):
        print "Actually ",FFTfor.__class__
        raise TypeError,"FFTfor MUST be a Python Obit FFT"
    if not Image.PIsA(inImage ):
        print "Actually ",inImage.__class__
        raise TypeError," inImage MUST be a Python Obit Image"
    if not FArray.PIsA(wtArray):
        print "Actually ",wtArray.__class__
        raise TypeError,"wtArray MUST be a Python Obit FArray"
    if not CArray.PIsA(accArray):
        print "Actually ",accArray.__class__
        raise TypeError,"accArray MUST be a Python Obit CArray"
    if not CArray.PIsA(workArray):
        print "Actually ",workArray.__class__
        raise TypeError,"workArray MUST be a Python Obit CArray"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Check compatability
    if not CArray.PIsCompatable (accArray, workArray):
        raise TypeError,"accArray and workArray incompatible"

    # Get array from image
    inArray = Image.PGetFArray(inImage)

    # Swaparoonie
    FArray.PCenter2D (inArray)
    
    # FFT
    PFFTR2C (FFTfor, inArray, workArray)

    # Multiply by weights
    CArray.PFMul(workArray, wtArray, workArray)

    # Scale by inverse of beam area to get units the same
    # Get image info from descriptor
    desc     = inImage.Desc
    descDict = desc.Dict
    beamMaj = 3600.0 * descDict["beamMaj"]
    beamMin = 3600.0 * descDict["beamMin"]
    cdelt   = descDict["cdelt"]
    factor = (abs(cdelt[1])/beamMaj) * (abs(cdelt[1])/beamMin)
    CArray.PSMul(workArray, factor)

    # Accumulate
    CArray.PAdd(accArray, workArray, accArray)

    # end PAccumImage
    
def PBackFFT(FFTrev, inArray, outArray, err):
    """ Back transform half plane complex to real

    inArray is FFTed (half plane complex - real) to outArray
    FFTref   = FFT object to FT inArray
    inArray  = CArray with image to be FFTed
               must be a size compatable with FFTrev
    outArray = FArray for output
               must be a size compatable with FT of inArray
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FFT.PIsA(FFTrev):
        print "Actually ",FFTrev.__class__
        raise TypeError,"FFTrev MUST be a Python Obit FFT"
    if not CArray.PIsA(inArray ):
        print "Actually ",inArray.__class__
        raise TypeError,"inArray MUST be a Python Obit CArray"
    if not FArray.PIsA(outArray ):
        print "Actually ",outArray.__class__
        raise TypeError,"outArray MUST be a Python Obit FArray"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # FFT
    PFFTC2R (FFTrev, inArray, outArray)

    FArray.PCenter2D (outArray)
  # end PBackFFT

def PInterpol(inImage, tmplImage, outImage, err):
    """ HGEOM-like operation (Before EWG got to it)

    outImage is inImage interpolated to the grid of tmplImage
    inImage  = Image to be interpolated
    tmplImage= Image whose geometry is to be used.
    outImage = values from inImage on grid of tmplImage
               undefined values set to zero to allow FFT.
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage ):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not Image.PIsA(tmplImage ):
        print "Actually ",tmplImage.__class__
        raise TypeError,"tmplImage MUST be a Python Obit Image"
    if not Image.PIsA(outImage ):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Clone, interpolate, deblank
    Image.PClone(tmplImage, outImage, err)
    #OErr.printErrMsg(err, "Error cloning padded image")

    # Copy selected header info
    descDict = inImage.Desc.Dict
    beamMaj = descDict["beamMaj"]
    beamMin = descDict["beamMin"]
    beamPA  = descDict["beamPA"]
    obsdat  = descDict["obsdat"]
    teles   = descDict["teles"]
    crval   = descDict["crval"]
    #Image.POpen(outImage, Image.READWRITE, err)
    outImage.Open(Image.READWRITE, err)
    #OErr.printErrMsg(err, "Error updating padded image")
    descDict = outImage.Desc.Dict
    descDict["beamMaj"] = beamMaj
    descDict["beamMin"] = beamMin
    descDict["beamPA"]  = beamPA
    descDict["obsdat"]  = obsdat
    descDict["teles"]   = teles
    outImage.Desc.Dict = descDict     # Update header
    Image.PDirty(outImage)            # Force header update
    outImage.Close(err)
    #OErr.printErrMsg(err, "Error updating padded image")

    # Interpolate image
    ImageUtil.PInterpolateImage(inImage, outImage, err)
    #OErr.printErrMsg(err, "Error writing padded image")

    # Deblank
    outImage.Open(Image.READWRITE, err)
    outImage.Read(err)
    PDeblank(outImage, 0.0)
    #outImage.Write(err)
    outImage.Close(err)
    #OErr.printErrMsg(err, "Error deblanking padded image")
   
    # reread image to have in FArray
    #outImage.Open(Image.READONLY, err)
    #outImage.Read(err)
    #outImage.Close(err)
    #OErr.printErrMsg(err, "Error rereading padded image")
   
    # end  PInterpol

def PSubImage(inImage, inArray, outImage, err):
    """ Extract the subimage in inAray corresponding to outImage

    This assumes that both input and output images have pixels
    on the same locations (i.e. one the padded version of the other)
    outImage is updated in permanent storage (disk)
    inImage  = Image describing inArray
    inArray  = Image to be extracted
    outImage = accepts values from inImage, must exist and be fully defined
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage ):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not FArray.PIsA(inArray ):
        print "Actually ",inArray.__class__
        raise TypeError,"inArray MUST be a Python Obit FArray"
    if not Image.PIsA(outImage ):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # get selected input header info
    descDict = inImage.Desc.Dict
    inCrpix   = descDict["crpix"]
    inInaxes  = descDict["inaxes"]
    # get selected output header info
    descDict = outImage.Desc.Dict
    outCrpix   = descDict["crpix"]
    outInaxes  = descDict["inaxes"]

    # determine window in image
    blc = [1,1,1,1,1,1,1]
    blc[0] = (int(inCrpix[0]) - int(outCrpix[0]))
    blc[1] = (int(inCrpix[1]) - int(outCrpix[1]))
    trc = inInaxes
    trc[0] = blc[0] + outInaxes[0] - 1
    trc[1] = blc[1] + outInaxes[1] - 1

    # Extract array
    outArray = FArray.PSubArr(inArray, blc, trc, err)
    #OErr.printErrMsg(err, "Error extracting sub-image")

    # rewrite image
    outImage.Open(Image.READWRITE, err)
    outImage.WriteFA(outArray, err)
    outImage.Close(err)
    #OErr.printErrMsg(err, "Error writing sub-image image")
   
    # end  PSubImage
