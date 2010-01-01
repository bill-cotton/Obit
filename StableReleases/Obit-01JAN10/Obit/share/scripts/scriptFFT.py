# python/Obit script to FFT an image producing amplitude and phase images
# Works only on FITS images
# Arguments:
# 1) Name of input image
# 2) Name of output amplitude image
# 3) Name of output phase image
#    outputs are cloned from input
#

import sys, Obit, Image, FArray, CArray, FFT, ImageUtil, FeatherUtil, OSystem, OErr, InfoList, History

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FFT", 1, 1, 0, ["Def"], 1, ["./"], True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

def FFTHeaderUpdate(inIm, naxis, err):
    """ Fix Image header for an image being FFTed

    Update first two axes for the effect of FFT
    inID  = image with descriptor to update
    naxis = dimensionality of array being FFTed (not size in inID)
    err   = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inIm):
        raise TypeError,"inIm MUST be a Python Obit Image"
    header = inIm.Desc.Dict
    # Image to uv plane
    if header["ctype"][0]=="RA---SIN":
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
    
# Get file names
inFile = sys.argv[1]
outAFile = sys.argv[2]
outPFile = sys.argv[3]
inDisk   = 1
outADisk = 1
outPDisk = 1

# Convert files into Images
inImage  = Image.newPFImage(inFile, inFile, inDisk, True, err)
outAImage  = Image.newPFImage(outAFile, outAFile, outADisk, False, err)
inImage.Clone(outAImage,err)
outPImage  = Image.newPFImage(outPFile, outPFile, outPDisk, False, err)
inImage.Clone(outPImage,err)
OErr.printErrMsg(err, "Error initializing images")

# Size of FFT
inImage.Open(Image.READONLY, err)
inImage.Read(err)
OErr.printErrMsg(err, "Error reading input")
inHead = inImage.Desc.Dict
FFTdim = [FFT.PSuggestSize(inHead["inaxes"][0]), FFT.PSuggestSize(inHead["inaxes"][1])]

# Create float arrays for FFT size
inFArray  = FArray.FArray("inF", naxis=FFTdim)
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
print "RMS",outFArray.RMS

# Extract output portion and write
outPImage.Open(Image.WRITEONLY,err)
outPImage.FArray = FeatherUtil.PExtract (FFTfor, outFArray, outPImage.FArray, err)
OErr.printErrMsg(err, "Error extracting output phase image")
print "RMS 2",outPImage.FArray.RMS
outPImage.WriteFA(outPImage.FArray, err)
# Fix header
FFTHeaderUpdate(outPImage, FFTdim, err)
outPImage.Close(err)
OErr.printErrMsg(err, "Error writing output phase image")

# Say something
print "FFTed",inFile,"to",outAFile,"and",outPFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

