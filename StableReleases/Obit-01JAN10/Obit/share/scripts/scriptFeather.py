# python/Obit script to feather together a series of images
# Arguments:
# 1) Name of output combined image
# 2-n) Names of the images to be combined, in order of DECREASING resolution
#      output is cloned from first named
#
# Python scriptFeather.py FeatherOut.fits  ComaA-P-C.fits ComaA-P-D.fits ComaAPGBT1950.fits
# Python scriptFeather.py FeatherOut.fits  ComaA-P-Bsub.fits ComaA-P-C.fits ComaA-P-D.fits ComaAPGBT1950.fits
# Python scriptFeather.py GCFeather.fits Gal.Cen.30s.Mos30.fits GCLBand.fits

import sys, Obit, Image, FArray, CArray, FFT, ImageUtil, FeatherUtil, OSystem, OErr, InfoList, History
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Feather", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Bomb()

# Get file names
outFile = sys.argv[1]
inFile = sys.argv[2:]
inDisk = 1
outDisk = 1

# Convert files into Images
inImage   = [Image.newPFImage(inFile[0], inFile[0], inDisk, 1, err)]
for x in inFile[1:]:
    inImage.append(Image.newPFImage("Input image", x, inDisk, 1, err))
OErr.printErrMsg(err, "Error initializing images")

# Create FFTs
FFTfor = FeatherUtil.PCreateFFT(inImage[0], 1)
FFTrev = FeatherUtil.PCreateFFT(inImage[0], 2)

# Create padded images for FFT size
print "********** Pad/interpolate Images to same grid *******"
padImage = []
i = 0;
# First image sets the grid for the following
print "Pad Loop",i,inFile[i]
name = "Pad"+inFile[i]
padIm = Image.newPFImage("Pad"+inFile[i], name, outDisk,  False, err)
OErr.printErrMsg(err, "Error creating padded image for "+inFile[i])
FeatherUtil.PPad(FFTfor, inImage[i], padIm, err)
OErr.printErrMsg(err, "Error padding image for "+inFile[i])
padImage.append(padIm)
i = i+1;
# Interpolate/pad others to same grid
for x in inImage[1:]:
    print "Pad Loop",i,inFile[i]
    name = "Pad"+inFile[i]
    padIm = Image.newPFImage("Pad"+inFile[i], name, outDisk,  False, err)
    OErr.printErrMsg(err, "Error creating padded image for "+inFile[i])
    FeatherUtil.PInterpol(inImage[i], padImage[0], padIm, err)
    OErr.printErrMsg(err, "Error writing padded image for "+inFile[i])
    padImage.append(padIm)
    i = i+1;

# Create masks in FArrays, first get weights from restoring beams/resolution
print "********** Create weighting masks *****************"
wtArrays = []
i = 0
for x in padImage:
    wtArray = FeatherUtil.PMakeBeamMask(x, FFTfor, err)
    OErr.printErrMsg(err, "Error creating weight array for "+inFile[i])
    wtArrays.append(wtArray)
    i = i + 1

# Create weight masks from FT of beams
# Weights are 1 with a Gaussian hole in the middle representing
# the uv coverage of the next smallest array/telescope
i = 0
for x in wtArrays:
    # replace with 1s
    FArray.PFill(x, 1.0)
    # Subtract next if it exists
    if i+1 < len(wtArrays):
        FArray.PSub(x, wtArrays[i+1], x)
    i = i + 1

# Make accumulation array and work array
print "********** Accumulate Weighted FFTs of images *****************"
accArray  = FeatherUtil.PCreateFFTArray(FFTfor)
CArray.PFill(accArray, [0.0,0.0])
workArray = FeatherUtil.PCreateFFTArray(FFTfor)

# Loop accumulating images
i = 0
for x in padImage:
    FeatherUtil.PAccumImage(FFTfor, x, wtArrays[i], accArray, workArray, err)
    OErr.printErrMsg(err, "Error accumulating array for "+inFile[i])
    i = i + 1

# FFT back to image domain
print "********** FFT back to image *****************"
resultArray = FArray.PClone(Image.PGetFArray(padImage[0]), err);
FeatherUtil.PBackFFT(FFTrev, accArray, resultArray, err)
OErr.printErrMsg(err, "Error back transforming to image plane ")

# Get normalization by repeating but using the padded images replaced by the beam
print "********** Get normalization using point models *****************"
CArray.PFill(accArray, [0.0,0.0])

# Loop accumulating normalization images
i = 0
for x in padImage:
    modelArray = Image.PGetFArray(x);  # get the FArray
    FeatherUtil.PCreateModel(x, modelArray)
    FeatherUtil.PAccumImage(FFTfor, x, wtArrays[i], accArray, workArray, err)
    OErr.printErrMsg(err, "Error accumulating array for "+inFile[i])
    i = i + 1

# FFT normalization image back to image domain
print "********** FFT back to image *****************"
workArray2 = FArray.PClone(resultArray, err)  # scratch FArray
OErr.printErrMsg(err, "Error cloning scratch array ")
FeatherUtil.PBackFFT(FFTrev, accArray, workArray2, err)
OErr.printErrMsg(err, "Error back transforming to image plane ")

# Do normalization from peak in workArray2
pos = [0, 0]
peak = FArray.PMax(workArray2, pos)
print "peak in normalization image",peak
norm = 1.0 / peak
FArray.PSMul(resultArray, norm)

# Extract to output from resultArray
outImage = Image.newPFImage(outFile,   outFile,  outDisk,  0, err)
Image.PClone(inImage[0], outImage, err)   # Same structure etc. as in 1
OErr.printErrMsg(err, "Error initializing output image")
FeatherUtil.PSubImage (padImage[0], resultArray, outImage, err)
OErr.printErrMsg(err, "Error writing output image")

# Do history to scratch image as table
inHistory  = History.History("history", inImage[0].List, err)
midHistory = History.History("history", padImage[0].List, err)
outHistory = History.History("history", outImage.List, err)
History.PCopyHeader(inHistory, midHistory, err)
# Add this programs history
midHistory.Open(History.READWRITE, err)
midHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
i = 0
for x in inFile:
    i = i + 1
    midHistory.WriteRec(-1,ObitSys.pgmName+" / input = "+x,err)
midHistory.Close(err)

# Copy to output header
History.PCopy2Header(midHistory, outHistory, err)
OErr.printErrMsg(err, "Error with history")

# Delete scratch files
for x in padImage:
    x.Zap(err)
    OErr.printErrMsg(err, "Error deleting padded images")


# Say something
print "Feathered",inFile,"to",outFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

