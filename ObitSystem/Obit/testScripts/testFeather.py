# python/Obit script to feather together a series of images
# The argument, if given, is the data directory, defaults to "../testIt"

import Image, FArray, CArray, FFT, ImageUtil, FeatherUtil, OSystem, OErr, InfoList, History, sys
from OErr import Bomb


if len(sys.argv)>=2:
    dataDir = sys.argv[1]
else:
    dataDir = "../testIt/"

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Feather", 1, 100, 1, ["../AIPSdata/"], 1, [dataDir], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#Bomb()

# Get file names
outFile = 'FeatherTestOut.fits'
inFile = ['FeatherTestIn1.fits','FeatherTestIn2.fits']
inDisk = 1
outDisk = 1
masterDisk = 1
masterFile  = 'FeatherTestMaster.fits'

# Convert files into Images
inImage   = [Image.newPFImage(inFile[0], inFile[0], inDisk, 1, err)]
for x in inFile[1:]:
    inImage.append(Image.newPFImage("Input image", x, inDisk, 1, err))
OErr.printErrMsg(err, "Error initializing images")

# Create FFTs
FFTfor = FeatherUtil.PCreateFFT(inImage[0], 1)
FFTrev = FeatherUtil.PCreateFFT(inImage[0], 2)

# Create padded images for FFT size
print "Pad/interpolate Images to same grid"
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
print "Create weighting masks"
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
print "Accumulate Weighted FFTs of images"
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
print "FFT back to image"
resultArray = FArray.PClone(Image.PGetFArray(padImage[0]), err);
FeatherUtil.PBackFFT(FFTrev, accArray, resultArray, err)
OErr.printErrMsg(err, "Error back transforming to image plane ")

# Get normalization by repeating but using the padded images replaced by the beam
print "Get normalization using point models"
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

# Generate scratch file from inImage[0]
tmpImage  = Image.PScratch(inImage[0], err)
Image.POpen(tmpImage, Image.WRITEONLY, err)   # Open
OErr.printErrMsg(err, "Error cloning template")

# Extract to tmpImage from resultArray
FeatherUtil.PSubImage (padImage[0], resultArray, tmpImage, err)
OErr.printErrMsg(err, "Error writing scratch output image")

# Do history to scratch image as table
inHistory  = History.History("history", inImage[0].List, err)
outHistory = History.History("history", tmpImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
i = 0
for x in inFile:
    i = i + 1
    outHistory.WriteRec(-1,ObitSys.pgmName+" / input = "+x,err)
outHistory.Close(err)
OErr.printErrMsg(err, "Error with history")

# Create output with same geometry as inImage[0]
outImage = Image.newPFImage(outFile, outFile,  outDisk,  0, err)
Image.PClone(inImage[0], outImage, err)
OErr.printErrMsg(err, "Error initializing output image")

# Copy to quantized integer image with history
print "Write output image"
inHistory  = History.History("history", tmpImage.List, err)
Image.PCopyQuantizeFITS (tmpImage, outImage, err, inHistory=inHistory)

# Compare with master lie [rms diff, max abs diff, max. master]
masterImage  = Image.newPFImage("Master image",   masterFile,  masterDisk,  True, err)
diff = Image.PCompare(masterImage, outImage, err);
print "Comparison, rel. max. residual",diff[1]/diff[0], " rel RMS residual",diff[2]/diff[0]

# Delete scratch files
for x in padImage:
    Image.PZap(x, err);
    OErr.printErrMsg(err, "Error deleting padded images")

# Say something
print "Feathered",inFile,"to",outFile

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)

