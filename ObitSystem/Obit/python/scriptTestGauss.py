# python/Obit equivalent of AIPSish IMMOD

import Obit, Image, FArray, ImageUtil, FeatherUtil, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("TestIMMOD", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (FITS)
inDisk = 1
inFile   = 'GaussTest.fits'
outDisk = 1
outFile  = '!TestIMMOD.fits'

# Model
amp = 1.0
Cen = [38.0, 41.0];
GauMod = [6.0, 3.0, 45.0]

# Set data
inImage   = Image.newPImage("Input image",  inFile,   inDisk,   1, err)
outImage  = Image.newPImage("Output image", outFile,  outDisk,  0, err)
Image.PClone(inImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

# For debugging
#Obit.Bomb()

# Get FArray
outArray = Image.PGetFArray (outImage)

# Add Model
#FArray.PEGauss2D(outArray, amp, Cen, GauMod)

# Replace image with resolution
FeatherUtil.PCreateModel(inImage, outArray)

# Write
Image.PWrite(outImage, err)
OErr.printErrMsg(err, "Error writing")

# Say something
print "Replaced",inFile,"with Gaussian model to",outFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

