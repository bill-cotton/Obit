# python/Obit equivalent of AIPSish 1.0/PBCor
# Arguments:
# 1) Name of input image
# 2) name of image from which pointing is to be obtained
# 3) Name of the output image
# example
# Python scriptPBApply.py ComaAGBT.fits ComaA-P-Bsub2.fits ComaAGBTPBcor.fits

import  sys, Obit, Image, ImageUtil, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("PBApply", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Obit.Bomb()

# Get file names
inFile = sys.argv[1]
pntFile = sys.argv[2]
outFile = sys.argv[3]
inDisk = 1
outDisk = 1

# Set data
inImage   = Image.newPImage("Input image",    inFile,   inDisk,   1, err)
pntImage  = Image.newPImage("Pointing image", pntFile, inDisk, 1, err)
outImage  = Image.newPImage("Output image",   outFile,  outDisk,  0, err)
Image.PClone(inImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

# do it - defaulting plane, antena size
ImageUtil.PPBApply(inImage, pntImage, outImage, err)
OErr.printErrMsg(err, "Error correcting image")

# Say something
print "PB applied to",inFile,"writing",outFile,", using pointing from",pntFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

