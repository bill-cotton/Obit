# python/Obit script to quantize a FITS image
# The output image is the input image written as 16 or 32 bit integers
# with a quantization level a given fraction of the minimum RMS in any plane.
# Arguments:
# 1) Name of input FITS image to be quantized
# 2) Name of output FITS image
# 3) The fraction of the RMS for the quantization [def. 0.25]
# All files in ./
# example:
# Python scriptQuantize.py HGeomTestOut.fits \!QuantTestOut.fits 0.25
#
import sys, Image, History, OSystem, OErr
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Quantize", 1, 100, 1, ["../AIPSdata/"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Bomb()

# Get Filenames (FITS)
inDisk   = 1
outDisk  = 1
inFile   = sys.argv[1]
outFile = sys.argv[2]
# Fraction of RMS if given
if (len(sys.argv)>=4):
    fract  = float(sys.argv[3])
else:
    fract = 0.25

# Set data
inImage   = Image.newPFImage("Input image",    inFile,   inDisk, True, err)
outImage  = Image.newPFImage("Output image",   outFile,  outDisk,  False, err)
#Image.PClone(inImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

# Copy to quantized integer image with history
print "Write quantized output image"
Image.PCopyQuantizeFITS (inImage, outImage, err, fract=fract)

# Copy History
inHistory  = History.History("inhistory",  inImage.List, err)
outHistory = History.History("outhistory", outImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" / Quantized at "+str(fract)+" RMS",err)
outHistory.Close(err)
# Copy back to header
inHistory  = History.History("inhistory",  outImage.List, err)
History.PCopy2Header (inHistory, outHistory, err)
OErr.printErrMsg(err, "Error with history")

# Say something
print "Quantized",inFile,"to",outFile,"with fraction",fract

# Shutdown Obit
OErr.printErr(err)
del ObitSys

