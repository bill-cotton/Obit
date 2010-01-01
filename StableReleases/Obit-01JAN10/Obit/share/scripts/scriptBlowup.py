# python/Obit script to expand (shrink) the scale of an image
# Arguments:
# 1) scale factor
# 2) Name of input FITS image to be resampled (assumed
# 3) Name of output FITS image
# All files in .
#
import sys, Obit, Image, ImageUtil, OSystem, OErr, History
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("HGeom", 1, 100, 0, ["None"], 1, ["./"], True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Bomb()

# Get Filenames (FITS)
inDisk   = 1
tmplDisk = 1
outDisk  = 1
scale    = float(sys.argv[1])
inFile   = sys.argv[2]
outFile  = sys.argv[3]

# Set data
inImage   = Image.newPFImage("Input image",    inFile,   inDisk,   1, err)
outImage  = Image.newPFImage("Output image",   outFile,  outDisk,  0, err)
OErr.printErrMsg(err, "Error initializing")

# Update header scale times bigger image with scale times smaller increment
inHead = inImage.Desc.Dict
inHead["cdelt"][0] /= scale
inHead["cdelt"][1] /= scale
inHead["crpix"][0]  = (scale*(inHead["crpix"][0]-1.0)) + 1.0
inHead["crpix"][1]  = (scale*(inHead["crpix"][1]-1.0)) + 1.0
inHead["inaxes"][0] = int(scale*(inHead["inaxes"][0]))
inHead["inaxes"][1] = int(scale*(inHead["inaxes"][1]))
outImage.Desc.Dict = inHead
outImage.UpdateDesc(err)
OErr.printErrMsg(err, "Error setting output header")

# Interpolate
ImageUtil.PInterpolateImage(inImage, outImage, err)
OErr.printErrMsg(err, "Error interpolating")

# copy history 
inHistory  = History.History("history", inImage.List, err)
outHistory = History.History("history", outImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" inFile = "+inFile,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" scale= "+str(scale),err)
outHistory.Close(err)
# Copy to output header - damn cfitsio
#inHistory = History.History("history", outImage.List, err)
#History.PCopy2Header(inHistory, outHistory, err)
#outHistory.Zap(err)   # Get rid of history table
OErr.printErrMsg(err, "Error with history")

# Say something
print "Interpolated",inFile,"to",outFile,"scaling by ",scale

# Shutdown Obit
OErr.printErr(err)
del ObitSys

