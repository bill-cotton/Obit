# python/Obit equivalent of AIPSish HGEOM for FITS images
# Arguments:
# 1) Name of input FITS image to be resampled (assumed
# 2) Name of template FITS image (defines grid for output)
# 3) Name of output FITS image
# All files in ../PythonData
# example:
# Python scriptHGeom.py 3C435BeditSep00.fits 3C435BeditJan02.fits \!HGeomOut.fits
#
import sys, Obit, Image, ImageUtil, OSystem, OErr, History
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("HGeom", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Bomb()

# Get Filenames (FITS)
inDisk   = 1
tmplDisk = 1
outDisk  = 1
inFile   = sys.argv[1]
tmplFile = sys.argv[2]
outFile  = sys.argv[3]

# Set data
inImage   = Image.newPFImage("Input image",    inFile,   inDisk,   1, err)
tmplImage = Image.newPFImage("Template image", tmplFile, tmplDisk, 1, err)
outImage  = Image.newPFImage("Output image",   outFile,  outDisk,  0, err)
Image.PClone(tmplImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

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
outHistory.WriteRec(-1,ObitSys.pgmName+" template= "+tmplFile,err)
outHistory.Close(err)
# Copy to output header - damn cfitsio
#inHistory = History.History("history", outImage.List, err)
#History.PCopy2Header(inHistory, outHistory, err)
#outHistory.Zap(err)   # Get rid of history table
OErr.printErrMsg(err, "Error with history")

# Say something
print "Interpolated",inFile,"to",outFile,"a clone of ",tmplFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

