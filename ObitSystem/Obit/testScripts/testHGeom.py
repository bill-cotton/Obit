# Test for python/Obit equivalent of AIPSish HGEOM

import Obit, Image, ImageUtil, OSystem, History, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("HGeom", 1, 100, 1, ["../AIPSdata/"], 1, ["../testIt/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (FITS)
inDisk = 1
inFile   = 'HGeomTestIn.fits'
tmplDisk = 1
tmplFile = 'HGeomTestTmpl.fits'
masterDisk = 1
masterFile  = 'HGeomTestMaster.fits'
outDisk = 1
outFile  = '!HGeomTestOut.fits'

# Set data
inImage   = Image.newPFImage("Input image",    inFile,   inDisk,   True, err)
tmplImage = Image.newPFImage("Template image", tmplFile, tmplDisk, True, err)
outImage  = Image.newPFImage("Output image",   outFile,  outDisk,  False, err)
Image.PClone(tmplImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

# For debugging
#Obit.Bomb()

# Generate scratch file from tmplFile
tmpImage  = Image.PScratch(tmplImage, err)
Image.POpen(tmpImage, Image.WRITEONLY, err)   # Open
OErr.printErrMsg(err, "Error cloning template")

# Interpolate
ImageUtil.PInterpolateImage(inImage, tmpImage, err)
OErr.printErrMsg(err, "Error interpolating")

# Do history to scratch image as table
inHistory  = History.History("history", inImage.List, err)
outHistory = History.History("history", tmpImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" / input = "+inFile,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" / template = "+tmplFile,err)
outHistory.Close(err)
OErr.printErrMsg(err, "Error with history")

# Copy to quantized integer image with history
print "Write output image"
inHistory  = History.History("history", tmpImage.List, err)
Image.PCopyQuantizeFITS (tmpImage, outImage, err, inHistory=inHistory)

# Compare with master lie [rms diff, max abs diff, max. master]
masterImage  = Image.newPFImage("Master image",   masterFile,  masterDisk,  True, err)
diff = Image.PCompare(outImage, masterImage, err);
print "Comparison, rel. max. residual",diff[1]/diff[0], " rel RMS residual",diff[2]/diff[0]

# Say something
print "Interpolated",inFile,"to",outFile,"a clone of ",tmplFile

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)

