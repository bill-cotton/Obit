# Test interferometric imaging script
# The argument, if given, is the data directory, defaults to "../testIt"

import Obit, UV, UVImager, Image, ImageMosaic, OSystem, History, OErr, sys
from OErr import Bomb

if len(sys.argv)>=2:
    dataDir = sys.argv[1]
else:
    dataDir = "../testIt/"

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Imager", 1, 100, 1, ["../AIPSdata/"], 1, [dataDir], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Allow multiple threads
OSystem.PAllowThreads(2)  # 2 threads

# Files (AIPS)  C346R422    .SPLIT .   1
inDisk = 1
inFile = 'UVImageTestIn.uvtab'
outDisk = 1
outFile = "UVImageTestOut.fits"
masterDisk = 1
masterFile  = 'UVImageTestMaster.fits'

# Imaging parameters
FOV = 0.1/60.0  # Smaller
FOV = 10.0/60.0 # larger 20 cm field
Stokes = 'I'
TimeRange = [0.0,10.0]
UVRange   = [0.0,0.0]
Robust    = 0.0
UVTaper   = [0.0,0.0]


# NVSS catalog
Catalog = 'NVSSVZ.FIT'
OutlierDist = 1.0  # Maximum distance to add outlyers (deg)
OutlierFlux = 0.01 # Minimum estimated outlier flux density (Jy)
OutlierSI   = -1.0 # Spectral index to estimate flux density
OutlierSize = 50   # Size of outlyer field

# If debugging
#Bomb()

# Set data
print "Set data"
inData  = UV.newPFUV("Input data", inFile, inDisk, True, err) 
OErr.printErrMsg(err, "Error initializing")

# Calibration/selection parameters
Input = UVImager.UVCreateImagerInput
# Data selection, cal etc.
Input['InData']      = inData
Input['doCalSelect'] = True
Input['Stokes']      = Stokes
Input['BChan']       = 0
Input['EChan']       = 0
Input['BIF']         = 0
Input['EIF']         = 0
Input['doCalib']     = -1
Input['gainUse']     = 1
Input['flagVer']     = 1
Input['timeRange']   = TimeRange
Input['UVRange']     = UVRange
# Imaging control
Input['doBeam'] = True
Input['Type']   = 0     # FITS
Input['Name']   = 'Image'  # Have to pretend it's AIPS
Input['Class']  = 'Class'
Input['Seq']    = 1
Input['Disk']   = outDisk
Input['FOV']    = 10.0/60.0
Input['FOV']    = FOV # Smaller
Input['xCells'] = 1.0
Input['yCells'] = 1.0
Input['doFull'] = True
Input['NField'] = 0
Input['RAShift'] = [0.0]
Input['DecShift'] = [0.0]
Input['Catalog'] = Catalog
Input['OutlierDist'] = OutlierDist
Input['OutlierFlux'] = OutlierFlux
Input['OutlierSI']   = OutlierSI
Input['OutlierSize'] = OutlierSize
Input['Robust']      =  Robust
Input['UVTaper']     =  UVTaper

print "Create UV Imager"
myImager = UVImager.PUVCreateImager(err, "myMosaic", Input)
OErr.printErrMsg(err, "Error creating image mosaic")

print "Weight uv data"
UVImager.PUVWeight(myImager, err)
OErr.printErrMsg(err, "Error creating image mosaic")

# Make Image
print "Make images"
UVImager.PUVImage(myImager, err, doWeight=False, doBeam=True)
OErr.printErrMsg(err, "Error making images")

# Flatten
print "Flatten images"
UVImager.PUVFlatten(myImager, err)
OErr.printErrMsg(err, "Error flattening image mosaic")

# Full field image for output
outMosaic = UVImager.PGetMosaic(myImager, err)
tmpImage  = ImageMosaic.PGetFullImage (outMosaic, err)

# Do history to scratch image as table
inHistory  = History.History("history", inData.List, err)
outHistory = History.History("history", tmpImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" inFile = "+inFile,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" FOV = "+str(FOV),err)
outHistory.WriteRec(-1,ObitSys.pgmName+" Stokes = "+Stokes,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" TimeRange = "+str(TimeRange),err)
outHistory.WriteRec(-1,ObitSys.pgmName+" UVRange = "+str(UVRange),err)
outHistory.WriteRec(-1,ObitSys.pgmName+" Robust = "+str(Robust),err)
outHistory.WriteRec(-1,ObitSys.pgmName+" UVTaper = "+str(UVTaper),err)
outHistory.Close(err)
OErr.printErrMsg(err, "Error with history")

# output image
outImage  = Image.newPFImage("Output image", outFile,  outDisk,  False, err)
Image.PClone(tmpImage, outImage, err)   # Same structure etc.

# Copy to quantized integer image with history
print "Write output image"
inHistory  = History.History("history", tmpImage.List, err)
Image.PCopyQuantizeFITS (tmpImage, outImage, err, inHistory=inHistory)

# Compare with master lie [rms diff, max abs diff, max. master]
masterImage  = Image.newPFImage("Master image",   masterFile,  masterDisk,  True, err)
diff = Image.PCompare(outImage, masterImage, err);
print "Comparison, rel. max. residual",diff[1]/diff[0], " rel RMS residual",diff[2]/diff[0]

# Cleanup
ImageMosaic.PZapImage(outMosaic, -1, err)      # delete mosaic images
tmpImage.Zap(err)                              # Full field too
OErr.printErrMsg(err, "Error cleaning up files")

# Say something
print "Imaged",inFile,"to", outFile

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
