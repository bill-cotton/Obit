#testbed for CleanVis
import Obit, OErr, OSystem, CleanVis, UV, Image, ImageMosaic, History

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("CleanVis", 1, 100, 1, ["../AIPSdata/"], 1, ["../testIt/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Obit.Bomb()

# Get file names
inDisk   = 1
inFile = "UVImageTestIn.uvtab"
inFile = 'UVImageTestIn.uvtab'
inFile = 'FLS3minBADASS.uvtab'
inFile = 'WideField20.uvtab'
#debug inFile = "PointModel.uvtab"
outDisk  = 1
outFile = "CleanVisTestOut.fits"
masterDisk = 1
masterFile  = 'CleanVisMaster.fits'

# Imaging parameters
#FOV = 1.0/60.0  # Smaller
FOV = 10.0/60.0 # larger 20 cm field
FOV = 25.0/60.0 # BADASS
Stokes = 'I'
TimeRange = [0.0,10.0]
UVRange   = [0.0,0.0]
Robust    = 0.0
UVTaper   = [0.0,0.0]


# NVSS catalog
Catalog = 'NVSSVZ.FIT'
OutlierDist = 1.0   # Maximum distance to add outlyers (deg)
OutlierFlux = 0.001 # Minimum estimated outlier flux density (Jy)
OutlierSI   = -0.7  # Spectral index to estimate flux density
OutlierSize = 50    # Size of outlyer field

# Convert files into uvdata
inUV  = UV.newPFUV("input data", inFile, inDisk,  1, err)
OErr.printErrMsg(err, "Error initializing uvdata")

# Set inputs
Input = CleanVis.CleanInput
#Input['Sources'] = ["C346R422"]
Input['Sources'] = ["C346R424"]
Input['doCalSelect'] = True;
Input['Niter'] = 1000
Input['minFlux'] = 0.0001
Input['minPatch'] = 200
Input['Gain'] = 0.1
Input['CCVer'] = 1
#Input['BMAJ'] = 4.0/3600.0
#Input['BMIN'] = 4.0/3600.0
#Input['BPA'] = 0.0
# Data selection, cal etc.
Input['doCalSelect'] = True
Input['autoWindow']  = True
Input['Stokes']      = Stokes
Input['BChan']       = 0
Input['EChan']       = 0
Input['BIF']         = 0
Input['EIF']         = 0
Input['doCalib']     = 2
Input['gainUse']     = 2
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
Input['FOV']    = FOV 
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
Input['PBCor']       =  True

# debug
#Input['Niter'] = 200
#Input['minFlux'] = 0.000001
#Input['Gain'] = 0.1
#Input['FOV']    = 10.0/3660.0

#  Make Clean object 
Clean = CleanVis.PCreate ("Clean Object", inUV, err)
OErr.printErrMsg(err, "Error making clean")

# Set default clean window - not with autoWindow
#CleanVis.PDefWindow (Clean, err)
OErr.printErrMsg(err, "Error setting default window")

#  Clean 
CleanVis.PClean (Clean, err)
OErr.printErrMsg(err, "Error CLEANing")

# Full field image for output
outMosaic = CleanVis.PGetMosaic(Clean)
tmpImage  = ImageMosaic.PGetFullImage (outMosaic, err)

# Do history to scratch image as table
inHistory  = History.History("history", inUV.List, err)
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

# Say something
print "CLEAN: Dirty Image ",inFile,"Clean",outFile

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
