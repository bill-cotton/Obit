# Test interferometric imaging script

import Obit, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Imager", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

import UV, UVImager, Image, ImageMosaic
from Obit import Bomb

# Files (AIPS)  C346R422    .UVMOD .   1
inDisk = 1
inName = '1331+305'
inName = 'C346R422'
inClass = 'UVWAX'
inClass = 'UVMOD'
inClass = 'SPLIT'
inSeq = 1
outName = "Pytho3C286"
outName = "C346R422Obit"
outClass = 'ID'
outDisk = 1
outSeq = 2  # Single
outSeq = 3  # Mosaic
outSeq = 4  # Robust
outSeq = 2666  # Single
outSeq = 1666  # Single

# NVSS catalog
Catalog = 'NVSSVZ.FIT'
OutlierDist = 1.0  # Maximum distance to add outlyers (deg)
OutlierFlux = 0.01 # Minimum estimated outlier flux density (Jy)
OutlierSI   = -1.0 # Spectral index to estimate flux density
OutlierSize = 50   # Size of outlyer field

# Bombs away
#Bomb()

# Set data
print "Set data"
inData  = UV.newPAUV("Input data", inName, inClass, inDisk, inSeq, True, err) 
OErr.printErrMsg(err, "Error initializing")

# Make scratch version
print "Make scratch file"
scrData = UV.PScratch (inData, err)
OErr.printErrMsg(err, "Error making scratch file")

# Calibration/selection parameters
Input = UVImager.UVCalInput
Input['InData']      = inData
Input['OutData']     = scrData
Input['doCalSelect'] = True
Input['Stokes']      = 'I'
Input['BChan']       = 0
Input['EChan']       = 0
Input['BIF']         = 0
Input['EIF']         = 0
Input['doCalib']     = -1
Input['gainUse']     = 1
Input['flagVer']     = 1
Input['timeRange']   = [0.0,10.0]
Input['UVRange']     = [0.0,0.0]

# Calibrate
print "Calibrate"
UVImager.UVCal(err, Input)
OErr.printErrMsg(err, "Error calibrating to scratch file")

# Generate output ImageMosaic
Input = UVImager.UVCreateImageInput
Input['InData'] = scrData
Input['doBeam'] = True
Input['Type']   = 1
Input['Name']   = outName
Input['Class']  = outClass
Input['Seq']    = outSeq
Input['Disk']   = outDisk
#Mosaic:
#Input['FOV']    = 15.0/60.0  # 3C286
#Input['FOV']    = 20.0/60.0  # Fullfield
Input['FOV']    = 500.0/3600.0  # Bigfield
Input['doFull'] = True
Input['NField'] = 0
#Single field:
#Input['NField'] = 1
#Input['xCells'] = 1.018
#Input['yCells'] = 1.018
#Input['xCells'] = 0.5
#Input['yCells'] = 0.5
Input['nx'] = [1024]
Input['ny'] = [1024]
Input['RAShift']  = [0.0/3600.0]
Input['DecShift'] = [0.0/3600.0]
Input['Catalog'] = Catalog
Input['OutlierDist'] = OutlierDist
Input['OutlierFlux'] = OutlierFlux
Input['OutlierSI']   = OutlierSI
Input['OutlierSize'] = OutlierSize

print "Create images"
outMosaic = UVImager.UVCreateImage(err, "myMosaic", Input)
OErr.printErrMsg(err, "Error creating image mosaic")

# Make Image
Input = UVImager.UVImageInput
Input['InData']    = scrData
Input['OutImages'] = outMosaic
Input['DoBeam']    =  True                                 
Input['DoWeight']  =  True                                 
Input['Robust']    =  7.0     # Natural weight
Input['Robust']    =  0.0                                
Input['UVTaper']   =  [0.0,0.0]
Input['Channel']   =  0
print "Make images"
UVImager.UVImage(err, Input)
OErr.printErrMsg(err, "Error making images")

# Flatten
print "Flatten images"
ImageMosaic.PFlatten(outMosaic, err)
OErr.printErrMsg(err, "Error flattening image mosaic")

# Say something
print "Imaged",inName,inClass,"to", outName,outClass

# Shutdown Obit
OErr.printErr(err)
del ObitSys
