# Test interferometric imaging script

import Obit, UV, UVImager, Image, ImageMosaic, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Imager", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (AIPS)  C346R422    .UVMOD .   1
inDisk = 1
inName = '1331+305'
inName = 'C346R422'
inClass = 'SPLIT'
inClass = 'UVSmal'
inSeq = 2
inClass = 'UVCOP'
inSeq = 1
inClass = 'UVMOD'
inSeq = 1
outName = "Pytho3C286"
outName = "PythonTest" # Bigfield
outClass = 'ID'
outDisk = 1
outSeq = 3  # Mosaic
outSeq = 4  # Robust
outSeq = 2000  # Single

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

# FUCK Python
#Obit.Bomb()
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
#Input['FOV']    = 500.0/3600.0  # Bigfield
#Input['doFull'] = True
#Input['NField'] = 0
#Single field:
Input['NField'] = 1
Input['xCells'] = 1.018
Input['yCells'] = 1.018
Input['nx'] = [64]
Input['ny'] = [64]
Input['RAShift']  = [300.0/3600.0]
Input['DecShift'] = [300.0/3600.0]
Input['RAShift']  = [0.0]
Input['DecShift'] = [0.0]

print "Create images"
outMosaic = UVImager.UVCreateImage(err, "myMosaic", Input)
OErr.printErrMsg(err, "Error creating image mosaic")

# Make Image
Input = UVImager.UVImageInput
Input['InData']    = scrData
Input['OutImages'] = outMosaic
Input['DoBeam']    =  True                                 
Input['DoWeight']  =  True                                 
Input['Robust']    =  0.0                                
Input['Robust']    =  7.0     # Natural weight
Input['WtSize'] = 512;          # Same weighting grid as AIPS
Input['UVTaper']   =  [0.0,0.0]
Input['Channel']   =  0
print "Make images"
UVImager.UVImage(err, Input)
OErr.printErrMsg(err, "Error making images")

# Flatten
#print "Flatten images"
#ImageMosaic.PFlatten(outMosaic, err)
OErr.printErrMsg(err, "Error flattening image mosaic")

# Say something
print "Imaged",inName,inClass,"to", outName,outClass

# Shutdown Obit
OErr.printErr(err)
del ObitSys
