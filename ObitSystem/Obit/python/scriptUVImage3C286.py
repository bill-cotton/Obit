# Test interferometric imaging script

import Obit, UV, UVImager, Image, ImageMosaic, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Imager", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (AIPS)  C346R422    .SPLIT .   1
inDisk = 1
inName = 'C346R422'
inName = '1331+305'
inClass = 'SPLIT'
inSeq = 1
outName = "PythonTest" # Bigfield
outName = "Pytho3C286"
outClass = 'ID'
outDisk = 1
outSeq = 1

# FUCK Python
#Obit.Bomb()

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
Input['FOV']    = 10.0/60.0 # C346R422
Input['FOV']    = 2.0/60.0  # 3C286
Input['FOV']    = 0.0
Input['doFull'] = True
Input['NField'] = 1
#Input['xCells'] = 0.01
#Input['yCells'] = 0.01
#Input['nx'] = [256]
#Input['ny'] = [256]
Input['RAShift'] = [60.0/3600.0]
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
