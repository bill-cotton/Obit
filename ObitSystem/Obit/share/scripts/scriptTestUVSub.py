# Test interferometric model subtraction script

import Obit, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("UVSub", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

import UV, UVImager, Image, ImageMosaic, SkyModel
from Obit import Bomb

# Files (FITS)
inDisk = 1
outDisk = 1
inFile  = 'PModel.uvtab'
inModel = 'PModelQ.fits'
outFile = 'UVSubTestOut.uvtab'

# Bombs away
#Bomb()

# Set data
print "Set data"
inData  = UV.newPFUV("Input uv data", inFile, inDisk, True, err)
inImage = Image.newPFImage("Input image",inModel, inDisk, True, err)
outData = UV.newPFUV("Output uv data", outFile, outDisk, True, err)
OErr.printErrMsg(err, "Error initializing")

# Make Mosaic
mosaic = ImageMosaic.newObit("Mosaic", 1, err)
OErr.printErrMsg(err, "Error making mosaic")
# Add image
ImageMosaic.PSetImage(mosaic, 0, inImage)

# Make SkyModel
model = SkyModel.newObit("SkyModel", err)
OErr.printErrMsg(err, "Error making SkyModel")
SkyModel.PSetMosaic(model, mosaic)

# control parameters
Input = SkyModel.PSubUVInput
Input['InData']      = inData
Input['SkyModel']    = model
Input['OutData']     = outData
Input['doCalSelect'] = False
Input['STOKES']      = 'I'
Input['BCHAN']       = 0
Input['ECHAN']       = 0
Input['BIF']         = 0
Input['EIF']         = 0
Input['DOCALIB']     = -1
Input['GAINUSE']     = 1
Input['FLAGVER']     = 1
Input['TIMERANGE']   = [0.0,10.0]
Input['UVRANGE']     = [0.0,0.0]

# Subtract
print "Subtract"
SkyModel.PSubUV(err, Input)
OErr.printErrMsg(err, "Error subtracting")

# Say something
print "Subtracted",inModel,"From", inFile,"to",outFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys
