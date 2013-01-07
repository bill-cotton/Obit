# Test interferometric model subtraction /division script
# The argument, if given, is the data directory, defaults to "../testIt"
# Output UVdata should have ~zero phase and ~unit amplitude.
# Need some automated uv data comparison

import Obit, OSystem, OErr, sys

if len(sys.argv)>=2:
    dataDir = sys.argv[1]
else:
    dataDir = "../testIt/"
# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("UVSub", 1, 100, 1, ["../AIPSdata/"], 1, [dataDir], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Allow multiple threads
OSystem.PAllowThreads(2)  # 2 threads

import UV, UVImager, Image, ImageMosaic, SkyModel
from Obit import Bomb

# Files (FITS)
inDisk = 1
outDisk = 1
inFile  = 'UVSubTestIn.uvtab'
inModel = 'UVSubTestModIn.fits'
outFile = 'UVSubTestOut.uvtab'
masterDisk = 1
masterFile  = 'UVSubTestMaster.uvtab'

# Bombs away
#Bomb()

# Set data
print "Set data"
inData  = UV.newPFUV("Input uv data", inFile, inDisk, True, err)
inImage = Image.newPFImage("Input image",inModel, inDisk, True, err)
outData = UV.newPFUV("Output uv data", outFile, outDisk, False, err)
OErr.printErrMsg(err, "Error initializing")

# Make Mosaic
mosaic = ImageMosaic.newObit("Mosaic", 2, err)
#mosaic = ImageMosaic.newObit("Mosaic", 1, err)
OErr.printErrMsg(err, "Error making mosaic")
# Add image
ImageMosaic.PSetImage(mosaic, 0, inImage)
ImageMosaic.PSetImage(mosaic, 1, inImage)

# Make SkyModel model
model = SkyModel.PCreate("SkyModel", mosaic)
OErr.printErrMsg(err, "Error making SkyModel")

# control parameters
Input = SkyModel.UVSubInput
Input['InData']      = inData
Input['SkyModel']    = model
Input['OutData']     = outData
Input['doCalSelect'] = False
Input['Stokes']      = '    '
Input['CCVer']       = [2]
Input['BChan']       = 0
Input['EChan']       = 0
Input['BIF']         = 0
Input['EIF']         = 0
Input['Factor']      = 1.0
Input['Factor']      = 0.5
Input['Mode']        = 0  # Fastest
Input['Mode']        = 1  # DFT
#Input['Mode']        = 2  # Grid
#Input['Type']        = 1  # Image
Input['Type']        = 0  # CC
#replace data with model?
#Input['REPLACE'] = True

# Subtract
#print "Subtract"
#print "Replace"
#SkyModel.PSubUV(err, Input)
#OErr.printErrMsg(err, "Error subtracting")

# Divide
print "Divide"
SkyModel.PDivUV(err, Input)
OErr.printErrMsg(err, "Error dividing")

# Compare with master lie [rms diff]
masterData  = UV.newPFUV("Master UVData",   masterFile,  masterDisk,  True, err)
diff = UV.PUtilVisCompare(outData, masterData, err);
print "Comparison with master lie, fractional RMS R,I difference",diff

# Say something
#print "Subtracted",inModel,"From", inFile,"to",outFile
print "Divided FT of",inModel,"into", inFile,"to",outFile

#def PUtilUVUtilVisCompare (in1UV, in2UV, err):
# returns RMS

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
