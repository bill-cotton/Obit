# Program to subtract an image from OTF data
import OTF, OTFUtil, OSystem, Image, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("OTFSub", 1, 103, 1, ["None"], 1, ["../FITSdata/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "GBTDaisyX2OTF.fits"         # input OTF data
imageFile= "PythonDaisyXClean.fits"     # Input image to subtract
imageFile= "PythonTDaisyX.fits"     # Input image to subtract
outFile  = "!GBTDaisyXSubOTF.fits"      # output OTF data
tmpFile  = "PythonScratch.fits"         # "scratch" file

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
outData = OTF.newPOTF("Output data", outFile, disk, 0, err)
OTF.POTFClone(inData, outData, err)     # Same structure etc
#tttData = OTF.newPOTF("tmp data",    "!"+tmpFile, disk, 0, err)
#OTF.POTFClone(inData, tttData, err)     # Same structure etc
modImage   = Image.newPImage("Model Image", imageFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Read image
imageData  = Image.PImageReadPlane (modImage, err)
imageDesc = Image.PImageGetDesc(modImage) 
OErr.printErrMsg(err, "Error reading image")

# Make scratch version
tmpData = OTF.POTFScratch (inData, err)
OErr.printErrMsg(err, "Error making scratch file")

# zero scratch data
scale = 0.0
offset = 0.0
OTFUtil.POTFUtilScale(inData, tmpData, scale, offset, err)
OErr.printErrMsg(err, "Error scaling OTF")

# Subtract image from scratch
OTFUtil.POTFUtilSubImage(tmpData, outData, imageData, imageDesc, err)
OErr.printErrMsg(err, "Error subtracting image")

# Say something
print "Subtracted",imageFile,"from",inFile," and wrote to",outFile

# Shutdown Obit
OErr.printErr(err)

