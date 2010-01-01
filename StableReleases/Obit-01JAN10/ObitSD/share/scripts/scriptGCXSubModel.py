# Program to subtract an image from GC Xband OTF data
import OTF, OTFUtil, OSystem, Image, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("OTFSub", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "GCXbandSelOTF.fits"         # input OTF data
imageFile= "GCXDirty.fits"            # Input image to subtract
outFile  = "!GCXbandCalSubOTF.fits"     # output OTF data

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
outData = OTF.newPOTF("Output data", outFile, disk, 0, err)
OTF.PClone(inData, outData, err)     # Same structure etc
modImage   = Image.newPImage("Model Image", imageFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Use prior calibration
gainuse = 0
flagver = 0
docal = gainuse >= 0
if docal:
    dcal = 1
else:
    dcal = 0
dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = OTF.PGetList(inData)
InfoList.PPutBoolean(inInfo, "doCalSelect", dim, [docal], err)
InfoList.PPutInt(inInfo, "DOCALIB",  dim, [dcal], err)
InfoList.PPutInt(inInfo, "GAINUSE",  dim, [gainuse], err)
InfoList.PPutInt(inInfo, "FLAGUSE",  dim, [flagver], err)

# Read image
imageData  = Image.PReadPlane (modImage, err)
imageDesc = Image.PGetDesc(modImage) 
OErr.printErrMsg(err, "Error reading image")

# Subtract image from input
OTFUtil.PSubImage(inData, outData, imageData, imageDesc, err)
OErr.printErrMsg(err, "Error subtracting image")

# Say something
print "Subtracted",imageFile,"from",inFile," and wrote to",outFile

# Shutdown Obit
OErr.printErr(err)

