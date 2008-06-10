# Program to subtract an image from OTF data
import OTF, OTFUtil, OSystem, Image, OErr, InfoList, FArray

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("OTFSub", 1, 103, 1, ["None"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
# Dirty
inFile   = "OTFDirtyFull.fits"        # input OTF data
imageFile= "DirtyDirty.fits"
outFile  = "!OTFDirtySub.fits"      # output OTF data
#Clean
#inFile   = "OTFCleanFull.fits"        # input OTF data
#imageFile= "CleanClean.fits"
#outFile  = "!OTFCleanSub.fits"      # output OTF data

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
outData = OTF.newPOTF("Output data", outFile, disk, 0, err)
OTF.POTFClone(inData, outData, err)     # Same structure etc
modImage   = Image.newPImage("Model Image", imageFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Read image
imageData  = Image.PImageReadPlane (modImage, err)
imageDesc  = Image.PImageGetDesc(modImage) 
OErr.printErrMsg(err, "Error reading image")

# Apply prior calibration as requested
flagver = -1
gainuse = 0
dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = OTF.POTFGetList(inData)
InfoList.PInfoListAlwaysPutInt (inInfo, "FLAGVER", dim, [flagver])
InfoList.PInfoListAlwaysPutInt (inInfo, "GAINUSE", dim, [gainuse])
itemp = gainuse >= 0
InfoList.PInfoListAlwaysPutInt (inInfo, "DOCALIB", dim, [itemp])
doCalSelect = (gainuse >= 0) or (flagver>0)
InfoList.PInfoListAlwaysPutBoolean (inInfo,"doCalSelect" , dim, [doCalSelect])

# clip image below minFlux
minFlux = 0.1
print "Clip image below ", minFlux
imageData = Image.PImageGetFArray(modImage)
FArray.PFArrayClip (imageData, minFlux, 1.0e20, 0.0)


# Subtract image from input
OTFUtil.POTFUtilSubImage(inData, outData, imageData, imageDesc, err)
OErr.printErrMsg(err, "Error subtracting image")

# Say something
print "Subtracted",imageFile,"from",inFile," and wrote to",outFile

# Shutdown Obit
OErr.printErr(err)

