# Program to subtract an image from another
import Image, OSystem, OErr, FArray

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("ImageSub", 1, 103, 1, ["None"], 1, ["../FITSdata/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "PythonDaisyXClean.fits"         # input image
in2File  = "PythonTDaisyX.fits"             # Input image to subtract
outFile  = "!PythonDaisyXCleanSub.fits"     # output image

# Set data
inImage  = Image.newPImage("Input image",  inFile, disk, 1, err)
in2Image = Image.newPImage("Input 2 image",  in2File, disk, 1, err)
outImage=  Image.newPImage("Output image", outFile, disk, 0, err)
Image.PImageClone(inImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

# Read images
imageData  = Image.PImageReadPlane (inImage, err)
image2Data = Image.PImageReadPlane (in2Image, err)

# Check compatability
if not FArray.PFArrayIsCompatable (imageData, image2Data):
    raise RuntimeError,'Images incompatable'

# Subtract images
FArray.PFArraySub(imageData, image2Data, imageData)    # result in imageData

# Tell RMS difference
print "RMS difference",FArray.PFArrayRMS(imageData)

# Write output
Image.PImageWritePlane (outImage, imageData, err)

# Say something
print "Subtracted",in2File,"from",inFile," and wrote to",outFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

