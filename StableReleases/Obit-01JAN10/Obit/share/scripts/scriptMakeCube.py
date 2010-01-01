# Add images into planes of a cube
# Arguments:
# 1) Name of output combined image
# 2-n) Names of the images to be combined, in order of increasing depth in the image
#      output is cloned from first named
# All files in ./
# Example: 
#Python MakeMovie.py \!NGC315CDirtyMovie.fits NGC315CDirtyAtm.fits NGC315CDirtyBase.fits NGC315CDirty12sec.fits NGC315CDirty6sec.fits NGC315CDirty3sec.fits

import sys, OSystem, OErr, InfoList, Image, ImageDesc, CleanOTF
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem  ("Image", 1, 100, 1, ["None"], 1, ["."], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Bomb if needed for debugging
#Bomb()

# Get file names
outFile = sys.argv[1]
inFile = sys.argv[2:]
inDisk = 1
outDisk = 1

# Convert files into Images
inImage   = [Image.newPFImage(inFile[0], inFile[0], inDisk, True, err)]
for x in inFile[1:]:
    inImage.append(Image.newPFImage("Input image", x, inDisk, True, err))
OErr.printErrMsg(err, "Error initializing images")

# output image
outImage = Image.newPFImage(outFile, outFile, outDisk, False, err)
OErr.printErrMsg(err, "Error initializing output image")

# Add number of planes
id = inImage[0].Desc
od = outImage.Desc
ImageDesc.PCopyDesc(id, od, err)
dict = od.Dict
inaxes = id.Dict["inaxes"]
inaxes[2] = len(inImage)
dict["naxis"]  = len(inaxes)
dict["inaxes"] = inaxes
dict["bitpix"] = -32
od.Dict = dict

# Glue 'em together
i = 1
outImage.Open(2,err)

# Loop over input images
for x in inImage:
    x.Open(1,err)
    fa = x.FArray
    x.Read(err)
    OErr.printErrMsg(err, "Error initializing output image")
    plane = [i,1,1,1,1]
    outImage.WriteFA(fa,err)
    OErr.printErrMsg(err, "Error initializing output image")
    x.Close(err)
    i=i+1
    OErr.printErrMsg(err, "Error initializing output image")
outImage.Close(err) # Close output

# Say something
print "Combined",inFile,"into",outFile

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
