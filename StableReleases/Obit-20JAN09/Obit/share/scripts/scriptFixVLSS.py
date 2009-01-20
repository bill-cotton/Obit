# python/Obit script to relabel the headers in VLSS images to J2000
# Arguments:
# 1) Name of Image to be updated

import sys, Obit, Image, ImageDesc, OSystem, OErr

# Init Obit
err=OErr.OErr()
#ObitSys=OSystem.OSystem ("Feather", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
ObitSys=OSystem.OSystem ("Feather", 1, 100, 1, ["../AIPSdata/"], 1, ["/home/nraoweb/cv/content/4mass/MAPS/work/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#Obit.Bomb()

# Get file names
inFile = sys.argv[1]
inDisk = 1
outDisk = 1
equinox = 2000.0

# Debug
print "input",inFile

# Convert file into Images
inImage = Image.newPImage(inFile, inFile, inDisk, 1, err)
OErr.printErrMsg(err, "Error initializing image")

# Open/get descriptor
Image.POpen(inImage, 3, err)
desc = Image.PGetDesc(inImage)

# update descriptor
descDict = ImageDesc.PGetDict(desc)       # Get descriptor as Python Dict
descDict["epoch"]   = equinox             # Update "epoch"
descDict["equinox"] = equinox             # Update equinox
descDict["origin"]  = "Obit to fix equinox"    # Update origin
ImageDesc.PSetDict(desc, descDict)        # update descriptor
Image.PClose(inImage, err)                # Close to update disk
OErr.printErrMsg(err, "Error writing updated header for "+Image.PGetName(inImage))

# Say something
print "Updated Equinox (EPOCH) in",inFile,"to",equinox

# Shutdown Obit
OErr.printErr(err)
del ObitSys

