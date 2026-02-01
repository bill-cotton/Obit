# python/Obit script to difference a pair of FITS images
# Arguments:
# 1) Name of output differenced FITS image
# 2-n) Names of the FITS images to be difference (first-second)
# Only does first plane
# All files are in the current working directory

import sys, Obit, Image, ImageDesc, FArray, ImageUtil, OSystem, OErr, InfoList, History

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("DiffImage", 1, 100, 1, ["None"], 1, ["./"], True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get file names
outFile = sys.argv[1]
inFile = sys.argv[2:]
inDisk = 1
outDisk = 1
print ("Files",outFile,inFile)

# Convert files into Images
inImage   = [Image.newPFImage(inFile[0], inFile[0], inDisk, 1, err)]
x=Image.newPFImage(inFile[0], inFile[0], inDisk, 1, err)
x.List.set("BLC",[1,1,1]) ; x.List.set("TRC",[0,0,1]) # Only fitst plane
for x in inFile[1:]:
    inImage.append(Image.newPFImage("Input image", x, inDisk, 1, err))
OErr.printErrMsg(err, "Error initializing images")

# Read accumulate images
accumArray = None
for x in inImage:
    # First
    x.Open(Image.READONLY,err)
    x.Read(err)
    fa = x.FArray
    # Copy first to accumArray
    if not accumArray:
        accumArray = FArray.PCopy(fa,err)
    else:
        # Subtract second
        FArray.PSub(accumArray, fa, accumArray)
    x.Close(err)
    OErr.printErrMsg(err, "Error accumulating images")
    del x

# Clone output from input
outImage=Image.PFArray2FITS(accumArray, outFile, err, outDisk=0)
ImageDesc.PCopyDesc(inImage[0].Desc, outImage.Desc, err)
outImage.UpdateDesc(err)
# Say something
print ("Differenced",inFile,"to",outFile)

inHistory  = History.History("history", inImage[0].List, err)
outHistory = History.History("history", outImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
i = 0
for x in inFile:
    i = i + 1
    outHistory.WriteRec(-1,ObitSys.pgmName+" / input = "+x,err)
outHistory.Close(err)
OErr.printErrMsg(err, "Error with history")

# Shutdown Obit
OErr.printErr(err)
del ObitSys

