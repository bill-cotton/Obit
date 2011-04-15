# python/Obit script to average a series of FITS images
# Arguments:
# 1) Name of output combined FITS image
# 2-n) Names of the FITS images to be combined
# All files are in the current working directory

import sys, Obit, Image, FArray, CArray, FFT, ImageUtil, OSystem, OErr, InfoList, History

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("AvgImage", 1, 100, 1, ["None"], 1, ["./"], True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get file names
outFile = sys.argv[1]
inFile = sys.argv[2:]
inDisk = 1
outDisk = 1

# Convert files into Images
inImage   = [Image.newPFImage(inFile[0], inFile[0], inDisk, 1, err)]
for x in inFile[1:]:
    inImage.append(Image.newPFImage("Input image", x, inDisk, 1, err))
OErr.printErrMsg(err, "Error initializing images")

# Read accumulate images
accumArray = None
countArray = None
for x in inImage:
    x.Open(Image.READONLY,err)
    x.Read(err)
    fa = x.FArray
    if not accumArray:
        accumArray = FArray.PCopy(fa,err)
        countArray = FArray.PClone(fa,err)   # Count of elements
        FArray.PFill(countArray, 1.0)
        FArray.PBlank(countArray, fa, countArray)
        FArray.PDeblank(countArray, 0.0)
    else:
        FArray.PSumArr(accumArray, fa, accumArray)
        cntTmp = FArray.PClone(fa,err)   # Count of elements
        FArray.PFill(cntTmp, 1.0)
        FArray.PBlank(cntTmp, fa, cntTmp)
        FArray.PDeblank(cntTmp, 0.0)
        FArray.PSumArr(countArray, cntTmp, countArray)
        FArray.PUnref(cntTmp)     # delete temporary counting array
    x.Close(err)
    OErr.printErrMsg(err, "Error accumulating images")
    del x

# Normalize
FArray.PDivClip (accumArray, countArray, 1.0, accumArray)

# Clone output from first input
outImage = Image.newPFImage(outFile,   outFile,  outDisk,  False, err)
Image.PClone(inImage[0], outImage, err)   # Same structure etc. as in 1
OErr.printErrMsg(err, "Error initializing output image")

# Write output
outImage.Open(Image.WRITEONLY,err)
outImage.WriteFA(accumArray,err)
outImage.Close(err)
OErr.printErrMsg(err, "Error writing "+outFile)
# Say something
print "Averaged",inFile,"to",outFile

# Do history 
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

