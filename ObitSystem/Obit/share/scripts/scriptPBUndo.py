# Remove primary beam correction
# python/Obit equivalent of AIPSish 1.0/PBCor
# Arguments:
# 1) Name of input FITS image
# 2) name of FITS image from which pointing is to be obtained
# 3) Name of the output FITS image

import  sys, Obit, Image, ImageUtil, OSystem, FArray, History, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("PBUndo", 1, 100, 1, ["../AIPSdata/"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Obit.Bomb()

# Get file names
inFile = sys.argv[1]
pntFile = sys.argv[2]
outFile = sys.argv[3]
inDisk = 1
outDisk = 1

# Set data
inImage   = Image.newPFImage("Input image",    inFile,   inDisk,  1, err)
pntImage  = Image.newPFImage("Pointing image", pntFile,  inDisk,  1, err)
outImage  = Image.newPFImage("Output image",   outFile,  outDisk, 0, err)
Image.PClone(inImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

# Make scratchfile for beam
beam = inImage.Scratch(err)
ImageUtil.PPBImage(pntImage, beam, err,antSize=25.,minGain=0.05)
plane = [1,1,1,1,1]
beam.GetPlane(None,plane,err)

# Multiply
nplane = inImage.Desc.Dict['inaxes'][2]
for iplane in range(0,nplane):
    plane = [iplane+1,1,1,1,1]
    inImage.GetPlane(None, plane, err)
    FArray.PMul(inImage.FArray, beam.FArray, inImage.FArray)
    outImage.PutPlane(inImage.FArray, plane, err)

beam.Zap(err)
OErr.printErrMsg(err, "Error correcting image")

# Copy History
inHistory  = History.History("inhistory",  inImage.List, err)
outHistory = History.History("outhistory", outImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" / Undo Primary beam correction",err)
outHistory.Close(err)

# Say something
print "PB removed from",inFile,"writing",outFile,", using pointing from",pntFile

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(err)

