# Test History

import Obit, Image, History, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("TestHist", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (FITS)
inDisk = 1
inFile   = 'HiTest.fits'
outFile  = 'GaussTest.fits'
# AIPS test AIPS image user 100 C346R422    .POLI  .   2,, disk 8/142
aname = "C346R422"
aclass = "POLI"
aseq = 2

# Set data
inImage   = Image.newPImage("Input image", inFile,  inDisk,   1, err)
outImage   = Image.newPAImage("Output image", aname, aclass, inDisk, aseq, 1, err)
OErr.printErrMsg(err, "Error initializing")

# For debugging
#Obit.Bomb()

# Make history
inInfo = Image.PGetList(inImage)
outInfo = Image.PGetList(outImage)
history = History.History("history", outInfo, err)
OErr.printErrMsg(err, "Error initializing history")
History.POpen(history, 3, err)
History.PTimeStamp(history,"Label",err)
History.PTimeStamp(history,"Something longer",err)
str = History.PReadRec(history,144,err)
print "Read ",str
History.PWriteRec(history,-1," Hello from Python",err)
History.PClose(history, err)
OErr.printErrMsg(err, "Error writing history")
# Copy header stuff
#history2 = History.History("history", outInfo, err)
#History.PCopy2Header(history,history2,err)
OErr.printErrMsg(err, "Error copying history")

# Say something
print "Added Timestamp to History",inFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

