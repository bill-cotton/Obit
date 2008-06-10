# Program to select a section of OTF data
import OTF, OTFUtil, OSystem, Image, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("OTFSelect", 1, 103, 1, ["None"], 1, ["../FITSdata/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "GCXbandCalOTF.fits"         # input OTF data
outFile  = "!GCXbandSelOTF.fits"      # output OTF data

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
outData = OTF.newPOTF("Output data", outFile, disk, 0, err)
OTF.POTFClone(inData, outData, err)     # Same structure etc
OErr.printErrMsg(err, "Error initializing")

# Set time range
timerange = [0.13,0.2]
scans = [70,71]
dim = OTF.dim
dim[0] = 2; dim[1] = 1
inInfo = OTF.POTFGetList(inData)
InfoList.PInfoListAlwaysPutFloat(inInfo, "TIMERANGE",  dim, timerange)
InfoList.PInfoListAlwaysPutInt(inInfo, "SCANS",  dim, scans)
dim[0] = 1
InfoList.PInfoListAlwaysPutBoolean (inInfo,"doCalSelect" , dim, [True])

# Subtract image from scratch
OTF.POTFCopy(inData, outData,  err)
OErr.printErrMsg(err, "Error selecting data")

# Say something
print "Selected timerange",timerange,"from",inFile,"and wrote to",outFile

# Shutdown Obit
OErr.printErr(err)

