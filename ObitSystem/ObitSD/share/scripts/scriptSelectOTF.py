# Program to select a section of OTF data
import OTF, OTFUtil, OSystem, Image, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("OTFSelect", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "GCLbandBothOTF.fits"         # input OTF data
outFile  = "!GCLbandSelOTF.fits"      # output OTF data

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
outData = OTF.newPOTF("Output data", outFile, disk, 0, err)
OTF.PClone(inData, outData, err)     # Same structure etc
OErr.printErrMsg(err, "Error initializing")

# Set time range
timerange = [0.0,1.0]
scans = [16,19]
dim = OTF.dim
dim[0] = 2; dim[1] = 1
inInfo = OTF.POTFGetList(inData)
InfoList.PAlwaysPutFloat(inInfo, "TIMERANGE",  dim, timerange)
InfoList.PAlwaysPutInt(inInfo, "SCANS",  dim, scans)
dim[0] = 1
InfoList.PAlwaysPutBoolean (inInfo,"doCalSelect" , dim, [True])
flagver=-1
gainuse=4
InfoList.PAlwaysPutInt (inInfo, "FLAGVER", dim, [flagver])
InfoList.PAlwaysPutInt (inInfo, "GAINUSE", dim, [gainuse])
itemp = -1
InfoList.PAlwaysPutInt (inInfo, "DOCALIB", dim, [itemp])
 
# Subtract image from scratch
OTF.PCopy(inData, outData,  err)
OErr.printErrMsg(err, "Error selecting data")

# Say something
print "Selected timerange",timerange,"from",inFile,"and wrote to",outFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys
