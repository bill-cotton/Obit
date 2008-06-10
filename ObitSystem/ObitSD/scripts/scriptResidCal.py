# Program to self calibrate OTF data
import Obit, OTF, Image, OSystem, OErr, OTFGetSoln, InfoList, Table

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Python", 1, 103, 1, ["None"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
# Dirty
inFullFile  = "OTFDirtyFull.fits"       # input Full OTF data
inSubFile   = "OTFDirtySub.fits"        # input Full OTF data
#Clean
#inFullFile  = "OTFCleanFull.fits"       # input Full OTF data
#inSubFile   = "OTFCleanSub.fits"        # input Full OTF data

# Set data
fullData = OTF.newPOTF("Input data", inFullFile, disk, 1, err)
subData  = OTF.newPOTF("Input data", inSubFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")

# Calibration parameters
calType = "Filter"
solint = 5.0 / 86400.0
minRMS = 0.0
minEl = 0.0
calJy = [1.0,1.0]
dim = OTF.dim
dim[0] = 1
inInfo = OTF.POTFGetList(subData) 
InfoList.PInfoListAlwaysPutFloat(inInfo, "SOLINT",  dim, [solint])
InfoList.PInfoListAlwaysPutFloat(inInfo, "MINRMS",  dim, [minRMS])
InfoList.PInfoListAlwaysPutFloat(inInfo, "MINEL",   dim, [minEl])
dim[0] = len(calJy)
InfoList.PInfoListAlwaysPutFloat(inInfo, "CALJY", dim, calJy)
dim[0] = len(calType)
InfoList.PInfoListAlwaysPutString(inInfo, "calType", dim, [calType])
dim[0] = 1
solnTable = OTFGetSoln.POTFGetSolnFilter (subData, fullData, err)
soln = Table.PTableGetVer(solnTable)

# Update Cal table
# Soln2Cal parameters (most defaulted)
OTF.Soln2CalInput["InData"]  = fullData
OTF.Soln2CalInput["soln"]  = soln
# Use highest extant Cal table as input
oldCal = Obit.OTFGetHighVer(fullData.me, "OTFCal")
if oldCal == 0:   # Must not be one
    oldCal = -1
OTF.Soln2CalInput["oldCal"]  = oldCal
OTF.Soln2CalInput["newCal"]  = 0
OTF.Soln2Cal(err,OTF.Soln2CalInput)

# Shutdown 
OErr.printErr(err)
print 'Done, calibrated',inFullFile
