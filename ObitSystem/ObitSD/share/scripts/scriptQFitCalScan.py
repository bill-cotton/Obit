# Program to select a section of OTF data
import OTF, OTFUtil, OSystem, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FitCal", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile = "QbandOTF.fits"

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Set scans
scan = [80,81,82,83] # 16Apr04 Q band test
tau0=0.03 #Qband - wild guess
ATemp = [0.735,0.722,0.722,0.722,0.722,0.722,0.722,0.722] # wild guess
dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = OTF.PGetList(inData)
InfoList.PAlwaysPutFloat(inInfo, "TAU0",  dim, [tau0])
dim[0] = len(ATemp)
InfoList.PAlwaysPutFloat (inInfo,"ATEMP" , dim, ATemp)
dim[0] = 4
InfoList.PAlwaysPutInt(inInfo, "Scan",  dim, scan)

# Do fitting
OTFUtil.PFitCal(inData, -1,  err)
OErr.printErrMsg(err, "Error fitting cal scan")

# Give results
print "Fitting scans",scan,"in",inFile
stuff = InfoList.PGet (inInfo, "TRX")
print "Average TRx = ",stuff[4]
stuff = InfoList.PGet (inInfo, "CALJY")
print "Average CalJy = ",stuff[4]
stuff = InfoList.PGet (inInfo, "Timeoff")
print "Average Timeoff = ",stuff[4]
print "Use opposite sign of Timeoff in DCR2OTR input"

# Shutdown Obit
OErr.printErr(err)

