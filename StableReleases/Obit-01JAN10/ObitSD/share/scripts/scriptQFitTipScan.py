# Program to fit tipping scan
import OTF, OTFUtil, OSystem, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FitTip", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile = "QbandOTF.fits"

# Set data
inData  = OTF.newPOTF(inFile,  inFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Set scans
scan = [15]
tsky  = [300.0]
minEl = [10.0]
tcal = [9.6, 9.8,5.65,5.55,9.6, 9.8,5.65,5.55]  # Qband
dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = OTF.PGetList(inData)
InfoList.PAlwaysPutInt(inInfo, "Scan",  dim, scan)
InfoList.PAlwaysPutFloat(inInfo, "TSKY",  dim, tsky)
InfoList.PAlwaysPutFloat(inInfo, "MINEL",  dim, minEl)
dim[0] = len(tcal)
InfoList.PAlwaysPutFloat (inInfo,"TCAL" , dim, tcal)

# Do fitting
OTFUtil.PFitTip(inData,  err)
OErr.printErrMsg(err, "Error fitting tipping scan")

# Give results
print "Fitting scan",scan[0],"in",inFile,"above elev",minEl[0]
stuff = InfoList.PGet (inInfo, "TAU0")
print "tau0 = ",stuff[4]
stuff = InfoList.PGet (inInfo, "TRX")
print "TRx = ",stuff[4],"cal units"
print "TRx = ",stuff[4][0]*tcal[0],stuff[4][1]*tcal[1],"K"
stuff = InfoList.PGet (inInfo, "ATEMP")
print "ATemp = ",stuff[4],"cal units per airmass"
print "ATemp = ",stuff[4][0]*tcal[0],stuff[4][1]*tcal[1],"K"
stuff = InfoList.PGet (inInfo, "TIPRMS")
print "RMS residuals (K) = ",stuff[4]

# Shutdown Obit
OErr.printErr(err)

