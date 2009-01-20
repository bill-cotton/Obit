# Program to fit tipping scan
import OTF, OTFUtil, OSystem, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FitTip", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
#inFile   = "GCXbandDay1OTF.fits"         # input OTF data
inFile   = "GCCbandDay2OTF.fits"         # input OTF data

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Set scans
#GCXbandDay1OTF.fits scan = [13]
scan = [2011]
tsky  = [300.0]
minEl = [10.0]
#XBand tcal = [3.24,3.35]
tcal = [3.93,2.88]  # Cband
dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = OTF.PGetList(inData)
InfoList.PAlwaysPutInt(inInfo, "Scan",  dim, scan)
InfoList.PAlwaysPutFloat(inInfo, "TSKY",  dim, tsky)
InfoList.PAlwaysPutFloat(inInfo, "MINEL",  dim, minEl)
dim[0] = 2
InfoList.PAlwaysPutFloat (inInfo,"TCAL" , dim, tcal)

# Do fitting
OTFUtil.POTFUtilFitTip(inData,  err)
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
stuff = InfoList.Pet (inInfo, "TIPRMS")
print "RMS residuals (K) = ",stuff[4]

# Shutdown Obit
OErr.printErr(err)

