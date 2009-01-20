# Program to process an On/Off calibrator pair of scans
import OTF, OTFUtil, OSystem, OErr, InfoList
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FitOnOff", 1, 100, 1, ["None"], 1, ["FITSdata/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "ComaPbandAvg2OTF.fits"      # Part 2
#inFile   = "ComaPbandAvg1OTF.fits"     # Part 1

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Set scans
scan = [126,128] # Part 2
#scan = [40,42] # Part 1, first
#scan = [77,75] # Part 1, second

dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = OTF.PGetList(inData)
dim[0] = 2
InfoList.PAlwaysPutInt(inInfo, "Scan",  dim, scan)
#Bomb()
# Do fitting
OTFUtil.PFitOnOff(inData, -1,  err)
OErr.printErrMsg(err, "Error fitting cal scan")

# Give results
print "Fitting scans",scan,"in",inFile
stuff = InfoList.PGet (inInfo, "TRX")
print "Average TRx = ",stuff[4]
stuff = InfoList.PGet (inInfo, "CALJY")
print "Average CalJy = ",stuff[4]

# Shutdown Obit
OErr.printErr(err)

