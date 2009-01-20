# Program to do test reading DCR data 
# of GBT Qband test and write calibrated output.
import OTF, OSystem, OErr, GBTDCROTF, InfoList
from Obit import Bomb

inDir  = "../../RPrestage/"     # Input data directory
outDir = "../PythonData/"      # Output data directory
outFile = "!QTestReadOTF.fits" # output OTF data
inDisk  = 1
outDisk = 2
scans = ["2004_04_16_01:20:04","2004_04_16_01:20:44","2004_04_16_01:21:56"]

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Python", 1, 103, 1, ["None"], 2, [inDir,outDir], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# debug
#Bomb()

# Set data
outOTF = OTF.newPOTF("Output data", outFile, outDisk, 0, err)
OErr.printErrMsg(err, "Error initializing output")

# Create converter
myGDO = GBTDCROTF.newPGBTDCROTF("DCR_OTF_Converter", outOTF, err)
OErr.printErrMsg(err, "Error initializing converter")

for x in scans:
    GBTDCROTF.PConvert(myGDO, inDisk, x, err)
    OErr.printErrMsg(err, "Error converting scan "+x)

# Shutdown 
OErr.printErr(err)
print 'Done: wrote data in',outFile
