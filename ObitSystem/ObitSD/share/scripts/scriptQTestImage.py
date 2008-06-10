# Program to make an image from OTF data and CLEAN it
# Test for Richard Prestages Q band data
import OTF, Image, OErr, OSystem, InfoList
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem  ("Python", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Bomb if needed for debugging
#Bomb()

# Files
disk = 1
inFile   = "QbandCalOTF.fits"        # input OTF data
dirtFile = "!QTest.fits"          # output dirty image file
dirtFile = "!QTestSecond.fits"          # output dirty image file
#beamFile = "DirtyBeam.fits"      # input dirty beam image
#cleanFile= "!GCXXClean.fits"     # output CLEAN image

# Set data
inData = OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")

# Select scans
dim = OTF.dim
dim[0] = 2; dim[1] = 1
scans = [108,133]
scans = [134,184] # Second scan
inInfo = OTF.PGetList(inData)
InfoList.PAlwaysPutInt(inInfo, "SCANS",  dim, scans)
dim[0] = 1
InfoList.PAlwaysPutBoolean (inInfo,"doCalSelect" , dim, [True])
print "inInfo",inInfo.me


# Imaging parameters
OTF.ImageInput["InData"]  = inData
OTF.ImageInput["disk"]    = disk
OTF.ImageInput["OutName"] = dirtFile
OTF.ImageInput["ra"]  = 194.04646          # Center RA
OTF.ImageInput["dec"] = -5.78892           # Center Dec
OTF.ImageInput["xCells"] = 4.0 / 3600.0    # "X" cell spacing, deg
OTF.ImageInput["yCells"] = 4.0 / 3600.0    # "Y" cell spacing, deg
OTF.ImageInput["nx"] = 100                 # number of cells in X
OTF.ImageInput["ny"] = 100                 # number of cells in Y
OTF.ImageInput["gainuse"] = 0              # Which cal table to apply, -1 = none
OTF.ImageInput["flagver"] = 1              # Which flag table to apply, -1 = none


# Make image
OTF.input(OTF.ImageInput)
Dirty = OTF.makeImage(err)

# Shutdown Obit
OErr.printErr(err)
