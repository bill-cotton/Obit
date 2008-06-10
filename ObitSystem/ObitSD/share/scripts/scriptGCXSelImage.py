# Program to image selected data
import OTF, Image, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Python", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "GCXbandSelOTF.fits"       # input OTF data
dirtFile = "!XDirty.fits"           # output dirty image file

# Set data
inData = OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")

# Imaging parameters
OTF.ImageInput["InData"]  = inData
OTF.ImageInput["disk"]    = disk
OTF.ImageInput["OutName"] = dirtFile
OTF.ImageInput["ra"]  = 266.2540           # Center RA
OTF.ImageInput["dec"] = -29.38028          # Center Dec
OTF.ImageInput["xCells"] = 20.0 / 3600.0   # "X" cell spacing, deg
OTF.ImageInput["yCells"] = 20.0 / 3600.0   # "Y" cell spacing, deg
OTF.ImageInput["nx"] = 500                # number of cells in X
OTF.ImageInput["ny"] = 500                # number of cells in X
OTF.ImageInput["gainuse"] = 0              # Which cal table to apply, -1 = none
OTF.ImageInput["flagver"] = -1             # Which flag table to apply, -1 = none

OTF.input(OTF.ImageInput)
Dirty = OTF.makeImage(err)

# Shutdown 
OErr.printErr(err)
print 'Imaged',inFile,'to',dirtFile
