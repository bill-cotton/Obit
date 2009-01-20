# Program to make an image from OTF data and CLEAN it
import OTF

# Init Obit
err=OTF.Obit.ObitErrCreate()
ObitSys=OTF.Obit.Startup ("Python", 1, 103, 1, ["None"], 1, ["../FITSdata/"], 1, 0, err)
OTF.printErrMsg(err, "Error with Obit startup")

# Imaging parameters
OTF.ImageInput["ra"] = 24.42216            # Center RA
OTF.ImageInput["dec"] = 33.15863           # Center Dec
OTF.ImageInput["xCells"] = 20.0 / 3600.0   # "X" cell spacing, deg
OTF.ImageInput["yCells"] = 20.0 / 3600.0   # "Y" cell spacing, deg
OTF.ImageInput["nx"] = 64                  # number of cells in X
OTF.ImageInput["ny"] = 64                  # number of cells in X
OTF.ImageInput["gainuse"] = -1             # Which cal table to apply, -1 = none
OTF.ImageInput["flagver"] = -1             # Which flag table to apply, -1 = none

# CLEAN parameters
OTF.CleanInput["niter"]    = 300    # Number of iterations
OTF.CleanInput["gain"]     = 0.05   # CLEAN loop gain
OTF.CleanInput["beamsize"] = 87.4   # CLEAN restoring beam size in asec
OTF.CleanInput["minFlux"]  = 0.3    # Minimum image brightness to CLEAN

# Files
disk = 1
inFile   = "GBTDaisyX2OTF.fits"        # input OTF data
dirtFile = "!PythonTDaisyX.fits"       # output dirty image file
beamFile = "DirtyBeam.fits"            # input dirty beam image
cleanFile= "!PythonDaisyXClean.fits"   # output CLEAN image

# Set data
inData =OTF.newOTF("Input data", inFile, disk, 1, err)
OTF.printErrMsg(err, "Error creating input data object")

# Get descriptor, test rewrite values
desc = OTF.Obit.OTFGetDesc(inData)
dic = OTF.Obit.OTFDescGetDict(desc)
OTF.input(dic)
dic["obsra"] = 1.234
dic["bunit"] = "noUnit"
OTF.Obit.OTFDescSetDict(desc,dic)
dic2 = OTF.Obit.OTFDescGetDict(desc)
OTF.input(dic2)
