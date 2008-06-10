# Program to make an image from OTF data and CLEAN it
import OTF, Image, OErr, OSystem

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem  ("Python", 1, 103, 1, ["None"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Imaging parameters
OTF.ImageInput["ra"] = 24.42216            # Center RA
OTF.ImageInput["dec"] = 33.15863           # Center Dec
OTF.ImageInput["xCells"] = 20.0 / 3600.0   # "X" cell spacing, deg
OTF.ImageInput["yCells"] = 20.0 / 3600.0   # "Y" cell spacing, deg
OTF.ImageInput["nx"] = 64                  # number of cells in X
OTF.ImageInput["ny"] = 64                  # number of cells in X
OTF.ImageInput["gainuse"] = 0             # Which cal table to apply, -1 = none
OTF.ImageInput["flagver"] = -1             # Which flag table to apply, -1 = none

# CLEAN parameters
OTF.CleanInput["niter"]    = 300    # Number of iterations
OTF.CleanInput["gain"]     = 0.05   # CLEAN loop gain
OTF.CleanInput["beamsize"] = 80.2   # CLEAN restoring beam size in asec
OTF.CleanInput["beamsize"] = 80.0   # CLEAN restoring beam size in asec
OTF.CleanInput["minFlux"]  = 0.3    # Minimum image brightness to CLEAN

# Files
disk = 1
# Dirty
inFile   = "OTFDirtyFull.fits"        # input OTF data
dirtFile = "!DirtyDirty.fits"       # output dirty image file
beamFile = "DirtyBeam.fits"            # input dirty beam image
cleanFile= "!DirtyClean.fits"   # output CLEAN image
#Clean
#inFile   = "OTFCleanFull.fits"        # input OTF data
#dirtFile = "!CleanDirty.fits"       # output dirty image file
#cleanFile= "!CleanClean.fits"   # output CLEAN image

# Set data
inData = OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")

# Make image
OTF.ImageInput["InData"]  = inData
OTF.ImageInput["disk"]    = disk
OTF.ImageInput["OutName"] = dirtFile
#OTF.input(OTF.ImageInput)
Dirty = OTF.makeImage(err)

# Images for cleaning
Beam  = Image.newPImage("Dirty Beam", beamFile, disk, 1, err)
Clean = Image.newPImage("Clean Image", cleanFile, disk, 0, err)
OTF.CleanInput["Dirty"] = Dirty
OTF.CleanInput["Beam"]  = Beam
OTF.CleanInput["Clean"] = Clean
OErr.printErrMsg(err, "Error creating clean images")

# Clean it
#OTF.input(OTF.CleanInput)
OTF.Clean (err, OTF.CleanInput)
OErr.printErrMsg(err, "Error Cleaning")

# Shutdown Obit
OErr.printErr(err)
del ObitSys
