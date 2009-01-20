# Program to self calibrate OTF data
import OTF, Image, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Python", 1, 103, 1, ["None"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "GBTDaisyX2OTF.fits"        # input OTF data
dirtFile = "!CombDirty.fits"           # output dirty image file
beamFile = "DirtyBeam.fits"            # input dirty beam image
cleanFile= "!CombClean.fits"           # output CLEAN image

# Images for cleaning
Beam  = Image.newPImage("Dirty Beam", beamFile, disk, 1, err)
Clean = Image.newPImage("Clean Image", cleanFile, disk, 0, err)

# Set data
inData =OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")

# delete any prior calibration tables
OTF.ClearCal(inData,err)

# Imaging parameters
OTF.ImageInput["InData"]  = inData
OTF.ImageInput["disk"]    = disk
OTF.ImageInput["OutName"] = dirtFile
OTF.ImageInput["ra"]  = 24.42216           # Center RA
OTF.ImageInput["dec"] = 33.15863           # Center Dec
OTF.ImageInput["xCells"] = 20.0 / 3600.0   # "X" cell spacing, deg
OTF.ImageInput["yCells"] = 20.0 / 3600.0   # "Y" cell spacing, deg
OTF.ImageInput["nx"] = 64                  # number of cells in X
OTF.ImageInput["ny"] = 64                  # number of cells in X
OTF.ImageInput["gainuse"] = 0              # Which cal table to apply, -1 = none
OTF.ImageInput["flagver"] = -1             # Which flag table to apply, -1 = none

# CLEAN parameters
OTF.CleanInput["Beam"]     = Beam          # Dirty beam
OTF.CleanInput["Clean"]    = Clean         # Clean image
OTF.CleanInput["niter"]    = 200           # Number of iterations
OTF.CleanInput["gain"]     = 0.05          # CLEAN loop gain
OTF.CleanInput["beamsize"] = 80.0          # CLEAN restoring beam size in asec
OTF.CleanInput["minFlux"]  = 0.3           # Minimum image brightness to CLEAN

# Calibration parameters
OTF.ResidCalInput["InData"]  = inData      # Input data object
OTF.ResidCalInput["solType"] = "Filter"    # Solution type
OTF.ResidCalInput["solint"]  = 500.0       # Solution interval (sec)
OTF.ResidCalInput["minFlux"] = 1.0         # Minimum image brightness to use in model
OTF.ResidCalInput["gainuse"] = 0           # Prior calibration, 0> highest

# Soln2Cal parameters (most defaulted)
OTF.Soln2CalInput["InData"]  = inData       # Input data object
OTF.Soln2CalInput["oldCal"]  = 0            # Use highest extant Cal table as input

# Solution intervals, minFlux pairs
soln = [(10.0,0.5)]
soln = [(10.0,0.45), (10.0,0.25), (5.0,0.1), (5.0,0.08), (5.0,0.05)]
#soln = [(40.0,0.5), (25.0,0.4), (10.0,0.3)]
#soln = [(40.0,0.5), (25.0,0.4)]
# Loop over self cal
count=0
for si,mf in soln:
    count = count+1
    print "\n *** Self calibration loop ",count
    OTF.ResidCalInput["solint"]  = si
    OTF.ResidCalInput["minFlux"] = mf
    #OTF.SelfCal(err, OTF.ImageInput, OTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
    OTF.SelfCal(err, OTF.ImageInput, None, OTF.ResidCalInput, OTF.Soln2CalInput)

print 'Finished with loop, final image and clean'

# Final image and Clean
image = OTF.makeImage(err,  OTF.ImageInput)
#OTF.Clean (err,  OTF.CleanInput)

# Shutdown 
OErr.printErr(err)
print 'Done, calibrated',inFile,'image',cleanFile
