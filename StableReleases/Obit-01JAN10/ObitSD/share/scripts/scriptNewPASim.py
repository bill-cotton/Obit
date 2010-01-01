# Program to self calibrate new Penn Array simulation data
# Also applies tmospheric calibration
# Takes input data in OTF FITS format
# The script iteratively inproves the derived image and calibration.
# Finally, writes a calibrated version of the data

import OTF, OTFGetSoln, Image, OSystem, InfoList, OErr
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("NewPASim", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "NewPASimOTF.fits"         # input OTF data
dirtFile = "!NewPASimDirty.fits"      # output dirty image file
outFile  = "!NewPASimOTFCal.fits"     # output OTF data

# Set data
inData = OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")

# select only RALongMap scan data
target = ["RALongMap"]
dim = OTF.dim
dim[0] = len(target[0]); dim[1]=1; dim[2]=1; dim[3]=1;
InfoList.PAlwaysPutString(inData.List, "TARGETS", dim, target)

# Bomb if needed for debugging
#Bomb()

# Imaging parameters
OTF.ImageInput["InData"]  = inData
OTF.ImageInput["disk"]    = disk
OTF.ImageInput["OutName"] = dirtFile
OTF.ImageInput["ra"]  = 189.00             # Center RA
OTF.ImageInput["ra"]  =  22.05             # Center RA
OTF.ImageInput["dec"] = 0.0                # Center Dec
OTF.ImageInput["xCells"] = 2.0 / 3600.0    # "X" cell spacing, deg
OTF.ImageInput["yCells"] = 2.0 / 3600.0    # "Y" cell spacing, deg
OTF.ImageInput["nx"] = 800                 # number of cells in X
OTF.ImageInput["ny"] = 800                 # number of cells in X
OTF.ImageInput["gainuse"] = 0              # Which cal table to apply, -1 = none
OTF.ImageInput["flagver"] = 1              # Which flag table to apply, -1 = none

# Atmospheric calibration parameters
OTF.AtmCalInput["InData"]  = inData # Input data object
OTF.AtmCalInput["solint"]  =  600.0 # Solution interval(sec)
OTF.AtmCalInput["tau0"]    =  0.15  # Zenith opacity
OTF.AtmCalInput["tau0"]    =  0.0   # Simulations have geometry all wrong
OTF.AtmCalInput["aTemp"]   =[290.0] # Atm. temperature (K) per detector
OTF.AtmCalInput["tRx"]     =[100.0] # Recvr. Temp. (K)
OTF.AtmCalInput["calJy"]   =[0.03]  # Cal. signal in Jy
OTF.AtmCalInput["minEl"]   = -90.   # Min elev. (deg)

# Calibration parameters (some reset in loop)
OTF.ResidCalInput["InData"]  = inData      # Input data object
OTF.ResidCalInput["solType"] = "MultiBeam" # Solution type
OTF.ResidCalInput["solint"]  = 500.0       # Solution interval (sec)
OTF.ResidCalInput["minFlux"] = 1.0         # Minimum image brightness to use in model
OTF.ResidCalInput["Clip"]    = 0.1         # Minimum image brightness to use in model
OTF.ResidCalInput["gainuse"] = 0           # Prior calibration, 0> highest
OTF.ResidCalInput["minEl"] = -90.0         # minimum elevation

# Soln2Cal parameters (most defaulted)
OTF.Soln2CalInput["InData"]  = inData       # Input data object
OTF.Soln2CalInput["oldCal"]  = 0            # Use highest extant Cal table as input

# Initialize calibtration tables
# delete any prior calibration tables
OTF.ClearCal(inData,err)

# Create an initial dummy table with a interval 1/4 of the shortest
# Filter type solution interval.
inter = 0.25
inter = 0.125
inter = 0.05
inter = 0.025
OTFGetSoln.POTFGetDummyCal (inData, inData, inter, 1, 1, err)

# Atmospheric calibration
print "Atmospheric calibration"
OTF.AtmCal(err, OTF.AtmCalInput)
OTF.Soln2Cal(err, OTF.Soln2CalInput)  # Apply

# Individual detector offsets
print "Individual detector offsets"
OTF.ResidCalInput["solType"] = "Filter" # Solution type
OTF.ResidCalInput["solint"]  = 60.0
OTF.ResidCalInput["solint"]  = 6.0
OTF.ResidCalInput["minFlux"] = 1000.0
OTF.ResidCalInput["Clip"] = 1000.0
OTF.SelfCal(err, OTF.ImageInput, "None", OTF.ResidCalInput, OTF.Soln2CalInput)

# Solution intervals, Residual clipping pairs
# The residual clipping is needed to suppress artifacts due to large
# residuals near bright point sources; it should start off large
# and decrease to several times the noise.
# It not needed, set to a large value (1.0e20)
# These define the number and parameters of the iterations
OTF.ResidCalInput["solType"] = "MultiBeam" # Solution type
soln = [(60.0,30.0), (45.0,5.0), (30.0,0.5), (20.0,0.2), (10.0,0.2), (5.0,0.2), (5.0,0.2), (5.0,0.2)]
soln = [(10.0,1000.0,100.0), (5.0,1000.0,5.0), (2.0,1000.0,5.0), (1.0,1000.0,5.0)]
soln = [(1.0,10.0,100.0), (0.5,1.0,5.0), (0.25,0.5,5.0), (0.1,0.1,1.0)]
soln = [(1.0,10.0,100.0), (0.5,1.0,5.0), (0.25,0.5,5.0)]
soln = [(0.50,1.0,0.1), (0.40,0.05,0.01), (0.20,0.01,0.003), (0.20,0.005,0.001)]

# Loop over self cal
count=0
for si,mf,fl in soln:
    count = count+1
    print "\n *** Self calibration loop ",count,"si=",si
    OTF.ResidCalInput["solint"]  = si*4.0
    OTF.ResidCalInput["minFlux"] = fl
    OTF.ResidCalInput["Clip"] = mf
    # Longer term individual detector
    OTF.ResidCalInput["solType"] = "Filter" # Solution type
    OTF.SelfCal(err, OTF.ImageInput, "None", OTF.ResidCalInput, OTF.Soln2CalInput)
    #OTF.SelfCal(err, OTF.ImageInput, OTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
    # multibeam atmospheric solution
    OTF.ResidCalInput["solint"]  = si
    OTF.ResidCalInput["solType"] = "MultiBeam" # Solution type
    OTF.SelfCal(err, OTF.ImageInput, "None", OTF.ResidCalInput, OTF.Soln2CalInput)

print 'Finished with loop, final image'

# Final image
image = OTF.makeImage(err,  OTF.ImageInput)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Write calibrated output
print "Write calibrated data "
outData = OTF.newPOTF("Output data", outFile, disk, False, err)
OTF.PClone(inData, outData, err)     # Same structure etc
OErr.printErrMsg(err, "Error initializing output")

# Set time range
dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = inData.List
#dim[0] = 2
#timerange=[0.0,1.0]
#scans=[1,10]
#InfoList.PAlwaysPutFloat(inInfo, "TIMERANGE",  dim, timerange)
#InfoList.PAlwaysPutInt(inInfo, "SCANS",  dim, scans)
dim[0] = 1
InfoList.PAlwaysPutBoolean (inInfo,"doCalSelect" , dim, [True])
flagver=-1
gainuse=0
InfoList.PAlwaysPutInt (inInfo, "FLAGVER", dim, [flagver])
InfoList.PAlwaysPutInt (inInfo, "GAINUSE", dim, [gainuse])
itemp = 1
InfoList.PAlwaysPutInt (inInfo, "DOCALIB", dim, [itemp])
 
# Copy/calibrate
OTF.PCopy(inData, outData,  err)
OErr.printErrMsg(err, "Error selecting data")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Shutdown 
OErr.printErr(err)
print 'Done, calibrated',inFile,'image',dirtFile
del ObitSys
