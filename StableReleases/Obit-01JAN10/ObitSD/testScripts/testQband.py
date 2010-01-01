# Test program to do Q Band basic atmospheric and instrumental calibration
# Should be run from ObitTalk
import OSystem, OErr, InfoList, Image, Table, TableUtil, History, ODisplay
import ObitTask, FITSDir
import OTF, OTFUtil, CleanOTF, OTFGetSoln, OTFGetSoln
import GBTUtil 

# Init Obit
err=OErr.OErr()
FITS = ["../testIt/"]
ObitSys=OSystem.OSystem  ("QBand", 1, 1, 1, ["None"], len(FITS), FITS, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")


#######################################################################################
# Define parameters

# Root of data directory
DataRoot = None  # To suppress attempt to update from archive

# Target to image
target = ["1256-057"]
scans = [108,133]  # First Obs, surface corr on


# Define data
# OTF file
inFile   = "testQBandOTF.fits"      # input OTF data
inDisk   = 1

# Use target name to define output files
outFile  = "!testQBandCalOTF.fits"           # Output calibrated data
dirtFile = "!testQBandOut.fits"              # scratch output dirty image
tmpFile  = "!tmpQBandOut.fits"               # working clean/model
cleanFile= "!testQBandClean.fits"           # Final output image (Gaussian beam)
BeamFile = "KBandGaussBeam.fits"             # Dirty Gaussian beam
masterFile = "testQBandMaster.fits"   # Standard lie

# List of feeds to use, empty list = all
feeds=[1,2,5,6]  # Sig feed
feeds=[3,4,7,8]  # Ref feed
feeds = []       # Both feeds

# Image Info
timerange=[0.0,1.0]               # time range in days, All times
cells = 4.0                       # cell spacing in asec
nx = 100                          # no. cells in x
ny = 100                          # no. cells in y

niter   = 1000                    # Number of iteration of CLEAN
gain    = 0.1                     # CLEAN loop gain
minFlux = 0.001                   # Minimum image brightness to CLEAN
flagver = 1                       # Flag table

# Calibration info
solInt = 5.0                    # Min solInt in seconds
BLInt  = 120.0                  # Baseline filter time in sec
AtmInt = 60.0                   # Atmospheric filter time in sec
tau0   = 0.1                    # Zenith opacity
# Air temperature (units of cal)
ATemp  = [1.643, 1.526, 2.999, 2.779, 1.686, 1.514, 2.964, 2.872]  # From tipping scan
# Rx Temp in K,  one for all or per detector
TRx    = [11.51, 8.601, 18.21, 14.31, 9.461, 8.278, 21.68, 15.38]  # From tipping scan
# Cal values in Jy, one for all or per detector
CalJy = [9.9, 10.1, 6.17, 5.97, 10.0, 11.18, 5.75, 5.55]  # From 3C286 very crude

# set of iterations (soln int, clip level for residuals, minimum flux in model)
soln = [(12*solInt,20.0,0.00,"Offset"),(9*solInt,10.0,0.00,"Offset"),\
        (2*solInt,2.0,0.00,"Offset"),(solInt,1.0,0.00,"Offset")]

#######################################################################################
#----------------------------------------------------------------
# Strategy
# 1) 
# 2) 
# 3) 
# 4) 
#----------------------------------------------------------------

# Set data
inOTF =OTF.newPOTF("Input data", inFile, inDisk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")
inInfo = inOTF.List

################################## Set parameters #####################################
# Get position from OTF
pos =  GBTUtil.GetTargetPos(inOTF, target[0], err)
ra  = pos[0]                      # ra of center
dec = pos[1]                      # dec of center

# dirty beam
dirtyBeam  = Image.newPFImage("dirty beam", BeamFile,  inDisk,  True, err)
OErr.printErrMsg(err, "Error initializing dirty beam")

# Image display
disp = ODisplay.ODisplay("ObitView", "ObitView", err)

# Imaging parameters
OTF.ImageInput["InData"]  = inOTF
OTF.ImageInput["disk"]    = inDisk
OTF.ImageInput["OutName"] = dirtFile
OTF.ImageInput["Beam"]    = dirtyBeam
OTF.ImageInput["ra"]  = ra                   # Center RA
OTF.ImageInput["dec"] = dec                  # Center Dec
OTF.ImageInput["xCells"] = cells             # "X" cell spacing
OTF.ImageInput["yCells"] = cells             # "Y" cell spacing
OTF.ImageInput["nx"] = nx                    # number of cells in X
OTF.ImageInput["ny"] = ny                    # number of cells in Y
OTF.ImageInput["gainUse"] = 0                # Which cal table to apply, -1 = none
OTF.ImageInput["flagVer"] = flagver          # Which flag table to apply, -1 = none
OTF.ImageInput["minWt"] = 1.0e-2             # Minimum weight in imaging - includes data weight
OTF.ImageInput["ConvType"] = 5               # Convolving fn = pillbox, 3 = Gaussian,
                                             # 4 = Exp*Sinc, 5 = Spherodial wave
# Calibration parameters (some reset in loop)
OTF.ResidCalInput["InData"]  = inOTF       # Input data object
OTF.ResidCalInput["solType"] = "MultiBeam" # Solution type
OTF.ResidCalInput["solInt"]  = 500.0       # Solution interval (sec)
OTF.ResidCalInput["minFlux"] = 1.0         # Minimum image brightness to CLEAN
OTF.ResidCalInput["Clip"]    = 0.1         # Minimum image brightness to use in model
OTF.ResidCalInput["gainUse"] = 4           # Prior calibration, 0> highest
OTF.ResidCalInput["minEl"]   = -90.0       # minimum elevation
OTF.ResidCalInput["flagVer"] = flagver     # Which flag table to apply, -1 = none

################################## Initialize calibration #############################
# delete any prior calibration tables
print "Remove previous calibration"
OTF.ClearCal(inOTF,err)

# Create an initial dummy table with a interval 1/4 of the shortest
# Filter type solution interval.
inter = solInt/4
OTFGetSoln.POTFGetDummyCal (inOTF, inOTF, inter, 1, 1, err)

################################## Gain/Weight calibration ############################
# Gain/Weight calibration
print "Gain/Weight  calibration"
inInfo = inOTF.List
inInfo.set("calJy", CalJy)
inInfo.set("doWate", True) # Do weight calibration


OTFGetSoln.POTFGetSolnPARGain(inOTF, inOTF, err)
OTF.Soln2CalInput["InData"]  = inOTF        # Input data object
OTF.Soln2CalInput["oldCal"]  = 1            # Use Gain cal output
OTF.Soln2CalInput["newCal"]  = 2            # (Re)Use
OTF.Soln2Cal(err, OTF.Soln2CalInput)        # Apply

################################## Target calibration #################################

inInfo.set("Targets", target)         # select only target data
inInfo.set("Stokes", "    ")          # Set Stokes
inInfo.set("timeRange", timerange)    # Set timerange
inInfo.set("Scans", scans)            # Select scans
inInfo.set("doCalSelect", True)       

################################## Baseline filter ####################################
print "Baseline Filter"
solint = BLInt/86400.0
inInfo.set("solInt", solint)
#DEBUG clip = 1.015 # DEBUG was 0.015
#DEBUG inInfo.set("Clip", clip)
inInfo.set("doCalSelect", True)
inInfo.set("flagVer", flagver)
gainuse=2
inInfo.set("gainUse", gainuse)
inInfo.set("doCalib", 1)
OTFGetSoln.PFilter(inOTF, inOTF, err)
# Soln2Cal parameters for filter cal (most defaulted)
OTF.Soln2CalInput["InData"]  = inOTF       # Input data object
OTF.Soln2CalInput["oldCal"]  = 2           # Use Atm cal output
OTF.Soln2CalInput["newCal"]  = 3           # (Re)Use 3
OTF.Soln2Cal(err, OTF.Soln2CalInput)  # Apply

############################## Common atmosphere + offset #############################
print "Common atmosphere removal"
inInfo.set("doCalSelect", True)
inInfo.set("flagVer", flagver)
gainuse=3
inInfo.set("gainUse", gainuse)
inInfo.set("doCalib", 1)
solint = AtmInt/86400.0
inInfo.set("solInt", solint)
clipsig = 5.0
inInfo.set("ClipSig", clipsig)
plotDet = -10
inInfo.set("plotDet", plotDet)
OTFGetSoln.PMBBase(inOTF, inOTF, err)
# Soln2Cal parameters for baseline cal (most defaulted)
OTF.Soln2CalInput["InData"]  = inOTF       # Input data object
OTF.Soln2CalInput["oldCal"]  = 3            # Use Atm cal output
OTF.Soln2CalInput["newCal"]  = 4           # (Re)Use 3
OTF.Soln2Cal(err, OTF.Soln2CalInput)  # Apply

############################## Create initial image ###################################
print "Make Dirty Image"

#gainuse=3
#inInfo.set("gainUse", gainuse)
#inInfo.set("doCalib", 1)
#OTF.ImageInput["gainUse"] = 4                # Which cal table to apply, -1 = none
if len(feeds)>0:
    inInfo.set("Feeds", feeds)                # Which feeds
Dirty = OTF.makeImage(err, OTF.ImageInput)
#ODisplay.PImage(disp, Dirty, err)
OErr.printErrMsg(err, "Error displaying initial image")

################################# CLEAN parameters ####################################
# Images for cleaning
CleanTmp = Image.newPFImage("Clean Image", tmpFile, inDisk, False, err)

# Create CleanOTF
CleanObj = CleanOTF.PCreate("Clean", Dirty, dirtyBeam, CleanTmp, err)
OErr.printErrMsg(err, "Error creating CLEAN object")

# CLEAN parameters
CleanOTF.CleanInput["CleanOTF"] = CleanObj     # Clean object
CleanOTF.CleanInput["Patch"]    = 40           # Beam patch
CleanOTF.CleanInput["Niter"]    = niter        # Number of iterations
CleanOTF.CleanInput["Gain"]     = gain         # CLEAN loop gain
CleanOTF.CleanInput["BeamSize"] = 8.0/3600.0   # CLEAN restoring beam size in deg
CleanOTF.CleanInput["minFlux"]  = minFlux      # Minimum image brightness to CLEAN
CleanOTF.CleanInput["CCVer"]    = 1            # Clean components table version
CleanOTF.CleanInput["Feeds"]    = feeds        # list of feeds to use
CleanOTF.CleanInput["noResid"]  = False        # Include residuals in calibration
CleanOTF.CleanInput["noResid"]  = True         # Don't include residuals in calibration

#  CLEAN window
CleanOTF.PAddWindow(CleanObj, [-1,6, nx/2+1, ny/2+1], err)
# Mars
# CleanOTF.PAddWindow(CleanObj, [-1,4, 127, 133], err)
# CleanOTF.PAddWindow(CleanObj, [-1,3, 134, 131], err)

window = CleanOTF.PGetWindow(CleanObj)
OErr.printErrMsg(err, "Error setting window")

# Reset Soln2Cal parameters for self cal
OTF.Soln2CalInput["InData"]  = inOTF       # Input data object
OTF.Soln2CalInput["oldCal"]  = 4           # Use bl cal
OTF.Soln2CalInput["newCal"]  = 5           # (Re)Use 5

# Solution intervals, Residual clipping pairs
# The residual clipping is needed to suppress artifacts due to large
# residuals near bright point sources; it should start off large
# and decrease to several times the noise.
# It not needed, set to a large value (1.0e20)
# These define the number and parameters of the iterations
OTF.ResidCalInput["solType"] = "Filter" # Solution type - may be flakey
OTF.ResidCalInput["solType"] = "MultiBeam" # Solution type
#OTF.ResidCalInput["solType"] = "Offset"    # Solution type

################################ Self calibration loop ################################
count=0
for si,mf,fl,ty in soln:
    count = count+1
    print "\n *** Self calibration loop ",count,"si=",si,ty
    # Edit window
    print "Display Dirty image for editing"
    ODisplay.PImage(disp, Dirty, err, window=window)
    OErr.printErrMsg(err, "Error editing CLEAN boxes")
    OTF.ResidCalInput["solInt"]  = si
    OTF.ResidCalInput["minFlux"] = fl
    OTF.ResidCalInput["Clip"]    = mf
    OTF.ResidCalInput["solType"] = ty
    OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
    #print "Display Clean image for editing"
    #ODisplay.PImage(disp, CleanTmp, err, window=window)
    OErr.printErrMsg(err, "Error editing CLEAN boxes")
    # Cleanup table
    #print "DEBUG don't zap Soln"
    inOTF.ZapTable("OTFSoln",4,err)

print 'Finished with loop, final image'

################################## Final image/CLEAN ##################################
# Final image
image = OTF.makeImage(err,  OTF.ImageInput)
OErr.printErrMsg(err, "Error in final dirty image")
# Edit window
print "Display Final Dirty image for editing"
ODisplay.PImage(disp, Dirty, err, window=window)
OErr.printErrMsg(err, "Error editing CLEAN boxes")

# Final Clean
resid = CleanObj.Clean                       # Copy image just produced
Image.PCopy(image, resid, err)               # to clean(resid)        
CleanOTF.CleanInput["noResid"]  = False      # include residuals in calibration
CleanOTF.PClean(err,  CleanOTF.CleanInput)
OErr.printErrMsg(err, "Error Cleaning")

############################### Write calibrated output ###############################
print "Write calibrated data "
outData = OTF.newPOTF("Output data", outFile, inDisk, False, err)
OTF.PClone(inOTF, outData, err)     # Same structure etc
OErr.printErrMsg(err, "Error initializing output")

# Set time range
inInfo.set("timeRange", timerange)
inInfo.set("doCalSelect", True)
flagver=1
inInfo.set("flagVer", flagver)
gainuse=0
# DEBUGgainuse=4 # DEBUG
inInfo.set("gainUse", gainuse)
inInfo.set("doCalib", 1)
 
# Copy/calibrate
OTF.PCopy(inOTF, outData,  err)
OErr.printErrMsg(err, "Error selecting data")


##################################### History #########################################
print "Copy history"
tmpImage = Image.newPFImage("temp output",  tmpFile, inDisk, True, err)
# Do history to scratch image as table
inHistory  = History.History("history", inOTF.List, err)
outHistory = History.History("history", tmpImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" inFile = "+inFile,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" target = "+str(target),err)
outHistory.WriteRec(-1,ObitSys.pgmName+" scans = "+str(scans),err)
outHistory.WriteRec(-1,ObitSys.pgmName+" timerange = "+str(timerange),err)
outHistory.Close(err)
OErr.printErrMsg(err, "Error with history")

############################ Write output, Quantized image ############################
# output image
outImage  = Image.newPFImage("Output image", cleanFile,  inDisk,  False, err)
Image.PClone(tmpImage, outImage, err)   # Same structure etc.

# Copy to quantized integer image with history
print "Write output image"
inHistory  = History.History("history", tmpImage.List, err)
Image.PCopyQuantizeFITS (Dirty, outImage, err, \
                         inHistory=inHistory, fract=0.05)
# Copy CC tables
Image.PCopyTables(CleanTmp, outImage, [], ["AIPS CC"], err)
OErr.printErrMsg(err, "Error copying CC tables")
# Compress CC tables
cctab = outImage.NewTable(Table.READWRITE, "AIPS CC",1,err)
TableUtil.PCCMerge(cctab, cctab, err)
OErr.printErrMsg(err, "Error merging CC table")
print "Display output image"
ODisplay.PImage(disp, outImage, err)
OErr.printErrMsg(err, "Error displaying output image")

# Cleanup
#tmpImage.Zap(err) # Scratch float
OErr.printErrMsg(err, "Error cleaning up files")

# Compare with master lie [rms diff, max abs diff, max. master]
masterImage  = Image.newPFImage("Master image",   masterFile,  inDisk,  True, err)
diff = Image.PCompare(outImage, masterImage, err);
print "Comparison, rel. max. residual",diff[1]/diff[0], " rel RMS residual",diff[2]/diff[0]

# Final error check
OErr.printErr(err)

