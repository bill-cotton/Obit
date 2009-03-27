# Template ObitTalk script for processing Mustang data
import OSystem, OErr, InfoList, Image, Table, TableUtil, History, ODisplay
import OTF, OTFUtil, CleanOTF, OTFGetSoln, OTF, OTFGetAtmCor
import GBTUtil, FITSDir
import PARCal

# Init Obit
err=OErr.OErr()
FITS = ["./FITSdata"]  # Where to put output/scratch files
ObitSys=OSystem.OSystem  ("Mustang", 1, 1, 1, ["None"], len(FITS), FITS, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Allow multiple threads
OSystem.PAllowThreads(2)  # 2 threads
#######################################################################################
# Define parameters

# Root of data directory
DataRoot="/home/gbtdata/TPAR21/"
DataRoot = None  # To suppress attempt to update from archive

# Target to image and range of scans
target=["casa"]
scans = [77,104]
target=["A1835"]   # Cluster
scans = [20,1000]
target=["CRL2688"]   # calibrator
scans = [65,65]      # In focus OOF
scans = [72,72]      # In focus OOF

# Define data
# OTF file
inFile   = "AGBT09A_052_01OTF.fits" # OTF data file
#           ^^^^^^^^^^^^^^^^^^^^^
inDisk   = 0                        # 0 means current working directory
outDisk  = 1                        # Where resultant files will go (./FITSdata)

# Use target name to define output files
session= "52_01"
#         ^^^^^
outFile  = target[0]+session+"OTFCal.fits"        # Output calibrated data
dirtFile = "!"+target[0]+session+"Dirty.fits"     # Final output dirty image
cleanFile= "!"+target[0]+session+"Clean.fits"     # Final CLEAN output image
wtFile   = "!"+target[0]+session+"Wt.fits"        # Weighting image
BeamFile = "PARGaussBeam.fits"                    # Dirty Gaussian beam
priorFile= target[0]+session+"Clean.fits"         # Prior model?
BeamSize = 8.0                            # Beam size in asec
doScale = True                            # Do scaling by beam ratios

# List of feeds to use, empty list = all
feeds=[]
#feeds=[16]

# Default Image Info
timerange=[0.0,10.0]               # time range in days, All times
cells = 2.0                       # cell spacing in asec
nx = 256                          # no. cells in x
ny = 256                          # no. cells in y
niter   = 1000                    # Number of iteration of CLEAN
gain    = 0.1                     # CLEAN loop gain
minFlux = 0.001                   # Minimum image brightness to CLEAN
minResidFlux = 0.0                # Minimum Model brightness to use
maxResidFlux = 1.0                # Maximum Model brightness to use
minWt   = 0.001                   # Minimum weight in imaging wrt maximum
CLEANbox=[[-1,10, nx/2+1,ny/2+1]] # Clean window, circular at center
convType = 4                      # Gridding fn  Exp*sinc
convParm = [5.0, 1.65, 2.52, 2.0, 0.,0.,  0.,0., 0.,0.] # Gridding fn  params

# Default Calibration info
CalJy    = [38.5]                  # Cal values in Jy, one for all or per detector
solInt    = 1.0                    # Min solint in seconds
BLInt     = 30.0                   # Baseline filter time in sec
AtmInt    = 20.0                   # Atmospheric filter time in sec
CommonInt = None                   # Common mode filtering? Time is sec if desired.
PriorInt  = None                   # Common mode filtering? Time is sec if desired.
tau0      = 0.1                    # Zenith opacity
AtmEm     = 290.0                  # Zenith atmosphetic temperature eq. in Jy

# Default editing info
flagver   = 1                   # Flag table
flagInt   = 1.0                 # Flagging interval (sec), <=0 -> no flagging
maxRMS    = 0.4                 # Maximum allowable detector residual RMS in Jy
maxRatio  = 6.0                 # Max. allowable ratio to equivalent model RMS

# Table of pointing offsets in time order [time(day) d Xel (asec), d el (asec)]
PointTab=[\
    [0.0, 0.0, 0.0], \
    [1.0, 0.0, 0.0]]
#    *********** SET THIS ********
#PointTab=[
#    [0.070726878941059113, -4.0502607253392489, 0.96176494510748145] ,\
#    [0.076489381492137909, -1.6885605689484702, 1.6298337194204839] ,\
#    [0.14786446094512939, -1.2084478759180266, 1.3333042846993715] ,\
#    [0.19242700934410095, -2.0431437801421537, 0.38002386797447746] ,\
#    ]

# Table of opacities in time order [time(day), zenith opacity(nepers)]
tau0 = [[0.0000, 0.100],       \
        [1.0000, 0.100]]

# set of time scales for iterations
soln = [(6*solInt), (3*solInt), (solInt), (solInt)]
doOffset = True    # Do Offset cal before each MultiBeam cal?
deMode   = True    # Subtract the mode of the image when forming

# The following resets parameters for particular objects
# Cluster
if target[0]=="A1835":
    BLInt  = 20.0                   # Baseline filter time in sec
    AtmInt = 20.0                   # Atmospheric filter time in sec
    CommonInt = 6.                  # Common mode
    nx = 400; ny=400
    solInt  = 0.5
    niter   = 50000
    gain    = 0.3
    minFlux = 0.0001                   # Minimum image brightness to CLEAN
    minWt   = 0.001                   # Minimum weight in imaging wrt maximum
    soln = [4*solInt, 3*solInt, 2*solInt, solInt]
    CLEANbox=[[-1,55,200,215]]
    minResidFlux = -0.015           # Minimum Model brightness to use
    maxResidFlux = 0.0002           # Maximum Model brightness to use
    convType = 3                    # Gridding fn  Gaussian
    convParm = [0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0] # Gridding fn  params
# Calibrator
if target[0]=="CRL2688":
    solInt  = 0.5
    convType = 3                    # Gridding fn  Gaussian
    convParm = [0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0] # Gridding fn  params
    CLEANbox=[[-1,6,128,130]]
# Cas A
if target[0]=="casa":
    minWt   = 0.001                 # Minimum weight in imaging wrt maximum
    BLInt  = 30.0                   # Baseline filter time in sec
    AtmInt = 20.0                   # Atmospheric filter time in sec
    nx = 400; ny=400
    solInt = 1.0
    solInt = 0.5
    soln = [12*solInt, 6*solInt, 3*solInt, solInt]
    CLEANbox=[[-1,108,137,151]]
    niter = 50000
    gain = 0.1
    maxRMS    = 0.3                 # Maximum allowable detector residual RMS in Jy
    maxRatio  = 5.0                 # Max. allowable ratio to equivalent model RMS
# Crab 
if target[0]=="crab":
    nx = 500; ny=500
    solInt = 1.0
    niter=50000
    gain=0.03
    soln = [6*solInt, 3*solInt, solInt, solInt]
    CLEANbox=[[-1,101,252,253]]
    minWt   = 0.00001  # Minimum weight in imaging wrt maximum
# 2253+1608
if target[0]=="2253+1608":
    nx = 200; ny=200
    solInt = 1.0
    niter=100
    gain=0.1
    soln = [3*solInt, solInt]
    CLEANbox=[[-1,3,95,102]]
    flagInt   = -1.0                 # Flagging interval (sec), <=0 -> no flagging
    minWt   = 0.001  # Minimum weight in imaging wrt maximum

################################# Update data #########################################
# Update OTF with any new scans in archive, 50 msec averaging
inOTF = GBTUtil.UpdateOTF ("PAROTF","Rcvr_PAR",inFile, inDisk, DataRoot, err, \
                           avgTime=0.05)
OErr.printErrMsg(err, "Error creating/updating input data object")
inInfo = inOTF.List

print "Processing", target, "scans", scans

########################### If prior model given #############################
# dirty beam
dirtyBeam  = Image.newPFImage("dirty beam", BeamFile,  inDisk,  True, err)
OErr.printErrMsg(err, "Error initializing dirty beam")

if priorFile and FITSDir.PExist(priorFile, inDisk, err):
    print "Using prior model in",priorFile
    prior  = Image.newPFImage("prior model", priorFile,  inDisk,  True, err)
    PSF    = dirtyBeam
else:
    prior = None
    PSF   = None

################################## Initial calibration #############################
PARCal.InitCal(inOTF, target, err,\
               flagver=flagver, CalJy=CalJy, BLInt=BLInt, AtmInt=AtmInt,tau0=tau0, \
               prior=prior, PSF=PSF, PointTab=PointTab)

############################### Write calibrated output ###############################
# Apply current calibration and use result for remaining calibration
print "Write calibrated data "
# Delete output if it exists */
if FITSDir.PExist (outFile, outDisk, err):
    zapOTF = OTF.newPOTF("Output data", outFile, outDisk, True, err)
    zapOTF.Zap(err)
outOTF = OTF.newPOTF("Output data", outFile, outDisk, False, err)
OTF.PClone(inOTF, outOTF, err)     # Same structure etc
OErr.printErrMsg(err, "Error initializing output")

# Set data selection
inInfo.set("Targets",target)
inInfo.set("Scans",scans)
inInfo.set("timeRange", timerange)
inInfo.set("doCalSelect", True)
inInfo.set("flagVer", flagver)
gainuse=0
inInfo.set("gainUse", gainuse)
inInfo.set("doCalib", 1)
 
# Copy/calibrate
OTF.PCopy(inOTF, outOTF,  err)
OTFUtil.PIndex(outOTF, err)
OErr.printErrMsg(err, "Error selecting data")
# Create an initial dummy table with a interval 1/4 of the shortest
# Filter type solution interval.
inter = solInt/4
OTFGetSoln.POTFGetDummyCal (outOTF, outOTF, inter, 1, 1, err)
inInfo = outOTF.List

############################ Initial Common Mode filter ########################
baseCal = 1
if (CommonInt):
    print "Common mode only calibration, si=",CommonInt
    ResidCalInput = OTF.ResidCalInput
    ResidCalInput["InData"]  = outOTF
    ResidCalInput["solInt"]  = CommonInt
    ResidCalInput["solType"] = "MultiBeam" # Common mode
    ResidCalInput["gainUse"] = 1
    ResidCalInput["flagVer"] = flagver     # Which flag table to apply, -1 = none  
    ResidCalInput["minFlux"] = 100.0       # Clipping level
    ResidCalInput["Model"]   = None
    OTF.ResidCal(err, input=ResidCalInput)
    OErr.printErrMsg(err, "Error in initial commonMode")
    # Apply calibration
    Soln2CalInput = OTF.Soln2CalInput
    Soln2CalInput["InData"]  = outOTF      # Input data object
    Soln2CalInput["oldCal"]  = 1           # Input cal table
    Soln2CalInput["newCal"]  = 2           # New cal table 
    Soln2CalInput["soln"]    = 0           # highest
    OTF.Soln2Cal(err, input=Soln2CalInput)
    OErr.printErrMsg(err, "Error in initial commonMode calibration")
    baseCal = 2
    # End common mode

############################ Initial prior residual filter ########################
if (priorInt):
    print "Initial prior calibration, si=",priorInt
    gainuse = 0
    inInfo.set("gainUse", gainuse)
    inInfo.set("doCalib", 1)
    inInfo.set("flagVer", flagver)
    if priorModel:
        # Use CC table
        resid = OTFUtil.PSubModel(outOTF, None, prior, PSF, err)
    else:
        # Use Image
        prior.Open(Image.READONLY,err)
        prior.Read(err)
        prior.Close(err)
        resid = OTF.PScratch (outOTF, err)
        OTFUtil.PSubImage(outOTF, resid, prior.FArray, prior.Desc, err)
        OErr.printErrMsg(err, "Error with residial data")
        
    # residual calibration
    scrInfo = resid.List
    scrInfo.set("solInt", priorInt/86400.0)
    scrInfo.set("calJy", CalJy)
    scrInfo.set("calType", "MultiBeam")   # Common mode
    solnTable = Obit.OTFGetSolnCal (resid.me, outOTF.me, err.me)
    #scrInfo.set("calType", "Offset")      # Detector offsets
    #solnTable = Obit.OTFGetSolnGain (resid.me, outOTF.me, err.me)
    #del resid
    OErr.printErrMsg(err, "Error in initial commonMode")
    # Apply calibration
    Soln2CalInput = OTF.Soln2CalInput
    Soln2CalInput["InData"]  = outOTF      # Input data object
    Soln2CalInput["oldCal"]  = baseCal     # Input cal table
    Soln2CalInput["newCal"]  = baseCal+1   # New cal table 
    Soln2CalInput["soln"]    = 0           # highest
    OTF.Soln2Cal(err, input=Soln2CalInput)
    OErr.printErrMsg(err, "Error in initial commonMode calibration")
    baseCal += 1
    # End prior cal

################################## Set parameters ##############################
# Get position from OTF
pos =  GBTUtil.GetTargetPos(inOTF, target[0], err)
ra  = pos[0]                      # ra of center
dec = pos[1]                      # dec of center

# Imaging parameters
OTF.ImageInput["InData"]   = outOTF
OTF.ImageInput["disk"]     = outDisk
OTF.ImageInput["OutName"]  = dirtFile
OTF.ImageInput["Beam"]     = dirtyBeam
OTF.ImageInput["OutWeight"]= wtFile
OTF.ImageInput["ra"]       = ra              # Center RA
OTF.ImageInput["dec"]      = dec             # Center Dec
OTF.ImageInput["xCells"]   = cells           # "X" cell spacing
OTF.ImageInput["yCells"]   = cells           # "Y" cell spacing
OTF.ImageInput["nx"]       = nx              # number of cells in X
OTF.ImageInput["ny"]       = ny              # number of cells in Y
OTF.ImageInput["gainUse"]  = 0               # Which cal table to apply, -1 = none
OTF.ImageInput["flagVer"]  = flagver         # Which flag table to apply, -1 = none
OTF.ImageInput["minWt"   ] = minWt           # Minimum weight in imaging - includes data weight
OTF.ImageInput["ConvType"] = convType        # Gridding fn 0= pillbox, 3=Gaussian, 4=Exp*Sinc
OTF.ImageInput["ConvParm"] = convParm        # Gridding fn parameters

# Beam scaling?
jnInfo = outOTF.List
if not doScale:
        jnInfo.set("beamNx",96)
        jnInfo.set("beamNy",96)
jnInfo.set("doScale",doScale)

# Calibration parameters (some reset in loop)
OTF.ResidCalInput["InData"]  = outOTF       # Input data object
OTF.ResidCalInput["solType"] = "MultiBeam"  # Solution type
OTF.ResidCalInput["solInt"]  = 500.0        # Solution interval (sec)
OTF.ResidCalInput["minFlux"] = minResidFlux # Minimum Model brightness to use      
OTF.ResidCalInput["maxFlux"] = maxResidFlux # Maximum Model brightness to use      
OTF.ResidCalInput["Clip"]    = 1000.0       # Minimum image brightness to use in model
OTF.ResidCalInput["gainUse"] = baseCal      # Prior calibration, 0-> highest
OTF.ResidCalInput["minEl"]   = -90.0        # minimum elevation
OTF.ResidCalInput["flagVer"] = flagver      # Which flag table to apply, -1 = none

############################## Create initial image ###################################
print "Make Dirty Image"

OTF.ImageInput["gainUse"] = 0               # Which cal table to
if len(feeds)>0:
    inInfo.set("Feeds", feeds)              # Which feeds
DirtyImg = OTF.makeImage(err, OTF.ImageInput)
OErr.printErrMsg(err, "Error making initial image")

################################# CLEAN parameters ####################################
# Image display
disp = ODisplay.ODisplay("ObitView", "ObitView", err)

# Image for cleaning
CleanImg = Image.newPFImage("Clean Image", cleanFile, outDisk, False, err)

# Create CleanOTF
CleanObj = CleanOTF.PCreate("Clean", DirtyImg, dirtyBeam, CleanImg, err)
OErr.printErrMsg(err, "Error creating CLEAN object")

# CLEAN parameters
CleanOTF.CleanInput["CleanOTF"] = CleanObj     # Clean object
CleanOTF.CleanInput["disp"]     = disp         # Image display object
CleanOTF.CleanInput["Patch"]    = 40           # Beam patch
CleanOTF.CleanInput["Niter"]    = niter        # Number of iterations
CleanOTF.CleanInput["Gain"]     = gain         # CLEAN loop gain
CleanOTF.CleanInput["BeamSize"] = BeamSize/3600.0   # CLEAN restoring beam size in deg
CleanOTF.CleanInput["minFlux"]  = minFlux      # Minimum image brightness to CLEAN
CleanOTF.CleanInput["CCVer"]    = 1            # Clean components table version
CleanOTF.CleanInput["Feeds"]    = feeds        # list of feeds to use
CleanOTF.CleanInput["noResid"]  = True         # Don't include residuals in calibration

#  CLEAN window
for win in CLEANbox:
    CleanOTF.PAddWindow(CleanObj, win, err)

# Reset Soln2Cal parameters for self cal
OTF.Soln2CalInput["InData"]  = outOTF      # Input data object
OTF.Soln2CalInput["oldCal"]  = baseCal     # Input cal table
OTF.Soln2CalInput["newCal"]  = baseCal+1   # New cal table

# Beam scaling?
CleanOTF.CleanInput["scale"]=doScale

################################ Self calibration loop ################################
# For each iteration specify the solution interval in soln
# Depending on doOffset, do pair of calibrations, the first with
#  solType="Offset" with interval 3*si followed by a "MultiBeam"
#  solution with interval si or simply 
# a "MultiBeam" solution with interval si

count=0
OTF.ImageInput["gainUse"]    = baseCal                       # Which cal table to apply
OTF.ResidCalInput["gainUse"] = OTF.ImageInput["gainUse"]     # Prior calibration,
inInfo.set("deMode", False)  # Don't remove mode first cycle
if doOffset:
    for si in soln:
        count = count+1
        print "\n *** Self calibration loop ",count,"si=",3*si,"Offset"
        # First calibration of pair
        OTF.ResidCalInput["solInt"]  = 3*si
        OTF.ResidCalInput["solType"] = "Offset"
        OTF.Soln2CalInput["oldCal"]  = baseCal
        OTF.Soln2CalInput["newCal"]  = baseCal+1
        OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
        OErr.printErrMsg(err, "Error in self cal")
        OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
        OTF.ResidCalInput["gainUse"] = OTF.Soln2CalInput["newCal"]   # Prior calibration
        # Second
        print "\n *** second calibration of loop ",count,"si=",si,"MultiBeam"
        OTF.ResidCalInput["solInt"]  = si
        OTF.ResidCalInput["solType"] = "MultiBeam"
        OTF.Soln2CalInput["oldCal"]  = baseCal+1
        OTF.Soln2CalInput["newCal"]  = baseCal+2 
        OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
        OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
        OTF.ResidCalInput["gainUse"] = baseCal   # Prior calibration for next cycle
        # Cleanup Soln tables
        outOTF.ZapTable("OTFSoln",-1,err)
        inInfo.set("deMode", deMode) # remove mode
else:   # Only MultiBeam
    for si in soln:
        count = count+1
        print "\n *** calibration loop ",count,"si=",si,"MultiBeam"
        OTF.ResidCalInput["solInt"]  = si
        OTF.ResidCalInput["solType"] = "MultiBeam"
        OTF.Soln2CalInput["oldCal"]  = baseCal
        OTF.Soln2CalInput["newCal"]  = baseCal+1
        OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
        OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
        OTF.ResidCalInput["gainUse"] = baseCal  # Prior calibration for next cycle
        # Cleanup Soln tables
        outOTF.ZapTable("OTFSoln",-1,err)
        inInfo.set("deMode", deMode) # remove mode

print 'Finished with loop, final image'

################################## Editing data #############################
if flagInt>0.0:
    # Editing
    print " Editing data"
    if flagver<=0:
        flagver = 1
    inInfo = outOTF.List
    inInfo.set("flagInt", flagInt/86400.)
    inInfo.set("maxRMS", maxRMS)
    inInfo.set("maxRatio", maxRatio)
    inInfo.set("flagVer", flagver)
    inInfo.set("FGVer", flagver)
    inInfo.set("doCalSelect", True)
    inInfo.set("doCalib", 1)
    inInfo.set("gainUse", 0)
    OTFGetSoln.PFlag(outOTF, CleanImg, outOTF, 1, err)
    OErr.printErr(err)

################################## Final image/CLEAN ##################################
# Final image
inInfo = outOTF.List
inInfo.set("doCalib", 1)
inInfo.set("gainUse", 0)
DirtyImg = OTF.makeImage(err,  OTF.ImageInput)
OErr.printErrMsg(err, "Error in final dirty image")

info = CleanObj.List
info.set("doRestore", True)
CleanOTF.CleanInput["doRestore"]= True            # Do restore CCs
CleanOTF.CleanInput["BeamSize"] = BeamSize/3600.0 # CLEAN restoring beam size in deg

# Final Clean
CleanOTF.CleanInput["noResid"]  = False      # include residuals
CleanOTF.CleanInput["doScale"]  = doScale    # Scale residuals
resid = CleanObj.Clean                       # Copy image just produced
Image.PCopy(DirtyImg, resid, err)            # to clean (residual)        
CleanOTF.PClean(err,  CleanOTF.CleanInput)
OErr.printErrMsg(err, "Error Cleaning")

# Compress CC tables
cctab = CleanImg.NewTable(Table.READWRITE, "AIPS CC", 1, err)
TableUtil.PCCMerge(cctab, cctab, err)
OErr.printErrMsg(err, "Error merging CC table")

##################################### History #########################################
print "Copy history"
# Loop over dirty, clean images, output data
for img in [DirtyImg, CleanImg, outOTF]:
    inHistory  = History.History("in history", inOTF.List, err)
    outHistory = History.History("out history", img.List, err)
    History.PCopy(inHistory, outHistory, err)
    
    # Add this programs history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" target = "+str(target),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" scans = "+str(scans),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" timerange = "+str(timerange),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" doScale = "+str(doScale),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" BeamFile = "+BeamFile,err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" BeamSize = "+str(BeamSize),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" doOffset = "+str(doOffset),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" deMode = "+str(deMode),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" minWt = "+str(minWt),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" convType = "+str(convType),err)
    line = ObitSys.pgmName+"convParm = [%f,%f,%f,%f]"%(convParm[0],convParm[1],convParm[2],convParm[3] )
    outHistory.WriteRec(-1,line,err)
    # Clean
    outHistory.WriteRec(-1,ObitSys.pgmName+" niter = "+str(niter),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" gain = "+str(gain),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" minFlux = "+str(minFlux),err)
    i = 0
    for box in CLEANbox:
        i += 1
        line = ObitSys.pgmName+"box[%d] = [%d,%d,%d,%d]"%(i, box[0],box[1],box[2],box[3])
        outHistory.WriteRec(-1,line,err)
        
    # Calibration
    outHistory.WriteRec(-1,ObitSys.pgmName+" solInt = "+str(solInt),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" minResidFlux = "+str(minResidFlux),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" maxResidFlux = "+str(maxResidFlux),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" BLInt = "+str(BLInt),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" AtmInt = "+str(AtmInt),err)
    if priorFile and FITSDir.PExist(priorFile, inDisk, err):
        outHistory.WriteRec(-1,ObitSys.pgmName+" priorFile = "+priorFile,err)
    if priorInt:
        outHistory.WriteRec(-1,ObitSys.pgmName+" priorInt = "+str(priorInt),err)
    if CommonInt:
        outHistory.WriteRec(-1,ObitSys.pgmName+" CommonInt(CM) = "+str(CommonInt),err)
    i = 0
    for si in soln:
        i += 1
        line = ObitSys.pgmName+"Solution interval[%d] = %f sec"%(i, si)
        outHistory.WriteRec(-1,line,err)
    # Editing
    outHistory.WriteRec(-1,ObitSys.pgmName+" flagver = "+str(flagver),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" flagInt = "+str(flagInt)+" \ <=0 -> no editing",err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" maxRMS = "+str(maxRMS),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" maxRatio = "+str(maxRatio),err)
    outHistory.Close(err)
    # Copy history to header
    inHistory  = History.History("in history",  img.List, err)
    outHistory = History.History("out history", img.List, err)
    History.PCopy2Header(inHistory, outHistory, err)
    OErr.printErrMsg(err, "Error with history")


###################### Display final CLEAN image, shutdown #######################
print "Display final CLEAN image"
ODisplay.PImage(disp, CleanImg, err)
OErr.printErrMsg(err, "Error displaying output image")

# Shutdown Obit
OErr.printErr(err)
del ObitSys
