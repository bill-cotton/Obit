# Template ObitTalk script for processing Mustang data
import OSystem, OErr, InfoList, Image, Table, TableUtil, History, ODisplay
import OTF, OTFUtil, CleanOTF, OTFGetSoln, OTF, OTFGetAtmCor
import GBTUtil, FITSDir
import PARCal
import math

# Init Obit
err=OErr.OErr()
FITS = ["./FITSdata"]  # Where to put putput/scratch files
ObitSys=OSystem.OSystem  ("Mustang", 1, 1, 1, ["None"], len(FITS), FITS, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Allow multiple threads
#OSystem.PAllowThreads(2)  # 2 threads
OSystem.PAllowThreads(1)  # 1 thread

#######################################################################################
# Define parameters

# Root of data directory
DataRoot="/home/gbtdata/AGBT08A_056_01/"
DataRoot = None  # To suppress attempt to update from archive

# Target to image and range of scans
target = ["MARS"]
scans = [41,48]                 # range of scans SET THIS
target = ["m82"]
scans = [109,112]              # range of scans SET THIS
target=["eskimo"]
scans = [55,72]                # range of scans SET THIS
target = ["m87"]
scans = [93,104]               # range of scans SET THIS
target=["bolorionN2","bolorionNbars","bolorionS"]
scans = [36,56]                 # range of scans SET THIS
target=["0530+135"]
scans = [27,28]

# Define data
# OTF file
inFile   = "AGBT08A_056_01OTF.fits" # OTF data file
inDisk   = 0                        # 0 means current working directory
outDisk  = 1                        # Where resultant files will go (./FITSdata)

# Use target name to define output files
outFile  = target[0]+"OTFCal.fits"        # Output calibrated data
dirtFile = "!"+target[0]+"Dirty.fits"     # Final output dirty image
cleanFile= "!"+target[0]+"Clean.fits"     # Final CLEAN output image
wtFile   = "!"+target[0]+"Wt.fits"        # Weighting image
BeamFile = "PAR1DBeam9.fits"              # Dirty beam from 3C279, 9 FWHM asec core
BeamSize = 9.0                            # Beam size in asec
doScale = False                           # No scaling by beam ratios
BeamFile = "PARGaussBeam.fits"            # Dirty Gaussian beam
BeamSize = 8.0                            # Beam size in asec
doScale = True                            # Do scaling by beam ratios
priorFile= target[0]+"Clean.fits"         # Prior model?

# List of feeds to use, empty list = all
feeds=[]

# Default Image Info
timerange=[0.0,1.0]               # time range in days, All times
cells = 2.0                       # cell spacing in asec
nx = 256                          # no. cells in x
ny = 256                          # no. cells in y
niter   = 1000                    # Number of iteration of CLEAN
gain    = 0.1                     # CLEAN loop gain
minFlux = 0.001                   # Minimum image brightness to CLEAN
minWt   = 0.00001                 # Minimum weight in imaging wrt maximum
CLEANbox=[[-1,10, nx/2+1,ny/2+1]] # Clean window, circular at center
convType = 4                      # Gridding fn  Exp*sinc
convParm = [0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0] # Gridding fn  params
pointOff= [0.0,0.0]               # RA, Dec Offset in cells to add to position
                                  # This is the uncorrected part of the pointing error

# Default Calibration info
CalJy = [38.5]                  # Cal values in Jy, one for all or per detector
solInt = 1.0                    # Min solint in seconds
BLInt  = 30.0                   # Baseline filter time in sec
AtmInt = 20.0                   # Atmospheric filter time in sec
tau0   = 0.1                    # Zenith opacity
AtmEm  = 290.0                  # Zenith atmosphetic temperature eq. in Jy

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
#    *********** SET THIS ********
# From Ron's weather server from CLEO 
#tau0 = [[0.0, 0.102453241896],        \
#        [0.0416667, 0.100114034451],  \
#        [0.0833333, 0.0999824523373], \
#        [0.125, 0.106008646735],      \
#        [0.1666667, 0.113737894416],  \
#        [0.2083333, 0.106998485604],  \
#        [0.25, 0.154141027655]] 


# set of time scales for iterations
soln = [(6*solInt), (3*solInt), (solInt), (solInt)]
doOffset = True    # Do Offset cal before each MultiBeam cal?
deMode   = True    # Subtract the mode of the image when forming

# The following resets parameters for particular objects
# Orion 
if target[0]=="bolorionN":
    minWt   = 1.0e-10               # Minimum weight in imaging wrt maximum
    BLInt  = 30.0                   # Baseline filter time in sec
    AtmInt = 20.0                   # Atmospheric filter time in sec
    nx = 500; ny=500
    solInt = 1.0
    soln = [6*solInt, 3*solInt, solInt, solInt]
    CLEANbox=[[177,49,332,333]]
    niter = 50000
    gain = 0.05
    # Crab 
if target[0]=="crab":
    nx = 500; ny=500
    solInt = 1.0
    niter=50000
    gain=0.03
    soln = [6*solInt, 3*solInt, solInt, solInt]
    CLEANbox=[[-1,101,252,253]]
    minWt   = 0.00001  # Minimum weight in imaging wrt maximum
# M82
if target[0]=="m82":
    inFile   = "M82RawOTF.fits"       # M82 data
    solInt = 1.0
    soln = [5*solInt, 3*solInt, 2*solInt, solInt]
    CLEANbox=[[119, 124, 147, 141]]
# Eskimo
if target[0]=="eskimo":
    inFile   = "EskimoRawOTF.fits"    # Eskimo nebula + cal only
    solInt = 2.0
    soln = [3*solInt, 2*solInt, solInt]
    CLEANbox=[[-1,16,129,126]]
    niter = 500
    doOffset = False
# M87
if target[0]=="m87":
    inFile   = "M87RawOTF.fits"       # M87
    niter = 5000
    solInt = 1.0
    soln = [(5*solInt),3*(solInt), (solInt), (solInt)]
    CLEANbox=[[-1,11,136,131],[-1,15,120,124],[-1,8,130,120],[-1,13,102,125]]

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
# Set number of records depending on number of threads
nrec = 5000*OSystem.PGetNoThreads()
outOTF = OTF.newPOTF("Output data", outFile, outDisk, False, err, nrec=nrec)
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
inInfo.set("nRecPIO", 1000)   # Size of read in records
 
# Copy/calibrate
OTF.PCopy(inOTF, outOTF,  err)
OTFUtil.PIndex(outOTF, err)
OErr.printErrMsg(err, "Error selecting data")
# Create an initial dummy table with a interval 1/4 of the shortest
# Filter type solution interval.
inter = solInt/4
OTFGetSoln.POTFGetDummyCal (outOTF, outOTF, inter, 1, 1, err)
inInfo = outOTF.List

################################## Atmospheric emission ##############################
print "Atmospheric emission"
solint = solInt/86400.0
inInfo.set("Tau0", 0.2)                   # Opacity table
inInfo.set("AtmEm", AtmEm)                # Zenith emission in Jy
inInfo.set("solInt", solint)
inInfo.set("doCalSelect", True)
inInfo.set("flagVer", flagver)
gainuse = 0
inInfo.set("gainUse", gainuse)
inInfo.set("doCalib", 1)
inInfo.set("doFilter", True)
OTFGetAtmCor.PAtmEm(outOTF, outOTF, err)
OErr.printErrMsg(err, "Error with Baseline calibration")

# Soln2Cal parameters for filter cal (most defaulted)
OTF.Soln2CalInput["InData"]  = outOTF      # Input data object
OTF.Soln2CalInput["oldCal"]  = 1           # Use gain cal output
OTF.Soln2CalInput["newCal"]  = 2           # New cal table
OTF.Soln2Cal(err, OTF.Soln2CalInput)  # Apply
OErr.printErrMsg(err, "Error updating Cal table with Soln")
    
################################## Set parameters ##############################
# Get position from OTF
pos =  GBTUtil.GetTargetPos(inOTF, target[0], err)
ra  = pos[0]                      # ra of center
dec = pos[1]                      # dec of center
# Fudge pointing if needed
ra  += math.cos(dec/57.296)*pointOff[0]*cells/3600.0
dec += pointOff[1]*cells/3600.0

# Imaging parameters
OTF.ImageInput["InData"]  = outOTF
OTF.ImageInput["disk"]    = outDisk
OTF.ImageInput["OutName"] = dirtFile
OTF.ImageInput["Beam"]    = dirtyBeam
OTF.ImageInput["OutWeight"]= wtFile
OTF.ImageInput["ra"]  = ra                   # Center RA
OTF.ImageInput["dec"] = dec                  # Center Dec
OTF.ImageInput["xCells"] = cells             # "X" cell spacing
OTF.ImageInput["yCells"] = cells             # "Y" cell spacing
OTF.ImageInput["nx"] = nx                    # number of cells in X
OTF.ImageInput["ny"] = ny                    # number of cells in Y
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
OTF.ResidCalInput["InData"]  = outOTF      # Input data object
OTF.ResidCalInput["solType"] = "MultiBeam" # Solution type
OTF.ResidCalInput["solInt"]  = 500.0       # Solution interval (sec)
OTF.ResidCalInput["minFlux"] = 0.0         # Minimum Model brightness to use
OTF.ResidCalInput["Clip"]    = 1000.0      # Minimum image brightness to use in model
OTF.ResidCalInput["gainUse"] = 1           # Prior calibration, 0-> highest
OTF.ResidCalInput["minEl"]   = -90.0       # minimum elevation
OTF.ResidCalInput["flagVer"] = flagver     # Which flag table to apply, -1 = none

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
OTF.Soln2CalInput["oldCal"]  = 2           # Input cal table
OTF.Soln2CalInput["newCal"]  = 3           # New cal table

# Beam scaling?
CleanOTF.CleanInput["scale"]=doScale

################################ Self calibration loop ################################
# For each iteration specify the solution interval in soln
# Depending on doOffset, do pair of calibrations, the first with
#  solType="Offset" with interval 3*si followed by a "MultiBeam"
#  solution with interval si or simply 
# a "MultiBeam" solution with interval si

count=0
OTF.ImageInput["gainUse"]    = 2                             # Which cal table to apply
OTF.ResidCalInput["gainUse"] = OTF.ImageInput["gainUse"]     # Prior calibration,
inInfo.set("deMode", False)  # Don't remove mode first cycle
if doOffset:
    for si in soln:
        count = count+1
        print "\n *** Self calibration loop ",count,"si=",3*si,"Offset"
        # First calibration of pair
        OTF.ResidCalInput["solInt"]  = 3*si
        OTF.ResidCalInput["solType"] = "Offset"
        OTF.Soln2CalInput["oldCal"]  = 2
        OTF.Soln2CalInput["newCal"]  = 3 
        OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
        OErr.printErrMsg(err, "Error in self cal")
        OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
        OTF.ResidCalInput["gainUse"] = OTF.Soln2CalInput["newCal"]   # Prior calibration
        # Second
        print "\n *** second calibration of loop ",count,"si=",si,"MultiBeam"
        OTF.ResidCalInput["solInt"]  = si
        OTF.ResidCalInput["solType"] = "MultiBeam"
        OTF.Soln2CalInput["oldCal"]  = 3
        OTF.Soln2CalInput["newCal"]  = 4 
        OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
        OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
        OTF.ResidCalInput["gainUse"] = 2   # Prior calibration for next cycle
        # Cleanup Soln tables
        outOTF.ZapTable("OTFSoln",2,err)
        outOTF.ZapTable("OTFSoln",3,err)
        inInfo.set("deMode", deMode) # remove mode
else:   # Only MultiBeam
    for si in soln:
        count = count+1
        print "\n *** calibration loop ",count,"si=",si,"MultiBeam"
        OTF.ResidCalInput["solInt"]  = si
        OTF.ResidCalInput["solType"] = "MultiBeam"
        OTF.Soln2CalInput["oldCal"]  = 2
        OTF.Soln2CalInput["newCal"]  = 3
        OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
        OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
        OTF.ResidCalInput["gainUse"] = 2   # Prior calibration for next cycle
        # Cleanup Soln tables
        outOTF.ZapTable("OTFSoln",2,err)
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
DirtyImg = OTF.makeImage(err,  OTF.ImageInput)
OErr.printErrMsg(err, "Error in final dirty image")

# Final Clean
resid = CleanObj.Clean                       # Copy image just produced
Image.PCopy(DirtyImg, resid, err)            # to clean (residual)        
CleanOTF.CleanInput["noResid"]  = False      # Include residuals
CleanOTF.CleanInput["doScale"]  = doScale    # Scale residuals
CleanOTF.PClean(err,  CleanOTF.CleanInput)
OErr.printErrMsg(err, "Error Cleaning")

# Compress CC tables
cctab = CleanImg.NewTable(Table.READWRITE, "AIPS CC", 1, err)
TableUtil.PCCMerge(cctab, cctab, err)
OErr.printErrMsg(err, "Error merging CC table")

############################## Remove pointing offset #################################
for img in [DirtyImg, CleanImg]:
    head = img.Desc.Dict
    head["crval"][0] -= math.cos(dec/57.296)*pointOff[0]*cells/3600.0
    head["crval"][1] -= pointOff[1]*cells/3600.0
    img.Desc.Dict = head
    img.UpdateDesc(err)
    OErr.printErrMsg(err, "Error updating header")

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
    outHistory.WriteRec(-1,ObitSys.pgmName+" inFile = "+inFile,err)
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
    outHistory.WriteRec(-1,ObitSys.pgmName+" pointOff = "+str(pointOff),err)
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
    outHistory.WriteRec(-1,ObitSys.pgmName+" BLInt = "+str(BLInt),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+" AtmInt = "+str(AtmInt),err)
    outHistory.WriteRec(-1,ObitSys.pgmName+"  = "+str(),err)
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
