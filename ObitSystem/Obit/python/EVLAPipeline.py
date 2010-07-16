# Pipeline processing 
# AIPS/FITS setup and Parameter file given as arguments:
# > python EVLAPipeline.py AIPSSetup.py parms.py
#
import sys
import OErr, OSystem, UV, AIPS, FITS
import ObitTalkUtil
from AIPS import AIPSDisk
from FITS import FITSDisk
from EVLACal import *

############################# Initialize OBIT ##########################################                                 
setup = sys.argv[1]
noScrat     = []    
execfile (setup)

############################# Default parameters ##########################################                                 
# Generic parameters
seq           = 1          # AIPS sequence number
gain          = 0.10       # CLEAN loop gain
doLoad        = True       # Load data from FITS?, else already in AIPS 
doLoadArchive = False      # Load from archive?
dataClass     = "UVData"   # AIPS class of uv data
Compress      = False      # Use compressed UV data?
archRoot      = "."        # Archive directory root
selBand       = " "        # Selected band, def = first  
selChan       = 0          # Selected number of channels, def = first  
calInt        = 0.5        # Calibration table interval in min.
check         = False      # Only check script, don't execute tasks
debug         = False      # run tasks debug

# Editing
doClearTab  = True         # Clear cal/edit tables
doGain      = True         # Clear SN and CL tables >1
doFlag      = True         # Clear FG tables > 1
doBP        = True         # Clear BP tables?
doCopyFG    = True         # Copy FG 1 to FG 2quack
Reason      = "Quack"      # Reason code for Quack
doQuack     = True         # Quack data?
quackBegDrop= 0.1          # Time to drop from start of each scan in min
quackEndDrop= 0.1          # Time to drop from end of each scan in min
quackReason   = "Quack"                 # Reason string
doAutoFlag    = True       # Autoflag editing?
RMSAvg      = 20.0         # AutoFlag Max RMS/Avg for time domain RMS filtering
IClip       = [1000.0,0.1] # AutoFlag Stokes I clipping
VClip       = [10.0,0.05]  # AutoFlag Stokes V clipping
timeAvg     = 2.0          # AutoFlag time averaging in min.
doAFFD      = False        # do AutoFlag frequency domain flag
FDmaxAmp    = 0.0          # Maximum average amplitude (Jy)
FDmaxV      = 0.0          # Maximum average VPol amp (Jy)
FDwidMW     = 5            # Width of the median window
FDmaxRMS    = 0.0          # Channel RMS limits (Jy)
FDmaxRes    = 6.0          # Max. residual flux in sigma
FDmaxResBL  = 6.0          # Max. baseline residual
FDbaseSel   = [0,0,0,0]    # Channels for baseline fit
doMedn      = True         # Median editing?
timeWind    = 2.0          # Median window width in min for median flagging
avgTime     = 10.0/60.     # Averaging time in min
avgFreq     = 0            # 1=>avg chAvg chans, 2=>avg all chan, 3=> avg chan and IFs
chAvg       = 1            # number of channels to average

# Special editing list
doEditList  = False        # Edit using editList?
editFG      = 1            # Table to apply edit list to
editList = [
#    {"timer":("0/06:09:0.0","0/06:13:0.0"),"Ant":[ 8,0],"IFs":[2,2],"Chans":[1,0],"Stokes":'1110',"Reason":"bad data"},
    ]

# Bandpass Calibration
doBP        = True         # Clear BP tables
doBPCal     = True         # Determine Bandpass calibration
BPCal       = None         # Bandpass calibrator
bpBChan1      = 1          # Low freq. channel,  initial cal
bpEChan1      = 0          # Highest freq channel, initial cal, 0=>all
bpDoCenter1   = None       # Fraction of  channels in 1st, overrides bpBChan1, bpEChan1
bpBChan2      = 1          # Low freq. channel for BP cal
bpEChan2      = 0          # Highest freq channel for BP cal,  0=>all 
bpChWid2      = 1          # Number of channels in running mean BP soln

# Amp/phase calibration
PCal          = None                    # Phase calibrator
ACal          = None                    # Amplitude calibrator
solint        = 0.0                     # Calibration solution time
solsmo        = 0.0                     # Smooth calibration
ampScalar     = False                   # Amp-scalar operation in Calib
AcalModel     = "3C286_K.MODEL"         # Flux calibrator model file name
AcalDisk      = 1                       # Flux calibrator model FITS disk
refAnt      = 0            # Reference antenna

# Imaging
doImage     = True         # Image targets
targets     = []           # targets
outIclass   = "IClean"     # Output image class
do3D        = False        # 2/3D imaging
prtLev      = 2            # Imager print level
BLFact      = 1.01         # Baseline dependent time averaging factor
BLchAvg     = False        # If True, and BLFact >1.0 then also average channels.
Robust      = 0.0          # Weighting robust parameter
FOV         = 1.0          # Field of view radius in deg.
Niter       = 100          # Max number of clean iterations
minFlux     = 0.0          # Minimum CLEAN flux density
xCells      = 0.0          # x cell spacing in asec
yCells      = 0.0          # y cell spacing in asec
doAmpPhaseCal = True       # Amplitude/phase calibration
solPType    = "L1"         # L1 solution for phase self cal
solPMode    = "DELA"       # Delay solution for phase self cal
avgPol      = True         # Average poln in self cal?
solAType    = "    "       # Type solution for amp+phase self cal
avgIF       = False        # Average IF in self cal?
maxPSCLoop  = 0            # Max. number of phase self cal loops
minFluxPSC  = 0.050        # Min flux density peak for phase self cal
solPInt     = 0.25         # phase self cal solution interval (min)
maxASCLoop  = 0            # Max. number of phase self cal loops
minFluxASC  = 0.2          # Min flux density peak for amp+phase self cal
solAInt     = 3.0          # amp+phase self cal solution interval (min)
OutlierDist = 0.0          # Distance (deg) to which  to add outlying fields
OutlierFlux = 0.0          # Minimum outlier flux density in Jy

# Final
doSaveRes     = True       # Save results
doCleanup     = True       # Destroy AIPS files
prtLv         = 2          # Amount of diagnostics

############################# Set Project Processing parameters ##################
parmFile = sys.argv[2]
execfile (parmFile)

################################## Process #####################################
# Logging directly to logFile
OErr.PInit(err, prtLv, logFile)
retCode = 0

mess = "Start project "+project+" with UV data "+inFile+" disk "+str(inDisk)
printMess(mess, logFile)
if debug:
    mess = "Using Debug mode "
    printMess(mess, logFile)
if check:
    mess = "Only checking script"
    printMess(mess, logFile)

# Load Data from FITS
uv = None
if doLoad:
    uv = EVLAUVLoadT(inFile, inDisk, project, dataClass, disk, seq, err, logfile=logFile, \
                         Compress=True, check=check, debug=debug)
    if uv==None:
        raise RuntimeError,"Cannot load "+inFile
# Load Data from Archive directory
if doLoadArchive:
    uv = EVLAUVLoadArch(archRoot, project, dataClass, disk, seq, err, logfile=logFile, \
                            selBand=selBand, selChan=selChan, calInt=calInt, \
                            Compress=True, check=check, debug=debug)
    if uv==None:
        raise RuntimeError,"Cannot load "+inFile
# Otherwise set uv
if uv==None:
    uv = UV.newPAUV("AIPS UV DATA", project, dataClass, disk, seq, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")

# Clear any old calibration/editing 
if doClearTab:
    EVLAClearCal(uv, err, doGain=doGain, doFlag=doFlag, doBP=doBP)
    OErr.printErrMsg(err, "Error resetting calibration")

# Copy FG 1 to FG 2
if doCopyFG:
    retCode = EVLACopyFG (uv, err, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error Copying FG table"

# Special editing
if doEditList:
    for edt in editList:
        UV.PFlag(uv,err,timeRange=[dhms2day(edt["timer"][0]),dhms2day(edt["timer"][1])], \
                     flagVer=editFG, Ants=edt["Ant"], Chans=edt["Chans"], IFs=edt["IFs"], \
                     Stokes=edt["Stokes"], Reason=edt["Reason"])
        OErr.printErrMsg(err, "Error Flagging")

# Quack to remove data from start and end of each scan
if doQuack:
    retCode = EVLAQuack (uv, err, begDrop=quackBegDrop, endDrop=quackEndDrop, Reason=quackReason, \
                             logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error Quacking data"

# Median editing
if doMedn:
    retCode = EVLAMedianFlag (uv, "    ", err, noScrat=noScrat, nThreads=nThreads, \
                                  avgTime=avgTime, avgFreq=avgFreq,  chAvg= chAvg, \
                                  timeWind=timeWind, flagVer=2, logfile=logFile, \
                                  check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in MednFlag"

# Bandpass calibration if needed
if doBPCal and BPCal:
    retCode = EVLABPCal(uv, BPCal, err, noScrat=noScrat, solInt1=bpsolint2, solInt2=bpsolint2, \
                            BChan1=bpBChan1, EChan1=bpEChan1, BChan2=bpBChan2, EChan2=bpEChan2, ChWid2=bpChWid2, \
                            doCenter1=bpDoCenter1, refAnt=refAnt, specIndex=specIndex, logfile=logFile, \
                            check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in Bandpass calibration"

# Amp & phase Calibrate
if doAmpPhaseCal:
    retCode = EVLACal (uv, targets, ACal, err, PCal=PCal, doBand=1, BPVer=1, \
                           calModel=AcalModel, calDisk=AcalDisk, nThreads=nThreads, \
                           solInt=solint, solSmo=solsmo, noScrat=noScrat, ampScalar=ampScalar, \
                           refAnt=refAnt, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error calibrating"

# More editing
if doAutoFlag:
    retCode = EVLAAutoFlag (uv, targets, err, RMSAvg=RMSAvg, flagVer=2, \
                                doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                                IClip=IClip, VClip=VClip, timeAvg=timeAvg, \
                                doFD=doAFFD, FDmaxAmp=FDmaxAmp, FDmaxV=FDmaxV, FDwidMW=5, \
                                FDmaxRes=FDmaxRes,  FDmaxResBL= FDmaxResBL, FDbaseSel=FDbaseSel, \
                                logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in AutoFlag"

# Image targets
if doImage:
    img = EVLASetImager(uv, targets, outIclass=outIclass, nThreads=nThreads, \
                            logfile=logFile)
    img.xCells      = xCells       # x cell spacing in asec
    img.yCells      = yCells       # x cell spacing in asec
    img.outSeq      = seq          # output sequence number
    img.out2Seq     = seq          # output sequence number
    img.Niter       = Niter        # Maximum number of CLEAN components
    img.Gain        = gain         # CLEAN loop gain
    img.minFlux     = minFlux      # Minimum CLEAN flux density
    img.Robust      = Robust       # Weighting robust parameter 
    img.FOV         = FOV          # Field of view radius in deg
    img.maxPSCLoop  = maxPSCLoop   # Max. number of phase self cal loops
    img.minFluxPSC  = minFluxPSC   # Min flux density peak for phase self cal
    img.solPType    = solPType     # L1 solution for phase self cal
    img.solPMode    = solPMode     # Solution mode for phase self cal
    img.solPInt     = solPInt      # Phase self cal solution interval
    img.avgPol      = avgPol       # Average poln in self cal
    img.maxASCLoop  = maxASCLoop   # Max. number of Amp&phase self cal loops
    img.minFluxASC  = minFluxASC   # Min flux density peak for amp+phase self cal
    img.solAType    = solAType     # Type solution for amp+phase self cal
    img.solAInt     = solAInt      # amp+phase self cal solution interval
    img.avgIF       = avgIF        # Don't average IF in self cal
    img.do3D        = do3D         # 2/3D imagingDon't average IF in self cal
    img.OutlierDist = OutlierDist  # Distance (deg) to which  outlying fields
    img.OutlierFlux = OutlierFlux  # Minimum outlier apparent flux density
    img.dispURL     = "None"       # No display
    img.prtLv       = prtLev       # Control amount of messages
    img.BLFact      = BLFact       # Baseline dependent time averaging factor
    img.BLchAvg     = BLchAvg      # If True, and BLFact >1.0 then also average channels.
    # Bandpass calibration?
    if BPCal:
        img.doBand=1; img.BPVer = 1;
    #img.debug = True
    img.noScrat = noScrat
    if debug:
        img.i
        img.debug = debug
    # Trap failure
    try:
        if not check:
            img.g                      # Image/deconvolve
    except:
        mess = "Imaging failed "
        printMess(mess, logFile)
    else:
        pass
    mess = "Imager Return code = "+str(img.retCode)
    printMess(mess, logFile)

# cleanup
mess ="Write results/Cleanup" 
printMess(mess, logFile)
if doCleanup:
    uv.Zap(err)  # UV data

# Write results, cleanup    
# Save UV data 
if doSaveRes:
    filename = project+"Cal.uvtab"
    fuv = EVLAUVFITS (uv, filename, 0, err, compress=Compress)
# Delete UV data
if doCleanup:
    uv.Zap(err)
# Imaging results
for target in targets:
    if doSaveRes:
        x = Image.newPAImage("out", target, outIclass, disk, seq, True, err)
        outfile = target+"."+outIclass+".fits"
        mess ="Write " +outfile+" on disk "+str(outDisk)
        printMess(mess, logFile)
        xf = EVLAImFITS (x, outfile, outDisk, err)
        # Statistics
        imstat(x, err, logfile=logFile)
    if doCleanup:
        x.Zap(err) # cleanup
        # Zap Imager work file
        u = UV.newPAUV("out", target, "Imager", disk, seq, True, err)
        u.Zap(err) # cleanup
# end writing loop
OErr.printErrMsg(err, "Writing output/cleanup")

# Shutdown
mess = "Finished project "+project
printMess(mess, logFile)
OErr.printErr(err)
OSystem.Shutdown(ObitSys)

