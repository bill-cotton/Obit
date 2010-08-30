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
dataClass     = "UVData"   # AIPS class of raw uv data
Compress      = False      # Use compressed UV data?
archRoot      = "."        # Archive directory root
selBand       = " "        # Selected band, def = first  
selNIF        = 0          # Selected number of IFs, def = first  
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
mednSigma   = 10.0         # Median sigma clipping level
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
bpsolMode     = 'A&P'      # Band pass type 'A&P', 'P', 'P!A'
bpsolint1     = 10.0/60.0  # BPass phase correction solution in min
bpsolint2     = 10.0       # BPass bandpass solution in min

# Amp/phase calibration
PCal          = None                    # Phase calibrator
ACal          = None                    # Amplitude calibrator
solint        = 0.0                     # Calibration solution time
solsmo        = 0.0                     # Smooth calibration
ampScalar     = False                   # Amp-scalar operation in Calib
AcalModel     = None                    # Flux calibrator model file name, None=>use point
AcalFlux      = None                    # Flux for amp calibrator, None=>use model or SetJy value
AcalDisk      = 1                       # Flux calibrator model FITS disk
refAnt        = 0                       # Reference antenna

# Apply calibration and average?
doCalAvg      = True                    # calibrate and average data
avgClass      = "UVAvg"                 # AIPS class of calibrated/averaged uv data
CalAvgTime    = 10.0                    # Time for averaging calibrated uv data (sec)
CABIF         = 1                       # First IF to copy
CAEIF         = 0                       # Highest IF to copy

# Poln  Cal
doPolCal      = False     # Do polarization calibration?
PCInsCals     = []        # List of instrumental poln calibrators
PCFixPoln     = False     # Fix polarization to value in source table?
PCAvgIF       = False     # Determine average calibration for all IFs?
PCSolInt      = 2.0       # Pcal solution averaging time in min.
PCRefAnt      = 0         # Poln cal reference antenna

# R-L phase/delay calibration
doRLCal       = False    # Determine R-L bandpass?
RLPCal        = None     # Polarization angle (R-L phase) calibrator
PCRLPhase     = 0.0      # R-L phase difference for RLPCal
RM            = 0.0      # rotation measure (rad/m^2) for RLPCal
rlBChan       = 1        # Low freq. channel
rlEChan       = 0        # High freq. channel
rlCalCode     = "    "   # Calcode to select
rlDoCal       = -1       # docalib
rlgainUse     = 0        # CL/SN table to apply
rltimerange   = [0.,0.,0.,0.,0.,0.,0.,0.]
rlDoBand      = -1       # Apply and type of bandpass
rlBPVer       = 0        # BP table to apply, 0=>highest
rlflagVer     = 0        # Flagging table to apply, 0=>highest
rlrefAnt      = 1        # Reference antenna REQUIRED
rldoPol       = False    # Apply polarization cal?
rlChWid       = 1        # Number of channels in running mean phase BP soln
rlBPSoln      = 0        # Number of output RL phase BP table
rlsolint1     = 10.0/60.0 # RL BPass phase correction solution in min
rlsolint2     = 10.0      # RL BPass bandpass solution in min

# Imaging
doImage     = True         # Image targets
targets     = []           # targets
outIclass   = "IClean"     # Output image class
Stokes      = "I"          # Stokes parameter(s) to image
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
doPol       = False        # Poln cal in imaging
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
doSaveUV      = True       # Save uv data
doSaveImg     = True       # Save images
doSaveTab     = True       # Save Tables
doCleanup     = True       # Destroy AIPS files
prtLv         = 2          # Amount of diagnostics

############################# Set Project Processing parameters ##################
parmFile = sys.argv[2]
execfile (parmFile)

################################## Process #####################################
# Logging directly to logFile
OErr.PInit(err, prtLv, logFile)
retCode = 0

mess = "Start project "+project
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
                            selBand=selBand, selChan=selChan, selNIF=selNIF, calInt=calInt, \
                            Compress=True, check=check, debug=debug)
    if uv==None and not check:
        raise RuntimeError,"Cannot load "+inFile
# Otherwise set uv
if uv==None and not check:
    uv = UV.newPAUV("AIPS UV DATA", project, dataClass, disk, seq, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")

# Clear any old calibration/editing 
if doClearTab:
    EVLAClearCal(uv, err, doGain=doGain, doFlag=doFlag, doBP=doBP, check=check)
    OErr.printErrMsg(err, "Error resetting calibration")

# Copy FG 1 to FG 2
if doCopyFG:
    retCode = EVLACopyFG (uv, err, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error Copying FG table"

# Special editing
if doEditList and not check:
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
                              timeWind=timeWind, flagVer=2,flagSig=mednSigma, logfile=logFile, \
                              check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in MednFlag"

# Bandpass calibration if needed
if doBPCal and BPCal:
    retCode = EVLABPCal(uv, BPCal, err, noScrat=noScrat, solInt1=bpsolint2, solInt2=bpsolint2, solMode=bpsolMode, \
                            BChan1=bpBChan1, EChan1=bpEChan1, BChan2=bpBChan2, EChan2=bpEChan2, ChWid2=bpChWid2, \
                            doCenter1=bpDoCenter1, refAnt=refAnt, specIndex=specIndex, flagVer=2, logfile=logFile, \
                            check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in Bandpass calibration"

# Amp & phase Calibrate
if doAmpPhaseCal:
    retCode = EVLACal (uv, targets, ACal, err, PCal=PCal, doBand=1, BPVer=1, flagVer=2, \
                           calModel=AcalModel, calDisk=AcalDisk, calFlux=AcalFlux, nThreads=nThreads, \
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

# Calibrate and average data
if doCalAvg:
#    retCode = EVLACalAvg (uv, avgClass, seq, CalAvgTime, err, \
#                          flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
#                          BIF=CABIF, EIF=CAEIF, \
#                          FOV=FOV, maxFact=1.004, Compress=Compress, \
#                          logfile=logFile, check=check, debug=debug)
    retCode = EVLACalAvg2(uv, avgClass, seq, CalAvgTime, err, \
                          flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                          BIF=CABIF, EIF=CAEIF, Compress=Compress, \
                          logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in CalAvg"
   
# Get calibrated/averaged data
if not check:
    uv = UV.newPAUV("AIPS UV DATA", project, avgClass, disk, seq, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating cal/avg AIPS data")


# R-L  delay calibration cal if needed, creates new CL table
# THIS NEEDS CLEANUP
if doRLCal:
    retCode = EVLARLCal(uv, err, RLPCal=RLPCal, RLPhase=PCRLPhase, RM=RM, \
                        BChan=rlBChan, EChan=rlEChan, ChWid=rlChWid, \
                        calcode=rlCalCode, doCalib=rlDoCal, gainUse=rlgainUse, \
                        timerange=rltimerange, \
                        doBand=rlDoBand, BPVer=rlBPVer, flagVer=rlflagVer, \
                        refAnt=rlrefAnt, doPol=rldoPol,  \
                        BPSoln=rlBPSoln, solInt1=rlsolint1, solInt2=rlsolint2, \
                        nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                        check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in RL phase spectrum calibration"

# Polarization calibration
if doPolCal:
    retCode = EVLAPolCal(uv, PCInsCals, err, \
                         doCalib=-1, gainUse=0, doBand=-1, BPVer=0, flagVer=0, \
                         fixPoln=PCFixPoln, avgIF=PCAvgIF, \
                         solInt=PCSolInt, refAnt=PCRefAnt, \
                         check=check, debug=debug, \
                         noScrat=noScrat, logfile =logFile)
    if retCode!=0 and (not check):
       raise  RuntimeError,"Error in polarization calibration: "+str(retCode)
    # end poln cal.


# Image targets
if doImage:
    img = EVLASetImager(uv, targets, outIclass=outIclass, nThreads=nThreads, \
                            check=check, logfile=logFile)
    img.Stokes      = Stokes       # Stokes parameters to image
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
    img.flagVer     = 2            # Flag table version
    img.doPol       = doPol        # Polcal?
    # Bandpass calibration?
    if BPCal:
        img.doBand=1; img.BPVer = 0;
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
# Save UV data? 
if doSaveUV and (not check):
    filename = project+"Cal.uvtab"
    fuv = EVLAUVFITS (uv, filename, 0, err, compress=Compress)
# Save UV data tables?
if doSaveTab and (not check):
    filename = project+"CalTab.uvtab"
    fuv = EVLAUVFITSTab (uv, filename, 0, err)
# Delete UV data
if doCleanup and (not check):
    uv.Zap(err)
# Imaging results
for target in targets:
    if doSaveImg and (not check):
        for s in Stokes:
            oclass = s+outIclass[1:]
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            outfile = target+"."+outIclass+".fits"
            mess ="Write " +outfile+" on disk "+str(outDisk)
            printMess(mess, logFile)
            xf = EVLAImFITS (x, outfile, outDisk, err)
            # Statistics
            zz=imstat(x, err, logfile=logFile)
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

