# Pipeline processing for calibrating and imaging EVLA data
# AIPS/FITS setup and Parameter file given as arguments:
# > ObitTalk EVLAPipeline.py AIPSSetup.py parms.py
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
project       = "Unspecified"   # Project name (12 char or less, used as AIPS Name)
session       = "?"             # Project session code
band          = "?"             # Observing band
seq           = 1               # AIPS sequence number
gain          = 0.10            # CLEAN loop gain
check         = False           # Only check script, don't execute tasks
debug         = False           # run Obit tasks debug

# Archive parameters
doLoadArchive = True       # Load from archive?
archRoot      = "."        # Archive directory root
selBand       = " "        # Selected band, def = first  
selConfig     = 0          # Selected frequency config, def = first  
selNIF        = 0          # Selected number of IFs, def = first  
selChan       = 0          # Selected number of channels, def = first  
calInt        = 0.5        # Calibration table interval in min.
doSwPwr       = False      # Make EVLA Switched power corr?

# Load from FITS
doLoad       = False       # Load FITS file?
FITSIn       = ""          # Input FITS file
FITSinDisk   = 0           # Input FITS disk

# General data parameters
dataClass     = "UVData"   # AIPS class of raw uv data
Compress      = False      # Use compressed UV data?

# Editing
doClearTab  = True         # Clear cal/edit tables at (re)start
doClearGain = True         # Clear SN and CL tables >1
doClearFlag = True         # Clear FG tables > 1
doClearBP   = True         # Clear BP tables?
doCopyFG    = True         # Copy FG 1 to FG 2 
Reason      = "Quack"      # Reason code for Quack
doQuack     = True         # Quack data?
quackBegDrop= 0.1          # Time to drop from start of each scan in min
quackEndDrop= 0.1          # Time to drop from end of each scan in min
quackReason   = "Quack"    # Reason string

doMedn      = True         # Median editing?
mednSigma   = 10.0         # Median sigma clipping level
timeWind    = 2.0          # Median window width in min for median flagging
avgTime     = 10.0/60.     # Averaging time in min
avgFreq     = 0            # 1=>avg chAvg chans, 2=>avg all chan, 3=> avg chan and IFs
chAvg       = 1            # number of channels to average

doFD1       = True         # Do initial frequency domain flagging
FD1widMW    = 15           # Width of the initial FD median window
FD1maxRes   = 10.0         # Clipping level in sigma
FD1TimeAvg  = 2.0          # Time averaging in min. for initial FD flagging

doRMSAvg    = True         # Edit calibrators by RMSAvg?
RMSAvg      = 5.0          # AutoFlag Max RMS/Avg for time domain RMS filtering
RMSTimeAvg  = 1.0          # AutoFlag time averaging in min.

doAutoFlag    = True       # Autoflag editing?
RMSAvg      = 20.0         # AutoFlag Max RMS/Avg for time domain RMS filtering
IClip       = [1000.0,0.1] # AutoFlag Stokes I clipping
minAmp      = 0.0001       # Min. allowable calibrated amplitude
VClip       = [10.0,0.05]  # AutoFlag Stokes V clipping
timeAvg     = 2.0          # AutoFlag time averaging in min.
doAFFD      = False        # do AutoFlag frequency domain flag
FDmaxAmp    = IClip[0]     # Maximum average amplitude (Jy)
FDmaxV      = VClip[0]     # Maximum average VPol amp (Jy)
FDwidMW     = 15           # Width of the median window
FDmaxRMS    = [6.0,0.1]    # Channel RMS limits (Jy)
FDmaxRes    = 10.0         # Max. residual flux in sigma
FDmaxResBL  = 10.0         # Max. baseline residual
FDbaseSel   = [0,0,0,0]    # Channels for baseline fit

# Special editing list
doEditList  = False        # Edit using editList?
editFG      = 2            # Table to apply edit list to
editList = [
#    {"timer":("0/06:09:0.0","0/06:13:0.0"),"Ant":[ 8,0],"IFs":[2,2],"Chans":[1,0],"Stokes":'1110',"Reason":"bad data"},
    ]

# Parallactic angle correction, need for PCAL w/ "ORI-"
doPACor     = True         # Make parallactic angle correction

# Bandpass Calibration
doBP        = True         # Clear BP tables
doBPCal     = True         # Determine Bandpass calibration
BPCal       = None         # Bandpass calibrator
bpSpecIndex = 0.0          # BPCal spectral index
bpBChan1      = 1          # Low freq. channel,  initial cal
bpEChan1      = 0          # Highest freq channel, initial cal, 0=>all
bpDoCenter1   = None       # Fraction of  channels in 1st, overrides bpBChan1, bpEChan1
bpBChan2      = 1          # Low freq. channel for BP cal
bpEChan2      = 0          # Highest freq channel for BP cal,  0=>all 
bpChWid2      = 1          # Number of channels in running mean BP soln
bpsolMode     = 'A&P'      # Band pass type 'A&P', 'P', 'P!A'
bpsolint1     = 10.0/60.0  # BPass phase correction solution in min
bpsolint2     = 10.0       # BPass bandpass solution in min

# Delay calibration from DCal (see below)
doDelayCal  = True         # Determine/apply delays from DCal
doDelayCal2 = True         # Determine/apply delays from DCal on averaged data
doTwo       = True         # Use 1 and 2 bl combinations?

# Amp/phase calibration
PCal          = None                    # Phase calibrator(s)
ACal          = None                    # Amplitude calibrator
DCal          = None                    # Delay calibrator(s)
solint        = 0.0                     # Calibration solution time
solsmo        = 0.0                     # Smooth calibration
ampScalar     = False                   # Amp-scalar operation in Calib
AcalModel     = None                    # Flux calibrator model file name, None=>use point
AcalFlux      = None                    # Flux for amp calibrator, None=>use model or SetJy value
AcalDisk      = 1                       # Flux calibrator model FITS disk
refAnt        = 0                       # Reference antenna

# Sample spectra
doRawSpecPlot = True        # Plot Raw spectrum
doSpecPlot    = True        # Plot spectrum at various stages of processing
plotSource    = None        # Source or None
plotTime      = [0.,1.]     # timerange

# Apply calibration and average?
doCalAvg      = True                    # calibrate and average data
avgClass      = "UVAvg"                 # AIPS class of calibrated/averaged uv data
CalAvgTime    = 10.0/60.0               # Time for averaging calibrated uv data (min)
CABIF         = 1                       # First IF to copy
CAEIF         = 0                       # Highest IF to copy
avgFreq       = 0                       # average channels, 0=no, 1=yes, 2=all, 3=also IFs
chAvg         = 1                       # Number of channels to average

# Poln  Cal
doPolCal      = False     # Do polarization calibration?
PCInsCals     = []        # List of instrumental poln calibrators
PCpmodel      = [1.0,0.,0.,0.,0.,0.]     # Polarization model
PCFixPoln     = False     # Fix polarization to value in source table?
PCAvgIF       = False     # Determine average calibration for all IFs?
PCSolInt      = 2.0       # Pcal solution averaging time in min.
PCRefAnt      = 0         # Poln cal reference antenna
PCSolType     ="ORI-"     # Solution type "ORI-", "RAPR"
doRecal       = True      # Redo calibration after editing

# R-L phase/delay calibration
doRLCal       = False    # Determine R-L bandpass prior to pcal?
doRLCal2      = False    # Determine R-L bandpass after pcal?
RLPCal        = None     # Polarization angle (R-L phase) calibrator
PCRLPhase     = 0.0      # R-L phase difference for RLPCal
RM            = 0.0      # rotation measure (rad/m^2) for RLPCal
RLDCal        = None     # R-L delay calibrator
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
rldoPol       = True     # Apply polarization cal?
rlChWid       = 1        # Number of channels in running mean phase BP soln
rlBPSoln      = 0        # Number of output RL phase BP table
rlsolint1     = 10.0/60.0 # RLPass phase correction solution in min
rlsolint2     = 10.0      # RLPass bandpass solution in min
rlUVRange     = [0.,0.]   # Range of baseline used in kilowavelengths in RL Delay cal.

# Imaging
doImage     = True         # Image targets
targets     = []           # targets, empty = all
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
solPMode    = "P"          # Delay solution for phase self cal
solAType    = "  "         # L1 solution for amp&phase self cal
solAMode    = "A&P"        # Delay solution for amp&phase self cal
avgPol      = False        # Average poln in self cal?
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
CleanRad    = None         # CLEAN radius about center or None=autoWin

# Final
doReport      = True       # Individual source report
doSaveUV      = True       # Save uv data
doSaveImg     = True       # Save images
doSaveTab     = True       # Save Tables
doCleanup     = True       # Destroy AIPS files
prtLv         = 2          # Amount of diagnostics
doMB          = True       # Use wideband imaging?
MBnorder      = 2          # Order of wideband imaging
MBmaxFBW      = 0.05       # Max. fractional IF center bandwidth

############################# Set Project Processing parameters ##################
parmFile = sys.argv[2]
execfile (parmFile)

################################## Process #####################################
logFile       = project+"_"+session+"_"+band+".log"  # Processing log file
# Logging directly to logFile
OErr.PInit(err, prtLv, logFile)
retCode = 0

mess = "Start project "+project+" session "+session+" "+band+" Band"+" AIPS user no. "+str(AIPS.userno)
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
    uv = EVLAUVLoadT(FITSIn, FITSinDisk, project+session, dataClass, disk, seq, err, logfile=logFile, \
                         Compress=True, check=check, debug=debug)
    if uv==None:
        raise RuntimeError,"Cannot load "+inFile
# Load Data from Archive directory
if doLoadArchive:
    uv = EVLAUVLoadArch(archRoot, project+session, dataClass, disk, seq, err, \
                            selConfig=selConfig, doSwPwr=doSwPwr, \
                            selBand=selBand, selChan=selChan, selNIF=selNIF, calInt=calInt, \
                            logfile=logFile, Compress=True, check=check, debug=debug)
    if uv==None and not check:
        raise RuntimeError,"Cannot load "+inFile
# Otherwise set uv
if uv==None and not check:
    uv = UV.newPAUV("AIPS UV DATA", project+session, dataClass, disk, seq, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")

# Clear any old calibration/editing 
if doClearTab:
    mess =  "Clear previous calibration"
    printMess(mess, logFile)
    EVLAClearCal(uv, err, doGain=doClearGain, doFlag=doClearFlag, doBP=doClearBP, check=check)
    OErr.printErrMsg(err, "Error resetting calibration")

# Copy FG 1 to FG 2
if doCopyFG:
    retCode = EVLACopyFG (uv, err, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error Copying FG table"

# Special editing
if doEditList and not check:
    mess =  "Special editing"
    printMess(mess, logFile)
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

# Median window time editing, for RFI impulsive in time
if doMedn:
    mess =  "Median window time editing, for RFI impulsive in time:"
    printMess(mess, logFile)
    retCode = EVLAMedianFlag (uv, "    ", err, noScrat=noScrat, nThreads=nThreads, \
                              avgTime=avgTime, avgFreq=avgFreq,  chAvg= chAvg, \
                              timeWind=timeWind, flagVer=2,flagSig=mednSigma, logfile=logFile, \
                              check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in MednFlag"

# Median window frequency editing, for RFI impulsive in frequency
if doFD1:
    mess =  "Median window frequency editing, for RFI impulsive in frequency:"
    printMess(mess, logFile)
    retCode = EVLAAutoFlag (uv, targets, err,  flagVer=2, doCalib=-1, doBand=-1,   \
                                timeAvg=FD1TimeAvg, \
                                doFD=True, FDmaxAmp=1.0e20, FDmaxV=1.0e20, FDwidMW=FD1widMW,  \
                                FDmaxRMS=[1.0e20,0.1], FDmaxRes=FD1maxRes,  FDmaxResBL= FD1maxRes,  \
                                logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in AutoFlag"


# RMS/Mean editing for calibrators
if doRMSAvg:
    mess =  "RMS/Mean editing for calibrators:"
    printMess(mess, logFile)
    clist = [ACal]   # Calibrator list
    for s in PCal:
        clist.append(s)
    retCode = EVLAAutoFlag (uv, clist, err,  flagVer=2, doCalib=-1, doBand=-1,   \
                                RMSAvg=RMSAvg, timeAvg=RMSTimeAvg, \
                                logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in AutoFlag"


# Plot Raw, edited data?
if doRawSpecPlot and plotSource:
    plotFile = "./"+project+"_"+session+"_"+band+"RawSpec.ps"
    retCode = EVLASpectrum(uv, plotSource, plotTime, plotFile, refAnt, err, \
                           Stokes=["RR","LL"], doband=-1,          \
                           check=check, debug=debug, logfile=logFile )
    if retCode!=0:
        raise  RuntimeError,"Error in Plotting spectrum"

# Parallactic angle correction?
if doPACor:
    retCode = EVLAPACor(uv, err, noScrat=noScrat, \
                            logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in Parallactic angle correction"

# Bandpass calibration
if doBPCal and BPCal:
    retCode = EVLABPCal(uv, BPCal, err, noScrat=noScrat, solInt1=bpsolint2, solInt2=bpsolint2, solMode=bpsolMode, \
                        BChan1=bpBChan1, EChan1=bpEChan1, BChan2=bpBChan2, EChan2=bpEChan2, ChWid2=bpChWid2, \
                        doCenter1=bpDoCenter1, refAnt=refAnt, specIndex=bpSpecIndex, \
                        doCalib=2, gainUse=0, flagVer=2, doPlot=False, \
                        logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in Bandpass calibration"
    
    # Plot corrected data?
    if doSpecPlot and plotSource:
        plotFile = "./"+project+"_"+session+"_"+band+"BPSpec.ps"
        retCode = EVLASpectrum(uv, plotSource, plotTime, plotFile, refAnt, err, \
                               Stokes=["RR","LL"], doband=1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"

# delay calibration
if doDelayCal and DCal and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"DelayCal.ps"
    retCode = EVLADelayCal(uv, err, calSou=DCal, CalModel=None, \
                           doCalib=2, flagVer=2, doBand=1, \
                           solInt=solint, smoTime=20.0/60.0,  \
                           refAnts=[refAnt], doTwo=doTwo, \
                           doPlot=doSNPlot, plotFile=plotFile, \
                           nThreads=nThreads, noScrat=noScrat, \
                           logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in delay calibration"
    
    # Plot corrected data?
    if doSpecPlot and plotSource:
        plotFile = "./"+project+"_"+session+"_"+band+"DelaySpec.ps"
        retCode = EVLASpectrum(uv, plotSource, plotTime, plotFile, refAnt, err, \
                               Stokes=["RR","LL"], doband=1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"

# Amp & phase Calibrate
if doAmpPhaseCal:
    plotFile = "./"+project+"_"+session+"_"+band+"APCal.ps"
    retCode = EVLACalAP (uv, targets, ACal, err, PCal=PCal, doCalib=2, doBand=1, BPVer=1, flagVer=2, \
                         calModel=AcalModel, calDisk=AcalDisk, calFlux=AcalFlux, nThreads=nThreads, \
                         solInt=solint, solSmo=solsmo, noScrat=noScrat, ampScalar=ampScalar, \
                         doPlot=doSNPlot, plotFile=plotFile, \
                         refAnt=refAnt, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error calibrating"

# More editing, with flux limits
if doAutoFlag:
    mess =  "Post calibration editing:"
    printMess(mess, logFile)
    retCode = EVLAAutoFlag (uv, targets, err, flagVer=2, \
                                doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                                IClip=IClip, minAmp=minAmp, VClip=VClip, timeAvg=timeAvg, \
                                doFD=doAFFD, FDmaxAmp=FDmaxAmp, FDmaxV=FDmaxV, \
                                FDwidMW=FDwidMW, FDmaxRMS=FDmaxRMS, \
                                FDmaxRes=FDmaxRes,  FDmaxResBL= FDmaxResBL, FDbaseSel=FDbaseSel, \
                                logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in AutoFlag"

# Redo the calibration using new flagging?
if doRecal:
    mess =  "Redo calibration:"
    printMess(mess, logFile)
    EVLAClearCal(uv, err, doGain=True, doFlag=False, doBP=True, check=check)
    OErr.printErrMsg(err, "Error resetting calibration")
    # Parallactic angle correction?
    if doPACor:
        retCode = EVLAPACor(uv, err, noScrat=noScrat, \
                            logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Parallactic angle correction"

    # Bandpass calibration
    if doBPCal and BPCal:
        retCode = EVLABPCal(uv, BPCal, err, noScrat=noScrat, \
                            solInt1=bpsolint2, solInt2=bpsolint2, solMode=bpsolMode, \
                            BChan1=bpBChan1, EChan1=bpEChan1, BChan2=bpBChan2, EChan2=bpEChan2, ChWid2=bpChWid2, \
                            doCenter1=bpDoCenter1, refAnt=refAnt, specIndex=bpSpecIndex, \
                            doCalib=2, gainUse=0, flagVer=2, doPlot=False, \
                            logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Bandpass calibration"


    # Amp & phase Calibrate
    if doAmpPhaseCal:
        plotFile = "./"+project+"_"+session+"_"+band+"APCal.ps"
        retCode = EVLACalAP (uv, targets, ACal, err, PCal=PCal, doCalib=2, doBand=1, BPVer=1, flagVer=2, \
                             calModel=AcalModel, calDisk=AcalDisk, calFlux=AcalFlux, nThreads=nThreads, \
                             solInt=solint, solSmo=solsmo, noScrat=noScrat, ampScalar=ampScalar, \
                             doPlot=doSNPlot, plotFile=plotFile, \
                             refAnt=refAnt, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error calibrating"

# end recal

# Calibrate and average data
# Note, PCAL doesn't handle flagging so this has to be here
if doCalAvg:
    retCode = EVLACalAvg (uv, avgClass, seq, CalAvgTime, err, \
                          flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1, doPol=False, \
                          avgFreq=avgFreq, chAvg=chAvg, \
                          BIF=CABIF, EIF=CAEIF, Compress=Compress, \
                          logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in CalAvg"
   
# Get calibrated/averaged data
if not check:
    uv = UV.newPAUV("AIPS UV DATA", project+session, avgClass, disk, seq, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating cal/avg AIPS data")


# delay calibration
if doDelayCal2 and DCal and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"DelayCal.ps"
    retCode = EVLADelayCal(uv, err, calSou=DCal, CalModel=None, \
                           doCalib=2, flagVer=2, doBand=-1, \
                           solInt=solint*2, smoTime=5.0/60.0,  \
                           refAnts=[refAnt], doTwo=doTwo, \
                           doPlot=doSNPlot, plotFile=plotFile, \
                           nThreads=nThreads, noScrat=noScrat, \
                           logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in delay calibration"

# Plot corrected data?
if doSpecPlot and plotSource:
    plotFile = "./"+project+"_"+session+"_"+band+"Spec.ps"
    retCode = EVLASpectrum(uv, plotSource, plotTime, plotFile, refAnt, err, \
                           Stokes=["RR","LL"], doband=-1,          \
                           check=check, debug=debug, logfile=logFile )
    if retCode!=0:
        raise  RuntimeError,"Error in Plotting spectrum"

# R-L  delay calibration cal if needed, creates new BP table
if doRLCal:
    plotFile = "./"+project+"_"+session+"_"+band+"RLSpec.ps"
    retCode = EVLARLCal(uv, err,\
                        RLDCal=RLDCal, BChan=rlBChan, EChan=rlEChan, UVRange=rlUVRange, \
                        ChWid2=rlChWid, solInt1=rlsolint1, solInt2=rlsolint2, \
                        RLPCal=None, RLPhase=PCRLPhase, RM=RM, CleanRad=CleanRad, \
                        calcode=rlCalCode, doCalib=rlDoCal, gainUse=rlgainUse, \
                        timerange=rltimerange, FOV=FOV, \
                        doBand=rlDoBand, BPVer=rlBPVer, flagVer=rlflagVer, \
                        refAnt=rlrefAnt, doPol=False,  \
                        doPlot=doSNPlot, plotFile=plotFile, \
                        nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                        check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in RL phase spectrum calibration"

# Polarization calibration 
if doPolCal:
    retCode = EVLAPolCal(uv, PCInsCals, err, \
                         doCalib=2, gainUse=0, doBand=1, BPVer=0, flagVer=0, \
                         fixPoln=PCFixPoln, pmodel=PCpmodel, avgIF=PCAvgIF, \
                         solInt=PCSolInt, refAnt=PCRefAnt, soltype=PCSolType, \
                         check=check, debug=debug, \
                         noScrat=noScrat, logfile=logFile)
    if retCode!=0 and (not check):
       raise  RuntimeError,"Error in polarization calibration: "+str(retCode)
    # end poln cal.


# R-L  delay calibration cal if needed, creates new BP table (again)
if doRLCal2:
    plotFile = "./"+project+"_"+session+"_"+band+"RLSpec2.ps"
    retCode = EVLARLCal(uv, err,\
                        RLDCal=RLDCal, BChan=rlBChan, EChan=rlEChan, UVRange=rlUVRange, \
                        ChWid2=rlChWid, solInt1=rlsolint1, solInt2=rlsolint2, \
                        RLPCal=RLPCal, RLPhase=PCRLPhase, RM=RM, CleanRad=CleanRad, \
                        calcode=rlCalCode, doCalib=rlDoCal, gainUse=rlgainUse, \
                        timerange=rltimerange, FOV=FOV, \
                        doBand=-1, BPVer=1, flagVer=rlflagVer, \
                        refAnt=rlrefAnt, doPol=True,  \
                        doPlot=doSNPlot, plotFile=plotFile, \
                        nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                        check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in RL phase spectrum calibration"

# Image targets
if doImage:
    # If targets not specified, image all
    if len(targets)<=0:
        slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    EVLAImageTargets (uv, err, Sources=slist, seq=seq, sclass=outIclass, \
                      doCalib=2, doBand=1,  flagVer=2, doPol=doPol,  \
                       Stokes=Stokes, FOV=FOV, Robust=Robust, Niter=Niter, \
                       CleanRad=CleanRad, \
                       maxPSCLoop=maxPSCLoop, minFluxPSC=minFluxPSC, \
                       solPInt=solPInt, solPMode=solPMode, solPType=solPType, \
                       maxASCLoop=maxASCLoop, minFluxASC=minFluxASC, \
                       solAInt=solAInt, solAMode=solAMode, solAType=solAType, \
                       avgPol=avgPol, avgIF=avgIF, minSNR = 3.0, refAnt=refAnt, \
                       do3D=do3D, BLFact=BLFact, BLchAvg=BLchAvg, \
                       doMB = doMB, norder=MBnorder, maxFBW=MBmaxFBW, \
                       nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                       check=check, debug=debug)                       
    # End image

# Get report on sources
if doReport:
    # If targets not specified, image all
    if len(targets)<=0:
        slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    Report = EVLAReportTargets(uv, err, Sources=slist, seq=seq, sclass=outIclass, \
                                   Stokes=Stokes, logfile=logFile, check=check, debug=debug)
    # Save to pickle jar
    ReportPicklefile = "./"+project+"_"+session+"_"+band+"Report.pickle"   # Where results saved
    SaveObject(Report, ReportPicklefile, True) 
   
# Write results, cleanup    
# Save UV data? 
if doSaveUV and (not check):
    filename = project+session+band+"Cal.uvtab"
    fuv = EVLAUVFITS (uv, filename, 0, err, compress=Compress)
# Save UV data tables?
if doSaveTab and (not check):
    filename = project+session+band+"CalTab.uvtab"
    fuv = EVLAUVFITSTab (uv, filename, 0, err, logfile=logFile)
# Imaging results
# If targets not specified, save all
if len(targets)<=0:
    slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
else:
    slist = targets
for target in slist:
    if doSaveImg and (not check):
        for s in Stokes:
            oclass = s+outIclass[1:]
            # Test if image exists
            cno = AIPSDir.PTestCNO(disk, user, target, oclass, "MA", seq, err)
            if cno <= 0 :
                continue
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            outfile = project+session+band+target+"."+oclass+".fits"
            xf = EVLAImFITS (x, outfile, 0, err, logfile=logFile)
            # Statistics
            zz=imstat(x, err, logfile=logFile)
# end writing loop

# cleanup
if doCleanup and (not check):
    if len(targets)<=0:
        slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    mess ="Cleanup AIPS Files" 
    printMess(mess, logFile)
    uv.Zap(err)  # UV data
    # Initial data 
    uv = UV.newPAUV("AIPS UV DATA", project+session, dataClass, disk, seq, True, err)
    uv.Zap(err)  # UV data
    # If targets not specified, save all
    for target in slist:
        for s in Stokes:
            oclass = s+outIclass[1:]
            # Test if image exists
            cno = AIPSDir.PTestCNO(disk, user, target, oclass, "MA", seq, err)
            if cno <= 0 :
                mess = "Image"+target+" "+oclass+" Not found"
                printMess(mess, logFile)
                continue
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            x.Zap(err) # cleanup
# end cleanup

OErr.printErrMsg(err, "Writing output/cleanup")

# Shutdown
mess = "Finished project "+project+" session "+session+" "+band+" Band"+" AIPS user no. "+str(AIPS.userno)
printMess(mess, logFile)
OErr.printErr(err)
OSystem.Shutdown(ObitSys)

