# VLBA Continuum Pipeline processing 
# AIPS/FITS setup and Parameter file given as arguments:
# > python VLBAContPipe.py AIPSSetup.py parms.py
#
import sys
import OErr, OSystem, UV, AIPS, FITS
import ObitTalkUtil
import pydoc
from AIPS import AIPSDisk
from FITS import FITSDisk
from VLBACal import *

############################# Initialize OBIT ##########################################                                 
setup = sys.argv[1]
noScrat     = []    
execfile (setup)

############################# Default parameters ##########################################                                 
# Generic parameters
project       = "Unspecified"               # Project name (12 char or less, used as AIPS Name)
session       = "?"                         # Project session code
band          = "?"                         # Observing band
logFile       = project+"_"+session+"_"+band+".log"  # Processing log file
seq           = 1          # AIPS sequence number
gain          = 0.10       # CLEAN loop gain
doLoadIDI     = True       # Load data from IDI FITS?
doLoadUVF     = False      # Load the "AIPS Friendly" (uvfits) FITS  version
dataInUVF     = None       # Input uvfits data file name
dataInIDI     = None       # Input FITS-IDI file or list
dataClass     = "Raw"      # AIPS class of raw uv data
Compress      = False      # Use compressed UV data?
calInt        = 0.15       # Calibration table interval in min.
wtThresh      = 0.8        # Data weight  threshold
check         = False      # Only check script, don't execute tasks
debug         = False      # run tasks debug

# Quantization correction
doQuantCor  = True         # Do quantization correction
QuantSmo    = 0.5          # Smoothing time (hr) for quantization corrections
QuantFlag   = 0.0          # If >0, flag solutions < QuantFlag (use 0.9 for 1 bit, 0.8 for 2 bit)

# Parallactic angle correction
doPACor     = True         # Make parallactic angle correction

# Opacity/Tsys correction
doOpacCor   = True         # Make Opacity/Tsys/gain correction?
OpacSmoo    = 0.25         # Smoothing time (hr) for opacity corrections

# Special editing list
doEditList  = False        # Edit using editList?
editFG      = 2            # Table to apply edit list to
editList = [
#    {"timer":("0/06:09:0.0","0/06:13:0.0"),"Ant":[ 8,0],"IFs":[2,2],"Chans":[1,0],"Stokes":'1110',"Reason":"bad data"},
    ]

# Apply phase cal corrections?
doPCcor  = True           # Apply PC table?
doPCPlot = True           # Plot results?

# 'Manual" phase cal - even to tweak up PCals
doManPCal      = True      # Determine and apply manual phase cals?
manPCsolInt    = 0.5       # Manual phase cal solution interval (min)
manPCSmoo      = 10.0      # Manual phase cal smoothing time (min)
doManPCalPlot  = True      # Plot the phase and delays from manual phase cal

# Bandpass Calibration?
doBPCal       = True       # Determine Bandpass calibration
bpBChan1      = 1          # Low freq. channel,  initial cal
bpEChan1      = 0          # Highest freq channel, initial cal, 0=>all
bpDoCenter1   = None       # Fraction of  channels in 1st, overrides bpBChan1, bpEChan1
bpBChan2      = 1          # Low freq. channel for BP cal
bpEChan2      = 0          # Highest freq channel for BP cal,  0=>all 
bpChWid2      = 1          # Number of channels in running mean BP soln
bpdoAuto      = False      # Use autocorrelations rather than cross?
bpsolMode     = 'A&P'      # Band pass type 'A&P', 'P', 'P!A'
bpsolint1     = 10.0/60.0  # BPass phase correction solution in min
bpsolint2     = 10.0       # BPass bandpass solution in min
specIndex     = 0.0        # Spectral index of BP Cal
doSpecPlot    = True       # Plot the amp. and phase across the spectrum

# Editing
doClearTab  = True         # Clear cal/edit tables
doGain      = True         # Clear SN and CL tables >1
doFlag      = True         # Clear FG tables > 1
doBP        = True         # Clear BP tables?
doCopyFG    = True         # Copy FG 1 to FG 2quack
doQuack     = True         # Quack data?
quackBegDrop= 0.1          # Time to drop from start of each scan in min
quackEndDrop= 0.0          # Time to drop from end of each scan in min
quackReason = "Quack"      # Reason string

# Amp/phase calibration parameters
refAnt        = 0                       # Reference antenna
refAnts       = [0]                     # List of Reference antenna for fringe fitting

# Imaging calibrators (contCals) and targets
doImgCal    = True         # Image calibrators
targets     = []           # targets
doImgTarget = True         # Image targets?
outCclass   = "ICalSC"     # Output calibrator image class
outTclass   = "IImgSC"     # Output target temporary image class
outIclass   = "IClean"     # Output target final image class
Robust      = 0.0          # Weighting robust parameter
FOV         = 0.1/3600     # Field of view radius in deg.
Niter       = 100          # Max number of clean iterations
minFlux     = 0.0          # Minimum CLEAN flux density
minSNR      = 5.0          # Minimum Allowed SNR
solMode     = "P"          # Delay solution for phase self cal
avgPol      = True         # Average poln in self cal?
avgIF       = False        # Average IF in self cal?
maxPSCLoop  = 6            # Max. number of phase self cal loops
minFluxPSC  = 0.1          # Min flux density peak for phase self cal
solPInt     = 0.25         # phase self cal solution interval (min)
maxASCLoop  = 1            # Max. number of phase self cal loops
minFluxASC  = 0.5          # Min flux density peak for amp+phase self cal
solAInt     = 3.0          # amp+phase self cal solution interval (min)

# Find good calibration data
doFindCal    = True         # Search for good calibration/reference antenna
findSolInt   = 0.5          # Solution interval (min) for Calib
findTimeInt  = 10.0         # Maximum timerange, large=>scan
contCals     = None         # Name or list of continuum cals
contCalModel = VLBAImageModel(contCals,outCclass, disk, seq, err)         # Check for any
targetModel  = None         # No target model yet

# Delay calibration
doDelayCal  = True         # Determine/apply delays from contCals

# Amplitude calibration
doAmpCal    = True         # Determine/smooth/apply amplitudes from contCals

# Apply calibration and average?
doCalAvg      = True                    # calibrate and average cont. calibrator data
avgClass      = "UVAvg"                 # AIPS class of calibrated/averaged uv data
CalAvgTime    = 10.0/60.0               # Time for averaging calibrated uv data (min)
CABIF         = 1                       # First IF to copy
CAEIF         = 0                       # Highest IF to copy
CABChan       = 1                       # First Channel to copy
CAEChan       = 0                       # Highest Channel to copy
chAvg         = 10000000                # Average all channels
avgFreq       = 1                       # Average all channels

# Phase calibration of all targets in averaged calibrated data
doPhaseCal    = True       # Phase calibrate all data with self-cal?

# Instrumental polarization cal?
doInstPol     = True      # determination instrumental polarization from instPolCal
instPolCal    = None      # Defaults to contCals

# Right-Left phase (EVPA) calibration 
doRLCal      = True       # Set RL phases from RLCal - also needs RLCal
RLCal        = None       # RL Calibrator source name, if given, a list of triplets, 
                          # (name, R-L phase(deg@1GHz), RM (rad/m^2))

# Final Image/Clean
doImgFullTarget = True    # Final Image/Clean/selfcal
Stokes          = "I"     # Stokes to image
doKntrPlots     = True    # Contour plots

# Final
outDisk       = 0          # FITS disk number for output (0=cwd)
doSaveUV      = True       # Save uv data
doSaveImg     = True       # Save images
doSaveTab     = True       # Save Tables
doCleanup     = True       # Destroy AIPS files

# diagnostics
doSNPlot      = True       # Plot SN tables etc
doDiagPlots   = True       # Plot single source diagnostics
prtLv         = 2          # Amount of task print diagnostics
doReport      = True       # Individual source report2

############################# Set Project Processing parameters ##################
parmFile = sys.argv[2]
execfile (parmFile)

################################## Process #####################################
# Init cal pickle jar
goodCalPicklefile = "./"+project+"_"+session+"_"+band+"GoodCal.pickle"   # Where results saved
# Default "best" calibration
goodCal = {"Source":"  ", "souID":0,"timeRange":(0.0,100.0), "Fract":0.0, "SNR":0.0, "bestRef":0}
SaveObject(goodCal, goodCalPicklefile, False)   # Save initial default

# Logging directly to logFile
OErr.PInit(err, prtLv, logFile)
retCode = 0

mess = "Start project "+project+" session "+session+" "+band+" Band"+" AIPS user no. "+str(AIPS.userno)
printMess(mess, logFile)
if debug:
    pydoc.ttypager = pydoc.plainpager # don't page task input displays
    mess = "Using Debug mode "
    printMess(mess, logFile)
if check:
    mess = "Only checking script"
    printMess(mess, logFile)

# Load Data from FITS
uv  = None   # Raw data
uvc = None   # Cal/averaged data
if doLoadIDI:
    if type(dataInIDI)==list:
        # Read list
        for dataIn in datainIDI:
            uv = VLBAIDILoad(dataIn, project, session, band, dataClass, disk, seq, err, logfile=logFile, \
                                 wtThresh=wtThresh, calInt=calInt, Compress=Compress, \
                                 check=check, debug=debug)
            if not UV.PIsA(uv):
                raise RuntimeError,"Cannot load "+dataIn
    else:
        # Single IDI file
        uv = VLBAIDILoad(dataInIDI, project, session, band, dataClass, disk, seq, err, logfile=logFile, \
                             wtThresh=wtThresh, calInt=calInt, Compress=Compress, \
                             check=check, debug=debug)
        if not UV.PIsA(uv):
            raise RuntimeError,"Cannot load "+dataInIDI
if doLoadUVF:
    # Single uv fits file
    uv = VLBAIDILoad(dataInUVF, project, session, band, dataClass, disk, seq, err, logfile=logFile, \
                         wtThresh=wtThresh, calInt=calInt, Compress=Compress, \
                         check=check, debug=debug)
    # Adding check condition to avoid error when checking
    if not UV.PIsA(uv) and not check:
        raise RuntimeError,"Cannot load "+dataInUVF
# Otherwise set uv
if uv==None and not check:
    Aname = VLBAAIPSName(project, session)
    uv = UV.newPAUV("AIPS UV DATA", Aname, dataClass, disk, seq, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")

# Clear any old calibration/editing 
if doClearTab:
    VLBAClearCal(uv, err, doGain=doGain, doFlag=doFlag, doBP=doBP, check=check, logfile=logFile)
    OErr.printErrMsg(err, "Error resetting calibration")

# Copy FG 1 to FG 2
if doCopyFG:
    retCode = VLBACopyFG (uv, err, logfile=logFile, check=check, debug=debug)
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
    retCode = VLBAQuack (uv, err, begDrop=quackBegDrop, endDrop=quackEndDrop, Reason=quackReason, \
                             logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error Quacking data"

# Quantization correction?
if doQuantCor:
    plotFile = "./"+project+"_"+session+"_"+band+"Quant.ps"
    retCode = VLBAQuantCor(uv, QuantSmo, QuantFlag, err, \
                               doSNPlot=doSNPlot, plotFile=plotFile, \
                               logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in quantization correcting/flagging"

# Parallactic angle correction?
if doPACor:
    retCode = VLBAPACor(uv, err, noScrat=noScrat, \
                            logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in quantization correcting/flagging"

# Opacity/Tsys/gain correction
if doOpacCor:
    plotFile = "./"+project+"_"+session+"_"+band+"Opacity.ps"
    retCode = VLBAOpacCor(uv, OpacSmoo, err,  \
                              doSNPlot=doSNPlot, plotFile=plotFile, \
                              logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in opacity/gain.Tsys correction"

# Find best calibration
if doFindCal:
    goodCal = VLBAGoodCal(uv,  err, solInt=findSolInt, timeInt=findTimeInt, \
                          calSou=contCals, \
                          #CalModel=contCalModel, \
                          doCalib=-1, flagVer=2, refAnts=refAnts, \
                          noScrat=noScrat, nThreads=nThreads, \
                          logfile=logFile, check=check, debug=debug)
    if not goodCal and not check:
        raise RuntimeError,"Error in finding best calibration data"
    # Save it to a pickle jar
    SaveObject(goodCal, goodCalPicklefile, True)
else:
    # Fetch from pickle
    goodCal = FetchObject(goodCalPicklefile)

# Apply Phase cals from PC table?
if doPCcor and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"PC.ps"
    retCode = VLBAPCcor(uv, err, calSou=goodCal["Source"], \
                        timeRange=goodCal["timeRange"], \
                        doCalib=-1, flagVer=2, solInt=manPCsolInt, \
                        PCin=1, SNout=0, refAnt=goodCal["bestRef"], \
                        doPCPlot=doPCPlot, plotFile=plotFile, \
                        noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in PC calibration"

# Plot amplitude and phase vs. frequency
if doSpecPlot:
    plotFile = "./"+project+session+band+".PCcor.spec.ps"
    VLBASpecPlot( uv, goodCal, err, doband=0, plotFile=plotFile, logfile=logFile )

# manual phase cal
if doManPCal and not check:
    plotFile = "./"+project+session+band+".ManPCal.ps"
    retCode = VLBAManPCal(uv, err, calSou=goodCal["Source"], \
                          #CalModel=contCalModel, \
                          timeRange=goodCal["timeRange"], \
                          solInt=manPCsolInt, smoTime=manPCSmoo,  \
                          refAnts=[goodCal["bestRef"]], doCalib=2, flagVer=2, 
                          doManPCalPlot=doManPCalPlot, plotFile=plotFile, noScrat=noScrat, \
                          nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in manual phase calibration"

# Plot amplitude and phase vs. frequency
if doSpecPlot:
    plotFile = "./"+project+session+band+".ManPCal.spec.ps"
    VLBASpecPlot( uv, goodCal, err, doband=0, plotFile=plotFile, logfile=logFile )

# Bandpass calibration if needed
if doBPCal and not check:
    retCode = VLBABPass(uv, goodCal["Source"], err, CalModel=contCalModel, \
                        timeRange=goodCal["timeRange"], doCalib=2, flagVer=2, \
                        noScrat=noScrat, solInt1=bpsolint1, solInt2=bpsolint2, solMode=bpsolMode, \
                        BChan1=bpBChan1, EChan1=bpEChan1, BChan2=bpBChan2, EChan2=bpEChan2, ChWid2=bpChWid2, \
                        doCenter1=bpDoCenter1, refAnt=goodCal["bestRef"], specIndex=specIndex, \
                        doAuto = bpdoAuto, \
                        nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in Bandpass calibration"

# Plot amplitude and phase vs. frequency
if doSpecPlot:
    plotFile = "./"+project+session+band+".spec.ps"
    VLBASpecPlot( uv, goodCal, err, doband=1, plotFile=plotFile, logfile=logFile )

# image cals
if doImgCal and not check:
    retCode = VLBAImageCals(uv, err, Sources=contCals, seq=seq, sclass=outCclass, \
                            doCalib=2, flagVer=2, doBand=1, \
                            FOV=FOV, Robust=Robust, \
                            maxPSCLoop=maxPSCLoop, minFluxPSC=minFluxPSC, solPInt=solPInt, solMode=solMode, \
                            maxASCLoop=maxASCLoop, minFluxASC=minFluxASC, solAInt=solAInt, \
                            avgPol=avgPol, avgIF=avgIF, minSNR=minSNR, refAnt=goodCal["bestRef"], \
                            nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in imaging calibrators"
    
# Check if calibrator models now available
contCalModel = VLBAImageModel(contCals, outCclass, disk, seq, err)

# delay calibration
if doDelayCal and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"DelayCal.ps"
    retCode = VLBADelayCal(uv, err, calSou=contCals, CalModel=contCalModel, \
                               doCalib=2, flagVer=2, doBand=1, \
                               solInt=manPCsolInt, smoTime=manPCSmoo,  \
                               refAnts=[goodCal["bestRef"]], \
                               doSNPlot=doSNPlot, plotFile=plotFile, \
                               nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in delay calibration"
    
# Amplitude calibration
if doAmpCal and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"AmpCal.ps"
    retCode = VLBAAmpCal(uv, err, calSou=contCals, CalModel=contCalModel, \
                         doCalib=2, flagVer=2, doBand=1, \
                         refAnt=goodCal["bestRef"], solInt=manPCsolInt, \
                         smoTimeA=1440.0, smoTimeP = 10.0, \
                         doSNPlot=doSNPlot, plotFile=plotFile, \
                         nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in amplitude calibration"
    
# Calibrate and average continuum calibrator data
if doCalAvg:
    retCode = VLBACalAvg (uv, avgClass, seq, CalAvgTime, err, \
                          flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                          BIF=CABIF, EIF=CAEIF, BChan=CABChan, EChan=CAEChan, \
                          chAvg=chAvg, avgFreq=avgFreq, Compress=Compress, \
                          logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in CalAvg"

# image targets phase only self-cal
if doImgTarget and not check:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    retCode = VLBAImageCals(uv, err, Sources=targets, seq=seq, sclass=outTclass, \
                            doCalib=2, flagVer=2, doBand=1, \
                            FOV=FOV, Robust=Robust, \
                            maxPSCLoop=maxPSCLoop, minFluxPSC=minFluxPSC, solPInt=solPInt, solMode=solMode, \
                            maxASCLoop=0, \
                            avgPol=avgPol, avgIF=avgIF, minSNR=minSNR, refAnt=goodCal["bestRef"], \
                            nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in imaging targets"
    
# Phase calibration using target models
if doPhaseCal:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    # Get list of all if no explicit list given
    if len(targets)<=0:
        slist = VLBAAllSource(uvc,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    targetModel  = VLBAImageModel(slist,outTclass,disk, seq, err)
    plotFile = "./"+project+"_"+session+"_"+band+"PhaseCal.ps"
    retCode = VLBAPhaseCal(uvc, err, calSou=slist, CalModel=targetModel, \
                         doCalib=-1, flagVer=0, doBand=-1, \
                         refAnt=goodCal["bestRef"], solInt=manPCsolInt, \
                         doSNPlot=doSNPlot, plotFile=plotFile, \
                         nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in phase calibration"
    
# Instrumental polarization calibration
if doInstPol:
    # calibrators defaults to strong calibrator list
    if not instPolCal:
        instPolCal = contCals
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    retCode = VLBAPolCal(uvc, instPolCal, err, \
                             doCalib=2, flagVer=0, doBand=-1, doSetJy=True, \
                             refAnt=goodCal["bestRef"], solInt=2.0, \
                             noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in instrumental poln calibration"
    
# RL Phase (EVPA) calibration as BP table
if doRLCal and RLCal:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    #retCode = VLBARLCal(uvc, err, RLPCal=RLCal,  \
    #                        doCalib=2, flagVer=0, doBand=-1, doPol=True, BPSoln=1, \
    #                        refAnt=goodCal["bestRef"], solInt1=1.0, \
    #                        nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    retCode = VLBARLCal2(uvc, err, RLPCal=RLCal, \
                            doCalib=2, gainUse=2, flagVer=0, doBand=-1, doPol=True,  \
                            refAnt=goodCal["bestRef"], niter=300, FOV=0.02/3600.0, \
                            nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in RL phase calibration"
    
# image targets possible with Stokes IQU
if doImgFullTarget:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    retCode = VLBAImageTargets(uvc, err, Sources=targets, seq=seq, sclass=outIclass, \
                  doCalib=2, flagVer=0, doBand=-1, \
                  Stokes=Stokes, FOV=FOV, Robust=Robust, \
                  maxPSCLoop=2, minFluxPSC=minFluxPSC, solPInt=solPInt, solMode=solMode, \
                  maxASCLoop=0, \
                  avgPol=avgPol, avgIF=avgIF, minSNR=minSNR, refAnt=goodCal["bestRef"], \
                  nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in imaging targets"

# Get report on sources
if doReport:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    Report = VLBAReportTargets(uvc, err, Sources=targets, seq=seq, sclass=outIclass, \
                                   Stokes=Stokes, logfile=logFile, check=check, debug=debug)
    # Save to pickle jar
    ReportPicklefile = "./"+project+"_"+session+"_"+band+"Report.pickle"   # Where results saved
    SaveObject(Report, ReportPicklefile, True) 
   
# cleanup
mess ="Write results/Cleanup" 
printMess(mess, logFile)

# Write results, cleanup    
# Save UV data? 
if doSaveUV and (not check):
    # Get calibrated/averaged data
    if not uvc:
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    # Write 
    filename = project+session+band+"CalAvg.uvtab"
    fuv = VLBAUVFITS (uvc, filename, 0, err, compress=Compress)
# Save UV data tables?
if doSaveTab and (not check):
    filename = project+session+band+"CalTab.uvtab"
    fuv = VLBAUVFITSTab (uv, filename, 0, err)
# Imaging results
if doSaveImg:
    # How many Stokes images
    nstok = len(Stokes)
    # Targets
    if len(targets)<=0:
        # fetch full list if needed
        targets = VLBAAllSource(uv, err, logfile=logFile,check=check,debug=debug)
    for target in targets:
        if not check:
            # intermediate images
            oclass = outTclass
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            if (not x.exist): 
                print target,"image not found. Skipping."
                continue
            outfile = project+session+band+target+"."+oclass+".fits"
            mess ="Write Intermediate target " +outfile+" on disk "+str(outDisk)
            printMess(mess, logFile)
            xf = VLBAImFITS (x, outfile, outDisk, err, fract=0.1)
            # Statistics
            zz=imstat(x, err, logfile=logFile)
            del x, xf
            # Final images
            for istok in range(0,nstok):
                oclass = Stokes[istok:istok+1]+outIclass[1:]
                x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
                outfile = project+session+band+target+"."+oclass+".fits"
                mess ="Write " +outfile+" on disk "+str(outDisk)
                printMess(mess, logFile)
                xf = VLBAImFITS (x, outfile, outDisk, err, fract=0.1)
                # Statistics
                zz=imstat(x, err, logfile=logFile)
                del x, xf
    # Calibrators
    for target in contCals:
        if not check:
            oclass = outCclass
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            if (not x.exist): 
                print target,"image not found. Skipping."
                continue
            outfile = project+session+band+target+"."+oclass+".fits"
            mess ="Write Calibrator " +outfile+" on disk "+str(outDisk)
            printMess(mess, logFile)
            xf = VLBAImFITS (x, outfile, outDisk, err, fract=0.1)
            # Statistics
            zz=imstat(x, err, logfile=logFile)
            del x, xf
    # end writing images loop

# Contour plots
if doKntrPlots:
    VLBAKntrPlots( err, project=project, session=session, band=band,
        debug=debug )
elif debug:
    print "Not creating contour plots ( doKntrPlots = ", doKntrPlots, " )"

# Diagnostic plots
if doDiagPlots:
    # Get the highest number avgClass catalog file
    Aname = VLBAAIPSName( project, session )
    uvc = None
    if not check:
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
    VLBADiagPlots( uvc, err, cleanUp=doCleanup, logfile=logFile, check=check, 
        debug=debug )
elif debug:
    print "Not creating diagnostic plots ( doDiagPlots = ", doDiagPlots, " )"

# Cleanup - delete AIPS files
if doCleanup and (not check):
    # Delete target images
    # How many Stokes images
    nstok = len(Stokes)
    for istok in range(0,nstok):
        oclass = Stokes[istok:istok+1]+outIclass[1:]
        AllDest(err, disk=disk,Aseq=seq,Aclass=oclass)
    
    # delete Calibrator images
    AllDest(err, disk=disk,Aseq=seq,Aclass=outCclass)
    
    # Delete intermediate target images
    AllDest(err, disk=disk,Aseq=seq,Aclass=outTclass)
    OErr.printErrMsg(err, "Deleting AIPS images")
    
    # Delete UV data
    uv.Zap(err)
    # Zap calibrated/averaged data
    if not uvc:
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    uvc.Zap(err)
    OErr.printErrMsg(err, "Writing output/cleanup")

# Shutdown
mess = "Finished project "+project
printMess(mess, logFile)
OErr.printErr(err)
OSystem.Shutdown(ObitSys)

