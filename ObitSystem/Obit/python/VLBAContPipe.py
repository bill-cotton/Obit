# VLBA Continuum Pipeline processing 
# AIPS/FITS setup and Parameter file given as arguments:
# > python VLBAContPipe.py AIPSSetup.py parms.py
#
import sys
import OErr, OSystem, UV, AIPS, FITS
import ObitTalkUtil
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

# Find good calibration data
doFindCal    = True         # Search for good calibration/reference antenna
findSolInt   = 0.5          # Solution interval (min) for Calib
findTimeInt  = 10.0         # Maximum timerange, large=>scan
contCals     = None         # Name or list of continuum cals
contCalModel = VLBAImageModel(contCals,"CalMod", disk, seq, err)         # Check for any

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
refAnts       = [0]                     # List of Reference antenna for fringe fitting

# Imaging calibrators (contCals)
doImgCal    = True         # Image calibrators
targets     = []           # targets
doImgTarget = True         # Image targets?
outCclass   = "ICalSC"     # Output calibrator image class
outIclass   = "IImgSC"     # Output target image class
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
    mess = "Using Debug mode "
    printMess(mess, logFile)
if check:
    mess = "Only checking script"
    printMess(mess, logFile)

# Load Data from FITS
uv = None
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
                        PCin=1, SNout=0, refAnt=goodCal["bestRef"], \
                        doPCPlot=doPCPlot, plotFile=plotFile, \
                        logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in PC calibration"

# manual phase cal
if doManPCal and not check:
    retCode = VLBAManPCal(uv, err, calSou=goodCal["Source"], \
                          #CalModel=contCalModel, \
                          timeRange=goodCal["timeRange"], \
                          solInt=manPCsolInt, smoTime=manPCSmoo,  \
                          refAnts=[goodCal["bestRef"]], doCalib=2, flagVer=2, noScrat=noScrat, \
                          nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in manual phase calibration"

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

# image targets
if doImgTarget and not check:
    retCode = VLBAImageCals(uv, err, Sources=targets, seq=seq, sclass=outIclass, \
                            doCalib=2, flagVer=2, doBand=1, \
                            FOV=FOV, Robust=Robust, \
                            maxPSCLoop=maxPSCLoop, minFluxPSC=minFluxPSC, solPInt=solPInt, solMode=solMode, \
                            maxASCLoop=maxASCLoop, minFluxASC=minFluxASC, solAInt=solAInt, \
                            avgPol=avgPol, avgIF=avgIF, minSNR=minSNR, refAnt=goodCal["bestRef"], \
                            nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in imaging targets"
    
# cleanup
mess ="Write results/Cleanup" 
printMess(mess, logFile)

# Write results, cleanup    
# Save UV data? 
if doSaveUV and (not check):
    # Get calibrated/averaged data
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
    # Targets
    if len(targets)<=0:
        # fetch full list if needed
        targets = VLBAAllSource(uv, err, logfile=logFile,check=check,debug=debug)
    for target in targets:
        if not check:
            oclass = outIclass
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            outfile = project+session+band+target+"."+outIclass+".fits"
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
            outfile = project+session+band+target+"."+outIclass+".fits"
            mess ="Write " +outfile+" on disk "+str(outDisk)
            printMess(mess, logFile)
            xf = VLBAImFITS (x, outfile, outDisk, err, fract=0.1)
            # Statistics
            zz=imstat(x, err, logfile=logFile)
            del x, xf
    # end writing images loop

# Diagnostic plots
if doDiagPlots:
    # Get the highest number avgClass catalog file
    Aname = VLBAAIPSName( project, session )
    uvc = None
    if not check:
        uvc = UV.newPAUV("AIPS UV DATA", Aname, avgClass, disk, seq, True, err)
    VLBADiagPlots( uvc, err, logfile=logFile, check=check, debug=debug )

# Cleanup - delete AIPS files
if doCleanup and (not check):
    # Delete target images
    if len(targets)<=0:
        # fetch full list if needed
        targets = VLBAAllSource(uv, err, logfile=logFile,check=check,debug=debug)
    for target in targets:
        oclass = outIclass
        x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
        if Image.PIsA(x):
            x.Zap(err) # cleanup
        del x
    # delete Calibrator images
    for target in contCals:
        oclass = outCclass
        x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
        if Image.PIsA(x):
            x.Zap(err) # cleanup
        del x
    # Delete UV data
    uv.Zap(err)
    # Zap calibrated/averaged data
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

