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
# Define data
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
prtLv         = 2          # Print level

# Initialize parameters
parms = VLBAInitContParms()

############################# Set Project Processing parameters ##################
parmFile = sys.argv[2]
execfile (parmFile)

################################## Process #####################################
# Init cal pickle jars
goodCalPicklefile = "./"+project+"_"+session+"_"+band+"GoodCal.pickle"   # Where results saved
# Default "best" calibration
goodCal = {"Source":"  ", "souID":0,"timeRange":(0.0,100.0), "Fract":0.0, "SNR":0.0, "bestRef":0}
SaveObject(goodCal, goodCalPicklefile, False)   # Save initial default
OKCalPicklefile = "./"+project+"_"+session+"_"+band+"OKCal.pickle"   # Where results saved
# Default OK calibrators
SaveObject(parms["contCals"], OKCalPicklefile, False)   # Save initial default

# Load the outputs pickle jar
VLBAFetchOutFiles()

# Logging directly to logFile
OErr.PInit(err, prtLv, logFile)
retCode = 0
VLBAAddOutFile( logFile, 'project', 'Pipeline log file' )

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
    uvname = project+"_"+session+"_"+band
    uv = UV.newPAUV(uvname, Aname, dataClass, disk, seq, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")

# frequency dependent default parameters
VLBAInitContFQParms(uv, parms, err, \
                        logfile=logFile, check=check, debug=debug)

# Clear any old calibration/editing 
if parms["doClearTab"]:
    VLBAClearCal(uv, err, doGain=parms["doGain"], doFlag=parms["doFlag"], doBP=parms["doBP"], \
                     check=check, logfile=logFile)
    OErr.printErrMsg(err, "Error resetting calibration")

# Copy FG 1 to FG 2
if parms["doCopyFG"]:
    retCode = VLBACopyFG (uv, err, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error Copying FG table"

# Special editing
if parms["doEditList"] and not check:
    for edt in parms["editList"]:
        UV.PFlag(uv,err,timeRange=[dhms2day(edt["timer"][0]),dhms2day(edt["timer"][1])], \
                     flagVer=editFG, Ants=edt["Ant"], Chans=edt["Chans"], IFs=edt["IFs"], \
                     Stokes=edt["Stokes"], Reason=edt["Reason"])
        OErr.printErrMsg(err, "Error Flagging")

# Quack to remove data from start and end of each scan
if parms["doQuack"]:
    retCode = VLBAQuack (uv, err, \
                             begDrop=parms["quackBegDrop"], endDrop=parms["quackEndDrop"], \
                             Reason=parms["quackReason"], \
                             logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error Quacking data"

# Quantization correction?
if parms["doQuantCor"]:
    plotFile = "./"+project+"_"+session+"_"+band+"Quant.ps"
    retCode = VLBAQuantCor(uv, parms["QuantSmo"], parms["QuantFlag"], err, \
                               doSNPlot=parms["doSNPlot"], plotFile=plotFile, \
                               logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in quantization correcting/flagging"

# Parallactic angle correction?
if parms["doPACor"]:
    retCode = VLBAPACor(uv, err, noScrat=noScrat, \
                            logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in quantization correcting/flagging"

# Opacity/Tsys/gain correction
if parms["doOpacCor"]:
    plotFile = "./"+project+"_"+session+"_"+band+"Opacity.ps"
    retCode = VLBAOpacCor(uv, parms["OpacSmoo"], err,  \
                              doSNPlot=parms["doSNPlot"], plotFile=plotFile, \
                              logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in opacity/gain.Tsys correction"
    
    VLBASaveOutFiles() # Save plot file in Outfiles

# Need to determine a list of calibrators?
if (parms["contCals"]==None) or (len(parms["contCals"])<=0):
    if parms["doFindOK"]:
        slist = VLBAAllSource(uv, err,logfile=logFile,check=check,debug=debug)
        parms["contCals"] = VLBAOKCal(uv, parms["minOKFract"], err, \
                                          solInt=parms["findSolInt"],  \
                                          calSou=slist, minSNR=parms["minOKSNR"], \
                                          doCalib=-1, flagVer=2, refAnts=parms["refAnts"], \
                                          noScrat=noScrat, nThreads=nThreads, \
                                          logfile=logFile, check=check, debug=debug)
        if not parms["contCals"] and not check:
            raise RuntimeError,"Error in finding acceptable calibrators"
    else:
        # Snatch from pickle jat
        parms["contCals"] = FetchObject(OKCalPicklefile)

# Save contCals to a pickle jar
SaveObject(parms["contCals"], OKCalPicklefile, True)

# Find best calibration source
if parms["doFindCal"]:
    goodCal = VLBAGoodCal(uv,  err, \
                              solInt=parms["findSolInt"], timeInt=parms["findTimeInt"], \
                              calSou=parms["contCals"], \
                              #CalModel=parms["contCalModel"], \
                              doCalib=-1, flagVer=2, refAnts=parms["refAnts"], \
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
if parms["doPCcor"] and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"PC.ps"
    retCode = VLBAPCcor(uv, err, calSou=goodCal["Source"], \
                        timeRange=goodCal["timeRange"], \
                        doCalib=-1, flagVer=2, solInt=parms["manPCsolInt"], \
                        PCin=1, SNout=0, refAnt=goodCal["bestRef"], \
                        doPCPlot=parms["doPCPlot"], plotFile=plotFile, \
                        noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in PC calibration"
    VLBASaveOutFiles() # Save plot file in Outfiles

# manual phase cal
if parms["doManPCal"] and not check:
    plotFile = "./"+project+session+band+".ManPCal.ps"
    retCode = VLBAManPCal(uv, err, calSou=goodCal["Source"], \
                              #CalModel=parms["contCalModel"], \
                              timeRange=goodCal["timeRange"], \
                              solInt=parms["manPCsolInt"], smoTime=parms["manPCSmoo"],  \
                              refAnts=[goodCal["bestRef"]], doCalib=2, flagVer=2, \
                              doManPCalPlot=parms["doManPCalPlot"], \
                              plotFile=plotFile, noScrat=noScrat, \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in manual phase calibration"

# Bandpass calibration if needed
if parms["doBPCal"] and not check:
    retCode = VLBABPass(uv, goodCal["Source"], err, CalModel=None, \
                            timeRange=goodCal["timeRange"], doCalib=2, flagVer=2, \
                            noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], BChan2=parms["bpBChan2"], \
                            EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=goodCal["bestRef"], specIndex=parms["specIndex"], \
                            doAuto = parms["bpdoAuto"], \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in Bandpass calibration"

# Plot amplitude and phase vs. frequency
if parms["doSpecPlot"]:
    plotFile = "./"+project+session+band+".spec.ps"
    VLBASpecPlot( uv, goodCal, err, doband=1, check=check, plotFile=plotFile, logfile=logFile )
    VLBASaveOutFiles() # Save plot file in Outfiles

# image cals
if parms["doImgCal"] and not check:
    retCode = VLBAImageCals(uv, err, Sources=parms["contCals"], seq=seq, sclass=parms["outCclass"], \
                                doCalib=2, flagVer=2, doBand=1, \
                                FOV=parms["FOV"], Robust=parms["Robust"], \
                                maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], \
                                solPInt=parms["solPInt"], solMode=parms["solMode"], \
                                maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], solAInt=parms["solAInt"], \
                                avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR=parms["minSNR"], refAnt=goodCal["bestRef"], \
                                nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in imaging calibrators"
    
# Check if calibrator models now available
parms["contCalModel"] = VLBAImageModel(parms["contCals"], parms["outCclass"], disk, seq, err)

# delay calibration
if parms["doDelayCal"] and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"DelayCal.ps"
    retCode = VLBADelayCal(uv, err, calSou=parms["contCals"], CalModel=parms["contCalModel"], \
                               doCalib=2, flagVer=2, doBand=1, \
                               solInt=parms["manPCsolInt"], smoTime=parms["delaySmoo"],  \
                               refAnts=[goodCal["bestRef"]], \
                               doSNPlot=parms["doSNPlot"], plotFile=plotFile, \
                               nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in delay calibration"
    VLBASaveOutFiles() # Save plot file in Outfiles
    
# Amplitude calibration
# NOTE: REALLY OUGHT TO HAVE MODEL FOR THIS
if parms["doAmpCal"] and not check:
    plotFile = "./"+project+"_"+session+"_"+band+"AmpCal.ps"
    retCode = VLBAAmpCal(uv, err, calSou=parms["contCals"], CalModel=parms["contCalModel"], \
                         doCalib=2, flagVer=2, doBand=1, \
                         refAnt=goodCal["bestRef"], solInt=parms["solAInt"], \
                         smoTimeA=24.0, smoTimeP = 10./60., \
                         doSNPlot=parms["doSNPlot"], plotFile=plotFile, \
                         nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in amplitude calibration"
    VLBASaveOutFiles() # Save plot file in Outfiles
    
# Calibrate and average  data
if parms["doCalAvg"]:
    retCode = VLBACalAvg (uv, parms["avgClass"], seq, parms["CalAvgTime"], err, \
                              flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                              BIF=parms["CABIF"], EIF=parms["CAEIF"], \
                              BChan=parms["CABChan"], EChan=parms["CAEChan"], \
                              chAvg=parms["chAvg"], avgFreq=parms["avgFreq"], Compress=Compress, \
                              logfile=logFile, check=check, debug=debug)
    if retCode!=0:
       raise  RuntimeError,"Error in CalAvg"

# image targets phase only self-cal
if parms["doImgTarget"] and not check:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    retCode = VLBAImageCals(uv, err, Sources=parms["targets"], seq=seq, sclass=parms["outTclass"], \
                                doCalib=2, flagVer=2, doBand=1, \
                                FOV=parms["FOV"], Robust=parms["Robust"], \
                                maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], \
                                solPInt=parms["solPInt"], solMode=parms["solMode"], \
                                maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], solAInt=parms["solAInt"], \
                                avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR=parms["minSNR"],\
                                refAnt=goodCal["bestRef"], \
                                nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in imaging targets"
    
# Phase calibration using target models
if parms["doPhaseCal"]:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    # Get list of all if no explicit list given
    if len(parms["targets"])<=0:
        slist = VLBAAllSource(uvc,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = parms["targets"]
    targetModel  = VLBAImageModel(slist,parms["outTclass"],disk, seq, err)
    plotFile = "./"+project+"_"+session+"_"+band+"PhaseCal.ps"
    retCode = VLBAPhaseCal(uvc, err, calSou=slist, CalModel=parms["targetModel"], \
                         doCalib=-1, flagVer=0, doBand=-1, \
                         refAnt=goodCal["bestRef"], solInt=parms["manPCsolInt"], \
                         doSNPlot=parms["doSNPlot"], plotFile=plotFile, \
                         nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in phase calibration"
    VLBASaveOutFiles() # Save plot file in Outfiles
    
# Instrumental polarization calibration
if parms["doInstPol"]:
    # calibrators defaults to strong calibrator list
    if not parms["instPolCal"]:
        instPolCal = contCals
    else:
        instPolCal = parms["instPolCal"]
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    retCode = VLBAPolCal(uvc, instPolCal, err, \
                             doCalib=2, flagVer=0, doBand=-1, doSetJy=True, \
                             refAnt=goodCal["bestRef"], solInt=2.0, \
                             noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in instrumental poln calibration"
    
# RL Phase (EVPA) calibration as BP table
if parms["doRLCal"] and parms["RLCal"]:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    retCode = VLBARLCal2(uvc, err, RLPCal=parms["RLCal"], \
                            doCalib=2, gainUse=2, flagVer=0, doBand=-1, doPol=True,  \
                            refAnt=goodCal["bestRef"], niter=300, FOV=0.02/3600.0, \
                            nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in RL phase calibration"
    
# image targets possible with Stokes I(QU)
if parms["doImgFullTarget"]:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    # Get list of all if no explicit list given
    if len(parms["targets"])<=0:
        slist = VLBAAllSource(uvc,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = parms["targets"]
    retCode = VLBAImageTargets(uvc, err, Sources=slist, seq=seq, sclass=parms["outIclass"], \
                                   doCalib=2, flagVer=0, doBand=-1, \
                                   Stokes=parms["Stokes"], FOV=parms["FOV"], Robust=parms["Robust"], \
                                   maxPSCLoop=2, minFluxPSC=parms["minFluxPSC"], solPInt=parms["solPInt"], solMode="P", \
                                   maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], solAInt=parms["solAInt"], \
                                   avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR=parms["minSNR"], refAnt=goodCal["bestRef"], \
                                   nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError,"Error in imaging targets"

# Save UV data? 
if parms["doSaveUV"] and (not check):
    mess ="Write calibrated and averaged UV data to disk" 
    printMess(mess, logFile)
    # Get calibrated/averaged data
    if not uvc:
        Aname = VLBAAIPSName(project, session)
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    # Write 
    filename = project+session+band+"CalAvg.uvtab"
    fuv = VLBAUVFITS (uvc, filename, 0, err, compress=Compress)
    VLBAAddOutFile( filename, 'project', "Calibrated Averaged UV data" )
    # Save list of output files
    VLBASaveOutFiles()

# Save UV data tables?
if parms["doSaveTab"] and (not check):
    filename = project+session+band+"CalTab.uvtab"
    fuv = VLBAUVFITSTab (uv, filename, 0, err)
    VLBAAddOutFile( filename, 'project', "Calibrated AIPS tables" )
    # Save list of output files
    VLBASaveOutFiles()
# Imaging results
outDisk = 0
if parms["doSaveImg"]:
    # How many Stokes images
    nstok = len(parms["Stokes"])
    # Targets
    if len(parms["targets"])<=0:
        # fetch full list if needed
        targets = VLBAAllSource(uv, err, logfile=logFile,check=check,debug=debug)
    for target in targets:
        if not check:
            # intermediate images
            oclass = parms["outTclass"]
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            if (not x.exist): 
                print target,"image not found. Skipping."
                continue
            outfile = project+session+band+target+"."+oclass+".fits"
            mess ="Write Intermediate target " +outfile+" on disk "+str(outDisk)
            VLBAAddOutFile( outfile, target, 'Intermediate target image' )
            printMess(mess, logFile)
            xf = VLBAImFITS (x, outfile, outDisk, err, fract=0.1)
            # Save list of output files
            VLBASaveOutFiles()
            # Statistics
            zz=imstat(x, err, logfile=logFile)
            del x, xf
            # Final images
            for istok in range(0,nstok):
                oclass = parms["Stokes"][istok:istok+1]+parms["outIclass"][1:]
                x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
                outfile = project+session+band+target+"."+oclass+".fits"
                mess ="Write " +outfile+" on disk "+str(outDisk)
                printMess(mess, logFile)
                xf = VLBAImFITS (x, outfile, outDisk, err, fract=0.1)
                VLBAAddOutFile( outfile, target, 'Image' )
                print "Writing file ", outfile                
                # Statistics
                zz=imstat(x, err, logfile=logFile)
                del x, xf
                # Save list of output files
                VLBASaveOutFiles()
    # Calibrators
    for target in parms["contCals"]:
        if not check:
            oclass = parms["outCclass"]
            x = Image.newPAImage("out", target, oclass, disk, seq, True, err)
            if (not x.exist): 
                print target,"image not found. Skipping."
                continue
            outfile = project+session+band+target+"."+oclass+".fits"
            mess ="Write Calibrator " +outfile+" on disk "+str(outDisk)
            printMess(mess, logFile)
            xf = VLBAImFITS (x, outfile, outDisk, err, fract=0.1)
            VLBAAddOutFile( outfile, target, 'Calibrator Image' )
            # Statistics
            zz=imstat(x, err, logfile=logFile)
            del x, xf
            # Save list of output files
            VLBASaveOutFiles()
    # end writing images loop

# Contour plots
if parms["doKntrPlots"]:
    VLBAKntrPlots( err, project=project, session=session, band=band,
        disk=disk, debug=debug )
    # Save list of output files
    VLBASaveOutFiles()
elif debug:
    print "Not creating contour plots ( doKntrPlots = ", parms["doKntrPlots"], " )"

# Source uv plane diagnostic plots
if parms["doDiagPlots"]:
    # Get the highest number avgClass catalog file
    Aname = VLBAAIPSName( project, session )
    uvc = None
    if not check:
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
    VLBADiagPlots( uvc, err, cleanUp=parms["doCleanup"], \
                       project=project, session=session, band=band, \
                       logfile=logFile, check=check, debug=debug )
    # Save list of output files
    VLBASaveOutFiles()
elif debug:
    print "Not creating diagnostic plots ( doDiagPlots = ", parms["doDiagPlots"], " )"

# Save metadata
srcMetadata = None
projMetadata = None
if parms["doMetadata"]:
    if not uvc:
        # Get calibrated/averaged data
        Aname = VLBAAIPSName(project, session)
        uvname = project+"_"+session+"_"+band+"_Cal"
        uvc = UV.newPAUV(uvname, Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")

    # Get source metadata; save to pickle file
    srcMetadata = VLBASrcMetadata( uvc, err, Sources=parms["targets"], seq=seq, 
        sclass=parms["outIclass"], Stokes=parms["Stokes"], logfile=logFile, check=check, 
        debug=debug )
    picklefile = "./"+project+"_"+session+"_"+band+"SrcReport.pickle" 
    SaveObject(srcMetadata, picklefile, True) 

    # Get project metadata; save to pickle file
    projMetadata = VLBAProjMetadata( uvc, AIPS_VERSION, err, contCals=parms["contCals"],
        goodCal = goodCal, project = project, session = session, band = band )
    picklefile = "./"+project+"_"+session+"_"+band+"ProjReport.pickle"
    SaveObject(srcMetadata, picklefile, True) 

# Write report
if parms["doHTML"]:
    VLBAHTMLReport( projMetadata, srcMetadata, \
                        outfile=project+"_"+session+"_"+band+"report.html", \
                        logFile=logFile )

# Save list of output files
VLBASaveOutFiles()

# Copy output files to specificed destination directory
if parms["copyDestDir"]:
    VLBACopyOutFiles( destDir=parms["copyDestDir"], logFile=logFile )

# Cleanup - delete AIPS files
if parms["doCleanup"] and (not check):
    # Delete target images
    # How many Stokes images
    nstok = len(parms["Stokes"])
    for istok in range(0,nstok):
        oclass = parms["Stokes"][istok:istok+1]+parms["outIclass"][1:]
        AllDest(err, disk=disk,Aseq=seq,Aclass=oclass)
    
    # delete Calibrator images
    AllDest(err, disk=disk,Aseq=seq,Aclass=parms["outCclass"])
    
    # Delete intermediate target images
    AllDest(err, disk=disk,Aseq=seq,Aclass=parms["outTclass"])
    OErr.printErrMsg(err, "Deleting AIPS images")
    
    # Delete UV data
    uv.Zap(err)
    # Zap calibrated/averaged data
    if not uvc:
        Aname = VLBAAIPSName(project, session)
        uvc = UV.newPAUV("AIPS UV DATA", Aname, parms["avgClass"], disk, seq, True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    uvc.Zap(err)
    OErr.printErrMsg(err, "Writing output/cleanup")

# Shutdown
mess = "Finished project "+project
printMess(mess, logFile)
OErr.printErr(err)
OSystem.Shutdown(ObitSys)

