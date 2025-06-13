##########################################################################
# Pipeline processing for calibrating and imaging MeerKAT data
# AIPS/FITS setup and Parameter file given as arguments:
# > ObitTalk MeerKATPipeline.py AIPSSetup.py parms.py
from __future__ import absolute_import
import sys, pydoc, math
import OErr, OSystem, UV, AIPS, FITS, OTObit
import ObitTalkUtil, AIPSTask, ObitTask
from AIPS import AIPSDisk
from FITS import FITSDisk
from PipeUtil import *
from MeerKATCal import *
import OPlot, Table, FArray
############################# Initialize OBIT ##########################################                                 
setup = sys.argv[1]
noScrat     = []    
#exec(compile(open(setup).read(), setup, 'exec'))
exec(open(setup).read())

############################# Set Project Processing parameters ########################
####### Initialize parameters dictionary ##### 
parms = MKInitContParms()
############################# Set Project Processing parameters ##################
parmFile = sys.argv[2]
#exec(compile(open(parmFile).read(), parmFile, 'exec'))
exec(open(parmFile).read())
MKAddOutFile(parmFile, 'project', 'Pipeline input parameters')

################################## Process #####################################
#logFile       = project+"_"+session+"_"+band+".log"  # Processing log file
# Logging directly to logFile
prtLv = 2
OErr.PInit(err, prtLv, logFile)
OSystem.PAllowThreads(nThreads)   # Allow threads in Obit/python
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


# General AIPS data parameters at script level
parms["data_class"]  = ("UVDa"+band)[0:6] # AIPS class of raw uv data
parms["delay_class"] = ("DELA"+band)[0:6] # AIPS class of DelayCal data
outIClass = parms["outIClass"] # image AIPS class
# Set uv  class, sequence number
data_class  = parms["data_class"] 
delay_class = parms["delay_class"] 
data_seq  = 1 
delay_seq = 1 

####################### Import data into AIPS from uvtab file via Hann  ##########################
# Is Hanned data there already?
exists = UV.AExist(MKAIPSName(project), data_class, disk, data_seq, err)
# Loading is done via Hann
if exists:
    uv = UV.newPAUV("AIPS UV DATA", MKAIPSName(project), data_class, disk, data_seq, True, err)
# Extract metadata from archive uvtab version
uf   = UV.newPFUV("Raw", parms["DataFile"], 0, True, err)
meta = MKGetMeta(uf, parms, logFile, err)

###################### Frequency dependent parameters #############################################
parms = MKInitContFQParms(meta,parms)
###################### Setup for polarization calibration #########################################
if parms["doPol"]:
    parms = MKPolSetup(parms, meta, logFile, err)
# Number of threads from AIPSSetup
parms["nThreads"] = nThreads
# Load the outputs pickle jar
MKFetchOutFiles()

retCode = 0
doBand = -1
BPVer = 0
maxgap = max(parms["CalAvgTime"], 160.*8.)/60. # for indexing
fileRoot =  './'+project+"_"+session+"_"+band  # root of file names
uvc           = None
avgClass      = ("UVAv"+band)[0:6]  # Averaged data AIPS class
outIClass     = parms["outIClass"]  # image AIPS class
    
# Log parameters
printMess("Parameter settings", logFile)
for p in parms:
    mess = "  "+p+": "+str(parms[p])
    printMess(mess, logFile)
clist = []
for DCal in parms["DCals"]:
    if DCal["Source"] not in clist:
        clist.append(DCal["Source"])
for PCal in parms["PCals"]:
    if PCal["Source"] not in clist:
        clist.append(PCal["Source"])
for ACal in parms["ACals"]:
    if ACal["Source"] not in clist:
        clist.append(ACal["Source"])
targets = parms['targets']

refAnt = parms['refAnt']

# Save parameters to pickle jar, manifest
ParmsPicklefile = fileRoot+".Parms.pickle"   # Where results saved
SaveObject(parms, ParmsPicklefile, True)
MKAddOutFile(os.path.basename(ParmsPicklefile), 'project', 'Processing parameters used' )

# Initial Checks
# If doPol, XYCals must be fully specified
if parms["doPol"]:
    badCal = not parms["XYCals"];
    if not badCal:
        for cal in parms["XYCals"]:
            badCal = badCal or (not cal[1]) or (not cal[2])
    if (badCal):
        mess = "Pol Cal wanted and XYDCal specified as "+str(parms["XYDCal"])
        printMess(mess, logFile)
        raise  RuntimeError("XY Delay calibrator(s) not fully specified")
    # Check that DCals and BPCals are in unpolarized list
    for c in parms["DCals"]:
        if not MKCheckUnpol(c["Source"]):
            raise  RuntimeError("Delay calibrator "+c["Source"]+" not known to be unpolarized")
    for c in parms["BPCals"]:
        if not MKCheckUnpol(c["Source"]):
            raise  RuntimeError("Bandpass calibrator "+c["Source"]+" not known to be unpolarized")
# end of doPol checks

# Hanning to load data if it doesn't already exist
doNDCal = len(parms["DCalFile"])>0  # Name given
if parms["doHann"]:
    uf = UV.newPFUV("Raw", parms["DataFile"], 0, True, err)
    # Create FG 1 for static flagging 
    MKStaticFlag(uf, 1, err)
    # Use Hann to trim outer channels leaving a multiple of 8
    # Tom's algorithm with update for effects of Hanning
    nchan = uf.Desc.Dict['inaxes'][uf.Desc.Dict['jlocf']]
    first_chan = int(nchan*parms["begChanFrac"])
    last_chan=nchan-int(nchan*parms["endChanFrac"])
    chan_after_ifs=(last_chan-first_chan +1)%16
    first_chan=first_chan+(chan_after_ifs//2)
    last_chan=last_chan-(chan_after_ifs-(chan_after_ifs//2))
    BChan = first_chan-1;  EChan = last_chan+1
    uv = MKHann(uf, MKAIPSName(project), data_class, disk, data_seq, err, \
                doDescm=parms["doDescm"], flagVer=1, BChan=BChan, EChan=EChan,\
                logfile=logFile, zapin=False, check=check, debug=debug)
    # Break into 8 IFs for calibration
    mess = "Divide into 8 IFs"; printMess(mess, logFile)
    import MakeIFs
    MakeIFs.UVMakeIF(uv,8,err)
    # Now Hanning the DelayCal data - do we have it
    doNDCal = len(parms["DCalFile"])>0  # Name given
    # Only if doPol and it doesn't already exist
    if parms["doPol"] and doNDCal:
        mess = "Hanning delay calibration scan"
        delay_uf = UV.newPFUV("Raw", parms["DCalFile"], 0, True, err)
        # Create FG 1 for static flagging 
        MKStaticFlag(delay_uf, 1, err)
        printMess(mess, logFile)
        delay_uv = MKHann(delay_uf, MKAIPSName(project), delay_class, disk, \
                          delay_seq, err, doDescm=parms["doDescm"], flagVer=1, 
                          BChan=BChan, EChan=EChan, logfile=logFile, zapin=False, 
                          check=check, debug=debug)
        parms["delay_uv"] = delay_uv
        MakeIFs.UVMakeIF(delay_uv,8,err)  # Break into 8 IFs
        
    if uv==None and not check:
        raise RuntimeError("Cannot Hann data ")

# Print the uv data header to screen.
uv = UV.newPAUV("AIPS UV DATA", MKAIPSName(project),data_class,disk,data_seq,True,err)
uv.Header(err)
OErr.printErrMsg(err, "Error Finding AIPS Data")

# Clear any old calibration/editing 
if parms["doClearTab"]:
    mess =  "Clear previous calibration"
    printMess(mess, logFile)
    MKClearCal(uv, err, doGain=parms["doClearGain"], doFlag=parms["doClearFlag"], \
               doBP=parms["doClearBP"], check=check)
    OErr.printErrMsg(err, "Error resetting calibration")

# Copy FG 1 to FG 2
if parms["doCopyFG"]:
    mess =  "Copy FG 1 to FG 2"
    printMess(mess, logFile)
    retCode = MKCopyFG(uv, err, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError("Error Copying FG table")

# Special editing
if parms["doEditList"] and not check:
    mess =  "Special editing"
    printMess(mess, logFile)
    for edt in parms["editList"]:
        print ("Debug",edt)
        UV.PFlag(uv,err,timeRange=[dhms2day(edt["timer"][0]),dhms2day(edt["timer"][1])], \
                     flagVer=parms["editFG"], Ants=edt["Ant"], Chans=edt["Chans"], IFs=edt["IFs"], \
                     Stokes=edt["Stokes"], Reason=edt["Reason"])
        OErr.printErrMsg(err, "Error Flagging")

# Flag antennas shadowed by others?
if parms["doShadow"]:
    retCode = MKShadow (uv, err, shadBl=parms["shadBl"], \
                        logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError("Error Shadow flagging data")

# Median window time editing, for RFI impulsive in time
if parms["doMedn"]:
    mess =  "Median window time editing, for RFI impulsive in time:"
    printMess(mess, logFile)
    retCode = MKMedianFlag (uv, clist, err, noScrat=noScrat, nThreads=nThreads, \
                            avgTime=parms["mednAvgTime"], avgFreq=parms["mednAvgFreq"],  \
                            chAvg= parms["mednChAvg"], \
                            timeWind=parms["mednTimeWind"],flagVer=2, flagTab=2, \
                            flagSig=parms["mednSigma"], \
                            logfile=logFile, check=check, debug=False)
    if retCode!=0:
        raise RuntimeError("Error in MednFlag")

# Median window frequency editing, for RFI impulsive in frequency
if parms["doFD1"]:
    mess =  "Median window frequency editing, for RFI impulsive in frequency:"
    printMess(mess, logFile)
    retCode = MKAutoFlag (uv, clist, err, flagVer=2, flagTab=2, doCalib=-1, doBand=-1,   \
                          timeAvg=parms["FD1TimeAvg"], \
                          doFD=True, FDmaxAmp=1.0e20, FDmaxV=1.0e20, FDwidMW=parms["FD1widMW"],  \
                          FDmaxRMS=[1.0e20,0.1], FDmaxRes=parms["FD1maxRes"],  \
                          FDmaxResBL= parms["FD1maxRes"],  FDbaseSel=parms["FD1baseSel"],\
                          nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise  RuntimeError("Error in AutoFlag")

# RMS/Mean editing for calibrators
if parms["doRMSAvg"]:
    mess =  "RMS/Mean editing for calibrators:"
    printMess(mess, logFile)
    clist = []   # Calibrator list
    for s in parms["ACals"]:
        if s['Source'] not in clist:
            clist.append(s['Source'])
    for s in parms["PCals"]:
        if s['Source'] not in clist:
            clist.append(s['Source'])
    for s in parms["DCals"]:
        if s['Source'] not in clist:
            clist.append(s['Source'])
    retCode = MKAutoFlag (uv, clist, err,  flagVer=2, doCalib=-1, doBand=-1,   \
                          RMSAvg=parms["RMSAvg"], timeAvg=parms["RMSTimeAvg"], \
                          nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise  RuntimeError("Error in AutoFlag")
    
# Need to find a reference antenna?  See if we have saved it?
if (not parms["refAnt"]):
    refAnt = FetchObject(fileRoot+".refAnt.pickle")
    if refAnt:
        parms["refAnt"] = refAnt

# Determine refAnt - Use bandpass calibrator and center half of each spectrum
if not parms["refAnt"]:
    mess = "Find best reference antenna: run Calib on BP Cal(s) "
    printMess(mess, logFile)
    parms["refAnt"] = MKGetRefAnt(uv, parms["BPCals"], err, flagVer=0, \
                                  solInt=parms["bpsolint1"], nThreads=nThreads, \
                                  logfile=logFile, check=check, debug=debug)
    if err.isErr:
        raise  RuntimeError("Error finding reference antenna")
    if parms["refAnts"][0]<=0:
        parms["refAnts"][0] = parms["refAnt"]
    mess = "Picked reference antenna "+str(parms["refAnt"])
    printMess(mess, logFile)
    # Save it
    #ParmsPicklefile = fileRoot+".Parms.pickle"   # Where results saved
    #FAILS SaveObject(parms, ParmsPicklefile, True)
    refAntPicklefile = fileRoot+".refAnt.pickle"   # Where results saved
    SaveObject(parms["refAnt"], refAntPicklefile, True)


# Plot Raw, edited data?
if parms["doRawSpecPlot"] and parms["plotSource"]:
    mess =  "Raw Spectral plot for: "+' '.join(parms["BPCal"])
    printMess(mess, logFile)
    plotFile = fileRoot+"_RawSpec.ps"
    retCode = MKSpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                         Stokes=["RR","LL"], doband=-1, flagVer=2,   \
                         check=check, debug=debug, logfile=logFile )
    if retCode!=0:
        raise  RuntimeError("Error in Plotting spectrum")
    MKAddOutFile(plotFile, 'project', 'Pipeline log file' )

if parms["doNDCal"] and parms["doPol"] and doNDCal:
    mess = "XYphase bandpass noise diode calibration"
    printMess(mess, logFile)
    # Delay cal data
    delay_uv = UV.newPAUV("AIPS UV DATA", MKAIPSName(project),delay_class,disk,delay_seq,True,err)
    retCode = MKXPhase(delay_uv, uv, err, logfile=logFile, check=check, debug=debug,
                       doCalib=-1, flagVer=0, doBand=-1, BPSoln=1, refAnt=parms['refAnt'])
    doBand = 1
    BPVer += 1
    if retCode!=0:
        raise RuntimeError("Error in XY phase calibration")
    # Plot table?
    if parms["doBPPlot"] and not check:
        # PLPlot doesn't like long names
        plotFile = project+"_XYPhaseBP.ps"
        MKPlotXYBPTab(uv, 1, plotFile, err, \
                      logfile=logFile, check=check, debug=debug)

# delay calibration
if parms["doDelayCal"] and parms["DCals"] and not check:
    plotFile = fileRoot+"_DelayCal.ps"
    retCode = MKDelayCal(uv, parms["DCals"], err,  \
                         BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                         doCalib=-1, flagVer=0, doBand=doBand, BPVer=BPVer, \
                         solInt=parms["delaySolInt"], smoTime=parms["delaySmoo"],  \
                         refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], 
                         doZeroPhs=parms["delayZeroPhs"], \
                         doAvgIF=parms["delayAvgIF"], doAvgPol=parms["delayAvgPol"], \
                         doPlot=parms["doSNPlot"], plotFile=plotFile, \
                         nThreads=nThreads, noScrat=noScrat, \
                         logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError("Error in delay calibration")
    
    # Plot corrected data?  No - turn off
    if parms["doSpecPlot"] and parms["plotSource"] and not True:
        plotFile = fileRoot+"_DelaySpec.ps"
        retCode = MKSpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, \
                             plotFile, parms["refAnt"], err, flagVer=2, \
                             Stokes=["RR","LL"], doband=doBand,          \
                             check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError("Error in Plotting spectrum")

# Bandpass calibration
if parms["doBPCal"] and parms["BPCals"]:
    retCode = MKBPCal(uv, parms["BPCals"], err, doBand=doBand, BPVer=BPVer, newBPVer=0,
                      noScrat=noScrat, solInt1=parms["bpsolint1"], \
                      solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                      BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                      BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                      doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                      UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=0, doPlot=False, \
                      nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError("Error in Bandpass calibration")
    # Plot corrected data? 
    if parms["doSpecPlot"] and  parms["plotSource"]:
        plotFile = fileRoot+"_BPSpec.ps"
        retCode = MKSpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, \
                             parms["refAnt"], err, Stokes=["RR","LL"], doband=1, flagVer=2, \
                             check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError("Error in Plotting spectrum")
 
# Amp & phase Calibrate
if parms["doAmpPhaseCal"]:
    plotFile = fileRoot+"_APCal.ps"
    retCode = MKCalAP (uv, [], parms["ACals"], parms["GCalList"], err, PCals=parms["PCals"],\
                       doCalib=2, doBand=1, BPVer=0, flagVer=2, \
                       BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                       solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                       doAmpEdit=parms["doAmpEdit"], ampSigma=parms["ampSigma"], \
                       ampEditFG=parms["ampEditFG"], avgPol=parms["doPol"], \
                       doPlot=parms["doSNPlot"], plotFile=plotFile,  refAnt=parms["refAnt"], \
                       nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
    
    if retCode!=0:
        raise RuntimeError("Error calibrating")
    # Plot corrected data?  NO, turn off
    if parms["doSpecPlot"] and  parms["plotSource"] and not True:
        plotFile = fileRoot+"_APSpec.ps"
        retCode = MKSpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, \
                             parms["refAnt"], err, Stokes=["RR","LL"], doband=1, flagVer=2, \
                             check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError("Error in Plotting spectrum")
        

# More editing
if parms["doAutoFlag"]:
    mess =  "Post calibration editing:"
    printMess(mess, logFile)
    # if going to redo then only calibrators
    if parms["doRecal"]:
        # Only calibrators
        clist = []
        for DCal in parms["DCals"]:
            if DCal["Source"] not in clist:
                clist.append(DCal["Source"])
        for PCal in parms["PCals"]:
            if PCal["Source"] not in clist:
                clist.append(PCal["Source"])
        for ACal in parms["ACals"]:
            if ACal["Source"] not in clist:
                clist.append(ACal["Source"])
    else:
        clist = []
        
        retCode = MKAutoFlag (uv, clist, err, flagVer=2, flagTab =2, \
                              doCalib=2, gainUse=0, doBand=1, BPVer=BPVer,  \
                              minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                              doFD=parms["doFirstAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                              FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                              FDmaxRes=parms["FDmaxRes"],  FDmaxResBL=parms["FDmaxResBL"], \
                              FDbaseSel=parms["FDbaseSel"], \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError("Error in AutoFlag")
    
# Calibrator Clipping if redoing the calibration
if parms["doClipCals"] and parms["doRecal"]:
    mess =  "Calibrator Clipping:"
    printMess(mess, logFile)
    # all calibrators should be in the PCals list
    clist = []
    for PCal in parms["PCals"]:
        if PCal["Source"] not in clist:
            clist.append(PCal["Source"])
    # Any polarized calibrators
    if parms["doXYDelay"]:
        for cal in parms["XYCals"]:
            if cal[0] not in clist:
                clist.append(cal[0])
    if meta['nstokes']>=4: # Full polarization?
        XClip = parms["XClip"]
    else:
        XClip = [0.0,0.0]
    mess =  "Clipping:"+str(clist)
    printMess(mess, logFile)
    # Clip
    retCode = MKAutoFlag (uv, clist, err, flagVer=2, flagTab =2, \
                          doCalib=2, gainUse=0, doBand=1, BPVer=BPVer, killAll=True, \
                          IClip=parms["IClip"], XClip=XClip, timeAvg=1.0/60.0, \
                          nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise  RuntimeError("Error in AutoFlag")

# Redo the calibration using new flagging?
if parms["doNDCal2"]==None:
    parms["doNDCal2"] = parms["doNDCal"]
if parms["doBPCal2"]==None:
    parms["doBPCal2"] = parms["doBPCal"]
if parms["doDelayCal2"]==None:
    parms["doDelayCal2"] = parms["doDelayCal2"]
if parms["doAmpPhaseCal2"]==None:
    parms["doAmpPhaseCal2"] = parms["doAmpPhaseCal"]
if parms["doAutoFlag2"]==None:
    parms["doAutoFlagCal2"] = parms["doAutoFlag"]
if parms["doRecal"]:
    mess =  "Redo calibration:"
    printMess(mess, logFile)
    MKClearCal(uv, err, doGain=True, doFlag=False, doBP=True, check=check, logfile=logFile)
    OErr.printErrMsg(err, "Error resetting calibration")
    BPVer = 0
    doBand = -1;
    # Run MKXPhase on delaycal data and attach BP table to UV data
    if parms["doNDCal2"] and parms["doPol"] and doNDCal:
        delay_uv = UV.newPAUV("AIPS UV DATA", MKAIPSName(project),delay_class,disk,delay_seq,True,err)
        mess = "XYphase bandpass noise diode calibration"
        printMess(mess, logFile)
        retCode = MKXPhase(delay_uv, uv, err, logfile=logFile, check=check, debug=debug,
                           doCalib=-1, flagVer=0, doBand=-1, refAnt=parms['refAnt'])
        BPVer += 1
        doBand = 1;
        if retCode!=0:
            raise RuntimeError("Error in Xphase calibration")
        # Plot table?
        if parms["doBPPlot"] and not check:
            # PLPlot doesn't like long names
            plotFile = project+"_XYPhaseBP2.ps"
            MKPlotXYBPTab(uv, 1, plotFile, err, \
                          logfile=logFile, check=check, debug=debug)
     
        
    # Delay recalibration
    if parms["doDelayCal2"] and parms["DCals"] and not check:
        plotFile = fileRoot+"_DelayCal2.ps"
        retCode = MKDelayCal(uv, parms["DCals"], err, \
                             BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                             doCalib=-1, flagVer=0, doBand=doBand, BPVer=BPVer, \
                             solInt=parms["delaySolInt"], smoTime=parms["delaySmoo"],  \
                             refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], \
                             doZeroPhs=parms["delayZeroPhs"], \
                             doAvgIF=parms["delayAvgIF"], doAvgPol=parms["delayAvgPol"], \
                             doPlot=parms["doSNPlot"], plotFile=plotFile, \
                             nThreads=nThreads, noScrat=noScrat, \
                             logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error in delay calibration")
        
        # Plot corrected data? - no turn off
        if parms["doSpecPlot"] and parms["plotSource"] and not True:
            plotFile = fileRoot+"_DelaySpec2.ps"
            retCode = MKSpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                                 Stokes=["RR","LL"], doband=doBand, flagVer=2,  \
                                 check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError("Error in Plotting spectrum")

    # Bandpass calibration
    if parms["doBPCal2"] and parms["BPCals"]:
        retCode = MKBPCal(uv, parms["BPCals"], err, doBand=doBand, BPVer=BPVer, newBPVer=0, \
                          noScrat=noScrat, solInt1=parms["bpsolint1"], \
                          solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                          BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                          BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                          doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                          UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=0, doPlot=False, \
                          nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error in Bandpass calibration")
        
    # Amp & phase Recalibrate
    if parms["doAmpPhaseCal2"]:
        plotFile = fileRoot+"_APCal2.ps"
        retCode = MKCalAP (uv, [], parms["ACals"],  parms["GCalList"], err, PCals=parms["PCals"], \
                           doCalib=2, doBand=1, BPVer=0, flagVer=0, \
                           BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                           solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                           doAmpEdit=True, ampSigma=parms["ampSigma"], \
                           ampEditFG=parms["ampEditFG"], avgPol=parms["doPol"], \
                           doPlot=parms["doSNPlot"], plotFile=plotFile, refAnt=parms["refAnt"], \
                           noScrat=noScrat, nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError("Error calibrating")
        
        # Plot corrected data?
        if parms["doSpecPlot"] and parms["plotSource"]:
            plotFile = fileRoot+"_APSpec2.ps"
            retCode = MKSpectrum(uv, parms["BPCal"], parms["plotTime"], maxgap, plotFile, parms["refAnt"], err, \
                                 Stokes=["RR","LL"], doband=1, flagVer=2,   \
                                 check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError("Error in Plotting spectrum")
    
    # More editing
    if parms["doAutoFlag2"]:
        mess =  "Post recalibration editing:"
        printMess(mess, logFile)
        retCode = MKAutoFlag (uv, [], err, flagVer=0, flagTab=2, \
                              doCalib=2, gainUse=0, doBand=1, BPVer=0,  \
                              IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                              doFD=parms["doSecAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                              FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                              FDmaxRes=parms["FDmaxRes"],  FDmaxResBL= parms["FDmaxResBL"], \
                              FDbaseSel=parms["FDbaseSel"], \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise  RuntimeError("Error in AutoFlag")
# end recal

# Calibrate and average data
if parms["doCalAvg"] == 'Splat':
    retCode = MKCalAvg (uv, avgClass, parms["seq"], parms["CalAvgTime"], err, \
                        flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=0, doPol=False, \
                        avgFreq=parms["avgFreq"], chAvg=parms["chAvg"], Stokes=parms["avgStokes"], \
                        BChan=1, EChan=0, doAuto=parms["doAuto"], \
                        BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                        nThreads=nThreads, logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise  RuntimeError("Error in CalAvg")
elif parms["doCalAvg"] == 'BL':
    retCode = MKBLCalAvg (uv, avgClass, parms["seq"], err, \
                          flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=0, doPol=False, \
                          avgFreq=parms["avgFreq"], chAvg=parms["chAvg"], FOV=parms['blFOV'], \
                          maxInt=min(parms["blMaxInt"],parms["solAInt"]), Stokes=parms["avgStokes"], \
                          BChan=1, EChan=0, timeAvg=parms["CalAvgTime"], \
                          BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                          logfile=logFile, check=check, debug=debug)
    if retCode!=0:
        raise  RuntimeError("Error in BLCalAvg")

# Get calibrated/averaged data
if not check:
    try:
        uvc = UV.newPAUV("AIPS UV DATA", MKAIPSName(project), avgClass[0:6], \
                         disk, parms["seq"], True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    except Exception as exception:
        print(exception)
        exit # Can't do much without it

# Phase calibrate poln calibrators?
gainUse=1
if parms["doPhsCal"] and uvc:
    mess =  "Poln. calibrator phase calibration:"
    printMess(mess, logFile)
    # Not sure this helps
    MKClearCal(uvc, err, doGain=True, check=check, logfile=logFile)
    clist = []
    for c in parms['UnPolCal']:
        if c not in clist:
            clist.append(c)
    for c in parms['XYDCal']:
        if c not in clist:
            clist.append(c)
    for c in parms['GCalList']:
        if c not in clist:
            clist.append(c)
    retCode = MKPhsCal (uvc, clist, parms["refAnt"], err,
                       doCalib=1, gainUse=1, flagVer=1, 
                       solnver=1, solInt=30.0/60.0, avgPol=True, avgIF=False, \
                       nThreads=nThreads,  noScrat=noScrat, \
                       #check=check, debug=True, logfile=logFile)
                       check=check, debug=debug, logfile=logFile)
    gainUse = 2
    if retCode!=0:
        raise  RuntimeError("Error in MKPhsCal")


#  Polarization calibration        
if parms["doPolCal"] and uvc:
    mess =  "Instrumental polarization calibration:"
    printMess(mess, logFile)
    # Use unpolarized and gain cals
    MKPolCal(uvc, parms['UnPolCal'], parms['GCalList'], err, \
             doCalib=1, gainUse=gainUse, doBand=-1, BPVer=0, flagVer=1, \
             solInt=parms["PCSolInt"], solnType=parms["PCSolType"], \
             refAnt=parms["PCRefAnt"], \
             ChInc= parms["PCChInc"], ChWid= parms["PCChWid"], \
             check=False, debug=debug, \
             nThreads=nThreads, noScrat=noScrat, logfile=logFile)
    if err.isErr:
        OErr.printErrMsg(err, "Error in instrumental polarization calibration.")

if parms["doXYDelay"] and uvc:
    mess =  "X-Y phase & delay calibration:"
    printMess(mess, logFile)
    retCode = MKXYDelay(uvc, err,\
                        XYCals=parms["XYCals"], 
                        UVRange=parms["xyUVRange"], timerange=parms["xytimerange"], \
                        doCalib=1, gainUse=gainUse, doPol=True, PDVer=1, doBand=-1, BPVer=0, \
                        flagVer=1, fitType=parms["xyFitType"], nChAvg=parms["xyChAvg"], \
                        timeAvg=parms["xyTimeAvg"], solInt=parms["xySolnInt"], \
                        refAnt=parms["xyrefAnt"], minSNR= parms["xyminSNR"], \
                        nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                        check=check, debug=debug)
    if retCode!=0:
        raise RuntimeError("Error in X-Y delay calibration")
    

if parms["doSaveTab"]:
    filename = project+".CalTab.uvtab"
    caltab = MKUVFITSTab (uv, filename, 0, err, logfile=logFile)
    if err.isErr:
        OErr.printErrMsg(err, "Error saving tables")

# Final Stokes I Spectrum
if parms["doSpecPlot"] and uvc:
    plotFile = fileRoot+"_Spec.ps"
    retCode = MKSpectrum(uvc, parms["BPCal"], parms["plotTime"], maxgap, \
                         plotFile, parms["refAnt"], err, flagVer=1, \
                         Stokes=["I"], doband=-1, docalib=-1,      \
                         check=check, debug=debug, logfile=logFile )
    if retCode!=0:
        raise  RuntimeError("Error in Plotting spectrum")

# Plot polarization spectra?
if parms["doPolSpecPlot"] and parms["XYDCal"] and uvc:
    mess =  "Polarized Spectral plot for: "+parms["XYDCal"][0]
    printMess(mess, logFile)
    plotFile = fileRoot+"_PolSpec.ps"
    # Lookup time range
    plotTime = [0.,0.]; plotSrc = parms["XYDCal"]
    for t in meta['targets']:
        if t[1]==plotSrc[0]:
            id = t[2]; break
    for s in meta['sched']:
        if id==s[1]:
            plotTime = [s[2],s[2]+s[3]]
            break
    # Have to mislabel to get POSSM to work.
    retCode = MKSpectrum(uvc, plotSrc, plotTime, maxgap, plotFile, parms["refAnt"], err, \
                         Stokes=["RL","LR"], doband=-1, flagVer=1,  docalib=1, doPol=True, PDVer=1,  \
                         check=check, debug=debug, logfile=logFile )
    if retCode!=0:
        raise  RuntimeError("Error in Plotting spectrum")
    MKAddOutFile(plotFile, 'project', 'Pipeline log file' )

# Save calibrated, averaged uv data
if parms['doSaveUV'] and uvc and not check:
    mess =  "Save UV data:"
    printMess(mess, logFile)
    MKUVFITab(uvc, fileRoot+'.uvtab', 0, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing uvtab data")
#Gzip the data?

# Image
if parms["doImage"] and uvc:
    mess =  "Image targets:"
    printMess(mess, logFile)
    # If targets not specified, image all
    if len(targets)<=0:
        slist = MKAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    MKImageTargets (uvc, err, Sources=slist, seq=parms["seq"], sclass=parms["outIClass"], \
                    doCalib=2, doBand=-1,  flagVer=1, doPol=parms["doPol"], PDVer=1,  \
                    Stokes=parms["Stokes"], FOV=parms["FOV"], Robust=parms["Robust"], \
                    Niter=parms["Niter"],  CleanRad=parms["CleanRad"], minFlux=parms["minFlux"], \
                    maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], \
                    solPInt=parms["solPInt"], solPMode=parms["solPMode"], solPType=parms["solPType"], \
                    maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], \
                    solAInt=parms["solAInt"], solAMode=parms["solAMode"], solAType=parms["solAType"], \
                    avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR = 5.0, refAnt=parms["refAnt"], \
                    do3D=False, BLFact=parms["BLFact"], BLchAvg=parms["BLchAvg"], \
                    doMB = True, norder=1, maxFBW=parms["MBmaxFBW"], \
                    nThreads=nThreads, doGPU=parms["doGPU"], doGPUGrid=parms["doGPUGrid"], \
                    noScrat=noScrat, logfile=logFile, check=check, debug=debug)                       
    if err.isErr:
        OErr.printErrMsg(err, "Error imaging targets")
    # End image

# Save images
if parms["doSaveImg"] and not check:
    mess =  "Save Target FITS images:"
    printMess(mess, logFile)
    if len(targets)<=0:
        slist = MKAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    for target in slist:
        for s in parms["Stokes"]:
            oclass = s+parms["outIClass"][1:]; outname = target
            # Test if image exists
            cno = AIPSDir.PTestCNO(disk, user, outname, oclass, "MA", parms["seq"], err)
            if cno <= 0 :
                mess = outname+"_"+oclass+"_"+disk+"_"+parms["seq"]+" Not Found"
                printMess(mess, logFile)
                continue
            x = Image.newPAImage("out", outname, oclass, disk, parms["seq"], True, err)
            outfile = project+session+"_"+band+"_"+target+"."+oclass+".fits"
            xf = MKImFITS (x, outfile, 0, err, logfile=logFile)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing FITS images")
    # end writing loop

# Get report on sources
if parms["doReport"] and not check:
    # If targets not specified, do all
    if len(targets)<=0:
        slist = MKAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    Report = MKReportTargets(uv, err, Sources=slist, seq=parms["seq"], sclass=parms["outIClass"], \
                             Stokes=parms["Stokes"], logfile=logFile, check=check, debug=debug)
    # Save to pickle jar
    ReportPicklefile = "./"+project+"_"+session+"_"+band+"Report.pickle"   # Where results saved
    SaveObject(Report, ReportPicklefile, True) 
   
# Cleanup
if parms["doCleanup"] and not check:
    mess =  "Cleanup AIPS directory:"
    printMess(mess, logFile)
    if len(targets)<=0:
        slist = MKAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = targets
    mess ="Cleanup AIPS Files" 
    printMess(mess, logFile)
    uv.Zap(err)   # Initial UV data
    uvc.Zap(err)  # Calibrated/averages UV data
    # Target list
    for target in slist:
        for s in parms["Stokes"]:
            oclass = s+parms["outIClass"][1:]
            # Test if image exists
            cno = AIPSDir.PTestCNO(disk, user, target, oclass, "MA", parms["seq"], err)
            if cno <= 0 :
                mess = "Image"+target+" "+oclass+" Not Found"
                printMess(mess, logFile)
                continue
            x = Image.newPAImage("out", target, oclass, disk, parms["seq"], True, err)
            x.Zap(err) # cleanup
    if err.isErr:
        OErr.printErrMsg(err, "Error cleaning up AIPS data")
# end cleanup

# Shutdown
mess = "Finished project "+project+" session "+session+" "+band+" Band"+" AIPS user no. "+str(AIPS.userno)
printMess(mess, logFile)
OErr.printErr(err)
OSystem.Shutdown(ObitSys)

class DataProductError(Exception):
    """ Exception for data product (output file) errors. """
    pass


