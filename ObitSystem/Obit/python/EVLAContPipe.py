#! /usr/bin/env ObitTalk
"""

The EVLA Continuum Pipeline.  The pipeline can be invoked from the command line 
as, ::

    $ EVLAContPipe AipsSetupScript PipelineParamScript

where the required arguments are

* *AipsSetupScript* = an AIPS setup script (an example of this file is stored in 
    ``Obit/share/scripts``)
* *PipelineParamScript* = the EVLA continuum pipeline input parameters script 
    (a template is in ``Obit/share/scripts``)
"""

import sys, pydoc 
from optparse import OptionParser
from ConfigParser import NoSectionError
import OErr, OSystem, UV, AIPS, FITS, OTObit
import ObitTalkUtil
from AIPS import AIPSDisk
from FITS import FITSDisk
from PipeUtil import *
from EVLACal import *

def pipeline( aipsSetup, parmFile):
    """
    EVLA Continuum pipeline.
    
    * *aipsSetup* = AIPS setup file
    * *parmFile* = pipeline input parameters file
    """
    ############################# Initialize OBIT ##########################################
    noScrat     = []    
    exec(open(aipsSetup).read())
    EVLAAddOutFile( aipsSetup, 'project', "Obit's AIPS setup file" )
    
    ############################# Default parameters ##########################################
    
    # Initialize parameters
    parms = EVLAInitContParms()
    
    ############################# Set Project Processing parameters ##################
    print "parmFile",parmFile
    exec(open(parmFile).read())
    EVLAAddOutFile( parmFile, 'project', 'Pipeline input parameters' )

    # frequency/configuration dependent default parameters
    EVLAInitContFQParms(parms)

    # General data parameters
    band      = parms["band"]           # Observing band
    dataClass = ("UVDa"+band)[0:6]      # AIPS class of raw uv data
    project   = parms["project"][0:12]  # Project name (12 char or less, used as AIPS Name)
    session   = parms["session"]        # Project session code

    ################################## Process #####################################
    fileRoot      = parms["project"]+"_"+parms["session"]+"_"+parms["band"] # root of file name
    logFile       = fileRoot +".log"   # Processing log file
    uv            = None
    uvc           = None
    avgClass      = ("UVAv"+band)[0:6]  # Averaged data AIPS class
    outIClass     =  parms["outIClass"] # image AIPS class

    # Load the outputs pickle jar
    EVLAFetchOutFiles()

    # Logging directly to logFile
    OErr.PInit(err, parms["prtLv"], logFile)
    OSystem.PAllowThreads(nThreads)   # Allow threads in Obit/oython
    retCode = 0
    EVLAAddOutFile( logFile, 'project', 'Pipeline log file' )
   
    mess = "Start project "+parms["project"]+" session "+parms["session"]+\
           " "+parms["band"]+" Band"+" AIPS user no. "+str(AIPS.userno)+\
           ", EVLA configuration "+parms["VLACfg"]
    printMess(mess, logFile)
    if debug:
        pydoc.ttypager = pydoc.plainpager # don't page task input displays
        mess = "Using Debug mode "
        printMess(mess, logFile)
    if check:
        mess = "Only checking script"
        printMess(mess, logFile)
    
    # Log parameters
    printMess("Parameter settings", logFile)
    for p in parms:
        mess = "  "+p+": "+str(parms[p])
        printMess(mess, logFile)

    # Save parameters to pickle jar, manifest
    ParmsPicklefile = project+"_"+session+"_"+band+".Parms.pickle"   # Where results saved
    SaveObject(parms, ParmsPicklefile, True)
    EVLAAddOutFile( ParmsPicklefile, 'project', 'Processing parameters used' )

    # Are we going to be doing Hanning?
    if parms["doHann"]:
        loadClass = parms["band"]+"Raw"
    else:
        loadClass = dataClass
    
    # Load Data from Archive directory
    if parms["doLoadArchive"]:
        uv = EVLAUVLoadArch(parms["archRoot"], EVLAAIPSName(project, session), loadClass, disk, parms["seq"], err, \
                            selConfig=parms["selConfig"], doSwPwr=parms["doSwPwr"], \
                            selBand=parms["selBand"], selChan=parms["selChan"], \
                            selNIF=parms["selNIF"], calInt=parms["calInt"], \
                            logfile=logFile, Compress=parms["Compress"], check=check, debug=debug)
        if uv==None and not check:
            raise RuntimeError,"Cannot load "+parms["DataRoot"]
    
    # Hanning
    if parms["doHann"]:
        # Set uv if not done
        if uv==None and not check:
            uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project, session), loadClass[0:6], disk, parms["seq"], True, err)
            if err.isErr:
                OErr.printErrMsg(err, "Error creating AIPS data")
    
        uv = EVLAHann(uv, EVLAAIPSName(project, session), dataClass, disk, parms["seq"], err, \
                      doDescm=parms["doDescm"], logfile=logFile, check=check, debug=debug)
        if uv==None and not check:
            raise RuntimeError,"Cannot Hann data "
    
    # Set uv is not done
    if uv==None and not check:
        uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project, session), dataClass[0:6], \
                        disk, parms["seq"], True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating AIPS data")
    
    # Clear any old calibration/editing 
    if parms["doClearTab"]:
        mess =  "Clear previous calibration"
        printMess(mess, logFile)
        EVLAClearCal(uv, err, doGain=parms["doClearGain"], doFlag=parms["doClearFlag"], doBP=parms["doClearBP"], check=check)
        OErr.printErrMsg(err, "Error resetting calibration")
    
    # Copy FG 1 to FG 2
    if parms["doCopyFG"]:
        mess =  "Copy FG 1 to FG 2"
        printMess(mess, logFile)
        retCode = EVLACopyFG (uv, err, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Copying FG table"

    # Drop end channels of spectra?  Only if new FG 2
    if parms["doCopyFG"] and (parms["BChDrop"]>0) or (parms["EChDrop"]>0):
        # Channels based on original number, reduced if Hanning
        nchan = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocf"]]
        fact = parms["selChan"]/nchan   # Hanning reduction factor
        BChDrop = parms["BChDrop"]/fact
        EChDrop = parms["EChDrop"]/fact
        mess =  "Trim %d channels from start and %d from end of each spectrum"%(BChDrop,EChDrop)
        printMess(mess, logFile)
        retCode = EVLADropChan (uv, BChDrop, EChDrop, err, flagVer=parms["editFG"], \
                                logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Copying FG table"
   
    # Special editing
    if parms["doEditList"] and not check:
        mess =  "Special editing"
        printMess(mess, logFile)
        for edt in parms["editList"]:
            UV.PFlag(uv,err,timeRange=[dhms2day(edt["timer"][0]),dhms2day(edt["timer"][1])], \
                         flagVer=parms["editFG"], Ants=edt["Ant"], Chans=edt["Chans"], IFs=edt["IFs"], \
                         Stokes=edt["Stokes"], Reason=edt["Reason"])
            OErr.printErrMsg(err, "Error Flagging")
    
    # Quack to remove data from start and end of each scan
    if parms["doQuack"]:
        retCode = EVLAQuack (uv, err, begDrop=parms["quackBegDrop"], endDrop=parms["quackEndDrop"], \
                             Reason=parms["quackReason"], \
                             logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Quacking data"
    
    # Flag antennas shadowed by others?
    if parms["doShad"]:
        retCode = EVLAShadow (uv, err, shadBl=parms["shadBl"], \
                              logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Shadow flagging data"
    
    # Median window time editing, for RFI impulsive in time
    if parms["doMedn"]:
        mess =  "Median window time editing, for RFI impulsive in time:"
        printMess(mess, logFile)
        retCode = EVLAMedianFlag (uv, "    ", err, noScrat=noScrat, nThreads=nThreads, \
                                  avgTime=parms["avgTime"], avgFreq=parms["avgFreq"],  chAvg= parms["chAvg"], \
                                  timeWind=parms["timeWind"], flagVer=2,flagSig=parms["mednSigma"], \
                                  logfile=logFile, check=check, debug=False)
        if retCode!=0:
            raise RuntimeError,"Error in MednFlag"
    
    # Median window frequency editing, for RFI impulsive in frequency
    if parms["doFD1"]:
        mess =  "Median window frequency editing, for RFI impulsive in frequency:"
        printMess(mess, logFile)
        retCode = EVLAAutoFlag (uv, "    ", err,  flagVer=2, doCalib=-1, doBand=-1,   \
                                timeAvg=parms["FD1TimeAvg"], \
                                doFD=True, FDmaxAmp=1.0e20, FDmaxV=1.0e20, FDwidMW=parms["FD1widMW"],  \
                                FDmaxRMS=[1.0e20,0.1], FDmaxRes=parms["FD1maxRes"],  \
                                FDmaxResBL= parms["FD1maxRes"],  \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in AutoFlag"
    
    
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
        retCode = EVLAAutoFlag (uv, clist, err,  flagVer=2, doCalib=-1, doBand=-1,   \
                                    RMSAvg=parms["RMSAvg"], timeAvg=parms["RMSTimeAvg"], \
                                    nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in AutoFlag"
    
    
    # Parallactic angle correction?
    if parms["doPACor"]:
        retCode = EVLAPACor(uv, err, \
                                logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Parallactic angle correction"
    
    # Need to find a reference antenna?  See if we have saved it?
    if (parms["refAnt"]<=0):
        refAnt = FetchObject(project+"_"+session+"_"+band+".refAnt.pickle")
        if refAnt:
            parms["refAnt"] = refAnt
    # Use bandpass calibrator and center half of each spectrum
    if parms["refAnt"]<=0:
        mess = "Find best reference antenna: run Calib on BP Cal(s) "
        printMess(mess, logFile)
        parms["refAnt"] = EVLAGetRefAnt(uv, parms["BPCals"], err, flagVer=2, \
                                        solInt=parms["bpsolint1"], nThreads=nThreads, \
                                        logfile=logFile, check=check, debug=debug)
        if err.isErr:
                raise  RuntimeError,"Error finding reference antenna"
        if parms["refAnts"][0]<=0:
            parms["refAnts"][0] = parms["refAnt"]
        mess = "Picked reference antenna "+str(parms["refAnt"])
        printMess(mess, logFile)
        # Save it
        ParmsPicklefile = project+"_"+session+"_"+band+".Parms.pickle"   # Where results saved
        SaveObject(parms, ParmsPicklefile, True)
        refAntPicklefile = project+"_"+session+"_"+band+".refAnt.pickle"   # Where results saved
        SaveObject(parms["refAnt"], refAntPicklefile, True)


    # Plot Raw, edited data?
    if parms["doRawSpecPlot"] and parms["plotSource"]:
        mess =  "Raw Spectral plot for: "+parms["plotSource"]
        printMess(mess, logFile)
        plotFile = "./"+fileRoot+"RawSpec.ps"
        retCode = EVLASpectrum(uv, parms["plotSource"], parms["plotTime"], plotFile, parms["refAnt"], err, \
                               Stokes=["RR","LL"], doband=-1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"
        EVLAAddOutFile( plotFile, 'project', 'Pipeline log file' )

    # delay calibration
    if parms["doDelayCal"] and parms["DCals"] and not check:
        plotFile = "./"+fileRoot+"DelayCal.ps"
        retCode = EVLADelayCal(uv, parms["DCals"], err,  \
                               BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                               doCalib=2, flagVer=2, doBand=-1, \
                               solInt=parms["solInt"], smoTime=1.0/60.0,  \
                               refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], 
                               doZeroPhs=parms["delayZeroPhs"], \
                               doPlot=parms["doSNPlot"], plotFile=plotFile, \
                               nThreads=nThreads, noScrat=noScrat, \
                               logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in delay calibration"
        
        # Plot corrected data?
        if parms["doSpecPlot"] and parms["plotSource"]:
            plotFile = "./"+fileRoot+"DelaySpec.ps"
            retCode = EVLASpectrum(uv, parms["plotSource"], parms["plotTime"], \
                                   plotFile, parms["refAnt"], err, \
                                   Stokes=["RR","LL"], doband=-1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"

    # Bandpass calibration
    if parms["doBPCal"] and parms["BPCals"]:
        retCode = EVLABPCal(uv, parms["BPCals"], err, noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                            BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                            UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=2, doPlot=False, \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Bandpass calibration"
        
        # Plot corrected data?
        if parms["doSpecPlot"] and  parms["plotSource"]:
            plotFile = "./"+fileRoot+"BPSpec.ps"
            retCode = EVLASpectrum(uv, parms["plotSource"], parms["plotTime"], plotFile, \
                                   parms["refAnt"], err, Stokes=["RR","LL"], doband=1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"

    # Amp & phase Calibrate
    if parms["doAmpPhaseCal"]:
        plotFile = "./"+fileRoot+"APCal.ps"
        retCode = EVLACalAP (uv, [], parms["ACals"], err, PCals=parms["PCals"], 
                             doCalib=2, doBand=1, BPVer=1, flagVer=2, \
                             BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                             solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                             doAmpEdit=parms["doAmpEdit"], ampSigma=parms["ampSigma"], \
                             ampEditFG=parms["ampEditFG"], \
                             doPlot=parms["doSNPlot"], plotFile=plotFile,  refAnt=parms["refAnt"], \
                             nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error calibrating"
    
    # More editing
    if parms["doAutoFlag"]:
        mess =  "Post calibration editing:"
        printMess(mess, logFile)
        # if going to redo then only calibrators
        if parms["doRecal"]:
            # Only calibrators
            clist = []
            for DCal in DCals:
                if DCal["Source"] not in clist:
                    clist.append(DCal["Source"])
            for PCal in PCals:
                if PCal["Source"] not in clist:
                    clist.append(DCal["Source"])
            for ACal in ACals:
                if ACal["Source"] not in clist:
                    clist.append(ACal["Source"])
        else:
            clist = []
    
        retCode = EVLAAutoFlag (uv, clist, err, flagVer=2, \
                                doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                                IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                                doFD=parms["doAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                                FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                                FDmaxRes=parms["FDmaxRes"],  FDmaxResBL=parms["FDmaxResBL"], \
                                FDbaseSel=parms["FDbaseSel"], \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in AutoFlag"
    
    # Redo the calibration using new flagging?
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
        EVLAClearCal(uv, err, doGain=True, doFlag=False, doBP=True, check=check, logfile=logFile)
        OErr.printErrMsg(err, "Error resetting calibration")
        # Parallactic angle correction?
        if parms["doPACor"]:
            retCode = EVLAPACor(uv, err, \
                                logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError,"Error in Parallactic angle correction"
    
        # Delay recalibration
        if parms["doDelayCal2"] and parms["DCals"] and not check:
            plotFile = "./"+fileRoot+"DelayCal2.ps"
            retCode = EVLADelayCal(uv, parms["DCals"], err, \
                                   BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                                   doCalib=2, flagVer=2, doBand=-1, \
                                   solInt=parms["solInt"], smoTime=1.0/60.0,  \
                                   refAnts=[parms["refAnt"]], doTwo=parms["doTwo"], \
                                   doZeroPhs=parms["delayZeroPhs"], \
                                   doPlot=parms["doSNPlot"], plotFile=plotFile, \
                                   nThreads=nThreads, noScrat=noScrat, \
                                   logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError,"Error in delay calibration"
                
            # Plot corrected data?
            if parms["doSpecPlot"] and parms["plotSource"]:
                plotFile = "./"+fileRoot+"DelaySpec2.ps"
                retCode = EVLASpectrum(uv, parms["plotSource"], parms["plotTime"], plotFile, parms["refAnt"], err, \
                                       Stokes=["RR","LL"], doband=-1,          \
                                       check=check, debug=debug, logfile=logFile )
                if retCode!=0:
                    raise  RuntimeError,"Error in Plotting spectrum"
    
    # Bandpass calibration
    if parms["doBPCal2"] and parms["BPCals"]:
        retCode = EVLABPCal(uv, parms["BPCals"], err, noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                            BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                            UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=2, doPlot=False, \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Bandpass calibration"
        
        # Plot corrected data?
        if parms["doSpecPlot"] and parms["plotSource"]:
            plotFile = "./"+fileRoot+"BPSpec2.ps"
            retCode = EVLASpectrum(uv, parms["plotSource"], parms["plotTime"], plotFile, parms["refAnt"], err, \
                                   Stokes=["RR","LL"], doband=1,          \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"
    
        # Amp & phase Recalibrate
        if parms["doAmpPhaseCal2"]:
            plotFile = "./"+fileRoot+"APCal2.ps"
            retCode = EVLACalAP (uv, [], parms["ACals"], err, PCals=parms["PCals"], \
                                 doCalib=2, doBand=1, BPVer=1, flagVer=2, \
                                 BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                                 solInt=parms["solInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
                                 doAmpEdit=parms["doAmpEdit"], ampSigma=parms["ampSigma"], \
                                 ampEditFG=parms["ampEditFG"], \
                                 doPlot=parms["doSNPlot"], plotFile=plotFile, refAnt=parms["refAnt"], \
                                 noScrat=noScrat, nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise RuntimeError,"Error calibrating"
    
        # More editing
        if parms["doAutoFlag2"]:
            mess =  "Post recalibration editing:"
            printMess(mess, logFile)
            retCode = EVLAAutoFlag (uv, [], err, flagVer=2, \
                                    doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                                    IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                                    doFD=parms["doAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                                    FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                                    FDmaxRes=parms["FDmaxRes"],  FDmaxResBL= parms["FDmaxResBL"], \
                                    FDbaseSel=parms["FDbaseSel"], \
                                    nThreads=nThreads, logfile=logFile, check=check, debug=debug)
            if retCode!=0:
                raise  RuntimeError,"Error in AutoFlag"
            
    # end recal
    
    
    
    # Calibrate and average data
    if parms["doCalAvg"]:
        retCode = EVLACalAvg (uv, avgClass, parms["seq"], parms["CalAvgTime"], err, \
                              flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1, doPol=False, \
                              avgFreq=parms["avgFreq"], chAvg=parms["chAvg"], \
                              BChan=parms["CABChan"], EChan=parms["CAEChan"], \
                              BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in CalAvg"
       
    # Get calibrated/averaged data
    if not check:
        uv = UV.newPAUV("AIPS UV DATA", EVLAAIPSName(project, session), avgClass[0:6], \
                        disk, parms["seq"], True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    
    # XClip
    if parms["XClip"] and parms["XClip"]>0.0:
        mess =  "Cross Pol clipping:"
        printMess(mess, logFile)
        retCode = EVLAAutoFlag (uv, [], err, flagVer=-1, flagTab=1, \
                                doCalib=2, gainUse=0, doBand=-1, maxBad=1.0,  \
                                XClip=parms["XClip"], timeAvg=1./60., \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise  RuntimeError,"Error in AutoFlag"
    
    # R-L  delay calibration cal if needed,
    if parms["doRLDelay"] and parms["RLDCal"][0][0]!=None:
        if parms["rlrefAnt"]<=0:
            parms["rlrefAnt"] =  parms["refAnt"]
        # parms["rlDoBand"] if before average, BPVer=parms["rlBPVer"], 
        retCode = EVLARLDelay(uv, err,\
                              RLDCal=parms["RLDCal"], BChan=parms["rlBChan"], \
                              EChan=parms["rlEChan"], UVRange=parms["rlUVRange"], \
                              soucode=parms["rlCalCode"], doCalib=parms["rlDoCal"], gainUse=parms["rlgainUse"], \
                              timerange=parms["rltimerange"], \
                              # NOT HERE doBand=parms["rlDoBand"], BPVer=parms["rlBPVer"],  \
                              flagVer=parms["rlflagVer"], \
                              refAnt=parms["rlrefAnt"], doPol=False,  \
                              nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                              check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in R-L delay calibration"
    
    # Polarization calibration
    if parms["doPolCal"]:
        if parms["PCRefAnt"]<=0:
            parms["PCRefAnt"] =  parms["refAnt"]
        retCode = EVLAPolCal(uv, parms["PCInsCals"], err, \
                             doCalib=2, gainUse=0, doBand=-1, flagVer=0, \
                             fixPoln=parms["PCFixPoln"], pmodel=parms["PCpmodel"], avgIF=parms["PCAvgIF"], \
                             solInt=parms["PCSolInt"], refAnt=parms["PCRefAnt"], solType=parms["PCSolType"], \
                             ChInc=parms["PCChInc"], ChWid=parms["PCChWid"], \
                             nThreads=nThreads, check=check, debug=debug, noScrat=noScrat, logfile=logFile)
        if retCode!=0 and (not check):
           raise  RuntimeError,"Error in polarization calibration: "+str(retCode)
        # end poln cal.
    
    
    # R-L phase calibration cal., creates new BP table
    if parms["doRLCal"] and parms["RLDCal"][0][0]!=None:
        plotFile = "./"+fileRoot+"RLSpec2.ps"
        if parms["rlrefAnt"]<=0:
            parms["rlrefAnt"] =  parms["refAnt"]
        retCode = EVLARLCal(uv, err,\
                            RLDCal=parms["RLDCal"], BChan=parms["rlBChan"],
                            EChan=parms["rlEChan"], UVRange=parms["rlUVRange"], \
                            ChWid2=parms["rlChWid"], solInt1=parms["rlsolint1"], solInt2=parms["rlsolint2"], \
                            RLPCal=parms["RLPCal"], RLPhase=parms["RLPhase"], \
                            RM=parms["RLRM"], CleanRad=parms["rlCleanRad"], \
                            calcode=parms["rlCalCode"], doCalib=parms["rlDoCal"], gainUse=parms["rlgainUse"], \
                            timerange=parms["rltimerange"], FOV=parms["rlFOV"], \
                            doBand=-1, BPVer=1, flagVer=parms["rlflagVer"], \
                            refAnt=parms["rlrefAnt"], doPol=parms["doPol"], PDVer=parms["PDVer"],  \
                            doPlot=parms["doSpecPlot"], plotFile=plotFile, \
                            nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                            check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in RL phase spectrum calibration"
    
    # VClip
    if parms["VClip"] and parms["XClip"]>0.0:
        mess =  "VPol clipping:"
        printMess(mess, logFile)
        retCode = EVLAAutoFlag (uv, [], err, flagVer=-1, flagTab=1, \
                                doCalib=2, gainUse=0, doBand=-1,  \
                                VClip=parms["VClip"], timeAvg=parms["timeAvg"], \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise  RuntimeError,"Error in AutoFlag VClip"
    
    # Plot corrected data?
    if parms["doSpecPlot"] and parms["plotSource"]:
        plotFile = "./"+fileRoot+"Spec.ps"
        retCode = EVLASpectrum(uv, parms["plotSource"], parms["plotTime"], \
                               plotFile, parms["refAnt"], err, \
                               Stokes=["RR","LL"], doband=-1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"
    
    # Image targets
    if parms["doImage"]:
        # If targets not specified, image all
        if len(parms["targets"])<=0:
            slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
        else:
            slist = parms["targets"]
        EVLAImageTargets (uv, err, Sources=slist, seq=parms["seq"], sclass=outIClass, \
                          doCalib=2, doBand=1,  flagVer=2, doPol=parms["doPol"], PDVer=parms["PDVer"],  \
                          Stokes=parms["Stokes"], FOV=parms["FOV"], Robust=parms["Robust"], Niter=parms["Niter"], \
                          CleanRad=parms["CleanRad"], minFlux=parms["minFlux"], \
                          maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], \
                          solPInt=parms["solPInt"], solPMode=parms["solPMode"], solPType=parms["solPType"], \
                          maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], \
                          solAInt=parms["solAInt"], solAMode=parms["solAMode"], solAType=parms["solAType"], \
                          avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR = 4.0, refAnt=parms["refAnt"], \
                          do3D=parms["do3D"], BLFact=parms["BLFact"], BLchAvg=parms["BLchAvg"], \
                          doMB=parms["doMB"], norder=parms["MBnorder"], maxFBW=parms["MBmaxFBW"], \
                          nTaper=parms["nTaper"], Tapers=parms["Tapers"], \
                          nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
        # End image
    
    # Get report on sources
    if parms["doReport"]:
        # If targets not specified, do all
        if len(parms["targets"])<=0:
            slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
        else:
            slist = parms["targets"]
        Report = EVLAReportTargets(uv, err, Sources=slist, seq=parms["seq"], sclass=outIClass, \
                                       Stokes=parms["Stokes"], logfile=logFile, check=check, debug=debug)
        # Save to pickle jar
        ReportPicklefile = "./"+fileRoot+"Report.pickle"   # Where results saved
        SaveObject(Report, ReportPicklefile, True) 
       
    # Write results, cleanup    
    # Save cal/average UV data? 
    if parms["doSaveUV"] and (not check):
        Aname = EVLAAIPSName(project, session)
        cno = AIPSDir.PTestCNO(disk, user, Aname, avgClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS CAL UV DATA", Aname, avgClass, disk, parms["seq"], True, err)
            filename = parms["project"]+parms["session"]+parms["band"]+"Cal.uvtab"
            fuv = EVLAUVFITS (uvt, filename, 0, err, compress=parms["Compress"], logfile=logFile)
            EVLAAddOutFile( filename, 'project', "Calibrated Averaged UV data" )
            # Save list of output files
            EVLASaveOutFiles()
            del uvt
    # Save raw UV data tables?
    if parms["doSaveTab"] and (not check):
        Aname = EVLAAIPSName(project, session)
        cno = AIPSDir.PTestCNO(disk, user, Aname, dataClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS RAW UV DATA", Aname, dataClass[0:6], disk, parms["seq"], True, err)
            filename = parms["project"]+parms["session"]+parms["band"]+"CalTab.uvtab"
            fuv = EVLAUVFITSTab (uvt, filename, 0, err, logfile=logFile)
            EVLAAddOutFile( filename, 'project', "Calibrated AIPS tables" )
            del uvt
            # Write History
            filename = project+'_'+session+'_'+band+".History.text"
            OTObit.PrintHistory(uv, file=filename)
            EVLAAddOutFile( filename, 'project', "Processing history of calibrated data" )
            # Save list of output files
            EVLASaveOutFiles()
    # Imaging results
    # If targets not specified, save all
    if len(parms["targets"])<=0:
        slist = EVLAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    else:
        slist = parms["targets"]
    for target in slist:
        if parms["doSaveImg"] and (not check):
            for s in parms["Stokes"]:
                oclass = s+outIClass[1:]
                outname = target
                # Test if image exists
                cno = AIPSDir.PTestCNO(disk, user, outname, oclass, "MA", parms["seq"], err)
                if cno <= 0 :
                    continue
                x = Image.newPAImage("out", outname, oclass, disk, parms["seq"], True, err)
                outfile = "./"+fileRoot+target+"."+oclass+".fits"
                xf = EVLAImFITS (x, outfile, 0, err, logfile=logFile)
                EVLAAddOutFile( outfile, target, 'Image of '+ target)
                # Statistics
                zz=imstat(x, err, logfile=logFile)
    # end writing loop
    
    # Save list of output files
    EVLASaveOutFiles()
    OErr.printErrMsg(err, "Writing output")
    
    # Contour plots
    if parms["doKntrPlots"]:
        mess = "INFO --> Contour plots (doKntrPlots)"
        printMess(mess, logFile)
        EVLAKntrPlots( err, imName=parms["targets"], project=project,
            session=session, band=band, disk=disk, debug=debug )
        # Save list of output files
        EVLASaveOutFiles()
    elif debug:
        mess = "Not creating contour plots ( doKntrPlots = "+str(parms["doKntrPlots"])+ " )"
        printMess(mess, logFile)
    
    # Source uv plane diagnostic plots
    if parms["doDiagPlots"]:
        mess = "INFO --> Diagnostic plots (doDiagPlots)"
        printMess(mess, logFile)
        # Get the highest number avgClass catalog file
        Aname = EVLAAIPSName( project, session )
        uvc = None
        if not check:
            uvname = project+"_"+session+"_"+band+"_Cal"
            uvc = UV.newPAUV(uvname, Aname, avgClass, disk, parms["seq"], True, err)
        EVLADiagPlots( uvc, err, cleanUp=parms["doCleanup"], \
                           project=project, session=session, band=band, \
                           logfile=logFile, check=check, debug=debug )
        # Save list of output files
        EVLASaveOutFiles()
    elif debug:
        mess = "Not creating diagnostic plots ( doDiagPlots = "+str(parms["doDiagPlots"])+ " )"
        printMess(mess, logFile)
    
    # Save metadata
    srcMetadata = None
    projMetadata = None
    if parms["doMetadata"]:
        mess = "INFO --> Save metadata (doMetadata)"
        printMess(mess, logFile)
        uvc = None
        if not uvc:
            # Get calibrated/averaged data
            Aname = EVLAAIPSName(project, session)
            uvname = project+"_"+session+"_"+band+"_Cal"
            uvc = UV.newPAUV(uvname, Aname, avgClass, disk, parms["seq"], True, err)
            if err.isErr:
                OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    
        # Get source metadata; save to pickle file
        srcMetadata = EVLASrcMetadata( uvc, err, Sources=parms["targets"], seq=parms["seq"], \
                                       sclass=outIClass, Stokes=parms["Stokes"],\
                                       logfile=logFile, check=check, debug=debug )
        picklefile = "./"+fileRoot+".SrcReport.pickle" 
        SaveObject( srcMetadata, picklefile, True ) 
        EVLAAddOutFile( picklefile, 'project', 'All source metadata' )
    
        # Get project metadata; save to pickle file
        projMetadata = EVLAProjMetadata( uvc, AIPS_VERSION, err, \
            PCals=parms["PCals"], ACals=parms["ACals"], \
            BPCals=parms["BPCals"], DCals=parms["DCals"], \
            project = project, session = session, band = band, \
            dataInUVF = parms["archRoot"], archFileID = 66666 )
        picklefile = "./"+fileRoot+".ProjReport.pickle"
        SaveObject(projMetadata, picklefile, True) 
        EVLAAddOutFile( picklefile, 'project', 'Project metadata' )
    else:
        # Fetch from pickle jar
         picklefile = "./"+fileRoot+".SrcReport.pickle"
         srcMetadata = FetchObject(picklefile)
         picklefile = "./"+fileRoot+".ProjReport.pickle"
         projMetadata = FetchObject(picklefile)
   
    # Write report
    if parms["doHTML"]:
        mess = "INFO --> Write HTML report (doHTML)"
        printMess(mess, logFile)
        EVLAHTMLReport( projMetadata, srcMetadata, \
                            outfile="./"+fileRoot+".report.html", \
                            logFile=logFile )
    
    # Write VOTable
    if parms["doVOTable"]:
        mess = "INFO --> Write VOTable (doVOTable)"
        printMess(mess, logFile)
        EVLAAddOutFile( 'VOTable.xml', 'project', 'VOTable report' ) 
        EVLAWriteVOTable( projMetadata, srcMetadata, filename='VOTable.xml' )
    
    # Save list of output files
    EVLASaveOutFiles()
    
    # Cleanup - delete AIPS files
    if parms["doCleanup"] and (not check):
        mess = "INFO --> Clean up (doCleanup)"
        printMess(mess, logFile)
        # Delete target images
        # How many Stokes images
        nstok = len(parms["Stokes"])
        for istok in range(0,nstok):
            oclass = parms["Stokes"][istok:istok+1]+outIClass[1:]
            AllDest(err, disk=disk,Aseq=parms["seq"],Aclass=oclass)
        
        # Delete initial UV data
        Aname = EVLAAIPSName(project, session)
        # Test if data exists
        cno = AIPSDir.PTestCNO(disk, user, Aname, dataClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS RAW UV DATA", Aname, dataClass[0:6], disk, parms["seq"], True, err)
            uvt.Zap(err)
            del uvt
            if err.isErr:
                OErr.printErrMsg(err, "Error deleting raw AIPS data")
        # Zap calibrated/averaged data
        # Test if data exists
        cno = AIPSDir.PTestCNO(disk, user, Aname, avgClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS CAL UV DATA", Aname, avgClass[0:6], disk, parms["seq"], True, err)
            uvt.Zap(err)
            del uvt
            if err.isErr:
                OErr.printErrMsg(err, "Error deleting cal/avg AIPS data")
        # Zap UnHanned data if present
        loadClass = parms["band"]+"Raw"
        # Test if image exists
        cno = AIPSDir.PTestCNO(disk, user, Aname, loadClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS CAL UV DATA", Aname, loadClass[0:6], disk, parms["seq"], True, err)
            uvt.Zap(err)
            del uvt
            if err.isErr:
                OErr.printErrMsg(err, "Error deleting cal/avg AIPS data")
        OErr.printErrMsg(err, "Writing output/cleanup")

    # Shutdown
    mess = "Finished project "+parms["project"]+" session "+parms["session"]+ \
    " "+parms["band"]+" Band"+" AIPS user no. "+str(AIPS.userno)
    printMess(mess, logFile)
    OErr.printErr(err)
    OSystem.Shutdown(ObitSys)
    # end pipeline

class DataProductError(Exception):
    """ Exception for data product (output file) errors. """
    pass

if __name__ == '__main__':
    usage = """usage: %prog [options] AIPSSetup PipelineParms
    AIPSSetup = pipeline AIPS setup file
    PipelineParms = pipeline input parameters file"""
    parser = OptionParser( usage=usage )
    (options, args) = parser.parse_args()
    # Unset LD_PRELOAD to avoid ld.so warnings from binary installation
    os.unsetenv('LD_PRELOAD')
    if len(args) != 2:
        parser.print_help()
        sys.exit()
    try:
        pipeline( args[0] , args[1])
    finally:
        pass
