#! /usr/bin/env ObitTalk
"""

The ALMA Pipeline.  The pipeline can be invoked from the command line 
as, ::

    ObitTalk ALMAPipe.py AipsSetupScript PipelineParamScript

where the required arguments are

* *AipsSetupScript* = an AIPS setup script (an example of this file is stored in 
    ``Obit/share/scripts``)
* *PipelineParamScript* = the ALMA continuum pipeline input parameters script 
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
from ALMACal import *

def pipeline( aipsSetup, parmFile):
    """
    ALMA Continuum pipeline.
    
    * *aipsSetup* = AIPS setup file
    * *parmFile* = pipeline input parameters file
    """
    ############################# Initialize OBIT ##########################################
    noScrat     = []    
    exec(open(aipsSetup).read())
    ALMAAddOutFile( aipsSetup, 'project', "Obit's AIPS setup file" )
    
    ############################# Default parameters ##########################################
    
    # Initialize parameters
    parms = ALMAInitContParms()
    
    ############################# Set Project Processing parameters ##################
    print "parmFile",parmFile
    exec(open(parmFile).read())
    ALMAAddOutFile( parmFile, 'project', 'Pipeline input parameters' )

    # frequency/configuration dependent default parameters
    ALMAInitContFQParms(parms)

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
    outCClass     =  parms["outCClass"] # image AIPS class

    # Load the outputs pickle jar
    ALMAFetchOutFiles()

    # Logging directly to logFile
    OErr.PInit(err, parms["prtLv"], logFile)
    OSystem.PAllowThreads(nThreads)   # Allow threads in Obit/oython
    retCode = 0
    ALMAAddOutFile( logFile, 'project', 'Pipeline log file' )
   
    mess = "Start project "+parms["project"]+" session "+parms["session"]+\
           " "+parms["band"]+" Band"+" AIPS user no. "+str(AIPS.userno)+\
           ", ALMA Max. baseline "+str(parms["ALMAMaxBl"])
    printMess(mess, logFile)
    if debug:
        pydoc.ttypager = pydoc.plainpager # don't page task input displays
        mess = "Using Debug mode "
        printMess(mess, logFile)
    if check:
        mess = "Only checking script"
        printMess(mess, logFile)
    
    # Log parameters
    if parms['doLogParms']:
        printMess("Parameter settings", logFile)
        for p in parms:
            mess = "  "+p+": "+str(parms[p])
            printMess(mess, logFile)

    # Save parameters to pickle jar, manifest
    ParmsPicklefile = project+"_"+session+"_"+band+".Parms.pickle"   # Where results saved
    SaveObject(parms, ParmsPicklefile, True)
    ALMAAddOutFile( ParmsPicklefile, 'project', 'Processing parameters used' )

    # Are we going to be doing Hanning?
    if parms["doHann"]:
        loadClass = parms["band"]+"Raw"
    else:
        loadClass = dataClass
    
    # Load Data from Archive directory
    if parms["doLoadArchive"]:
        uv = ALMAUVLoadArch(parms["archRoots"], ALMAAIPSName(project, session), loadClass, disk, parms["seq"], err, \
                            selConfig=parms["selConfig"], selBand=parms["selBand"], selChan=parms["selChan"], \
                            selChBW=parms["selChBW"], selNIF=parms["selNIF"], calInt=parms["calInt"], \
                            logfile=logFile, Compress=parms["Compress"], check=check, debug=debug)
        if uv==None and not check:
            raise RuntimeError,"Cannot load "+parms["archRoots"]
    
    # Hanning
    if parms["doHann"]:
        # Set uv if not done
        if uv==None and not check:
            uv = UV.newPAUV("AIPS UV DATA", ALMAAIPSName(project, session), loadClass[0:6], disk, parms["seq"], True, err)
            if err.isErr:
                OErr.printErrMsg(err, "Error creating AIPS data")
    
        uv = ALMAHann(uv, ALMAAIPSName(project, session), dataClass, disk, parms["seq"], err, \
                      doDescm=parms["doDescm"], logfile=logFile, check=check, debug=debug)
        if uv==None and not check:
            raise RuntimeError,"Cannot Hann data "
    
    # Set uv if not done
    if uv==None and not check:
        uv = UV.newPAUV("AIPS UV DATA", ALMAAIPSName(project, session), dataClass[0:6], \
                        disk, parms["seq"], True, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating AIPS data")
    
    # Clear any old calibration/editing 
    if parms["doClearTab"]:
        mess =  "Clear previous calibration"
        printMess(mess, logFile)
        ALMAClearCal(uv, err, doGain=parms["doClearGain"], doFlag=parms["doClearFlag"], doBP=parms["doClearBP"], check=check)
        OErr.printErrMsg(err, "Error resetting calibration")
    
    # Copy FG 1 to FG 2
    if parms["doCopyFG"]:
        mess =  "Copy FG 1 to FG 2"
        printMess(mess, logFile)
        retCode = ALMACopyFG (uv, err, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Copying FG table"

    # Drop end channels of spectra?  Only if new FG 2
    if parms["doCopyFG"] and (parms["BChDrop"]>0) or (parms["EChDrop"]>0):
        # Channels based on original number, reduced if Hanning
        nchan = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocf"]]
        fact = max (1,parms["selChan"]/nchan)   # Hanning reduction factor
        BChDrop = parms["BChDrop"]/fact
        EChDrop = parms["EChDrop"]/fact
        mess =  "Trim %d channels from start and %d from end of each spectrum"%(BChDrop,EChDrop)
        printMess(mess, logFile)
        retCode = ALMADropChan (uv, BChDrop, EChDrop, err, flagVer=parms["editFG"], \
                                logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Copying FG table"

    # Get list of source in data
    AllSource = ALMAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
    # Edit source lists to remove sources not present
    for c in parms["PCInsCals"]:
        if c not in AllSource:
             parms["PCInsCals"].remove(c)
    for t in parms["targets"]:
        if t not in AllSource:
            parms["targets"].remove(t)
    if parms["XYGainSource"] not in AllSource:
        parms["XYGainSource"] = None
    if parms["XYDelaySource"] not in AllSource:
        parms["XYDelaySource"] = None
    for c in parms["DCals"]:
        if c['Source'] not in AllSource:
            c['Source'] = None
    for c in parms["BPCals"]:
        if c['Source'] not in AllSource:
            c['Source'] = None
    for c in parms["PCals"]:
        if c['Source'] not in AllSource:
            c['Source'] = None
    for c in parms["ACals"]:
        if c['Source'] not in AllSource:
            c['Source'] = None

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
        retCode = ALMAQuack (uv, err, begDrop=parms["quackBegDrop"], endDrop=parms["quackEndDrop"], \
                             Reason=parms["quackReason"], \
                             logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Quacking data"
    
    # Flag antennas shadowed by others?
    if parms["doShad"]:
        retCode = ALMAShadow (uv, err, shadBl=parms["shadBl"], \
                              logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error Shadow flagging data"
    
    # Apply online calibration
    if parms["doOnlineCal"]:
        retCode = ALMAOnlineCal (uv,err, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error applying online calibration"
    
    # Median window time editing, for RFI impulsive in time
    if parms["doMedn"]:
        mess =  "Median window time editing, for RFI impulsive in time:"
        printMess(mess, logFile)
        retCode = ALMAMedianFlag (uv, "    ", err, noScrat=noScrat, nThreads=nThreads, \
                                  avgTime=parms["avgTime"], avgFreq=parms["avgFreq"],  chAvg= parms["chAvg"], \
                                  timeWind=parms["timeWind"], flagVer=2,flagSig=parms["mednSigma"], \
                                  logfile=logFile, check=check, debug=False)
        if retCode!=0:
            raise RuntimeError,"Error in MednFlag"
    
    # Median window frequency editing, for RFI impulsive in frequency
    if parms["doFD1"]:
        mess =  "Median window frequency editing, for RFI impulsive in frequency:"
        printMess(mess, logFile)
        retCode = ALMAAutoFlag (uv, "    ", err,  flagVer=2, doCalib=-1, doBand=-1,   \
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
        clist =  ALMACombineCals(parms["ACals"], parms["PCals"],  parms["DCals"])   # Calibrator list
        retCode = ALMAAutoFlag (uv, clist, err,  flagVer=2, doCalib=-1, doBand=-1,   \
                                    RMSAvg=parms["RMSAvg"], timeAvg=parms["RMSTimeAvg"], \
                                    nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in AutoFlag"

    # Need to find a reference antenna?  See if we have saved it?
    if (parms["refAnt"]<=0):
        refAnt = FetchObject(project+"_"+session+"_"+band+".refAnt.pickle")
        if refAnt:
            parms["refAnt"] = refAnt
    # Use bandpass calibrator and center half of each spectrum
    if parms["refAnt"]<=0:
        mess = "Find best reference antenna: run Calib on BP Cal(s) "
        printMess(mess, logFile)
        parms["refAnt"] = ALMAGetRefAnt(uv, parms["BPCals"], err, flagVer=2, \
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
        retCode = ALMASpectrum(uv, parms["plotSource"], parms["plotTime"], plotFile, parms["refAnt"], err, \
                               flagVer=2, Stokes=["XX","YY"], doband=-1,          \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"
        ALMAAddOutFile( plotFile, 'project', 'Pipeline log file' )

    # delay calibration
    if parms["doDelayCal"] and parms["DCals"] and not check:
        plotFile = "./"+fileRoot+"DelayCal.ps"
        retCode = ALMADelayCal(uv, parms["DCals"], err,  \
                               BChan=parms["delayBChan"], EChan=parms["delayEChan"], \
                               doCalib=2, flagVer=2, doBand=-1, \
                               solInt=parms["delaySolInt"], smoTime=parms["delaySmoo"],  \
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
            retCode = ALMASpectrum(uv, parms["plotSource"], parms["plotTime"], \
                                   plotFile, parms["refAnt"], err, \
                                   flagVer=2, Stokes=["XX","YY"], doband=-1,  \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"

    # Bandpass calibration
    if parms["doBPCal"] and parms["BPCals"]:
        plotFile = "./"+fileRoot+"BPCal.ps"
        retCode = ALMABPCal(uv, parms["BPCals"], err, noScrat=noScrat, solInt1=parms["bpsolint1"], \
                            solInt2=parms["bpsolint2"], solMode=parms["bpsolMode"], \
                            BChan1=parms["bpBChan1"], EChan1=parms["bpEChan1"], \
                            BChan2=parms["bpBChan2"], EChan2=parms["bpEChan2"], ChWid2=parms["bpChWid2"], \
                            doCenter1=parms["bpDoCenter1"], refAnt=parms["refAnt"], \
                            UVRange=parms["bpUVRange"], doCalib=2, gainUse=0, flagVer=2, doPlot=False, \
                            doAmpEdit=parms["doAmpEdit"], ampSigma=parms["ampSigma"], \
                            doBPPlot=parms["doBPPlot"], plotBPFile=plotFile, \
                            nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in Bandpass calibration"
        
        # Plot corrected data?
        if parms["doSpecPlot"] and  parms["plotSource"]:
            plotFile = "./"+fileRoot+"BPSpec.ps"
            retCode = ALMASpectrum(uv, parms["plotSource"], parms["plotTime"], plotFile, \
                                   parms["refAnt"], err, flagVer=2, Stokes=["XX","YY"], doband=1, \
                                   check=check, debug=debug, logfile=logFile )
            if retCode!=0:
                raise  RuntimeError,"Error in Plotting spectrum"

    # set X/Y gains and initial calibration
    if parms["doXYFixGain"] and parms["XYGainSource"] and not check:
        mess =  "Fix X/Y gain ratios"
        printMess(mess, logFile)
        retCode = ALMAXYGain(uv, err, \
                             XYCal=parms["XYGainSource"],timerange=parms["XYGainTime"],  \
                             doCalib=2, gainUse=0, doBand=1, flagVer=2, refAnt=parms["refAnt"],  \
                             nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                             check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in X-Y gain fix"
    
    # Self calibrate calibrators
    if parms["doImgCal"] and not check:
        mess =  "SelfCalibrate/Image calibrators"
        printMess(mess, logFile)
        src =  ALMACombineCals(parms["ACals"], parms["PCals"], parms["APCals"])   # Source list
        ALMAImageCals(uv, err, Sources=src, seq=parms["seq"], \
            sclass=parms["outCClass"], doCalib=2, flagVer=2, doBand=1, FOV=parms["CalFOV"], \
            maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], solPInt=parms["solPInt"], \
            maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"],\
            solAInt=parms["solAInt"], avgPol=parms["avgPol"],\
            avgIF=parms["avgIF"], minSNR=parms["minSNR"],\
            refAnt=parms["refAnt"], nThreads=nThreads, noScrat=noScrat,\
            logfile=logFile, check=check, debug=debug)

    # Self calibrated models now available
    ALMAImageModel(parms["ACals"],  parms["outCClass"], disk, parms["seq"], err)
    ALMAImageModel(parms["PCals"],  parms["outCClass"], disk, parms["seq"], err)
    ALMAImageModel(parms["APCals"], parms["outCClass"], disk, parms["seq"], err)
    
    # Phase Calibrate
    if parms["doPhaseCal"]:
        plotFile = "./"+fileRoot+"PhaseCal.ps"
        retCode = ALMAPhaseCal (uv, parms["PCals"], err,  ACals=parms["ACals"], \
                                doCalib=2, doBand=1, BPVer=1, flagVer=2, refAnt=parms["refAnt"], \
                                BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                                solInt=parms["solPInt"],  avgIF = parms['avgIF'], \
                                ampScalar=parms["ampScalar"], doPlot=parms["doSNPlot"], plotFile=plotFile,  \
                                nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error calibrating"
    
    # Amp & phase Calibrate
    if parms["doAmpPhaseCal"]:
        plotFile = "./"+fileRoot+"APCal.ps"
        retCode = ALMACalAP (uv, [], parms["ACals"], err, PCals=parms["APCals"], \
                             doCalib=2, doBand=1, BPVer=1, flagVer=2, \
                             BChan=parms["ampBChan"], EChan=parms["ampEChan"], \
                             solInt=parms["solAInt"], solSmo=parms["solSmo"], ampScalar=parms["ampScalar"], \
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
        clist = []
        retCode = ALMAAutoFlag (uv, clist, err, flagVer=2, \
                                doCalib=2, gainUse=0, doBand=1, BPVer=1,  \
                                IClip=parms["IClip"], minAmp=parms["minAmp"], timeAvg=parms["timeAvg"], \
                                doFD=parms["doAFFD"], FDmaxAmp=parms["FDmaxAmp"], FDmaxV=parms["FDmaxV"], \
                                FDwidMW=parms["FDwidMW"], FDmaxRMS=parms["FDmaxRMS"], \
                                FDmaxRes=parms["FDmaxRes"],  FDmaxResBL=parms["FDmaxResBL"], \
                                FDbaseSel=parms["FDbaseSel"], \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in AutoFlag"
    
    
    # Calibrate and average data
    if parms["doCalAvg"]:
        retCode = ALMACalAvg (uv, avgClass, parms["seq"], parms["CalAvgTime"], err, \
                              flagVer=2, doCalib=2, gainUse=0, doBand=1, BPVer=1, doPol=False, \
                              avgFreq=parms["CAavgFreq"], chAvg=parms["CAchAvg"], \
                              BChan=parms["CABChan"], EChan=parms["CAEChan"], \
                              BIF=parms["CABIF"], EIF=parms["CAEIF"], Compress=parms["Compress"], \
                              nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
           raise  RuntimeError,"Error in CalAvg"
       
    # Get calibrated/averaged data
    if not check:
        try:
            uv = UV.newPAUV("AIPS UV DATA", ALMAAIPSName(project, session), avgClass[0:6], \
                                disk, parms["seq"], True, err)
        except Exception, exception:
            mess = "No Averaged/calibrated datafile"
            printMess(mess, logFile)
            return
    
    # XClip
    if parms["XClip"] and parms["XClip"]>0.0:
        mess =  "Cross Pol clipping:"
        printMess(mess, logFile)
        retCode = ALMAAutoFlag (uv, [], err, flagVer=-1, flagTab=1, \
                                doCalib=2, gainUse=0, doBand=-1, maxBad=1.0,  \
                                XClip=parms["XClip"], timeAvg=1./60., \
                                nThreads=nThreads, logfile=logFile, check=check, debug=debug)
        if retCode!=0:
            raise  RuntimeError,"Error in AutoFlag"
    
    # X/Y Delay calibration
    if parms["doXYDelay"] and parms["XYDelaySource"]!=None:
        retCode = ALMAXYDelay(uv, err, \
                              XYDCal=parms["XYDelaySource"], timerange=parms["XYDelayTime"], \
                              BChan=parms["xyBChan"], EChan=parms["xyEChan"], UVRange=parms["xyUVRange"], \
                              doCalib=parms["xyDoCal"], gainUse=parms["xygainUse"], numIFs=parms["xynumIFs"], \
                              flagVer=parms["xyflagVer"], refAnt=parms["refAnt"], doPol=False,  \
                              nThreads=nThreads, noScrat=noScrat, logfile=logFile, \
                              check=check, debug=debug)
        if retCode!=0:
            raise RuntimeError,"Error in X-Y delay calibration"
    
    # Instrumental polarization calibration
    if parms["doPolCal"]:
        if parms["PCRefAnt"]==-10:
            parms["PCRefAnt"] =  parms["refAnt"]
        retCode = ALMAPolCal(uv, parms["PCInsCals"], err, InsCalPoln=parms["PCCalPoln"], \
                             doCalib=2, gainUse=0, doBand=-1, flagVer=0, doFitXY=parms["doFitXY"], \
                             solInt=parms["PCSolInt"], refAnt=parms["PCRefAnt"], solType=parms["PCSolType"], \
                             ChInc=parms["PCChInc"], ChWid=parms["PCChWid"], doFitOri=parms["doFitOri"], \
                             nThreads=nThreads, check=check, debug=debug, noScrat=noScrat, logfile=logFile)
        if retCode!=0 and (not check):
           raise  RuntimeError,"Error in polarization calibration: "+str(retCode)
        # end poln cal.
    
    
    # Plot corrected data?
    if parms["doSpecPlot"] and parms["plotSource"]:
        plotFile = "./"+fileRoot+"Spec.ps"
        retCode = ALMASpectrum(uv, parms["plotSource"], parms["plotTime"], \
                               plotFile, parms["refAnt"], err, \
                               flagVer=1, Stokes=["XX","YY"], doband=-1,    \
                               check=check, debug=debug, logfile=logFile )
        if retCode!=0:
            raise  RuntimeError,"Error in Plotting spectrum"
    
    # Image targets
    if parms["doImage"]:
        # If targets not specified, image all
        if len(parms["targets"])<=0:
            slist = ALMAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
        else:
            slist = parms["targets"]
        ALMAImageTargets (uv, err, Sources=slist, seq=parms["seq"], sclass=outIClass, \
                          doCalib=2, doBand=1,  flagVer=1, doPol=parms["doPol"], PDVer=parms["PDVer"],  \
                          Stokes=parms["Stokes"], FOV=parms["FOV"], Robust=parms["Robust"], Niter=parms["Niter"], \
                          CleanRad=parms["CleanRad"], minFlux=parms["minFlux"], UVRange=parms["UVRange"], \
                          maxPSCLoop=parms["maxPSCLoop"], minFluxPSC=parms["minFluxPSC"], \
                          solPInt=parms["solPInt"], solPMode=parms["solPMode"], solPType=parms["solPType"], \
                          maxASCLoop=parms["maxASCLoop"], minFluxASC=parms["minFluxASC"], \
                          solAInt=parms["solAInt"], solAMode=parms["solAMode"], solAType=parms["solAType"], \
                          avgPol=parms["avgPol"], avgIF=parms["avgIF"], minSNR = 4.0, refAnt=parms["refAnt"], \
                          do3D=parms["do3D"], BLFact=parms["BLFact"], BLchAvg=parms["BLchAvg"], \
                          doMB=parms["doMB"], norder=parms["MBnorder"], maxFBW=parms["MBmaxFBW"], \
                          PBCor=parms["PBCor"],antSize=parms["antSize"], \
                          nTaper=parms["nTaper"], Tapers=parms["Tapers"], \
                          doOutlier=parms["doOutlier"], OutlierDist=parms["OutlierDist"], OutlierFlux=parms["OutlierFlux"], \
                          nThreads=nThreads, noScrat=noScrat, logfile=logFile, check=check, debug=debug)
        # End image
    
    # Get report on sources
    if parms["doReport"]:
        # If targets not specified, do all
        if len(parms["targets"])<=0:
            slist = ALMAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
        else:
            slist = parms["targets"]
        Report = ALMAReportTargets(uv, err, Sources=slist, seq=parms["seq"], sclass=outIClass, \
                                       Stokes=parms["Stokes"], logfile=logFile, check=check, debug=debug)
        # Save to pickle jar
        ReportPicklefile = "./"+fileRoot+"Report.pickle"   # Where results saved
        SaveObject(Report, ReportPicklefile, True) 
       
    # Write results, cleanup    
    # Save cal/average UV data? 
    if parms["doSaveUV"] and (not check):
        Aname = ALMAAIPSName(project, session)
        cno = AIPSDir.PTestCNO(disk, user, Aname, avgClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS CAL UV DATA", Aname, avgClass, disk, parms["seq"], True, err)
            filename = parms["project"]+parms["session"]+parms["band"]+"Cal.uvtab"
            fuv = ALMAUVFITS (uvt, filename, 0, err, compress=parms["Compress"], logfile=logFile)
            ALMAAddOutFile( filename, 'project', "Calibrated Averaged UV data" )
            # Save list of output files
            ALMASaveOutFiles()
            del uvt
    # Save raw UV data tables?
    if parms["doSaveTab"] and (not check):
        Aname = ALMAAIPSName(project, session)
        cno = AIPSDir.PTestCNO(disk, user, Aname, dataClass[0:6], "UV", parms["seq"], err)
        if cno>0:
            uvt = UV.newPAUV("AIPS RAW UV DATA", Aname, dataClass[0:6], disk, parms["seq"], True, err)
            filename = parms["project"]+parms["session"]+parms["band"]+"CalTab.uvtab"
            fuv = ALMAUVFITSTab (uvt, filename, 0, err, logfile=logFile)
            ALMAAddOutFile( filename, 'project', "Calibrated AIPS tables" )
            del uvt
            # Write History
            filename = project+'_'+session+'_'+band+".History.text"
            OTObit.PrintHistory(uv, file=filename)
            ALMAAddOutFile( filename, 'project', "Processing history of calibrated data" )
            # Save list of output files
            ALMASaveOutFiles()
    # Imaging results
    # If targets not specified, save all
    if len(parms["targets"])<=0:
        slist = ALMAAllSource(uv,err,logfile=logFile,check=check,debug=debug)
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
                xf = ALMAImFITS (x, outfile, 0, err, logfile=logFile)
                ALMAAddOutFile( outfile, target, 'Image of '+ target)
                # Statistics
                zz=imstat(x, err, logfile=logFile)
    # end writing loop
    
    # Save list of output files
    ALMASaveOutFiles()
    OErr.printErrMsg(err, "Writing output")
    
    # Contour plots
    if parms["doKntrPlots"]:
        mess = "INFO --> Contour plots (doKntrPlots)"
        printMess(mess, logFile)
        ALMAKntrPlots( err, imName=parms["targets"], project=project,
            session=session, band=band, disk=disk, debug=debug )
        # Save list of output files
        ALMASaveOutFiles()
    elif debug:
        mess = "Not creating contour plots ( doKntrPlots = "+str(parms["doKntrPlots"])+ " )"
        printMess(mess, logFile)
    
    # Source uv plane diagnostic plots
    if parms["doDiagPlots"]:
        mess = "INFO --> Diagnostic plots (doDiagPlots)"
        printMess(mess, logFile)
        # Get the highest number avgClass catalog file
        Aname = ALMAAIPSName( project, session )
        uvc = None
        if not check:
            uvname = project+"_"+session+"_"+band+"_Cal"
            uvc = UV.newPAUV(uvname, Aname, avgClass, disk, parms["seq"], True, err)
        ALMADiagPlots( uvc, err, cleanUp=parms["doCleanup"], \
                           project=project, session=session, band=band, \
                           logfile=logFile, check=check, debug=debug )
        # Save list of output files
        ALMASaveOutFiles()
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
            Aname = ALMAAIPSName(project, session)
            uvname = project+"_"+session+"_"+band+"_Cal"
            uvc = UV.newPAUV(uvname, Aname, avgClass, disk, parms["seq"], True, err)
            if err.isErr:
                OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    
        # Get source metadata; save to pickle file
        srcMetadata = ALMASrcMetadata( uvc, err, Sources=parms["targets"], seq=parms["seq"], \
                                       sclass=outIClass, Stokes=parms["Stokes"],\
                                       logfile=logFile, check=check, debug=debug )
        picklefile = "./"+fileRoot+".SrcReport.pickle" 
        SaveObject( srcMetadata, picklefile, True ) 
        ALMAAddOutFile( picklefile, 'project', 'All source metadata' )
    
        # Get project metadata; save to pickle file
        projMetadata = ALMAProjMetadata( uvc, AIPS_VERSION, err, \
            PCals=parms["PCals"], ACals=parms["ACals"], \
            BPCals=parms["BPCals"], DCals=parms["DCals"], \
            project = project, session = session, band = band, \
            dataInUVF = parms["archRoots"][0], archFileID = 66666 )
        picklefile = "./"+fileRoot+".ProjReport.pickle"
        SaveObject(projMetadata, picklefile, True) 
        ALMAAddOutFile( picklefile, 'project', 'Project metadata' )
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
        ALMAHTMLReport( projMetadata, srcMetadata, \
                            outfile="./"+fileRoot+".report.html", \
                            logFile=logFile )
    
    # Write VOTable
    if parms["doVOTable"]:
        mess = "INFO --> Write VOTable (doVOTable)"
        printMess(mess, logFile)
        ALMAAddOutFile( 'VOTable.xml', 'project', 'VOTable report' ) 
        ALMAWriteVOTable( projMetadata, srcMetadata, filename='VOTable.xml' )
    
    # Save list of output files
    ALMASaveOutFiles()
    
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
            oclass = parms["Stokes"][istok:istok+1]+outCClass[1:]
            AllDest(err, disk=disk,Aseq=parms["seq"],Aclass=oclass)
        
        # Delete initial UV data
        Aname = ALMAAIPSName(project, session)
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
