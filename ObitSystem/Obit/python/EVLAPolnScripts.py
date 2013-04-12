"""
Scripts for polarization calibration for EVLA-like arrays with circular feeds
"""

import ObitTask, AIPSDir, AIPS
import EVLACal
import OSystem
from PipeUtil import *

def EVLAPolnFlag (uv, souModels, err, \
                  flagTab=1, timeAvg=0.0, maxAmp=1.0, \
                  doCalib=0, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
                  doPol=False, PDVer=1, nThreads=1, \
                  check=False, debug = False, noScrat=[], logfile = ""):
    """
    Subtract a polarization model from a data set and flag by residuals
    Data flagged using AutoFlag with points larger than maxAmp

    Returns task error code, 0=OK, else failed
    * uv       = UV data object to calibrate
    * souModels= List of Source Models as created by EVLACal.EVLACalModel
    * err      = Obit error/message stack
    * flagTab  = Output Flagging table version
    * timeAvg  = Time over which data is vector averaged
    * maxAmp   = maximum allowed residual, single or list per source
    #* flagSig  = Flagging level (sigma)
    #* alpha    = Smoothing parameter
    #* timeWind = MWF width (min)
    #* avgTime  = pre average time (min)
    #* avgFreq  = Averaging in frequency
    #             0 = none
    #             1 = average chAvg channels
    #             2 = average all channels
    #             3 = average all channels, IFs
    * doCalib  = Apply calibration table
    * chAvg    = Number of channels to average
    * gainUse  = CL/SN table to apply
    * doBand   = If >0.5 apply bandpass cal.
    * BPVer    = Bandpass table version
    * doPol    = If True apply poln cal
    * PDVer    = AIPS PD polarization cal table to apply
    * flagVer  = Input Flagging table version
    * check    = Only check script, don't execute tasks
    * debug    = Run tasks debug, show input
    * nThreads = Number of threads to use
    * noScrat  = list of disks to avoid for scratch files
    * logfile  = Log file for task
    """
    if len(souModels)<=0:    # Anything to work on?
        return 0
    mess =  "Flag residuals from models "
    printMess(mess, logfile)

    # Setup UVSub
    uvsub = ObitTask.ObitTask("UVSub")
    if not check:
        setname(uv,uvsub)
    uvsub.user     = AIPS.userno   # This gets lost somehow
    uvsub.taskLog  = logfile
    uvsub.nThreads = nThreads
    uvsub.noScrat  = noScrat
    uvsub.outDisk  = uv.Disk
    uvsub.flagVer  = flagVer
    uvsub.doCalib  = doCalib
    uvsub.gainUse  = gainUse
    uvsub.doBand   = doBand
    uvsub.BPVer    = BPVer
    if "PDVer" in uvsub.__dict__:
        uvsub.doPol    = doPol
        uvsub.PDVer    = PDVer
    # setup for flagging
    #flag = ObitTask.ObitTask("MednFlag")
    flag = ObitTask.ObitTask("AutoFlag")
    flag.user = AIPS.userno   # This gets lost somehow
    if not check:
        setoname(uv,flag)
    flag.taskLog  = logfile
    flag.nThreads = nThreads
    flag.flagTab  = flagTab
    flag.timeAvg  = timeAvg
    if type(maxAmp)!=list:
        flag.IClip    = [maxAmp, 0.1]
    #flag.flagSig  = flagSig
    #flag.alpha    = alpha
    #flag.timeWind = timeWind
    #flag.avgTime  = avgTime
    #flag.avgFreq  = avgFreq
    #flag.chAvg    = chAvg
    # Loop over sources
    OK = False
    i = -1
    for model in souModels:
        i += 1
        uvsub.Sources[0]= model["Source"]
        if "DataType2" in uvsub.__dict__:
            uvsub.DataType2 = model["CalDataType"]
        uvsub.in2File   = model["CalFile"]
        uvsub.in2Name   = model["CalName"][0:12]
        uvsub.in2Class  = model["CalClass"]
        uvsub.in2Seq    = model["CalSeq"] 
        uvsub.in2Disk   = model["CalDisk"]
        if "nmaps" in uvsub.__dict__:
            uvsub.nmaps     = model["CalNfield"]
        else:
            uvsub.nfield    = model["CalNfield"]
        uvsub.CCVer     = model["CalCCVer"]
        uvsub.BComp     = model["CalBComp"]
        uvsub.EComp     = model["CalEComp"]
        uvsub.Cmethod   = model["CalCmethod"]
        uvsub.Cmodel    = model["CalCmodel"]
        uvsub.Flux      = model["CalFlux"]
        if "Alpha" in uvsub.__dict__:
            uvsub.Alpha     = model["CalModelSI"]
        uvsub.modelFlux = model["CalModelFlux"]
        uvsub.modelPos  = model["CalModelPos"]
        uvsub.modelParm = model["CalModelParm"]
        # Setup scratch output file
        uvsub.outName  = model["Source"][0:12]
        uvsub.outClass = "PolFlg"
        uvsub.outSeq   = 6666
        if debug:
            uvsub.i
            uvsub.debug = debug
        # Trap failure
        try:
            mess = "Run UVSub on "+uvsub.Sources[0]+", then edit residuals above "+str(maxAmp[i])
            printMess(mess, logfile)
            if not check:
                uvsub.g
        except Exception, exception:
            print exception
            mess = "UVSub Failed retCode= "+str(uvsub.retCode)+" Source "+uvsub.Sources[0]
            printMess(mess, logfile)
            return uvsub.retCode
        else:
            OK = True
        # Now flag
        cno = AIPSDir.PTestCNO(uvsub.outDisk, AIPS.userno, \
                               uvsub.outName, uvsub.outClass, "UV", uvsub.outSeq, err)
        if cno>0:
            uvscr = UV.newPAUV("zap", uvsub.outName, uvsub.outClass, uvsub.outDisk, uvsub.outSeq, False, err)
        else:
            mess = "No subtracted data created for Source "+uvsub.Sources[0]
            printMess(mess, logfile)
            return 1
        # Flag using scratch writing input FG
        flag.Sources[0]= uvsub.Sources[0]
        if type(maxAmp)==list:
            flag.IClip    = [maxAmp[i], 0.1]
        if not check:
            setname(uvscr, flag)
        if debug:
            flag.i
            flag.debug = debug
        # Trap failure
        try:
            if not check:
                flag.g
        except Exception, exception:
            print exception
            mess = "AutoFlag Failed retCode= "+str(flag.retCode)+" Source "+flag.Sources[0]
            printMess(mess, logfile)
            return flag.retCode
        # Zap scratch
        uvscr.Zap(err); del uvscr
        
        # End loop over models
    return 0;  # Most be OK if it gets here
# end EVLAPolnFlag

def EVLAPolnSelfCal(uv, Cals, err, \
                    doCalib=-1, gainUse=0, doBand=0, BPVer=0, flagVer=-1,
                    doPol=False, PDVer=0, avgPol=False, avgIF=False, \
                    FQid=0, BChan=1, EChan=0, solMode="A&P", solType="L1", minSNR=5., \
                    solnver=0, solInt=10.0/60.0, refAnt=0, ampScalar=False, \
                    calIn=0, calOut=0, FOV=0.0, \
                    nThreads=1, check=False, debug = False, noScrat=[], logfile = ""):
    """
    Self Calibration using a set of models
    
    Returns task error code, 0=OK, else failed

    * uv       = UV data object to calibrate
    * Cals     = List of calibrators as created by EVLACal.EVLACalMode
                 If the Cals entry contains the following, it overrides the defaults
                 solInt  = solution interval (min)
                 solMode = solution mode "P", "DELA", "A&P"
                 minSNR  = minimum SNR
                 avgPol  = Average polns?
    * err      = Obit error/message stack
    * FQid     = Frequency Id to process, 0=>any
    * BChan    = First (1-rel channel to include
    * EChan    = Highest channel to include
    * doCalib  = Apply calibration table, positive=>calibrate
    * gainUse  = CL/SN table to apply
    * doBand   = If >0.5 apply previous bandpass cal.
    * BPVer    = Bandpass table (BP) version
    * flagVer  = Flagging table to apply
    * doPol    = Calibrate polarization?
    * PDVer    = Poln cal table (PD) version
    * solnver  = output SN table version (+1 if smooth), 0=>new
    * solInt   = solution interval (min)
    * solMode  = solution mode
    * solType  = solution type
    * refAnt   = Reference antenna
    * avgPol   = Average poln in SC?
    * avgIF    = Average IFs in SC?
    * ampScalar= If true, scalar average data in calibration?
    * calIn    = Input CL table for CLCal, 0=> gainUse or highest if 0
    * calOut   = Output CL table for CLCal, 0 => new
    * FOV      = If >0.0 then do baseline dependent averaging before running Calib
                 where FOV is the desired undistorted FOV (deg)
    * check    = Only check script, don't execute tasks
    * nThreads = Number of threads to use
    * debug    = Run tasks debug, show input
    * noScrat  = list of disks to avoid for scratch files
    * logfile  = Log file for tasks
    """
    ################################################################
    mess =  "Self Calibration with model"
    printMess(mess, logfile)
    solnVer2 = None
    tmpfile  = None  # Averaged data file, if used

    # output SN version
    if solnver<=0:
        # Open and close image to sync with disk 
        uv.Open(UV.READONLY, err)
        uv.Close(err)
        solnVer  = max(1,uv.GetHighVer("AIPS SN")+1)
    else:
        solnVer  = solnver

    # Setup Calib 
    calib = ObitTask.ObitTask("Calib")
    calib.user     = AIPS.userno   # This gets lost somehow
    calib.taskLog  = logfile
    if not check:
        setname(uv,calib)
    calib.flagVer  = flagVer
    calib.doCalib  = doCalib
    calib.gainUse  = gainUse
    calib.doBand   = doBand
    calib.BPVer    = BPVer
    calib.doPol    = doPol
    calib.avgIF    = avgIF
    calib.avgPol   = avgPol
    calib.minSNR   = minSNR
    calib.PDVer    = PDVer
    calib.ampScalar= ampScalar
    calib.solType  = solType
    calib.nThreads = nThreads
    calib.refAnts  = [refAnt]
    calib.noScrat  = noScrat
    calib.solnVer  = solnVer
    calib.prtLv    = 2
    OK      = False   # Must have some work
    OKCals2 = []      # List of ones that worked

    # Setup CLCal
    # Default CL tables
    if not check:
        hiCL   = uv.GetHighVer("AIPS CL")
        if calIn<=0:
            calIn = hiCL
        if calOut<=0:
            calOut = hiCL+1
    clcal = ObitTask.ObitTask("CLCal")
    clcal.user     = AIPS.userno   # This gets lost somehow
    clcal.taskLog  = logfile
    ical = 0
    if not check:
        setname(uv,clcal)
    clcal.solnVer   = solnVer
    clcal.calIn     = calIn
    clcal.calOut    = calOut
    clcal.interMode = "SELF" 
    clcal.FreqID    = FQid
    clcal.refAnt    = refAnt

    # Doing baseline average?
    if FOV>0.0:
        avg = ObitTask.ObitTask("UVBlAvg")
        avg.user       = AIPS.userno   # This gets lost somehow
        avg.taskLog    = logfile
        if not check:
            setname(uv,avg)
        avg.flagVer    = flagVer
        avg.gainUse    = gainUse
        avg.doBand     = doBand
        avg.BPVer      = BPVer
        if "PDVer" in avg.__dict__:
            avg.doPol  = doPol
            avg.PDVer  = PDVer
        avg.FOV        = FOV
        avg.maxInt     = calib.solInt
        avg.maxFact    = 1.01
        avg.outFile    = ("TMP"+calib.inFile).strip()
        avg.outName    = avg.inName
        avg.outClass   = "Temp"
        avg.outSeq     = 9999
        avg.outDisk    = calib.inDisk
        # Add sources
        avg.Sources = []
        for Cal in Cals:
            avg.Sources.append(Cal["Source"])
        # Trap failure
        try:
            mess = "Run UVBlAvg on "+avg.Sources[0]
            printMess(mess, logfile)
            if not check:
                avg.g
        except Exception, exception:
            print exception
            mess = "UVBlAvg Failed retCode= "+str(avg.retCode)+" Source "+avg.Sources[0]
            printMess(mess, logfile)
            tmpfile = None
            # Keep on Truckin'
        else:
            # Worked - use it
            if calib.DataType=="AIPS":
                tmpfile = UV.newPAUV("Avg", avg.outName, avg.outClass, avg.outDisk, avg.outSeq, True, err)
            else:
                tmpfile = UV.newPFUV("Avg", avg.outFile, avg.outDisk, True, err)
                
            if err.isErr:
                OErr.printErrMsg(err, "Error with UVBlAvg output")
            # No edit/cal in Calib
            if not check:
                setname(tmpfile, calib)
            calib.flagVer  = -1
            calib.doCalib  = -1
            calib.doBand   = -1
            calib.doPol    = False
            calib.solnVer  = 1
    # End if baseline dependent averaging
    
        
    # Loop over calibrators
    for Cal in Cals:
        calib.Sources[0]= Cal["Source"]
        calib.DataType2 = Cal["CalDataType"]
        calib.in2File   = Cal["CalFile"]
        calib.in2Name   = Cal["CalName"][0:12]
        calib.in2Class  = Cal["CalClass"]
        calib.in2Seq    = Cal["CalSeq"] 
        calib.in2Disk   = Cal["CalDisk"]
        calib.nfield    = Cal["CalNfield"]
        calib.CCVer     = Cal["CalCCVer"]
        calib.BComp     = Cal["CalBComp"]
        calib.EComp     = Cal["CalEComp"]
        calib.Cmethod   = Cal["CalCmethod"]
        calib.Cmodel    = Cal["CalCmodel"]
        calib.Flux      = Cal["CalFlux"]
        calib.Alpha     = Cal["CalModelSI"]
        calib.modelFlux = Cal["CalModelFlux"]
        calib.modelPos  = Cal["CalModelPos"]
        calib.modelParm = Cal["CalModelParm"]
        calib.solnVer   = solnVer
        # Source dependent overrides?
        if "solMode" in Cal:
            calib.solMode  = Cal["solMode"]
        else:
            calib.solMode  = solMode
        if "solInt" in Cal:
            calib.solInt   = Cal["solInt"]
        else: 
            calib.solInt   = solInt
        if "minSNR" in Cal:
            calib.minSNR   = Cal["minSNR"]
        else: 
            calib.minSNR   = minSNR
        if "avgPol" in Cal:
            calib.avgPol   = Cal["avgPol"]
        else: 
            calib.avgPol   = avgPol
        if debug:
            calib.i
            calib.debug = debug
        #calib.prtLv = 5
        # Trap failure
        try:
            mess = "Run Calib on "+calib.Sources[0]
            printMess(mess, logfile)
            if not check:
                calib.g
        except Exception, exception:
            print exception
            mess = "Calib Failed retCode= "+str(calib.retCode)+" Source "+calib.Sources[0]
            printMess(mess, logfile)
            #return 1  # allow some failures
        else:
            OK = True
            OKCals2.append(calib.Sources[0])
            # Were averaged data used?  If so,copy SN table contents
            if FOV>0.0 and tmpfile:
                tmpfile.Open(UV.READONLY,err)  # Update header
                tmpfile.Close(err)
                EVLAPolnSNAppend (tmpfile, solnVer, uv, solnVer, err)
                # Delete SN table on averaged data
                tmpfile.ZapTable("AIPS SN", solnVer, err)
            
            # Run ClCal for this source
            clcal.Sources[0] = calib.Sources[0]
            clcal.calSour[0] = calib.Sources[0]
            if debug:
                clcal.i
                clcal.debug = debug
            #calib.prtLv = 5
            # Trap failure
            try:
                if not check:
                    clcal.g
            except Exception, exception:
                print exception
                mess = "CLCal Failed retCode= "+str(clcal.retCode)+" Source "+clcal.Sources[0]
                printMess(mess, logfile)
                return 1;
            
    # end calibration loop
    # Were averaged data used?  Zap temp file
    if FOV>0.0 and tmpfile:
        tmpfile.Zap(err); del tmpfile; tmpfile = None  # Delete scratch
    return 0  # Must be OK
# end EVLAPolnSelfCal

def EVLAPolnPCal(uv, InsCals, err, \
                 doCalib=2, gainUse=0, doBand=-1, BPVer=0, BPSoln=0, flagVer=1, \
                 solType="    ", solInt=0.0, refAnt=0, polVer = 1, ChInc=1, ChWid=1, \
                 RLPhase=[-999.,], RM=[0.0], PPol=[0.], \
                 check=False, debug = False, \
                 nThreads=1, noScrat=[], logfile = ""):
    """
    Instrumental Polarization 
    
    Do Instrumental
    Instrumental cal uses PCal
    Returns task error code, 0=OK, else failed

    * uv       = UV data object to calibrate
    * InsCals  = Instrumental poln calibrators, name or list of names
      If None no instrumental cal
    * err      = Obit error/message stack
    * doCalib  = Apply prior calibration table, positive=>calibrate
    * gainUse  = CL/SN table to apply
    * doBand   = >0 => apply bandpass calibration
    * BPVer    = AIPS BP table to apply
    * BPSoln   = Output AIPS BP table, 0=>new
    * flagVer  = Input Flagging table version
    * solType  = solution type, "    ", "LM  "
    * PPol     = list of R-L fractional polarizations for fixed source linear poln.
    * RLPhase  = list of R-L phases (deg) for fixed source linear poln.
    * RM       = NYI Rotation measure of fixed source
    * solInt   = instrumental solution interval (min)
    * refAnt   = Reference antenna
    * polVer   = version number of PD, CP (BP) tables
    * ChInc    = channel increment for solutions
    * ChWid    = number of channels to average for solution.
    * nThreads = Number of threads to use in imaging
    * check    = Only check script, don't execute tasks
    * debug    = Run tasks debug, show input
    * noScrat  = list of disks to avoid for scratch files
    * logfile  = Log file for task
    """
    ################################################################
    # Don't bother if not full polarization
    d     = uv.Desc.Dict
    nstoke = d["inaxes"][d["jlocs"]]
    if nstoke<4:
        mess = "Skip Instrumental polarization corrections - not full stokes"
        printMess(mess, logfile)
        return 0
    mess =  "Instrumental polarization calibration "
    printMess(mess, logfile)
    # Instrumental calibration
    if InsCals!=None:
        pcal         = ObitTask.ObitTask("PCal")
        pcal.user    = AIPS.userno   # This gets lost somehow
        pcal.logFile = logfile
        if not check:
            setname(uv,pcal)
        if type(InsCals)==list:
            pcal.Sources = InsCals
        else:
            pcal.Sources = [InsCals]
        pcal.doCalib  = doCalib
        pcal.gainUse  = gainUse
        pcal.doBand   = doBand
        pcal.BPVer    = BPVer
        pcal.flagVer  = flagVer
        pcal.solnType = solType
        pcal.solInt   = solInt
        pcal.PPol     = PPol
        # If any set, fit R-L phase, not fit lin poln
        i=0
        for p in PPol:
            pcal.doFitI[i]    = True
            if p>0.0:
                pcal.doFitRL = True
                pcal.doFitPol[i] = False
            i += 1
        pcal.RLPhase  = RLPhase
        pcal.ChInc    = ChInc
        pcal.ChWid    = ChWid
        pcal.CPSoln   = polVer
        pcal.PDSoln   = polVer
        pcal.BPSoln   = BPSoln
        pcal.refAnt   = refAnt
        pcal.prtLv    = 2
        pcal.nThreads = nThreads
        for i in range(0,len(pcal.doFitI)):
            pcal.doFitI[i]   = True
        pcal.taskLog  = logfile
        i = 1;
        for d in noScrat:
            pcal.noScrat[i] = d
            i += 1
        if debug:
            pcal.i
            pcal.debug = debug
        # Trap failure
        try:
            if not check:
                pcal.g
        except Exception, exception:
            print exception
            mess = "PCal Failed retCode="+str(pcal.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
    # end instrumental poln cal
    
    return 0
    # End EVLAPolnPCal

def EVLAPolnRL(uv, err, \
              RLCal=None, BChan=1, EChan = 0, ChWid2=1, solInt1=1./6, solInt2=10., \
              UVRange=[0.0,0.0], timerange = [0.0,1000.0],  \
              FQid=0,  doCalib=-1, gainUse=0, \
              doBand=0, BPVer=0, BPSoln=0, flagVer=-1, refAnt=0, doPol=False, PDVer=1, 
              nThreads=1, noScrat=[], logfile = "",check=False, debug = False):
    """
    Determine R-L  phase calibration as BP table
    
    Returns task error code, 0=OK, else failed
    R-L Delay calibration using new BP table, if R-L phase (& RM) known for

    * uv       = UV data object to clear
    * err      = Obit error/message stack
    * RLCal    = list of R-L (polarization angle) calibrator as
                 (name, R-L phase (deg at 1 GHz), RM (rad/m**2))
    * solInt1  = first solution interval (min), 0=> scan average
    * solInt2  = second solution interval (min)
    * BChan    = First (1-rel) channel number
    * EChan    = Highest channel number. 0=> high in data.
    * ChWid2   = Number of channels in running mean RL BP soln, 
    * UVRange  = Range of baseline used in kilowavelengths
    * FQid     = Frequency Id to process
    * doCalib  = Apply calibration table, positive=>calibrate
    * gainUse  = CL/SN table to apply
    * timerange= time range of data (days)
    * doBand   = If >0.5 apply previous bandpass cal.
    * BPVer    = previous Bandpass table (BP) version
    * BPSoln   = Output Bandpass table (BP) version
    * flagVer  = Flagging table to apply
    * refAnt   = Reference antenna REQUIRED
    * doPol    = Apply polarization cal?
    * PDVer    = PD version for pol cal, -1=>use IF
    * noScrat  = list of AIPS disks to avoid for scratch files
    * nThreads = Number of threads to use in imaging
    * logfile  = Log file for task
    * check    = Only check script, don't execute tasks
    * debug    = Run tasks debug, show input
    """
    ################################################################
    # Don't bother if not full polarization
    d     = uv.Desc.Dict
    nstoke = d["inaxes"][d["jlocs"]]
    if nstoke<4:
        mess = "Skip R-L polarization corrections - not full stokes"
        printMess(mess, logfile)
        return 0
    mess =  "R-L polarization calibration "
    printMess(mess, logfile)

    lbpver = BPVer  # default bandpass in imaging
    # Want R-L phase cal ?
    if RLCal!=None:
        ncal = len(RLCal)  # How many calibrators?
        rlpass=ObitTask.ObitTask("RLPass")
        rlpass.user    = AIPS.userno   # This gets lost somehow
        rlpass.taskLog = logfile
        if not check:
            setname(uv,rlpass)
            #if Antennas:
            #    i = 0
            #    for a in Antennas:
            #        rlpass.Antennas[i] = a; i  += 1
        rlpass.timeRange[0] = timerange[0]
        rlpass.timeRange[1] = timerange[1]
        rlpass.BChan1  = BChan
        rlpass.EChan1  = EChan
        rlpass.BChan2  = BChan
        rlpass.EChan2  = EChan
        rlpass.ChWid2  = ChWid2
        rlpass.UVRange[0] = UVRange[0];
        rlpass.UVRange[1] = UVRange[1];
        rlpass.doCalib = doCalib
        rlpass.gainUse = gainUse
        rlpass.flagVer = flagVer
        rlpass.FreqID  = FQid
        rlpass.doPol   = doPol
        if "PDVer" in rlpass.__dict__:
            rlpass.PDVer = PDVer
        rlpass.doBand  = doBand
        rlpass.BPVer   = BPVer
        rlpass.refAnt  = refAnt
        rlpass.solInt1 = solInt1
        rlpass.solInt2 = solInt2
        rlpass.BPSoln  = 0
        rlpass.prtLv   = 1
        rlpass.nThreads = nThreads
        rlpass.noScrat  = noScrat
        # Loop over calibrators
        for ical in range (0,ncal):
            rlpass.Sources[0]= RLCal[ical][0]
            rlpass.RLPhase   = RLCal[ical][1]
            rlpass.RM        = RLCal[ical][2]
            mess =  "R-L channel phase calibration using "+rlpass.Sources[0]
            printMess(mess, logfile)
            if debug:
                print "timerange", rlpass.timerang
                rlpass.i
                rlpass.debug = True
            # Trap failure
            try:
                if not check:
                    rlpass.g
            except Exception, exception:
                print exception
                mess = "rlpass Failed retCode="+str(rlpass.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
            # end loop over calibrators
        # Get output BP table
        # Open and close image to sync with disk 
        uv.Open(UV.READONLY, err)
        uv.Close(err)
    # end R-L delay cal
    return 0
# end EVLAPolnRL

def EVLAPolnImage(uv, err, Sources=None,  FreqID=1, seq=1, sclass="IClean", band="", \
                  doCalib=-1, gainUse=0, doBand=-1, BPVer=0,  flagVer=-1, maxPixel=50000, \
                  doPol=False, PDVer=-1,  minFlux=0.0, BIF=1, EIF=0, \
                  Stokes="I", FOV=20.0/3600.0, Robust=0, Niter=300, CleanRad=None, \
                  maxPSCLoop=0, minFluxPSC=0.1, solPInt=20.0/60.,\
                  solPMode="P", solPType= "L1", Gain=0.1, \
                  maxASCLoop=0, minFluxASC=0.5, solAInt=2.0, \
                  solAMode="A&P", solAType= "L1",  Beam=[0.0,0.0,0.0], \
                  avgPol=False, avgIF=False, minSNR = 5.0, refAnt=0, \
                  do3D=False, BLFact=0.999, BLchAvg=False, doOutlier=None, \
                  doMB=False, norder=1, maxFBW=0.05, doComRes=False, \
                  nTaper=0, Tapers=[20.0], dispURL="None", \
                  nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """
    Image a list of sources with optional selfcal
    
    Uses Imager or MFImage to image a list of sources.
    Data must be at least approximately calibrated
    Returns list of source models, None on failure

    * uv         = UV data object
    * err        = Python Obit Error/message stack
    * Sources    = Source name or list of names to use
      If an empty list all sources in uv are included
    * seq        = sequence number of output
    * sclass     = Image output class
    * band       = project band - appended to name
    * BIF        = First IF
    * CIF        = Highest IF
    * FreqID     = Frequency group identifier
    * doCalib    = Apply calibration table
    * gainUse    = CL/SN table to apply
    * doBand     = If >0.5 apply bandpass cal.
    * BPVer      = Bandpass table version
    * flagVer    = Input Flagging table version
    * doPol      = Apply polarization cal?
    * PDVer      = PD version for pol cal, -1=>use IF
    * minFlux    = minimum flux density for initial CLEAN
    * Stokes     = Stokes parameters to image
    * FOV        = Field of view to image in deg
    * Robust     = Weighting robustness parameter
    * Niter      = max no. iterations
    * Gain       = CLEAN loop gain
    * maxPixel   = Maximum number of residuals in inner CLEAN
    * CleanRad   = CLEAN radius about center or None=autoWin
    * maxPSCLoop = max. number of phase sc loops
    * minFluxPSC = Trip level for phase self cal (Jy)
    * solPInt    = Phase solution interval (min)
    * solPMode   = Phase soln mode "P", "DELA"
    * solPType   = Phase soln type
    * maxASCLoop = max. number of amp&phase sc loops
    * minFluxASC = Trip level for amp&phase self cal (Jy)
    * solAInt    = Amp&phase solution interval (min)
    * solAMode   = Amp&Phase soln mode "A&P", "P", "DELA"
    * solAType   = Amp&PPhase soln type
    * Beam       = CLEAN restoring beam
    * avgPol     = Average poln in SC?
    * avgIF      = Average IFs in SC?
    * minSNR     = minimum acceptable SNR in SC
    * refAnt     = Reference antenna
    * do3D       = Use 3D imaging?
    * doComRes   = Force common resolution in frequency? (MF)
    * BLFact     = Baseline dependent averaging factor
    * BLchAvg    = If True and BLFact>=1.0 also average channels
    * doOutlier  = Outliers from NVSS?  Yes=> 4*FOV, 1 mJy
                   None = Default, Yes if freq<6 GHz
    * doMB       = If True is wideband imaging
    * norder     = order on wideband imaging
    * maxFBW     = max. fractional wideband imaging
    * nTaper     = number of (additional) multi resolution tapers
    * Tapers     = Sizes of additional tapers in pixels
    * dispURL    = Name of image display
    * nThreads   = Max. number of threads to use
    * noScrat    = list of disks to avoid for scratch files
    * logfile    = logfile for messages
    * check      = Only check script, don't execute tasks
    * debug      = show input
    """
    ################################################################
    mess = "Image a list of sources "
    printMess(mess, logfile)

    # Tolerate missing BP table
    # Open and close image to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)
    hiBP = uv.GetHighVer("AIPS BP")
    if hiBP<=0:
        doBand = -1
    # get reference Freq
    refFreq = uv.Desc.Dict["crval"][uv.Desc.Dict["jlocf"]]
    # If list empty get all sources
    if type(Sources)==list:
        sl = Sources
    else:
        sl = [Sources]

    if len(sl)<=0:
        slist = EVLAAllSource(uv,err,logfile=logfile,check=check,debug=debug)
    else:
        slist = sl
    if doMB:
        imager = ObitTask.ObitTask("MFImage")
        imager.norder = norder
        imager.maxFBW = maxFBW
        imager.prtLv = 2
    else:
        imager = ObitTask.ObitTask("Imager")
        imager.prtLv = 2
    imager.user = AIPS.userno   # This gets lost somehow
    imager.taskLog  = logfile
    if not check:
        setname(uv,imager)
    imager.outDisk     = imager.inDisk
    imager.outName     = "_"+band
    imager.out2Name    = "_"+band
    imager.out2Disk    = imager.inDisk
    imager.outSeq      = seq
    imager.out2Seq     = seq
    imager.outClass    = sclass
    imager.BIF         = BIF
    imager.EIF         = EIF
    imager.BLFact      = BLFact
    imager.BLchAvg     = BLchAvg
    imager.flagVer     = flagVer
    imager.doCalib     = doCalib
    imager.gainUse     = gainUse
    imager.doBand      = doBand
    imager.BPVer       = BPVer
    imager.doPol       = doPol
    if "PDVer" in imager.__dict__:
        imager.PDVer = PDVer
    imager.Stokes      = Stokes
    imager.FOV         = FOV
    imager.Robust      = Robust
    imager.Niter       = Niter
    imager.Gain        = Gain
    imager.Beam        = Beam
    imager.maxPixel    = maxPixel
    imager.minFlux     = minFlux
    imager.maxPixel    = maxPixel
    imager.maxPSCLoop  = maxPSCLoop
    imager.minFluxPSC  = minFluxPSC
    imager.solPInt     = solPInt
    imager.solPMode    = solPMode
    imager.solPType    = solPType
    imager.maxASCLoop  = maxASCLoop
    imager.minFluxASC  = minFluxASC
    imager.solAInt     = solAInt
    imager.solAMode    = solAMode
    imager.solAType    = solAType
    imager.avgPol      = avgPol
    imager.avgIF       = avgIF
    imager.refAnt      = refAnt
    imager.minSNR      = minSNR
    imager.do3D        = do3D
    imager.dispURL     = dispURL
    imager.nTaper      = nTaper
    imager.Tapers      = Tapers
    if doOutlier or ((doOutlier==None) and refFreq<6.0e9):
        FWHM = (45.0 /(refFreq*1.0e-9) ) / 60.   # FWHM in deg
        imager.OutlierDist = FWHM*2.0   # Outliers from NVSS for lower frequencies
        imager.OutlierFlux = 0.01
    # Auto window or centered box
    if CleanRad:
        imager.CLEANBox=[-1,CleanRad,0,0]
    else:
        imager.autoWindow  = True
    if "doComRes" in imager.__dict__:
        imager.doComRes  = doComRes
    imager.noScrat     = noScrat
    imager.nThreads    = nThreads
    if debug:
        imager.prtLv = 5
        imager.i
        imager.debug = debug
    OK = False   # Some must work
    modelList = []   # Init output model list
    # Loop over slist
    for sou in slist:
        imager.Sources[0] = sou
        mess = "Image "+sou
        printMess(mess, logfile)
        # Trap failure
        try:
            if not check:
                imager.g
        except Exception, exception:
            print exception
            mess = "Imager Failed retCode= "+str(imager.retCode)
            printMess(mess, logfile)
            #return 1  Allow some failures
            # Cleanup image mess
            AllDest(err,Atype="MA",Aname=imager.Sources[0], disk=imager.outDisk, Aseq=imager.outSeq);
        else:
            OK = True
        # delete Imager file if not debug
        if not debug:
            out2Name = imager.Sources[0].strip()+"_"+band
            out2Name = out2Name[0:12]
            if doMB:
                out2Class = "MFImage"
            else:
                out2Class = "Imager"
            # Tolerate failures
            try:
                # Test if file exists
                cno = AIPSDir.PTestCNO(imager.out2Disk, AIPS.userno, \
                                       out2Name, out2Class, "UV", imager.out2Seq, err)
                if cno>0:
                    u = UV.newPAUV("zap", out2Name, out2Class, imager.out2Disk, imager.out2Seq, False, err)
                    if UV.PIsA(u):
                        u.Zap(err) # cleanup
                    if err.isErr:
                        mess = "Error deleting Imager work file"
                        printMess(mess, logfile)
                        #return 1
                    del u
                # Create model, add to modelList
                modelList.append(EVLACal.EVLACalModel(sou, CalName=imager.Sources[0]+imager.outName, \
                                                      CalClass=imager.outClass,        \
                                                      CalSeq=imager.outSeq,             \
                                                      CalCmethod="DFT", CalDataType="AIPS",   \
                                                      CalDisk=imager.outDisk, CalNfield=1))
            except Exception, exception:
                print exception
                mess = "Imager Cleanup Failed source= "+imager.Sources[0].strip()+"_"+band
                printMess(mess, logfile)
                OErr.PClear(err)     # Clear any message/error
                #return 1  Allow some failures
            else:
                pass
    # end loop over sources
    # Something work?
    if not OK:
        printMess("All images failed", logfile)
        return None

    return modelList
# end EVLAPolnImage

def EVLAPolnGetFluxes(Sources, err, \
                      logfile = "",check=False, debug = False):
    """
    Determine I,Q,U flux densities for a set of calibrator models
    
    Sums CLEAN components
    Each entry in Sources contains"
    "IFlux", "QFlux", "UFlux"
    Returns 0 on success otherwise failure
    
    * Sources  = List of Source Models as created by EVLACal.EVLACalModel
    * err      = Obit error/message stack
    * logfile  = Log file for task
    * check    = Only check script
    * debug    = Run tasks debug, show input
    """
    ################################################################
    mess = "Get Source summed CLEAN components"
    printMess(mess, logfile)
    for sou in Sources:
        if sou["CalDataType"]=="AIPS":
            klass = "I"+sou["CalClass"][1:]
            img = Image.newPAImage("img", sou["CalName"][0:12], klass, \
                                   sou["CalDisk"], sou["CalSeq"], True, err)
        else:
            file = "I"+sou["CalFile"][1:]
            img = Image.newPFImage("img", file, sou["CalDisk"], True, err)
        if err.isErr:
            mess = "Error finding I image for "+sou["Source"]
            printMess(mess, logfile)
            return 1
        #img.Header(err) # DEBUG
        IFlux = EVLACal.EVLAGetSumCC(img, err, logfile=logfile)
        del img

        if sou["CalDataType"]=="AIPS":
            klass = "Q"+sou["CalClass"][1:]
            img = Image.newPAImage("img", sou["CalName"][0:12], klass, \
                                   sou["CalDisk"], sou["CalSeq"], True, err)
        else:
            file = "Q"+sou["CalFile"][1:]
            img = Image.newPFImage("img", file, sou["CalDisk"], True, err)
        if err.isErr:
            mess = "Error finding Q image for"+sou["Source"]
            printMess(mess, logfile)
            return 1
        QFlux = EVLACal.EVLAGetSumCC(img, err, logfile=logfile)
        del img
        
        if sou["CalDataType"]=="AIPS":
            klass = "U"+sou["CalClass"][1:]
            img = Image.newPAImage("img", sou["CalName"][0:12], klass, \
                                   sou["CalDisk"], sou["CalSeq"], True, err)
        else:
            file = "U"+sou["CalFile"][1:]
            img = Image.newPFImage("img", file, sou["CalDisk"], True, err)
        if err.isErr:
            mess = "Error finding U image for"+sou["Source"]
            printMess(mess, logfile)
            return 1
        UFlux = EVLACal.EVLAGetSumCC(img, err, logfile=logfile)
        del img
        sou["IFlux"] = IFlux
        sou["QFlux"] = QFlux
        sou["UFlux"] = UFlux
        mess = "Source "+sou["Source"]+" I,Q,U ="+str(IFlux)+","+str(QFlux)+","+str(UFlux)
        printMess(mess, logfile)
    # end loop over sources
    return 0
# end EVLAPolnGetFluxes

def EVLAPolnToFITS(uv, outFile, outDisk, Sources, session, err, \
                      logfile = "",check=False, debug = False):
    """
    Copy new calibration and editing tables from uv to the input FITS file

    Any SN, CL, BP, PD, CP tables not on outFile are copied from uv
    uv FG table 1 is copied to outFile table 2
    Also writes AIPS calibrator images (I, Q, U)
    Returns 0 on success otherwise failure

    * uv       = UV data with calibration and editing
    * outFile  = Name of FITS file to write tables to
    * outDisk  = FITS disk for inFile, 0=>cwd
    * Sources  = List of Source Models as created by EVLACal.EVLACalModel
    * session = Name of session, used for FITS images
    * err      = Obit error/message stack
    * logfile  = Log file for task
    * check    = Only check script
    * debug    = Run tasks debug, show input
    """
    ################################################################
    mess = "Save resultant tables, calibration images"
    printMess(mess, logfile)

    # Copy tables not on inFile
    # UV Data for infile
    outUV = UV.newPFUV("FITS UV DATA", outFile, outDisk, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with FITS data")
        
    # Use TabCopy
    taco = ObitTask.ObitTask("TabCopy")
    taco.user = AIPS.userno   # This gets lost somehow
    setname(uv, taco)
    taco.taskLog  = logfile
    taco.DataType = "AIPS"
    taco.outDType = "FITS"
    taco.outFile  = outFile
    taco.outDisk  = outDisk

    # FG
    inHi  = UV.PGetHighVer(uv, "AIPS FG")
    outHi = UV.PGetHighVer(outUV, "AIPS FG")
    if outHi<2:
        # copy uv FG table 1 to outUV FG 2 (new table)
        mess = "Copy FG 1 to FG 2"
        printMess(mess, logfile)
        taco.inTab= "AIPS FG"; taco.inVer = 1; taco.outVer = 2
        # Trap failure
        try:
            if not check:
                taco.g
        except Exception, exception:
            print exception
            mess = "TabCopy Failed retCode= "+str(taco.retCode)
            printMess(mess, logfile)
            return 1
    
    # Other tables
    tabs = ["BP", "CL", "SN", "PD", "CP"]
    for tt in tabs:
        ttype = "AIPS "+tt
        inHi  = UV.PGetHighVer(uv, ttype)
        outHi = UV.PGetHighVer(outUV, ttype)
        mess = "Copy "+tt+" tables "+str(outHi+1)+" to "+str(inHi)
        printMess(mess, logfile)
        for ver in range (outHi+1, inHi+1):
            taco.inTab=ttype; taco.inVer = ver; taco.outVer = ver
            # Trap failure
            try:
                if not check:
                    taco.g
            except Exception, exception:
                print exception
                mess = "TabCopy Failed retCode= "+str(taco.retCode)
                printMess(mess, logfile)
                return 1
    # end loop over tables

    # Write calibration images
    stokes = ["I","Q","U"]
    for sou in Sources:
        # Ignore FITS images
        if sou["CalDataType"]=="FITS":
            continue
        for s in stokes:
            klass = s+sou["CalClass"][1:]
            # Test if file exists
            cno = AIPSDir.PTestCNO(sou["CalDisk"], AIPS.userno, \
                                   sou["CalName"][0:12], klass, "MA", sou["CalSeq"], err)
            if cno>0:
                img = Image.newPAImage("img", sou["CalName"][0:12], klass, \
                                       sou["CalDisk"], sou["CalSeq"], True, err)
                if err.isErr:
                    mess = "Error finding "+s+" image for "+sou["Source"]
                    printMess(mess, logfile)
                    return 1
                imgfile = sou["Source"].strip()+session+s+"Pol.fits"
                EVLACal. EVLAImFITS (img, imgfile, outDisk, err, logfile=logfile)
                if err.isErr:
                    OErr.printErrMsg(err, "Error Writing image")
                    return 1
            else:
                mess = "Skip writing "+s+" image for "+sou["Source"]
                printMess(mess, logfile)
           # end if exist
    # end loops over Stokes and source
    return 0
# end EVLAPolnToFITS

def EVLAPolnCleanup(uv, Sources, err, \
                    logfile = "",check=False, debug = False):
    """
    Zap uv data and any calibration images in AIPS

    Returns 0 on success otherwise failure

    * uv       = UV data with calibration and editing
    * Sources  = List of Source Models as created by EVLACal.EVLACalModel
    * err      = Obit error/message stack
    * logfile  = Log file for task
    * check    = Only check script
    * debug    = Run tasks debug, show input
    """
    ################################################################
    mess = "Cleanup, Zap AIPS files, check = "+str(check)
    printMess(mess, logfile)

    # Zap data
    if UV.PIsA(uv):
        mess = "Zap uv data"
        printMess(mess, logfile)
        if not check:
            uv.Zap(err) # cleanup
            if err.isErr:
                mess = "Error Zapping uv file"
                printMess(mess, logfile)
                return 1
    
    # Zap calibration images
    stokes = ["I","Q","U"]
    for sou in Sources:
        # Ignore FITS images
        if sou["CalDataType"]=="FITS":
            continue
        for s in stokes:
            klass = s+sou["CalClass"][1:]
            # Test if file exists
            cno = AIPSDir.PTestCNO(sou["CalDisk"], AIPS.userno, \
                                   sou["CalName"][0:12], klass, "MA", sou["CalSeq"], err)
            if cno>0:
                img = Image.newPAImage("img", sou["CalName"][0:12], klass, \
                                       sou["CalDisk"], sou["CalSeq"], True, err)
                if err.isErr:
                    mess = "Error finding "+s+" image for "+sou["Source"]
                    printMess(mess, logfile)
                    return 1
                if Image.PIsA(img):
                    mess = "Zap "+sou["Source"].strip()+" "+s+"Pol"
                    printMess(mess, logfile)
                    if not check:
                        img.Zap(err) # cleanup
                        del img
                    if err.isErr:
                        mess = "Error Zapping image"
                        printMess(mess, logfile)
                        return 1
            # end if exist
    # end source/stokes loop
    return 0
# end EVLAPolnCleanup

def EVLAPolnSNAppend (inUV, inSNVer, outUV, outSNVer, err):
    """
    Append contents of one SN table to another

    * inUV      = Input UV data
    * inSNVer   = Input SN version number
    * outUV     = Output UV data
    * outSNVer  = Output SN version number
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Get tables
    inUV.Open(UV.READONLY,err)
    inUV.Close(err)
    #print "inSNVer",inSNVer
    inSN  = inUV.NewTable(Table.READONLY, "AIPS SN", inSNVer, err)
    numIF  = inSN.Desc.List.Dict["NO_IF"][2][0]
    numPol = inSN.Desc.List.Dict["NO_POL"][2][0]
    outSN = outUV.NewTable(Table.READWRITE, "AIPS SN", outSNVer, err, \
                           numIF=numIF, numPol=numPol)
    # Copy if new
    if outSN.Desc.Dict["nrow"]<=0:
        Table.PCopy(inSN, outSN,err)
    else:
        # Concatenate
        Table.PConcat(inSN, outSN,err)
    if err.isErr:
        return
    # end EVLAPolnSNAppend 

