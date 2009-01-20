""" Functions for calibrating and editing VLA data
"""
import UV, UVDesc, ObitTask, AIPSTask, OErr
from AIPS import AIPS
from FITS import FITS
from AIPSDir import AIPSdisks, nAIPS

def setname (inn, out):
    """ Copy file definition from inn to out as in...
    
    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    inn  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask.ObitTask:
        out.DataType = inn.FileType
        out.inDisk   = int(inn.Disk)
        if inn.FileType == 'FITS':
            out.inFile = inn.Fname
        else:   # AIPS
            out.inName  = inn.Aname
            out.inClass = inn.Aclass
            out.inSeq   = int(inn.Aseq)
    else:  # AIPS
        out.inname  = inn.Aname
        out.inclass = inn.Aclass
        out.inseq   = inn.Aseq
        out.indisk  = inn.Disk
    # end setname
    
def VLAClearCal(uv, err, doGain=True, doBP=False, doFlag=False):
    """ Clear previous calibration

    Delete all SN tables, all CL but CL 1
    uv       = UV data object to clear
    err      = Obit error/message stack
    doGain   = If True, delete SN and CL tables
    doBP     = If True, delete BP tables
    doFlag   = If True, delete FG tables
    """
    ################################################################
    # Gain tables
    if doGain:
        uv.ZapTable("AIPS SN",-1,err)
        ver = uv.GetHighVer("AIPS CL")
        while (ver>1):
            uv.ZapTable ('AIPS CL', ver, err)
            ver = ver-1
    # Bandpass
    if doBP:
        uv.ZapTable("AIPS BP",-1,err)
    # Flag tables
    if doFlag:
        uv.ZapTable("AIPS FG",-1,err)
    OErr.printErrMsg(err, "VLAClearCal: Error reseting calibration")
    # end VLAClearCal

def VLAMedianFlag(uv, target, err, \
                  flagTab=1, flagSig=10.0, alpha=0.5, timeWind=2.0,
                  doCalib=0, gainUse=0, doBand=0, BPVer=0, flagVer=-1):
    """ Does Median window flagging

    Flag data based on deviations from a running median
    See documentation for task MednFlag for details
    uv       = UV data object to flag
    target   = Target source name or list of names, blank = all
    err      = Obit error/message stack
    flagTab  = Output Flagging table version
    flagSig  = Flagging level (sigma)
    alpha    = Smoothing parameter
    timeWind = Averaging window (min)
    doCalib  = Apply calibration table
    gainUse  = CL/SN table to apply
    doBand   = If >0.5 apply bandpass cal.
    BPVer    = Bandpass table version
    flagVer  = Input Flagging table version
    """
    ################################################################
    medn=ObitTask.ObitTask("MednFlag")
    setname(uv,medn)
    if type(target)==list:
        medn.Sources=target
    else:
        medn.Sources=[target]
    medn.flagTab  = flagTab
    medn.flagSig  = flagSig
    medn.alpha    = alpha
    medn.timeWind = timeWind
    medn.doCalib  = doCalib
    medn.gainUse  = gainUse
    medn.doBand   = doBand
    medn.BPVer    = BPVer
    medn.flagVer  = flagVer
    medn.g
    # end VLAMedianFlag
    
def VLAAutoFlag(uv, target, err, \
                doCalib=0, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
                flagTab=1, VClip=[0.0,0.0], IClip=[0.0,0.0], RMSClip=[0.0,0.0], \
                RMSAvg=0.0, maxBad=0.25 ,timeAvg=1.0, \
                doFD=False, FDmaxAmp=0.0, FDmaxV=0.0, FDwidMW=5, FDmaxRMS=[0.0,0.0], \
                FDmaxRes=6.0, FDmaxResBL=6.0,  FDbaseSel=[0, 0, 0, 0]):
    """ Does Automated flagging

    Flag data based on any of a number of criteria
    See documentation for task AutoFlag for details
    uv         = UV data object to flag
    target     = Target source name or list of names, blank = all
    err        = Obit error/message stack
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    flagTab    = Output Flagging table version
    VClip      = If > 0.0 VPol clipping level
    IClip      = If > 0.0 IPol clipping level
    RMSClip    = Abs and fractional clip levels for
                 Time domain RMS filtering
    RMSAvg     = Max RMS/Avg for time domain RMS filtering
    maxBad     = Maximum fraction of baselines for
                 correlator or antenna to be
                 flagged before all are flagged
    timeAvg    = Flagging interval (min)
    doFD       = do frequency domain editing?
    FDmaxAmp   = Maximum average amplitude
    FDmaxV     = Maximum average VPol amp
    FDwidMW    = Width of the median window
    FDmaxRMS   = Channel RMS limits
    FDmaxRes   = Max. residual flux in sigma
    FDmaxResBL = Max. baseline residual
    FDbaseSel  =  Channels for baseline fit (start, end, iincrement,IF)
    """
    ################################################################
    af=ObitTask.ObitTask("AutoFlag")
    setname(uv,af)
    if type(target)==list:
        af.Sources=target
    else:
        af.Sources=[target]
    af.flagTab    = flagTab
    af.flagVer    = flagVer
    af.doCalib    = doCalib
    af.gainUse    = gainUse
    af.doBand     = doBand
    af.BPVer      = BPVer
    af.VClip      = VClip
    af.IClip      = IClip
    af.RMSClip    = RMSClip
    af.RMSAvg     = RMSAvg
    af.maxBad     = maxBad
    af.timeAvg    = timeAvg 
    af.doFD       = doFD
    af.FDmaxAmp   = FDmaxAmp
    af.FDmaxV     = FDmaxV
    af.FDwidMW    = FDwidMW
    af.FDmaxRMS   = FDmaxRMS
    af.FDmaxRes   = FDmaxRes 
    af.FDmaxResBL = FDmaxResBL
    af.FDbaseSel  = FDbaseSel
    af.g
    # end VLAAutoFlag

def VLACal(uv, target, ACal, err, \
           PCal=None, FQid=1, calFlux=None, \
           calModel=None, calDisk=0, \
           solnVer=1, solInt=10.0, nThreads=1, refAnt=0, ampScalar=False):
    """ Basic Amplitude and phase cal for VLA data

    Amplitude calibration can be based either on a point flux
    density or a calibrator model.
    If neither calFlux nor calModel is given, an attempt is made
    to use the setjy.OPType="CALC" option.
    uv       = UV data object to calibrate
    target   = Target source name or list of names
    ACal     = Amp calibrator
    err      = Obit error/message stack
    PCal     = if given, the phase calibrator name
    FQid     = Frequency Id to process
    calFlux  = ACal point flux density if given
    calModel = Amp. calibration model FITS file
               Has priority over calFlux
    calDisk  = FITS disk for calModel
    solnVer  = output SN table version
    solInt   = solution interval (min)
    nThreads = Number of threads to use
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    """
    ################################################################

    # Run SetJy
    setjy = ObitTask.ObitTask("SetJy")
    setname(uv,setjy)
    setjy.Sources=[ACal]
    if FQid:
        setjy.FreqID=FQid
    if calFlux:
        setjy.ZeroFlux=[calFlux]
    elif calModel==None:
        setjy.OPType="CALC"
    setjy.g
    if PCal:
        setjy.ZeroFlux=[0.0,0.0,0.0,0.0]
        if type(PCal)==list:
            setjy.Sources=PCal
        else:
            setjy.Sources=[PCal]
        setjy.OPType="REJY"
        setjy.g
    # Calib on Amp cal
    calib = ObitTask.ObitTask("Calib")
    setname(uv,calib)
    calib.Sources  = [ACal]
    calib.flagVer  = 1
    calib.ampScalar= ampScalar
    calib.doCalib  = 2
    calib.solMode  ="A&P"
    calib.solnVer  = solnVer
    calib.nThreads = nThreads
    calib.solInt   = solInt
    calib.refAnt   = refAnt
    # Given model?
    if calModel:
        calib.DataType2 = "FITS"
        calib.in2File   = calModel
        calib.in2Disk   = calDisk
        calib.nfield    = 1
        calib.CCVer     = 0
        calib.Cmethod   = "DFT"
        calib.Cmodel    = "COMP"
    calib.g

    # ClCal CL1 + SN1 => Cl2
    clcal=ObitTask.ObitTask("CLCal")
    setname(uv,clcal)
    clcal.solnVer = uv.GetHighVer("AIPS SN")
    clcal.calSour = [ACal]
    clcal.calIn   = uv.GetHighVer("AIPS CL")
    clcal.calOut  = clcal.calIn+1
    clcal.interMode="2PT"
    clcal.FreqID = FQid
    # Accumulate full list of sources to calibrate
    if type(target)==list:
        tarlist=[]
        for tar in target:
            tarlist.append(tar)
    else:
        tarlist=[target]
    tarlist.append(ACal)

    # Calib on phase reference if given
    if PCal:
        if type(PCal)==list:
            calib.Sources=PCal
        else:
            calib.Sources=[PCal]
        calib.in2File   = "    "
        calib.nfield    = 0
        calib.flagVer   = 1
        calib.ampScalar = True
        calib.doCalib   = 2
        calib.solMode   = "A&P"
        calib.Cmethod   = "DFT"
        calib.Cmodel    = " "
        calib.g
        
        # GetJy to set flux density scale
        getjy = ObitTask.ObitTask("GetJy")
        setname(uv,getjy)
        getjy.calSour=[ACal]
        getjy.Sources=[PCal]
        getjy.FreqID = FQid
        getjy.g

        # Set up for CLCal
        tarlist.append(PCal)
        clcal.calSour = [PCal]
        
    clcal.Sources = tarlist
    print "Apply calibration for",tarlist
    clcal.g
    
    # end PCal calibration
    # end VLACal

def VLABPCal(uv, BPCal, err, newBPVer=1,
             doCalib=2, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
             solInt=0.0, refAnt=0, ampScalar=False, specIndex=0.0):
    """ Bandbass calibration

    Do bandbass calibration, write BP table
    uv       = UV data object to calibrate
    BPCal    = Bandpass calibrator, name or list of names
    err      = Obit error/message stack
    newBPVer = output BP table
    doCalib  = Apply calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    doBand   = If >0.5 apply previous bandpass cal.
    BPVer    = previous Bandpass table (BP) version
    flagVer  = Input Flagging table version
    solInt   = solution interval (min), 0=> scan average
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    specIndex= spectral index of calibrator (steep=-0.70)
    """
    ################################################################
    bpass = AIPSTask.AIPSTask("bpass")
    setname(uv,bpass)
    if type(BPCal)==list:
        bpass.calsour[1:] = BPCal
    else:
        bpass.calsour[1:] = [BPCal]
    bpass.bpver   = BPVer
    bpass.outver  = newBPVer
    bpass.docalib = doCalib
    bpass.gainuse = gainUse
    bpass.doband  = doBand
    bpass.flagver = flagVer
    bpass.solint  = solInt
    bpass.specindx= specIndex
    bpass.refant  = refAnt
    if ampScalar:
        bpass.bpassprm[8] = 1.0
    bpass.g
    # End VLABPCal

def VLASplit(uv, target, err, FQid=1, outClass="      "):
    """ Write calibrated data

    uv       = UV data object to clear
    target   = Target source name source name or list of names
    err      = Obit error/message stack
    FQid     = Frequency Id to process
    """
    ################################################################
    split=ObitTask.ObitTask("Split")
    setname(uv,split)
    if type(target)==list:
        split.Sources=target
    else:
        split.Sources=[target]
    split.doCalib = 2
    split.gainUse = 0
    split.flagVer = 1
    split.FreqID = FQid
    split.outClass = outClass
    split.outDisk  = split.inDisk
    split.g
    # end VLAsplit

def VLASetImager (uv, target, outIclass="", nThreads=1):
    """ Setup to run Imager

    return Imager task interface object
    uv       = UV data object to image
    target   = Target source name or list of names
    outIclass= output class
    FQid     = Frequency Id to process
    nThreads = Number of threads to use
    """
    ################################################################
    img = ObitTask.ObitTask("Imager")
    setname(uv,img)
    img.outDisk  = img.inDisk
    img.out2Disk = img.inDisk
    if type(target)==list:
        img.Sources=target
    else:
        img.Sources=[target]
    img.outClass=outIclass
    img.doCalib = 2
    img.doBand = 1
    img.UVTaper=[0.0, 0.0, 0.0]
    img.UVRange=[0.0,0.0]
    img.FOV=0.05
    img.autoWindow=True
    img.Niter=5000
    img.Gain=0.10
    img.maxPSCLoop=3
    img.minFluxPSC=0.5
    img.solPInt=10.0/60.
    img.solPType="L1"
    img.maxASCLoop=1
    img.minFluxPSC=1.5
    img.solPInt=10.0/60.0
    img.minSNR=3.0
    img.avgPol=True
    img.avgIF=True
    img.nThreads = nThreads
    return img
# end VLASetImager
