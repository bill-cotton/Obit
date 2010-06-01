""" Functions for calibrating and editing VLA data
"""
import UV, UVDesc, Image, ImageDesc, FArray, ObitTask, AIPSTask, OErr, History
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
    
def imstat (inImage, err, blc=[1,1,1,1,1], trc=[0,0,0,0,0]):
    """ Get statistics in a specified region of an image plane

    Returns dictionary with statistics of selected region with entries:
        Mean    = Mean value
        RMSHist = RMS value from a histogram analysis
        RMS     = Simple RMS value
        Max     = maximum value
        MaxPos  = pixel of maximum value
        Min     = minimum value
        MinPos  = pixel of minimum value
        Flux    = Flux density if CLEAN beam given, else -1
        BeamArea= CLEAN Beam area in pixels
    inImage   = Python Image object, created with getname, getFITS
    err      = Obit error/message stack
    blc      = bottom left corner pixel (1-rel)
    trc      = top right corner pixel (1-rel)
    """
    ################################################################
    # Read plane
    p    = Image.PReadPlane(inImage,err,blc=blc,trc=trc)
    OErr.printErrMsg(err, "Error with input image")
    head = inImage.Desc.Dict  # Header

    # Get statistics
    Mean = p.Mean
    RMS  = p.RMS
    RawRMS  = p.RawRMS
    MaxPos=[0,0]
    Max = FArray.PMax(p, MaxPos)
    MaxPos[0] = MaxPos[0]+blc[0]
    MaxPos[1] = MaxPos[1]+blc[1]
    MinPos=[0,0]
    Min = FArray.PMin(p, MinPos)
    MinPos[0] = MinPos[0]+blc[0]
    MinPos[1] = MinPos[1]+blc[1]
    # Integrated flux density
    Flux = -1.0
    beamarea = 1.0
    if (head["beamMaj"]>0.0) :
        beamarea = 1.1331*(head["beamMaj"]/abs(head["cdelt"][0])) * \
                   (head["beamMin"]/abs(head["cdelt"][1]))
        Flux = p.Sum/beamarea
    print "Region Mean %g, RMSHist %g RMS %g" % (Mean, RMS, RawRMS)
    print "  Max %g @ pixel " % Max, MaxPos
    print "  Min %g @ pixel " % Min, MinPos
    if (head["beamMaj"]>0.0) :
        print "  Integrated Flux density %g, beam area = %7.1f pixels" % (Flux, beamarea)
   
    # Reset BLC, TRC
    blc = [1,1,1,1,1]
    trc = [0,0,0,0,0]
    Image.POpen(inImage, Image.READONLY, err, blc=blc, trc=trc)
    Image.PClose (inImage, err)
    OErr.printErrMsg(err, "Error with input image")
    
    del p, blc, trc
    return {"Mean":Mean,"RMSHist":RMS,"RMS":RawRMS,"Max":Max, \
            "MaxPos":MaxPos,"Min":Min,"MinPos":MinPos,"Flux":Flux,
            "BeamArea":beamarea}
    # end imstat
   
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

def VLAUVLoad(filename, inDisk, Aname, Aclass, Adisk, Aseq, err):
    """ Read FITS uvtab file into AIPS

    Read a UVTAB FITS UV data file and write an AIPS data set
    filename   = name of FITS file
    inDisk     = FITS directory number
    Aname      = AIPS name of file
    Aclass     = AIPS class of file
    Aseq       = AIPS sequence number of file, 0=> create new
    Adisk      = FITS directory number
    err        = Python Obit Error/message stack
    returns AIPS UV data object
    """
    ################################################################
    #
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Get input
    inUV = UV.newPFUV("FITS UV DATA", filename, inDisk, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with FITS data")
    # Get output, create new if seq=0
    if Aseq<1:
        OErr.printErr(err)   # Print any outstanding messages
        user = OSystem.PGetAIPSuser()
        Aseq=AIPSDir.PHiSeq(Adisk,user,Aname,Aclass,"MA",err)
        # If it already exists, increment seq
        if AIPSDir.PTestCNO(Adisk,user,Aname,Aclass,"MA",Aseq,err)>0:
            Aseq = Aseq+1
        OErr.PClear(err)     # Clear any message/error
    print "Creating AIPS UV file",Aname,".",Aclass,".",Aseq,"on disk",Adisk
    outUV = UV.newPAUV("AIPS UV DATA", Aname, Aclass, Adisk, Aseq, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")
    # Copy
    UV.PCopy (inUV, outUV, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying UV data to AIPS")
    # Copy History
    inHistory  = History.History("inhistory",  inUV.List, err)
    outHistory = History.History("outhistory", outUV.List, err)
    History.PCopyHeader(inHistory, outHistory, err)
    # Add history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit uvlod",err)
    outHistory.WriteRec(-1,"uvlod   / FITS file "+filename+" disk "+str(inDisk),err)
    outHistory.Close(err)
   #
    # Copy Tables
    exclude=["AIPS HI", "AIPS AN", "AIPS FQ", "AIPS SL", "AIPS PL", "History"]
    include=[]
    UV.PCopyTables (inUV, outUV, exclude, include, err)
    return outUV  # return new object
    # end VLAUVLoad

def VLAUVLoadT(filename, disk, Aname, Aclass, Adisk, Aseq, err, \
              Compress=False):
    """ Read FITS file into AIPS

    Read input uvtab FITS file, write AIPS
    Returns Obit uv object
    Filename = input FITS uvtab format file
    disk     = input FITS file disk number
    Aname    = output AIPS file name
    Aclass   = output AIPS file class
    Adisk    = output AIPS file disk
    Aseq     = output AIPS file sequence
    err      = Obit error/message stack
    Compress = Write AIPS data in compressed form?
    """
    ################################################################
    #
    uvc = ObitTask.ObitTask("UVCopy")
    uvc.DataType = "FITS"
    uvc.inFile   = filename
    uvc.inDisk   = disk
    uvc.outDType = "AIPS"
    uvc.outName  = Aname
    uvc.outClass = Aclass
    uvc.outSeq   = Aseq
    uvc.outDisk  = Adisk
    uvc.Compress = Compress
    uvc.g

    # Get output
    outuv = UV.newPAUV("UVdata", Aname, Aclass, Adisk, Aseq, True, err)
    return outuv
    # end VLAUVLoadT

def VLAImFITS(inImage, filename, outDisk, err, fract=None, quant=None, \
          exclude=["AIPS HI","AIPS PL","AIPS SL"], include=["AIPS CC"],
          headHi=False):
    """ Write AIPS image as FITS

    Write a Image data set as a FITAB format file
    History also copied
    inImage    = Image data to copy
    filename   = name of FITS file
    outDisk     = FITS directory number
    err        = Python Obit Error/message stack
    fract      = Fraction of RMS to quantize
    quant      = quantization level in image units, has precedence over fract
                 None or <= 0 => use fract.
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy
    headHi     = if True move history to header, else leave in History table
    """
    ################################################################
    #
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Set output
    outImage = Image.newPFImage("FITS Image DATA", filename, outDisk, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating FITS data")
    # Check for valid pixels
    if inImage.Desc.Dict["maxval"]<=inImage.Desc.Dict["minval"]:
        fract=None; quant=None
    # Copy
    if fract or quant:
        Image.PCopyQuantizeFITS (inImage, outImage, err, fract=fract, quant=quant)
    else:
        Image.PCopy (inImage, outImage, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying Image data to FITS")
    # Copy History
    inHistory  = History.History("inhistory",  inImage.List, err)
    outHistory = History.History("outhistory", outImage.List, err)
    History.PCopy(inHistory, outHistory, err)
    # Add this programs history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit imtab",err)
    if fract:
        outHistory.WriteRec(-1,"imtab   / Quantized at "+str(fract)+" RMS",err)
    outHistory.WriteRec(-1,"imtab   / FITS file "+filename+", disk "+str(outDisk),err)
    outHistory.Close(err)
    # History in header?
    if headHi:
        # Copy back to header
        inHistory  = History.History("inhistory",  outImage.List, err)
        History.PCopy2Header (inHistory, outHistory, err)
        # zap table
        outHistory.Zap(err)
    OErr.printErrMsg(err, "Error with history")
    # Copy Tables
    Image.PCopyTables (inImage, outImage, exclude, include, err)
    # end VLAImFITS

def VLAMedianFlag(uv, target, err, \
                  flagTab=1, flagSig=10.0, alpha=0.5, timeWind=2.0,
                  doCalib=0, gainUse=0, doBand=0, BPVer=0, flagVer=-1,
                  nThreads=1, noScrat=[]):
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
    nThreads = Number of threads to use
    noScrat  = list of disks to avoid for scratch files
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
    medn.nThreads = nThreads
    medn.noScrat  = noScrat
    medn.g
    # end VLAMedianFlag
    
def VLAQuack(uv, err, \
             Stokes = " ", BIF=1, EIF=0, Sources=["  "], FreqID=0, \
             subA=0, timeRange=[0.,0.], Antennas=[0], flagVer=1, \
             begDrop=0.0, endDrop=0.0, Reason="Quack"):
    """ Flags beginning and end of each scan

    Trim start and end of each selected scan,
    nothing done if begDrop=endDrop=0.0
    See documentation for task Quack for details
    uv       = UV data object to flag
    err      = Obit error/message stack
    Stokes   = Limit flagging by Stokes
    BIF      = Limit flagging to BIF-EIF
    EIF      = Limit flagging
    Sources  = Sources selected
    subA     = Subarray number 0=>all
    FreqID   = Freq. ID to flag. -1=>all
    timeRange= Time range to process
    Antennas = List of antennas to include
    flagVer  = Flag table version, 0 => highest
    begDrop  = Time (min) to drop from beginning
    endDrop  = Time (min) to drop from end
    Reason   = Reason (max 24 char.)
    """
    ################################################################
    # Anything to do?
    if (begDrop<=0) and (endDrop<=0):
        return
    
    quack=ObitTask.ObitTask("Quack")
    
    setname(uv, quack)
    quack.Stokes    = Stokes
    quack.BIF       = BIF
    quack.EIF       = EIF
    quack.Sources   = Sources
    quack.subA      = subA
    quack.FreqID    = FreqID
    quack.timeRange = timeRange
    quack.Antennas  = Antennas
    quack.flagVer   = flagVer
    quack.begDrop   = begDrop
    quack.endDrop   = endDrop
    quack.Reason    = Reason
    quack.g
    # end VLAQuack
    
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
           PCal=None, FQid=0, calFlux=None, \
           doCalib=-1, gainUse=0, doBand=0, BPVer=0, flagVer=-1, 
           calModel=None, calDisk=0, \
           solnVer=1, solInt=10.0, nThreads=1, refAnt=0, ampScalar=False,
           noScrat=[]):
    """ Basic Amplitude and phase cal for VLA data

    Amplitude calibration can be based either on a point flux
    density or a calibrator model.
    If neither calFlux nor calModel is given, an attempt is made
    to use the setjy.OPType="CALC" option.
    uv       = UV data object to calibrate
    target   = Target source name or list of names to calibrate
    ACal     = Amp calibrator
    err      = Obit error/message stack
    PCal     = if given, the phase calibrator name
    FQid     = Frequency Id to process, 0=>any
    calFlux  = ACal point flux density if given
    calModel = Amp. calibration model FITS file
               Has priority over calFlux
    calDisk  = FITS disk for calModel
    doCalib  = Apply calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    doBand   = If >0.5 apply previous bandpass cal.
    BPVer    = previous Bandpass table (BP) version
    solnVer  = output SN table version
    solInt   = solution interval (min)
    nThreads = Number of threads to use
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    noScrat  = list of disks to avoid for scratch files
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
    else:
        setjy.ZeroFlux=[1.0,0.0,0.0,0.0]
    setjy.g
    if PCal:
        setjy.ZeroFlux=[1.0,0.0,0.0,0.0]
        if type(PCal)==list:
            setjy.Sources=PCal
        else:
            setjy.Sources=[PCal]
        setjy.OPType="REJY"
        #setjy.debug = True # DEBUG
        if setjy.Sources[0]!=ACal:
            setjy.g
    # Calib on Amp cal if not in PCal
    calib = ObitTask.ObitTask("Calib")
    setname(uv,calib)
    calib.Sources  = [ACal]
    calib.flagVer  = 1
    calib.ampScalar= ampScalar
    calib.doCalib  = doCalib
    calib.gainUse  = gainUse
    calib.doBand   = doBand
    calib.BPVer    = BPVer
    calib.solMode  ="A&P"
    calib.solnVer  = solnVer
    calib.nThreads = nThreads
    calib.solInt   = solInt
    calib.refAnt   = refAnt
    calib.noScrat  = noScrat
    # Given model?
    if calModel:
        calib.DataType2 = "FITS"
        calib.in2File   = calModel
        calib.in2Disk   = calDisk
        calib.nfield    = 1
        calib.CCVer     = 0
        calib.Cmethod   = "DFT"
        calib.Cmodel    = "COMP"
    #calib.prtLv = 5  # DEBUG
    #calib.debug    = True  # DEBUG
    # Run if Amp cal if not in PCal
    if not ACal in PCal:
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

    # Calib on phase reference if given
    if PCal:
        if type(PCal)==list:
            calib.Sources=PCal
        else:
            calib.Sources=[PCal]
        calib.in2File   = "    "
        calib.nfield    = 0
        calib.flagVer   = 1
        calib.ampScalar = ampScalar
        calib.doCalib   = 2
        calib.solMode   = "A&P"
        calib.Cmethod   = "DFT"
        calib.Cmodel    = " "
        calib.g
        
        # GetJy to set flux density scale if ACal not in PCal
        if not ACal in PCal:
            getjy = ObitTask.ObitTask("GetJy")
            setname(uv,getjy)
            getjy.calSour=[ACal]
            if type(PCal)==list:
                getjy.Sources=PCal
            else:
                getjy.Sources=[PCal]
                getjy.FreqID = FQid
            #getjy.debug = True # DEBUG
            getjy.g

        # Set up for CLCal - only use phase calibrators
        if type(PCal)==list:
            clcal.calSour=PCal
        else:
            clcal.calSour=[PCal]
        
    if type(target)==list:
        clcal.Sources=target
    else:
        clcal.Sources=[target]
    print "Apply calibration for",target
    #clcal.debug=True
    clcal.g
    
    # end PCal calibration
    # end VLACal

def VLABPCal(uv, BPCal, err, newBPVer=1, timerange=[0.,0.], \
             doCalib=2, gainUse=0, doBand=0, BPVer=0, flagVer=-1,  \
             CalDataType="    ", CalFile="  ", CalName="  ", CalClass="  ", \
             CalSeq=0, CalDisk=0, CalNfield=0, CalCCVer=0, \
             CalBComp=[0], CalEComp=[0], CalCmethod=" ", CalCmodel=" ", CalFlux=0.0, \
             modelFlux=0.0, modelPos=[0.0], modelParm=[0.0], \
             doCenter1=None, BChan1=1, EChan1=0, \
             BChan2=1, EChan2=0, ChWid2=1, \
             solInt1=0.0, solInt2=0.0, refAnt=0, \
             ampScalar=False, specIndex=0.0, \
             doAuto=False, doPol=False, avgPol=False, avgIF=False, \
             nThreads=1, noScrat=[]):
    """ Bandbass calibration

    Do bandbass calibration, write BP table
    Calibration is done in two passes
    1)  First a wideband phase only calibration using channels
    BChan1 to EChan1 or the central doCenter1 fraction of the band
    using a solution interval of solInt1.  This solution is applied
    to all selected data and used in the second pass.
    2)  Second channels in the range BChan2 to EChan2 averaging blocks
    of ChWid2 are calibrated using solTYpe and solMode for solInt2 and
    the results written as the output BP table.
       The Calibrator model may be given as either and Image with CC table,
    a parameterized model or a point source with the flux density in 
    the SU table.
    See BPass documentation for details
    
    uv       = UV data object to calibrate
    BPCal    = Bandpass calibrator, name or list of names
    err      = Obit error/message stack
    newBPVer = output BP table
    doCalib  = Apply calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    doBand   = If >0.5 apply previous bandpass cal.
    BPVer    = previous Bandpass table (BP) version
    flagVer  = Input Flagging table version
    timerange= timerange in days to use
    CalDataType =  Calibrator model file data type (AIPS,FITS)
    CalFile  = Calibrator model FITS input image if Type=='FITS'
    CalName  = Calibrator model Cleaned map name 
    CalClass = Calibrator model Cleaned map class
    CalSeq   = Calibrator model Cleaned map seq
    CalDisk  = Calibrator model Cleaned map disk
    CalNfield= Calibrator model  No. maps to use for model
    CalCCVer = Calibrator model CC file version
    CalBComp = Calibrator model First CLEAN comp to use, 1/field
    CalEComp = Calibrator model  Last CLEAN comp to use, 0=>all
    CalCmethod= Calibrator model Modeling method, 'DFT','GRID','    '
    CalCmodel= Calibrator model Model type: 'COMP','IMAG'
    CalFlux  = Calibrator model  Lowest CC component used
    modelFlux= Parameterized model flux density (Jy)
    modelPos = Parameterized model Model position offset (asec)
    modelParm= Parameterized model Model parameters (maj, min, pa, type)
    doCenter1= If defined, the center fraction of the bandpass to use first pass
    solInt1  = first solution interval (min), 0=> scan average
    solInt2  = second solution interval (min)
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    specIndex= spectral index of calibrator (steep=-0.70)
    doAuto   = Use autocorrelation spectra? Else, crosscorrelation
    doPol    = Apply polarization cal?
    avgPol   = Avg. poln. in solutions?
    avgIF    = Avg. IFs. in solutions?
    nThreads = Number of threads to use
    noScrat  = list of disks to avoid for scratch files
    """
    ################################################################
    bpass = ObitTask.ObitTask("BPass")
    setname(uv,bpass)
    if type(BPCal)==list:
        bpass.Sources = BPCal
    else:
        bpass.Sources = [BPCal]
    bpass.doBand    = doBand
    bpass.BPVer     = BPVer
    bpass.BPSoln    = newBPVer
    bpass.doCalib   = doCalib
    bpass.gainUse   = gainUse
    bpass.flagVer   = flagVer
    bpass.doPol     = doPol
    bpass.solInt1   = solInt1
    bpass.solInt2   = solInt2
    bpass.Alpha     = specIndex
    bpass.refAnt    = refAnt
    bpass.timeRange = timerange
    bpass.DataType2 = CalDataType
    bpass.in2File   = CalFile
    bpass.in2Name   = CalName
    bpass.in2Class  = CalClass
    bpass.in2Seq    = CalSeq 
    bpass.in2Disk   = CalDisk
    bpass.nfield    = CalNfield
    bpass.CCVer     = CalCCVer
    bpass.BComp     = CalBComp
    bpass.EComp     = CalEComp
    bpass.Cmethod   = CalCmethod
    bpass.Cmodel    = CalCmodel
    bpass.Flux      = CalFlux
    bpass.modelFlux = modelFlux
    bpass.modelPos  = modelPos
    bpass.modelParm = modelParm
    bpass.doAuto    = doAuto
    bpass.avgPol    = avgPol
    bpass.avgIF     = avgIF
    bpass.ampScalar = ampScalar
    bpass.noScrat   = noScrat
    bpass.nThreads  = nThreads

    # Channel selection
    d     = uv.Desc.Dict
    nchan = d["inaxes"][d["jlocf"]]
    # Center fraction requested?
    if doCenter1:
        # Center doCenter1 fraction of channels for first cal
        mchan = int(nchan*doCenter1)
        bpass.BChan1 = max(1, (nchan/2)-(mchan/2))
        bpass.EChan1 = min(nchan, (nchan/2)+(mchan/2))
    else:
        bpass.BChan1 = BChan1
        bpass.EChan1 = EChan1
    bpass.BChan2 = BChan2
    bpass.EChan2 = EChan2
    if bpass.EChan2<=0:
        bpass.EChan2 = nchan
    
    #bpass.i
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

def VLASetImager (uv, target, outIclass="", nThreads=1, noScrat=[]):
    """ Setup to run Imager

    return Imager task interface object
    uv       = UV data object to image
    target   = Target source name or list of names
    outIclass= output class
    FQid     = Frequency Id to process
    nThreads = Number of threads to use
    noScrat  = list of disks to avoid for scratch files
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
    img.outClass   = outIclass
    img.doCalib    = 2
    img.doBand     = 1
    img.UVTaper    = [0.0, 0.0, 0.0]
    img.UVRange    = [0.0,0.0]
    img.FOV        = 0.05
    img.autoWindow = True
    img.Niter      = 5000
    img.Gain       = 0.10
    img.maxPSCLoop = 3
    img.minFluxPSC= 0.5
    img.solPInt   = 10.0/60.
    img.solPType  = "L1"
    img.maxASCLoop= 1
    img.minFluxPSC= 1.5
    img.solAInt   = 1.0
    img.minSNR    = 3.0
    img.avgPol    = True
    img.avgIF     = True
    img.nThreads  = nThreads
    img.noScrat   = noScrat
    return img
# end VLASetImager

def VLAPolCal(uv, InsCals, RLCal, RLPhase, err, RM=0.0, \
              doCalib=2, gainUse=0, flagVer=-1, \
              soltype="APPR", fixPoln=False, avgIF=False, \
              solInt=0.0, refAnt=0, \
              pmodel=[0.0,0.0,0.0,0.0,0.0,0.0,0.0], \
              FOV=0.05, niter = 100, \
              nThreads=1, noScrat=[]):
    """ Polarization calibration, both instrumental and orientation

    Do Instrumental and R-L calibration
    Instrumental cal uses PCAL, R-L cal is done by imaging each IF in Q and U
    and summing the CLEAN components.
    uv       = UV data object to calibrate
    InsCals  = Instrumental poln calibrators, name or list of names
               If None no instrumental cal
    RLCal    = R-L (polarization angle) calibrator,
               If None no R-L cal
    RLPhase  = R-L phase of RLCal (deg) at 1 GHz
    err      = Obit error/message stack
    RM       = Rotation measure of RLCal
    doCalib  = Apply prior calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    flagVer  = Input Flagging table version
    soltype  = solution type
    fixPoln  = if True, don't solve for source polarization in ins. cal
    avgIF    = if True, average IFs in ins. cal.
    solInt   = instrumental solution interval (min), 0=> scan average
    refAnt   = Reference antenna
    pmodel   = Instrumental poln cal source poln model.
               pmodel[0] = I flux density (Jy)
               pmodel[1] = Q flux density (Jy)
               pmodel[2] = U flux density (Jy)
               pmodel[3] = V flux density (Jy)
               pmodel[4] = X offset in sky (arcsec)
               pmodel[5] = Y offset in sky (arcsec)
    FOV      = field of view radius (deg) needed to image RLCal
    niter    = Number  of iterations of CLEAN in R-L cal
    nThreads = Number of threads to use in imaging
    noScrat  = list of disks to avoid for scratch files
    """
    ################################################################
    # Instrumental calibrtation
    if InsCals!=None:
        pcal = AIPSTask.AIPSTask("pcal")
        setname(uv,pcal)
        if type(InsCals)==list:
            pcal.calsour[1:] = InsCals
        else:
            pcal.calsour[1:] = [InsCals]
        pcal.docalib = doCalib
        pcal.gainuse = gainUse
        pcal.flagver = flagVer
        pcal.soltype = soltype
        pcal.solint  = solInt
        pcal.refant  = refAnt
        if fixPoln:
            pcal.bparm[10]=1.0
        if avgIF:
            pcal.cparm[1]=1.0
        pcal.pmodel[1:]  = pmodel
        pcal.i
        pcal.g
        # end instrumental poln cal

    # R-L phase cal
    if RLCal!=None:
        img = ObitTask.ObitTask("Imager")
        setname(uv,img)
        img.doCalib    = doCalib
        img.gainUse    = gainUse
        img.flagVer    = flagVer
        img.doPol      = True
        img.Sources[0] = RLCal
        img.Stokes     = "IQU"
        img.FOV        = FOV
        img.Niter      = niter
        img.autoWindow = True
        img.dispURL    = "None"
        img.Catalog    = "None"
        img.nThreads   = nThreads
        img.noScrat    = noScrat
        # Temporary output files
        if img.DataType=="AIPS":
            img.outName = "TEMP"
            img.outClass= "IPOLCL"
            img.outDisk = img.inDisk
            img.outSeq  = 6666
            img.out2Name = "TEMP"
            img.out2Class= "IPOLCL"
            img.out2Disk = img.inDisk
            img.out2Seq  = 7777
        elif img.DataType=="FITS":
            img.outFile  = "TEMPPOLCAL.fits"
            img.outDisk  = img.inDisk
            img.out2File = "TEMPPOLCAL2.uvtab"
            img.out2Disk = img.inDisk
       # How many IFs?
        h = uv.Desc.Dict
        if h["jlocif"]>=0:
            nif = h["inaxes"][h["jlocif"]]
        else:
            nif = 1

        # Lists of flux densities and RMSes
        IFlux = []
        IRMS  = []
        QFlux = []
        QRMS  = []
        UFlux = []
        URMS  = []
        
        # Loop over IF imaging I,Q, U
        for iif in range (1, nif+1):
            img.BIF = iif
            img.EIF = iif
            #img.dispURL    = "ObitView"  # DEBUG
            #img.debug=True               # DEBUG
            img.g
            
            # Get fluxes from inner quarter of images
            if img.DataType=="AIPS":
                outName = (img.Sources[0].strip()+"TEMP")[0:12]
                outDisk = img.outDisk
                outSeq  = 6666
                # Stokes I
                outClass="IPOLCL"
                x =  Image.newPAImage("I",outName, outClass, outDisk,outSeq,True,err)
                h = x.Desc.Dict
                blc = [h["inaxes"][0]/4,h["inaxes"][1]/4]
                trc = [3*h["inaxes"][0]/4,3*h["inaxes"][1]/4]
                stat = imstat(x, err, blc=blc,trc=trc)
                IFlux.append(stat["Flux"])
                IRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes Q
                outClass="QPOLCL"
                x =  Image.newPAImage("Q",outName, outClass, outDisk,outSeq,True,err)
                stat = imstat(x, err, blc=blc,trc=trc)
                QFlux.append(stat["Flux"])
                QRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes U
                outClass="UPOLCL"
                x =  Image.newPAImage("U",outName, outClass, outDisk,outSeq,True,err)
                stat = imstat(x, err, blc=blc,trc=trc)
                UFlux.append(stat["Flux"])
                URMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Delete UV output
                out2Name = (img.Sources[0].strip()+"TEMP")[0:12]
                out2Class="IPOLCL"
                out2Disk = img.inDisk
                out2Seq  = 7777
                u =  UV.newPAUV("UV",out2Name,out2Class,out2Disk,out2Seq,True,err)
                u.Zap(err)
                del u
            elif img.DataType=="FITS":
                # Stokes I
                outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                x =  Image.newPFImage("I",outFile,img.outDisk,True,err)
                h = x.Desc.Dict
                blc = [h["inaxes"][0]/4,h["inaxes"][1]/4]
                trc = [3*h["inaxes"][0]/4,3*h["inaxes"][1]/4]
                stat = imstat(x, err, blc=blc,trc=trc)
                IFlux.append(stat["Flux"])
                IRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes Q
                outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                stat = imstat(x, err, blc=blc,trc=trc)
                IFlux.append(stat["Flux"])
                IRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes U
                outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                stat = imstat(x, err, blc=blc,trc=trc)
                IFlux.append(stat["Flux"])
                IRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                out2File = img.Sources[0].strip()+"TEMPPOLCAL2.uvtab"
                u =  UV.newPFUV("UV",outFile,img.outDisk,True,err)
                u.Zap(err)
                del u
           # End accumulate statistics by file type
        # End loop over IF

        # Give results, compute R-L correction
        RLCor = []
        import math
        print " IF     IFlux    IRMS    QFlux   QRMS    UFlux  URMS  R-L Corr"
        for i in range (0,len(IFlux)):
            # REALLY NEED RM Correction!!!!!
            cor = RLPhase - 57.296 * math.atan2(UFlux[i],QFlux[i])
            RLCor.append(cor)
            print "%3d  %8.3f %8.3f %7.3f %7.3f %7.3f %7.3f %7.3f "%\
                  (i+1, IFlux[i], IRMS[i], QFlux[i], QRMS[i], UFlux[i], URMS[i], cor)
        
        # Copy highest CL table
        hiCL = uv.GetHighVer("AIPS CL")

        # Apply R-L phase corrections
        clcor = AIPSTask.AIPSTask("clcor")
        setname(uv,clcor)
        clcor.opcode   = "POLR"
        clcor.gainver  = hiCL
        clcor.gainuse  = hiCL+1
        clcor.clcorprm[1:] = RLCor
        clcor.g
        # end R-L Cal
    # End VLAPolCal

def unique (inn):
    """ Removes duplicate entries from an array of strings
    
    Returns an array of strings, also removes null and blank strings
    as well as leading or trailing blanks
    inn  = list of strings with possible redundancies
    """
    # Make local working copy and blank redundant entries
    linn = []
    for item in inn:
        sitem = item.strip()
        if len(sitem)>0:
            linn.append(sitem)
    # Remove duplicates from the end of the list
    n = len(linn)
    for jj in range(0,n):
        j = n-jj-1
        for i in range (0,j):
            if (linn[j]==linn[i]):
                linn[j] = ""
                break;
    # end loops
    # Copy to output string
    outl = []
    for item in linn:
        if len(item)>0:
            outl.append(item)
    return outl
# end unique

