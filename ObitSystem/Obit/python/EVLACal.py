""" Functions for calibrating and editing EVLA data
"""
import UV, UVDesc, Image, ImageDesc, FArray, ObitTask, AIPSTask, AIPSDir, OErr, History
import InfoList
from AIPS import AIPS
from FITS import FITS
from AIPSDir import AIPSdisks, nAIPS
import re

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
    
def imstat (inImage, err, blc=[1,1,1,1,1], trc=[0,0,0,0,0], logfile=None):
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
    logfile  = file to write results to
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
    mess =  "Image statistics:  Region Mean %g, RMSHist %g RMS %g" % (Mean, RMS, RawRMS)
    printMess(mess, logfile)
    mess =  "  Max %g @ pixel %s" % (Max, str(MaxPos))
    printMess(mess, logfile)
    mess = "  Min %g @ pixel %s" % (Min,  str(MinPos))
    printMess(mess, logfile)
    if (head["beamMaj"]>0.0) :
        mess = "  Integrated Flux density %g, beam area = %7.1f pixels" % (Flux, beamarea)
        printMess(mess, logfile)
   
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

def AllDest (err, disk=None, Atype="  ", Aname="            ", Aclass="      ", Aseq=0):
    """ Delete AIPS files matching a pattern

    Strings use AIPS wild cards:
        blank => any
        '?'   => one of any character
        "*"   => arbitrary string
    disk      = AIPS disk number, 0=>all
    Atype     = AIPS entry type, 'MA' or 'UV'; '  => all
    Aname     = desired AIPS name, using AIPS wildcards, None -> don't check
    Aclass    = desired AIPS class, using AIPS wildcards, None -> don't check
    Aseq      = desired AIPS sequence, 0=> any
    """
    ################################################################
    global Adisk
    if disk==None:
        disk = Adisk
    else:
        if disk>0:
            Adisk = disk
    # NO Confirm
    #prompt = "Do you really want to delete all AIPS "+Atype+" files on disk(s) "\
    #         +str(disk)+"\nwith names matching "+Aname+"."+Aclass+"."+str(Aseq)+ \
    #         " y/n "
    #ans = raw_input(prompt)
    #if ans.startswith('y'):
    AIPSDir.PAllDest (disk, err, Atype=Atype, Aname=Aname, Aclass=Aclass,
                      Aseq=Aseq)
    #else:
    #    print "Not confirmed"
    OErr.printErrMsg(err, "Error with destroying AIPS enreies")
    # end AllDest

def printMess (message, logfile=''):
    """ Print message, optionally in logfile

        
        message = message to print
        logfile = logfile for message
    """
    print message
    if len(logfile)>0:
        f = file(logfile,'a')
        f.writelines(message+"\n")
        f.close()
    # end printMess

def dhms2day(st):
    """ convert a time string in d/hh:mm:ss.s to days

    Returns time in days
    st        time string as "d/hh:mm:ss.s"
    """
    ################################################################
    stt = st
    if st.__contains__("/"):
        pp=stt.split("/")
        day = int(pp[0])
        stt = pp[1]
    else:
        day = 0
    pp=stt.split(":")
    if len(pp)>0:
        hour = int(pp[0])
    else:
        hour = 0
    if len(pp)>1:
        min = int(pp[1])
    else:
        min = 0
    if len(pp)>2:
        ssec = float(pp[2])
    else:
        ssec = 0.0
    tim = day + hour/24.0 + min/1440.0 + ssec/86400.0
    return tim
    # end dhms2day

def day2dhms(tim):
    """ convert a time in days to a string as d/hh:mm:ss.s

    Returns time as string:  "d/hh:mm:ss.s"
    tim       time in days
    """
    ################################################################
    day=int(tim)
    ttim = 24.0*(tim - day)
    thour = min (int(ttim), 23)
    ttim = 60.0*(ttim - thour)
    tmin = min (int(ttim), 59)
    ssec = 60.0*(ttim - tmin)
    return str(day)+"/"+str(thour).zfill(2)+":"+str(tmin).zfill(2)+\
           ":"+str(ssec)
    # end day2dhms

def EVLAClearCal(uv, err, doGain=True, doBP=False, doFlag=False):
    """ Clear previous calibration

    Delete all SN tables, all CL but CL 1
    uv       = UV data object to clear
    err      = Obit error/message stack
    doGain   = If True, delete SN and CL tables
    doBP     = If True, delete BP tables
    doFlag   = If True, delete FG tables except FG=1
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
        ver = uv.GetHighVer("AIPS FG")
        while (ver>1):
            uv.ZapTable ('AIPS FG', ver, err)
            ver = ver-1
    OErr.printErrMsg(err, "EVLAClearCal: Error reseting calibration")
    # end EVLAClearCal

def EVLACopyFG(uv, err, logfile='', check=False, debug = False):
    """ Copy AIPS FG table from 1 to 2

    Returns task error code, 0=OK, else failed
    uv       = UV data object to copy
    err      = Obit error/message stack
    logfile  = logfile for messages
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    taco = ObitTask.ObitTask("TabCopy")
    setname(uv, taco)
    taco.outName  = taco.inName
    taco.outClass = taco.inClass
    taco.outDisk  = taco.inDisk
    taco.outSeq   = taco.inSeq
    taco.inTab    = "AIPS FG"
    taco.inVer    = 1
    taco.outVer   = 2
    taco.logFile = logfile
    if debug:
        taco.debug = debug
        taco.i
    # Trap failure
    try:
        if not check:
            taco.g
    except:
        mess = "Copy of FG table Failed retCode="+str(taco.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLACopyFG

def EVLAUVLoad(filename, inDisk, Aname, Aclass, Adisk, Aseq, err, logfile=''):
    """ Read FITS uvtab file into AIPS

    Returns task error code, 0=OK, else failed
    Read a UVTAB FITS UV data file and write an AIPS data set
    filename   = name of FITS file
    inDisk     = FITS directory number
    Aname      = AIPS name of file
    Aclass     = AIPS class of file
    Aseq       = AIPS sequence number of file, 0=> create new
    Adisk      = FITS directory number
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
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
    mess = "Creating AIPS UV file "+Aname+"."+Aclass+"."+str(Aseq)+" on disk "+str(Adisk)
    printMess(mess, logfile)
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
    # end EVLAUVLoad

def EVLAUVLoadT(filename, disk, Aname, Aclass, Adisk, Aseq, err, logfile="  ", \
                    check=False, debug = False,  Compress=False):
    """ Read FITS file into AIPS

    Read input uvtab FITS file, write AIPS
    Returns Obit uv object, None on failure
    Filename = input FITS uvtab format file
    disk     = input FITS file disk number
    Aname    = output AIPS file name
    Aclass   = output AIPS file class
    Adisk    = output AIPS file disk
    Aseq     = output AIPS file sequence
    err      = Obit error/message stack
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
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
    uvc.logFile  = logfile
    if debug:
        uvc.i
        uvc.debug = debug
    # Trap failure
    try:
        if not check:
            uvc.g
    except:
        mess = "UVData load Failed "
        printMess(mess, logfile)
    else:
        pass
    
    # Get output
    if uvc.retCode==0:
        outuv = UV.newPAUV("UVdata", Aname, Aclass, Adisk, Aseq, True, err)
    else:
        outUV = None
    return outuv
    # end EVLAUVLoadT

def EVLAUVLoadArch(dataroot, Aname, Aclass, Adisk, Aseq, err, logfile="  ", \
                       selBand="", selChan=0, dropZero=True, calInt=0.5,  Compress=False, \
                       check=False, debug = False):
    """ Read EVLA archive into AIPS

    Read EVLA archive file, write AIPS
    Returns Obit uv object, None on failure
    dataroot = root of archive directory structure
    Aname    = output AIPS file name
    Aclass   = output AIPS file class
    Adisk    = output AIPS file disk
    Aseq     = output AIPS file sequence
    err      = Obit error/message stack
    selBand  = Selected band, def first
    selChan  = Selected number of channels, def first
    dropZero = If true drop records with all zeroes
    calInt   = CL table interval
    Compress = Write AIPS data in compressed form?
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    #
    outuv        = None
    bdf = ObitTask.ObitTask("BDFIn")
    bdf.DataRoot = dataroot
    bdf.DataType = "AIPS"
    bdf.outName  = Aname
    bdf.outClass = Aclass
    bdf.outSeq   = Aseq
    bdf.outDisk  = Adisk
    bdf.selBand  = selBand
    bdf.selChan  = selChan
    bdf.dropZero = dropZero
    bdf.calInt   = calInt
    bdf.Compress = Compress
    bdf.logFile  = logfile
    if debug:
        bdf.i
        bdf.debug = debug
    # Trap failure
    try:
        if not check:
            bdf.g
    except:
        mess = "UVData load Failed "
        printMess(mess, logfile)
    else:
        pass
    
    # Get output
    if bdf.retCode==0:
        outuv = UV.newPAUV("UVdata", Aname, Aclass, Adisk, Aseq, True, err)
    else:
        outUV = None
    
    # Dummy entry to ensure FG table 1
    UV.PFlag (outuv, err, timeRange=[1.0e20,1.0e21], Ants=[999,0], Reason="Dummy flag")

    return outuv
    # end EVLAUVLoadArch

def EVLAImFITS(inImage, filename, outDisk, err, fract=None, quant=None, \
          exclude=["AIPS HI","AIPS PL","AIPS SL"], include=["AIPS CC"],
          headHi=False):
    """ Write AIPS image as FITS

    Write a Image data set as a FITAB format file
    History also copied
    inImage    = Image data to copy
    filename   = name of FITS file, any whitespace characters replaced with underscore
    outDisk    = FITS directory number
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
    # Deblank filename
    fn = re.sub('\s','_',filename)
    # Set output
    outImage = Image.newPFImage("FITS Image DATA", fn, outDisk, False, err)
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
    outHistory.WriteRec(-1,"imtab   / FITS file "+fn+", disk "+str(outDisk),err)
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
    # end EVLAImFITS

def EVLAUVFITS(inUV, filename, outDisk, err, compress=False, \
              exclude=["AIPS HI", "AIPS AN", "AIPS FQ", "AIPS SL", "AIPS PL"], \
                  include=[], headHi=False):
    """ Write UV data as FITS file
    
    Write a UV data set as a FITAB format file
    History written to header
    inUV       = UV data to copy
    filename   = name of FITS file, any whitespace characters replaced with underscore 
    inDisk     = FITS directory number
    err        = Python Obit Error/message stack
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy (FQ, AN always done )
                 Exclude has presidence over include
    headHi     = if True move history to header, else leave in History table
    returns FITS UV data object
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Deblank filename
    fn = re.sub('\s','_',filename)
    # Set output
    outUV = UV.newPFUV("FITS UV DATA", fn, outDisk, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating FITS data")
    #Compressed?
    if compress:
        inInfo = UV.PGetList(outUV)    # 
        dim = [1,1,1,1,1]
        InfoList.PAlwaysPutBoolean (inInfo, "Compress", dim, [True])        
    # Copy
    UV.PCopy (inUV, outUV, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying UV data to FITS")
    # History
    inHistory  = History.History("inhistory",  outUV.List, err)
    outHistory = History.History("outhistory", outUV.List, err)
    # Add history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit uvtab",err)
    outHistory.WriteRec(-1,"uvtab   / FITS file "+fn+" disk "+str(outDisk),err)
    outHistory.Close(err)
    # History in header?
    if headHi:
        History.PCopy2Header (inHistory, outHistory, err)
        OErr.printErrMsg(err, "Error with history")
        # zap table
        outHistory.Zap(err)
    # Copy Tables
    UV.PCopyTables (inUV, outUV, exclude, include, err)
    return outUV  # return new object
    # end uvtab

def EVLAMedianFlag(uv, target, err, \
                       flagTab=2, flagSig=10.0, alpha=0.5, timeWind=2.0, \
                       doCalib=0, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
                       avgTime=0, avgFreq=0, chAvg=1, \
                       check=False, debug = False, \
                       nThreads=1, noScrat=[], logfile = ""):
    """ Does Median window flagging

    Flag data based on deviations from a running median
    See documentation for task MednFlag for details
    Returns task error code, 0=OK, else failed
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
    avgTime  = preaveraging time (min)
    avgFreq  = 1=>avg chAvg chans, 2=>avg all chan, 3=> avg chan and IFs
    chAvg    = number of channels to average
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    nThreads = Number of threads to use
    noScrat  = list of disks to avoid for scratch files
    logfile  = Log file for task
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
    medn.avgTime  = avgTime
    medn.avgFreq  = avgFreq
    medn.chAvg    = chAvg
    medn.flagVer  = flagVer
    medn.nThreads = nThreads
    medn.logFile  = logfile
    medn.noScrat  = noScrat
    medn.debug    = debug
    if debug:
        medn.i
    # Trap failure
    try:
        if not check:
            medn.g
    except:
        mess = "Median flagging Failed retCode="+str(medn.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAMedianFlag
    
def EVLAQuack(uv, err, \
                  Stokes = " ", BIF=1, EIF=0, Sources=["  "], FreqID=0, \
                  subA=0, timeRange=[0.,0.], Antennas=[0], flagVer=2, \
                  check=False, debug = False, \
                  begDrop=0.0, endDrop=0.0, Reason="Quack", logfile = ""):
    """ Flags beginning and end of each scan

    Trim start and end of each selected scan,
    nothing done if begDrop=endDrop=0.0
    See documentation for task Quack for details
    Returns task error code, 0=OK, else failed
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
    logfile  = Log file for task
    """
    ################################################################
    # Anything to do?
    if (begDrop<=0) and (endDrop<=0):
        return 0
    
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
    quack.logFile   = logfile
    if debug:
        quack.i
        quack.debug = debug
    # Trap failure
    try:
        if not check:
            quack.g
    except:
        mess = "Quack Failed retCode= "+str(quack.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAQuack
    
def EVLAAutoFlag(uv, target, err, \
                     doCalib=0, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
                     flagTab=2, VClip=[0.0,0.0], IClip=[0.0,0.0], RMSClip=[0.0,0.0], \
                     RMSAvg=0.0, maxBad=0.25 ,timeAvg=1.0, \
                     doFD=False, FDmaxAmp=0.0, FDmaxV=0.0, FDwidMW=5, FDmaxRMS=[0.0,0.0], \
                     FDmaxRes=6.0, FDmaxResBL=6.0,  FDbaseSel=[0, 0, 0, 0], \
                     check=False, debug = False, logfile = ""):
    """ Does Automated flagging

    Flag data based on any of a number of criteria
    See documentation for task AutoFlag for details
    Returns task error code, 0=OK, else failed
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
    FDbaseSel  =  Channels for baseline fit (start, end, increment, IF)
    check      = Only check script, don't execute tasks
    debug      = Run tasks debug, show input
    logfile    = Log file for task
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
    af.logFile    = logfile
    if debug:
        af.i
        af.debug = debug
    # Trap failure
    try:
        if not check:
            af.g
    except:
        mess = "AutoFlag Failed retCode="+str(af.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAAutoFlag

def EVLACal(uv, target, ACal, err, \
                PCal=None, FQid=0, calFlux=None, \
                doCalib=-1, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
                calModel=None, calDisk=0, \
                solnVer=1, solInt=10.0/60.0, solSmo=0.0, nThreads=1, refAnt=0, ampScalar=False, \
                check=False, debug = False, \
                noScrat=[], logfile = ""):
    """ Basic Amplitude and phase cal for EVLA data

    Amplitude calibration can be based either on a point flux
    density or a calibrator model.
    If neither calFlux nor calModel is given, an attempt is made
    to use the setjy.OPType="CALC" option.
    Returns task error code, 0=OK, else failed
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
    solnVer  = output SN table version (+1 if smooth)
    solInt   = solution interval (min)
    solSmo   = if solSmo<solInt smooth solutions to solSmo
    nThreads = Number of threads to use
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    noScrat  = list of disks to avoid for scratch files
    logfile  = Log file for tasks
     """
    ################################################################

    # Run SetJy
    setjy = ObitTask.ObitTask("SetJy")
    setjy.logFile  = logfile
    setname(uv,setjy)
    setjy.Sources=[ACal]
    if FQid:
        setjy.FreqID=FQid
    if calFlux:
        setjy.ZeroFlux=[calFlux]
    else:
        setjy.OPType="CALC"
        setjy.ZeroFlux=[1.0,0.0,0.0,0.0]
    if debug:
        setjy.i
        setjy.debug = debug
    # Trap failure
    try:
        if not check:
            setjy.g
    except:
        mess = "SetJy Failed retCode="+str(setjy.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    if PCal:
        setjy.ZeroFlux=[1.0,0.0,0.0,0.0]
        if type(PCal)==list:
            setjy.Sources=PCal
        else:
            setjy.Sources=[PCal]
        setjy.OPType="REJY"
        #setjy.debug = True # DEBUG
        if setjy.Sources[0]!=ACal:
            if debug:
                setjy.i
                setjy.debug = debug
            # Trap failure
            try:
                if not check:
                    setjy.g
            except:
                mess = "SetJy Failed retCode="+str(setjy.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
    # Calib on Amp cal if not in PCal
    calib = ObitTask.ObitTask("Calib")
    calib.logFile  = logfile
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
        if debug:
            calib.i
            calib.debug = debug
        # Trap failure
        try:
            if not check:
                calib.g
        except:
            mess = "Calib Failed retCode= "+str(calib.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
    # Run if Amp cal if not in PCal
    if not ACal in PCal:
        # ClCal CL1 + SN1 => Cl2
        clcal=ObitTask.ObitTask("CLCal")
        clcal.logFile  = logfile
        setname(uv,clcal)
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
            if debug:
                calib.i
                calib.debug = debug
            # Trap failure
            try:
                if not check:
                    calib.g
            except:
                mess = "Calib Failed retCode= "+str(calib.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
            
            # Smoothing?
            if solSmo>solInt:
                solnVer2 = solnVer+1
                snsmo=ObitTask.ObitTask("SNSmo")
                snsmo.logFile  = logfile
                setname(uv,snsmo)
                snsmo.solnIn  = solnVer
                snsmo.solnOut = solnVer2
                snsmo.smoType = "BOTH"
                snsmo.smoParm = [solSmo*60., solSmo*60.]
                snsmo.smoType = "MWF"
                snsmo.refAnt  = refAnt
                if debug:
                    snsmo.i
                    snsmo.debug = debug
                # Trap failure
                try:
                    if not check:
                        snsmo.g
                except:
                    mess = "SNSmo Failed retCode="+str(snsmo.retCode)
                    printMess(mess, logfile)
                    return 1
                else:
                    pass
            else:
                solnVer2 = solnVer
        # GetJy to set flux density scale if ACal not in PCal
        if not ACal in PCal:
            getjy = ObitTask.ObitTask("GetJy")
            getjy.logFile  = logfile
            setname(uv,getjy)
            getjy.calSour=[ACal]
            getjy.solnVer = solnVer2
            if type(PCal)==list:
                getjy.Sources=PCal
            else:
                getjy.Sources=[PCal]
                getjy.FreqID = FQid
            #getjy.debug = True # DEBUG
            if debug:
                getjy.i
                getjy.debug = debug
            # Trap failure
            try:
                if not check:
                    getjy.g
            except:
                mess = "GetJy Failed retCode="+str(getjy.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
        
        # Set up for CLCal - only use phase calibrators
        clcal.solnVer = uv.GetHighVer("AIPS SN")
        if type(PCal)==list:
            clcal.calSour=PCal
        else:
            clcal.calSour=[PCal]
        
    if type(target)==list:
        clcal.Sources=target
    else:
        clcal.Sources=[target]
    mess = "Apply calibration for "+str(target)
    printMess(mess, logfile)
    #clcal.debug=True
    if debug:
        clcal.i
        clcal.debug = debug
    # Trap failure
    try:
        if not check:
            clcal.g
    except:
        mess = "clcal Failed retCode="+str(clcal.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLACal

def EVLABPCal(uv, BPCal, err, newBPVer=1, timerange=[0.,0.], \
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
                  check=False, debug = False, \
                  nThreads=1, noScrat=[], logfile = ""):
    """ Bandbass calibration

    Do bandbass calibration, write BP table
    Returns task error code, 0=OK, else failed
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
    BChan1   = Low freq. channel,  initial cal
    EChan1   = Highest freq channel, initial cal
    BChan2   = Low freq. channel for BP cal
    EChan2   = Highest freq channel for BP cal
    ChWid2   = Number of channels in running mean BP soln, 
    solInt1  = first solution interval (min), 0=> scan average
    solInt2  = second solution interval (min)
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    specIndex= spectral index of calibrator (steep=-0.70)
    doAuto   = Use autocorrelation spectra? Else, crosscorrelation
    doPol    = Apply polarization cal?
    avgPol   = Avg. poln. in solutions?
    avgIF    = Avg. IFs. in solutions?
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    nThreads = Number of threads to use
    noScrat  = list of disks to avoid for scratch files
    logfile  = Log file for task
    """
    ################################################################
    bpass = ObitTask.ObitTask("BPass")
    bpass.logFile = logfile
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
    
    if debug:
        bpass.i
        bpass.debug = debug
   # Trap failure
    try:
        if not check:
            bpass.g
    except:
        mess = "BPass Failed retCode="+str(bpass.retCode)
        return 1
        printMess(mess, logfile)
    else:
        pass
    return 0
    # End EVLABPCal

def EVLASplit(uv, target, err, FQid=1, outClass="      ", logfile = "", \
                  check=False, debug = False):
    """ Write calibrated data

    Returns task error code, 0=OK, else failed
    uv       = UV data object to clear
    target   = Target source name source name or list of names
    err      = Obit error/message stack
    FQid     = Frequency Id to process
    logfile  = Log file for task
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    split=ObitTask.ObitTask("Split")
    split.logFile = logfile
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
    if debug:
        split.i
        split.debug = debug
    # Trap failure
    try:
        if not check:
            split.g
    except:
        mess = "??? Failed retCode="+str(split.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAsplit

def EVLASetImager (uv, target, outIclass="", nThreads=1, noScrat=[], logfile = ""):
    """ Setup to run Imager

    return MFImage task interface object
    uv       = UV data object to image
    target   = Target source name or list of names
    outIclass= output class
    FQid     = Frequency Id to process
    nThreads = Number of threads to use
    noScrat  = list of disks to avoid for scratch files
    logfile  = Log file for task
    """
    ################################################################
    img = ObitTask.ObitTask("MFImage")
    img.logFile = logfile
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
    img.BLFact     = 1.01
    img.BLchAvg    = True
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
# end EVLASetImager

def EVLAPolCal(uv, InsCals, RLCal, RLPhase, err, RM=0.0, \
                   doCalib=2, gainUse=0, flagVer=-1, \
                   soltype="APPR", fixPoln=False, avgIF=False, \
                   solInt=0.0, refAnt=0, \
                   pmodel=[0.0,0.0,0.0,0.0,0.0,0.0,0.0], \
                   FOV=0.05, niter = 100, \
                   check=False, debug = False, \
                   nThreads=1, noScrat=[], logfile = ""):
    """ Polarization calibration, both instrumental and orientation

    Do Instrumental and R-L calibration
    Instrumental cal uses PCAL, R-L cal is done by imaging each IF in Q and U
    and summing the CLEAN components.
    Returns task error code, 0=OK, else failed
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
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    nThreads = Number of threads to use in imaging
    noScrat  = list of disks to avoid for scratch files
    logfile  = Log file for task
    """
    ################################################################
    # Instrumental calibrtation
    if InsCals!=None:
        pcal = AIPSTask.AIPSTask("pcal")
        pcal.logFile = logfile
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
        if debug:
            pcal.i
            pcal.debug = debug
        # Trap failure
        try:
            if not check:
                pcal.g
        except:
            mess = "PCAL Failed retCode="+str(pcal.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
        # end instrumental poln cal
    
    # R-L phase cal
    if RLCal!=None:
        img = ObitTask.ObitTask("Imager")
        img.logFile    = logfile
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
            if debug:
                img.i
                img.debug = debug
            # Trap failure
            try:
                if not check:
                    img.g
            except:
                mess = "Imager Failed retCode="+str(img.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
            
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
        mess = " IF     IFlux    IRMS    QFlux   QRMS    UFlux  URMS  R-L Corr"
        printMess(mess, logfile)
        for i in range (0,len(IFlux)):
            # REALLY NEED RM Correction!!!!!
            cor = RLPhase - 57.296 * math.atan2(UFlux[i],QFlux[i])
            RLCor.append(cor)
            mess = "%3d  %8.3f %8.3f %7.3f %7.3f %7.3f %7.3f %7.3f "%\
                (i+1, IFlux[i], IRMS[i], QFlux[i], QRMS[i], UFlux[i], URMS[i], cor)
            printMess(mess, logfile)
        # Copy highest CL table
        hiCL = uv.GetHighVer("AIPS CL")

        # Apply R-L phase corrections
        clcor = AIPSTask.AIPSTask("clcor")
        clcor.logFile  = logfile
        setname(uv,clcor)
        clcor.opcode   = "POLR"
        clcor.gainver  = hiCL
        clcor.gainuse  = hiCL+1
        clcor.clcorprm[1:] = RLCor
        if debug:
            clcor.i
            clcor.debug = debug
        # Trap failure
        try:
            if not check:
                clcor.g
        except:
            mess = "CLCOR Failed retCode="+str(clcor.retCode)
            printMess(mess, logfile)
        else:
            pass
        return 1
        # end R-L Cal
    return 0
    # End EVLAPolCal

