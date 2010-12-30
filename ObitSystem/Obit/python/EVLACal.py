""" Functions for calibrating and editing EVLA data
"""
import UV, UVDesc, Image, ImageDesc, FArray, ObitTask, AIPSTask, AIPSDir, OErr, History
import InfoList, Table, OSystem
from AIPS import AIPS
from FITS import FITS
from AIPSDir import AIPSdisks, nAIPS
import re, pickle, math

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
    
def setoname (inn, out):
    """ Copy file definition from inn to out as out...
    
    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    inn  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask.ObitTask:
        out.DataType = inn.FileType
        out.outDisk   = int(inn.Disk)
        if inn.FileType == 'FITS':
            out.outFile = inn.Fname
        else:   # AIPS
            out.outName  = inn.Aname
            out.outClass = inn.Aclass
            out.outSeq   = int(inn.Aseq)
    else:  # AIPS
        out.outname  = inn.Aname
        out.outclass = inn.Aclass
        out.outseq   = inn.Aseq
        out.outdisk  = inn.Disk
    # end setoname
    
def set2name (in2, out):
    """ Copy file definition from in2 to out as in2...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    in2  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask:
        out.DataType  = in2.FileType
        out.in2Disk   = int(in2.Disk)
        if in2.FileType == 'FITS':
            out.in2File = in2.Fname
        else:   # AIPS
            out.in2Name  = in2.Aname
            out.in2Class = in2.Aclass
            out.in2Seq   = int(in2.Aseq)
    else: # AIPS
        out.in2Name  = in2.Aname
        out.in2Class = in2.Aclass
        out.in2Seq   = in2.Aseq
        out.in2Disk  = in2.Disk
    # end set2name
   
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
    logfile  = file to write results to, if None don't print
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
    if logfile:
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

def SaveObject (pyobj, file, update):
    """ Save python object to a pickle file

    pyobj    = python object to save
    file     = pickle file name
    update   = If True update, otherwise only if file doesn't already exist
    """
    ################################################################
    # Does file exist?, only do this is not or update
    if update or not os.path.isfile(file):
        fd = open(file, "w")
        pickle.dump(pyobj, fd)
        fd.close()
    # end SaveObject
   
def FetchObject (file):
    """ Fetch python object from a pickle file

    returns python object
    file     = pickle file name
    """
    ################################################################
    # unpickle file
    fd = open(file, "r")
    pyobj = pickle.load(fd)
    fd.close()
    return pyobj
    # end FetchObject
   
def EVLAClearCal(uv, err, doGain=True, doBP=False, doFlag=False, check=False):
    """ Clear previous calibration

    Delete all SN tables, all CL but CL 1
    uv       = UV data object to clear
    err      = Obit error/message stack
    doGain   = If True, delete SN and CL tables
    doBP     = If True, delete BP tables
    doFlag   = If True, delete FG tables except FG=1
    check    = Only check script, don't execute tasks
    """
    ################################################################
    # Only checking?
    if check:
        return
    # Open and close image to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)
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
    if not check:
        setname(uv, taco)
    taco.outName  = taco.inName
    taco.outClass = taco.inClass
    taco.outDisk  = taco.inDisk
    taco.outSeq   = taco.inSeq
    taco.inTab    = "AIPS FG"
    taco.inVer    = 1
    taco.outVer   = 2
    taco.taskLog = logfile
    if debug:
        taco.debug = debug
        taco.i
    # Trap failure
    try:
        if not check:
            taco.g
    except Exception, exception:
        print exception
        mess = "Copy of FG table Failed retCode="+str(taco.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLACopyFG

def EVLACopyTable(inObj, outObj, inTab, err, inVer=1, outVer=0,
                  logfile='', check=False, debug = False):
    """ Copy AIPS Table

    Returns task error code, 0=OK, else failed
    inObj    = Input Object (UV or Image)
    outObj   = Output object
    inTab    = Table type, e.g. "AIPS AN"
    err      = Obit error/message stack
    inVer    = intput version
    outVer   = output version
    logfile  = logfile for messages
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    mess =  "Copy "+inTab+" Table "+str(inVer)+" to "+str(outVer)
    printMess(mess, logfile)
    taco = ObitTask.ObitTask("TabCopy")
    if not check:
        setname(inObj, taco)
        setoname(outObj, taco)
    taco.inTab    = inTab
    taco.inVer    = inVer
    taco.outVer   = outVer
    taco.taskLog  = logfile
    if debug:
        taco.debug = debug
        taco.i
    # Trap failure
    try:
        if not check:
            taco.g
    except Exception, exception:
        print exception
        mess = "Copy of "+inTab+" table Failed retCode="+str(taco.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLACopyTable

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
    mess =  "Load FITS uvtab file into AIPS"
    printMess(mess, logfile)
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
    mess =  "Load FITS uvtab file into AIPS"
    printMess(mess, logfile)
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
    uvc.taskLog  = logfile
    if debug:
        uvc.i
        uvc.debug = debug
    # Trap failure
    try:
        if not check:
            uvc.g
    except Exception, exception:
        print exception
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

def EVLAUVLoadArch(dataroot, Aname, Aclass, Adisk, Aseq, err, \
                   selConfig=-1, selBand="", selChan=0, selNIF=0, \
                   dropZero=True, calInt=0.5,  doSwPwr=False, Compress=False, \
                   logfile = "", check=False, debug = False):
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
    selNIF   = Selected number of IFs, def first
    dropZero = If true drop records with all zeroes
    calInt   = CL table interval
    doSwPwr  = Make EVLA Switched power corr?
    Compress = Write AIPS data in compressed form?
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    mess =  "Load archive file into AIPS"
    printMess(mess, logfile)
    #
    outuv        = None
    bdf = ObitTask.ObitTask("BDFIn")
    bdf.DataRoot = dataroot
    bdf.DataType = "AIPS"
    bdf.outName  = Aname
    bdf.outClass = Aclass
    bdf.outSeq   = Aseq
    bdf.outDisk  = Adisk
    bdf.selConfig= selConfig
    bdf.selBand  = selBand
    bdf.selChan  = selChan
    bdf.selIF    = selNIF
    bdf.dropZero = dropZero
    bdf.calInt   = calInt
    bdf.doSwPwr  = doSwPwr
    bdf.Compress = Compress
    bdf.taskLog  = logfile
    if debug:
        bdf.i
        bdf.debug = debug
    # Trap failure
    try:
        if not check:
            bdf.g
    except Exception, exception:
        print exception
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
    if not check:
        UV.PFlag (outuv, err, timeRange=[1.0e20,1.0e21], Ants=[999,0], Reason="Dummy flag")

    return outuv
    # end EVLAUVLoadArch

def EVLAImFITS(inImage, filename, outDisk, err, fract=None, quant=None, \
          exclude=["AIPS HI","AIPS PL","AIPS SL"], include=["AIPS CC"],
          headHi=False, logfile=""):
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
    mess =  "Write Image to FITS "+filename+" on disk "+str(outDisk)
    printMess(mess, logfile)
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
                  include=[], headHi=False, logfile=""):
    """ Write UV data as FITS file
    
    Write a UV data set as a FITAB format file
    History written to header
    inUV       = UV data to copy
    filename   = name of FITS file, any whitespace characters replaced with underscore 
    outDisk    = FITS directory number
    err        = Python Obit Error/message stack
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy (FQ, AN always done )
                 Exclude has presidence over include
    headHi     = if True move history to header, else leave in History table
    returns FITS UV data object
    """
    ################################################################
    mess =  "Write Data to FITS UV data "+filename+" on disk "+str(outDisk)
    printMess(mess, logfile)

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
    # end EVLAUVFITS

def EVLAUVFITSTab(inUV, filename, outDisk, err, \
              exclude=["AIPS HI", "AIPS AN", "AIPS FQ", "AIPS SL", "AIPS PL"], \
                  include=[], logfile=""):
    """ Write Tables on UV data as FITS file
    
    Write Tables from a UV data set (but no data) as a FITAB format file
    History written to header
    inUV       = UV data to copy
    filename   = name of FITS file, any whitespace characters replaced with underscore 
    outDisk    = FITS directory number
    err        = Python Obit Error/message stack
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy (FQ, AN always done )
                 Exclude has presidence over include
    returns FITS UV data object
    """
    ################################################################
    mess =  "Write Tables to FITS UV data "+filename+" on disk "+str(outDisk)
    printMess(mess, logfile)
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
    # Clone
    UV.PClone (inUV, outUV, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error cloning UV data to FITS")
    # Copy back to header
    inHistory  = History.History("inhistory",  outUV.List, err)
    outHistory = History.History("outhistory", outUV.List, err)
    # Add history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit uvTabSave",err)
    outHistory.WriteRec(-1,"uvTabSave / FITS file "+fn+" disk "+str(outDisk),err)
    outHistory.Close(err)
    History.PCopy2Header (inHistory, outHistory, err)
    OErr.printErrMsg(err, "Error with history")
    # zap table
    outHistory.Zap(err)
    # Copy Tables
    UV.PCopyTables (inUV, outUV, exclude, include, err)
    return outUV  # return new object
    # end EVLAUVFITSTab

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
    mess =  "Median Window flagging"
    printMess(mess, logfile)
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
    medn.taskLog  = logfile
    medn.noScrat  = noScrat
    medn.debug    = debug
    #bombaroonie = BombsAwayWithCurtisLemay # DEBUG
    if debug:
        medn.i
    # Trap failure
    try:
        if not check:
            medn.g
    except Exception, exception:
        print exception
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
    
    mess =  "Quack data"
    printMess(mess, logfile)
    quack=ObitTask.ObitTask("Quack")
    
    if not check:
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
    quack.taskLog   = logfile
    if debug:
        quack.i
        quack.debug = debug
    # Trap failure
    try:
        if not check:
            quack.g
    except Exception, exception:
        print exception
        mess = "Quack Failed retCode= "+str(quack.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAQuack
    
def EVLAAutoFlag(uv, target, err, \
                     doCalib=0, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
                     flagTab=2, VClip=[0.0,0.0], IClip=[0.0,0.0], minAmp=0.0, \
                     RMSClip=[0.0,0.0], RMSAvg=0.0, maxBad=0.25 ,timeAvg=1.0, \
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
    minAmp     = Min flux for IClip flagging
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
    mess =  "AutoFlag data"
    printMess(mess, logfile)
    af=ObitTask.ObitTask("AutoFlag")
    if not check:
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
    af.minAmp     = minAmp
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
    af.taskLog    = logfile
    if debug:
        af.i
        af.debug = debug
    # Trap failure
    try:
        if not check:
            af.g
    except Exception, exception:
        print exception
        mess = "AutoFlag Failed retCode="+str(af.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAAutoFlag

def EVLAPACor(uv, err, CLver=0, FreqID=1,\
                  noScrat=[], logfile='', check=False, debug=False):
    """ Make parallactic angle correction

    Updates CL CLver, if only one, a new one (CL 2) is copied
    Returns task error code, 0=OK, else failed
    uv         = UV data object
    err        = Python Obit Error/message stack
    CLver      = Cl version to update, 0=> highest
    FreqID     = Frequency group identifier
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################

    # Which CL?
    if CLver<=0 and not check:
        # Open and close image to sync with disk 
        uv.Open(UV.READONLY, err)
        uv.Close(err)
        CLver = uv.GetHighVer("AIPS CL")
    # If CLver==1, copy to 2 first
    if CLver==1 and not check:
        EVLACopyTable (uv, uv, "AIPS CL", err, inVer=CLver, outVer=CLver+1, \
                       logfile=logfile, check=check, debug=debug)
        if err.isErr:
            print  "Error copying CL Table"
            return 1
        CLver += 1  # Use copy
        
    mess = "Parallactic angle corrections made to CL "+str(CLver)
    printMess(mess, logfile)

    clcor = AIPSTask.AIPSTask("clcor")
    if not check:
        setname(uv,clcor)
    clcor.opcode      = "PANG"
    clcor.gainver     = CLver
    clcor.gainuse     = CLver
    clcor.clcorprm[1] = 1.0
    clcor.freqid      = FreqID
    clcor.logFile     = logfile
    i = 1;
    for d in noScrat:
        clcor.baddisk[i] = d
        i += 1
    if debug:
        clcor.i
    # Trap failure
    try:
        if not check:
            clcor.g
    except Exception, exception:
        print exception
        mess = "CLCOR Failed "
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAPACor

def EVLADelayCal(uv, err, solInt=0.5, smoTime=10.0, calSou=None,  CalModel=None, \
                     timeRange=[0.,0.], FreqID=1, doCalib=-1, gainUse=0, minSNR=5.0, \
                     refAnts=[0], doBand=-1, BPVer=0, flagVer=-1, doTwo=True, \
                     doPlot=False, plotFile="./DelayCal.ps", \
                     nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Group delay calibration

    Determine delay corrections from a list of calibrators
    Solutions optionally smoothed to smoTime
    Apply this SN table to the highest CL table writing a new CL table (Obit/CLCal)
    Returns task error code, 0=OK, else failed
    err        = Python Obit Error/message stack
    calSou     = Source name or list of names to use
    CalModel = python dict with image model dicts keyed on calibrator name
               image object = "Image"
               also optional
               "nfield",    Calibrator model  No. maps to use for model
               "CCver",     Calibrator model CC file version
               "BComp",     Calibrator model First CLEAN comp to use, 1/field
               "EComp"      Calibrator model  Last CLEAN comp to use, 0=>all
               "Cmethod"    Calibrator model Modeling method, 'DFT','GRID','    '
               "CModel"     Calibrator model Model type: 'COMP','IMAG'
               "CalFlux"    Calibrator model  Lowest CC component used
               "modelFlux"  if ["Image"]=None, Parameterized model flux density (Jy)
               "modepPos"   if ["Image"]=None, Parameterized model Model position offset (asec)
               "modelParm"  if ["Image"]=None, Parameterized model Model parameters
                            (maj, min, pa, type)
    timeRange  = timerange of data to use
    solInt     = Calib solution interval (min)
    smoTime    = Smoothing time applied to SN table (hr) if >0.0
    FreqID     = Frequency group identifier
    minSNR     = minimum acceptable SNR in Calib
    refAnts    = List of reference antennas
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    doTwo      = If True, use one and two baseline combinations
                 for delay calibration, else only one baseline
    doPlot     = If True make plots of SN gains
    plotFile   = Name of postscript file for plots
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine parallel hand group delays"
    printMess(mess, logfile)

    # Open and close image to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)

    # Set output (new) SN table
    SNver = uv.GetHighVer("AIPS SN")+1

    calib = ObitTask.ObitTask("Calib")
    calib.taskLog  = logfile
    if not check:
        setname(uv,calib)
    calib.flagVer   = flagVer
    calib.timeRange = timeRange
    calib.doCalib   = doCalib
    calib.gainUse   = gainUse
    calib.doBand    = doBand
    calib.BPVer     = BPVer
    calib.solMode   = "DELA"
    calib.solType   = "  "
    calib.solInt    = solInt
    calib.minSNR    = minSNR
    calib.refAnts   = refAnts
    calib.solnVer   = SNver
    calib.noScrat   = noScrat
    calib.nThreads  = nThreads
    calib.doTwo     = doTwo
    #  If multiple sources given with models - loop
    if (type(calSou)==list) and (len(calSou)>1):
        ncloop = len(calSou)
        calsou = calSou[0] # setup for first
    else:
        ncloop = 1;
        calsou = calSou
    # Loop
    for ical in range(0,ncloop):
        if type(calsou)==list:
            calib.Sources[0] = calsou[ical]
        else:
            calib.Sources = [calsou]
        # Get model details
        if CalModel and calib.Sources[0] in CalModel:
            Model = CalModel[calib.Sources[0]]
            if Model["image"]:
                if not check:
                    set2name(Model["image"], calib)
            if "nfield" in Model:
                calib.nfield    = Model["nfield"]
            if "CCVer" in Model:
                calib.CCVer     = Model["CCVer"]
            if "BComp" in Model:
                calib.BComp     = Model["BComp"]
            if "EComp" in Model:
                calib.EComp     = Model["EComp"]
            if "Cmethod" in Model:
                calib.Cmethod   = Model["Cmethod"]
            if "Cmodel" in Model:
                calib.Cmodel    = Model["Cmodel"]
            if "Flux" in Model:
                calib.Flux      = Model["Flux"]
            if "modelFlux" in Model:
                calib.modelFlux = Model["modelFlux"]
            if "modelPos" in Model:
                calib.modelPos  = Model["modelPos"]
            if "modelParm" in Model:
                calib.modelParm = Model["modelParm"]
        if debug:
            calib.prtLv =6
            calib.i
            calib.debug = debug
        # Trap failure
        try:
            if not check:
                calib.g
        except Exception, exception:
            print exception
            mess = "Calib Failed retCode= "+str(calib.retCode)
            printMess(mess, logfile)
            return None
        else:
            pass
        # Setup for next if looping
        if ical<(ncloop-1):
            calsou = calSou[ical+1]
    # End loop over calibrators
    
    # Open/close UV to update header
    uv.Open(UV.READONLY,err)
    uv.Close(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1

    # Smooth if requested
    if smoTime>0.0:
        SNver = uv.GetHighVer("AIPS SN")
        snsmo = ObitTask.ObitTask("SNSmo")
        if not check:
            setname(uv, snsmo)
        snsmo.solnIn     = SNver
        snsmo.solnOut    = SNver+1
        snsmo.smoType    = 'VLBI'
        snsmo.smoFunc    = 'MWF '
        snsmo.smoParm    = [smoTime,smoTime,smoTime,smoTime,smoTime]
        snsmo.doBlank    = True
        snsmo.refAnt     = refAnts[0]
        snsmo.taskLog    = logfile
        snsmo.debug      = debug
        #snsmo.debug     = True  # DEBUG
        #bombaroonie     = BombsAwayWithCurtisLemay # DEBUG
        if debug:
            snsmo.i
        mess = "EVLADelayCal: SNSmo SN "+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
        printMess(mess, logfile)
        # Trap failure
        try:
            if not check:
                snsmo.g
        except Exception, exception:
            print exception
            mess = "SNSmo Failed retCode="+str(snsmo.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
    # End SNSmo

    # Plot fits?
    if doPlot:
        retCode = EVLAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="DELA", \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
  
        retCode = EVLAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="PHAS", \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        retCode = EVLAWritePlots (uv, 1, 0, plotFile, err, \
                                      logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        
    # end SN table plot
    # Apply to CL table
    retCode = EVLAApplyCal(uv, err, maxInter=600.0, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode
    
    # Open and close image to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)
    return 0
    # end EVLADelayCal

def EVLACalAP(uv, target, ACal, err, \
              PCal=None, FQid=0, calFlux=None, \
              doCalib=-1, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
              calModel=None, calDisk=0, \
              solnver=0, solInt=10.0/60.0, solSmo=0.0, nThreads=1, refAnt=0, ampScalar=False, \
              doPlot=False, plotFile="./APCal.ps", \
              check=False, debug = False, noScrat=[], logfile = ""):
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
    flagVer  = Flagging table to apply
    solnver  = output SN table version (+1 if smooth), 0=>new
    solInt   = solution interval (min)
    solSmo   = if solSmo<solInt smooth solutions to solSmo
    nThreads = Number of threads to use
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    doPlot   = If True make plots of solutions
    plotFile = Name of postscript file for plots
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    noScrat  = list of disks to avoid for scratch files
    logfile  = Log file for tasks
    """
    ################################################################
    mess =  "Amplitude and phase calibration"
    printMess(mess, logfile)

    # Run SetJy
    setjy = ObitTask.ObitTask("SetJy")
    setjy.taskLog  = logfile
    if not check:
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
    except Exception, exception:
        print exception
        mess = "SetJy Failed retCode="+str(setjy.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass

    # output SN version
    if solnver<=0:
        # Open and close image to sync with disk 
        uv.Open(UV.READONLY, err)
        uv.Close(err)
        solnVer  = max(1,uv.GetHighVer("AIPS SN")+1)
    else:
        solnVer  = solnver

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
            except Exception, exception:
                print exception
                mess = "SetJy Failed retCode="+str(setjy.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
    # Calib on Amp cal if not in PCal
    calib = ObitTask.ObitTask("Calib")
    calib.taskLog  = logfile
    if not check:
        setname(uv,calib)
    calib.Sources  = [ACal]
    calib.flagVer  = flagVer
    calib.ampScalar= ampScalar
    calib.doCalib  = doCalib
    calib.gainUse  = gainUse
    calib.doBand   = doBand
    calib.BPVer    = BPVer
    calib.solMode  = "A&P"
    calib.nThreads = nThreads
    calib.solInt   = solInt
    calib.refAnts  = [refAnt]
    calib.noScrat  = noScrat
    calib.solnVer  = solnVer
    # Given model?
    if calModel:
        calib.DataType2 = "FITS"
        calib.in2File   = calModel
        calib.in2Disk   = calDisk
        calib.nfield    = 1
        calib.CCVer     = 0
        calib.Cmethod   = "DFT"
        calib.Cmodel    = "COMP"
    if debug:
        calib.i
        calib.debug = debug
        #calib.prtLv = 5
    # Trap failure
    try:
        if not check:
            calib.g
    except Exception, exception:
        print exception
        mess = "Calib Failed retCode= "+str(calib.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    # Run if Amp cal if not in PCal
    if not ACal in PCal:
        # ClCal CL1 + SN1 => Cl2
        clcal=ObitTask.ObitTask("CLCal")
        clcal.taskLog  = logfile
        if not check:
            setname(uv,clcal)
        clcal.calSour = [ACal]
        if not check:
            # Open and close image to sync with disk 
            uv.Open(UV.READONLY, err)
            uv.Close(err)
            clcal.calIn   = uv.GetHighVer("AIPS CL")
        else:
            clcal.calIn   = 1
        clcal.calOut    = clcal.calIn+1
        clcal.interMode = "2PT"
        clcal.interMode = "SELF"  # DEBUG
        clcal.FreqID    = FQid
        
    # Calib on phase reference if given
        if PCal:
            if type(PCal)==list:
                slist = PCal
            else:
                slist = [PCal]
            calib.in2File   = "    "
            calib.nfield    = 0
            calib.flagVer   = flagVer
            calib.ampScalar = ampScalar
            calib.doCalib   = 2
            calib.solMode   = "A&P"
            calib.Cmethod   = "DFT"
            calib.Cmodel    = " "
            if debug:
                calib.i
                calib.debug = debug
            # Loop over sources
            for s in slist:
                calib.Sources = [s]
                # Trap failure
                try:
                    if not check:
                        calib.g
                except Exception, exception:
                    print exception
                    mess = "Calib Failed retCode= "+str(calib.retCode)
                    printMess(mess, logfile)
                    return 1
                else:
                    pass
            # end source list

        solnVer2 = calib.solnVer
        # Smoothing?   
        if solSmo>solInt:
            # Get list of calibrators
            smoList = [ACal]
            if PCal:
                if type(PCal)==list:
                    for s in PCal:
                        smoList.append(s)
                else:
                    smoList.append(PCal)
            solnVer2 = solnVer+1
            snsmo=ObitTask.ObitTask("SNSmo")
            snsmo.taskLog  = logfile
            if not check:
                setname(uv,snsmo)
            snsmo.solnIn  = solnVer
            snsmo.solnOut = solnVer2
            snsmo.smoType = "BOTH"
            snsmo.smoType = "MWF"
            snsmo.refAnt  = refAnt
            snsmo.clipSmo = [24.]  # Clip wild amplitudes
            snsmo.clipParm= [100.0]
            mess = "Smooth SN "+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
            printMess(mess, logfile)
            if debug:
                snsmo.i
                snsmo.debug = debug
            # run on all sources for clipping then on individual cal.
            # Trap failure
            try:
                if not check:
                    snsmo.g
            except Exception, exception:
                print exception
                mess = "SNSmo Failed retCode="+str(snsmo.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
            snsmo.clipSmo = [solSmo/60.0] 
            snsmo.clipParm= [0.5]
            # Next time smooth
            snsmo.smoParm = [solSmo/60., solSmo/60.]
            snsmo.solnIn  = solnVer2
            snsmo.solnOut = solnVer2+1
            solnVer2     +=1
            # Loop over sources
            for s in smoList:
                snsmo.Sources[0] = s
                # Trap failure
                try:
                    if not check:
                        snsmo.g
                except Exception, exception:
                    print exception
                    mess = "SNSmo Failed retCode="+str(snsmo.retCode)
                    printMess(mess, logfile)
                    return 1
                else:
                    pass
        else:   # No smoothing
            solnVer2 = solnVer
        # GetJy to set flux density scale if ACal not in PCal
        if not ACal in PCal:
            getjy = ObitTask.ObitTask("GetJy")
            getjy.taskLog  = logfile
            if not check:
                setname(uv,getjy)
            getjy.calSour=[ACal]
            getjy.solnVer = solnVer2
            if type(PCal)==list:
                getjy.Sources=PCal
            else:
                getjy.Sources=[PCal]
                getjy.FreqID = FQid
            #getjy.debug = True # DEBUG
            #bombaroonie = BombsAwayWithCurtisLemay # DEBUG
            if debug:
                getjy.i
                getjy.debug = debug
            # Trap failure
            try:
                if not check:
                    getjy.g
            except Exception, exception:
                print exception
                mess = "GetJy Failed retCode="+str(getjy.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
       
    # Plot gain corrections?
    if doPlot:
        # Amplitude corrections
        retCode = EVLAPlotTab(uv, "SN", solnVer2, err, nplots=6, optype="AMP ", \
                              logfile=logfile, check=check, debug=debug)
        # Phase corrections
        retCode = EVLAPlotTab(uv, "SN", solnVer2, err, nplots=6, optype="PHAS", \
                              logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        # R-L phase corrections
        retCode = EVLAPlotTab(uv, "SN", solnVer2, err, nplots=6, optype="PHAS", stokes="DIFF", \
                              logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        retCode = EVLAWritePlots (uv, 1, 0, plotFile, err, \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
    # end SN table plot

    # Set up for CLCal - only use phase calibrators
    if not check:
        # Open and close image to sync with disk 
        uv.Open(UV.READONLY, err)
        uv.Close(err)
        clcal.solnVer = uv.GetHighVer("AIPS SN")
    else:
        clcal.solnVer = 1
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
    mess = "Update CL "+str(clcal.calIn)+" with SN "+str(clcal.solnVer)+" to CL "+str(clcal.calOut)
    printMess(mess, logfile)

    if debug:
        clcal.i
        clcal.debug = debug
    # Trap failure
    try:
        if not check:
            clcal.g
    except Exception, exception:
        print exception
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
                  solInt1=0.0, solInt2=0.0, solMode="A&P", refAnt=0, \
                  ampScalar=False, specIndex=0.0, \
                  doAuto=False, doPol=False, avgPol=False, avgIF=False, \
                  doPlot=False, plotFile="./BPCal.ps", \
                  check=False, debug = False, nThreads=1, noScrat=[], logfile = ""):
    """ Bandbass calibration

    Do bandbass calibration, write BP table
    Returns task error code, 0=OK, else failed
    Calibration is done in two passes
    1)  First a wideband phase only calibration using channels
    BChan1 to EChan1 or the central doCenter1 fraction of the band
    using a solution interval of solInt1.  This solution is applied
    to all selected data and used in the second pass.
    2)  Second channels in the range BChan2 to EChan2 averaging blocks
    of ChWid2 are calibrated using solType and solMode for solInt2 and
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
    solMode  = solution mode 'A&P', 'P', 'P!A'
    refAnt   = Reference antenna
    ampScalar= If true, scalar average data in calibration
    specIndex= spectral index of calibrator (steep=-0.70)
    doAuto   = Use autocorrelation spectra? Else, crosscorrelation
    doPol    = Apply polarization cal?
    avgPol   = Avg. poln. in solutions?
    avgIF    = Avg. IFs. in solutions?
    doPlot   = If True make plots of corrected data
    plotFile = Name of postscript file for plots
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    nThreads = Number of threads to use
    noScrat  = list of AIPS disks to avoid for scratch files
    logfile  = Log file for task
    """
    ################################################################
    mess =  "Bandpass calibrate data"
    printMess(mess, logfile)
    bpass = ObitTask.ObitTask("BPass")
    bpass.taskLog = logfile
    if not check:
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
    bpass.solMode   = solMode
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
    if not check:
        d     = uv.Desc.Dict
        nchan = d["inaxes"][d["jlocf"]]
    else:
        nchan = 1
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
    except Exception, exception:
        print exception
        mess = "BPass Failed retCode="+str(bpass.retCode)
        return 1
        printMess(mess, logfile)
    else:
        pass
    # Plot corrected data?
    if doPlot:
        scr = EVLASpecPlot( uv, bpass.Sources[0], timerange, refAnt, err, \
                            Stokes=["RR","LL"], doband=1,          \
                            plotFile=plotFile, check=check, logfile=logfile )
        if not UV.PIsA(scr):
            return 1
        retCode = EVLAWritePlots (scr, 1, 0, plotFile, err, \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        scr.Zap(err)
        # end plots
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
    split.taskLog = logfile
    if not check:
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
    except Exception, exception:
        print exception
        mess = "split Failed retCode="+str(split.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end EVLAsplit

def EVLACalAvg(uv, avgClass, avgSeq, CalAvgTime,  err, \
               FQid=0, \
               flagVer=0, doCalib=2, gainUse=0, doBand=1, BPVer=0,  doPol=False, \
               BIF=1, EIF=0, BChan=1, EChan=0, \
               avgFreq=0, chAvg=1, Compress=False, \
               logfile = "", check=False, debug=False):
    """ Calibrate, select and/or average data to a multisource file

    Returns task error code, 0=OK, else failed
    Generates NX and initial dummy CL table if needed
    uv         = UV data object to clear
    avgClass   = Class name of averaged data
    avgSeq     = Sequence number of averaged data
    CalAvgTime = Averaging time in sec
    err        = Obit error/message stack
    FQid       = Frequency Id to process, 0=>all
    doCalib    = Apply calibration table, positive=>calibrate
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply previous bandpass cal.
    BPVer      = previous Bandpass table (BP) version
    doPol      = Calibrate polarization?
    BIF        = first IF to copy
    EIF        = highest IF to copy
    BChan      = first channel to copy
    EChan      = highest channel to copy
    flagVer    = Input Flagging table version
    avgFreq    = If 0 < avgFreq <= 1 then average channels
    chAvg      = Number of channels to average
    Compress   = Write "Compressed" data?
    logfile    = Log file for task
    check      = Only check script, don't execute tasks
    debug      = Run tasks debug, show input
    """
    ################################################################
    mess =  "Average/calibrate data"
    printMess(mess, logfile)
    splat=ObitTask.ObitTask("Splat")
    splat.taskLog = logfile
    if not check:
        setname(uv,splat)
    splat.doCalib  = doCalib
    splat.gainUse  = gainUse
    splat.doBand   = doBand
    splat.BPVer    = BPVer
    splat.doPol    = doPol
    splat.BIF      = BIF
    splat.EIF      = EIF
    splat.BChan    = BChan
    splat.EChan    = EChan
    splat.flagVer  = flagVer
    splat.FreqID   = FQid
    splat.timeAvg  = CalAvgTime
    splat.avgFreq  = avgFreq
    splat.chAvg    = chAvg
    splat.Compress = Compress
    splat.outClass = avgClass
    splat.outDisk  = splat.inDisk
    splat.outSeq   = avgSeq
    if debug:
        splat.i
        splat.debug = debug
    # Trap failure
    try:
        if not check:
            splat.g
    except Exception, exception:
        print exception
        mess = "Splat Failed retCode="+str(splat.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    # end average

    # Get calibrated/averaged data, index and make CL table 1 if doCalib>0
    if not check:
        try:
            uvc = UV.newPAUV("AIPS UV DATA", splat.inName, avgClass, splat.inDisk, avgSeq, True, err)
            if err.isErr:
                print "Error creating cal/avg AIPS data"
                OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
            # Dummy CL table
            solint = splat.timeAvg * 2   # CL table interval twice averaging
            # Open and close image to sync with disk 
            uvc.Open(UV.READONLY, err)
            uvc.Close(err)
            hiver = uvc.GetHighVer("AIPS CL")
            if (doCalib>0) or (hiver<=0):
                UV.PTableCLGetDummy(uvc, uvc, 0, err, solInt=solint)
            if err.isErr:
                print "Error creating cal/avg AIPS data CL table"
                OErr.printErrMsg(err, "Error creating cal/avg AIPS data CL table")
            # Index
            UV.PUtilIndex (uvc, err)
            if err.isErr:
                print  "Error indexing cal/avg AIPS data"
                OErr.printErrMsg(err, "Error indexing cal/avg AIPS data")
        except Exception, exception:
            print exception
            OErr.printErr(err)
            mess = "Indexing or creating CL table failed"
            printMess(mess, logfile)
            return 1
        else:
            pass
    return 0
    # end EVLACalAvg
    
def EVLACalAvg2(uv, avgClass, avgSeq, CalAvgTime,  err,  FQid=0, \
                flagVer=0, doCalib=2, gainUse=0, doBand=1, BPVer=0, doPol=False,  \
                BIF=1, EIF=0, BChan=1, EChan=0, chAvg=1, Compress=False, \
                logfile = "", check=False, debug=False):
    """ Calibrate and average data to a multisource file

    Returns task error code, 0=OK, else failed
    Generates NX and initial dummy CL table
    uv         = UV data object to clear
    avgClass   = Class name of averaged data
    avgSeq     = Sequence number of averaged data
    CalAvgTime = Averaging time in sec
    err        = Obit error/message stack
    FQid       = Frequency Id to process, 0=>all
    doPol      = Calibrate polarization?
    doCalib    = Apply calibration table, positive=>calibrate
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply previous bandpass cal.
    BPVer      = previous Bandpass table (BP) version
    BIF        = first IF to copy
    EIF        = highest IF to copy
    BChan      = first channel to copy
    EChan      = highest channel to copy
    flagVer    = Input Flagging table version
    Compress   = Write "Compressed" data?
    logfile    = Log file for task
    check      = Only check script, don't execute tasks
    debug      = Run tasks debug, show input
    """
    ################################################################
    mess =  "Average/calibrate calibrate data"
    printMess(mess, logfile)
    # Create output
    if not check:
        # Set calibration, editing and selection
        info = uv.List
        info.set("doCalSelect",  True)
        info.set("FreqID",  FQid)
        info.set("doPol",   doPol)
        info.set("doCalib", doCalib)
        info.set("gainUse", gainUse)
        info.set("doBand",  doBand)
        info.set("doPol",   doPol)
        info.set("BPVer",   BPVer)
        info.set("BIF",     BIF)
        info.set("EIF",     EIF)
        info.set("BChan",   BChan)
        info.set("EChan",   EChan)
        info.set("flagVer", flagVer)
        info.set("Compress", Compress,)
        #print "info", info.Dict  # DEBUG
        # Open and close to set
        uv.Open(UV.READCAL, err)
        outuv = UV.newPAUV("CalAvg", uv.Aname, avgClass, uv.Disk, avgSeq, False, err)
        uv.Clone (outuv, err)
        uv.Close(err)
        #outuv.Header(err) # debug
        if err.isErr:
            print "Error creating cal/avg AIPS uv data"
            OErr.printErrMsg(err, "Error creating cal/avg AIPS data")
    
    # Average
    if not check:
        try:
            mess = "Copy/average/calibrate data to "+\
                   outuv.Aname+" . "+outuv.Aclass+" . "+str(outuv.Disk)+ \
                   " . "+str(outuv.Aseq)+" cno: "+str(outuv.Acno)
            printMess(mess, logfile)
            info = outuv.List
            info.set("Compress", Compress,)
            UV.PUtilAvgT (uv, outuv, err, timeAvg=CalAvgTime/60.)
            if err.isErr:
                print "Error cal/avg AIPS uv data"
                OErr.printErrMsg(err, "Error cal/avg AIPS data")
                

            # Do History - previous already copied
            outuv = UV.newPAUV("CalAvg", uv.Aname, avgClass, uv.Disk, avgSeq, True, err)
            #print "DEBUG Copy history"
            inHistory  = History.History("inhistory",  uv.List, err)
            outHistory = History.History("outhistory", outuv.List, err)

            # Add history
            #print "DEBUG Add history"
            outHistory.Open(History.READWRITE, err)
            outHistory.TimeStamp(" Start Obit CalAvg",err)
            outHistory.WriteRec(-1,"CalAvg  CalAvgTime = "+str(CalAvgTime),err)
            outHistory.WriteRec(-1,"CalAvg  inName = "+uv.Aname,  err)
            outHistory.WriteRec(-1,"CalAvg  inClass = "+uv.Aclass, err)
            outHistory.WriteRec(-1,"CalAvg  inDisk = " +str(uv.Disk),err)
            outHistory.WriteRec(-1,"CalAvg  inSeq = " +str(uv.Aseq),err)
            outHistory.WriteRec(-1,"CalAvg  FreqID = "+str(FQid),err)
            outHistory.WriteRec(-1,"CalAvg  doPol = "+str(doPol),err)
            outHistory.WriteRec(-1,"CalAvg  doCalib = "+str(doCalib),err)
            outHistory.WriteRec(-1,"CalAvg  gainUse = "+str(gainUse),err)
            outHistory.WriteRec(-1,"CalAvg  doBand = "+str(doBand),err)
            outHistory.WriteRec(-1,"CalAvg  BPVer = "+str(BPVer),err)
            outHistory.WriteRec(-1,"CalAvg  BIF = "+str(BIF),err)
            outHistory.WriteRec(-1,"CalAvg  EIF = "+str(EIF),err)
            outHistory.WriteRec(-1,"CalAvg  BChan = "+str(BChan),err)
            outHistory.WriteRec(-1,"CalAvg  EChan = "+str(EChan),err)
            outHistory.WriteRec(-1,"CalAvg  flagVer = "+str(flagVer),err)
            outHistory.WriteRec(-1,"CalAvg  Compress = "+str(Compress),err)
            outHistory.Close(err)
            #print "DEBUG Copy history done"
            if err.isErr:
                print "Error cal/avg History"
                OErr.printErrMsg(err, "Error cal/avg History")
                # end copy+history
        except Exception, exception:
            print exception
            OErr.printErr(err)
            mess = "Calibrate and average uv data failed"
            printMess(mess, logfile)
            return 1
        else:
            pass
   
    # Index and make CL table
    if not check:
        try:
            # Dummy CL table
            solint = 2 * CalAvgTime/60.   # CL table interval twice averaging
            UV.PTableCLGetDummy(outuv, outuv, 0, err, solInt=solint)
            if err.isErr:
                print "Error creating cal/avg AIPS data CL table"
                OErr.printErrMsg(err, "Error creating cal/avg AIPS data CL table")
            # Index
            UV.PUtilIndex (outuv, err)
            if err.isErr:
                print  "Error indexing cal/avg AIPS data"
                OErr.printErrMsg(err, "Error indexing cal/avg AIPS data")
        except Exception, exception:
            print exception
            OErr.printErr(err)
            mess = "Indexing or creating CL table failed"
            printMess(mess, logfile)
            return 1
        else:
            pass
    return 0
    # end EVLACalAvg2
    
def EVLASetImager (uv, target, outIclass="", nThreads=1, noScrat=[], logfile = ""):
    """ Setup to run Imager or MFImage

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
    img.taskLog = logfile
    if not check:
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

def EVLAPolCal(uv, InsCals, err, RM=0.0, \
               doCalib=2, gainUse=0, doBand=1, BPVer=0, flagVer=-1, \
               soltype="ORI-", fixPoln=False, avgIF=False, \
               solInt=0.0, refAnt=0, \
               pmodel=[0.0,0.0,0.0,0.0,0.0,0.0,0.0], \
               check=False, debug = False, \
               noScrat=[], logfile = ""):
    """ Instrumental Polarization 

    Do Instrumental
    Instrumental cal uses PCAL, R-L cal is done by imaging each IF in Q and U
    and summing the CLEAN components.
    Returns task error code, 0=OK, else failed
    uv       = UV data object to calibrate
    InsCals  = Instrumental poln calibrators, name or list of names
               If None no instrumental cal
    err      = Obit error/message stack
    RM       = Rotation measure of RLCal
    doCalib  = Apply prior calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    doBand   = >0 => apply bandpass calibration
    BPVer    = AIPS BP table to apply
    flagVer  = Input Flagging table version
    soltype  = solution type, "ORI-", "RAPR"
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
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    noScrat  = list of disks to avoid for scratch files
    logfile  = Log file for task
    """
    ################################################################
    mess =  "Instrumental polarization calibration "
    printMess(mess, logfile)
    # Instrumental calibration
    if InsCals!=None:
        pcal = AIPSTask.AIPSTask("pcal")
        pcal.logFile = logfile
        if not check:
            setname(uv,pcal)
        if type(InsCals)==list:
            pcal.calsour[1:] = InsCals
        else:
            pcal.calsour[1:] = [InsCals]
        pcal.docalib = doCalib
        pcal.gainuse = gainUse
        pcal.doband  = doBand
        pcal.bpver   = BPVer
        pcal.flagver = flagVer
        pcal.soltype = soltype
        pcal.solint  = solInt
        pcal.refant  = refAnt
        pcal.prtlev  = 1
        i = 1;
        for d in noScrat:
            pcal.baddisk[i] = d
            i += 1
        if fixPoln:
            # pcal.domodel = 1. # Does not work
            pcal.bparm[10]=1.0
        if avgIF:
            pcal.cparm[1]=1.0
        pcal.pmodel[1:]  = pmodel
        if debug:
            pcal.i
        # Trap failure
        try:
            if not check:
                pcal.g
        except Exception, exception:
            print exception
            mess = "PCAL Failed retCode="+str(pcal.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
        # end instrumental poln cal
    
    return 0
    # End EVLAPolCal

def EVLARLCal(uv, err, \
              RLDCal=None, BChan=1, EChan = 0, ChWid2=1, solInt1=1./6, solInt2=10., \
              RLPCal=None, RLPhase=0.0, RM=0.0, UVRange=[0.0,0.0], \
              timerange = [0.,0.,0.,0.,0.,0.,0.,0.], \
              FQid=0, calcode="    ", doCalib=-1, gainUse=0, \
              doBand=0, BPVer=0, flagVer=-1,  \
              refAnt=0, doPol=-1, FOV=0.05, niter = 100, CleanRad=None, \
              doPlot=False, plotFile="./BPCal.ps", \
              nThreads=1, noScrat=[], logfile = "",check=False, debug = False):
    """ Determine R-L delay and/or phase calibration

    Returns task error code, 0=OK, else failed
    R-L Delay calibration using new BP table, if R-L phase (& RM) known for
    calibrator(s), this also does the R-L phase calibration
    R-L Phase Calibration applies to (new) highest numbered CL table on uv
    uv       = UV data object to clear
    err      = Obit error/message stack
    RLPCal   = R-L (polarization angle) calibrator,
               If None no R-L cal
    RLPhase  = R-L phase of RLPCal (deg) at 1 GHz
    RM       = R-L phase RM (NYI)
    RLDCal   = An array of triplets with R-L calibrators:
               (name, R-L phase (deg at 1 GHz), RM (rad/m**2))
               If None no R-L delay cal
    solInt1  = first solution interval (min), 0=> scan average
    solInt2  = second solution interval (min)
    RLDPhase = R-L phase of RLPCal (deg) at 1 GHz
    BChan    = First (1-rel) channel number
    EChan    = Highest channel number. 0=> high in data.
    ChWid2   = Number of channels in running mean RL BP soln, 
    UVRange  = Range of baseline used in kilowavelengths
    FQid     = Frequency Id to process
    calcode  = Calibrator code
    doCalib  = Apply calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    timerange= time range of data (aips format)
    doBand   = If >0.5 apply previous bandpass cal.
    BPVer    = previous Bandpass table (BP) version
    flagVer  = Flagging table to apply
    refAnt   = Reference antenna REQUIRED
    doPol    = Apply polarization cal?
    FOV      = field of view radius (deg) needed to image RLPCal
    niter    = Number  of iterations of CLEAN in R-L cal
    CleanRad = CLEAN radius about center or None=autoWin
    doPlot   = If True make plots of corrected data
    plotFile = Name of postscript file for plots
    noScrat  = list of AIPS disks to avoid for scratch files
    nThreads = Number of threads to use in imaging
    logfile  = Log file for task
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    mess =  "R-L polarization calibration "
    printMess(mess, logfile)

    lbpver = BPVer  # default bandpass in imaging
    # Want R-L delay cal?
    if RLDCal!=None:
        ncal = len(RLDCal)  # How many calibrators?
        rlpass=ObitTask.ObitTask("RLPass")
        rlpass.taskLog = logfile
        if not check:
            setname(uv,rlpass)
            #if Antennas:
            #    i = 0
            #    for a in Antennas:
            #        rlpass.Antennas[i] = a; i  += 1
            rlpass.timeRange[0] = timerange[0]+(timerange[1]+(timerange[2]+timerange[3]/60)/60)/24
            rlpass.timeRange[1] = timerange[4]+(timerange[5]+(timerange[6]+timerange[7]/60)/60)/24
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
            rlpass.doBand  = doBand
            rlpass.BPVer   = BPVer
            rlpass.refAnt  = refAnt
            rlpass.solInt1 = solInt1
            rlpass.solInt2 = solInt2
            rlpass.BPSoln  = 0
            rlpass.prtLv   = 1
            rlpass.nThreads = nThreads
            # Loop over calibrators
            for ical in range (0,ncal):
                rlpass.Sources[0]= RLDCal[ical][0]
                rlpass.RLPhase   = RLDCal[ical][1]
                rlpass.RM        = RLDCal[ical][2]
                mess =  "R-L delay calibration using "+rlpass.Sources[0]
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
        lbpver = uv.GetHighVer("AIPS BP")
    # end R-L delay cal
    
    # R-L phase cal
    if RLPCal!=None:
        mess =  "R-L phase calibration using "+RLPCal
        printMess(mess, logfile)
        img = ObitTask.ObitTask("Imager")
        img.taskLog    = logfile
        if not check:
            setname(uv,img)
        if RLDCal==None:
            img.doBand = doBand
            img.BPVer  = lbpver
        else:
            img.doBand = 1
            img.BPVer  = lbpver  # Just created one
        img.doCalib    = doCalib
        img.gainUse    = gainUse
        img.flagVer    = flagVer
        img.doPol      = True
        img.Sources[0] = RLPCal
        img.Stokes     = "IQU"
        img.FOV        = FOV
        img.Niter      = niter
        # Auto window or centered box
        if CleanRad:
            img.CLEANBox=[-1,CleanRad,0,0]
        else:
            img.autoWindow  = True
        img.dispURL    = "None"
        img.BLFact     = 1.004
        img.Catalog    = "None"
        img.nThreads   = nThreads
        img.maxPSCLoop = 2
        img.minFluxPSC = 0.05
        img.solPInt    = solInt1
        img.solPType   = "L1"
        img.prtLv      = 2
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
        if not check:
            h = uv.Desc.Dict
            if h["jlocif"]>=0:
                nif = h["inaxes"][h["jlocif"]]
            else:
                nif = 1
        else:
            nif = 1

        # Lists of flux densities and RMSes
        IFlux = []
        IRMS  = []
        QFlux = []
        QRMS  = []
        UFlux = []
        URMS  = []
        
        # Loop over IF imaging I,Q, U, allow failure
        for iif in range (1, nif+1):
            img.BIF = iif
            img.EIF = iif
            #img.dispURL    = "ObitView"  # DEBUG
            #img.debug=True               # DEBUG
            if debug:
                img.i
                img.debug = debug
            # Trap failure
            failed = False
            try:
                if not check:
                    img.g
            except Exception, exception:
                print exception
                mess = "Imager Failed IF "+str(iif)+" retCode="+str(img.retCode)
                printMess(mess, logfile)
                failed = True
            else:
                pass
            # Stub if failed
            if failed:
                IFlux.append(-1.0)
                IRMS.append(-1.0)
                QFlux.append(-1.0)
                QRMS.append(-1.0)
                UFlux.append(-1.0)
                URMS.append(-1.0)
                continue

            if check:      # Don't bother if only checking 
                continue
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
                stat = imstat(x, err, blc=blc,trc=trc,logfile=None)
                IFlux.append(stat["Flux"])
                IRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes Q
                outClass="QPOLCL"
                x =  Image.newPAImage("Q",outName, outClass, outDisk,outSeq,True,err)
                stat = imstat(x, err, blc=blc,trc=trc,logfile=None)
                QFlux.append(stat["Flux"])
                QRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes U
                outClass="UPOLCL"
                x =  Image.newPAImage("U",outName, outClass, outDisk,outSeq,True,err)
                stat = imstat(x, err, blc=blc,trc=trc,logfile=None)
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
                stat = imstat(x, err, blc=blc,trc=trc,logfile=None)
                IFlux.append(stat["Flux"])
                IRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes Q
                outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                stat = imstat(x, err, blc=blc,trc=trc,logfile=None)
                IFlux.append(stat["Flux"])
                IRMS.append(stat["RMSHist"])
                x.Zap(err)  # Cleanup
                del x
                # Stokes U
                outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                stat = imstat(x, err, blc=blc,trc=trc,logfile=None)
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
            mess = "%3d  %8.3f %8.4f %7.3f %7.4f %7.3f %7.4f %7.2f "%\
                (i+1, IFlux[i], IRMS[i], QFlux[i], QRMS[i], UFlux[i], URMS[i], cor)
            printMess(mess, logfile)
        # Copy gainUse to new highest CL table
        if not check:
            # Open and close image to sync with disk 
            uv.Open(UV.READONLY, err)
            uv.Close(err)
            hiCL = uv.GetHighVer("AIPS CL")
        else:
            hiCL = 1

        # Copy CL table to be modified
        if not check:
            EVLACopyTable (uv, uv, "AIPS CL", err, inVer=hiCL, outVer=hiCL+1, \
                           logfile=logfile, check=check, debug=debug)
        if err.isErr:
            print  "Error copying CL Table"
            return 1
        
        # Apply R-L phase corrections
        clcor = AIPSTask.AIPSTask("clcor")
        clcor.logFile  = logfile
        if not check:
            setname(uv,clcor)
        clcor.opcode   = "POLR"
        clcor.gainver  = hiCL+1
        clcor.gainuse  = hiCL+1
        clcor.clcorprm[1:] = RLCor
        if debug:
            clcor.i
            clcor.debug = debug
        # Trap failure
        try:
            if not check:
                clcor.g
        except Exception, exception:
            print exception
            mess = "CLCOR Failed retCode="+str(clcor.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
        # end R-L Cal
    # Plot corrected data?
    if doPlot:
        if RLPCal:
            pSou = RLPCal
        else:
            pSou = RLDCal[0][0]
        scr = EVLASpecPlot( uv, pSou, timerange, refAnt, err,    \
                                Stokes=["RL","LR"], doband=1,   \
                                plotFile=plotFile, \
                                check=check, debug=debug, logfile=logfile )
        if not scr.UVIsA():
            return 1
        retCode = EVLAWritePlots (scr, 1, 0, plotFile, err, \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        scr.Zap(err)
        # end plots
    return 0
    # end EVLARLCal

def EVLARLCal2(uv, err, uv2 = None, \
               RLDCal=None, BChan=1, EChan = 0,  \
               FQid=0, calcode="    ", doCalib=-1, gainUse=0, \
               timerange = [0.,0.,0.,0.,0.,0.,0.,0.], \
               doBand=0, BPVer=0, flagVer=-1, \
               refAnt=0, doPol=-1, smooth=[0.,0.,0.], dataInt=0., \
               RLPCal=None,  FOV=0.05, niter = 100, \
               nThreads=1, noScrat=[], logfile = "",check=False, debug = False):
    """ Determine R-L delay and phase calibration

    Returns task error code, 0=OK, else failed
    Calibration applies to (new) highest numbered CL table on uv
    uv       = UV data object to clear
    err      = Obit error/message stack
    uv2      = If gives, then copy AN table from uv to uv2 and apply same
               calibration (intended to calibrate CVel data)
    RLPCal   = An array of triplets with R-L calibrators:
               (name, R-L phase (deg at 1 GHz), RM (rad/m**2))
               If None no R-L cal
    RLDCal   = R-L delay calibrator name or list
               If None no R-L delay cal
    BChan    = First (1-rel) channel number
    EChan    = Highest channel number. 0=> high in data.
    FQid     = Frequency Id to process
    calcode  = Calibrator code
    doCalib  = Apply calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    timerange= time range of data (aips format)
    doBand   = If >0.5 apply previous bandpass cal.
    BPVer    = previous Bandpass table (BP) version
    flagVer  = Flagging table to apply
    refAnt   = Reference antenna REQUIRED
    doPol    = Apply polarization cal?
    smooth   = Channel smoothing function
    dataInt  = Data integration time (sec)
    FOV      = field of view radius (deg) needed to image RLPCal
    niter    = Number  of iterations of CLEAN in R-L cal
    noScrat  = list of AIPS disks to avoid for scratch files
    nThreads = Number of threads to use in imaging
    logfile  = Log file for task
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    mess =  "R-L polarization calibration "
    printMess(mess, logfile)
    # Want R-L delay cal?
    if RLDCal!=None:
        rldly=ObitTask.AIPSTask("rldly")
        rldly.logFile = logfile
        if not check:
            setname(uv,rldly)
        if type(RLDCal)!=list:
            rldly.calsour[1]=RLDCal
        else:
            i = 1
            for t in RLDCal:
                rldly.calsour[i] = t
                i  += 1
        i = 1
        for t in timerange:
            rldly.timerang[i] = t
            i  += 1
        rldly.bchan   = BChan
        rldly.echan   = EChan
        rldly.docalib = doCalib
        rldly.gainuse = gainUse
        rldly.flagver = flagVer
        rldly.freqid  = FQid
        rldly.calcode = calcode
        rldly.dopol   = doPol
        rldly.smooth[1]=smooth[0]; rldly.smooth[2]=smooth[1];rldly.smooth[3]=smooth[2];
        rldly.doband  = doBand
        rldly.bpver   = BPVer
        rldly.flagver = flagVer
        rldly.refant  = refAnt
        rldly.solint  = dataInt
        if debug:
            print "timerange", rldly.timerang
            rldly.i
        # Trap failure
        try:
            if not check:
                rldly.g
        except Exception, exception:
            print exception
            mess = "rldly Failed retCode="+str(rldly.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
        # Get new CL table number
        if not check:
            # Open and close image to sync with disk 
            uv.Open(UV.READONLY, err)
            uv.Close(err)
            gainUse = uv.GetHighVer("AIPS CL")
            
    # end R-L delay cal
    
    # R-L phase cal
    if RLPCal!=None:
        ncal = len(RLPCal)  # How many calibrators? 
        img = ObitTask.ObitTask("Imager")
        img.taskLog    = logfile
        if not check:
            setname(uv,img)
        img.doCalib    = doCalib
        img.gainUse    = gainUse
        img.flagVer    = flagVer
        img.doPol      = True
        img.Stokes     = "IQU"
        img.FOV        = FOV
        img.Niter      = niter
        img.autoWindow = True
        img.dispURL    = "None"
        img.Catalog    = "None"
        img.nThreads   = nThreads
        img.noScrat    = noScrat
        img.prtLv      = 2
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
        if not check:
            h = uv.Desc.Dict
            if h["jlocif"]>=0:
                nif = h["inaxes"][h["jlocif"]]
            else:
                nif = 1
        else:
            nif = 1

        
        # Loop over calibrators
        SouCal = []
        for ical in range (0,ncal):
            img.Sources[0]= RLPCal[ical][0]
            #rlpass.RLPhase   = RLPCal[ical][1]
            #rlpass.RM        = RLPCal[ical][2]
            # Loop over IF imaging I,Q, U
            # Lists of flux densities and RMSes
            IFlux = []
            IRMS  = []
            QFlux = []
            QRMS  = []
            UFlux = []
            URMS  = []
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
                except Exception, exception:
                    print exception
                    mess = "Imager Failed retCode="+str(img.retCode)
                    printMess(mess, logfile)
                    return 1
                else:
                    pass
    
                if check:      # Don't bother if only checking 
                    continue
                # Get fluxes from Summed CCs, RMS from inner quarter of images
                if img.DataType=="AIPS":
                    outName = (img.Sources[0].strip()+"TEMP")[0:12]
                    outDisk = img.outDisk
                    outSeq  = 6666
                    # Stokes I
                    outClass="IPOLCL"
                    x =  Image.newPAImage("I",outName, outClass, outDisk,outSeq,True,err)
                    SumCC = EVLAGetSumCC (x,err)
                    h = x.Desc.Dict
                    blc = [h["inaxes"][0]/4,h["inaxes"][1]/4]
                    trc = [3*h["inaxes"][0]/4,3*h["inaxes"][1]/4]
                    stat = imstat(x, err, blc=blc,trc=trc, logfile=logfile)
                    IFlux.append(SumCC)
                    IRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes Q
                    outClass="QPOLCL"
                    x =  Image.newPAImage("Q",outName, outClass, outDisk,outSeq,True,err)
                    SumCC = EVLAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc, logfile=logfile)
                    QFlux.append(SumCC)
                    QRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes U
                    outClass="UPOLCL"
                    x =  Image.newPAImage("U",outName, outClass, outDisk,outSeq,True,err)
                    SumCC = EVLAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc, logfile=logfile)
                    UFlux.append(SumCC)
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
                    SumCC = EVLAGetSumCC (x,err)
                    h = x.Desc.Dict
                    blc = [h["inaxes"][0]/4,h["inaxes"][1]/4]
                    trc = [3*h["inaxes"][0]/4,3*h["inaxes"][1]/4]
                    stat = imstat(x, err, blc=blc,trc=trc, logfile=logfile)
                    IFlux.append(SumCC)
                    IRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes Q
                    outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                    x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                    SumCC = EVLAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc, logfile=logfile)
                    QFlux.append(SumCC)
                    QRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes U
                    outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                    x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                    SumCC = EVLAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc, logfile=logfile)
                    UFlux.append(SumCC)
                    URMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    out2File = img.Sources[0].strip()+"TEMPPOLCAL2.uvtab"
                    u =  UV.newPFUV("UV",outFile,img.outDisk,True,err)
                    u.Zap(err)
                    del u
               # End accumulate statistics by file type
            # End loop over IF
            # Save source info
            SouCal.append({"name":img.Sources[0],"Phase":RLPCal[ical][1],"RM":RLPCal[ical][2], \
                               "IFlux":IFlux, "IRMS":IRMS, "QFlux":QFlux, "QRMS":QRMS, \
                               "UFlux":UFlux, "URMS":URMS})
        # end loop over calibrators

        # Give results, weighted compute R-L correction
        import math
        mess = '\n R-L Phase calibration results'
        printMess(mess, logfile)
        RLCor = []
        RLCorRSum = []
        RLCorISum = []
        RLCorWt   = []
        # Zero accumulators
        for i in range (0,len(IFlux)):
            RLCorRSum.append(0.0)
            RLCorISum.append(0.0)
            RLCorWt.append(0.0)
        
        for ical in range (0,ncal):
            IFlux   = SouCal[ical]["IFlux"]
            IRMS    = SouCal[ical]["IRMS"]
            QFlux   = SouCal[ical]["QFlux"]
            QRMS    = SouCal[ical]["QRMS"]
            UFlux   = SouCal[ical]["UFlux"]
            URMS    = SouCal[ical]["URMS"]
            RLPhase = SouCal[ical]["Phase"]
            RM      = SouCal[ical]["RM"]
            mess = SouCal[ical]["name"]+"\n IF     IFlux    IRMS    QFlux   QRMS    UFlux  URMS   R-L Corr     Wt"
            printMess(mess, logfile)
            for i in range (0,len(IFlux)):
                # REALLY NEED RM Correction!!!!!
                cor = RLPhase - 57.296 * math.atan2(UFlux[i],QFlux[i])
                if cor>180:
                    cor -= 360.0
                if cor<-180:
                    cor += 360.0
                wt = (QFlux[i]**2 + UFlux[i]**2) /(QRMS[i]**2 + URMS[i]**2)   # weight from SNR
                RLCorRSum[i] += (math.cos(cor/57.296)*wt)
                RLCorISum[i] += (math.sin(cor/57.296)*wt)
                RLCorWt[i]   += wt
                mess = "%3d  %8.3f %8.3f %7.3f %7.3f %7.3f %7.3f %8.3f %7.1f "% \
                    (i+1, IFlux[i], IRMS[i], QFlux[i], QRMS[i], UFlux[i], URMS[i], cor, wt)
                printMess(mess, logfile)
            # Copy gainUse to new highest CL table
            if not check:
                # Open and close image to sync with disk 
                uv.Open(UV.READONLY, err)
                uv.Close(err)
                hiCL = uv.GetHighVer("AIPS CL")
            else:
                hiCL = 1
        # end loop over calibrators

        # Loop making weighted average correction
        mess = "\n\n Weighted average corrections\n IF  R-L Corr"
        printMess(mess, logfile)
        for i in range (0,len(IFlux)):
            if RLCorWt[i]>0.0:
                corr = RLCorRSum[i]
                cori = RLCorISum[i]
                cor = math.atan2(cori,corr)*57.296
            else:
                cor = 0.0
            mess = "%3d  %7.3f "% (i+1, cor)
            printMess(mess, logfile)
            RLCor.append(cor)
        # end loop making weighted average

        # If calibrating second uv data, copy AN table 1
        if uv2:
            z = uv2.ZapTable("AIPS AN",1,err)
            EVLACopyTable (uv, uv2, "AIPS AN", err, \
                           logfile=logfile, check=check, debug=debug)
            if err.isErr:
                print  "Error copying AN Table"
                return 1
            
        # Copy CL table to be modified (CLCOR buggy)
        if not check:
            EVLACopyTable (uv, uv, "AIPS CL", err, inVer=hiCL, outVer=hiCL+1, \
                           logfile=logfile, check=check, debug=debug)
        if err.isErr:
            print  "Error copying CL Table"
            return 1
        
        # Apply R-L phase corrections
        clcor = AIPSTask.AIPSTask("clcor")
        clcor.logFile  = logfile
        if not check:
            setname(uv,clcor)
        clcor.opcode   = "POLR"
        clcor.gainver  = hiCL+1
        clcor.gainuse  = hiCL+1
        clcor.clcorprm[1:] = RLCor
        if debug:
            clcor.i
            clcor.debug = debug
        # Trap failure
        try:
            if not check:
                clcor.g
        except Exception, exception:
            print exception
            mess = "CLCOR Failed retCode="+str(clcor.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
        # If calibrating second uv data, run clcor
        if uv2:
            mess = "Also calibrate Secondary UV data"
            printMess(mess, logfile)
            if not check:
                setname(uv2,clcor)
                # Open and close image to sync with disk 
                uv2.Open(UV.READONLY, err)
                uv2.Close(err)
                hiCL = uv2.GetHighVer("AIPS CL")
                # Copy CL table to be modified (CLCOR buggy)
                EVLACopyTable (uv, uv, "AIPS CL", err, inVer=hiCL, outVer=hiCL+1, \
                               logfile=logfile, check=check, debug=debug)
                if err.isErr:
                    print  "Error copying CL Table"
                    return 1
                
                clcor.gainver  = hiCL+1
                clcor.gainuse  = hiCL+1
            if debug:
                clcor.i
                clcor.debug = debug
            # Trap failure
            try:
                if not check:
                    clcor.g
            except Exception, exception:
                print exception
                mess = "CLCOR Failed retCode="+str(clcor.retCode)
                printMess(mess, logfile)
                return 1
            else:
                pass
        # end R-L Cal
    return 0
    # end EVLARLCal2

def EVLAReportTargets(uv, err,  FreqID=1, Sources=None, seq=1, sclass="IClean", \
                          Stokes="I", logfile='', check=False, debug=False):
    """ Generate report info for a list of targets in AIPS files

    Returns a report which is a list of dicts, each of which contains
        "Source":    Source name
        "haveImage": True if images were made, 
        "ObsDate":   Observing date as "yyyy-mm-dd"
        "numVis":    Number of visibilities (ignoring flagging)
        "Exposure":  Total integration time (day)
        "RA"    :    Source RA (deg) at standard equinox
        "Dec"   :    Source Dec (deg) at standard equinox
        "IFlux" :    Source Table IPol flux per IF
        "QFlux" :    Source Table QPol flux per IF
        "UFlux" :    Source Table UPol flux per IF
        "VFlux" :    Source Table VPol flux per IF
    following present if haveImage True
        "RAPnt" :   Antenna pointing RA (deg) at standard equinox
        "DecPnt":   Antenna pointing Dec (deg) at standard equinox
        "Freq"  :   Reference frequency (Hz)
        "BW"    :   Image bandwidth (Hz)
        "Size"  :   Width of image in deg (From Stokes I)
        "Cells" :   Cell spacing in deg (From Stokes I)
        for each s in Stokes:
            "sSum" : Sum of clean components in Jy
            "sPeak": Peak pixel brightness in Jy
            "sRMS" : RMS noise in inner quarter (Jy)
            "sBeam": Beam (maj, min, PA) (deg)
    
    uv         = UV data object
    err        = Python Obit Error/message stack
    Sources    = Source name or list of names to use
                 If an empty list all sources in uv are included
    seq        = sequence number of images
    sclass     = Image class, first character replaced with char in Stokes
    FreqID     = Frequency group identifier
    Stokes     = Stokes parameters of images
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Generate source statistics "
    printMess(mess, logfile)

    # If list empty get all sources
    if type(Sources)==list:
        sl = Sources
    else:
        sl = [Sources]

    if len(sl)<=0:
        slist = EVLAAllSource(uv,err,logfile=logfile,check=check,debug=debug)
    else:
        slist = sl
    
    # Init output
    Report = []

    # Image disk assumed same as uv
    disk = uv.Disk
    user = OSystem.PGetAIPSuser()
    
    # Loop over slist
    for sou in slist:
        sdict = {"Source":sou, "haveImage":False}  # Init source structure
        sdict["ObsDate"]  = uv.Desc.Dict["obsdat"]
        # Observing stats
        obstat = EVLAGetTimes (uv, sou, err, logfile=logfile, check=check, debug=debug)
        sdict["numVis"]   = obstat["numVis"]
        sdict["Exposure"] = obstat["Exposure"]
        sdict["RA"]       = obstat["RA"]
        sdict["Dec"]      = obstat["Dec"]
        sdict["IFlux"]    = obstat["IFlux"]
        sdict["QFlux"]    = obstat["QFlux"]
        sdict["UFlux"]    = obstat["UFlux"]
        sdict["VFlux"]    = obstat["VFlux"]
        # Test if image exists
        cno = AIPSDir.PTestCNO(disk, user, sou, Stokes[0:1]+sclass[1:], "MA", seq, err)
        if cno <= 0 :
            Report.append(sdict)  # Save source info
            continue
        # Image statistics, loop over Stokes
        for s in Stokes:
            klass = s+sclass[1:]
            x = Image.newPAImage(s, sou, klass, disk, seq, True, err)
            hd = x.Desc.Dict
            sdict[s+"Beam"] = (hd["beamMaj"],hd["beamMin"],hd["beamPA"])
            # Some from Stokes I only
            if s == 'I':
                sdict["haveImage"] = True
                sdict["Size"]    = hd["inaxes"][1]*hd["cdelt"][1]
                sdict["Cells"]   = hd["cdelt"][1]
                sdict["RAPnt"]   = hd["obsra"]
                sdict["DecPnt"]  = hd["obsdec"]
                sdict["Freq"]    = hd["crval"][hd["jlocf"]]
                sdict["BW"]      = hd["cdelt"][hd["jlocf"]]
            blc = [hd["inaxes"][0]/4,hd["inaxes"][1]/4]
            trc = [3*hd["inaxes"][0]/4,3*hd["inaxes"][1]/4]
            stat = imstat(x,err,blc=blc,trc=trc)  # Image statistics inner quarter
            if abs(stat["Max"]) >  abs(stat["Min"]):
                sdict[s+"Peak"] = stat["Max"]
            else:
                sdict[s+"Peak"] = stat["Min"]
            sdict[s+"RMS"]  = stat["RMSHist"]
            sdict[s+"Sum"]  = EVLAGetSumCC(x, err, logfile=logfile, check=check, debug=debug)
        # End stokes image loop
        Report.append(sdict)  # Save source info
    # end loop over sources

    # Give terse listing
    mess = "\n Summary at frequency = "+"%8.3f"%(hd["crval"][hd["jlocf"]]*1.0e-9)+" GHz on "+ \
           uv.Desc.Dict["obsdat"]
    printMess(mess, logfile)
    for sdict in Report:
        mess = "\n Source = "+sdict["Source"]+", Exposure="+"%5.3f"%(sdict["Exposure"]*24.)+" hr"
        printMess(mess, logfile)
        if "IBeam" in sdict:
            mess = "IPol Beam = ("+"%8.3f"%(sdict["IBeam"][0]*3600.0)+", %8.3f"%(sdict["IBeam"][1]*3600.0)+ \
                   ", %6.1f"%(sdict["IBeam"][2])+") asec, asec, deg"
            printMess(mess, logfile)
        else:
            continue   # Nothing to report
        # Source table flux densities
        if "IFlux" in sdict:
            n = len(sdict["IFlux"])
            for i in range(0,n):
                mess = "IF "+str(i+1)+" IPol="+"%8.4f"%(sdict["IFlux"][i])+ \
                       ", QPol="+"%8.4f"%(sdict["QFlux"][i])+ \
                       ", UPol="+"%8.4f"%(sdict["UFlux"][i])+ \
                       ", VPol="+"%8.4f"%(sdict["VFlux"][i])
                printMess(mess, logfile)
        for s in Stokes:
            mess = "Stokes "+s+" Sum CC="+"%8.4f"%(sdict[s+"Sum"])+", Peak="+"%8.4f"%(sdict[s+"Peak"])+ \
                ", RMS="+"%8.5f"%(sdict[s+"RMS"])+" Jy"
            printMess(mess, logfile)
        # Polarization
        if Stokes=="IQU":
            ppolSum  = (sdict["QSum"]**2  + sdict["USum"]**2)**0.5
            ppolPeak = (sdict["QPeak"]**2 + sdict["UPeak"]**2)**0.5
            RLSum    = 57.296*math.atan2(sdict["USum"], sdict["QSum"])
            RLPeak   = 57.296*math.atan2(sdict["UPeak"],sdict["QPeak"])
            mess = "Sum CC PPol="+"%8.4f"%(ppolSum)+", R=L Phase="+"%8.2f"%(RLSum)+ \
                   "; Peak PPol="+"%8.4f"%(ppolPeak)+", R=L Phase="+"%8.2f"%(RLPeak)
            printMess(mess, logfile)
    # End terse listing
    return Report
    # end EVLAReportTargets

def EVLAGetSumCC(image, err, CCver=1,
                 logfile='', check=False, debug=False):
    """ Sum fluxes in a CC table
    
    Sums the flux densities in a CC Table on an image
    Returns sum
    Returns with err set on error
    image      = Image with CC table
    err        = Python Obit Error/message stack
    CCver      = CC table to sum
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
    if check:
        return 0.0
    if debug:
        return 0.0
    # Open and close image to sync with disk 
    image.Open(Image.READONLY, err)
    image.Close(err)
    
    CCTab = image.NewTable(Table.READONLY, "AIPS CC", CCver, err)
    if err.isErr:
        return 0.0
    # Open
    CCTab.Open(Table.READONLY, err)
    if err.isErr:
        return 0.0
    # Number of rows
    nrow    = CCTab.Desc.Dict["nrow"]
    sum     = 0.0
    # Loop over table
    for irow in range (1, nrow+1):
        CCrow = CCTab.ReadRow(irow, err)
        if err.isErr:
            return sum
        sum += CCrow["FLUX"][0]
    # End loop over table
    # Close table
    CCTab.Close(err)
    if err.isErr:
        return sum
    return sum
    # end EVLAGetSumCC

def EVLAGetTimes(uv, Source, err, 
                 logfile='', check=False, debug=False):
    """ Lookup observing times and number of visibilities for a source, other info
    
    Return dict {"numVis":no vis, "Exposure":Total integration time (day),
                 "RA": RA@equinox, "Dec" Dec@equinox,
                 "IFlux":IFlux, "QFlux":QFlux, "UFlux":UFlux, "VFlux":VFlux}
    uv         = UV data with AIPS SU and AIPS NX tables
    Source     = Source to lookup
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
    if check:
        return {"numVis":0, "Exposure":0.0, "RA":0.0, "Dec":0.0}
    # Open and close uv to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)
    
    # Lookup Source ID (SouID)
    SUtab = uv.NewTable(Table.READONLY, "AIPS SU", 1, err)
    SUtab.Open(Table.READONLY, err)
    if err.isErr:
        return  {"numVis":0, "Exposure":0.0, "RA":0.0, "Dec":0.0}
    # Number of rows
    nrow =  SUtab.Desc.Dict["nrow"]
    for i in range (0,nrow):    # Loop over rows
        SUrow = SUtab.ReadRow(i+1, err)
        if err.isErr:
            return  {"numVis":0, "Exposure":0.0, "RA":0.0, "Dec":0.0}
        SouID = SUrow["ID. NO."][0]
        RA    = SUrow["RAEPO"][0]
        Dec   = SUrow["DECEPO"][0]
        IFlux = SUrow["IFLUX"]
        QFlux = SUrow["QFLUX"]
        UFlux = SUrow["UFLUX"]
        VFlux = SUrow["VFLUX"]
        #if debug:
        #    mess="Source "+Source+" test "+SUrow["SOURCE"][0]+" ID ="+str(SouID)+ \
        #        " match="+str(SUrow["SOURCE"][0].rstrip()==Source.rstrip())
        #    printMess(mess, logfile)
        if SUrow["SOURCE"][0].rstrip()==Source.rstrip():   # This it?
            break;
    SUtab.Close(err)
    if err.isErr:
        return  {"numVis":0, "Exposure":0.0, "RA":RA, "Dec":Dec, \
                 "IFlux":IFlux, "QFlux":QFlux, "UFlux":UFlux, "VFlux":VFlux}
    
    # get observing stats from AIPS NX table
    cntVis  = 0
    sumTime = 0.0
    NXTab = uv.NewTable(Table.READONLY, "AIPS NX", 1, err)
    if err.isErr:
        return {"numVis":0, "Exposure":0.0, "RA":RA, "Dec":Dec, \
                 "IFlux":IFlux, "QFlux":QFlux, "UFlux":UFlux, "VFlux":VFlux}
    # Open
    NXTab.Open(Table.READONLY, err)
    if err.isErr:
        return {"numVis":cntVis, "Exposure":sumTime, "RA":RA, "Dec":Dec, \
                 "IFlux":IFlux, "QFlux":QFlux, "UFlux":UFlux, "VFlux":VFlux}
    # Number of rows
    nrow    = NXTab.Desc.Dict["nrow"]
    # Loop over table
    for irow in range (1, nrow+1):
        NXrow = NXTab.ReadRow(irow, err)
        if err.isErr:
            return {"numVis":cntVis, "Exposure":sumTime, "RA":RA, "Dec":Dec, \
                 "IFlux":IFlux, "QFlux":QFlux, "UFlux":UFlux, "VFlux":VFlux}
        #  Is this the desired source?
        if NXrow["SOURCE ID"][0]==SouID:
            sumTime += NXrow["TIME INTERVAL"][0]
            cntVis  += NXrow["END VIS"][0] - NXrow["START VIS"][0] + 1
    # End loop over table
    # Close table
    NXTab.Close(err)
    if err.isErr:
        return {"numVis":cntVis, "Exposure":sumTime, "RA":RA, "Dec":Dec, \
                 "IFlux":IFlux, "QFlux":QFlux, "UFlux":UFlux, "VFlux":VFlux}

    if debug:
        mess="EVLAGetTimes: Source "+Source+"="+str(SouID)+" numVis="+str(cntVis)+ \
            " Integration time = "+"%5.3f"%(sumTime*24.)+" hr"
        printMess(mess, logfile)
 
    return {"numVis":cntVis, "Exposure":sumTime, "RA":RA, "Dec":Dec, \
                 "IFlux":IFlux, "QFlux":QFlux, "UFlux":UFlux, "VFlux":VFlux}
    # end EVLAGetTimes

def EVLAImageTargets(uv, err, Sources=None,  FreqID=1, seq=1, sclass="IClean", band="", \
                     doCalib=-1, gainUse=0, doBand=-1, BPVer=0,  flagVer=-1,  doPol=False, \
                     Stokes="I", FOV=0.1/3600.0, Robust=0, Niter=300, CleanRad=None, \
                     maxPSCLoop=0, minFluxPSC=0.1, solPInt=20.0/60., \
                     solPMode="P", solPType= "  ", \
                     maxASCLoop=0, minFluxASC=0.5, solAInt=2.0, \
                     solAMode="A&P", solAType= "  ", \
                     avgPol=False, avgIF=False, minSNR = 5.0, refAnt=0, \
                     do3D=True, BLFact=0.999, BLchAvg=False, \
                     doMB=False, norder=2, maxFBW=0.05, \
                     nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Image a list of sources with optional selfcal

    Uses Imager or MFImage to image a list of sources.
    Data must be at least approximately calibrated
    Returns task error code, 0=OK, else failed
    
    uv         = UV data object
    err        = Python Obit Error/message stack
    Sources    = Source name or list of names to use
                 If an empty list all sources in uv are included
    seq        = sequence number of output
    sclass     = Image output class
    band       = project band - appended to name
    FreqID     = Frequency group identifier
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    doPol      = Apply polarization cal?
    Stokes     = Stokes parameters to image
    FOV        = Field of view to image in deg
    Robust     = Weighting robustness parameter
    Niter      = max no. iterations
    CleanRad   = CLEAN radius about center or None=autoWin
    maxPSCLoop = max. number of phase sc loops
    minFluxPSC = Trip level for phase self cal (Jy)
    solPInt    = Phase solution interval (min)
    solPMode   = Phase soln mode "P", "DELA"
    solPType   = Phase soln type
    maxASCLoop = max. number of amp&phase sc loops
    minFluxASC = Trip level for amp&phase self cal (Jy)
    solAInt    = Amp&phase solution interval (min)
    solAMode   = Amp&Phase soln mode "A&P", "P", "DELA"
    solAType   = Amp&PPhase soln type
    avgPol     = Average poln in SC?
    avgIF      = Average IFs in SC?
    minSNR     = minimum acceptable SNR in SC
    refAnt     = Reference antenna
    do3D       = Use 3D imaging?
    BLFact     = Baseline dependent averaging factor
    BLchAvg    = If True and BLFact>=1.0 also average channels
    doMB       = If True is wideband imaging
    norder     = order on wideband imaging
    maxFBW     = max. fractional wideband imaging
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
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
    imager.BLFact      = 1.004
    imager.BLchAvg     = BLchAvg
    imager.flagVer     = flagVer
    imager.doCalib     = doCalib
    imager.gainUse     = gainUse
    imager.doBand      = doBand
    imager.BPVer       = BPVer
    imager.doPol       = doPol
    imager.Stokes      = Stokes
    imager.FOV         = FOV
    imager.Robust      = Robust
    imager.Niter       = Niter
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
    imager.dispURL     = "None"
    # Auto window or centered box
    if CleanRad:
        imager.CLEANBox=[-1,CleanRad,0,0]
    else:
        imager.autoWindow  = True
    imager.noScrat     = noScrat
    imager.nThreads    = nThreads
    if debug:
        imager.prtLv = 5
        imager.i
        imager.debug = debug
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
            return 1
        else:
            pass
        # delete Imager file if not debug
        if not debug:
            if doMB:
                u = UV.newPAUV("zap", imager.Sources[0], "MFImage", imager.out2Disk, imager.out2Seq, True, err)
            else:
                u = UV.newPAUV("zap", imager.Sources[0], "Imager", imager.out2Disk, imager.out2Seq, True, err)
            if UV.PIsA(u):
                u.Zap(err) # cleanup
                if err.isErr:
                    mess = "Error deleting Imager work file"
                    printMess(mess, logfile)
                    return 1
            del u
    # end loop over sources

    return 0
    # end EVLAImageTargets

def EVLAAllSource(uv, err, \
               logfile='', check=False, debug=False):
    """ Flag SN table entries with real < Flag

    Returns List of all sources, empty if no SU table
    uv         = UV data object
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
     # Open and close uv to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)

    allSou = []
    if check or debug:
        return allSou
    # Is there an SU table?
    hiver = uv.GetHighVer("AIPS SU")
    if hiver<=0:
        return allSou
    mess = "List of sources in database"  
    printMess(mess, logfile)
    SUTab = uv.NewTable(Table.READONLY, "AIPS SU", 1, err)
    if err.isErr:
        return allSou
    # Open
    SUTab.Open(Table.READONLY, err)
    if err.isErr:
        return allSou
    # Number of rows
    nrow =  SUTab.Desc.Dict["nrow"]
    for i in range (0,nrow):    # Loop over rows
        SUrow = SUTab.ReadRow(i+1, err)
        if err.isErr:
            return
        allSou.append(SUrow["SOURCE"][0].strip())
        mess = "Source("+str(i+1)+") = "+SUrow["SOURCE"][0]  
        printMess(mess, logfile)
    # end loop over rows
            
    # Close table
    SUTab.Close(err)
    if err.isErr:
        return allSou
    return allSou
    # end EVLAAllSource

def EVLAPlotTab(uv, inext, invers, err, \
                source=None, timerang=[0.,0.,0.,0.,0.,0.,0.,0.], \
                stokes="HALF", optype="    ", opcode="    ", nplots=1,  \
                logfile=None, check=False, debug=False):
    """ Makes AIPS plots of tables
    
    Returns task error code, 0=OK, else failed
    uv       = UV data object to plot
    inext    = AIPS Table ("SN", "CL", "TY", "PC")
    inver    = version number, 0-> highest
    source   = if given the name of the source
    timerang = timerange to plot.
    stokes   = Stokes type to plot
    optype   = Data to be plotted (see help snplt)
    opcode   = Plot type (see help snplt)
    nplots   = number of plots per page
    err      = Obit error/message stack
    logfile  = logfile for messages
    check    = Only check script, don't execute tasks
    debug    = show input
    """
    ################################################################
    snplt = AIPSTask.AIPSTask("snplt")
    if not check:
        setname(uv,snplt)
    snplt.inext   = inext
    snplt.invers  = invers
    snplt.optype  = optype
    snplt.opcode  = opcode
    snplt.nplots  = nplots
    snplt.stokes  = stokes
    snplt.msgkill = 5        # Suppress blather
    i = 1
    for t in timerang:
        snplt.timerang[i] = t
        i  += 1
    snplt.logFile = logfile
    if debug:
        snplt.i
    # Trap failure
    try:
        if not check:
            snplt.g
    except Exception, exception:
        print exception
        mess = "SNPLT Failed "
        printMess(mess, logfile)
        return 1
    else:
        pass

    # Open/close UV to update header
    if not check:
        uv.Open(UV.READONLY,err)
        uv.Close(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1
    
    return 0
    # end EVLAPlotTab

def EVLAWritePlots(uv, loPL, hiPL, plotFile, err, \
                     logfile=None, check=False, debug=False):
    """ Writes plots to an external postscript file

    All Plots deleted from AIPS
    Returns task error code, 0=OK, else failed
    uv       = UV data object to plot
    loPL     = Lowest (1-rel) plot number
    hiPL     = Highest PL version number (0->all)
    plotFile = plot file
    err      = Obit error/message stack
    logfile  = logfile for messages
    check    = Only check script, don't execute tasks
    debug    = show input
    """
    ################################################################
    # Open/close UV to update header
    if not check:
        uv.Open(UV.READONLY,err)
        uv.Close(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1
    if hiPL<=0 and not check:
        hiPL = uv.GetHighVer("AIPS PL")

    lwpla = AIPSTask.AIPSTask("lwpla")
    if not check:
        setname(uv,lwpla)
    lwpla.plver   = max(1, loPL)
    lwpla.invers  = hiPL
    lwpla.outfile = plotFile
    lwpla.logFile = logfile
    lwpla.msgkill = 5         # Suppress blather
    if debug:
        lwpla.i
    # Trap failure
    try:
        if not check:
            lwpla.g
    except Exception, exception:
        print exception
        mess = "Lwpla Failed - continuing anyway"
        printMess(mess, logfile)
        # return 1  # Continue in spite of lwpla failure
    else:
        pass
    
    # Delete plot files
    if not check:
        zz=uv.ZapTable("AIPS PL", -1,err)
    
    return 0
    # end EVLAWritePlots

def EVLASpecPlot(uv, Source, timerange, refAnt, err, Stokes=["RR","LL"], \
                 doband=0, plotFile="./spec.ps",
                 check=False, debug=False, logfile = ""):
    """
    Plot amplitude and phase across the spectrum.

    returns scratch file with plot
    Note: possm can't apply flags so data copied to scratch file
    Returns task error code, 0=OK, else failed
    uv        = uv data object
    Source    = Name of source to plot
    timerange = timerange (Obit form) to plot
    refAnt    = ref. Ant, only baselines to this antenna plotted
    err       = Obit error object
    Stokes    = List of stokes types to plot
    doband    = do bandpass calibration before plotting (requires BP table)
    plotFile  = name of output PS file
    check     = Only check script, don't execute tasks
    debug     = Run tasks debug, show input
    logfile   = Log file for task
    """
    # Calibrate and edit data
    scr  = uv.Scratch(err)
    info = uv.List
    info.set("doCalSelect",True)
    info.set("doCalib",2)
    info.set("gainUse",0)
    info.set("doBand",doband)
    info.set("BPVer",0)
    info.set("flagVer",2)
    info.set("Sources",[Source])
    info.set("Stokes","    ")
    info.set("timeRange",timerange)
    #uv.Header(err) # DEBUG 
    uv.Copy(scr, err)
    scr.Info(err)     # Get file information
    
    # Setup and run POSSM
    possm = AIPSTask.AIPSTask("possm")
    setname(scr, possm)
    source = [ Source ] # get BP calibration source, in list format
    possm.sources= AIPSTask.AIPSList( source )
    timerang = [ timerange[0],  0,   0,   0, timerange[1],  0,   0,   0 ]
    possm.timerang = AIPSTask.AIPSList( timerang )
    solint          = (timerange[1]-timerange[0])*1440.0
    possm.baseline[1] = refAnt
    possm.flagver  = -1           # POSSM can't flag
    possm.aparm    = AIPSTask.AIPSList( [0] * 10 ) # initialize with zeroes
    possm.aparm[1] = -1           # scalar average
    possm.aparm[9] = 3            # all IFs and pols in same frame
    possm.nplots   = 6            # plot each baseline in seperate frame on page
    possm.ltype    = 3            # include all labels
    possm.solint   = solint       # time interval of plot
    possm.logFile  = logfile
    possm.msgkill  = 5    # Suppress blather
    # Loop over Stokes
    for s in Stokes:
        possm.stokes   = s
        # Trap failure
        try:
            if not check:
                possm.g
        except Exception, exception:
            print exception
            mess = "POSSM Failed - continue anyway"
            printMess(mess, logfile)
            # return 1
        else:
            pass
        # End Stokes loop
    return scr
# end EVLASpecPlot

def EVLAApplyCal(uv, err, SNver=0, CLin=0, CLout=0, maxInter=240.0, \
                     logfile=None, check=False, debug=False):
    """ Applies an SN table to a CL table and writes another

    Returns task error code, 0=OK, else failed
    uv       = UV data object to clear
    err      = Obit error/message stack
    SNver    = SN table to apply, 0=>highest
    CLin     = input CL table, 0=>highest
    CLout    = output CL table, 0=>create new
    maxInter = Max time (min) over which to interpolate
    logfile  = logfile for messages
    check    = Only check script, don't execute tasks
    debug    = show input, ObitTasks debug
    """
    ################################################################
    # Open/close UV to update header
    if not check:
        uv.Open(UV.READONLY,err)
        uv.Close(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1
    if not check:
        if SNver<=0:
            SNver = uv.GetHighVer("AIPS SN")
        if CLin<=0:
            CLin = uv.GetHighVer("AIPS CL")
        if CLout<=0:
            CLout = uv.GetHighVer("AIPS CL")+1

    if CLin<1:
        mess = "No input CL table to update"
        printMess(mess, logfile)
        uv.Header(err)
        return 1
    mess = "Update CL "+str(CLin)+" with SN "+str(SNver)+" to CL "+str(CLout)
    printMess(mess, logfile)

    clcal = ObitTask.ObitTask("CLCal")
    if not check:
        setname(uv,clcal)
    clcal.solnVer  = SNver
    clcal.calIn    = CLin
    clcal.calOut   = CLout
    clcal.maxInter = maxInter
    clcal.taskLog  = logfile
    clcal.debug    = debug
    if debug:
        clcal.i
    # Trap failure
    try:
        if not check:
            clcal.g
    except Exception, exception:
        print exception
        mess = "CLCal Failed retCode="+str(clcal.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    # End CLCal

    # Open/close UV to update header
    if not check:
        uv.Open(UV.READONLY,err)
        uv.Close(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1
    
    return 0
    # end EVLAApplyCal

def EVLASpectrum(uv, plotSource, plotTime, plotFile, refAnt, err, \
                 Stokes=["RR","LL"], doband=-1,                   \
                 logfile=None, check=False, debug=False):
    """ Spectrum plot of selected data

    Returns task error code, 0=OK, else failed
    uv         = UV data object to clear
    plotSource = Name of source to plot
    plotTime   = timerange (Obit form) to plot
    plotFile   = name of output PS file
    refAnt     = ref. Ant, only baselines to this antenna plotted
    err        = Obit error/message stack
    Stokes     = List of stokes types to plot
    doband     = do bandpass calibration before plotting (requires BP table)
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input, ObitTasks debug
    """
    ################################################################
    # POSSM can't apply flags so write scratch file and plot
    scr = EVLASpecPlot( uv, plotSource,  plotTime, refAnt, err, \
                        Stokes=Stokes, doband=doband,          \
                        plotFile=plotFile, check=check, logfile=logfile )
    if scr.UVIsA():
        retCode = EVLAWritePlots (scr, 1, 0, plotFile, err, \
                                  logfile=logfile, check=check, debug=debug)
    if retCode!=0:
        return 1
    scr.Zap(err)
    return 0
    # end EVLASpectrum
