""" Functions for calibrating and editing VLBA data
"""
import UV, UVDesc, Image, ImageDesc, FArray, ObitTask, AIPSTask, AIPSDir, OErr, History
import InfoList, Table, AIPSDir, OSystem
import os, pickle, math
from AIPS import AIPS
from FITS import FITS
from AIPSDir import AIPSdisks, nAIPS
from OTObit import Acat, AMcat, getname, zap, imhead, tabdest
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
    if err.isErr:   # Ignore if error condition
        return
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
    if logfile and len(logfile) > 0:
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
    tmp = int(ssec*100 + 0.5)  # truncate
    ssec = tmp*0.01
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
   
def VLBAAIPSName( project, session):
    """ Derive AIPS Name

    AIPS file name will be project+session with project truncated to fit in 12 characters
    project    = project name
    session    = session code
    """
    ################################################################
    Aname = Aname=project.strip()[0:12-len(session)]+session
    return Aname
    # end VLBAAIPSName

def VLBAClearCal(uv, err, doGain=True, doBP=False, doFlag=False, logfile=None, check=False):
    """ Clear previous calibration

    Delete all SN tables, all CL but CL 1
    uv       = UV data object to clear
    err      = Obit error/message stack
    doGain   = If True, delete SN and CL tables
    doBP     = If True, delete BP tables
    doFlag   = If True, delete FG tables except FG=1
    logfile  = logfile for messages
    check    = Only check script, don't execute tasks
    """
    ################################################################
    # Only checking?
    if check:
        return
    # Gain tables
    if doGain:
        mess = "Clearing old SN/CL tables"
        printMess(mess, logfile)
        uv.ZapTable("AIPS SN",-1,err)
        ver = uv.GetHighVer("AIPS CL")
        while (ver>1):
            uv.ZapTable ('AIPS CL', ver, err)
            ver = ver-1
    # Bandpass
    if doBP:
        mess = "Clearing old BP tables"
        printMess(mess, logfile)
        uv.ZapTable("AIPS BP",-1,err)
    # Flag tables
    if doFlag:
        mess = "Clearing old FG tables > 1"
        printMess(mess, logfile)
        ver = uv.GetHighVer("AIPS FG")
        while (ver>1):
            uv.ZapTable ('AIPS FG', ver, err)
            ver = ver-1
    OErr.printErrMsg(err, "VLBAClearCal: Error reseting calibration")
    # end VLBAClearCal

def VLBACopyFG(uv, err, logfile='', check=False, debug = False):
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
    taco.logFile = logfile
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
    # end VLBACopyFG

def VLBACopyTable(inObj, outObj, inTab, err, inVer=1, outVer=0,
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
    taco = ObitTask.ObitTask("TabCopy")
    if not check:
        setname(inObj, taco)
        setoname(outObj, taco)
    taco.inTab    = inTab
    taco.inVer    = inVer
    taco.outVer   = outVer
    taco.logFile  = logfile
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
    # end VLBACopyTable

def VLBAIDILoad(filename, project, session, band, Aclass, Adisk, Aseq, err, \
                    wtThresh=0.8,calInt=1.0,logfile='',Compress=False, \
                    check=False, debug=False):
    """ Read FITS IDI file into AIPS

    Read a IDI FITS UV data file and write an AIPS data set
    AIPS file name will be project+session with project truncated to fit in 12 characters
    Returns AIPS data file
    filename   = name of FITS file
    project    = project name
    session    = session code
    band       = observing band code
    Aclass     = AIPS class of file
    Aseq       = AIPS sequence number of file, 0=> create new
    Adisk      = FITS directory number
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
    wtThresh   = Data weight threshold
    calInt     = CL entry interval (min)
    Compress   = Write AIPS data in compressed form?
    check      = Only check script, don't execute tasks
    debug      = show input
    returns AIPS UV data object
    """
    ################################################################
    outuv = None   # In case of failure
    #  Get AIPS Name
    Aname = VLBAAIPSName(project,session)
    fitld = AIPSTask.AIPSTask("fitld")
    fitld.datain   = filename
    fitld.outname  = Aname
    fitld.outclass = Aclass
    fitld.outseq   = Aseq
    fitld.outdisk  = Adisk
    if Compress:
        fitld.douvcomp = 1.0
    fitld.doconcat = 1.0
    fitld.clint    = calInt
    fitld.wtthresh = wtThresh
    fitld.logFile  = logfile
    if debug:
        fitld.i
    # Trap failure
    try:
        if not check:
            fitld.g
    except Exception, exception:
        print exception
        mess = "FITLD load Failed "
        printMess(mess, logfile)
    else:
        pass
    
    # Get output
    if not check:
        outuv = UV.newPAUV("UVdata", Aname, Aclass, Adisk, Aseq, True, err)

    return outuv
    # end VLBAIDILoad

def VLBAFreqInfo(uv, restFreq, srcVel, err, FreqID=1, \
                    logfile='',Compress=False,  check=False, debug=False):
    """ Set lost frequency/velocity infortmation

    Update AIPS SU table, spectrum assumes centered on velocities given
    Returns task error code, 0=OK, else failed
    uv         = UV data object
    restFreq   = list of rest frequencies per IF (Hz)
    srcVel     = list of src name plus velocities (Km/s) and Flux densities (Jy), 
                 eg [("myStar", (0.0,0.0), (1.0,1.0))]  
    err        = Python Obit Error/message stack
    FreqID     = Frequency group identifier
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Set rest frequency/velocity/flux density info into SU table"
    printMess(mess, logfile)

    # Number of channels
    nchan = 1
    setjy=ObitTask.ObitTask("SetJy")
    if not check:
        nchan = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocf"]]
        setname(uv,setjy)
    setjy.Parms[0] = nchan/2+1    # Central channel
    for s in srcVel:              # Loop over source
        nif = len(restFreq)       # Number of IFs
        for iif in range(0,nif):  # Loop over IF
            setjy.Sources[0]  = s[0]
            setjy.BIF         = iif+1
            setjy.EIF         = iif+1
            setjy.RestFreq    = restFreq[iif]
            setjy.SysVel      = s[1][iif]
            setjy.ZeroFlux[0] = s[2][iif]
            setjy.FreqID      = FreqID
            setjy.logFile     = logfile
            setjy.debug       = debug
            if debug:
                setjy.i
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
        # end IF loop
    # end source loop
    return 0
    # end VLBAFreqInfo

def VLBAApplyCal(uv, err, SNver=0, CLin=0, CLout=0, maxInter=240.0, \
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
    # end VLBAApplyCal

def VLBAPlotTab(uv, inext, invers, err, \
                    source=None, timerang=[0.,0.,0.,0.,0.,0.,0.,0.], \
                    optype="    ", opcode="    ", nplots=1,  \
                    logfile=None, check=False, debug=False):
    """ Makes AIPS plots of tables

    Returns task error code, 0=OK, else failed
    uv       = UV data object to plot
    inext    = AIPS Table ("SN", "CL", "TY", "PC")
    inver    = version number, 0-> highest
    source   = if given the name of the source
    timerang = timerange to plot.
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
    # end VLBAPlotTab

def VLBAWritePlots(uv, loPL, hiPL, plotFile, err, \
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
    # end VLBAWritePlots

def VLBAQuantCor(uv, QuantSmoo, QuantFlag, err, \
                     doSNPlot=False, plotFile="./Quant.ps", \
                     logfile='', check=False, debug=False):
    """ Determine/apply quantization correstion, possibly with flagging

    Determine quantization corrections using AIPS/ACCOR to new SN table
    Smooth derived SN table (Obit/SNSmo) to QuamntSmoo
    If QuantFlag, then the SN table is flagged with values below QuantFlag
    Apply this SN table to the highest CL table writing a new CL table (Obit/CLCal)
    Returns task error code, 0=OK, else failed
    uv         = UV data object
    QuantSmoo  = Smoothing time in Hr
    QuantFlag  = if >0.0 flag SN table entries with real<QuantFlag
    err        = Python Obit Error/message stack
    doSNPlot   = If True make plots of SN gains
    plotFile   = Name of postscript file for plots
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine/apply quantization corrections"
    printMess(mess, logfile)

    accor = AIPSTask.AIPSTask("accor")
    if not check:
        setname(uv,accor)
    accor.logFile  = logfile
    if debug:
        accor.i
    # Trap failure
    try:
        if not check:
            accor.g
    except Exception, exception:
        print exception
        mess = "ACCOR Failed "
        printMess(mess, logfile)
        return 1
    else:
        pass

    # Open/close UV to update header
    if not check:
        uv.Open(UV.READONLY,err)
        uv.Close(err)
    #uv.UpdateTables(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1

    # Smooth
    SNver = 0
    if not check:
        SNver = uv.GetHighVer("AIPS SN")
    snsmo = ObitTask.ObitTask("SNSmo")
    if not check:
        setname(uv, snsmo)
    snsmo.solnIn     = SNver
    snsmo.solnOut    = SNver+1
    snsmo.smoType    = 'AMPL'
    snsmo.smoParm[0] = QuantSmoo
    snsmo.taskLog    = logfile
    snsmo.debug      = debug
    if debug:
        snsmo.i
    mess = "VLBAQuantCor: SNSmo SN"+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
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
    
    # Flagging if requested
    if QuantFlag>0.0:
        # Open/close UV to update header
        uv.Open(UV.READONLY,err)
        uv.Close(err)
        if err.isErr:
            OErr.printErr(err)
            mess = "Update UV header failed"
            printMess(mess, logfile)
            return 1
        #print "DEBUG version",snsmo.solnOut
        VLBAClipSN (uv, snsmo.solnOut, QuantFlag, err, \
                    logfile=logfile, check=check, debug=debug)
    if err.isErr:
        OErr.printErr(err)
        mess = "Flagging of ACCOR soln failed"
        printMess(mess, logfile)
        return 1
    # end flagging

    # Plot opacity fits?
    if doSNPlot:
        retCode = VLBAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="REAL", \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
  
        retCode = VLBAWritePlots (uv, 1, 0, plotFile, err, \
                                      logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        
    # end SN table plot
    # Apply to CL table
    retCode = VLBAApplyCal(uv, err, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode
    return 0
    # end VLBAQuantCor

def VLBAPACor(uv, err, CLver=0, FreqID=1,\
                  noScrat=[], logfile='', check=False, debug=False):
    """ Make parallactic angle correction

    Updates CL CLver
    Returns task error code, 0=OK, else failed
    uv         = UV data object
    err        = Python Obit Error/message stack
    CLver      = Cl version to update, 0=> highest
    FreqID     = Frequency group identifier
    noScrat  = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Make parallactic angle corrections"
    printMess(mess, logfile)

    # Which CL?
    if CLver<=0 and not check:
        CLver = uv.GetHighVer("AIPS CL")
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
    # end VLBAPACor

def VLBAOpacCor(uv, OpacSmoo, err, FreqID=1, WXver=0, TYver=0, GCver=0, \
                    doSNPlot=False, plotFile="./Opacity.ps", \
                    logfile='', check=False, debug=False):
    """ Make Opacity/gain correction from TY/GC table

    Creates new SN table, smooths, applies to highest numbered Cl 
    and writes new CL
    Returns task error code, 0=OK, else failed
    uv         = UV data object
    OpacSmoo   = Smoothing time (hr) for SN table
    err        = Python Obit Error/message stack
    FreqID     = Frequency group identifier
    TYver      = TY version to use, 0=> highest
    GCver      = GC version to use, 0=> highest
    WXver      = WX version to use, 0=> highest, -1=>none
    doSNPlot   = If True make plots of SN gains
    plotFile   = Name of postscript file for plots
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine/apply opaticity/gin/Tsys corrections"
    printMess(mess, logfile)

    # Open/close UV to update header
    if not check:
        uv.Open(UV.READONLY,err)
        uv.Close(err)
    #uv.UpdateTables(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1
    # Input tables
    if not check:
        if WXver==0:
            WXver = uv.GetHighVer("AIPS WX")
        if TYver==0:
            TYver = uv.GetHighVer("AIPS TY")
        if GCver==0:
            GCver = uv.GetHighVer("AIPS GC")

    apcal = AIPSTask.AIPSTask("apcal")
    if not check:
        setname(uv,apcal)
    apcal.opcode      = "GRID"
    apcal.inver       = WXver
    apcal.tyver       = TYver
    apcal.gcver       = GCver
    apcal.aparm[5]    = 4.0     # Max airmass to include in fit
    apcal.freqid      = FreqID
    apcal.ltype       = 3       # 
    apcal.logFile     = logfile
    for i in range(1,len(apcal.dofit)):  # Fit all antennas
        apcal.dofit[i] = 1.0
    if debug:
        apcal.i
    # Trap failure
    try:
        if not check:
            apcal.g
    except Exception, exception:
        print exception
        mess = "APCAL Failed "
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

    # Smooth
    SNver = 0
    if not check:
        SNver = uv.GetHighVer("AIPS SN")
    snsmo = ObitTask.ObitTask("SNSmo")
    if not check:
        setname(uv, snsmo)
    snsmo.solnIn     = SNver
    snsmo.solnOut    = SNver+1
    snsmo.smoType    = 'AMPL'
    snsmo.smoParm[0] = OpacSmoo
    snsmo.taskLog    = logfile
    snsmo.debug      = debug
    if debug:
        snsmo.i
    mess = "VLBAOpacCor: SNSmo SN"+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
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
    
    # Plot opacity fits?
    if doSNPlot:
        retCode = VLBAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="REAL", \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
  
        retCode = VLBAWritePlots (uv, 1, 0, plotFile, err, \
                                      logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        
    # end SN table plot
        

    # Delete any plot files left
    if not check: 
        zz=uv.ZapTable("AIPS PL", -1,err)

    # Apply latest SN to CL table
    retCode = VLBAApplyCal(uv, err, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode

    return 0
    # end VLBAOpacCor

def VLBAGoodCal(uv, err, solInt=0.5, timeInt=100., FreqID=1, \
                calSou=None, CalModel=None, \
                doCalib=-1, gainUse=0, minSNR = 10.0, refAnts=[0], \
                doBand=-1, BPVer=0,  \
                flagVer=-1,  \
                nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Find best calibration source and timerange and best refAnt

    Determines good calibrator source and timerange and best
    reference antenna by running a phase Calib and looking at the
    SNRs of the solutions.
    Return dict {"Source":source, "souID":souID, "timeRange":(tbeg,tend), \
                "Fract":fract_OK, "SNR":avg_SNR,"bestRef":refAnt}
                or None on failure
    uv         = UV data object
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
    solInt     = Calib solution interval (min)
    timeInt    = interval width (min) to use for timerange (large=>scan)
    FreqID     = Frequency group identifier
    minSNR     = minimum acceptable SNR in Calib
    refAnts    = List of reference antennas in priority order,
                 Also, if values given then list of acceptable ref. ants
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Find best calibrator source, timerange, refAnt"
    printMess(mess, logfile)

    # Set output (new) SN table
    SNver = 0
    if not check:
        SNver = uv.GetHighVer("AIPS SN")+1

    calib = ObitTask.ObitTask("Calib")
    calib.taskLog  = logfile
    if not check:
        setname(uv,calib)
    calib.flagVer   = flagVer
    calib.doCalib   = doCalib
    calib.gainUse   = gainUse
    calib.doBand    = doBand
    calib.BPVer     = BPVer
    calib.solMode   = "DELA"
    calib.solType   = "  "
    calib.solInt    = solInt
    calib.refAnts   = refAnts
    calib.minSNR    = minSNR
    calib.solnVer   = SNver
    calib.nThreads  = nThreads
    calib.noScrat   = noScrat
    #  If multiple sources given with models - loop
    if (type(calSou)==list) and (len(calSou)>1) and CalModel:
        ncloop = len(calSou)
        calsou = calSou[0] # setup for first
    else:
        ncloop = 1;
        calsou = calSou
    # Loop
    for ical in range(0,ncloop):
        if type(calsou)==list:
            calib.Sources = calsou
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
            calib.prtLv = 2
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
    # no further for check only
    if check:
        return

    # Open/close UV to update header
    uv.Open(UV.READONLY,err)
    uv.Close(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return None

    # Get SN table version
    SNver = uv.GetHighVer("AIPS SN")

    # Get statistics
    out = VLBASNStats(uv, SNver, timeInt, err, refAnts=refAnts, \
                      logfile=logfile, check=check, debug=debug)

    # Tell results
    mess = "Best Calibrator is "+out["Source"]+" during "+\
        day2dhms(+out["timeRange"][0])+" - "+day2dhms(+out["timeRange"][1])
    printMess(mess, logfile)
    mess = "Best reference antenna "+str(out["bestRef"])+" Fract "+\
        str(out["Fract"])+" Avg SNR "+str(out["SNR"])
    printMess(mess, logfile)
 
    return out
    # end VLBAGoodCal

def VLBAPCcor(uv, err, calSou=None,  timeRange=[0.,0.], FreqID=1, \
                  doCalib=-1, gainUse=0, doBand=-1, BPVer=0,  \
                  flagVer=-1,  solInt = 0.5, \
                  PCin=1, SNout=0, refAnt=1, doPCPlot=False, plotFile="./PC.ps", \
                  noScrat=[], logfile='', check=False, debug=False):
    """ Apply corrections from PC table

    Use a short section of calibrator data to resolve ambiguities in tghe PC
    measurements and apply the calibration derived from a PC file and write an
    SN table.  This SN table is then applied to the previous highest CL table
    creating a new CL table.
    Returns task error code, 0=OK, else failed
    err        = Python Obit Error/message stack
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    calSou     = Source name to use
    timeRange  = timerange of data to use
    FreqID     = Frequency group identifier
    PCin       = Input PC table input number
    SNout      = Output SN table, 0=> create new
    refAnt     = Reference antenna
    solInt     = solution interval for cal data (min)
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine Phase cal (PC) corrections"
    printMess(mess, logfile)

    # Set output (new) SN table
    if SNout<=0:
        SNver = uv.GetHighVer("AIPS SN")+1
    else:
        SNver = SNout

    pccor = ObitTask.ObitTask("PCCor")
    if not check:
        setname(uv,pccor)
    pccor.flagVer      = flagVer
    pccor.doCalib      = doCalib
    pccor.gainUse      = gainUse
    pccor.doBand       = doBand
    pccor.BPVer        = BPVer
    pccor.solnVer      = SNver
    pccor.PCVer        = PCin
    pccor.refAnt       = refAnt
    pccor.calSour      = calSou
    pccor.timeRange[0] = timeRange[0]
    pccor.timeRange[1] = timeRange[1]
    pccor.solInt       = solInt
    pccor.prtLv        = 1
    pccor.taskLog      = logfile
    pccor.noScrat      = noScrat
    pccor.debug        = debug
    pccor.debug = True
    if debug:
        pccor.i
    # Trap failure
    try:
        if not check:
            pccor.g
    except Exception, exception:
        print exception
        mess = "PCCor Failed "
        printMess(mess, logfile)
        return 1
    else:
        pass

    # Open/close UV to update header
    uv.Open(UV.READONLY,err)
    uv.Close(err)
    if err.isErr:
        OErr.printErr(err)
        mess = "Update UV header failed"
        printMess(mess, logfile)
        return 1

    # Plot
    SNver = uv.GetHighVer("AIPS SN")
    # Plot PC corrections?
    if doPCPlot:
        # Phase corrections
        retCode = VLBAPlotTab(uv, "SN", SNver, err, nplots=6, optype="PHAS", \
                              logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        # Delay corrections
        retCode = VLBAPlotTab(uv, "SN", SNver, err, nplots=6, optype="DELA", \
                              logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
  
        retCode = VLBAWritePlots (uv, 1, 0, plotFile, err, \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
    # end SN table plot

    # Apply to CL table - one scan for whole dataset
    retCode = VLBAApplyCal(uv, err, maxInter=2880.0, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode
    return 0
    # end VLBAPCcor

def VLBAManPCal(uv, err, solInt=0.5, smoTime=10.0, calSou=None, CalModel=None, 
                timeRange=[0.,0.], FreqID=1, doCalib=-1, gainUse=0, 
                minSNR = 10.0, refAnts=[0], doBand=-1, BPVer=0, flagVer=-1,  
                doManPCalPlot=True, plotFile='ManPCal.ps', 
                nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Manual phase cal

    Determine phase corrections from a short section of strong source data.
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
    smoTime    = Smoothing time applied to SN table (min)
    FreqID     = Frequency group identifier
    minSNR     = minimum acceptable SNR in Calib
    refAnt   = Reference antenna
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine manual phase corrections"
    printMess(mess, logfile)

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
    #  If multiple sources given with models - loop
    if (type(calSou)==list) and (len(calSou)>1) and CalModel:
        ncloop = len(calSou)
        calsou = calSou[0] # setup for first
    else:
        ncloop = 1;
        calsou = calSou
    # Loop
    for ical in range(0,ncloop):
        if type(calsou)==list:
            calib.Sources = calsou
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
            calib.prtLv = 2
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

    # Smooth
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
    snsmo.taskLog    = logfile
    snsmo.debug      = debug
    if debug:
        snsmo.i
    mess = "VLBAManPCal: SNSmo SN"+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
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

    # Filter to remove phase fluctuations, average delays
    #?? VLBARefMB(uv, SNver+1, err, logfile=logfile, check=check,debug=debug)
    #?? if err.isErr:
    #??     mess = "VLBARefMB Failed"
    #??     printMess(mess, logfile)
    #??     return 1
    
    # Plot man p-cal corrections?
    if doManPCalPlot:
        # Phase corrections
        retCode = VLBAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="PHAS", \
                              logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        # Delay corrections
        retCode = VLBAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="DELA", \
                              logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        retCode = VLBAWritePlots (uv, 1, 0, plotFile, err, \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
    # end SN table plot

    # Apply to CL table - one scan for whole dataset
    retCode = VLBAApplyCal(uv, err, maxInter=2880.0, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode
    return 0
    # end VLBAManPCal

def VLBABPass(uv, BPCal, err, CalModel=None, newBPVer=1, timeRange=[0.,0.], \
                  doCalib=2, gainUse=0, doBand=0, BPVer=0, flagVer=-1,  \
                  doCenter1=None, BChan1=1, EChan1=0, \
                  BChan2=1, EChan2=0, ChWid2=1, \
                  solInt1=0.0, solInt2=0.0, solMode="A&P", refAnt=0, \
                  ampScalar=False, specIndex=0.0, \
                  doAuto=False, doPol=False, avgPol=False, avgIF=False, \
                  check=False, debug = False, \
                  nThreads=1, noScrat=[], logfile = ""):
    """ Bandpass calibration

    Determine bandpass from a short section of strong source data.
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
    Note: this makes no Doppler corrections - may have limited accuracy
    Returns task error code, 0=OK, else failed
    uv       = UV data object to calibrate
    BPCal    = Bandpass calibrator, name or list of names
    err      = Obit error/message stack
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
    newBPVer = output BP table, should be set for multiple calibrators w/ model
    doCalib  = Apply calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    doBand   = If >0.5 apply previous bandpass cal.
    BPVer    = previous Bandpass table (BP) version
    flagVer  = Input Flagging table version
    timerange= timerange in days to use
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
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    nThreads = Number of threads to use
    noScrat  = list of AIPS disks to avoid for scratch files
    logfile  = Log file for task
    """
    ################################################################
    mess = "Determine bandpass corrections"
    printMess(mess, logfile)

    bpass = ObitTask.ObitTask("BPass")
    bpass.taskLog = logfile
    if not check:
        setname(uv,bpass)
    bpass.doBand    = doBand
    bpass.BPVer     = BPVer
    bpass.BPSoln    = newBPVer
    bpass.doCalib   = doCalib
    bpass.gainUse   = gainUse
    bpass.flagVer   = flagVer
    bpass.doPol     = doPol
    bpass.BChan1    = BChan1
    bpass.EChan1    = EChan1
    bpass.BChan2    = BChan2
    bpass.EChan2    = EChan2
    bpass.ChWid2    = ChWid2
    bpass.solInt1   = solInt1
    bpass.solInt2   = solInt2
    bpass.solMode   = solMode
    bpass.Alpha     = specIndex
    bpass.refAnt    = refAnt
    bpass.timeRange = timeRange
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
    
    #  If multiple sources given with models - loop
    if (type(BPCal)==list) and (len(BPCal)>1) and CalModel:
        ncloop = len(BPCal)
        calsou = BPCal[0] # setup for first
    else:
        ncloop = 1;
        calsou = BPCal

    # Loop
    for ical in range(0,ncloop):
        if type(calsou)==list:
            bpass.Sources = calsou
        else:
            bpass.Sources = [calsou]
        # Get model details
        if CalModel and bpass.Sources[0] in CalModel:
            Model = CalModel[bpass.Sources[0]]
            if Model["image"]:
                if not check:
                    set2name(Model["image"], bpass)
            if "nfield" in Model:
                bpass.nfield    = Model["nfield"]
            if "CCVer" in Model:
                bpass.CCVer     = Model["CCVer"]
            if "BComp" in Model:
                bpass.BComp     = Model["BComp"]
            if "EComp" in Model:
                bpass.EComp     = Model["EComp"]
            if "Cmethod" in Model:
                bpass.Cmethod   = Model["Cmethod"]
            if "Cmodel" in Model:
                bpass.Cmodel    = Model["Cmodel"]
            if "Flux" in Model:
                bpass.Flux      = Model["Flux"]
            if "modelFlux" in Model:
                bpass.modelFlux = Model["modelFlux"]
            if "modelPos" in Model:
                bpass.modelPos  = Model["modelPos"]
            if "modelParm" in Model:
                bpass.modelParm = Model["modelParm"]

        if debug:
            bpass.i
            bpass.prtLv = 5
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
        # Setup for next if looping
        if ical<(ncloop-1):
            calsou = BPCal[ical+1]
    # End calibrator loop
    return 0
# end VLBABPass

def VLBASpecPlot(uv, goodCal, err, doband=0, plotFile="./spec.ps", logfile = ""):
    """
    Plot amplitude and phase across the spectrum.

    uv = uv data object
    goodCal = good calibrator info; dictionary from VLBAGoodCal
    err = Obit error object
    doband = do bandpass calibration before plotting (requires BP table)
    plotFile = name of output PS file
    logfile  = Log file for task
    """
    # Remove any pre-existing PL tables
    tabdest(uv, "AIPS PL", -1)

    # Setup and run POSSM
    possm = AIPSTask.AIPSTask("possm")
    setname(uv, possm)
    source = [ goodCal["Source"] ] # get BP calibration source, in list format
    possm.sources= AIPSTask.AIPSList( source )
    tr = goodCal["timeRange"] # get BP calibration time range
    #          [   day, hr, min, sec,   day, hr, min, sec ] <- possm format
    timerang = [ tr[0],  0,   0,   0, tr[1],  0,   0,   0 ]
    possm.timerang = AIPSTask.AIPSList( timerang )
    possm.stokes = 'HALF' # plot only RR & LL
    possm.docalib = 1 # calibrate data & weights
    possm.doband = doband # calibrate bandpass
    possm.bpver = 0 # apply highest BP table
    possm.aparm = AIPSTask.AIPSList( [0] * 10 ) # initialize with zeroes
    possm.aparm[9] = 3 # all IFs and pols in same frame
    possm.nplots = 2 # plot each baseline in seperate frame on page
    possm.ltype = 3 # include all labels
    possm.logFile = logfile
    possm.go()

    # Setup and run LWPLA
    lwpla = AIPSTask.AIPSTask("lwpla")
    setname(uv,lwpla)
    lwpla.plver = 1
    print "PL high ver = ", uv.GetHighVer("AIPS PL")
    uv.Header(err) # this is required in order for GetHighVer to work
    print "PL high ver = ", uv.GetHighVer("AIPS PL")
    lwpla.invers = uv.GetHighVer("AIPS PL")
    lwpla.outfile = plotFile
    lwpla.logFile = logfile
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

    # Remove any tables
    tabdest(uv, "AIPS PL", -1)

# end VLBASpecPlot

def VLBAImageCals(uv, err,  FreqID=1, Sources=None, seq=1, sclass="ImgSC", \
                  doCalib=-1, gainUse=0, doBand=-1, BPVer=0,  flagVer=-1,  \
                  FOV=0.1/3600.0, Robust=0, Niter=300, \
                  maxPSCLoop=6, minFluxPSC=0.1, solPInt=20.0/60., solMode="P", \
                  maxASCLoop=1, minFluxASC=0.5, solAInt=2.0, \
                  avgPol=False, avgIF=False, minSNR = 5.0, refAnt=0, \
                  nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Image a list of sources

    Uses SCMap to image a list of sources.
    Returns task error code, 0=OK, else failed
    
    uv         = UV data object
    err        = Python Obit Error/message stack
    Sources    = Source name or list of names to use
                 If an empty list all sources in uv are included
    seq        = sequence number of output
    sclass     = Image output class
    FreqID     = Frequency group identifier
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    FOV        = Field of view to image in deg
    Robust     = Weighting robustness parameter
    Niter      = max no. iterations
    maxPSCLoop = max. number of phase sc loops
    minFluxPSC = Trip level for phase self cal (Jy)
    solPInt    = Phase solution interval (min)
    solMode    = Phase soln mode "P", "DELA"
    maxASCLoop = max. number of amp&phase sc loops
    minFluxASC = Trip level for amp&phase self cal (Jy)
    solAInt    = Amp&phase solution interval (min)
    avgPol     = Average poln in SC?
    avgIF      = Average IFs in SC?
    minSNR     = minimum acceptable SNR in SC
    refAnt     = Reference antenna
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Image a list of sources "
    printMess(mess, logfile)

    # If list empty get all sources
    if type(Sources)==list:
        sl = Sources
    else:
        sl = [Sources]

    if len(sl)<=0:
        slist = VLBAAllSource(uv,err,logfile=logfile,check=check,debug=debug)
    else:
        slist = sl
    
    scmap = ObitTask.ObitTask("SCMap")
    scmap.taskLog  = logfile
    if not check:
        setname(uv,scmap)
    scmap.outDisk     = scmap.inDisk
    scmap.out2Disk    = scmap.inDisk
    scmap.outSeq      = seq
    scmap.out2Seq     = seq
    scmap.outClass    = sclass
    scmap.BLFact      = 1.004
    scmap.BLchAvg     = True
    scmap.flagVer     = flagVer
    scmap.doCalib     = doCalib
    scmap.gainUse     = gainUse
    scmap.doBand      = doBand
    scmap.BPVer       = BPVer
    scmap.FOV         = FOV
    scmap.Robust      = Robust
    scmap.Niter       = Niter
    scmap.maxPSCLoop  = maxPSCLoop
    scmap.minFluxPSC  = minFluxPSC
    scmap.solPInt     = solPInt
    scmap.solPMode    = solMode
    scmap.maxASCLoop  = maxASCLoop
    scmap.minFluxASC  = minFluxASC
    scmap.solAInt     = solAInt
    scmap.avgPol      = avgPol
    scmap.avgIF       = avgIF
    scmap.refAnt      = refAnt
    scmap.minSNR      = minSNR
    scmap.dispURL     = "None"
    scmap.autoWindow  = True
    scmap.modelFlux   = 1.0
    scmap.noScrat     = noScrat
    scmap.nThreads    = nThreads
    scmap.prtLv       = 3
    if debug:
        scmap.prtLv = 5
        scmap.i
        scmap.debug = debug
    # Loop over slist
    for sou in slist:
        scmap.Sources[0] = sou
        mess = "Image/selfcal "+sou
        printMess(mess, logfile)
        # Trap failure
        try:
            if not check:
                scmap.g
        except Exception, exception:
            print exception
            mess = "SCMap Failed retCode= "+str(scmap.retCode)
            printMess(mess, logfile)
            return 1
        else:
            pass
        # Save SN table and delete SCMap file if not debug
        if not debug:
            u = UV.newPAUV("zap", scmap.Sources[0], "SCMap", scmap.out2Disk, scmap.out2Seq, True, err)
            if UV.PIsA(u):
                u.Zap(err) # cleanup
                if err.isErr:
                    mess = "Error deleting SCMap work file"
                    printMess(mess, logfile)
                    return 1
            del u
    # end loop over sources

    return 0
    # end VLBAImageCals

def VLBAImageTargets(uv, err,  FreqID=1, Sources=None, seq=1, sclass="IClean", \
                         doCalib=-1, gainUse=0, doBand=-1, BPVer=0,  flagVer=-1,  \
                         Stokes="I", FOV=0.1/3600.0, Robust=0, Niter=300, \
                         maxPSCLoop=0, minFluxPSC=0.1, solPInt=20.0/60., solMode="P", \
                         maxASCLoop=0, minFluxASC=0.5, solAInt=2.0, \
                         avgPol=False, avgIF=False, minSNR = 5.0, refAnt=0, \
                         nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Image a list of sources with optional selfcal

    Uses Imager to image a list of sources.
    Data must be at least approximately calibrated
    Returns task error code, 0=OK, else failed
    
    uv         = UV data object
    err        = Python Obit Error/message stack
    Sources    = Source name or list of names to use
                 If an empty list all sources in uv are included
    seq        = sequence number of output
    sclass     = Image output class
    FreqID     = Frequency group identifier
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    Stokes     = Stokes parameters to image
    FOV        = Field of view to image in deg
    Robust     = Weighting robustness parameter
    Niter      = max no. iterations
    maxPSCLoop = max. number of phase sc loops
    minFluxPSC = Trip level for phase self cal (Jy)
    solPInt    = Phase solution interval (min)
    solMode    = Phase soln mode "P", "DELA"
    maxASCLoop = max. number of amp&phase sc loops
    minFluxASC = Trip level for amp&phase self cal (Jy)
    solAInt    = Amp&phase solution interval (min)
    avgPol     = Average poln in SC?
    avgIF      = Average IFs in SC?
    minSNR     = minimum acceptable SNR in SC
    refAnt     = Reference antenna
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Image a list of sources "
    printMess(mess, logfile)

    # If list empty get all sources
    if type(Sources)==list:
        sl = Sources
    else:
        sl = [Sources]

    if len(sl)<=0:
        slist = VLBAAllSource(uv,err,logfile=logfile,check=check,debug=debug)
    else:
        slist = sl
    imager = ObitTask.ObitTask("Imager")
    imager.taskLog  = logfile
    if not check:
        setname(uv,imager)
    imager.outDisk     = imager.inDisk
    imager.out2Disk    = imager.inDisk
    imager.outSeq      = seq
    imager.out2Seq     = seq
    imager.outClass    = sclass
    imager.BLFact      = 1.004
    imager.BLchAvg     = True
    imager.flagVer     = flagVer
    imager.doCalib     = doCalib
    imager.gainUse     = gainUse
    imager.doBand      = doBand
    imager.BPVer       = BPVer
    imager.Stokes      = Stokes
    imager.FOV         = FOV
    imager.Robust      = Robust
    imager.Niter       = Niter
    imager.maxPSCLoop  = maxPSCLoop
    imager.minFluxPSC  = minFluxPSC
    imager.solPInt     = solPInt
    imager.solPMode    = solMode
    imager.maxASCLoop  = maxASCLoop
    imager.minFluxASC  = minFluxASC
    imager.solAInt     = solAInt
    imager.avgPol      = avgPol
    imager.avgIF       = avgIF
    imager.refAnt      = refAnt
    imager.minSNR      = minSNR
    imager.dispURL     = "None"
    imager.autoWindow  = True
    imager.noScrat     = noScrat
    imager.nThreads    = nThreads
    imager.prtLv       = 3
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
    # end VLBAImageTargets

def VLBADelayCal(uv, err, solInt=0.5, smoTime=10.0, calSou=None,  CalModel=None, \
                     timeRange=[0.,0.], FreqID=1, doCalib=-1, gainUse=0, minSNR=5.0, \
                     refAnts=[0], doBand=-1, BPVer=0, flagVer=-1,  \
                     doSNPlot=False, plotFile="./DelayCal.ps", \
                     nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Group delay calibration

    Determine delay corrections from a list of calibrators
    Solutions smoothed to smoTime
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
    smoTime    = Smoothing time applied to SN table (min)
    FreqID     = Frequency group identifier
    minSNR     = minimum acceptable SNR in Calib
    refAnts    = List of reference antennas
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    doSNPlot   = If True make plots of SN gains
    plotFile   = Name of postscript file for plots
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine group delays"
    printMess(mess, logfile)

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
    #  If multiple sources given with models - loop
    if (type(calSou)==list) and (len(calSou)>1) and CalModel:
        ncloop = len(calSou)
        calsou = calSou[0] # setup for first
    else:
        ncloop = 1;
        calsou = calSou
    # Loop
    for ical in range(0,ncloop):
        if type(calsou)==list:
            calib.Sources = calsou
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
            calib.prtLv = 2
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

    # Smooth
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
    if debug:
        snsmo.i
    mess = "VLBADelayCal: SNSmo SN "+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
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

    # Filter to remove phase fluctuations, average delays
    VLBARefMB(uv, SNver+1, err, smoTime=smoTime,
              logfile=logfile, check=check,debug=debug)
    if err.isErr:
        mess = "VLBARefMB Failed"
        printMess(mess, logfile)
        return 1
    
    # Plot fits?
    if doSNPlot:
        retCode = VLBAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="DELA", \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
  
        retCode = VLBAWritePlots (uv, 1, 0, plotFile, err, \
                                      logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        
    # end SN table plot
    # Apply to CL table
    retCode = VLBAApplyCal(uv, err, maxInter=600.0, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode
    return 0
    # end VLBADelayCal

def VLBAAmpCal(uv, err, solInt=0.5, smoTimeA=1440.0, smoTimeP=10.0, \
               calSou=None,  CalModel=None, \
               timeRange=[0.,0.], FreqID=1, doCalib=-1, gainUse=0, minSNR=5.0, \
               refAnt=0, doBand=-1, BPVer=0, flagVer=-1,  \
               doSNPlot=False, plotFile="./AmpCal.ps", \
               nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Amplitude calibration

    Determine amplitude calibratiion from a list of calibrators
    Solutions smoothed to smoTime
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
    smoTimeA   = Amplitude Smoothing time applied to SN table (min)
    smoTimeP   = Phase Smoothing time applied to SN table (min)
    FreqID     = Frequency group identifier
    minSNR     = minimum acceptable SNR in Calib
    refAnt     = Reference antenna
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    doSNPlot   = If True make plots of SN gains
    plotFile   = Name of postscript file for plots
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine amplitude calibration"
    printMess(mess, logfile)

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
    calib.solMode   = "A&P"
    calib.solType   = "  "
    calib.solInt    = solInt
    calib.minSNR    = minSNR
    calib.refAnts   = [refAnt]
    calib.prtLv     = 3
    calib.solnVer   = SNver
    calib.nThreads  = nThreads
    calib.noScrat   = noScrat
    #  If multiple sources given with models - loop
    if (type(calSou)==list) and (len(calSou)>1) and CalModel:
        ncloop = len(calSou)
        calsou = calSou[0] # setup for first
    else:
        ncloop = 1;
        calsou = calSou
    # Loop
    for ical in range(0,ncloop):
        if type(calsou)==list:
            calib.Sources = calsou
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
            calib.nfield    = max (1, calib.nfield)
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
            calib.prtLv = 2
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

    # Smooth
    SNver = uv.GetHighVer("AIPS SN")
    snsmo = ObitTask.ObitTask("SNSmo")
    if not check:
        setname(uv, snsmo)
    snsmo.solnIn     = SNver
    snsmo.solnOut    = SNver+1
    snsmo.smoType    = 'AMPL'
    snsmo.smoFunc    = 'MWF '
    snsmo.smoParm    = [smoTimeA,smoTimeP]
    snsmo.clipSmo    = [smoTimeA,smoTimeP]
    snsmo.clipParm   = [2.0]   # May be too conservative
    snsmo.doBlank    = True
    snsmo.refAnt     = refAnt
    snsmo.taskLog    = logfile
    snsmo.debug      = debug
    if debug:
        snsmo.i
    mess = "VLBAAmpCal: SNSmo SN "+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
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
    if doSNPlot:
        retCode = VLBAPlotTab(uv, "SN", SNver+1, err, nplots=6, optype="AMP", \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
  
        retCode = VLBAWritePlots (uv, 1, 0, plotFile, err, \
                                      logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        
    # Apply to CL table
    retCode = VLBAApplyCal(uv, err, maxInter=600.0, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode
    return 0
    # end VLBAAmpCal

def VLBAPhaseCal(uv, err, solInt=0.5, smoTimeA=1440.0, smoTimeP=10.0/60.0, doSmooth=False, \
               calSou=None,  CalModel=None, \
               timeRange=[0.,0.], FreqID=1, doCalib=-1, gainUse=0, minSNR=5.0, \
               refAnt=0, doBand=-1, BPVer=0, flagVer=-1,  \
               doSNPlot=False, plotFile="./PhaseCal.ps", \
               nThreads=1, noScrat=[], logfile='', check=False, debug=False):
    """ Phase calibration

    Determine phase calibration from a list of calibrators
    Solutions smoothed to smoTime
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
    doSmooth   = If true Smooth table
    smoTimeA   = Amplitude Smoothing time applied to SN table (min)
    smoTimeP   = Phase Smoothing time applied to SN table (min)
    FreqID     = Frequency group identifier
    minSNR     = minimum acceptable SNR in Calib
    refAnt     = Reference antenna
    doCalib    = Apply calibration table
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    doSNPlot   = If True make plots of SN gains
    plotFile   = Name of postscript file for plots
    nThreads   = Max. number of threads to use
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Determine phase calibration"
    printMess(mess, logfile)

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
    calib.solMode   = "P"
    calib.solType   = "L1"
    calib.solInt    = solInt
    calib.minSNR    = minSNR
    calib.refAnts   = [refAnt]
    calib.prtLv     = 3
    calib.solnVer   = SNver
    calib.nThreads  = nThreads
    calib.noScrat   = noScrat
    #  If multiple sources given with models - loop
    if (type(calSou)==list) and (len(calSou)>1) and CalModel:
        ncloop = len(calSou)
        calsou = calSou[0] # setup for first
    else:
        ncloop = 1;
        calsou = calSou
    # Loop
    for ical in range(0,ncloop):
        if type(calsou)==list:
            calib.Sources = calsou
        else:
            calib.Sources = [calsou]
        mess = "Calibrate "+ str(calsou)
        printMess(mess, logfile)

        # Get model details
        if CalModel and calib.Sources[0] in CalModel:
            Model = CalModel[calib.Sources[0]]
            if Model["image"]:
                # Make sure there is a model
                CChi = Model["image"].GetHighVer("AIPS CC")
                if CChi<=0:
                    mess = "No model - skipping "
                    printMess(mess, logfile)
                    # Setup for next if looping
                    if ical<(ncloop-1):
                        calsou = calSou[ical+1]
                    continue
                if not check:
                    set2name(Model["image"], calib)
            if "nfield" in Model:
                calib.nfield    = Model["nfield"]
            calib.nfield    = max (1, calib.nfield)
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
            calib.prtLv = 2
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

    # Smooth
    SNver = uv.GetHighVer("AIPS SN")
    if doSmooth:
        snsmo = ObitTask.ObitTask("SNSmo")
        if not check:
            setname(uv, snsmo)
        snsmo.solnIn     = SNver
        snsmo.solnOut    = SNver+1
        snsmo.smoType    = 'AMPL'
        snsmo.smoFunc    = 'MWF '
        snsmo.smoParm    = [smoTimeA,smoTimeP]
        snsmo.clipSmo    = [smoTimeA,smoTimeP]
        snsmo.doBlank    = True
        snsmo.refAnt     = refAnt
        snsmo.taskLog    = logfile
        snsmo.debug      = debug
        if debug:
            snsmo.i
        mess = "VLBAPhaseCal: SNSmo SN "+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
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
        SNver +=1    
    # End SNSmo

    # Plot fits?
    if doSNPlot:
        retCode = VLBAPlotTab(uv, "SN", SNver, err, nplots=6, optype="PHAS", \
                                  logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
  
        retCode = VLBAWritePlots (uv, 1, 0, plotFile, err, \
                                      logfile=logfile, check=check, debug=debug)
        if retCode!=0:
            return retCode
        
    # Apply to CL table
    retCode = VLBAApplyCal(uv, err, maxInter=600.0, logfile=logfile, check=check,debug=debug)
    if retCode!=0:
        return retCode
    return 0
    # end VLBAPhaseCal

def VLBADoppler(uv, err, Sources=None, timeRange=[0.,0.,0.,0., 0.,0.,0.,0], FreqID=1, \
                    doBand=-1, BPVer=0, flagVer=-1,  \
                    noScrat=[], logfile='', check=False, debug=False):
    """ Diurnal Doppler corrections

    Determine and apply Doppler corrections producing a new output file, class = "CVel"
    Returns task error code, 0=OK, else failed
    err        = Python Obit Error/message stack
    Sources    = if given, the list of sources to shift
    timeRange  = (AIPS style) timerange of data to use
    FreqID     = Frequency group identifier
    doBand     = If >0.5 apply bandpass cal.
    BPVer      = Bandpass table version
    flagVer    = Input Flagging table version
    noScrat    = list of disks to avoid for scratch files
    logfile    = logfile for messages
    check      = Only check script, don't execute tasks
    debug      = show input
    """
    ################################################################
    mess = "Doppler correction"
    printMess(mess, logfile)

    # Set output (new) SN table
    SNver = uv.GetHighVer("AIPS SN")+1

    cvel = AIPSTask.AIPSTask("cvel")
    cvel.logFile   = logfile
    if not check:
        setname(uv,cvel)
    i = 1;
    cvel.outclass  = "CVel"
    cvel.outdisk   = cvel.indisk
    cvel.outseq    = cvel.inseq
    cvel.flagver   = flagVer
    cvel.freqid    = FreqID
    cvel.doband    = doBand
    cvel.bpver     = BPVer
    i = 1;
    for t in timeRange:
        cvel.timerang[i] = t
        i += 1
    i = 1;
    if Sources:
        for s in Sources:
            cvel.sources[i] = s
            i += 1
    i = 1;
    for d in noScrat:
        cvel.baddisk[i] = d
        i += 1
    if debug:
        cvel.prtLv = 2
        cvel.i
        cvel.debug = debug
    # Trap failure
    try:
        if not check:
            cvel.g
    except Exception, exception:
        print exception
        mess = "CVEL Failed retCode= "+str(cvel.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end VLBADoppler

def VLBAClipSN(uv, SNver, Flag, err, \
               logfile='', check=False, debug=False):
    """ Flag SN table entries with real < Flag

    Returns with err set on error
    uv         = UV data object
    SNver      = SN table to flag
    Flag       = if >0.0 flag SN table entries with real<Flag
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
    # Number of IFs
    nif   = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocif"]]
    # Number of Stokes
    npoln = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocs"]]
    SNTab = uv.NewTable(Table.READONLY, "AIPS SN", SNver, err, \
                        numIF=nif, numPol=npoln)
    if err.isErr:
        return
    # Open
    SNTab.Open(Table.READWRITE, err)
    if err.isErr:
        return
    # Number of rows
    nrow =  SNTab.Desc.Dict["nrow"]
    count = 0
    for i in range (0,nrow):    # Loop over rows
        SNrow = SNTab.ReadRow(i+1, err)
        if err.isErr:
            return
        dirty = False
        # Loop over IF
        for iif in range (0, nif):
            if SNrow["REAL1"][iif]<Flag:
                SNrow["WEIGHT 1"][iif] = 0.0
                count += 1
            # Second Poln
            if npoln>1:
                if SNrow["REAL2"][iif]<Flag:
                    SNrow["WEIGHT 2"][iif] = 0.0
                    count += 1
       # Rewrite if modified
        if dirty:
            SNTab.WriteRow(i+1, SNRow, err)
            if err.isErr:
                return
    # end loop over rows
            
    # Close table
    SNTab.Close(err)
    if err.isErr:
        return
    
    mess = "Flagged "+str(count)+" Gain entries from ACCOR"
    printMess(mess, logfile)
    # end VLBAClipSN
    
def VLBASNStats(uv, SNver, solInt, err, refAnts=[0], logfile='', check=False, debug=False):
    """ Find good timerange on the basis of an SN table

    Returns with err set on error
    Return dict {"Source":source, "souID":souID, "timeRange":(tbeg,tend), \
                  "Fract":fract_OK, "SNR":avg_SNR,"bestRef":refAnt}
    If there is no SU table or source ID not in table, source name is blank
    uv         = UV data object
    SNver      = SN table to flag
    solInt     = statistics interval (min)
    refAnts    = If values given, then list of acceptable ref. ants
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
    # Number of IFs
    nif   = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocif"]]
    # Number of Stokes
    npoln = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocs"]]
    npoln = min(2, npoln)   # No more than 2
    SNtab = uv.NewTable(Table.READONLY, "AIPS SN", SNver, err, \
                        numIF=nif, numPol=npoln)
    if err.isErr:
        return
    # Make sure sorted
    Table.PSort(SNtab, "TIME", False, err)
    if err.isErr:
        return
    # Number of antennas
    nant = SNtab.Desc.List.Dict["NO_ANT"][2][0]
    # Open
    SNtab.Open(Table.READONLY, err)
    if err.isErr:
        return
    # Number of rows
    nrow =  SNtab.Desc.Dict["nrow"]
    # Initialize
    solIntD = solInt/1440.
    time0   = -1.0e20       # Beginning time of current interval
    timeE   = -1.0e20       # End time of current interval
    tlast   = -1.0e20       # Last time
    souId   = -10           # Current source ID
    accum   = []            # solution period statistics array
    times   = None
    SNR1    = None
    SNR2    = None
    antCnt  = None
    fract   = 0.0
    avgSNR  = 0.0
    
    totAnt  = []  # Total valid IF/poln count per antenna
    snrAnt  = []  # Total valid IF/poln SNR per antenna
    for i in range(0,nant):
        totAnt.append(0)
        snrAnt.append(0.0)
    
    # For each interval collect an accum entry containing
    # 0) source ID
    # 1) (beginning_time, end_time)
    # 2) Fraction of total ant/IF/poln occuring
    # 3) Average SNR
    # 4) [antenna_occurance_count]
    # 5) [[avg_SNR_per_ant/IF]]  Poln 1
    # 6) [[avg_SNR_per_ant/IF]]  Poln 2
    for i in range (0,nrow):    # Loop over rows
        SNrow = SNtab.ReadRow(i+1, err)
        if err.isErr:
            return
        time     = SNrow["TIME"][0]
        curSouID = SNrow["SOURCE ID"][0]
        # New interval?
        if (time>timeE) or (souID!=curSouID):
            # Save any current values to accum
            if times:
                times[1] = tlast    # actual end time
                # Normalize accumulations by counts, overall statistics
                acnt   = 0
                sum    = 0.0
                fract  = 0.0
                avgSNR = 0.0
                for i in range(0,nant):
                    for j in range(0,nif):
                        if CNT1[i][j]>0:
                            totAnt[i]  += CNT1[i][j]
                            snrAnt[i]  += SNR1[i][j]
                            SNR1[i][j] /= CNT1[i][j]
                            acnt += 1
                            sum  += SNR1[i][j]
                        if (npoln>1) and  CNT2[i][j]>0:
                            snrAnt[i]  += SNR2[i][j]
                            totAnt[i]  += CNT2[i][j]
                            SNR2[i][j] /= CNT2[i][j]
                            acnt += 1
                            sum  += SNR1[i][j]
                if acnt>0:
                    avgSNR = sum / acnt
                fract = float(acnt) / float(nant*nif*npoln)
                pastSI = [souID, times, fract, avgSNR, antCnt, SNR1, SNR2]
                accum.append(pastSI)
            # Build new accumulators
            times = [time,time+solIntD]
            antCnt = []            # Occurences of this antenna
            SNR1   = []            # Antenna array of sums of poln1
            CNT1   = []            # Antenna array of counts of poln1
            for i in range(0,nant):
                antCnt.append(0)
                # per antenna
                snr1 = []
                cnt1 = []
                for j in range(0,nif):
                    snr1.append(0.0)
                    cnt1.append(0)
                SNR1.append(snr1)
                CNT1.append(cnt1)
            # Second poln?
            if npoln>1:
                SNR2 = []          # Antenna array of sums of poln2
                CNT2 = []          # Antenna array of sums of poln2
                for i in range(0,nant):
                    snr2 = []
                    cnt2 = []
                    for j in range(0,nif):
                        snr2.append(0.0)
                        cnt2.append(0)
                    SNR2.append(snr2)
                    CNT2.append(cnt2)
            # end build accumulators
        timeE = time + solIntD
        souID = curSouID
        tlast = time
        # end new period

        # Accumulate
        tlast = time  # Save last time
        iant = SNrow["ANTENNA NO."][0] - 1   # 0-rel antenna no.
        antCnt[iant] += 1
        # Loop over IF
        for iif in range (0, nif):
            if SNrow["WEIGHT 1"][iif]>0.0:
                SNR1[iant][iif] += SNrow["WEIGHT 1"][iif];
                CNT1[iant][iif] += 1;
            # Second Poln
            if npoln>1:
                if SNrow["WEIGHT 1"][iif]>0.0:
                    SNR2[iant][iif] += SNrow["WEIGHT 2"][iif];
                    CNT2[iant][iif] += 1;

    # end loop over rows
            
    # Close table
    SNtab.Close(err)
    if err.isErr:
        return

    # Find highest fraction
    hiFract = 0.0
    for s in accum:
        hiFract = max (hiFract, s[2])

    # eliminate (negate avg SNR) entries with lower fract
    for s in accum:
        if s[2]<0.99*hiFract:
            s[3] = -s[3]

    # Find highest avg SNR
    hiSNR = 0.0
    hi = None
    for s in accum:
        if s[3]>hiSNR:
            hiSNR = s[3]
            hi = s

    # Normalize antenna average SNRs
    for i in range (0,nant):
        if totAnt[i]>0:
            snrAnt[i] /= totAnt[i]

    # deselect antennas not in refAnts (if any non zero)
    if refAnts[0]>0:
        for i in range (0,nant):
            drop = True
            for ra in refAnts:
                if ra==i+1:
                    drop = False
            # found it?
            if drop:
                snrAnt[i] = 0.0
                totAnt[i] = 0
    
    # Find best refant count - one with most valid occurences
    bestCnt = 0
    for i in range (0,nant):
        if totAnt[i]>bestCnt:
            bestCnt = totAnt[i]

    # Find antenna with count equal to bestCnt with highest SNR
    bestRef = 0
    hiSNR = 0.0
    for i in range (0,nant):
        if (totAnt[i]>=bestCnt) and (snrAnt[i]>hiSNR):
            bestRef = i+1
            hiSNR   = snrAnt[i]

    # Lookup source name if SU table present
    hiSU = uv.GetHighVer("AIPS SU")
    souName = "            "    # default
    if hiSU>= 1:
        SUtab = uv.NewTable(Table.READONLY, "AIPS SU", 1, err, \
                        numIF=nif,)
        SUtab.Open(Table.READONLY, err)
        if err.isErr:
            return
        # Number of rows
        nrow =  SUtab.Desc.Dict["nrow"]
        for i in range (0,nrow):    # Loop over rows
            SUrow = SUtab.ReadRow(i+1, err)
            if err.isErr:
                return
            curSouID = SUrow["ID. NO."][0]
            if curSouID==hi[0]:   # This it?
                souName = SUrow["SOURCE"][0]
                break;
        SUtab.Close(err)
        if err.isErr:
            return
    
    if debug:
        print totAnt,"\n", snrAnt,"\n"
        for s in accum:
            print s[0],s[1],s[2],s[3]

    # Create output structure
    out = {"Source":souName, "souID":hi[0],"timeRange":hi[1], "Fract":hi[2], "SNR":hi[3], "bestRef":bestRef}
    return out
    # end VLBSNStats

def VLBAImFITS(inImage, filename, outDisk, err, fract=None, quant=None, \
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
    # Open and close input to update
    inImage.Open(Image.READONLY,err)
    inImage.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error Opening image")
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
    # end VLBAImFITS

def VLBAUVFITS(inUV, filename, outDisk, err, compress=False, \
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
    inInfo = UV.PGetList(outUV)    # 
    dim = [1,1,1,1,1]
    #InfoList.PAlwaysPutInt (inInfo, "corrType", dim, [1])  
    inInfo.set("corrType", [1]) 
    #Compressed?
    if compress:
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
    # end VLBAUVFITS

def VLBAUVFITSTab(inUV, filename, outDisk, err, \
              exclude=["AIPS HI", "AIPS AN", "AIPS FQ", "AIPS SL", "AIPS PL"], \
                  include=[]):
    """ Write Tables on UV data as FITS file
    
    Write Tables from a UV data set (but no data) as a FITAB format file
    History written to header
    inUV       = UV data to copy
    filename   = name of FITS file, any whitespace characters replaced with underscore 
    inDisk     = FITS directory number
    err        = Python Obit Error/message stack
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy (FQ, AN always done )
                 Exclude has presidence over include
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
    # end VLBAUVFITSTab

def VLBAMedianFlag(uv, target, err, \
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
    except Exception, exception:
        print exception
        mess = "Median flagging Failed retCode="+str(medn.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end VLBAMedianFlag
    
def VLBAQuack(uv, err, \
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
    quack.logFile   = logfile
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
    # end VLBAQuack
    
def VLBAAutoFlag(uv, target, err, \
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
    except Exception, exception:
        print exception
        mess = "AutoFlag Failed retCode="+str(af.retCode)
        printMess(mess, logfile)
        return 1
    else:
        pass
    return 0
    # end VLBAAutoFlag

def VLBACal(uv, target, ACal, err, \
                PCal=None, FQid=0, calFlux=None, \
                doCalib=-1, gainUse=0, doBand=0, BPVer=0, flagVer=-1, \
                calModel=None, calDisk=0, \
                solnVer=1, solInt=10.0/60.0, solSmo=0.0, nThreads=1, refAnt=0, ampScalar=False, \
                check=False, debug = False, \
                noScrat=[], logfile = ""):
    """ Basic Amplitude and phase cal for VLBA data

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
    calib.logFile  = logfile
    if not check:
        setname(uv,calib)
    calib.Sources  = [ACal]
    calib.flagVer  = flagVer
    calib.ampScalar= ampScalar
    calib.doCalib  = doCalib
    calib.gainUse  = gainUse
    calib.doBand   = doBand
    calib.BPVer    = BPVer
    calib.solMode  ="A&P"
    calib.solnVer  = solnVer
    calib.nThreads = nThreads
    calib.solInt   = solInt
    calib.refAnts  = [refAnt]
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
        clcal.logFile  = logfile
        if not check:
            setname(uv,clcal)
        clcal.calSour = [ACal]
        if not check:
            clcal.calIn   = uv.GetHighVer("AIPS CL")
        else:
            clcal.calIn   = 1
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
                calib.flagVer   = flagVer
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
            except Exception, exception:
                print exception
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
            if not check:
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
            mess = "VLBACal: SNSmo SN"+str(snsmo.solnIn)+" to "+str(snsmo.solnOut)
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
        else:
            solnVer2 = solnVer
        # GetJy to set flux density scale if ACal not in PCal
        if not ACal in PCal:
            getjy = ObitTask.ObitTask("GetJy")
            getjy.logFile  = logfile
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
        
    # Set up for CLCal - only use phase calibrators
    if not check:
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
    #clcal.debug=True
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
    # end VLBACal

def VLBASplit(uv, target, err, FQid=1, outClass="      ", logfile = "", \
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
    # end VLBAsplit

def VLBACalAvg(uv, avgClass, avgSeq, CalAvgTime,  err, \
               FQid=0, \
               flagVer=0, doCalib=2, gainUse=0, doBand=1, BPVer=0,  \
               BIF=1, EIF=0, BChan=1, EChan=0, \
               avgFreq=0, chAvg=1, Compress=False, \
               logfile = "", check=False, debug=False):
    """ Calibrate, select and/or average data to a multisource file

    Returns task error code, 0=OK, else failed
    Generates NX and initial dummy CL table if needed
    uv         = UV data object to clear
    avgClass   = Class name of averaged data
    avgSeq     = Sequence number of averaged data. 0 => highest unique
    CalAvgTime = Averaging time in sec
    err        = Obit error/message stack
    FQid       = Frequency Id to process, 0=>all
    doCalib    = Apply calibration table, positive=>calibrate
    gainUse    = CL/SN table to apply
    doBand     = If >0.5 apply previous bandpass cal.
    BPVer      = previous Bandpass table (BP) version
    BIF        = first IF to copy
    EIF        = highest IF to copy
    BChan      = first channel to copy
    EChan      = highest channel to copy
    flagVer    = Input Flagging table version
    avgFreq    = If 0 < avgFreq <= 1 then avg chans (see splat for details)
    chAvg      = Number of channels to average
    Compress   = Write "Compressed" data?
    logfile    = Log file for task
    check      = Only check script, don't execute tasks
    debug      = Run tasks debug, show input
    """
    ################################################################
    splat=ObitTask.ObitTask("Splat")
    splat.logFile = logfile
    if not check:
        setname(uv,splat)
    splat.doCalib  = doCalib
    splat.gainUse  = gainUse
    splat.doBand   = doBand
    splat.BPVer    = BPVer
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
    # end VLBACalAvg
    
def VLBACalAvg2(uv, avgClass, avgSeq, CalAvgTime,  err, \
                Source=None, FQid=0, doPol=False, \
                flagVer=0, doCalib=2, gainUse=0, doBand=1, BPVer=0,  \
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
    Source     = Selected source or list of sources
    FQid       = Frequency Id to process, 0=>all
    doPol      = Calibratio polarization?
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
    check      = Only check script, dont execute tasks
    debug      = Run tasks debug, show input
    """
    ################################################################
    # Create output
    if not check:
        # Set calibration, editing and selection
        info = uv.List
        info.set("doCalSelect",  True)
        if type(Source)==list:
            info.set("Sources",    Source)
        elif calSou:
             info.set("Sources",    [Source])
        info.set("FreqID",     FQid)
        info.set("doPol",      doPol)
        info.set("doCalib",    doCalib)
        info.set("gainUse",    gainUse)
        info.set("copyCalTab", True)   # Copy CL/SN tables if not calibrating
        info.set("doBand",     doBand)
        info.set("BPVer",      BPVer)
        info.set("BIF",        BIF)
        info.set("EIF",        EIF)
        info.set("BChan",      BChan)
        info.set("EChan",      EChan)
        info.set("Stokes",     "    ")
        info.set("flagVer",    flagVer)
        #print "info", info.Dict  # DEBUG
        # Open and close to set
        uv.Open(UV.READCAL, err)
        outuv = UV.newPAUV("CalAvg", uv.Aname, avgClass, uv.Disk, avgSeq, False, err)
        info = outuv.List
        info.set("Compress", Compress,)
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
        except Exception, exception:
            print exception
            OErr.printErr(err)
            mess = "Calibrate and average uv data failed"
            printMess(mess, logfile)
            return 1
        else:
            pass
              

    if not check:
        try:
            # Do History - previous already copied
            ooutuv = UV.newPAUV("CalAvg", uv.Aname, avgClass, uv.Disk, avgSeq, True, err)
            inHistory  = History.History("inhistory",  uv.List, err)
            outHistory = History.History("outhistory", ooutuv.List, err)

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
                OErr.printErrMsg(err,"History Error")
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
            UV.PTableCLGetDummy(ooutuv, ooutuv, 0, err, solInt=solint)
            if err.isErr:
                print "Error creating cal/avg AIPS data CL table"
                OErr.printErrMsg(err,"Error in dummy CL table")
            # Index
            UV.PUtilIndex (ooutuv, err)
            if err.isErr:
                print  "Error indexing cal/avg AIPS data"
                OErr.printErrMsg(err,"Error indexing")
        except Exception, exception:
            print exception
            OErr.printErr(err)
            mess = "Indexing or creating CL table failed"
            printMess(mess, logfile)
            return 1
        else:
            pass
    return 0
    # end VLBACalAvg2
    
def VLBASetImager (uv, target, outIclass="", check=False, nThreads=1, \
                   noScrat=[], logfile = ""):
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
# end VLBASetImager

def VLBAPolCal(uv, InsCals, err, \
               doCalib=2, gainUse=0, doBand=1, BPVer=0, flagVer=-1, \
               soltype="ORI-", fixPoln=False, avgIF=False, \
               solInt=0.0, refAnt=0, doSetJy=False, \
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
    doCalib  = Apply prior calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    doBand   = >0 => apply bandpass calibration
    BPVer    = AIPS BP table to apply
    flagVer  = Input Flagging table version
    soltype  = solution type
    fixPoln  = if True, don't solve for source polarization in ins. cal
    avgIF    = if True, average IFs in ins. cal.
    solInt   = instrumental solution interval (min), 0=> scan average
    refAnt   = Reference antenna
    doSetJy  = If True, use SetJy to set flux densities of calibrators to 1
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
    # Need to set flux densities?
    if doSetJy:
        setjy=ObitTask.ObitTask("SetJy")
        if not check:
            setname(uv,setjy)
        if type(InsCals)==list:
            setjy.Sources=InsCals
        else:
            setjy.Sources[0]=InsCals
        setjy.ZeroFlux[0] = 1.0
        setjy.logFile     = logfile
        setjy.debug       = debug
        if debug:
            setjy.i
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
    # end doSetJy

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
    # End VLBAPolCal

def VLBARLCal(uv, err, RLPCal=None, \
              BChan=1, EChan = 0, ChWid=1, solInt1=1./6, solInt2=10., \
              UVRange=[0.,0.], WtUV=0.0, \
              FQid=0, calcode="    ", doCalib=-1, gainUse=0, \
              timerange = [0.,0.], Antennas=None, \
              doBand=0, BPVer=0, flagVer=-1, BPSoln=0, \
              refAnt=0, doPol=-1, smooth=[0.,0.,0.], \
              nThreads=1, noScrat=[], logfile = "",check=False, debug = False):
    """ Determine R-L phase bandpass

    Determines new BP table by looping over entries in RLPCal
    Returns task error code, 0=OK, else failed
    uv       = UV data object to clear
    err      = Obit error/message stack
    RLPCal   = An array of triplets with R-L calibrators:
               (name, R-L phase (deg at 1 GHz), RM (rad/m**2))
    BChan    = First (1-rel) channel number
    EChan    = Highest channel number. 0=> high in data.
    ChWid    = Number of channels to average to determine soln.
    solInt1  = first solution interval (min), 0=> scan average
    solInt2  = second solution interval (min)
    UVRange  = Range of baseline used in kilowavelengths
    WTUV     = Weight to use outside of UVRange
    FQid     = Frequency Id to process
    calcode  = Calibrator code
    doCalib  = Apply calibration table, positive=>calibrate
    gainUse  = CL/SN table to apply
    timerange= time range of data (days)
    Antennas = if given, a list of antennas to use
    doBand   = If >0.5 apply previous bandpass cal.
    BPVer    = previous Bandpass table (BP) version
    BPSoln   = Output BP table, 0=>new
    flagVer  = Flagging table to apply
    refAnt   = Reference antenna REQUIRED
    doPol    = Apply polarization cal?
    smooth   = Channel smoothing function
    noScrat  = list of AIPS disks to avoid for scratch files
    nThreads = Number of threads to use in imaging
    logfile  = Log file for task
    check    = Only check script, don't execute tasks
    debug    = Run tasks debug, show input
    """
    ################################################################
    # Anything to do?
    if (not RLPCal) or (len(RLPCal)<=0):
        return 0
    ncal = len(RLPCal)  # How many calibrators?
    rlpass=ObitTask.ObitTask("RLPass")
    rlpass.taskFile = logfile
    if not check:
        setname(uv,rlpass)
    if Antennas:
        i = 0
        for a in Antennas:
            rlpass.Antennas[i] = a; i  += 1
    rlpass.timeRange[0] = timerange[0];rlpass.timeRange[1] = timerange[1];
    rlpass.BChan1  = BChan
    rlpass.EChan1  = EChan
    rlpass.BChan2  = BChan
    rlpass.EChan2  = EChan
    rlpass.ChWid2  = ChWid
    rlpass.UVRange[0] = UVRange[0]; rlpass.UVRange[1] = UVRange[1];
    rlpass.doCalib = doCalib
    rlpass.gainUse = gainUse
    rlpass.flagVer = flagVer
    rlpass.FreqID  = FQid
    rlpass.souCode = calcode
    rlpass.doPol   = doPol
    rlpass.doBand  = doBand
    rlpass.BPVer   = BPVer
    rlpass.refAnt  = refAnt
    rlpass.solInt1 = solInt1
    rlpass.solInt2 = solInt2
    rlpass.BPSoln  = BPSoln
    rlpass.prtLv   = 1
    rlpass.nThreads = nThreads
    # Loop over calibrators
    for ical in range (0,ncal):
        rlpass.Sources[0]= RLPCal[ical][0]
        rlpass.RLPhase   = RLPCal[ical][1]
        rlpass.RM        = RLPCal[ical][2]
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
    return 0
    # end VLBARLCal

def VLBARLCal2(uv, err, uv2 = None, \
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
                    SumCC = VLBAGetSumCC (x,err)
                    h = x.Desc.Dict
                    blc = [h["inaxes"][0]/4,h["inaxes"][1]/4]
                    trc = [3*h["inaxes"][0]/4,3*h["inaxes"][1]/4]
                    stat = imstat(x, err, blc=blc,trc=trc)
                    IFlux.append(SumCC)
                    IRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes Q
                    outClass="QPOLCL"
                    x =  Image.newPAImage("Q",outName, outClass, outDisk,outSeq,True,err)
                    SumCC = VLBAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc)
                    QFlux.append(SumCC)
                    QRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes U
                    outClass="UPOLCL"
                    x =  Image.newPAImage("U",outName, outClass, outDisk,outSeq,True,err)
                    SumCC = VLBAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc)
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
                    SumCC = VLBAGetSumCC (x,err)
                    h = x.Desc.Dict
                    blc = [h["inaxes"][0]/4,h["inaxes"][1]/4]
                    trc = [3*h["inaxes"][0]/4,3*h["inaxes"][1]/4]
                    stat = imstat(x, err, blc=blc,trc=trc)
                    IFlux.append(SumCC)
                    IRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes Q
                    outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                    x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                    SumCC = VLBAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc)
                    QFlux.append(SumCC)
                    QRMS.append(stat["RMSHist"])
                    x.Zap(err)  # Cleanup
                    del x
                    # Stokes U
                    outFile  = img.Sources[0].strip()+"ITEMPPOLCAL.fits"
                    x =  Image.newPFImage("Q",outFile,img.outDisk,True,err)
                    SumCC = VLBAGetSumCC (x,err)
                    stat = imstat(x, err, blc=blc,trc=trc)
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
            VLBACopyTable (uv, uv2, "AIPS AN", err, \
                           logfile=logfile, check=check, debug=debug)
            if err.isErr:
                print  "Error copying AN Table"
                return 1
        # Copy CL table to be modified
        VLBACopyTable (uv, uv, "AIPS CL", err, inVer=gainUse, outVer=hiCL+1, \
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
            #clcor.debug = debug
        # Trap failure
        try:
            if not check:
                #clcor.i
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
                hiCL = uv2.GetHighVer("AIPS CL")
                # Copy CL table to be modified
                VLBACopyTable (uv, uv, "AIPS CL", err, inVer=hiCL, outVer=hiCL+1, \
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
    # end VLBARLCal2

def VLBAReportTargets(uv, err,  FreqID=1, Sources=None, seq=1, sclass="IClean", \
                          Stokes="I", logfile='', check=False, debug=False):
    """ Generate report info for a list of targets in AIPS files

    Returns a report which is a list of dicts, each of which contains
        "Source":   Source name
        "ObsDate":  Observing date as "yyyy-mm-dd"
        "RA"    :   Source RA (deg0 at standard equinox
        "Dec"   :   Source Dec (deg) at standard equinox
        "RAPnt" :   Antenna pointing RA (deg) at standard equinox
        "DecPnt":   Antenna pointing Dec (deg) at standard equinox
        "Freq" :    Reference frequency (Hz)
        "numVis":   Number of visibilities (ignoring flagging)
        "Size"  :   Width of image in deg (From Stokes I)
        "Cells" :   Cell spacing in deg (From Stokes I)
        "Exposure": Total integration time (day)
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
        slist = VLBAAllSource(uv,err,logfile=logfile,check=check,debug=debug)
    else:
        slist = sl
    
    # Init output
    Report = []

    # Image disk assumed same as uv
    disk = uv.Disk
    
    # Loop over slist
    for sou in slist:
        sdict = {"Source":sou}  # Init source structure
        # Image statistics, loop over Stokes
        for s in Stokes:
            klass = s+sclass[1:]
            x = Image.newPAImage(s, sou, klass, disk, seq, True, err)
            hd = x.Desc.Dict
            sdict[s+"Beam"] = (hd["beamMaj"],hd["beamMin"],hd["beamPA"])
            # Some from Stokes I only
            if s == 'I':
                sdict["Size"]    = hd["inaxes"][1]*hd["cdelt"][1]
                sdict["Cells"]   = hd["cdelt"][1]
                sdict["RA"]      = hd["crval"][0]
                sdict["Dec"]     = hd["crval"][1]
                sdict["RAPnt"]   = hd["obsra"]
                sdict["DecPnt"]  = hd["obsdec"]
                sdict["ObsDate"] = hd["obsdat"]
                sdict["Freq"]    = hd["crval"][hd["jlocf"]]
            blc = [hd["inaxes"][0]/4,hd["inaxes"][1]/4]
            trc = [3*hd["inaxes"][0]/4,3*hd["inaxes"][1]/4]
            stat = imstat(x,err,blc=blc,trc=trc)  # Image statistics inner quarter
            sdict[s+"Peak"] = stat["Max"]
            sdict[s+"RMS"]  = stat["RMSHist"]
            sdict[s+"Sum"]  = VLBAGetSumCC(x, err, logfile=logfile, check=check, debug=debug)
        # End stokes image loop
        # Observing stats
        obstat = VLBAGetTimes (uv, sou, err, logfile=logfile, check=check, debug=debug)
        sdict["numVis"]   = obstat["numVis"]
        sdict["Exposure"] = obstat["Exposure"]
        Report.append(sdict)  # Save source info
    # end loop over sources

    # Give terse listing
    for sdict in Report:
        mess = "\n Source = "+sdict["Source"]+", Exposure="+"%5.3f"%(sdict["Exposure"]*24.)+" hr"
        printMess(mess, logfile)
        mess = "IPol Beam = ("+"%8.3f"%(sdict["IBeam"][0]*3600000.0)+", %8.3f"%(sdict["IBeam"][1]*3600000.0)+ \
            ", %6.1f"%(sdict["IBeam"][2])+") mas, mas, deg"
        printMess(mess, logfile)
        for s in Stokes:
            mess = "Stokes "+s+" Sum CC="+"%8.3f"%(sdict[s+"Sum"])+", Peak="+"%8.3f"%(sdict[s+"Peak"])+ \
                ", RMS="+"%8.5f"%(sdict[s+"RMS"])+" Jy"
            printMess(mess, logfile)
    # End terse listing
    return Report
    # end VLBAReportTargets

def VLBAGetSumCC(image, err, CCver=1,
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
    # end VLBAGetSumCC

def VLBAGetTimes(uv, Source, err, 
                 logfile='', check=False, debug=False):
    """ Lookup observing times and number of visibilities for a source
    
    Return dict {"numVis":no vis, "Exposure":Total integration time (day)}
    uv         = UV data with AIPS SU and AIPS NX tables
    Source     = Source to lookup
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
    if check:
        return {"numVis":0, "Exposure":0.0}
    # Open and close uv to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)
    
    # Lookup Source ID (SouID)
    SUtab = uv.NewTable(Table.READONLY, "AIPS SU", 1, err)
    SUtab.Open(Table.READONLY, err)
    if err.isErr:
        return  {"numVis":0, "Exposure":0.0}
    # Number of rows
    nrow =  SUtab.Desc.Dict["nrow"]
    for i in range (0,nrow):    # Loop over rows
        SUrow = SUtab.ReadRow(i+1, err)
        if err.isErr:
            return  {"numVis":0, "Exposure":0.0}
        SouID = SUrow["ID. NO."][0]
        #if debug:
        #    mess="Source "+Source+" test "+SUrow["SOURCE"][0]+" ID ="+str(SouID)+ \
        #        " match="+str(SUrow["SOURCE"][0].rstrip()==Source.rstrip())
        #    printMess(mess, logfile)
        if SUrow["SOURCE"][0].rstrip()==Source.rstrip():   # This it?
            break;
    SUtab.Close(err)
    if err.isErr:
        return  {"numVis":0, "Exposure":0.0}
    
    # get observing stats from AIPS NX table
    cntVis  = 0
    sumTime = 0.0
    NXTab = uv.NewTable(Table.READONLY, "AIPS NX", 1, err)
    if err.isErr:
        return 0.0
    # Open
    NXTab.Open(Table.READONLY, err)
    if err.isErr:
        return {"numVis":cntVis, "Exposure":sumTime}
    # Number of rows
    nrow    = NXTab.Desc.Dict["nrow"]
    # Loop over table
    for irow in range (1, nrow+1):
        NXrow = NXTab.ReadRow(irow, err)
        if err.isErr:
            return {"numVis":cntVis, "Exposure":sumTime}
        #  Is this the desired source?
        if NXrow["SOURCE ID"][0]==SouID:
            sumTime += NXrow["TIME INTERVAL"][0]
            cntVis  += NXrow["END VIS"][0] - NXrow["START VIS"][0] + 1
    # End loop over table
    # Close table
    NXTab.Close(err)
    if err.isErr:
        return {"numVis":cntVis, "Exposure":sumTime}

    if debug:
        mess="VLBAGetTimes: Source "+Source+"="+str(SouID)+" numVis="+str(cntVis)+ \
            " Integration time = "+"%5.3f"%(sumTime*24.)+" hr"
        printMess(mess, logfile)
 
    return {"numVis":cntVis, "Exposure":sumTime}
    # end VLBAGetTimes

def VLBARefMB(uv, SNver, err, smoTime=1.0e20, \
              logfile='', check=False, debug=False):
    """ Average delays and reference phases to first IF
    
    Averages all delays per antenna/Poln/IF in an SN table and 
    referes all phases to IF 1 (or first IF with valid data).
    This preserves multiband delays but removes phase fluctuations in time.
    Intended for use for a single scan SN table
    Returns with err set on error
    uv         = UV data object
    SNver      = SN table to modify
    err        = Python Obit Error/message stack
    smoTime    = Averaging interval (min)
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
    mess = "VLBARefMB: Filter SN "+str(SNver)+" to remove phase fluctuations, average delays"
    printMess(mess, logfile)
    # Open and close uv to sync with disk 
    uv.Open(UV.READONLY, err)
    uv.Close(err)

    fblank = FArray.PGetBlank()  
    # Number of IFs
    nif   = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocif"]]
    # Number of Stokes
    npoln = uv.Desc.Dict["inaxes"][uv.Desc.Dict["jlocs"]]
    SNTab = uv.NewTable(Table.READONLY, "AIPS SN", SNver, err, \
                        numIF=nif, numPol=npoln)
    if err.isErr:
        return

    nant = SNTab.Desc.List.Dict["NO_ANT"][2][0]  # Number of antennas
    # Create storage arrays
    bigList  = []
    foundAnt = []
    for iant in range(0,nant):
        foundAnt.append(False)
    
    # Open
    SNTab.Open(Table.READWRITE, err)
    if err.isErr:
        return
    # Number of rows
    nrow    =  SNTab.Desc.Dict["nrow"]
    count   = 0
    irow    = 1
    hiRow   = 1
    loRow   = 1
    done    = False
    tstart  = -1.0e20
    smoDay = smoTime/1440.0   # Interval in days
    # Loop over time segments
    while not done:
        while irow<=nrow:  # Loop over rows
            #for i in range (lastRow,nrow):    # Loop over rows
            SNrow = SNTab.ReadRow(irow, err)
            if err.isErr:
                return
            if debug:
                mess = "VLBARefMB: read irow "+str(irow)+" time="+str( SNrow["TIME"][0])
                printMess(mess, logfile)
            # Time
            if tstart<-1.0e10:
                tstart = SNrow["TIME"][0]
            # End of interval
            if SNrow["TIME"][0]>tstart+smoDay:
                if debug:
                    mess = "VLBARefMB: Finished segment, irow "+str(irow)+" tstart="+str(tstart)
                    printMess(mess, logfile)
                break
            # Track antennas found
            ant = SNrow["ANTENNA NO."][0]
            foundAnt[ant-1] = True
            # Start row structure 
            datum = {"row":irow, "ant":ant, "suba":SNrow["SUBARRAY"][0], \
                     "MBDelay1":SNrow["MBDELAY1"][0]}
            if SNrow["MBDELAY1"][0]!=fblank:
                datum["MBDelay1"] = SNrow["MBDELAY1"][0]
            else:
                datum["MBDelay1"] = 0.0
            if (npoln>1) and SNrow["MBDELAY2"][0]!=fblank:
                datum["MBDelay2"] = SNrow["MBDELAY2"][0]
            else:
                datum["MBDelay2"] = 0.0
            # Phase of first IF with data
            phaseR0 =  -1.0e20
            phaseL0 =  -1.0e20
            # First poln
            if1 = []
            for iif in range (0, nif):
                if SNrow["WEIGHT 1"][iif]>0.0:
                    ph = math.atan2(SNrow["IMAG1"][iif],SNrow["REAL1"][iif])
                    if phaseR0<-1.0e19:
                        phaseR0 = ph
                    phase =  ph - phaseR0
                    if1.append([SNrow["WEIGHT 1"][iif], phase, SNrow["DELAY 1"][iif]])
                else:
                    if1.append([0.0,0.0,0.0])
            datum["IF1"] = if1   # save
            # Second Poln
            if npoln>1:
                if2 = []
                for iif in range (0, nif):
                    if SNrow["WEIGHT 2"][iif]>0.0:
                        ph = math.atan2(SNrow["IMAG2"][iif],SNrow["REAL2"][iif])
                        if phaseL0<-1.0e19:
                            phaseL0 = ph
                        phase =  ph - phaseL0
                        if2.append([SNrow["WEIGHT 2"][iif], phase, SNrow["DELAY 2"][iif]])
                    else:
                        if2.append([0.0,0.0,0.0])
                datum["IF2"] = if2   # save
            # Add datum to bigList
            bigList.append(datum)
            hiRow = irow
            irow += 1;  # Next row
        # end loop over rows reading
    
        # Loop over antennas averaging
        for ant in range (1,nant+1):
            if not foundAnt[ant-1]:
                continue
            # Loop over IFs
            for iif in range  (0, nif):
                # loop over data entries
                sumMB1 = 0.0; cntMB1 = 0
                sumDly = 0.0; sumRe  = 0.0; sumIm  = 0.0; sumWt  = 0.0
                for s in bigList:
                    # Want antenna?
                    if s["ant"] != ant:
                        continue
                    sumMB1 += s["MBDelay1"]
                    cntMB1 += 1
                    # Valid
                    d = s["IF1"][iif]
                    if d[0]>0:
                        sumWt  += d[0]
                        sumRe  += d[0]*math.cos(d[1])
                        sumIm  += d[0]*math.sin(d[1])
                        sumDly += d[0]*d[2]
                # end loop over lists
                # average
                if cntMB1>0:
                    sumMB1 /= cntMB1
                if sumWt>0.0:
                    sumDly /= sumWt
                phase = math.atan2(sumIm, sumRe)
                # Save averages
                for s in bigList:
                    # Want antenna?
                    if s["ant"] != ant:
                        continue
                    s["MBDELAY1"] = sumMB1
                    d = s["IF1"][iif]
                    d[1] = phase
                    d[2] = sumDly
                # end loop over lists
                # Second Poln
                if npoln>1:
                    # loop over data entries
                    sumMB1 = 0.0; cntMB1 = 0
                    sumDly = 0.0; sumRe  = 0.0; sumIm  = 0.0; sumWt  = 0.0
                    for s in bigList:
                        # Want antenna?
                        if s["ant"] != ant:
                            continue
                        #print "DEBUG2", s, iif
                        sumMB1 += s["MBDelay2"]
                        cntMB1 += 1
                        # Valid
                        d = s["IF2"][iif]
                        if d[0]>0:
                            sumWt  += d[0]
                            sumRe  += d[0]*math.cos(d[1])
                            sumIm  += d[0]*math.sin(d[1])
                            sumDly += d[0]*d[2]
                    # end loop over lists
                    # average
                    if cntMB1>0:
                        sumMB1 /= cntMB1
                    if sumWt>0.0:
                        sumDly /= sumWt
                    phase = math.atan2(sumIm, sumRe)
                    # Save averages
                    for s in bigList:
                        # Want antenna?
                        if s["ant"] != ant:
                            continue
                        s["MBDELAY2"] = sumMB1
                        d = s["IF2"][iif]
                        d[1] = phase
                        d[2] = sumDly
                    # end loop over lists
                # end IF loop
            # end antenna loop
        
        # Loop over entries rewriting
        for s in bigList:
            i = s["row"]
            SNrow = SNTab.ReadRow(i, err)
            if err.isErr:
                return
            # Loop over IF
            for iif in range (0, nif):
                d = s["IF1"][iif]
                if (d[0]>0.0) and (SNrow["WEIGHT 1"][iif]>0.0):
                    SNrow["REAL1"][iif]  = math.cos(d[1])
                    SNrow["IMAG1"][iif]  = math.sin(d[1])
                    SNrow["DELAY 1"][iif] = d[2]
                # Second Poln
                if npoln>1:
                    d = s["IF2"][iif]
                    if (d[0]>0.0) and (SNrow["WEIGHT 2"][iif]>0.0):
                        SNrow["REAL2"][iif]  = math.cos(d[1])
                        SNrow["IMAG2"][iif]  = math.sin(d[1])
                        SNrow["DELAY 2"][iif] = d[2]
            # Rewrite
            SNrow = SNTab.WriteRow(i, SNrow, err)
            if err.isErr:
                return
        # end loop over rows rewriting
        done = irow>=nrow   # Finished?
        if debug:
            mess = "VLBARefMB: Finished pass, irow "+str(irow)+" tstart="+str(tstart)
            printMess(mess, logfile)
        # Clear accumulators for next period
        tstart  = -1.0e20
        bigList  = []
        foundAnt = []
        for iant in range(0,nant):
            foundAnt.append(False)
    # end loop over time segments
           
    # Close table
    SNTab.Close(err)
    if err.isErr:
        return
    
    mess = "VLBARefMB: Processed SN table "+str(SNver)
    printMess(mess, logfile)
    # end VLBARefMB
    
def VLBAImageModel(snames, sclass, sdisk, sseq, err, \
              logfile='', check=False, debug=False):
    """ Generate Dict of image Models from images

    Create a source list dict with source models for AIPS images.
    Creates an image object for each image found matching criteria
    Return python dict with image model dicts keyed on source name,
    image object = "Image"
    names      = source name or list of names, None=noop
    sclass     = image class
    sdisk      = image disk
    sseq       = image sequence
    err        = Python Obit Error/message stack
    logfile    = logfile for messages
    check      = Only check script
    debug      = Only debug - no effect
    """
    ################################################################
    out = None
    if check or (not snames):
        return out
    if type(snames)==list:  # list
        for name in snames:
            try:
                cno = AIPSDir.PFindCNO(sdisk, OSystem.PGetAIPSuser(), name, sclass, "MA", sseq, err)
                if cno>0:
                    x =  Image.newPAImage(name, name, sclass, sdisk, sseq, True, err)
                else:
                    x = None
            except Exception, exception:
                err.Clear()
            else:
                if x:
                    if not out:  # Need to create empty dist?
                        out = {}
                    out[name] = {"image":x}  # add entry
    else:   # Single
        try:
            cno = AIPSDir.PFindCNO(sdisk, OSystem.PGetAIPSuser(), name, sclass, "MA", sseq, err)
            if cno>0:
                x =  Image.newPAImage(snames, snames, sclass, sdisk, sseq, True, err)
            else:
                x = None
        except Exception, exception:
            err.Clear()
        else:
            if x:
                out = {sname:{"image":x}}
    return out
    # end VLBAImageModel
    
def VLBAAllSource(uv, err, logfile='', check=False, debug=False):
    """ 
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
        printMess("Debug mode: returning empty source list")
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
    # end VLBAAllSource
    
def VLBASNAppend (inSN, outSN, err):
    """  Append contents of one SN table to another
    
    inSN      = Input Obit SN Table
    outSN     = Output Obit SN Table
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Open
    inSN.Open(Table.READONLY, err)
    outSN.Open(Table.READWRITE, err)
    if err.isErr:
        return
    print "Input",inSN.Desc.Dict
    print "Output",outSN.Desc.Dict
    # Number of rows
    nrow =  inSN.Desc.Dict["nrow"]
    for i in range (0,nrow):    # Loop over rows
        SNrow = inSN.ReadRow(i+1, err)
        outSN.WriteRow(-1, SNrow, err)
        if err.isErr:
            return
    # end loop over rows
            
    # Close table
    inSN.Close(err)
    outSN.Close(err)
    if err.isErr:
        return
    
    # end VLBASNAppend 

def VLBADiagPlots( uv, err, cleanUp=True, JPEG=True, sources=None, project='',
    session='', band='', logfile=None, check=False, debug=False ):
    """
    Generate single source diagnostic plots.

    This method uses the averaged, calibrated data generated by the
    pipeline to produced single source diagnostics. The diagnostics
    can be used to assess the quality of the visibility data underlying
    the pipeline image maps and to assess the quality of the pipeline
    algorithms.

    uv = UV data object to plot
    err = Python Obit Error/message stack
    cleanUp = clean up temporary files when finished
    JPEG = if True, make JPEG plots; else make PS 
    sources = list of sources; None = make plots for each source
    logfile =  logfile for messages
    check = Only check script
    debug = Turn on debug mode
    """

    mess = "Generating diagnostic plots for each source"
    printMess(mess, logfile)

    avgClass = 'UVAvgT' # uv average temporary data
    avgSeq = 1
    
    # Average data over: 1 sec, all IFs, all channels
    calAvgTime = 1 # temporal averaging (sec)
    printMess("Averaging: "+str(calAvgTime)+" sec interval, all IFs, all channels",
        logfile = logfile)
    rtn = VLBACalAvg( uv, avgClass=avgClass, avgSeq=avgSeq, err = err, 
        logfile = logfile, check=check, debug = debug, CalAvgTime = calAvgTime, 
        avgFreq = 3, # avg all IFs
        chAvg = 0, # avg all channels (should already have been done)
        doCalib = 2, # apply calibration
        doBand = 0 # do not calibrate bandpass; already calibrated
        )
    if rtn != 0:
        mess = "Error averaging data. VLBACalAvg returned: " + str(rtn)
        printMess(mess, logfile)
        return rtn

    # Get the newly averaged data set: most recent file with class UVAvg
    uvAvg = None
    if not check:
        uvAvg = UV.newPAUV("AIPS UV DATA", uv.Aname, avgClass, uv.Disk, avgSeq, 
            True, err)

    # Put source list into slist
    if not sources:
        slist = VLBAAllSource(uvAvg,err,logfile=logfile,check=check,debug=debug)
    else:
        slist = sources
    if not type(slist) == list:
        slist = [slist]

    # Setup UVPLT
    uvplt = AIPSTask.AIPSTask("uvplt")
    if not check:
        setname(uvAvg, uvplt)
    uvplt.stokes = 'I' # unpolarized
    uvplt.ltype = -3 # Omit PL number and creation time
    printMess("Plotting stokes "+uvplt.stokes, logfile=logfile)

    # Setup LWPLA
    lwpla = AIPSTask.AIPSTask("lwpla")
    if not check:
        setname(uvAvg, lwpla)
    
    # Define plots: file => filename string, bparm => UVPLT 
    plotTypes =  ( { 'file' : 'amp', 'bparm' : [3,1]   },  # Amp vs. uv Dist
                   { 'file' : 'uv' , 'bparm' : [7,6,2] },  # u vs. v 
                   { 'file' : 'ri' , 'bparm' : [10,9]  } ) # Re vs. Im

    # Loop over sources
    for (i,s) in enumerate(slist):
        print "Generatig diagnostic plots for source "+s+ \
            " ("+str(i+1)+"/"+str(len(slist))+")"
        uvplt.sources[1] = s
        # Loop over plot types
        for plot in plotTypes:
            uvplt.bparm = AIPSTask.AIPSList( plot['bparm'] )
            if not check:
                uvplt.go()

            # Create output file name
            outfile = project+session+band+s+'.'+plot['file']+'.ps'
            lwpla.outfile = 'PWD:'+outfile # output to current directory

            # Remove preexisting file
            if os.path.exists(outfile): os.remove(outfile) 
    
            if not check:
                # Trap failure
                try:
                    if not check:
                        lwpla.g
                except Exception, exception:
                    print exception
                    mess = "Lwpla Failed - continuing anyway"
                    printMess(mess, logfile)
                else:
                    if JPEG:
                        # Convert PS to JPG
                        jpg = outfile[:-3]+'.jpg'
                        printMess('Converting '+outfile+' -> '+jpg,logfile)
                        cmd = 'convert -depth 96 '+outfile+' '+outfile[:-3]+'.jpg'
                        rtn = os.system(cmd)
                        if rtn == 0: 
                            if cleanUp: 
                                os.remove(outfile) # Remove the PS file
                            else:
                                # Print error message and leave the PS file
                                mess="Error occurred while converting PS to JPG"
                                printMess(mess,logfile)
    
    if cleanUp:
        printMess('Cleaning up temporary data',logfile)
        if not check:
            zap(uvAvg)

    # end VLBADiagPlot

def VLBAKntrPlots( err, catNos=[], imClass='?Clean', imName=[], project='tProj', 
    session='tSes', band='tB', debug=False ):
    """
    Create contour plots for the specified images. Image selection is made
    based on the input catalog numbers (catNos), or, if catalog numbers are not
    given, based on a pattern match to the image name and class. Pattern
    matching follows the rules of function AMcat(). One PS file is generated
    for each unique image name. Multiple images with the same name will be added
    to the same file on different pages. Arugments project, session, and band
    are used only in creating file names.

    err = Python Obit Error/message stack
    catNos = catalog numbers of images
    imClass = class of images to plot (used only if catNos is empty)
    imName = name of images to plot; None = make plots for each source (used 
             only if catNos is empty)
    project = project name
    session = project session
    band = project receiver band code
    debug = Turn on debug mode
    """
    # Setup AIPS task KNTR
    kntr = AIPSTask.AIPSTask("kntr")
    kntr.dogrey = 0
    kntr.dovect = 0
    kntr.ltype = -2 # Show border and labels w/o creation date and PL version
    kntr.cbplot = -18 # half-power beam in bottom right; no contour overlay
    # Set contour levels in units of cntr.clev (defined below). Contours begin 
    #   with -2, -2^0.5, 2^0.5, and then increase as powers of root two.
    levs = [ -2, -2**(0.5), 2**(0.5) ]
    for i in range(27):
        l = levs[-1] * 2.0**( 0.5 )
        levs = levs + [ l ]
    kntr.levs = AIPSTask.AIPSList( levs )

    # Instantiate AIPS task LWPLA
    lwpla = AIPSTask.AIPSTask("lwpla")    
    
    # If catalog numbers not given, get all images matching class imClass
    # and with names in list imName.
    if (not catNos):
        if not imName:
            imName = '*'
        elif not type(imName) == list:
            imName = [ imName ]
        for n in imName:
            catNos += AMcat( Aname=n, Aclass=imClass, giveList=True )
    for cno in catNos: # loop over images
        image = getname(cno)

        # Run CNTR to make plot
        setname(image, kntr)
        # Contour level unit = 2 * RMS noise
        stats = imstat(image, err)
        kntr.clev = 2 * stats['RMSHist']
        kntr.go()

        # Run LWPLA to make PS file
        setname(image, lwpla)
        name = image.Aname.rstrip() # Image name w/o whitespace
        outfile = project+session+band+name+'.cntr.ps'
        lwpla.outfile = './'+outfile # output to current directory
        # Trap failure
        try:
            if not check:
                lwpla.g
        except Exception, exception:
            print exception
            mess = "Lwpla Failed - continuing anyway"
            printMess(mess, logfile)
        else:
            pass
        # Delete plot files
        if not check:
            zz=image.ZapTable("AIPS PL", -1,err)
    # end VLBACntrPlots
