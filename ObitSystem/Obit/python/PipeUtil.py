import ObitTask, Image, AIPSDir, OErr, FArray
import urllib, urllib2, os.path, pickle

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
    """ 
    Print message, optionally in logfile
        
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
   
def QueryArchive(startTime, endTime):
    """
    Query the NRAO Archive for data files.

    startTime = Start of query time range ( YYYY-mmm-DD [ HH:MM:SS ] )
    endTime = End of query time range
    """
    freqRange = '1000 - 16000' # frequency range (MHz) (L to U bands)
    format = 'FITS-AIPS' # Archive file format (FITS-AIPS good prior to ? 2010)
    # HTTP POST parameters
    dataList = [ ('PROTOCOL','TEXT-stream'),
                 ('TELESCOPE','VLBA'),
                 ('QUERYTYPE','ARCHIVE'),
                 ('OBSFREQ1', freqRange ),
                 ('ARCHFORMAT', format),
                 ('TIMERANGE1', startTime),
                 ('TIMERANGE2', endTime) ]
    data = urllib.urlencode( dataList )
    url = 'https://archive.nrao.edu/archive/ArchiveQuery'
    # Archive is known to occasionally send a null response. This is an error.
    # At the least, a header line should be returned. Repeat the query 
    # until a non-null response is obtained. If the response contains a header
    # only, repeat the query once to verify.
    lineCount = 0
    while 1:
        response = urllib2.urlopen( url, data ) # Submit query
        lines = response.readlines() # Extract response into a list of lines
        if len(lines) == 0: # null response
            print "Archive response is null. Repeating query."
            lineCount = len(lines)
        else: 
            lines.pop() # remove tailing blank line
        if len(lines) == 1: # Header only, no data
            if lineCount != 1:
                print "Archive response contains no files. " + \
                    "Repeating query to verify."
                lineCount = len(lines)
            else:
                print "Verified: Archive response contains no files."
                break
        elif len(lines) > 1: # Data is present
            break
    return lines

def ParseArchiveResponse( responseLines ):
    """
    Parse the archive response returned by QueryArchive and return a list
    containing one dictionary for each file. The values in each dictionary will
    correspond to items in each row of the query response.  The keys will 
    correspond to the items in the response header row, with additional column
    headers added where needed. Redundant columns are not included in the 
    dictionaries.

    responseLines = Archive query response returned by QueryArchive
    """
    head = responseLines[0] # header line
    headers = head.split(' ')
    headers.pop(0) # remove leading '#'
    headers.pop(1) # remove extra space
    # Add headers to columns where there are none
    headers[14:] = ['DATE', 'RESPONSE_ROW', '?', 'FILESIZE_UNIT']
    fileList = []
    for line in responseLines[1:]:
        metadata = line.split(',')
        # Remove redundant metadata
        metadata[14:19] = [] 
        metadata.pop(17)
        dict = {}
        for i,d in enumerate(metadata): # iterate over every column
            dict[ headers[i] ] = d.strip()
        fileList.append( dict )
    return fileList

def DownloadArchiveFile( fileDict, email, destination ):
    """
    Download a file from the archive.

    fileDict = a single file dictionary returned in the response from 
               ParseArchiveResponse
    email = Email address to receive download confirmation
    destination = Download destination directory
    """
    FTPCHECKED=fileDict['file_root']+','+ \
               fileDict['logical_file']+','+ \
               fileDict['raw_project_code']+','+ \
               fileDict['segment']+','+ \
               fileDict['lock_status']+','+ \
               fileDict['DATE']+','+ \
               fileDict['RESPONSE_ROW']+','+ \
               fileDict['?']+','+ \
               fileDict['format']+','+ \
               fileDict['FILESIZE_UNIT']
    # This IP lookup method may be host dependent. Since the archive accepts
    # the request without parameter REMOTEADDR, leave this off for now.
    # localhostIP =  \
    #     commands.getoutput("/sbin/ifconfig").split("\n")[1].split()[1][5:]
    dataList = [
        ('EMAILADDR', email ),
        ('FTPCHECKED', FTPCHECKED ),
        ('COPYFILEROOT', destination ),
        ('ARCHIVEROOT', '/home/archive/e2e/archive'),
        ('CONVERT2FORMAT', 'UVFITS'),
        ('DOWNLOADFTPCHK', 'Download Checked Files'),
        ('FILENAMESTYLE', 'ANALYSTS'),
        ('FLAGGING', 'FLAG'),
        ('PROTOCOL', 'HTML'),
        ('SPECTAVG', 'x1'),
        ('TIMEAVG','0s'),
        ('USERSTRING', '')
    ]
    data = urllib.urlencode( dataList )
    url = 'https://archive.nrao.edu/cgi-bin/e2eftp.cgi'
    response = urllib2.urlopen( url, data )
    lines = response.readlines()
    print "--- Archive Download Response ---"
    for ln in lines: print ln,
    print "--- --- "

def SummarizeArchiveResponse( fileList ):
    """
    Print a table summary of the archive response.

    fileList = List of file dictionaries returned by ParseArchiveResponse
    """
    formatHead = "%-2s %-6s %-1s %-18s %-18s %-6s"
    print formatHead % \
        ( "#-", "PCODE-", "S", "STARTTIME---------", "STOPTIME----------", 
        "SIZE--" )
    for i,file in enumerate(fileList):
        formatStr = "%2d %6s %1s %18s %18s %6s" 
        print formatStr % ( i, file['project_code'], 
            file['segment'], file['starttime'], file['stoptime'], 
            file['FILESIZE_UNIT'] )
