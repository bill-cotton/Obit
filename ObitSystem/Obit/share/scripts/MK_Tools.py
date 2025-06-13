# Package of tools for processing MeerKAT data
#exec(open('MK_Tools.py').read())
SetMFImage=None; ImageA2F=None; RMFitQU=None; MakeMask=None
HugoPerley_3C286_l2 =None; HugoPerley_3C286=None; PlotRM=None; PlotSpec=None
linRegW =None; RMClip=None

import UV, Image, OErr, RMFit, OSystem, Table, ObitTask
import ImageDesc, FArray, OPlot
from OTObit import setname, addParam
import math
from math import isnan

del SetMFImage
def SetMFImage(uv, src, gainUse=0, PDVer=-1, outDisk=1, outSeq=1,
               Niter=1000, minFlux=50.0e-6, ref=0, nThreads=16):
    """
    Setup Parameters for MFImage, generalized defaults for MeerKAT

    Parameters may need adjusting before execution
    * uv       = Python Obit UV object
    * src      = Source to be imaged
    * gainUse  = if >= 0 apply gain cal with gainUse (0->hi)
    * PDVer    = if >  0 apply poln cal with PDVer
    * outDisk  = AIPS disk number for images
    * outSeq   = AIPS sequence number for images
    * Niter    = Maximum number of iterations
    * minFlux  = Minimum flux density for first CLEAN
    * ref      = Reference antenna
    * nThreads = number of threads to use
    returns MFImage task object
    """
    mf=ObitTask.ObitTask('MFImage'); setname(uv,mf)
    mf.Sources=[src]; mf.flagVer=1; mf.Stokes='IQUV'
    mf.outDisk=outDisk; mf.outSeq=outSeq; mf.outClass='IPol'
    mf.out2Disk=mf.inDisk; mf.out2Seq=mf.inSeq; 
    mf.Niter=Niter; mf.Gain=0.05; mf.Robust=-1.5; 
    mf.nThreads=nThreads; 
    mf.prtLv=2; mf.logFile="MFImage_"+src+"."+str(outSeq)+".log"
    # Calibration
    if gainUse>=0:
        mf.doCalib=2; mf.gainUse=gainUse
    if PDVer>0:
        mf.doPol=True; mf.PDVer=PDVer; mf.keepLin=True
    # Others
    mf.BLFact=1.01; mf.BLchAvg=True  # Baseline dependent averaging
    mf.doFit= False;  mf.PBCor=False; mf.ccfLim=0.5
    mf.minFlux=minFlux
    # Selfcal
    mf.refAnt=ref; mf.avgPol=True; 
    mf.maxPSCLoop=2; mf.minFluxPSC=50.0e-6; mf.solPInt=0.5; mf.solPType='L1'; mf.solPMode='P'
    mf.maxASCLoop=1; mf.minFluxASC=0.5; mf.solAInt=10.0; mf.solAType='L1'; mf.solAMode='A&P'
    
    # Get frequency to set FOV, outliers
    mf.Catalog='AllSkyVZ.FIT'; mf.OutlierSize=520
    d=uv.Desc.Dict
    nu_0 = d['crval'][d['jlocf']]
    UHF=False; LBand=False; SBand=False
    UHF=nu_0<700e6
    LBand = (nu_0>700e6) and (nu_0<1400e6)
    SBand = nu_0>1400e6
    if UHF:
        mf.FOV=1.5; mf.OutlierDist=2.0
    if LBand:
        mf.FOV=1.0; mf.OutlierDist=1.5
    if SBand:
        mf.FOV=0.75; mf.OutlierDist=1.0

    # Additions to input
    NiterQU = Niter//2
    addParam(mf,"NiterQU", paramVal=NiterQU, shortHelp="Niter for Q,U", \
             longHelp="  NiterQU.....Niter for Q,U\n")
    minFluxQU = 20.0e-6
    addParam(mf,"minFluxQU", paramVal=minFluxQU, \
             shortHelp="minFlux for Q,U", \
             longHelp="  minFluxQU...minFlux for Q,U\n")

    NiterV=1000; minFluxV=0.00005 # Don't need much V CLEANing
    addParam(mf,"NiterV", paramVal=NiterV, shortHelp="Niter for V", \
             longHelp="  NiterV.....Niter for V\n")
    addParam(mf,"minFluxV", paramVal=minFluxV, \
             shortHelp="minFlux for V", \
             longHelp="  minFluxV...minFlux for V\n")
    maxAWLoop = 1 # Seems needed for MeerKAT
    addParam(mf,"maxAWLoop", paramVal=maxAWLoop, shortHelp="Max. middle CLEAN loop", \
             longHelp="  maxAWLoop....Max. middle CLEAN Loop count\n"+ \
             "              Override the default behavior for the autoWin middle loop\n"+
             "              if > 0.\n")
    minFList = [0.0001,0.00005]  # 0.0001 Jy after 1st, 0.00005 after 2nd
    addParam(mf,"minFList", paramVal=minFList, \
             shortHelp="minFlux list after SC", \
             longHelp="  minFList....minFluxes to use in IPol after selfcals\n")
    return mf


# end SetMFImage

del ImageA2F
def ImageA2F (img, err, dir="./"):
    """
    Write an AIPS Image to a disk file

    * img      = Obit/python AIPS image 
    * err      = Obit error/message stack
    * dir      = Directory name for FITS image
    returns obit Image object to FITS file
    """
    from OTObit import imtab
    f = dir+'/'+img.Aname.strip()+"_"+img.Aclass.strip()+ \
        "_"+str(img.Aseq)+".fits"
    out = imtab(img,f,0,err)
    OErr.printErr(err)
    return out

# end ImageA2F
  
del RMFitQU
def RMFitQU (name, imQ, imU, err, doRMSyn=True, nThreads=16, maxRM=100.):
    """
    Faraday synthesis of a Q/U image cube pair

    See help(RMFit.Cube) for details
    * name     = root of output file name
    * imQ      = Obit/python Q cube
    * imU      = Obit/python U cube
    * err      = Obit error/message stack
    * doRMSyn  = If True use peak of RM spectrum, 
                 else least squares
    * nThreads = Number of threads to use
    * maxRM    = search +/- maxRM rad/m^2
    returns RM Obit Image object
    """
    OSystem.PAllowThreads(nThreads)  # Threading
    OErr.printErr(err)
    fit = RMFit.PCreate('Fitter')
    fit.List.set('doError',True);    fit.List.set('doRMSyn',doRMSyn)
    fit.List.set('minRMSyn',-maxRM); fit.List.set('maxRMSyn',+maxRM)
    fit.List.set('delRMSyn',0.5);    fit.List.set('maxChi2',10000.0)
    fit.List.set('minQUSNR',1.0);    fit.List.set('minFrac',0.25);   
    fit.List.set('refLamb2',1.0e-6);  
    # which method?
    if doRMSyn:
        outRM = Image.newPFImage('RM', name+'_RMSyn.fits',0,False,err)
    else:
        outRM = Image.newPFImage('RM', oname+'_RMLSQ.fits',0,False,err)
    fit.Cube(imQ, imU, outRM, err)
    OErr.printErr(err)
    del fit
    return outRM

# end RMFit
  
del RMClip
def RMClip (imRM, err, minPPol=50.0e-6 ):
    """
    Clip RM and EVLA planes in RMFit cube by ppol

    WARNING: This routine modifies the input image;
    ONLY RUN IT ONCE.
    * imRM     = Obit/python Image as produced by RMFitQU
    * err      = Obit error/message stack
    * minPPol  = min PPol (inRM plane 3) in Jy
    """
    # Check
    if not imRM.ImageIsA():
        raise TypeError("input MUST be a Python Obit Image")
    # Image info
    d = imRM.Desc.Dict
    (nx,ny) = d['inaxes'][0:2]
    # Read planes
    RMPix = FArray.FArray('RM', [nx,ny])
    imRM.GetPlane(RMPix, [1,1,1,1,1], err)
    OErr.printErrMsg(err, "Error reading image")
    EVPAPix = FArray.FArray('EVPA', [nx,ny])
    imRM.GetPlane(EVPAPix, [2,1,1,1,1], err)
    PPolPix = FArray.FArray('PPol', [nx,ny])
    imRM.GetPlane(PPolPix, [3,1,1,1,1], err)
    OErr.printErrMsg(err, "Error reading image")
    # Make mask in PPolPix
    FArray.PInClip(PPolPix, -1.0e6, minPPol, FArray.fblank)
    # Blank
    FArray.PBlank(RMPix,   PPolPix, RMPix)
    FArray.PBlank(EVPAPix, PPolPix, EVPAPix)
    # Rewrite
    imRM.PutPlane(RMPix,   [1,1,1,1,1], err)
    imRM.PutPlane(EVPAPix, [2,1,1,1,1], err)
    OErr.printErrMsg(err, "Error writing image")

    # Add history
    outHistory = History.History("history", imRM.List, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit RMClip",err)
    outHistory.WriteRec(-1,"RMClip"+" minPPol = "+str(minPPol),err)
    outHistory.Close(err)
# end RMClip
  
del MakeMask
def MakeMask (img, err, name=None, sigma=10.0, nThreads=10, maxWid=60., minSNR=5.):
    """
    Make a CLEAN mask from an image

    Generates a <?>.mask file in cwd which can be used as 
    CLEANFile in MFImage or Imager.  Also, a test *.cat file
    Uses Obit/FndSou for source finding.
    NB: Any previous MF or VL tables deleted
    * img      = Obit/python Stokes I image, AIPS or FITS
    * err      = Obit error/message stack
    * name     = if given, the root mask, <name>.mask
                 if none AIPS or FITS file name
    * sigma    = multiple of RMS to search for sources.
    * nThreads = Number of threads to use
    * maxWid   = Maximum component size in asec
    * minSNR   = Minimum SNR of component
    """
    # Output file root
    ftype = img.FileType # Input file type
    if name:
        froot = name
    elif ftype=='FITS':
        froot = img.FileName
        # Remove ".fits" if possible
        e = froot.find(".fits")
        if e>0:
            froot=froot[:e]
    else:
        froot = img.Aname.strip()
    print(" Root of output files names is ",froot)
    OSystem.PAllowThreads(nThreads)    # Threading
    s = imstat(img); rms=s['RMSHist']  # Image statistics
    print ("Searching to ",rms*sigma)
    # Clear any previous MF, VL tables
    # If FITS need to refresh header
    if ftype=='FITS':
        img = Image.newPFImage("cat",img.FileName, img.Disk, True, err)
    img.ZapTable("AIPS MF",-1,err); img.ZapTable("AIPS VL",-1,err); 
    # Generate source list
    fs=ObitTask.ObitTask("FndSou"); setname(img,fs)
    fs.doVL=True; fs.NGauss=10000; fs.NPass=1; fs.RMSsize=50
    fs.prtLv=1; fs.nThreads=nThreads
    fs.CutOff=0.9*rms*sigma; fs.doMult=True; fs.doWidth=True
    fs.Parms=[sigma*rms, maxWid, 0., 3., 0.] 
    fs.OutPrint=froot+".cat"
    fs.g
    # If FITS need to refresh header
    if ftype=='FITS':
        img2 = Image.newPFImage("cat",img.FileName, img.Disk, True, err)
    # Generate mask from VL table
    maskFile = froot+".mask"; vlver = 1
    fd = open(maskFile,'w')
    vltab=img2.NewTable(Table.READONLY,'AIPS VL',vlver,err)
    maxscale = 1.0/abs(img2.Desc.Dict['cdelt'][0])
    nrow=vltab.Desc.Dict['nrow']
    vltab.Open(Table.READONLY,err)
    OErr.printErr(err)
    for irow in range(1,nrow+1):
        row = vltab.ReadRow(irow,err)
        if row['PEAK INT'][0]/row['I RMS'][0]>minSNR:
            ra =  ImageDesc.PRA2HMS(row['RA(2000)'][0])
            dec = ImageDesc.PDec2DMS(row['DEC(2000)'][0])
            size = min(50, int(row['MAJOR AX'][0]*maxscale+0.5))
            line = "%s %s %d\n"%(ra,dec,size)
            fd.write(line)
    
    vltab.Close(err)
    OErr.printErr(err)
    fd.close()
# end MakeMask

import math
from math import isnan
del HugoPerley_3C286
def HugoPerley_3C286 (nu):
    """
    Calculate polarization model 
    
    from SARAO memo SSA-0004E-001 B
    returns (P [fract lin. pol.], EVPA[deg])
    nu = frequency (Hz)
    """
    clight = 2.997924562e8   # Speed of light m/s
    lamb = (clight/nu); lamb2 = lamb**2   # Wavelength and squared
    nug = nu*1.0e-9          # GHz
    P = 10.0; EVPA = 33.0  # Default
    # EVPA By frequency range
    if (nug>=1.7) and (nug<=12.0):
        EVPA = 32.64 + 85.37*lamb2
    elif nug<1.7:
        c1 = math.log10(nug)**3
        EVPA = 29.53 + lamb2 *  (4005.88 *c1 - 39.38)
    # P By frequency range
    c2 =  math.log10(lamb2)
    if (nug>=1.1) and (nug<=12.0):
        P = 0.080 - 0.053*lamb2 - 0.015*c2
    elif nug<1.1:
        P = 0.029- 0.172*lamb2 - 0.067*c2
    return (P,EVPA)
# end HugoPerley_3C286 

del HugoPerley_3C286_l2 
def HugoPerley_3C286_l2 (lamb2):
    """
    Calculate polarization model 
    
    from SARAO memo SSA-0004E-001 B
    returns (P [fract. lin. pol.], EVPA[deg])
    lamb2 = wavelength^2 in m^2
    """
    lamb = lamb2**0.5
    clight = 2.997924562e8   # Speed of light m/s
    nu = clight/lamb
    nug = nu*1.0e-9          # GHz
    P = 10.0; EVPA = 33.0  # Default
    # EVPA By frequency range
    if (nug>=1.7) and (nug<=12.0):
        EVPA = 32.64 + 85.37*lamb2
    elif nug<1.7:
        c1 = math.log10(nug)**3
        EVPA = 29.53 + lamb2 *  (4005.88 *c1 - 39.38)
    # P By frequency range
    c2 =  math.log10(lamb2)
    if (nug>=1.1) and (nug<=12.0):
        P = 0.080 - 0.053*lamb2 - 0.015*c2
    elif nug<1.1:
        P = 0.029- 0.172*lamb2 - 0.067*c2
    return (P,EVPA)
# end HugoPerley_3C286_l2

del PlotRM
def PlotRM (src, icube, qcube, ucube, pos, plotfile, err,
           label=None,  nThreads=1):
    """
    Make rotation measure plot, EVPA v lambda^2

    Make RM and fractional poln. plots at position pos
    * src     = Source name
    * icube   = Obit IPol image cube
    * qcube   = Obit QPol image cube
    * ucube   = Obit UPol image cube
    * pos     = (RA, dec) in deg. Uses nearest pixel.
    * label   = label for plot, defaults to src
    * plotfile= name of postscript plot file
    * err     = Python Obit Error/message stack
    * nthreads= the number of threads to be used in calculations
    """
    OSystem.PAllowThreads(nThreads)  # threading
    HPAlias = ["3C286","J1331+3030"] # Hugo/Perley Aliases
    fblank = FArray.fblank
    if label==None:
        label = src

    # Perley/Butler values (3c286 wrong - use Hugo/Perley)
    PB={"3C48": [(1.050, 0.003, 25.),(1.45, 0.005, 140.), (1.64, 0.007, -5.0)], \
        "3C138":[(1.050, 0.056, -14.),(1.45, 0.075, -11.), (1.64, 0.084, -10.0)], \
        "3C286":[(1.050, 0.086, 33.),(1.45, 0.095, 33.), (1.64, 0.099, 33.)], \
        "J0521+1638_S":[(1.050, 0.056, -14.),(1.45, 0.075, -11.), (1.64, 0.084, -10.0),
                        (1.95,9.0,-10.), (2.45,10.4,-9.), (2.95,10.7,-10.), (3.25,10.0,-10.)], \
        "J1331+3030_S":[(1.050, 0.086, 33.),(1.45, 0.095, 33.), (1.64, 0.099, 33.),
                  (1.95,10.1,33.), (2.45,10.5,33.), (2.95,10.8,-33.), (3.25,10.9,-33.)], \
    }
    
    # Convert position to pixel
    pp   = ImageDesc.PGetPixel(icube.Desc, pos, err)
    OErr.printErrMsg(err,message='Determining pixel')
    pixel = [int(pp[0]+0.5), int(pp[1]+0.5)]  # nearest pixel
    # get lambda sq
    d=icube.Desc.List.Dict
    clight=2.997924562e8
    nterm = d['NTERM'][2][0]
    nspec = d['NSPEC'][2][0]
    lamb2 = []
    for ip in range(nterm+1,nterm+nspec+1):
        key = "FREQ%4.4d"%(ip-nterm)
        frq = d[key][2][0]
        lamb2.append((clight/frq)**2)
        #print ip,frq
    #
    frq0 = d['FREL0001'][2][0] # ref freq
    lamb20 = (clight/frq0)**2
    #print('lambdas',lamb2[0],lamb20)
    #pixel = s[1]  # Pixel to plot
    icube.Open(Image.READONLY,err)
    qcube.Open(Image.READONLY,err)
    ucube.Open(Image.READONLY,err)
    iflux=[]; iRMS=[]  # Ipol flux, RMS in plane
    qflux=[]; qRMS=[]  # Qpol flux, RMS in plane
    uflux=[]; uRMS=[]  # Upol flux, RMS in plane
    llamb2 = []; nplot = 0
    plane = [1,1,1,1,1]
    for ip in range(nterm+1,nterm+nspec+1):
        plane[0] = ip
        icube.GetPlane(None, plane, err)
        if icube.FArray.get(pixel[0]-1, pixel[1]-1)!=0.0:
            llamb2.append(lamb2[ip-nterm-1])
            nplot += 1  # How many valid data points
            iRMS.append(icube.FArray.RMS)
            iflux.append(icube.FArray.get(pixel[0]-1, pixel[1]-1))
            qcube.GetPlane(None, plane, err)
            qRMS.append(qcube.FArray.RMS)
            qflux.append(qcube.FArray.get(pixel[0]-1, pixel[1]-1))
            ucube.GetPlane(None, plane, err)
            uRMS.append(ucube.FArray.RMS)
            uflux.append(ucube.FArray.get(pixel[0]-1, pixel[1]-1))
    
    icube.Close(err)
    qcube.Close(err)
    ucube.Close(err)
    OErr.printErrMsg(err,message='Error reading data')

    # Get quanties to plot and errors
    phs=[];  ephs=[]; l2=[]; nsamp = 0
    fpol=[]; efpol=[]; qqflux=[]; uuflux=[]; qqRMS=[]; uuRMS=[]
    for ip in range(0,nplot):
        # Datum OK?
        OK = (not isnan(iflux[ip])) and (not isnan(qflux[ip])) and (not isnan(uflux[ip])) and \
             (not (iflux[ip]==fblank) and (not (qflux[ip]==fblank)) and (not (uflux[ip]==fblank)) \
              and (not (iflux[ip]==0.0)) and (not (qflux[ip]==0.0)))
        if OK:
            nsamp += 1
            l2.append(lamb2[ip])
            phs.append(0.5*math.degrees(math.atan2(uflux[ip],qflux[ip])))
            arg    = abs(uflux[ip]/qflux[ip])
            sigarg = arg*((qRMS[ip]**2/qflux[ip]**2)+(uRMS[ip]**2/uflux[ip]**2))**0.5
            ephs.append(abs(math.degrees(sigarg/(1+arg**2))))
            arg = ((uflux[ip]**2+qflux[ip]**2)**0.5)/iflux[ip]
            fpol.append(100.*arg)
            vararg  = (qRMS[ip]**2 +uRMS[ip]**2)
            sigarg  = arg*((vararg/arg**2) + iRMS[ip]**2/iflux[ip]**2)**0.5
            efpol.append(100.*sigarg)
            ip, iflux[ip], (uflux[ip]**2+qflux[ip]**2)**0.5, 0.5*math.degrees(math.atan2(uflux[ip],qflux[ip]))
            qqflux.append(qflux[ip]); qqRMS.append(qRMS[ip])
            uuflux.append(uflux[ip]); uuRMS.append(uRMS[ip])
    
    # Remove phase wraps
    for ip in range(1,nsamp):
        if phs[ip]-phs[ip-1]>130.:
            phs[ip] -=180.
        if phs[ip]-phs[ip-1]<-130.:
            phs[ip] +=180.
        
    # Fit RM
    refLamb2 = lamb2[nspec//2]
    RMParms = RMFit.PSingle(refLamb2, l2, qqflux, qqRMS, uuflux, uuRMS, err)
    evpa0 =  math.degrees(RMParms[1]+RMParms[0]*(lamb20-refLamb2)) # EVPA @ ref. freq
    eevpa0 = math.degrees(RMParms[3])
    print ("RMParms",RMParms)
    #print (src, math.degrees(RMParms[1]), lamb20, evpa0, eevpa0)
    #print ('phs',phs)
    # Guess offset (180 ambiguity) from difference in evpa and phs[0]
    pdif =  phs[0]-math.degrees(RMParms[1])
    evpa0 += round(pdif/180)*180
    #print (math.degrees(RMParms[1]),phs[0],pdif)

    plot=OPlot.newOPlot("Plot", err,output=plotfile+".ps/ps",ny=2)
    # In Hugo & Perley 2024 (3C286)?
    HPP = []; HPEVPA = []
    if src in HPAlias:
        for ll2 in l2:
            (hpP,hpEVPA) = HugoPerley_3C286_l2 (ll2)
            HPP.append(100*hpP);  HPEVPA.append(hpEVPA); 
            
    # Plot EVPA
    plot.List.set("TITLE"," %s: RM=%6.2f (%5.3f), EVPA#d_ref#u =%7.2f (%5.2f)"\
                  %(label, RMParms[0], RMParms[2],evpa0,eevpa0))
    plot.List.set("YLABEL","EVPA (deg)")
    plot.List.set("XLABEL","#gl#u2#d (m#u2#d)")
    plot.List.set("XMAX",1.05*lamb20)
    plot.List.set("XMIN",0.95*lamb2[nspec-1])
    avgP = sum(phs)/len(phs)
    #print ("avgP",avgP)
    plot.List.set("YMIN",(min(min(phs)-15,avgP-20.)))
    plot.List.set("YMAX",(max(max(phs)+15,avgP+20.)))
    plot.List.set("CSIZE",1)
    plot.List.set("SSIZE",3)
    plot.List.set("LWIDTH",3)
    OPlot.PXYErr(plot, 2, l2, phs, ephs, err)
    #print ("Phs", phs); print ("ephs",ephs)
    # In Hugo& Perley 2024
    if src in HPAlias:
        OPlot.PXYOver(plot, 8, l2, HPEVPA, err)
        #print (HPEVPA); print (HPP)
        #print ("HP Phs", phs); print ("HP FPol",fpol)
    # In Perley&Butler?
    elif src in PB:
        for dd in PB[srs]:
            x1 = ((clight/(dd[0]*1.0e9))**2)
            y1 = dd[2];
            #print "PB", dd,x1,y1
            OPlot.PDrawSymbol(plot,x1,y1,8,err)
    # Plot RM fit - a few times
    noff = 7
    offs = [0.0, 180, -180., 360., -360, 540. -540.]
    for off in offs:
        x1 = 1.1*lamb2[0] - refLamb2
        y1 = math.degrees(RMParms[1]+RMParms[0]*x1)+off
        x2 = 0.9*lamb2[nspec-1] - refLamb2
        y2 = math.degrees(RMParms[1]+RMParms[0]*x2)+off
        OPlot.PDrawLine(plot, 1.1*lamb2[0], y1, 0.9*lamb2[nspec-1], y2, err)
    OPlot.PSetLineStyle(plot, 3, err)
    OPlot.PDrawLine(plot, lamb20, -2000., lamb20, 2000., err)
    OPlot.PSetLineStyle(plot, 1, err)
    # Plot fractional poln
    plot.List.set("TITLE"," ")
    plot.List.set("XLABEL","#gl#u2#d (m#u2#d)")
    plot.List.set("YLABEL","frac. pol. (%)")
    plot.List.set("YMIN",0.0)
    plot.List.set("YMAX",min(1.20*max(fpol),50.0))
    OPlot.PXYErr(plot, 2, lamb2, fpol, efpol, err)
    #print ("fpol", fpol); print ("efpol",efpol)
    # In Perley&Butler?
    # In Hugo& Perley 2024?
    if src in HPAlias:
        OPlot.PXYOver(plot, 8, l2, HPP, err)
    # Else In Perley&Butler?
    elif src in PB:
        for dd in PB[src]:
            x1 = ((clight/(dd[0]*1.0e9))**2)
            y1 = 100.*dd[1];
            OPlot.PDrawSymbol(plot,x1,y1,8,err)
    OPlot.PShow(plot,err)
    OErr.printErrMsg(err,message='Plotting Rotation Measure')

# end PlotRM

del linRegW 
def linRegW (nu, s, e):
    """
    Weighted linear regression for 2 term spectral fitting
    
    returns (flux@freq[0], spectral_index)
    * nu    = list of frequencies 
    * s     = list of flux densities
    * e     = list of flux density errors
    """
    from math import log, exp
    s_w=0; s_x=0.; s_y=0.; s_xx=0.; s_yy=0.; s_xy=0.
    nu_0 = nu[0]; s_0 = s[0]
    lnu = []; ls = []; w = []
    for i in range(0,len(nu)):
        if s[i] and s[i]>0:
            lnu.append(log(nu[i]/nu_0))
            ls.append(log(s[i]))
            w.append(1/e[i])
    n = len (lnu)
    for j in range(0,n):
        s_x += w[j]*lnu[j]; s_y += w[j]*ls[j]; s_xy += w[j]*lnu[j]*ls[j]
        s_xx += w[j]*lnu[j]**2; s_yy += w[j]*ls[j]**2; s_w+=w[j]
    a = (s_w*s_xy-s_x*s_y)/(s_w*s_xx-s_x*s_x)
    b = (s_y-a*s_x)/s_w
    #print (a, b, "S_0=", exp(b), "alpha=", a)
    return (exp(b),a)
# end linRegW

del PlotSpec
def PlotSpec (src, icube, pos, plotfile, err,
              label=None,  nThreads=1):
    """
    Plot the spectrum at a location in an image
    
    * src     = Source name
    * icube   = Obit IPol image cube
    * pos     = (RA, dec) in deg.
    * label   = label for plot, defaults to src
    * plotfile= name of postscript plot file
    * err     = Python Obit Error/message stack
    * nthreads= the number of threads to be used in calculations
    """
    OSystem.PAllowThreads(nThreads)  # threading
    import math
    from math import isnan, log
    import FITS, OPlot, Image, ImageDesc, ImageMF, FArray, FInterpolate, ImageUtil, SpectrumFit
    fblank = FArray.fblank
    # Image info
    nterm = icube.Desc.List.Dict['NTERM'][2][0]
    nspec = icube.Desc.List.Dict['NSPEC'][2][0]
    # Interpolator
    fi = FInterpolate.FInterpolate('FI', icube.FArray, icube.Desc, 2)
    # Get pixel location
    pixel = ImageDesc.PGetPixel(icube.Desc, pos, err)
    # end loop
    OErr.printErrMsg(err,message='Finding positions')
    # Loop over planes interpolating values
    freqs = []; rms = []; vals = []
    for i in range(0,nspec):
        # Read plane
        plane=[i+nterm+1,1,1,1,1]
        Image.PGetPlane (icube, None, plane, err)
        rmss = icube.FArray.RMS
        if not isnan(rmss): # Check validity
            # Interpolate position
            val = FInterpolate.PPixel(fi, pixel, err)
            #print (i,j,val,rms[i],freqs[i])
            if val>0.0:
                vals.append(val)
                key = 'FREQ%4.4d'%(i+1)
                freqs.append(1.0e-9*icube.Desc.List.Dict[key][2][0]) # GHz
                rms.append(rmss) # Plane RMS
                
                
    # end loop over planes
    OErr.printErrMsg(err,message='Interpolating values')

    # Plot
    plot=OPlot.newOPlot("Plot", err,output=plotfile+".ps/ps")
    # Fit Spectrum
    nterm=2; refFreq = freqs[0]
    # Inconsistent results
    #FitParms = SpectrumFit.PSingle(nterm,refFreq, freqs, vals, rms ,err)
    #if FitParms[1]==fblank:  # fit failed?
    #    FitParms[1] = 0.0
    FitParms = linRegW (freqs, vals, rms)
    #print ("\nFit",FitParms)
    #print ("vals",vals)
    #print ("rms",rms)
    ras = ImageDesc.PRA2HMS(pos[0])[:13]; decs = ImageDesc.PDec2DMS(pos[1])[:11];
    #label = "%s %s %s s=%6.4f(%8.5f) #ga=%5.2f(%4.2f)"\
    #        %(src, ras[1:], decs, FitParms[0], FitParms[2],FitParms[1], FitParms[3])
    # PLPlot can't handle the longer label
    label = "%s %s %s s=%6.4f #ga=%5.2f"\
            %(src, ras[1:], decs, FitParms[0], FitParms[1])
    plot.List.set("TITLE",label)
    plot.List.set("YLABEL","Flux density (Jy)")
    plot.List.set("XLABEL","Frequency (GHz)")
    plot.List.set("XOPT","LBCN")
    plot.List.set("YOPT","LBCN")
    ymn = 0.9*10**int(-0.9+math.log10(min(vals)))
    plot.List.set("YMIN",ymn)
    plot.List.set("CSIZE",1)
    plot.List.set("SSIZE",3)
    plot.List.set("LWIDTH",3)
    OPlot.PXYErr(plot, 2, freqs, vals, rms, err)
    # Plot fitted spectrum
    x0 = min(freqs); y0 = FitParms[0]*(x0/refFreq)**FitParms[1]
    x1 = max(freqs); y1 = FitParms[0]*(x1/refFreq)**FitParms[1]
    OPlot.PDrawLine(plot, x0, y0, x1, y1, err)
    OErr.printErrMsg(err,message='Plotting Spectrum')
    OPlot.PShow(plot,err)
    OErr.printErrMsg(err,message='Plotting Spectrum')

# end PlotSpec
