# Do instrumental polarization cal for MeerKAT L band Datasets
# Uses PCal w/ upolarized calibrators and
# XYDly for known polarized calibrators
# Generates PD and SN (+CL) tables
# write uvtab file with BP and PD tables
# Verbose log file from PCal and XYDly in 'PolnCal_'+uv.{AF}name.strip()+'.log'
#exec(open('MK_PolCal_XY.py').read())
MKPolCalXY = None; SwapPDXY=None; SwapSNXY=None; PDNeg=None; CorrPD=None
import ObitTask, OErr
err = OErr.OErr()
from OTObit import setname, imhead

del SwapPDXY
def SwapPDXY(uv, err, inVer=1, outVer=2):
    """
    Copy PD table inVer to outVer swapping X and Y solutions

    * uv       = Python Obit UV object, AIPS or FITS
    * err      = Python Obit Error/message stack
    * inVer    = input PD version
    * outVer   = output PD version
    """
    from math import pi
    import FArray, Table, History
    fblank=FArray.fblank
    uv.Header(err)
    # Copy table
    print ("SwapPDXY inVer=",inVer,"outVer=", outVer)
    pdtabi=uv.NewTable(Table.READONLY,"AIPS PD",inVer,err)
    pdtabi.Open(Table.READONLY,err)
    numIF=pdtabi.keys['numIF']; numPol=pdtabi.keys['numPol'];
    numChan=pdtabi.keys['numChan']; 
    pdtabo=Table.PClone(pdtabi,None)
    pdtabo=uv.NewTable(Table.WRITEONLY,"AIPS PD",outVer,err,
                       numIF=numIF,numPol=numPol,numChan=numChan)
    pdtabo.Open(Table.WRITEONLY,err)
    # Patch things in header which may get lost
    pdtabo.keys['numAnt']=pdtabi.keys['numAnt']
    pdtabo.keys['polType']=pdtabi.keys['polType']
    # 'REAL 1', 'IMAG 1','REAL 2', 'IMAG 2', 'P_DIFF'
    inrow = pdtabi.Desc.Dict['nrow']; nchan= pdtabi.Desc.Dict['repeat'][5];
    for irow in range(1,inrow+1):
        row=pdtabi.ReadRow(irow,err)
        for i in range(0,nchan):
            # Swap
            tre1 = row['REAL 1'][i]; tim1 = row['IMAG 1'][i]
            tre2 = row['REAL 2'][i]; tim2 = row['IMAG 2'][i]
            row['REAL 1'][i]=tre2;  row['IMAG 1'][i]=tim2
            row['REAL 2'][i]=tre1;  row['IMAG 2'][i]=tim1
            if  row['P_DIFF'][i]!=fblank:
                row['P_DIFF'][i] = -row['P_DIFF'][i]
        pdtabo.WriteRow (irow, row, err)
    pdtabi.Close(err); pdtabo.Close(err)
    OErr.printErrMsg(err, "Error Swapping PD tables")
    # History entry
    outHistory = History.History("history", uv.List, err)
    z=History.POpen(outHistory, History.READWRITE, err)
    z=History.PTimeStamp(outHistory," Start Obit SwapPDXY",err)
    z=History.PWriteRec(outHistory,-1,"SwapPDXY / inPD = "+str(inVer),err)
    z=History.PWriteRec(outHistory,-1,"SwapPDXY / outPD = "+str(outVer),err)
    z=History.PClose(outHistory, err)
    OErr.printErrMsg(err, "Error with history")

# end SwapPDXY

del SwapSNXY
def SwapSNXY(uv, err, inVer=1, outVer=2):
    """
    Copy SN table inVer to outVer swapping X and Y solutions

    * uv       = Python Obit UV object, AIPS or FITS
    * err      = Python Obit Error/message stack
    * inVer    = input PD version
    * outVer   = output PD version
    """
    from math import pi
    import FArray, Table, History
    fblank=FArray.fblank
    uv.Header(err)
    # Copy table
    print ("SwapSNXY inVer=",inVer,"outVer=", outVer)
    sntabi=uv.NewTable(Table.READONLY,"AIPS SN",inVer,err)
    sntabi.Open(Table.READONLY,err)
    numIF=sntabi.keys['numIF']; numPol=sntabi.keys['numPol'];
    sntabo=Table.PClone(sntabi,None)
    sntabo=uv.NewTable(Table.WRITEONLY,"AIPS SN",outVer,err,
                       numIF=numIF,numPol=numPol)
    sntabo.Open(Table.WRITEONLY,err)
    # Patch things in header which may get lost
    sntabo.keys['numAnt']=sntabi.keys['numAnt']
    # 'REAL 1', 'IMAG 1','REAL 2', 'IMAG 2', 
    inrow = sntabi.Desc.Dict['nrow']; 
    for irow in range(1,inrow+1):
        row=sntabi.ReadRow(irow,err)
        # Swap
        for i in range(0,numIF):
            tre1 = row['REAL1'][i]; tim1 = row['IMAG1'][i]
            tre2 = row['REAL2'][i]; tim2 = row['IMAG2'][i]
            twt1 = row['WEIGHT 1'][i]; twt2 = row['WEIGHT 2'][i]
            row['REAL1'][i]=tre2;  row['IMAG1'][i]=tim2
            row['REAL2'][i]=tre1;  row['IMAG2'][i]=tim1
            row['WEIGHT 1'][i]=twt2; row['WEIGHT 2'][i]=twt1;  
        sntabo.WriteRow (irow, row, err)
    sntabi.Close(err); sntabo.Close(err)
    OErr.printErrMsg(err, "Error Swapping SN tables")
    # History entry
    outHistory = History.History("history", uv.List, err)
    z=History.POpen(outHistory, History.READWRITE, err)
    z=History.PTimeStamp(outHistory," Start Obit SwapSNXY",err)
    z=History.PWriteRec(outHistory,-1,"SwapSNXY / inSN = "+str(inVer),err)
    z=History.PWriteRec(outHistory,-1,"SwapSNXY / outSN = "+str(outVer),err)
    z=History.PClose(outHistory, err)
    OErr.printErrMsg(err, "Error with history")

# end SwapSNXY

del PDNeg
def PDNeg(uv, err, inVer=1, outVer=2):
    """
    Copy PD table inVer to outVer negating X and Y

    * uv       = Python Obit UV object, AIPS or FITS
    * err      = Python Obit Error/message stack
    * inVer    = input PD version
    * outVer   = output PD version
    """
    from math import pi
    import FArray, Table, History
    fblank=FArray.fblank
    uv.Header(err)
    # Copy table
    print ("PDNeg inVer=",inVer,"outVer=", outVer)
    pdtabi=uv.NewTable(Table.READONLY,"AIPS PD",inVer,err)
    pdtabi.Open(Table.READONLY,err)
    numIF=pdtabi.keys['numIF']; numPol=pdtabi.keys['numPol'];
    numChan=pdtabi.keys['numChan']; 
    pdtabo=Table.PClone(pdtabi,None)
    pdtabo=uv.NewTable(Table.WRITEONLY,"AIPS PD",outVer,err,
                       numIF=numIF,numPol=numPol,numChan=numChan)
    pdtabo.Open(Table.WRITEONLY,err)
    # Patch things in header which may get lost
    pdtabo.keys['numAnt']=pdtabi.keys['numAnt']
    pdtabo.keys['polType']=pdtabi.keys['polType']
    # 'REAL 1', 'IMAG 1','REAL 2', 'IMAG 2', 'P_DIFF'
    inrow = pdtabi.Desc.Dict['nrow']; nchan= pdtabi.Desc.Dict['repeat'][5];
    for irow in range(1,inrow+1):
        row=pdtabi.ReadRow(irow,err)
        for i in range(0,nchan):
            # Negate
            tre1 = row['REAL 1'][i]; tim1 = row['IMAG 1'][i]
            tre2 = row['REAL 2'][i]; tim2 = row['IMAG 2'][i]
            row['REAL 1'][i]=-tre1;  row['IMAG 1'][i]=-tim1
            row['REAL 2'][i]=-tre2;  row['IMAG 2'][i]=-tim2
            if  row['P_DIFF'][i]!=fblank:
                row['P_DIFF'][i] = -row['P_DIFF'][i]
            pdtabo.WriteRow (irow, row, err)
    pdtabi.Close(err); pdtabo.Close(err)
    OErr.printErrMsg(err, "Error Negating PD tables")
    # History entry
    outHistory = History.History("history", uv.List, err)
    z=History.POpen(outHistory, History.READWRITE, err)
    z=History.PTimeStamp(outHistory," Start Obit PDNeg",err)
    z=History.PWriteRec(outHistory,-1,"PDNeg / inPD = "+str(inVer),err)
    z=History.PWriteRec(outHistory,-1,"PDNeg / outPD = "+str(outVer),err)
    z=History.PClose(outHistory, err)
    OErr.printErrMsg(err, "Error with history")
# end PDNeg

del MKPolCalXY
def MKPolCalXY (uv, cals, err, refAnt=59, nthreads=1, 
                doCalib=True, doClip=False, doPCal=True, 
                doSwap=False, doNeg=False, doXYDly=True,
                ChWid=17, timeAvg=2.0, chAvg=8, solInt=1440, fitType=0,
                debug = False, noScrat=[0,0,0]):
    """
    Polarization calibrate a MeerKAT dataset using PCal and XYDly

    A list of "standard" calibrators has known MeerKAT L band values.
    Log file given in file 'PolnCal_'+uv.Aname.strip()+'.log'
    * uv       = Python Obit UV object, AIPS or FITS
    * cals     = list of calibrator source names, recognizes known calibrators:
                 unpolarized: 1934-638, J1939-6342, 0408-65, 
                 Polarized: J0108+0134, 3C138, 3C286, J1130-1449, J2329-4730
                 Others will be included in PCal as sources with joint
                 source & instrumental polarization.
    * err      = Python Obit Error/message stack
    * refAnt   = the reference antenna used for parallel hand calibration
                 NB: this is NOT used as the reference antenna in PCal.
    * nthreads = the number of threads to be used in calculations
    * doClip   = If True, clip excessively high values in data.
    * doCalib  = If True, do point source calibration for sources in cals writing a new CL table
                 All SN tables and CL versions>1 are deleted
    * doPCal   = If True, run PCal, old PD tables are removed
    * doSwap   = If True swap X and Y solutions if PD table from PCal, create new PD
    * doNeg    = If True negate X and Y solutions if PD table from PCal, create new PD
    * doXYDly  = If True, run XYDly
    * ChWid    = width of block of channels for PCal
    * timeAvg  = Data Averaging time (min) for PCal, XYDly
    * chAvg    = Number of channels to average in XYDly
    * solInt   = Solution interval for XYDly
    * fitType  = fitting type in XYDly, 0=>both, 1=>XY, 2=>YX
    * debug    = if True, save input files
    * noScrat  = List of AIPS disks to avoid for scratch files
    """
    import math
    # Control parameters
    # Known L Band calibrators, EVPA expressed as R-L phase (EVPA*2)
    MKCals = {}
    MKCals['1934-638']   = {"PPol":0.0,  "dPPol":0.0,  "RLPhase":0.0,  "RM":0.0, "Clip":20.}
    MKCals['J1939-6342'] = {"PPol":0.0,  "dPPol":0.0,  "RLPhase":0.0, "RM":0.0,  "Clip":20.}
    MKCals['0408-65']    = {"PPol":0.0,  "dPPol":0.0,  "RLPhase":0.0, "RM":0.0,  "Clip":25.}
    MKCals['3C138']      = {"PPol":0.06, "dPPol":0.03, "RLPhase":-61.,"RM":-3.8, "Clip":12.}
    MKCals['J0521+1638'] = {"PPol":0.06, "dPPol":0.03, "RLPhase":-61.,"RM":-3.8, "Clip":12.}
    MKCals['3C286']      = {"PPol":0.09, "dPPol":0.025,"RLPhase":47., "RM":-1.3, "Clip":20.}
    MKCals['J1331+3030'] = {"PPol":0.09, "dPPol":0.025,"RLPhase":47., "RM":-1.3, "Clip":20.}
    MKCals['J1130-1449'] = {"PPol":0.03 ,"dPPol":0.035,"RLPhase":57., "RM":43.2, "Clip":10.}
    MKCals['J2329-4730'] = {"PPol":0.025,"dPPol":0.03, "RLPhase":-6., "RM":16.0, "Clip":10.}
    MKCals['J0108+0134'] = {"PPol":0.035,"dPPol":0.0,  "RLPhase":89.2, "RM":-13.3,"Clip":10.}
    MKCals['Other']      = {"PPol":0.0,  "dPPol":0.0,  "RLPhase":-999.,"RM":0.0, "Clip":20.}

    # Clip
    if doClip:
        af=ObitTask('AutoFlag'); setname(uv,af)
        af.minAmp=1.0
        af.nThreads=nthreads; af.flagVer=1; af.flagTab=1;
        i=0
        for c in cals:
            if c in MKCals:
                af.Sources[0] = c; af.IClip=[MKCals[c]["Clip"],0.1]; 
                af.XClip=[max(2.0,5*MKCals[c]["Clip"]*MKCals[c]["PPol"]),0.1]; 
            print ("Clip",c, af.IClip,af.XClip)
            af.debug = debug
            af.g
    
    # end Clip loop
    # Calibrate?
    temp_CL = -1; temp_SN = -1           # No Temporary calibration tables yet
    aclver = uv.GetHighVer("AIPS CL")    # default initial CL table
    snver  = uv.GetHighVer("AIPS SN")+1  # default initial old SN table
    if doCalib:
        # Clear any old tables
        z=uv.ZapTable('AIPS CL',2,err);z=uv.ZapTable('AIPS CL',3,err);z=uv.ZapTable('AIPS CL',4,err)
        z=uv.ZapTable('AIPS SN',-1,err)
        snver = uv.GetHighVer("AIPS SN")+1  # SN table to write
        clver = uv.GetHighVer("AIPS CL")    # Initial CL table to apply
        fgver = uv.GetHighVer("AIPS FG")    # FlaG table to apply

        calib=ObitTask.ObitTask('Calib'); setname(uv,calib)
        calib.solnVer=snver
        calib.doCalib=1;calib.gainUse=clver;calib.flagVer=fgver;calib.doPol=False
        calib.modelFlux=1.0; calib.refAnt=[refAnt,56]; calib.prtLv=1
        calib.solInt=0.5; calib.solType='L1'; calib.solMode='P'
        calib.avgPol=True; calib.avgIF=False; calib.nThreads=nthreads; 
        calib.Sources = cals; calib.noScrat = noScrat
        calib.debug = debug
        calib.g
        
        clcal=ObitTask.ObitTask('CLCal'); setname(uv,clcal)
        clcal.interMode='SELF'; clcal.maxInter=1440.0
        clcal.solnVer=snver; clcal.calIn=clver; clcal.calOut=clver+1; clcal.refAnt=refAnt
        clcal.debug = debug
        clcal.g
        aclver = clver+1                   # Version of  CL to apply in PCal, XYDly
    
    # end doCalib
  
    # Run PCal in blocks of ChWid channels
    imhead(uv) 
    outSoln = uv.GetHighVer("AIPS PD")+1    # extant high PD table version
    fgver   = uv.GetHighVer("AIPS FG")      # FlaG table to apply
    pcal=ObitTask.ObitTask('PCal'); setname(uv,pcal)
    pcal.UVRange=[3.0,200.]  # Ignore very short
    pcal.doCalib=2;pcal.gainUse=aclver; pcal.flagVer=fgver; pcal.solnType='  '
    pcal.ChWid=ChWid;pcal.ChInc=ChWid; 
    pcal.CPSoln=outSoln; pcal.BPSoln=outSoln; pcal.PDSoln=outSoln
    pcal.refAnt=-1;pcal.nThreads=nthreads; pcal.prtLv=2; pcal.solInt=timeAvg
    pcal.doFitRL=False; pcal.doFitGain=False; pcal.doFitOri=True
    pcal.doFitI=[True,True,True,True,True,True]
    pcal.noScrat = noScrat
    if (uv.FileType=="AIPS"):
        pcal.taskLog='PolnCal_'+uv.Aname.strip()+'.log'
    elif (uv.FileType=="FITS"):
        pcal.taskLog='PolnCal_'+uv.Fname.strip()+'.log'
    else:
        pcal.taskLog='PolnCal.log'
    i=0
    # Unpolarized calibrators
    srcs = []
    for c in cals:
        if (c in MKCals) and (MKCals[c]["PPol"]<=0.0):
            srcs.append(c)
            pcal.doFitPol[i] = False
            pcal.PPol[i]     = 0.0
            pcal.dPPol[i]    = 0.0
            pcal.RLPhase[i]  = 0.0
            pcal.RM[i]       = 0.0; i+=1
    
    # end loop
    # Calibrator to fit - not in MKCals
    for c in cals:
        if not (c in MKCals):
            srcs.append(c)
            pcal.doFitPol[i]=True; pcal.doFitI[i]=True
            pcal.RLPhase[1]=MKCals["Other"]["RLPhase"]; 
            pcal.RM[i]=MKCals["Other"]["RM"]; 
            pcal.PPol[i]=MKCals["Other"]["PPol"]; 
            pcal.dPPol[i]=MKCals["Other"]["dPPol"]; 
            i+=1
    
    pcal.Sources = srcs
    pcal.debug = debug
    if doPCal:
        # Clear any old PD tables
        z=uv.ZapTable('AIPS PD',-1,err)
        pcal.PDSoln=1; 
 
        pcal.i
        pcal.g
        if doSwap:
            SwapPDXY(uv, err, pcal.PDSoln, pcal.PDSoln+1)
        elif doNeg:
            PDNeg(uv, err, pcal.PDSoln, pcal.PDSoln+1)
    # end if doPCal
    
    # X-Y Phase and delay using XYDly
    imhead(uv)
    outSoln = uv.GetHighVer("AIPS SN")+1  # Solution tables to write
    pdver   = uv.GetHighVer("AIPS PD")    # PD table to apply
    xydly=ObitTask.ObitTask('XYDly'); setname(uv,xydly)
    xydly.UVR_Full=[3.0,200.]; xydly.WtUV=0.1  # Downweight very short
    xydly.doCalib=2;xydly.gainUse=aclver;xydly.flagVer=fgver; 
    xydly.doPol=True; xydly.PDVer=pdver; xydly.keepLin=True
    xydly.SNSoln=outSoln; xydly.minSNR=5.; xydly.fitType=fitType
    xydly.refAnt=refAnt;xydly.nThreads=nthreads; xydly.prtLv=2
    xydly.timeAvg=timeAvg; xydly.chAvg=max(1,chAvg); xydly.solInt=solInt
    xydly.noScrat = noScrat
    xydly.taskLog = pcal.taskLog
    i=0
    # Polarized calibrators
    srcs = []
    for c in cals:
        if (c in MKCals) and (MKCals[c]["PPol"]>0.0):
            srcs.append(c)
            xydly.EVPA[i]  = MKCals[c]["RLPhase"]/2.
            xydly.RM[i]    = MKCals[c]["RM"]; i+=1
    
    # end loop
    xydly.Sources = srcs
    xydly.debug = debug
    if doXYDly:
        xydly.i
        xydly.g
    
    # Apply to CL table
    nclver = uv.GetHighVer("AIPS CL")    # Current high CL version
    clcal=ObitTask.ObitTask('CLCal'); setname(uv,clcal)
    clcal.maxInter=1440.0; clcal.refAnt=-1
    clcal.solnVer=xydly.SNSoln; clcal.calIn=aclver; clcal.calOut=nclver+1; 
    clcal.debug = debug
    if doXYDly:
        clcal.g
        print ("Calibration tables: PD",pdver, ", CL",clcal.calOut)
        print ("Log file for PCal,XYDly in",pcal.taskLog)
    
   
 # end MKPolCalXY

del CorrPD
def CorrPD(uv, err, inPD=1, outPD=2, RLCorr=0.0, RLRM=0.0):
    """
    Copy PD table inVer to outVer applying EVPA/RM corrections

    * uv       = Python Obit UV object, AIPS or FITS
    * err      = Python Obit Error/message stack
    * inPD    = input PD version
    * outPD    = output PD version
    * RLCorr   = EVPA correction at ref. freq (deg)
    * RLRM     = Rotation measure correction (rad/m^2)
    """
    import UV, Table, math, FArray, OErr, UVPolnUtil
    clight = 2.997924562e8   # Speed of light m/s
    print ("CorrPD inPD=",inPD,"outPD=", outPD,"RLCorr=",RLCorr, "RLRM=",RLRM)
    # Open and close to update uv descriptor
    err = OErr.OErr()
    uv.Open(UV.READONLY,err)
    uv.Close(err)
    refFreq = uv.Desc.Dict['crval'][uv.Desc.Dict['jlocf']]
    refLamb2 = (clight/refFreq)**2   # Reference lambda**2 for RL correction m**2
    lamb2 = UVPolnUtil.GetLamb2(uv,err)
    
    inpdtab  = uv.NewTable(Table.READONLY,"AIPS PD",inPD,err)
    nrow  = inpdtab.Desc.Dict['nrow']
    nif   = inpdtab.Desc.List.Dict['NO_IF'][2][0];   nchan = inpdtab.Desc.List.Dict['NO_CHAN'][2][0]
    nant  = inpdtab.Desc.List.Dict['NO_ANT'][2][0];  npol  = inpdtab.Desc.List.Dict['NO_POL'][2][0]; 
    ncorr = nif*nchan
    tPD   = uv.NewTable(Table.READWRITE,"AIPS PD",outPD,err,numIF=nif,numChan=nchan, numPol=npol)
    outpdtab = Table.PClone(inpdtab,tPD)
    
    fblank = FArray.fblank
    inpdtab.Open(Table.READONLY, err)
    outpdtab.Open(Table.READWRITE, err); orow = 1
    RLCor_r = math.radians(RLCorr)
    for irow in range (1,nrow+1):
        row = inpdtab.ReadRow(irow, err)
        for i in range(0,ncorr):
            corr = RLCor_r+(lamb2[i]-refLamb2)*RLRM # Phase for frequency
            if row['IMAG 1'][i]!=fblank:
                row['IMAG 1'][i] += corr
                if row['IMAG 2'][i]!=fblank:
                    row['IMAG 2'][i] += corr
                    
            outpdtab.WriteRow(orow, row,err); orow+=1
        
    # Other header stuff
    outpdtab.keys['numAnt']  = inpdtab.keys['numAnt']
    outpdtab.keys['polType'] = inpdtab.keys['polType']
    Table.PDirty(outpdtab)
    # Seems OK in spite of what this says print ("NO_ANT",outpdtab.Desc.List.Dict['NO_ANT'][2][0])
    outpdtab.Close(err)
    inpdtab.Close(err)
        
    # History entry
    outHistory = History.History("history", uv.List, err)
    z=History.POpen(outHistory, History.READWRITE, err)
    z=History.PTimeStamp(outHistory," Start Obit CorrPD",err)
    z=History.PWriteRec(outHistory,-1,"CorrPD / inPD = "+str(inPD),err)
    z=History.PWriteRec(outHistory,-1,"CorrPD / outPD = "+str(outPD),err)
    z=History.PWriteRec(outHistory,-1,"CorrPD / RLCorr = "+str(RLCorr),err)
    z=History.PWriteRec(outHistory,-1,"CorrPD / RLRM = "+str(RLRM),err)
    z=History.PClose(outHistory, err)
    OErr.printErrMsg(err, "Error with history")

# end CorrPD
