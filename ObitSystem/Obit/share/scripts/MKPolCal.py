# do instrumental polarization cal for MeerKAT L band Datasets
# write uvtab file with BP and PD tables

#exec(open('MKPolCal.py').read())
#exec(open('../../PolCal/MKPolCal.py').read())
MKPolCal = None; del MKPolCal
def MKPolCal (uv, cals, err, refAnt=59, nthreads=1, doCalib=True, doEVPA=False, doClip=False):
    """
    Polarization calibrate a MeerKAT dataset

    A list of "standard" calibrators has known MeerKAT L band values, others will be solved for
    but need observations over a range of parallactic angle.
    Log file from PCal given in file 'PCal_'+uv.Aname.strip()+'.log'
    * uv       = Python Obit UV object, AIPS or FITS
    * cals     = list of calibrator source names, recognizes known calibrators
                 1934-638, J1939-6342, 0408-65, 3C138, J0521+1638, 3C286, J1130-4049, 
                 J2329-4730, J2253+1608 (3C454.3)
                 any other source will be treated as an unknown and be solved for
    * err      = Python Obit Error/message stack
    * refAnt   = the reference antenna used for parallel hand calibration
                 nb: this is NOT used as the reference antenna in PCal.
    * nthreads = the number of threads to be used in model calculation for Calib
    * doCalib  = If True, do point source calibration for sources in cals writing a new CL table
    * doEVPA   = If True, correct EVPA for ad hoc correction,  If true, use PD table 2, else 1.
                 DON'T USE!!!!
    * doClip   = If True,  clip excessively high values in data.
    """
    import math
    # Control parameters
    inTab = 1; outTab = 2      # PD tables for EVPA correction
    RLCorr = math.radians(7.7) # R-L phase correcton in degrees GP82
    RLRM = 1.0                 # R-L phase correction RM rad/m**2
    # Known calibrators
    MKCals = {}
    MKCals['1934-638'] = {"doFitPol":False,"PPol":0.0,"dPPol":0.0,"RLPhase":0.0,"RM":0.0,"Clip":20.}
    MKCals['J1939-6342'] = {"doFitPol":False,"PPol":0.0,"dPPol":0.0,"RLPhase":0.0,"RM":0.0,"Clip":20.}
    MKCals['0408-65']  = {"doFitPol":False,"PPol":0.0,"dPPol":0.0,"RLPhase":0.0,"RM":0.0,"Clip":25.}
    MKCals['3C138']    = {"doFitPol":False,"PPol":0.05,"dPPol":0.03,"RLPhase":-36.,"RM":-1.3,"Clip":12.}
    MKCals['J0521+1638']= {"doFitPol":False,"PPol":0.05,"dPPol":0.03,"RLPhase":-36.,"RM":-1.3,"Clip":12.}
    MKCals['3C286']    = {"doFitPol":False,"PPol":0.09,"dPPol":0.025,"RLPhase":66.,"RM":0.0,"Clip":20.}
    MKCals['J1331+3030'] = {"doFitPol":False,"PPol":0.09,"dPPol":0.025,"RLPhase":66.,"RM":0.0,"Clip":20.}
    MKCals['J1130-4049'] = {"doFitPol":False,"PPol":0.025,"dPPol":0.035,"RLPhase":78.,"RM":34.3,"Clip":10.}
    MKCals['J2329-4730'] = {"doFitPol":False,"PPol":0.025,"dPPol":0.03,"RLPhase":-6.,"RM":16.0,"Clip":10.}
    #MKCals['3C454.3'] = {"doFitPol":False,"PPol":0.06,"dPPol":0.02,"RLPhase":2*28.1,"RM":-58.,"Clip":10.}
    #MKCals['J2253+1608'] = {"doFitPol":False,"PPol":0.06,"dPPol":0.02,"RLPhase":2*28.1,"RM":-58.,"Clip":10.}
    MKCals['3C454.3'] = {"doFitPol":True,"PPol":0.06,"dPPol":0.02,"RLPhase":-999.,"RM":-58.,"Clip":10.}
    MKCals['J2253+1608'] = {"doFitPol":True,"PPol":0.06,"dPPol":0.02,"RLPhase":-999.,"RM":-58.,"Clip":10.}
    MKCals['Other']    = {"doFitPol":True,"PPol":0.0,"dPPol":0.0,"RLPhase":-999.,"RM":0.0,"Clip":20.}

    # Clip
    if doClip:
        af=ObitTask('AutoFlag'); setname(uv,af)
        af.minAmp=1.0
        af.nThreads=nthreads; af.flagVer=1; af.flagTab=1;
        i=0
        for c in cals:
            if c in MKCals:
                print ("Known calibrator",c)
                af.Sources[0] = c; af.IClip=[MKCals[c]["Clip"],0.1]; 
                af.XClip=[max(2.0,5*MKCals[c]["Clip"]*MKCals[c]["PPol"]),0.1]; 
            else:
                print ("Unknown calibrator",c)
                af.Sources[0] = c; af.IClip=[MKCals['Other']["Clip"],0.1]; 
                af.XClip=[max(2.0,5*MKCals['Other']["Clip"]*MKCals['Other']["PPol"]),0.1]; 
                print ("Clip",c, af.IClip,af.XClip)
            af.g
    
    # end loop
    # Calibrate?
    if doCalib:
        z=uv.ZapTable('AIPS CL',2,err);z=uv.ZapTable('AIPS CL',3,err);z=uv.ZapTable('AIPS CL',4,err)
        z=uv.ZapTable('AIPS SN',-1,err)
        snver = uv.GetHighVer("AIPS SN")+1  # SN table to write
        clver = uv.GetHighVer("AIPS CL")    # CL table to apply

        calib=ObitTask('Calib'); setname(uv,calib)
        calib.solnVer=snver
        calib.doCalib=1;calib.gainUse=clver;calib.flagVer=1;calib.doPol=False
        calib.modelFlux=1.0; calib.refAnt=[refAnt,56]; calib.prtLv=1
        calib.solInt=0.5; calib.solType='L1'; calib.solMode='P'
        calib.avgPol=True; calib.avgIF=False; calib.nThreads=nthreads; 
        calib.Sources = cals
        calib.g
        
        clcal=ObitTask('CLCal'); setname(uv,clcal)
        clcal.interMode='SELF'; clcal.maxInter=1440.0
        clcal.solnVer=snver; clcal.calIn=clver; clcal.calOut=clver+1; clcal.refAnt=refAnt
        clcal.g
    
    # end doCalib
    # Run Polcal in blocks of 17 channels
    pcal=ObitTask('PCal'); setname(uv,pcal)
    pcal.UVRange=[3.0,200.]  # Ignore very short
    pcal.doCalib=2;pcal.gainUse=2;pcal.flagVer=1; pcal.solnType='  '
    pcal.ChWid=17;pcal.ChInc=17; pcal.CPSoln=1; pcal.BPSoln=1; pcal.PDSoln=1
    pcal.refAnt=-1;pcal.nThreads=nthreads; pcal.prtLv=2
    pcal.doFitRL=True; pcal.doFitGain=False; pcal.doFitOri=True
    pcal.doFitI=[True,True,True,True,True,True]
    pcal.Sources = cals; 
    if "Aname" in uv.__dict__:
        pcal.taskLog='PCal_'+uv.Aname.strip()+'.log'
    else:
        pcal.taskLog='PCal_'+uv.Fname.strip()+'.log' # Better be FITS
    i=0
    for c in cals:
        if c in MKCals:
            print ("Known calibrator",c)
            pcal.doFitPol[i] = MKCals[c]["doFitPol"]
            pcal.PPol[i]     = MKCals[c]["PPol"]
            pcal.dPPol[i]    = MKCals[c]["dPPol"]
            pcal.RLPhase[i]  = MKCals[c]["RLPhase"]
            pcal.RM[i]       = MKCals[c]["RM"]; i+=1
        else:
            print ("Unknown calibrator",c)
            pcal.doFitPol[i] = MKCals['Other']["doFitPol"]
            pcal.PPol[i]     = MKCals['Other']["PPol"]
            pcal.dPPol[i]    = MKCals['Other']["dPPol"]
            pcal.RLPhase[i]  = MKCals['Other']["RLPhase"]
            pcal.RM[i]       = MKCals['Other']["RM"]; i+=1
    
    # end loop

    pcal.g
    
    # Add correction to MK PCal PD table to correct EVPA based on 3C286
    if (doEVPA):
        import UV, Table, math, FArray, OErr, UVPolnUtil
        # Open and close to update uv descriptor
        err = OErr.OErr()
        uv.Open(UV.READONLY,err)
        uv.Close(err)
        refLamb2 = 0.05699140      # Reference lambda**2 for RL correction m**2
        lamb2 = UVPolnUtil.GetLamb2(uv,err)
        
        inPD  = uv.NewTable(Table.READONLY,"AIPS PD",inTab,err)
        nrow = inPD.Desc.Dict['nrow']
        nif = inPD.Desc.List.Dict['NO_IF'][2][0];   nchan = inPD.Desc.List.Dict['NO_CHAN'][2][0]
        nant = inPD.Desc.List.Dict['NO_ANT'][2][0];  npol = inPD.Desc.List.Dict['NO_POL'][2][0]; 
        ncorr = nif*nchan
        tPD   = uv.NewTable(Table.READWRITE,"AIPS PD",outTab,err,numIF=nif,numChan=nchan, numPol=npol)
        outPD = Table.PClone(inPD,tPD)
        
        fblank = FArray.fblank
        inPD.Open(Table.READONLY, err)
        outPD.Open(Table.READWRITE, err); orow = 1
        for irow in range (1,nrow+1):
            row = inPD.ReadRow(irow, err)
            for i in range(0,ncorr):
                corr = RLCorr+(lamb2[i]-refLamb2)*RLRM # Phase for frequency
                if row['IMAG 1'][i]!=fblank:
                    row['IMAG 1'][i] += corr
                    if row['IMAG 2'][i]!=fblank:
                        row['IMAG 2'][i] += corr
                        
            outPD.WriteRow(orow, row,err); orow+=1
        
        # Other header stuff
        outPD.keys['numAnt']  = inPD.keys['numAnt']
        outPD.keys['polType'] = inPD.keys['polType']
        Table.PDirty(outPD)
        print (outPD.Desc.List.Dict['NO_ANT'])
        outPD.Close(err)
        inPD.Close(err)
        
    
    # end doEVPA
    return pcal   # Debug
# end MKPolCal
