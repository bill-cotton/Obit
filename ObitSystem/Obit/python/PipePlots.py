# Package of diagnostic plots for Obit pipelines
# Uses matplotlib when possible, otherwise ObitPlot or AIPS
#exec(open('PipePlots.py').read())
MKPlotXYBPTab=None; MKPlotBPTab=None; MKPlotPDTab=None;  VisPlot=None;
PlotSNTab=None; PlotSNDlyTab=None; isCirc=None

import UV, Image, OErr, RMFit, OSystem, Table, ObitTask
import ImageDesc, FArray, OPlot
from OTObit import setname, addParam
import math
from math import isnan, degrees, radians, atan2, cos, sin
from UVPolnUtil import GetFreqArr
from PipeUtil import ProjMetadata
from PipeUtil import printMess

del MKPlotXYBPTab
def MKPlotXYBPTab(uv, BPVer, plotfile, err, nplot=4,
                logfile=None, check=False, debug=False):
    """
    * Plot cross-hand phase from a BP table, phase of second poln.
    * tolerates failure
    * uv       = UV data object to plot
    * BPVer    = BP table version number, 0-> highest
    * PlotFile = root of plot file, ".pdf" or ".ps" added
    * err      = Obit error/message stack
    * nplot    = number of plots per page
    * logfile  = logfile for messages
    * check    = Only check script, don't execute tasks
    * debug    = show input
    """
    ################################################################
    from OTObit import day2dhms
    mess = "MKPlotXYBPTab with plotfile "+plotfile+" BPVer "+str(BPVer)
    print (mess)
    if check:
        return 0
    fblank = FArray.fblank
    try:
        bptab = uv.NewTable(Table.READONLY, "AIPS BP",BPVer, err)
        bptab.Open(Table.READONLY, err)
        # This work?
    except Exception as exception:
        print(exception)
        mess = "MKPlotXYBPTab Failed "
        printMess(mess, logfile)
        return 0
    # Get frequencies
    freqs = GetFreqArr(uv, err)
    nchan = len(freqs)
    nrow  = bptab.Desc.Dict['nrow']
    for i in range(0,len(freqs)):
        freqs[i] *= 1.0e-6 # to MHz
    xlabel = "Frequency (MHz)"; ylabel = "X-Y Phase (deg.)"
    # Get antenna names
    meta = ProjMetadata( uv, "AIPSVERSION", err)

    # Plot - first try matplotlib
    nplt = 0
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.backends.backend_pgf import PdfPages
        # Loop over pages 
        with PdfPages(plotfile+'.pdf') as pdf:
            # loop over row
            for ir in range(1,nrow+1):
                r = bptab.ReadRow(ir, err)
                phz = []; x = []
                for j in range(0,nchan):
                    if r['REAL 2'][j] != fblank:
                        p = degrees(atan2(r['IMAG 2'][j],r['REAL 2'][j]))
                        phz .append(p); x.append(freqs[j])
                aant = max(1,r['ANTENNA'][0])
                if ((aant-1)>=len(meta['anNames'])):  #Shouldn't happen
                    continue
                plt.figure()
                plt.title("X-Y Phase Antenna "+str(aant)+" ("+meta['anNames'][aant-1]+")")
                plt.xlabel(xlabel); plt.ylabel(ylabel)
                plt.scatter(x, phz, marker="+")
                pdf.savefig(); plt.close()  # add to file
                nplt += 1  # Count
                OErr.printErrMsg(err, "Error plotting BP Table")
        print ("MKPlotXYBPTab: Wrote",nplt,"pages")
        return
    except Exception as exception:
        print(exception)
        print ("Will try PLPlot")
        
    # PLPlot verson
    try:
        plot = OPlot.newOPlot("plot", err, output=plotfile+"/ps", ny=nplot)
        plot.List.set("XLABEL",xlabel);    
        plot.List.set("YLABEL",ylabel)    
        # loop over row
        for ir in range(1,nrow+1):
            r = bptab.ReadRow(ir, err)
            aant = max(1,r['ANTENNA'][0])
            if ((aant-1)>=len(meta['anNames'])):  #Shouldn't happen
                continue
            plot.List.set("TITLE","X-Y Phase Antenna "+str(aant)+" ("+meta['anNames'][aant-1]+")")
            phz = []; x = []
            for j in range(0,nchan):
                if r['REAL 2'][j] != fblank:
                    p = degrees(atan2(r['IMAG 2'][j],r['REAL 2'][j]))
                    phz.append(p); x.append(freqs[j])
            OPlot.PXYPlot(plot, 2, x, phz, err)
            OPlot.PShow(plot,err)
            OErr.printErrMsg(err, "Error plotting BP Table")
    except Exception as exception:
        print(exception)
        mess = "MKPlotXYBPTab Failed "
        printMess(mess, logfile)
        # Allow failure
        bptab.Close(err)
        return 0
    else:
        pass
    bptab.Close(err)
    OErr.printErrMsg(err, "Error plotting BP Table")
    return 0
# end MKPlotXYBPTab

del MKPlotBPTab
def MKPlotBPTab(uv, BPVer, plotfile, err, 
                logfile=None, check=False, debug=False):
    """
    * Plot a BP table, tolerates failure
    * uv       = UV data object to plot
    * BPVer    = BP table version number, 0-> highest
    * PlotFile = root of plot file, ".pdf" or ".ps" added
    * err      = Obit error/message stack
    * logfile  = logfile for messages
    * check    = Only check script, don't execute tasks
    * debug    = show input
    """
    ################################################################
    from OTObit import day2dhms
    mess = "MKPlotBPTab with plotfile "+plotfile+" BPVer "+str(BPVer)
    print (mess)
    if check:
        return 0
    fblank = FArray.fblank
    try:
        bptab = uv.NewTable(Table.READONLY, "AIPS BP",BPVer, err)
        bptab.Open(Table.READONLY, err)
        # This work?
    except Exception as exception:
        print(exception)
        mess = "MKPlotBPTab Failed "
        printMess(mess, logfile)
        return 0
    if 'REAL 2' in bptab.Desc.Dict['FieldName']:
        npol = 2;
    else:
        npol = 1
    # Get frequencies
    freqs =  GetFreqArr(uv, err)
    nchan = len(freqs)
    nrow  = bptab.Desc.Dict['nrow']
    for i in range(0,len(freqs)):
        freqs[i] *= 1.0e-6  # in MHz
    xlabel = "Frequency (MHz)"; 
    # Get antenna names
    meta = ProjMetadata( uv, "AIPSVERSION", err)
    # Circular or linear feeds?
    if isCirc(uv):
        RX = "R"; LY = "L"
    else:
        RX = "X"; LY = "Y"
    # Plot - first try matplotlib
    nplot = 0
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.backends.backend_pgf import PdfPages
        freq = np.array(freqs)
        # Loop over pages 
        with PdfPages(plotfile+'.pdf') as pdf:
            # loop over row
            for ir in range(1,nrow+1):
                r = bptab.ReadRow(ir, err)
                timeStr = day2dhms(r['TIME'][0])[0:12]
                aant = max(1,r['ANTENNA'][0])
                if ((aant-1)>=len(meta['anNames'])):  #Shouldn't happen
                    continue
                re1 = np.array(r['REAL 1']); im1 = np.array(r['IMAG 1']);
                x  = np.where(re1==fblank,np.nan, freq)
                # Any Valid data?
                if np.sum(x, where=~np.isnan(x))<=0.0:
                    continue
                re = np.where(re1==fblank,np.nan, re1)
                im = np.where(im1==fblank,np.nan, im1)
                amp = (re*re+im*im)**0.5; phs = np.degrees(np.arctan2(im,re))
                fig, ax = plt.subplots(nrows=2*npol,sharex='col')
                z=ax[0].set_title("BP Antenna "+str(aant)+" ("+meta['anNames'][aant-1]+") "\
                                  +timeStr)
                z=ax[0].set_ylabel(RX+" Amp")
                z=ax[0].scatter(x, amp, marker="+")
                if (npol==1):
                    z=ax[1].set_xlabel(xlabel);
                z=ax[1].set_ylabel(RX+" Phase ($^\circ$)")
                z=ax[1].scatter(x, phs, marker="+")
                if npol==2:
                    re2 = np.array(r['REAL 2']); im2 = np.array(r['IMAG 2']);
                    x  = np.where(re2==fblank,np.nan, freq)
                    re = np.where(re2==fblank,np.nan, re2)
                    im = np.where(im2==fblank,np.nan, im2)
                    amp = (re*re+im*im)**0.5; phs = np.degrees(np.arctan2(im,re))
                    z=ax[2].set_ylabel(LY+" Amp")
                    z=ax[2].scatter(x, amp, marker="+")
                    z=ax[3].set_xlabel(xlabel); ax[3].set_ylabel(LY+" Phase ($^\circ$)")
                    z=ax[3].scatter(x, phs, marker="+")
                pdf.savefig(); plt.close(fig=fig)  # add to file
                nplot += 1  # Count
                OErr.printErrMsg(err, "Error plotting BP Table")
        print ("MKPlotBPTab: Wrote",nplot,"pages")
        return
    except Exception as exception:
        print(exception)
        print ("Will try PLPlot")
        
    # PLPlot verson
    try:
        plot = OPlot.newOPlot("plot", err, output=plotfile+"/ps",bgcolor=OPlot.WHEAT, ny=npol*2)
        plot.List.set("XLABEL",xlabel)    
        # loop over row
        for ir in range(1,nrow+1):
            r = bptab.ReadRow(ir, err)
            aant = max(1,r['ANTENNA'][0])
            if ((aant-1)>=len(meta['anNames'])):  #Shouldn't happen
                break
            plot.List.set("TITLE","BP Antenna "+str(r['ANTENNA'][0])+\
                          " ("+meta['anNames'][aant-1]+")")
            amp1 = []; amp2 = []; phs1 = []; phs2 = []; x1 = []; x2 = []
            for j in range(0,nchan):
                if r['REAL 1'][j] != fblank:
                    amp = (r['REAL 1'][j]**2 + r['IMAG 1'][j]**2)**0.5
                    amp1.append(amp); x1.append(freqs[j])
                    ph = degrees(atan2(r['IMAG 1'][j],r['REAL 1'][j]))
                    ph1.append(ph)
                if (npol==2) and (r['REAL 2'][j] != fblank):
                    amp = (r['REAL 2'][j]**2 + r['IMAG 2'][j]**2)**0.5
                    amp2.append(amp); x2.append(freqs[j])
                    ph = degrees(atan2(r['IMAG 2'][j],r['REAL 2'][j]))
                    ph2.append(ph)
            plot.List.set("YLABEL",RX+" Amp")    
            OPlot.PXYPlot(plot, 2, x1, amp1, err)
            plot.List.set("TITLE","")
            plot.List.set("YLABEL",RX+" Phase")    
            OPlot.PXYPlot(plot, 2, x1, ph1, err)
            if npol==2:
                plot.List.set("YLABEL",LY+" Amp")    
                OPlot.PXYPlot(plot, 2, x2, amp2, err)
                plot.List.set("YLABEL",LY+" Phase")    
                OPlot.PXYPlot(plot, 2, x2, ph2, err)
        OPlot.PShow(plot,err)
        OErr.printErrMsg(err, "Error plotting BP Table")
    except Exception as exception:
        print(exception)
        mess = "MKPlotAmpBPTab Failed "
        printMess(mess, logfile)
        # Allow failure
        bptab.Close(err)
        return
    else:
        pass
    bptab.Close(err)
    OErr.printErrMsg(err, "Error plotting BP Table")
    return
# end MKPlotBPTab

del MKPlotPDTab
def MKPlotPDTab(uv, PDVer, plotfile, err, 
                logfile=None, check=False, debug=False):
    """
    * Plot a PD table, tolerates failure
    * uv       = UV data object to plot
    * PDVer    = PD table version number, 0-> highest
    * PlotFile = root of plot file, ".pdf" or ".ps" added
    * err      = Obit error/message stack
    * logfile  = logfile for messages
    * check    = Only check script, don't execute tasks
    * debug    = show input
    """
    ################################################################
    from OTObit import day2dhms
    mess = "MKPlotPDTab with plotfile "+plotfile+" PDVer "+str(PDVer)
    print (mess)
    if check:
        return 0
    fblank = FArray.fblank
    try:
        pdtab = uv.NewTable(Table.READONLY, "AIPS PD",PDVer, err)
        pdtab.Open(Table.READONLY, err)
        # This work?
    except Exception as exception:
        print(exception)
        mess = "MKPlotPDTab Failed "
        printMess(mess, logfile)
        return 0
    npol = 2;
    # Get frequencies
    freqs =  GetFreqArr(uv, err)
    nchan = len(freqs)
    nrow  = pdtab.Desc.Dict['nrow']
    for i in range(0,len(freqs)):
        freqs[i] *= 1.0e-6  # in MHz
    xlabel = "Frequency (MHz)"; 
    # Get antenna names
    meta = ProjMetadata( uv, "AIPSVERSION", err)
    # Circular or linear feeds?
    if isCirc(uv):
        RX = "R"; LY = "L"
    else:
        RX = "X"; LY = "Y"
    # Plot - first try matplotlib
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.backends.backend_pgf import PdfPages
        freq = np.array(freqs)
        # Loop over pages 
        with PdfPages(plotfile+'.pdf') as pdf:
            # loop over row
            for ir in range(1,nrow+1):
                r = pdtab.ReadRow(ir, err)
                aant = max(1,r['ANTENNA'][0])
                if ((aant-1)>=len(meta['anNames'])):  #Shouldn't happen
                    continue
                re1 = np.array(r['REAL 1']); im1 = np.array(r['IMAG 1']);
                x  = np.where(re1==fblank,np.nan, freq)
                # Any Valid data?
                if np.sum(x, where=~np.isnan(x))<=0.0:
                    continue
                re = np.where(re1==fblank,np.nan, re1)
                im = np.degrees(np.where(im1==fblank,np.nan, im1))
                fig, ax = plt.subplots(nrows=2*npol,sharex='col')
                z=ax[0].set_title("PD Antenna "+str(aant)+" ("+meta['anNames'][aant-1]+")")
                z=ax[0].set_ylabel(RX+" Elip")
                z=ax[0].scatter(x, re, marker="+")
                z=ax[1].set_ylabel(RX+" Ori ($^\circ$)")
                z=ax[1].scatter(x, im, marker="+")
                re2 = np.array(r['REAL 2']); im2 = np.array(r['IMAG 2']);
                x  = np.where(re2==fblank,np.nan, freq)
                re = np.where(re2==fblank,np.nan, re2)
                im = np.degrees(np.where(im2==fblank,np.nan, im2))
                z=ax[2].set_ylabel(LY+" Elip")
                z=ax[2].scatter(x, re, marker="+")
                z=ax[3].set_xlabel(xlabel); ax[3].set_ylabel(LY+" Ori ($^\circ$)")
                z=ax[3].scatter(x, im, marker="+")
                pdf.savefig(); plt.close(fig=fig)  # add to file
                OErr.printErrMsg(err, "Error plotting PD Table")
        return
    except Exception as exception:
        print(exception)
        print ("Will try PLPlot")
        
    # PLPlot verson
    try:
        plot = OPlot.newOPlot("plot", err, output=plotfile+"/ps",bgcolor=OPlot.WHEAT, ny=npol*2)
        plot.List.set("XLABEL",xlabel)    
        # loop over row
        for ir in range(1,nrow+1):
            r = pdtab.ReadRow(ir, err)
            aant = max(1,r['ANTENNA'][0])
            if ((aant-1)>=len(meta['anNames'])):  #Shouldn't happen
                break
            plot.List.set("TITLE","PD Antenna "+str(r['ANTENNA'][0])+\
                          " ("+meta['anNames'][aant-1]+")")
            elp1 = []; elp2 = []; ori1 = []; ori2 = []; x1 = []; x2 = []
            for j in range(0,nchan):
                if r['REAL 1'][j] != fblank:
                    elp1.append(r['REAL 1'][j]); ori1.append(degrees(r['IMAG 1'][j]))
                    x1.append(freqs[j])
                if (npol==2) and (r['REAL 2'][j] != fblank):
                    elp1.append(r['REAL 2'][j]); ori1.append(degrees(r['IMAG 2'][j]))
                    x2.append(freqs[j])
            plot.List.set("YLABEL",RX+" Elip")    
            OPlot.PXYPlot(plot, 2, x1, elp1, err)
            plot.List.set("TITLE","")
            plot.List.set("YLABEL",RX+" Ori")    
            OPlot.PXYPlot(plot, 2, x1, ori1, err)
            if npol==2:
                plot.List.set("YLABEL",LY+" Elip")    
                OPlot.PXYPlot(plot, 2, x2, elp2, err)
                plot.List.set("YLABEL",lY+" Ori")    
                OPlot.PXYPlot(plot, 2, x2, ori2, err)
        OPlot.PShow(plot,err)
        OErr.printErrMsg(err, "Error plotting PD Table")
    except Exception as exception:
        print(exception)
        mess = "MKPlotAmpPDTab Failed "
        printMess(mess, logfile)
        # Allow failure return 1
        pdtab.Close(err)
        return
    else:
        pass
    pdtab.Close(err)
    OErr.printErrMsg(err, "Error plotting PD Table")
    return
# end MKPlotPDTab

del VisPlot
def VisPlot(uv, label, plotfile, err, selAnt=1, XPol=False,
            logfile=None, check=False, debug=False):
    """
    * Plot Visibilities to a selected antenna
    * Only matplotlib supported
    * uv       = UV data object to plot
    * label    = string to plot titles
    * PlotFile = root of plot file, ".pdf" added
    * err      = Obit error/message stack
    * selAnt   = plot baselines involving ant selAnt, 1 -rel
    * XPol     = If True plot cross pol vis, else parallel hand
    * logfile  = logfile for messages
    * check    = Only check script, don't execute tasks
    * debug    = show input
    * Returns True if successful, else failed
    """
    ################################################################
    #selAnt=35; # XPol=False # DEBUG
    #uv.Header(err) #DEBUG
    from OTObit import day2dhms
    mess = "VisPlot to plotfile "+plotfile
    print (mess)
    if check:
        return True
    fblank = FArray.fblank
    if uv.GetHighVer("AIPS NX")<1:
        UV.PUtilIndex(uv, err)  # Needs index
    uv.List.set("nVisPIO", 1)  # 1 vis per call
    z      = uv.Open(UV.READONLY, err)
    idd    = uv.Desc.Dict
    nvis   = idd['nvis']
    nrparm = idd['nrparm']
    lrec   = idd["inaxes"][idd['jlocc']] * idd["inaxes"][idd['jlocs']]
    nchan  = idd["inaxes"][idd['jlocf']]
    nif    = idd["inaxes"][idd['jlocif']]; nchif = nchan*nif
    ns     = idd["inaxes"][idd['jlocs']]
    if ns>=2:
        nstok = 2;
    else:
        nstok = 1
    # Get frequencies
    freqs =  GetFreqArr(uv, err)
    for i in range(0,len(freqs)):
        freqs[i] *= 1.0e-6  # in MHz
    xlabel = "Frequency (MHz)"; 
    # Get antenna names
    meta = ProjMetadata( uv, "AIPSVERSION", err)
    # Circular or linear feeds?
    if isCirc(uv):
        if XPol:  # Cross?
            P = "RL"; Q = "LR"
        else:     # Parallel?
            P = "RR"; Q = "LL"
    else:
       if XPol:  # Cross?
           P = "XY"; Q = "YX"
       else:     # Parallel?
           P = "XX"; Q = "YY";
    if nstok==1:
        P = "I"  # Only Stokes I
    # plot if have matplotlib and numpy
    nplot = 0
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.backends.backend_pgf import PdfPages
        # Plot - first try matplotlib
        freq = np.array(freqs)     # Frequencies as numpy array
        # Loop over pages
        with PdfPages(plotfile+'.pdf') as pdf:
            # loop over row
            for iv in range(1,nvis+1):
                uv.Read(err,firstVis=iv)
                buff = np.array(uv.VisBuf) # I/O buffer as numpy array
                rparm = np.frombuffer(buff,count=nrparm,dtype=np.float32) # Random parms
                # Want this one?
                a1 = int(rparm[idd['ilocb']]//256); a2 = int(rparm[idd['ilocb']]-256*a1)
                if (a1!=selAnt) and (a2!=selAnt):
                    continue
                timeStr = day2dhms(rparm[idd['iloct']])[0:12]
                vis= np.frombuffer(buff,offset=nrparm*4,dtype=np.float32)
                if (XPol):  # Cross hand?
                    P_Real=vis[6::lrec]; P_Imag=vis[7::lrec];  P_Wt=vis[8::lrec]
                    Q_Real=vis[9::lrec]; Q_Imag=vis[10::lrec]; Q_Wt=vis[11::lrec]
                else:      # Parallel hand
                    P_Real=vis[0::lrec]; P_Imag=vis[1::lrec]; P_Wt=vis[2::lrec]
                    if (nstok>1):
                        Q_Real=vis[3::lrec]; Q_Imag=vis[4::lrec]; Q_Wt=vis[5::lrec]
                f = np.where(P_Wt<=0.0,np.nan, freq)
                if np.isnan(np.sum(f, where=~np.isnan(f))):  # Anything?
                    continue
                re = np.where(P_Wt<=0.0,np.nan, P_Real)
                im = np.where(P_Wt<=0.0,np.nan, P_Imag)
                amp = (re*re+im*im)**0.5; phs = np.degrees(np.arctan2(im,re))
                fig, ax = plt.subplots(nrows=2*nstok,sharex='col')
                ax[0].scatter(f, amp, marker="+")
                ax[1].scatter(f, phs, marker="+")
                ax[0].set_title(label+" bl "+str(a1)+"-"+str(a2)+\
                                " ("+meta['anNames'][a1-1]+"-"+meta['anNames'][a2-1]+") "+
                                timeStr)
                ax[0].set_ylabel(P+" Amp (Jy)")
                ax[1].set_ylabel(P+" Phase($^\circ$)")
                if nstok==2:
                    f = np.where(Q_Wt<=0.0,np.nan, freq)
                    re = np.where(Q_Wt<=0.0,np.nan, Q_Real)
                    im = np.where(Q_Wt<=0.0,np.nan, Q_Imag)
                    amp = (re*re+im*im)**0.5; phs = np.degrees(np.arctan2(im,re))
                    ax[2].scatter(f, amp, marker="+")
                    ax[3].scatter(f, phs, marker="+")
                    ax[2].set_ylabel(Q+" Amp (Jy)")
                    ax[3].set_ylabel(Q+" Phase($^\circ$)")
                    ax[3].set_xlabel("Freq (MHZ)")
                else:
                    ax[1].set_xlabel("Freq (MHZ)")
                pdf.savefig(); plt.close(fig=fig)  # add to file
                nplot += 1  # Count
                OErr.printErrMsg(err, "Error plotting Visibility")
        print ("VisPlot: Wrote",nplot,"Spectra")
        return True
    except Exception as exception:
        print(exception)
        OErr.printErrMsg(err, "Error plotting visibilities")
        return False
# end VisPlot

del PlotSNTab
def PlotSNTab(uv, SNVer, plotfile, err, 
                logfile=None, check=False, debug=False):
    """
    * Plot Amp & Phase from SN table, tolerates failure
    * Only matplotlib supported
    * uv       = UV data object to plot
    * SNVer    = SN table version number, 0-> highest
    * PlotFile = root of plot file, ".pdf" or ".ps" added
    * err      = Obit error/message stack
    * logfile  = logfile for messages
    * check    = Only check script, don't execute tasks
    * debug    = show input
    * returns True if successful, else False
    """
    ################################################################
    from OTObit import day2dhms
    mess = "MKPlotSNTab with plotfile "+plotfile+" SNVer "+str(SNVer)
    print (mess)
    if check:
        return True
    # try matplotlib
    nplot = 0
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.backends.backend_pgf import PdfPages
        fblank = FArray.fblank
        try:
            sntab = uv.NewTable(Table.READONLY, "AIPS SN",SNVer, err)
            sntab.Open(Table.READONLY, err)
            # This work?
        except Exception as exception:
            print(exception)
            mess = "MKPlotSNTab Failed "
            printMess(mess, logfile)
            return True
        nrow = sntab.Desc.Dict['nrow']
        nant = sntab.Desc.List.Dict['NO_ANT'][2][0]
        npol = sntab.Desc.List.Dict['NO_POL'][2][0]
        nif  = sntab.Desc.List.Dict['NO_IF'][2][0]
        xlabel = "Time (Hr)"; 
        # Get antenna names
        meta = ProjMetadata( uv, "AIPSVERSION", err)
        # Circular or linear feeds?
        if isCirc(uv):
            P = "R"; Q = "L"
        else:
            P = "X"; Q = "Y"
        # Create dictionary of arrays with time, p1re, p1in, q1re,q1im...pnre,pnim,qnre,qnim
        data ={}
        for j in range(0,nant):
            l=[];
            for i in range(0,1+npol*nif*2):
                l.append([])
            data[j]=l
        
        # Suck in data
        for ir in range(1,nrow+1):
            r = sntab.ReadRow(ir,err)
            ant  =r['ANTENNA NO.'][0]-1 # 0 rel
            ll = data[ant]
            ll[0].append(24.0*r['TIME'][0])  # Hours
            for iif in range(0,nif):
                off = iif*2*npol
                ll[1+off].append(r['REAL1'][iif])
                ll[2+off].append(r['IMAG1'][iif])
                if (npol>1):
                    ll[3+off].append(r['REAL2'][iif])
                    ll[4+off].append(r['IMAG2'][iif])
            # end IF loop
        # end loop over data
        sntab.Close(err)
        OErr.printErrMsg(err, "Error reading SN Table")

        # Loop over pages 
        with PdfPages(plotfile+'.pdf') as pdf:
            # loop over antenna
            for ant in range(0,nant):
                ll = data[ant]  # Fetch data
                # loop over IF
                for iif in range(0,nif):
                    off = iif*2*npol
                    t   = np.array(ll[0])
                    rep = np.array(ll[1+off]); imp =  np.array(ll[2+off]);
                    if(npol==2):
                        req = np.array(ll[3+off]); imq =  np.array(ll[4+off]);
                    
                    fig, ax = plt.subplots(nrows=2*npol,sharex='col')
                    x   = np.where(rep==fblank,np.nan, t)
                    # Any data?
                    if np.sum(x, where=~np.isnan(x))<=0.0:
                        continue
                    re  = np.where(rep==fblank,np.nan, rep)
                    im  = np.where(rep==fblank,np.nan, imp)
                    amp = (re*re+im*im)**0.5; phs = np.degrees(np.arctan2(im,re))
                    z = ax[0].plot(x, amp, "b+")
                    z = ax[1].plot(x, phs, "b+")
                    z = ax[0].set_title("SN Antenna "+str(ant+1)+" ("+meta['anNames'][ant]+") "\
                                  +"IF "+str(iif+1))
                    z = ax[0].set_ylabel(P+" Amp")
                    if (not np.isnan(sum(amp))) and (len(amp)>1):
                        z = ax[0].set_ylim([0.95*min(amp),1.05*max(amp)])
                    z = ax[1].set_ylabel(P+" Phase($^\circ$)")
                    if npol==2:
                        x   = np.where(req==fblank,np.nan, t)
                        re  = np.where(req==fblank,np.nan, req)
                        im  = np.where(req==fblank,np.nan, imq)
                        amp = (re*re+im*im)**0.5; phs = np.degrees(np.arctan2(im,re))
                        z = ax[2].plot(x, amp, "b+")
                        z = ax[3].plot(x, phs, "b+")
                        z = ax[2].set_ylabel(Q+" Amp")
                        if (not np.isnan(sum(amp))) and (len(amp)>1):
                            z = ax[2].set_ylim([0.95*min(amp),1.05*max(amp)])
                        z = ax[3].set_ylabel(Q+" Phase($^\circ$)")
                        z = ax[3].set_xlabel(xlabel)
                    else:
                        z = ax[1].set_xlabel(xlabel)
                    pdf.savefig(); plt.close(fig="all")  # add to file
                    nplot += 1  # Count
        print ("MKPlotSNTab: Wrote",nplot,"pages")
        return True
    except Exception as exception:
        print(exception)
        return False
# end PlotSNTab

del PlotSNDlyTab
def PlotSNDlyTab(uv, SNVer, plotfile, err, 
                logfile=None, check=False, debug=False):
    """
    * Plot Delay from SN table, tolerates failure
    * Only matplotlib supported
    * uv       = UV data object to plot
    * SNVer    = SN table version number, 0-> highest
    * PlotFile = root of plot file, ".pdf" or ".ps" added
    * err      = Obit error/message stack
    * logfile  = logfile for messages
    * check    = Only check script, don't execute tasks
    * debug    = show input
    * returns True if successful, else False
    """
    ################################################################
    from OTObit import day2dhms
    mess = "MKPlotSNDlyTab with plotfile "+plotfile+" SNVer "+str(SNVer)
    print (mess)
    if check:
        return True
    # try matplotlib
    nplot = 0
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.backends.backend_pgf import PdfPages
        fblank = FArray.fblank
        try:
            sntab = uv.NewTable(Table.READONLY, "AIPS SN",SNVer, err)
            sntab.Open(Table.READONLY, err)
            # This work?
        except Exception as exception:
            print(exception)
            mess = "MKPlotSNDlyTab Failed "
            printMess(mess, logfile)
            return True
        nrow = sntab.Desc.Dict['nrow']
        nant = sntab.Desc.List.Dict['NO_ANT'][2][0]
        npol = sntab.Desc.List.Dict['NO_POL'][2][0]
        nif  = sntab.Desc.List.Dict['NO_IF'][2][0]
        xlabel = "Time (Hr)"; 
        # Get antenna names
        meta = ProjMetadata( uv, "AIPSVERSION", err)
        # Circular or linear feeds?
        if isCirc(uv):
            P = "R"; Q = "L"
        else:
            P = "X"; Q = "Y"
        # Create dictionary of arrays with time, p1re, p1in, q1re,q1im...pnre,pnim,qnre,qnim
        data ={}
        for j in range(0,nant):
            l=[];
            for i in range(0,1+npol*nif):
                l.append([])
            data[j]=l
        
        # Suck in data
        for ir in range(1,nrow+1):
            r = sntab.ReadRow(ir,err)
            ant  =r['ANTENNA NO.'][0]-1 # 0 rel
            ll = data[ant]
            ll[0].append(24.0*r['TIME'][0])  # Hours
            for iif in range(0,nif):
                off = iif*npol
                ll[1+off].append(r['DELAY 1'][iif])
                if (npol>1):
                    ll[2+off].append(r['DELAY 2'][iif])
            # end IF loop
        # end loop over data
        sntab.Close(err)
        OErr.printErrMsg(err, "Error reading SN Table")

        # Loop over pages 
        with PdfPages(plotfile+'.pdf') as pdf:
            # loop over antenna
            for ant in range(0,nant):
                ll = data[ant]  # Fetch data
                # loop over IF
                for iif in range(0,nif):
                    off = iif*npol
                    t   = np.array(ll[0])
                    dlyp = np.array(ll[1+off]); 
                    if(npol==2):
                        dlyq = np.array(ll[2+off]);
                    
                    fig, ax = plt.subplots(nrows=npol,sharex='col')
                    x   = np.where(dlyp==fblank,np.nan, t)
                    # Any data?
                    if np.sum(x, where=~np.isnan(x))<=0.0:
                        continue
                    dly  = 1.0e9*np.where(dlyp==fblank,np.nan, dlyp) # to nsec
                    z = ax[0].plot(x, dly, "b+")
                    z = ax[0].set_title("SN Antenna "+str(ant+1)+" ("+meta['anNames'][ant]+") "\
                                  +"IF "+str(iif+1))
                    z = ax[0].set_ylabel(P+" Delay (nsec)")
                    if npol==2:
                        x    = np.where(dlyq==fblank,np.nan, t)
                        dly  = 1.0e-9* np.where(dlyq==fblank,np.nan, dlyq) # to nsec
                        z = ax[1].plot(x, dly, "b+")
                        z = ax[1].set_ylabel(Q+"  Delay (nsec)")
                        z = ax[1].set_xlabel(xlabel)
                    else:
                        z = ax[0].set_xlabel(xlabel)
                    pdf.savefig(); plt.close(fig="all")  # add to file
                    nplot += 1  # Count
        print ("MKPlotSNDlyTab: Wrote",nplot,"pages")
        return True
    except Exception as exception:
        print(exception)
        return False
# end PlotSNDlyTab

del isCirc
def isCirc(uv):
    """
    * Does the file uv use circular feeds?lot a PD table, tolerates failure
    * uv       = UV data object to plot
    * Returns True or False
    """
    ################################################################
    return uv.Desc.Dict['crval'][uv.Desc.Dict['jlocs']]>-4
# end isCirc
