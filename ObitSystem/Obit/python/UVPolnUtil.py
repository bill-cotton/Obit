"""
  Python utility package for fitting & plotting source RM and 
  Instrumental poln from a CP or PD table
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2013
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

import Obit, UV, Table, RMFit, OPlot, OErr
from FArray import fblank
import math
clight = 2.997924562e8   # Speed of light m/s

def FitSource(inUV, souID, err, CPVer=1, inc=1):
    """
    Fit the RM of a source in a CP table
    
    Return a dict with:
        "RM"       = fitted RM (rad/m**2)
        "RMErr"    = fitted RM error
        "EVPA0"    = EVPA at "RefLamb2 (rad)
        "EVPAErr"  = error in EVPA
        "RefLamb2" = reference Lambda^2
        "Lamb2"    = array of lamb2 (m**2)
        "FPol"  = array of fractional poln
        "EVPA   = array of EVPA (rad) of poln
    * inUV     = Obit UV data with tables
    * souID    = source ID
    * err      = Python Obit Error/message stack
    * CPVer    = CP table version
    * inc      = incremant between channels to use
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        print "Actually ",inUV.__class__
        raise TypeError,"inUV MUST be a Python Obit UV"

    # Get table
    cptab = inUV.NewTable(Table.READONLY,"AIPS CP",CPVer,err)
    cptab.Open(Table.READONLY,err)
    nrow = cptab.Desc.Dict["nrow"]
    # Find souID
    IPol = None
    for irow in range (1,nrow+1):
        cprow=cptab.ReadRow(irow,err)
        if souID == cprow["SOURCE ID"][0]:
            IPol  = cprow["I"]
            QPol  = cprow["Q"]
            UPol  = cprow["U"]
            break
    # end loop looking for souID
    # Find it?
    if IPol==None:
        raise  RuntimeError,"No CP entry for source ID "+str(souID)
    # Get lambda^2 array
    lamb2 = GetLamb2(inUV, err)
    # Check consistency
    if len(IPol)!=len(lamb2):
        raise  RuntimeError,"Size of CP entries ("+str(len(IPol))+") NOT number of channels ("+str(len(lamb2))+")"
    # Get fitting parameters with descimation by inc
    l2 = [];  QObs = []; Qsig = []; UObs = []; Usig = [];
    FPol = []; EVPA = []
    for i in range (0, len(lamb2), inc):
        if IPol[i]>0.0:
            l2.append(lamb2[i])
            QObs.append(QPol[i]); Qsig.append(0.001)
            UObs.append(UPol[i]); Usig.append(0.001)
            EVPA.append(0.5*math.atan2(UPol[i], QPol[i]))
            FPol.append(((QPol[i]**2+UPol[i]**2)**0.5)/IPol[i])
    # end loop getting data arrays
    # Fit RM
    refLamb2 = l2[len(l2)/2]
    RMParms = RMFit.PSingle(refLamb2, l2, QObs, Qsig, UObs, Usig, err)
    print "debug RMParms",RMParms
    print "Fitted RM=%f (%f), EVPA=%f(%f)"\
        %(RMParms[0],RMParms[2],math.degrees(RMParms[1]),math.degrees(RMParms[3]))
    return{"RM":RMParms[0], "RMErr":RMParms[2], "EVPA0":RMParms[1], "EVPAErr":RMParms[3], \
               "RefLamb2":refLamb2 ,"Lamb2":l2, "FPol":FPol, "EVPA":EVPA}
    # end FitSource

def GetLamb2(inUV, err, inc = 1):
    """
    Get Array of lambda^2 in data
    
    Return list of channel lambda^2 (m^2)
    * inUV     = Obit UV data with tables
    * err      = Python Obit Error/message stack
    * inc      = increment in channel
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        print "Actually ",inUV.__class__
        raise TypeError,"inUV MUST be a Python Obit UV"
    d = inUV.Desc.Dict # UV data descriptor dictionary
    refFreq = d["crval"][ d["jlocf"]]    # reference frequency
    nchan   = d["inaxes"][ d["jlocf"]]   # Number of channels per IF
    nif     = d["inaxes"][ d["jlocif"]]  # Number of IFs
    chwid   = abs(d["cdelt"][ d["jlocf"]]) # channel width from header
    fqtab = inUV.NewTable(Table.READONLY,"AIPS FQ",1,err)
    fqtab.Open(Table.READONLY,err)
    fqrow = fqtab.ReadRow(1, err)   # Only bother with first
    OErr.printErrMsg(err) # catch table errors
    IFfreq   = fqrow['IF FREQ']     # IF offset
    ChWidth  = fqrow['CH WIDTH']    # IF channel width
    sideband = fqrow['SIDEBAND']    # IF sizeband, 1=USB, -1=LSB
    lamb2 = []
    for iif in range(0,nif):
        for ifq in range(0,nchan,inc):
            frq =  refFreq + IFfreq[iif] + sideband[iif] * ifq * min (chwid, ChWidth[iif])
            lamb2.append((clight/frq)**2)
    # End loops
    return lamb2
    # end GetLamb2

def GetFreqArr(inUV, err, inc = 1):
    """
    Get Array of channel frequencies in data
    
    Return list of channel Frequencies (Hz)
    * inUV     = Obit UV data with tables
    * err      = Python Obit Error/message stack
    * inc      = increment in channel
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        print "Actually ",inUV.__class__
        raise TypeError,"inUV MUST be a Python Obit UV"
    d = inUV.Desc.Dict # UV data descriptor dictionary
    refFreq = d["crval"][ d["jlocf"]]    # reference frequency
    nchan   = d["inaxes"][ d["jlocf"]]   # Number of channels per IF
    nif     = d["inaxes"][ d["jlocif"]]  # Number of IFs
    chwid   = abs(d["cdelt"][ d["jlocf"]]) # channel width from header
    fqtab = inUV.NewTable(Table.READONLY,"AIPS FQ",1,err)
    fqtab.Open(Table.READONLY,err)
    fqrow = fqtab.ReadRow(1, err)   # Only bother with first
    OErr.printErrMsg(err) # catch table errors
    IFfreq   = fqrow['IF FREQ']     # IF offset
    ChWidth  = fqrow['CH WIDTH']    # IF channel width
    sideband = fqrow['SIDEBAND']    # IF sizeband, 1=USB, -1=LSB
    freqs = []
    for iif in range(0,nif):
        for ifq in range(0,nchan,inc):
            frq =  refFreq + IFfreq[iif] + sideband[iif] * ifq * min (chwid, ChWidth[iif])
            freqs.append(frq)
    # End loops
    return freqs
    # end GetFreqArr

def PlotRM(inUV, fitDict, plotFile, title, err):
    """
    Plot data and model
    
    Return list of channel lambda^2 (m^2)
    * inUV     = Obit UV data with tables
    * fitDict  = dict returned by PFitSource
    * plotFile = output ps plotfile name
    * title    = plot title
    * err      = Python Obit Error/message stack
    """
    ################################################################
    plot=OPlot.newOPlot("Plot", err,output=plotFile+"/ps",ny=2)
    # Extract data, fits
    RM    = fitDict["RM"];     RMErr    = fitDict["RMErr"];
    EVPA0 = fitDict["EVPA0"];  EVPAErr  = fitDict["EVPAErr"];
    lamb2 = fitDict["Lamb2"];  refLamb2 = fitDict["RefLamb2"];
    EVPA  = fitDict["EVPA"];   fpol     = fitDict["FPol"];
    # angles to deg
    EVPA0 = math.degrees(EVPA0); EVPAErr = math.degrees(EVPAErr); 
    phs = []
    for angle in EVPA :
        phs.append(math.degrees(angle))
    
    # Plot EVPA
    plot.List.set("TITLE","%s RM=%6.2f (%5.2f) EVPA=%6.2f (%5.2f) "\
                      %(title, RM, RMErr, EVPA0, EVPAErr))
    plot.List.set("YLABEL","EVPA (deg)")
    plot.List.set("CSIZE",1)
    plot.List.set("SSIZE",3)
    plot.List.set("LWIDTH",3)
    OPlot.PXYPlot(plot, 2, lamb2, phs, err)
    # Plot RM fit
    nspec = len(lamb2)
    x1 = lamb2[0] - refLamb2
    y1 = EVPA0 + math.degrees(RM*x1)
    x2 = lamb2[nspec-1] - refLamb2
    y2 = EVPA0 + math.degrees(RM*x2)
    OPlot.PDrawLine(plot, lamb2[0], y1, lamb2[nspec-1], y2, err)
    # Plot fractional poln
    plot.List.set("TITLE"," ")
    plot.List.set("XLABEL","#gl#u2#d (m#u2#d)")
    plot.List.set("YLABEL","frac. pol.")
    plot.List.set("XMAX",lamb2[0])
    plot.List.set("XMIN",lamb2[nspec-1])
    plot.List.set("YMIN",0.0)
    plot.List.set("YMAX",0.1)
    OPlot.PXYPlot(plot, 2, lamb2, fpol, err)
    OPlot.PShow(plot,err)
    # end PlotRM

# Extract/Plot values from PD table
def ExtractPD (inUV, ant, err, PDVer=1, inc=1):
    """ Extract data from PD table

    return dict:
     'Elip1' array of elipticities poln 1 in deg
     'Ori1'  array of orientations poln 1 in deg
     'Elip2' array of elipticities poln 2 in deg
     'Ori2'  array of orientations poln 2 in deg
     'Freqs' array of frequencies corresponding to data
     'refant' reference antenna
     'ant'    antenna
    * inUV  = Obit UV data with tables
    * ant   = antenna number
    * err   = Obit error/message stack
    * PDVer = PD table version
    * inc   = increment in spectrum
    """
      ################################################################
    # Checks
    if not UV.PIsA(inUV):
        print "Actually ",inUV.__class__
        raise TypeError,"inUV MUST be a Python Obit UV"
    
    e1 = []; e2 = []; o1 = []; o2 = []
    # Get table
    PDTab = inUV.NewTable(Table.READONLY,"AIPS PD",PDVer,err)
    # Open table, find antenna
    PDTab.Open(Table.READONLY, err)
    nrow = PDTab.Desc.Dict['nrow']
    for irow in range (1,nrow+1):
        row = PDTab.ReadRow(irow, err)
        if row['ANTENNA'][0]==ant:
            break
    for i in range(0,len(row['REAL 1']),inc):
        e = row['REAL 1'][i]
        # Bad soln?
        if e==fblank:
            e1.append(fblank); e2.append(fblank); 
            o1.append(fblank); o2.append(fblank); 
        else:
            # Fold elipticities to range [-45, 45]
            e = math.degrees(e)
            o = math.degrees(row['IMAG 1'][i])
            if e>45.0:
                e = 90.0 - e
                o = -o
            if o>180.0:
                o -= 360.
            if o<-180.0:
                o += 360
            # Default values probably bad
            if e>=44.999 and o==0.0:
                e = fblank; o = fblank
            e1.append(e)
            o1.append(o)
            e = math.degrees(row['REAL 2'][i])
            o = math.degrees(row['IMAG 2'][i])
            if e<-45.0:
                e = -(90.0 + e)
                o = -o
            if o>180.0:
                o -= 360.
            if o<-180.0:
                o += 360
            if e<=-44.999 and o==0.0:
                e = fblank; o = fblank
            e2.append(e)
            o2.append(o)
    refant = row['REFANT'][0]
    ant    = row['ANTENNA'][0]
    del row
    PDTab.Close(err)
    # Get frequencies
    f = GetFreqArr(inUV, err)
    # select by inc in GHz
    frqs = []
    for i in range(0,len(f),inc):
        frqs.append(f[i]*1.0e-9)
    return {'Elip1':e1, 'Ori1':o1, 'Elip2':e2, 'Ori2':o2, 'refant':refant, 'ant':ant, 'Freqs':frqs}
# end ExtractPD 

def PlotPD (inUV, PDDict, plotfile, err, title="plot"):
    """
    Plot array of PD elipticities and orientations

    * inUV     = Obit UV data with tables
    * fitDict  = dict returned by ExtractPD with entries
     'Elip1' array of elipticities poln 1 in deg
     'Ori1'  array of orientations poln 1 in deg
     'Elip2' array of elipticities poln 2 in deg
     'Ori2'  array of orientations poln 2 in deg
     'Freqs' array of frequencies corresponding to data
     'refant' reference antenna
     'ant'    antenna
    * plotfile = name of postscript plot file
    * err      = Obit Error/message stack
    * title    = title for plot
    """
    Earr = []; Oarr = []
    Earr.append(PDDict['Elip1']); Oarr.append(PDDict['Ori1']);
    Earr.append(PDDict['Elip2']); Oarr.append(PDDict['Ori2']);
    xa = PDDict['Freqs']
    plot = OPlot.newOPlot("plot", err, output=plotfile+"/ps", nx=2, ny=2)
    nplot = len(Earr[0])
    if xa==None:
        x = []
        for i in range(0, len(Earr[0])):
            x.append(float(i))
    else:
        x = xa
    # Loop over poln
    pol = [" RCP"," LCP"]
    for ipol in range (0,2):
        symb = 2
        # Plot elipticities
        plot.List.set("TITLE",title+pol[ipol])
        if (max(Earr[ipol])>0. and max(Earr[ipol])<1000.) or  min(Earr[ipol])>0.:
            plot.List.set("YMIN",35.0); plot.List.set("YMAX",46.0);
            #print "RCP data", min(Earr[ipol])
        else:
            plot.List.set("YMIN",-46.0); plot.List.set("YMAX",-35.0);
            #print "LCP data", min(Earr[ipol])
        if xa==None:
            plot.List.set("XLABEL","Channel")
        else:
            plot.List.set("XLABEL","Frequency (GHz)")
        plot.List.set("YLABEL","Ellipticity (deg)")
        OPlot.PXYPlot(plot, symb, x, Earr[ipol], err)
        # Plot orientations
        symb = 2
        plot.List.set("TITLE","    ")
        plot.List.set("YMIN",-180.0); plot.List.set("YMAX",180.0);
        plot.List.set("YLABEL","Orientation (deg)")
        OPlot.PXYPlot(plot, symb, x, Oarr[ipol], err)
    # end poln loop
    OPlot.PShow(plot,err)
# end PlotPD

