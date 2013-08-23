""" GBT VEGAS Utility module
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

# Python ObitOTF utilities
import Obit, OTF,OErr, FInterpolate, FArray, Image, ImageDesc, InfoList
import Table

def PAverage(inOTF, outOTF, chAvg, err):
    """ Frequency average an OTF

    inOTF   = input Python Obit OTF
    outOTF  = output Python Obit OTF, must be previously defined
    chAvg   = Number of channels to average, -1 => all
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OTF.PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OTF.PIsA(outOTF):
        raise TypeError,"outOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    Obit.VEGASAverage(inOTF.me, outOTF.me, chAvg, err.me)
    # end PAverage

import OPlot, OTFRec
def PlotSpectrum (inData, targets, scans, feed, Stokes, channs, err, \
              output="None", bgcolor=0):
    """ Plot selected data

    Plot data in inData selected by targets and scans
    inData  = OTF data set to be plotted, any calibration and editing will be applied
    targets = list of target names, empty list = all
    scans   = Range of scan number, 0's => all
    feed    = feed to plot
    Stokes  = Stokes product 0=XX, 1=YY, 2=Real(XY), 3=Imag(XY)
    channs  = [first, last] channels to plot, 1-rel
    err     = Python Obit Error/message stack
    output  = name and type of output device:
              "None"  interactive prompt
              "xwin"  X-Window (Xlib)
              "gcw"   Gnome Canvas Widget (interacts with ObitTalk)
              "ps"    PostScript File (monochrome)
              "psc"   PostScript File (color)
              "xfig"  Fig file
              "png"   PNG file
              "jpeg"  JPEG file
              "gif"   GIF file
              "null"  Null device
    bgcolor   = background color index (1-15), symbolic names:
                BLACK, RED(default), YELLOW, GREEN, 
                AQUAMARINE, PINK, WHEAT, GRAY, BROWN,
                BLUE, BLUEVIOLET, CYAN, TURQUOISE,
                MAGENTA, SALMON, WHITE
    """
    # Set selection
    inInfo = inData.List
    inInfo.set("doCalSelect",True)       # Do selection
    if len(targets)>0:
        inInfo.set("Targets", targets)   # select only target data
    if scans[0]>0 and scans[1]>0:
        lscans = scans
    else:
        lscans = [1,10000000]
    inInfo.set("Scans",lscans)
    inInfo.set("Feeds",[feed])
    inInfo.set("BChan",channs[0])
    inInfo.set("EChan",channs[1])
  
    inData.Open(OTF.READONLY,err)
    d = inData.Desc.Dict
    nrec  = d["nrecord"]
    nchan = d["inaxes"][d["jlocf"]]
    incs  = d["incs"]/2
    incf  = d["incf"]/2
    off = Stokes*incs
    # Init accumulators
    spec = []; wt = []
    for i in range(0,nchan):
        spec.append(0.0)
        wt.append(0.0)
    # Loop over data
    for i in range(1,nrec+1):
        rec = inData.ReadRec (err)
        # eod of file?
        eof = 'EOF' in rec.__dict__
        if eof:
            break
        OErr.printErrMsg(err, "Error reading OTF data")
        for j in range(0,nchan):
            indx = off+j*incf
            if (rec.data[indx][1]!=0.0) and (rec.data[indx][0]!=FArray.fblank):
                spec[j] += rec.data[indx][0] * rec.data[indx][1]
                wt[j]   += rec.data[indx][1]
    # end loop over data
    inData.Close(err)
    
    ch = []
    for j in range(0,nchan):
        ch.append(float(j+channs[0]))

        if wt[j] > 0.0:
            spec[j] /=  wt[j]
            #spec[j] = min(spec[j],1000.0)
        else:   # Hide at 0
            spec[j] = FArray.fblank
            ch[j]   = ch[0]
    # Anything selected
    if len(spec)<=0:
        raise RuntimeError,'No data selected to plot'
    # plot
    plot = OPlot.newOPlot("plot", err, output=output, bgcolor=bgcolor)
    plotInfo = plot.List
    plotInfo.set("XLABEL","Channel")
    plotInfo.set("YLABEL","Flux density")
    plotInfo.set("XMIN",1)
    #plotInfo.set("YMIN",0.0)
    #plotInfo.set("YMAX",1000.0)
    plotInfo.set("TITLE","Poln "+str(Stokes))
    #print spec
    OPlot.PXYPlot(plot, 2, ch, spec, err)
    OPlot.PShow(plot,err)
    OErr.printErrMsg(err, "Error plotting OTF/VEGAS data")
   # end PlotSpectrum

