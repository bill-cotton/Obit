# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2024,2025
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
 
# Python class for image slices
import Obit, OErr, Image, FArray, FArrayUtil, FInterpolate, OPlot, OWindow, ODisplay, ImageDesc

class Slice(Obit.Source):
    """ Python Obit Image slice class
    
    """
    def __init__(self, name, img):
        self.__dict__ = {"Name":name, "Image":img, "Pos1":[0,0], "Pos2":[0,0], 
                         "incr":0.3, "Slice":[]}; 

    def set_corn (self, disp, err):
        """
        Set ends of slice using TV window edit
        Only seems to work on port 8765
        """
        win = OWindow.PCreate1('win',self.Image.Desc.Dict['inaxes'][0:2], err)
        print("Set rectangular box for BLC,TRC as slice ends")
        ODisplay.PImage(disp, self.Image, err, window=win)
        OErr.printErr(err)
        lll=OWindow.PGetList(win,1,err)
        if len(lll[0])<6:
            print ('editing failed, window:',lll)
            raise RuntimeError("Window editing failed")
        self.Pos1 = lll[0][2:4]; self.Pos2 = lll[0][4:]
    # end set_corn
    
    def fill (self, err, plane=1, scale=1.0):
        """
        Calculate slice values interpolating image between Pos1 and Pos2
        plane = plane in image
        scale = scaling factor
        """
        delta_x = (self.Pos2[0]-self.Pos1[0])
        delta_y = (self.Pos2[1]-self.Pos1[1])
        delta = (((delta_x)**2 + (delta_y)**2)**0.5)
        npoint = int (1+delta/self.incr)
        # Actual increments
        delta_x /= npoint; delta_y /= npoint; 
        #print (npoint,delta,delta_x,delta_y)
        img = self.Image
        img.GetPlane(None,[plane,1,1,1,1],err)
        fi = FInterpolate.FInterpolate("FI", img.FArray, img.Desc, 2)
        self.Slice = []
        for i in range(0,npoint):
            pixel = [self.Pos1[0]+i*delta_x,self.Pos1[1]+i*delta_y]
            val = scale*FInterpolate.PPixel(fi, pixel, err)
            #print (i,pixel,val)
            self.Slice.append(val)
    # end fill
    
    def fill_FA (self, err, plane=1, scale=1.0):
        """
        Calculate slice values interpolating image between Pos1 and Pos2
        And return as an FArray
        plane = plane in image
        scale = scaling factor
        """
        # Fill slice
        self.fill(err,plane,scale)
        npos = [len(self.Slice)]
        # create FArray
        fa = FArray.FArray("Slice", npos)
        for i in range(0,npos[0]):
            fa.set(self.Slice[i],i)
        return fa
    # end fill_FA
    
    def fit_Gauss (self, fa, err, FWHM=[0.], center=[0.], peak=[0.]):
        """
        Fits a Gaussian to a slice in an FArray
        fa     = FArray with pixel data
        returns  list:
        [0]= Full width Half Max (pixels) of fitted Gaussian
        [1]= peak value in fitted Gaussian
        [2]= x pixel (0-rel) coordinates of peak in pixels
        [3]= 0th order baseline term
        [4]= 1st order baseline term
        [5]= RMS residual        """
        # Fit
        res = FArrayUtil.PFit1DGauss2(fa,len(FWHM),err,FWHM,center,peak)
        return res
    # end fit_Gauss
    
    def clip (self, ref_slice, prange, val, err):
        """
        Replace slice values with corresponding ref_slice values
        out of prange with val.
        ref_slice = Slice to use to clip self.Slice
        prange    = [min, max] values to accept
        val       = value for replacement
        """
        # Check compatibility
        ncheck = len(ref_slice.Slice);
        if len(self.Slice) != ncheck:
            print ('Slices are incompatible: %d %d', \
                   len(self.Slice), ncheck)
            raise RuntimeError("Incompatible slices")
        for i in range(0,ncheck):
            if (ref_slice.Slice[i]<prange[0]) or (ref_slice.Slice[i]>prange[1]):
                self.Slice[i] = val
    # end clip
    
    def plot (self, pname, err, label='',prange=None, ylab = None):
        """
        Plot to postscript file pname
        pname  = name of postscript file
        label  = prefix for title
        prange = [min, max] values to plot
        ylab   = Y axis label
        """
        delta_x = self.incr * self.Image.Desc.Dict['cdelt'][1]*3600
        x = []; y = []
        # Scaling
        peak = max(max(self.Slice),abs(min(self.Slice)))
        scale = 1.0
        if peak<1.0e-1:
            scale = 1.0e3; pre = "m"
        if peak<1.0e-4:
            scale = 1.0e6; pre = "#gm"
        for i in range(0,len(self.Slice)):
            x.append((i-len(self.Slice)/2)*delta_x)
            y.append(scale*self.Slice[i])
        plot = OPlot.newOPlot('plt', err, output=pname+'/ps')
        plot.List.set("XLABEL","Position along slice (asec)")
        if ylab:
            plot.List.set("YLABEL",ylab)
        else:
            plot.List.set("YLABEL","Flux Density ("+pre+"Jy)")
        # Title giving positions
        pix1 = [float(self.Pos1[0]),float(self.Pos1[1])]
        pos1 = ImageDesc.PGetPos(self.Image.Desc,pix1,err)
        pix2 = [float(self.Pos2[0]),float(self.Pos2[1])]
        # Center to pos 2
        pix2 = [(pix1[0]+pix2[0])*0.5,(pix1[1]+pix2[1])*0.5]
        pos2 = ImageDesc.PGetPos(self.Image.Desc,pix2,err)
        ra2 = ImageDesc.PRA2HMS(pos2[0])[2:13]+' '; 
        dec2 = ImageDesc.PDec2DMS(pos2[1])[0:11]+' '
        title = label+ra2+dec2+str(self.Pos1)+'-'+str(self.Pos2) # Make fit
        plot.List.set("TITLE",title)
        if prange:
            plot.List.set("YMIN",scale*prange[0])
            plot.List.set("YMAX",scale*prange[1])
            for i in range(0,len(y)):
                y[i] = max(y[i], scale*prange[0])
                y[i] = min(y[i], scale*prange[1])
        OPlot.PXYPlot(plot,0, x,y, err)
        OPlot.PShow(plot,err)
# end class
