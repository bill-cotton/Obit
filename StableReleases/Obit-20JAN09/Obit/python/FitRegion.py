""" Python Obit FitModel class

This class contains info about an image fitting region

FitModel Members with python interfaces:
name      An optional name for the object.
corner    bottom left corner in selected region of image (0-rel)
dim       dimension of region
peak      peak in region 
peakResid peak in region residual after model subtraction
RMSResid  RMS residual
fluxResid Sum of pixel values in residual
nmodel    Number of models
models    Array of FitModels
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007
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

# Obit FitRegion
import Obit, OErr, FitModel, Image, ODisplay, OWindow, FInterpolate
import pydoc

# Python shadow class to ObitFitRegion class

# class name in C
myClass = "ObitFitRegion"
 
    
class FitRegionPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.FitRegionUnref(Obit.FitRegion_me_get(self.this))
            # In with the new
            Obit.FitRegion_me_set(self.this,value)
            return
        # members
        if name=="corner":
            Obit.FitRegionSetCorner(self.me, value)
            return
        if name=="dim":
            Obit.FitRegionSetDim(self.me, value)
            return
        if name=="peak":
            Obit.FitRegionSetPeak(self.me, value)
            return
        if name=="peakResid":
            Obit.FitRegionSetPeakResid(self.me, value)
            return
        if name=="RMSResid":
            Obit.FitRegionSetRMSResid(self.me, value)
            return
        if name=="fluxResid":
            Obit.FitRegionSetFluxResid(self.me, value)
            return
        if name=="nmodel":
            Obit.FitRegionSetNmodel(self.me, value)
            return
        if name=="models":
            nmodel = Obit.FitRegionGetNmodel(self.me)
            # Resize if necessary
            if nmodel!=len(value):
                self.nmodel = len(value)
            # Pass one at a time
            i=0
            for x in value:
                Obit.FitRegionSetModels(self.me, x.me, i)
                i+= 1;
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != FitRegion:
            return
        if name == "me" : 
            return Obit.FitRegion_me_get(self.this)
        # members
        if name=="corner":
            return Obit.FitRegionGetCorner(self.me)
        if name=="dim":
            return Obit.FitRegionGetDim(self.me)
        if name=="peak":
            return Obit.FitRegionGetPeak(self.me)
        if name=="peakResid":
           return  Obit.FitRegionGetPeakResid(self.me)
        if name=="RMSResid":
            return Obit.FitRegionGetRMSResid(self.me)
        if name=="fluxResid":
            return Obit.FitRegionGetFluxResid(self.me)
        if name=="nmodel":
            return Obit.FitRegionGetNmodel(self.me)
        if name=="models":
            # Make into array of python structures
            out = []
            nmodel = Obit.FitRegionGetNmodel(self.me)
            for i in range(0,nmodel):
                tout = FitModel.FitModel("model")
                tout.me = Obit.FitRegionGetModels(self.me, i)
                out.append(tout)
            return out
    def __repr__(self):
        if self.__class__ != FitRegion:
            return
        return "<C FitRegion instance> " + Obit.FitRegionGetName(self.me)
#
class FitRegion(FitRegionPtr):
    """ Python Obit FitRegion class
    
    This class contains info about an image fitting region
    Interactive definitions of a region using the Image viewer
    is possible using PSetup to create the Fit Region.
    
    FitRegion Members with python interfaces:
    """
    def __init__(self, name="no_name", corner=[0,0], dim=[0,0], peak=0.0, \
                 peakResid=0.0, RMSResid=0.0, fluxResid=0.0, models=[]):
        self.this = Obit.new_FitRegion(name, corner, dim, peak, peakResid, \
                                       RMSResid, fluxResid)
        if len(models)>0:  # Set models if given
            self.nmodel = len(models)
            self.models = models
        self.myClass = myClass
    def __del__(self):
        if Obit!=None:
            Obit.delete_FitRegion(self.this)

    def Print (self, ImDesc, file=None):
        """ Display human readable contents
        
        self     = object with Model to display
        ImDesc   = Image Descriptor with Beam, etc.
        file     = if present, the name of a file into which to write
        the information rather than displaying it on the screen
        """
        ################################################################
        # Start output string
        id = ImDesc.Dict
        modelInfo = "\nModel fit for "+id["object"]+"\n"
        
        # Loop over models
        corner = self.corner
        for imod in self.models:
            modelInfo += imod.Print(ImDesc, corner)

        # Display or log
        if file:
            fd = open(file, "a")
            fd.write(modelInfo)
            fd.close()
        else:
            #pydoc.ttypager(modelInfo)
            print modelInfo
        del modelInfo
        # end Print
        
        
    # end class FitRegion
    
def PSetup (inImage, disp, err):
    """ Interactive initial definition of fitting region

    Interactively allows the user to set the region of the image
    to be fitted and the initial model.
    The fitting region is first specified with a rectangular window
    and then the initial models to be fitted with circular windows.
    Returns FitRegion, leaves image pixel array on inImage
    image  = image to be fitted
    disp   = image display to use
    err    = Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not ODisplay.PIsA(disp):
        print "Actually ",disp.__class__
        raise TypeError,"disp MUST be a Python Obit Display"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    id = inImage.Desc.Dict
    naxis = id["inaxes"][0:2]
    window = OWindow.PCreate1("window", naxis, err)
    # Give user instructions
    print "Specify the region to fit with a rectangular box."
    print "Followed by circular boxes to mark Gaussian component"
    print "initial locations and initial sizes."
    boxes = []
    while len(boxes)<2:
        ODisplay.PImage(disp, inImage, err, window)
        if err.isErr:
            printErrMsg(err, "Error with interactive display")
            
        boxes = OWindow.PGetList(window, 1, err)
        if err.isErr:
            printErrMsg(err, "Error getting boxes")
            print "DEBUG BOXES",boxes
            
        # Checks
        if len(boxes)<2:
            print "You MUST specity fitting region and at least one Gaussian"
        # End loop 'til user gets it right
    # Fitting window
    box = boxes[0]
    if box[1]!=OWindow.RectangleType:
        raise RuntimeError,"Fitting region NOT rectangular"
    corner = [min(box[2],box[4]), min(box[3],box[5])]
    dim    = [abs(box[4]-box[2])+1, abs(box[5]-box[3])+1]
    # Image interpolator
    inImage.Open(Image.READONLY,err)
    inImage.Read(err)
    fi = FInterpolate.FInterpolate("inter", inImage.FArray, inImage.Desc, 2)
    inImage.Close(err)
    if err.isErr:
        printErrMsg(err, "Error reading image for interpolator")
    # Models
    models = []
    for box in boxes[1:]:
        if box[1]!=OWindow.RoundType:
            raise RuntimeError,"Indicate Gaussians with Round boxes"
        type = FitModel.GaussMod
        x = box[3] - corner[0]
        y = box[4] - corner[1]
        pixel = [float(box[3]), float(box[4])]
        parms = [float(box[2]), float(box[2]), 0.0]
        s = FInterpolate.PPixel(fi, pixel,err)
        models.append(FitModel.FitModel(type=type, Peak=s, DeltaX=x, DeltaY=y, parms=parms))

    # create output 
    out = FitRegion(name=id["object"], corner=corner, dim=dim, models=models)
    return out
    # end PSetup

def PIsA (inFitRegion):
    """ Tells if input really a Python Obit FitRegion

    return True, False (1,0)
    inFitRegion   = Python FitRegion object
    """
    ################################################################
    if inFitRegion.__class__ != FitRegion:
        print "Actually is",inFitRegion.__class__ 
        return 0
    return Obit.FitRegionIsA(inFitRegion.me)
    # end PIsA
