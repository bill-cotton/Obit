""" Python Obit FitModel class

This class contains a parameterized image model component

FitModel Members with python interfaces:
name   An optional name for the object.
type   Model type of the model component:
   PointMod    Point
   GaussMod    Eliptical Gaussian
   USphereMod  Uniform optically thin sphere
   Background  Background wedge
Peak   Peak density
DeltaX "X" (RA) offset (deg) of center from reference position
DeltaY "Y" (Dec) offset (deg) of center from reference position
nparm  Number of parameters
parms  Model parameters, type dependent 
ePeak   Error in Peak density
eDeltaX Error in "X" (RA) offset (deg) of center from reference position
eDeltaY Error in "Y" (Dec) offset (deg) of center from reference position
eparms  Error in Model parameters, type dependent 
"""
# $Id: FitModel.py,v 1.4 2007/10/11 13:38:33 bcotton Exp $
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

# Obit FitModel
import Obit, OErr, ImageMosaic, InfoList, UV, ImageDesc, SkyGeom

# Python shadow class to ObitFitModel class

# class name in C
myClass = "ObitFitModel"
 
# Class data model type codes
PointMod   = 0   # Point
GaussMod   = 1   # Eliptical Gaussian
USphereMod = 2   # Uniform optically thin sphere
Background = 3   # Background wedge
    
class FitModelPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.FitModelUnref(Obit.FitModel_me_get(self.this))
            # In with the new
            Obit.FitModel_me_set(self.this,value)
            return
        # members
        if name=="type":
            Obit.FitModelSetType(self.me,value)
            return
        if name=="Peak":
            Obit.FitModelSetPeak(self.me,value)
            return
        if name=="DeltaX":
            Obit.FitModelSetDeltaX(self.me,value)
            return
        if name=="DeltaY":
            Obit.FitModelSetDeltaY(self.me,value)
            return
        if name=="nparm":
            Obit.FitModelSetNparm(self.me,value)
            return
        if name=="parms":
            Obit.FitModelSetNparm(self.me,len(value))
            Obit.FitModelSetParms(self.me,value)
            return
        if name=="ePeak":
            Obit.FitModelSetePeak(self.me,value)
            return
        if name=="eDeltaX":
            Obit.FitModelSeteDeltaX(self.me,value)
            return
        if name=="eDeltaY;":
            Obit.FitModelSeteDeltaY(self.me,value)
            return
        if name=="eparms":
            Obit.FitModelSeteParms(self.me,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != FitModel:
            return
        if name == "me" : 
            return Obit.FitModel_me_get(self.this)
        # members
        if name=="type":
            return Obit.FitModelGetType(self.me)
        if name=="Peak":
            return Obit.FitModelGetPeak(self.me)
        if name=="DeltaX":
            return Obit.FitModelGetDeltaX(self.me)
        if name=="DeltaY":
            return Obit.FitModelGetDeltaY(self.me)
        if name=="nparm":
            return Obit.FitModelGetNparm(self.me)
        if name=="parms":
            return Obit.FitModelGetParms(self.me)
        if name=="ePeak":
            return Obit.FitModelGetePeak(self.me)
        if name=="eDeltaX":
            return Obit.FitModelGeteDeltaX(self.me)
        if name=="eDeltaY":
            return Obit.FitModelGeteDeltaY(self.me)
        if name=="eparms":
            return Obit.FitModelGeteParms(self.me)
    def __repr__(self):
        if self.__class__ != FitModel:
            return
        return "<C FitModel instance> " + Obit.FitModelGetName(self.me)
#
class FitModel(FitModelPtr):
    """ Python Obit FitModel class
    
    This class contains a parameterized image model
    
    FitModel Members with python interfaces:
    """
    def __init__(self, name="no_name", type=PointMod, Peak=0.0, DeltaX=0.0, DeltaY=0.0, parms=[]) :
        self.this = Obit.new_FitModel(name, type, Peak, DeltaX, DeltaY, len(parms), parms)
        self.myClass = myClass
    def __del__(self):
        if Obit!=None:
            Obit.delete_FitModel(self.this)
    def cast(self, toClass):
        """ Casts object pointer to specified class
        
        self     = object whose cast pointer is desired
        toClass  = Class string to cast to
        """
        ################################################################
        # Get pointer with type of this class
        out =  self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
    
    def DeconGau (self, ImDesc):
        """ Deconvolves Beam on an image descriptor from a Gaussian model component
        
        Returns a FitModel with the deconvolved values
        self     = object with Gaussian to be desired
        ImDesc   = Image Descriptor with Beam
        """
        ################################################################
        try:
            out = FitModel("Deconvolved")
            out.nparm = 3
            Obit.DeconGau(self.me, out.me, ImDesc.me)
            #print "dconv OK"
        except:
            #print "dconv failed"
            out = None
        else:
            #print "dconv failed2"
            #out = None
            pass
        return out
    # end DeconGau
    
    def Print (self, ImDesc, corner, file=None):
        """ Prepare human readable contents

        Returns string with description of model
        self     = object with Model to display
        ImDesc   = Image Descriptor with Beam, etc.
        corner   = bottom left corner in selected region of image (0-rel)
        """
        ################################################################
        # Start output string
        id = ImDesc.Dict
        modelInfo = ""
        
        # Collect info
        type    = self.type
        parms   = self.parms
        eparms  = self.eparms

        # Gaussian
        if type==GaussMod:
            # Get celestial position
            xpix = corner[0]+self.DeltaX; ypix = corner[1]+self.DeltaY
            pos=SkyGeom.PWorldPos (xpix, ypix, id["crval"][0], id["crval"][1], \
                                   id["crpix"][0], id["crpix"][1], \
                                   id["cdelt"][0], id["cdelt"][1], id["crota"][1], \
                                   id["ctype"][0][4:8], 0.0, 0.0)
            rast  = ImageDesc.PRA2HMS(pos[1])
            decst = ImageDesc.PDec2DMS(pos[2])
            # Position errors
            era  = self.eDeltaX * abs(id["cdelt"][0])*3600.0
            edec = self.eDeltaY * abs(id["cdelt"][1])*3600.0
            modelInfo += "RA  "+rast+" (%8.3g asec),  pixel %8.3f (%8.3g)\n" % (era, xpix, self.eDeltaX)
            modelInfo += "Dec  "+decst+" (%8.3g asec),  pixel %8.3f (%8.3g)\n" % (edec, ypix, self.eDeltaY)
            
            modelInfo += "Peak Flux density %8.3g (%8.3g) %s\n" % (self.Peak, self.ePeak, id["bunit"])
            # Ratio of Gaussian beams if beam in ImDesc
            if (id["beamMaj"]>0.0) and (id["beamMin"]>0.0):
                ratio = (parms[0] * abs(id["cdelt"][0])*parms[1] * abs(id["cdelt"][1])) / \
                        (id["beamMaj"]*id["beamMin"])
                modelInfo += "Integrated Flux density %8.3g (%8.3g) %s\n" % \
                             (self.Peak*ratio, self.ePeak*ratio, "Jy")
                
                modelInfo += "Fitted Major axis %8.3f (%8.3g) asec, %8.3f (%8.3g) pixels\n" % \
                             (parms[0]*abs(id["cdelt"][0])*3600.0, eparms[0]*abs(id["cdelt"][0])*3600.0, \
                              parms[0], eparms[1])
                modelInfo += "Fitted Minor axis %8.3f (%8.3g) asec, %8.3f (%8.3g) pixels\n" % \
                             (parms[1]*abs(id["cdelt"][1])*3600.0, eparms[1]*abs(id["cdelt"][1])*3600.0, \
                              parms[1], eparms[1])
                modelInfo += "Fitted Position angle %8.5g (%8.3g) deg\n" % \
                             (parms[2]*57.296, eparms[2]*57.296)
                # Deconvolve
                deconMod = self.DeconGau(ImDesc)
                if PIsA(deconMod) and deconMod.type>0:
                    modelInfo += "\nDeconvolved model\n"
                    dparms  = deconMod.parms
                    deparms = deconMod.eparms
                    modelInfo += "Deconvolved Major axis %8.3g (%8.3g) asec, %8.3f (%8.3g) pixels\n" % \
                                 (dparms[0]*abs(id["cdelt"][0])*3600.0, deparms[0]*abs(id["cdelt"][0])*3600.0, \
                                  dparms[0], deparms[1])
                    modelInfo += "Deconvolved Minor axis %8.3g (%8.3g) asec, %8.3f (%8.3g) pixels\n" % \
                                 (dparms[1]*abs(id["cdelt"][1])*3600.0, deparms[1]*abs(id["cdelt"][1])*3600.0, \
                                  dparms[1], deparms[1])
                    modelInfo += "Deconvolved Position angle %8.5g (%8.3g) deg\n" % \
                                 (dparms[2]*57.296, deparms[2]*57.296)
           
                else:
                    modelInfo += "\nDeconvolution failed\n"
                    # end deconvolved
        # end Gaussian
        # done
        return modelInfo
        # end Print

    # end class FitModel
    
def PIsA (inFitModel):
    """ Tells if input really a Python Obit FitModel

    return True, False (1,0)
    inFitModel   = Python FitModel object
    """
    ################################################################
    if inFitModel.__class__ != FitModel:
        print "Actually",inFitModel.__class__
        return 0
    # Checks - allow inheritence
    return Obit.FitModelIsA(inFitModel.me)
    # end PIsA
