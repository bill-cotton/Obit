""" Python Obit FullBeam class

This class provides values of the beam shape derived from an image

"""
# $Id: FullBeam.py 2 2008-06-10 15:32:27Z bill.cotton $
#-----------------------------------------------------------------------
#  Copyright (C) 2009
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

# Obit FullBeam
import Obit, OErr, Image, InfoList

# Python shadow class to ObitFullBeam class

# class name in C
myClass = "ObitFullBeam"
 
    
class FullBeamPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.FullBeamUnref(Obit.FullBeam_me_get(self.this))
            # In with the new
            Obit.FullBeam_me_set(self.this,value)
            return
        # members
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != FullBeam:
            return
        if name == "me" : 
            return Obit.FullBeam_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != FullBeam:
            return
        return "<C FullBeam instance> " + Obit.FullBeamGetName(self.me)
#
class FullBeam(FullBeamPtr):
    """ Python Obit FullBeam class
    
    This class provides values of the beam shape derived from an image
    Gains at specified offsets from the beam center at giv3en parallactic
    angles are interpolated from a Full Beam image.
    
    FullBeam Members with python interfaces:
    """
    def __init__(self, name="no_name", image=None, err=None) :
        self.this = Obit.new_FullBeam(name, image, err)
        self.myClass = myClass
    def __del__(self):
        if Obit!=None:
            Obit.delete_FullBeam(self.this)
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
    
    def Gain (self, dra, ddec, parAng, plane, err):
        """ 
        
        Returns Gain
        self     = the FullBeam object
        ra       = RA (deg) offset of direction for gain
        dec      = Dec (deg) offset of direction for gain
        parAng   = Parallactic angle (deg)
        plane    = plane in FullBeam (from FindPlane)
        err      = Obit error/message stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError,"self MUST be a Python Obit FullBeam"
        return Obit.FullBeamValue(self.me, dra, ddec, parAng, plane, err.me)
    # end Gain
    
    def FindPlane(self, freq):
        """ 
        
        Returns nearest plane to a given frequency
        self     = the FullBeam object
        freq     = frequency (Hz)
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError,"self MUST be a Python Obit FullBeam"
        return Obit.FullBeamFindPlane(self.me, freq)
    # end FindPlane
    
    def FullBeamIsA (self):
        """ Tells if input really a Python Obit FullBeam
        
        return true, false (1,0)
        self   = Python FullBeam object
        """
        ################################################################
        # Allow derived types
        return Obit.FullBeamIsA(self.cast(myClass))
        # end FullBeamIsA 
    # end class FullBeam
    
def PIsA (inFullBeam):
    """ Tells if input really a Python Obit FullBeam

    return True, False (1,0)
    inFullBeam   = Python FullBeam object
    """
    ################################################################
    if inFullBeam.__class__ != FullBeam:
        print "Actually",inFullBeam.__class__
        return 0
    # Checks - allow inheritence
    return Obit.FullBeamIsA(inFullBeam.me)
    # end PIsA

def PCreate (name, image, err):
    """ Create the underlying structures of a FullBeam

    return object created.
    name      = Name to be given to object
    image     = Python Image for which beam shape is desired
    err       = Obit error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(image):
        raise TypeError,"image MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    out = FullBeam("None",image=image.me,err=err.me);
    out.me = Obit.FullBeamCreate(name, image.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating beam object")
    return out;
    # end PCreate

