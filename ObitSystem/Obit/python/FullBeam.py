""" Python Obit FullBeam class

This class provides values of the beam shape derived from an image
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2009-2020
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
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, OErr, Image, InfoList

# Python shadow class to ObitFullBeam class

# class name in C
myClass = "ObitFullBeam"
     
#
class FullBeam(Obit.FullBeam):
    """
    Python Obit FullBeam class
    
    This class provides values of the beam shape derived from an image
    Gains at specified offsets from the beam center at giv3en parallactic
    angles are interpolated from a Full Beam image.
    
    FullBeam Members with python interfaces:
    """
    def __init__(self, name="no_name", image=None, err=None) :
        super(FullBeam, self).__init__()
        Obit.CreateFullBeam(self.this, name, image, err)
        self.myClass = myClass
    def __del__(self, DeleteFullBeam=_Obit.DeleteFullBeam):
        if _Obit!=None:
            DeleteFullBeam(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            if self.this!=None:
                Obit.FullBeamUnref(Obit.FullBeam_Get_me(self.this))
            # In with the new
            Obit.FullBeam_Set_me(self.this,value)
            return
        # members
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, FullBeam):
            return "Bogus dude "+str(self.__class__)
        if name == "me" : 
            return Obit.FullBeam_Get_me(self.this)
        raise AttributeError(name)
    def __repr__(self):
        if not isinstance(self, FullBeam):
            return "Bogus dude "+str(self.__class__)
        return "<C FullBeam instance> " + Obit.FullBeamGetName(self.me)
    
    def Gain (self, dra, ddec, parAng, plane, err):
        """
        Returns Gain

        * self     = the FullBeam object
        * ra       = RA (deg) offset of direction for gain
        * dec      = Dec (deg) offset of direction for gain
        * parAng   = Parallactic angle (deg)
        * plane    = plane in FullBeam (from FindPlane)
        * err      = Obit error/message stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit FullBeam")
        return Obit.FullBeamValue(self.me, dra, ddec, parAng, plane, err.me)
    # end Gain
    
    def FindPlane(self, freq):
        """
        Returns nearest plane to a given frequency

        * self     = the FullBeam object
        * freq     = frequency (Hz)
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit FullBeam")
        return Obit.FullBeamFindPlane(self.me, freq)
    # end FindPlane
    
    def FullBeamIsA (self):
        """
        Tells if input really a Python Obit FullBeam
        
        return True, False

        * self   = Python FullBeam object
        """
        ################################################################
        # Allow derived types
        return Obit.FullBeamIsA(self.me)!=0
        # end FullBeamIsA 
    # end class FullBeam
    
def PIsA (inFullBeam):
    """
    Tells if input really a Python Obit FullBeam
    
    return True, False

    * inFullBeam   = Python FullBeam object
    """
    ################################################################
    if not isinstance(inFullBeam, FullBeam):
        print("Actually",inFullBeam.__class__)
        return False
    # Checks - allow inheritence
    return Obit.FullBeamIsA(inFullBeam.me)!=0
    # end PIsA

def PCreate (name, image, err):
    """
    Create the underlying structures of a FullBeam
    
    return object created.

    * name      = Name to be given to object
    * image     = Python Image for which beam shape is desired
    * err       = Obit error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(image):
        raise TypeError("image MUST be a Python Obit Image")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    #
    out = FullBeam("None",image=image.me,err=err.me);
    out.me = Obit.FullBeamCreate(name, image.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating beam object")
    return out;
    # end PCreate

