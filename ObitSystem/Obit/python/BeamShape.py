""" 
Python Obit BeamShape class

This class provides estimates of the beam shape
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2008
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

# Obit BeamShape
import Obit, OErr, Image, InfoList

# Python shadow class to ObitBeamShape class

# class name in C
myClass = "ObitBeamShape"
 
    
class BeamShapePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.BeamShapeUnref(Obit.BeamShape_me_get(self.this))
            # In with the new
            Obit.BeamShape_me_set(self.this,value)
            return
        # members
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != BeamShape:
            return
        if name == "me" : 
            return Obit.BeamShape_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != BeamShape:
            return
        return "<C BeamShape instance> " + Obit.BeamShapeGetName(self.me)
#
class BeamShape(BeamShapePtr):
    """
    Python Obit BeamShape class
    
    This class provides estimates of the beam shape
    
    BeamShape Members with python interfaces:
    """
    def __init__(self, name="no_name", image=None, pbmin=0.05, antSize=25, doGain=True) :
        self.this = Obit.new_BeamShape(name, image, pbmin, antSize, doGain)
        self.myClass = myClass
    def __del__(self):
        if Obit!=None:
            Obit.delete_BeamShape(self.this)
    def cast(self, toClass):
        """ Casts object pointer to specified class
        
        * self     = object whose cast pointer is desired
        * toClass  = Class string to cast to
        """
        ################################################################
        # Get pointer with type of this class
        out =  self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
    
    def Gain (self, ra, dec, parAng=0.0):
        """
        Returns Gain

        * self     = the BeamShape object
        * ra       = RA (deg) of direction for gain
        * dec      = RA (deg) of direction for gain
        * parAng   = Parallactic angle (rad) NYI
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError,"self MUST be a Python Obit BeamShape"
        return Obit.BeamShapeGain(self.me, ra, dec, parAng)
    # end Gain
    
    def GainSym (self, Angle):
        """
        Calculate gain in a given offset from a symmetric beam shape.
        
        Returns Gain, Simple function of distance from pointing center.

        * self     = the BeamShape object
        * Angle    = Angular distance (deg) from pointing center
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError,"self MUST be a Python Obit BeamShape"
        return Obit.BeamShapeGainSym(self.me, Angle)
    # end GainSym
    
    def BeamShapeIsA (self):
        """
        Tells if input really a Python Obit BeamShape
        
        return true, false (1,0)

        * self   = Python BeamShape object
        """
        ################################################################
        # Allow derived types
        return Obit.BeamShapeIsA(self.cast(myClass))
        # end BeamShapeIsA 
    # end class BeamShape
    
def PIsA (inBeamShape):
    """
    Tells if input really a Python Obit BeamShape
    
    return True, False (1,0)

    * inBeamShape   = Python BeamShape object
    """
    ################################################################
    if inBeamShape.__class__ != BeamShape:
        print "Actually",inBeamShape.__class__
        return 0
    # Checks - allow inheritence
    return Obit.BeamShapeIsA(inBeamShape.me)
    # end PIsA

def PCreate (name, image, pbmin, antSize, doGain):
    """
    Create the underlying structures of a BeamShape
    
    return object created.

    * name      = Name to be given to object
    * image     = Python Image for which beam shape is desired
    * pbmin     = Minimum gain, lower values will be clipped at this value
    * antSize   = Size of Antenna in (m)
    * doGain    = If true gain wanted, else gain set to 1.0
    """
    ################################################################
    # Checks
    if not Image.PIsA(image):
        raise TypeError,"uvData MUST be a Python Obit UV"
    #
    out = BeamShape("None",image=image.me,pbmin=pbmin,antSize=antSize,doGain=doGain);
    out.me = Obit.BeamShapeCreate(name, image.me, pbmin, antSize, doGain)
    return out;
    # end PCreate
