""" Python Obit ImageInterp class

This class provides values of the beam shape derived from an image
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2009,2019
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

# Obit ImageInterp
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, OErr, Image, InfoList

# Python shadow class to ObitImageInterp class

# class name in C
myClass = "ObitImageInterp"

class ImageInterp(Obit.ImageInterp):
    """
    Python Obit ImageInterp class
    
    This class interpolates pixel values in an image, possibly selecting
    by frequency and optionally rotating by a "parallactic angle".
    
    ImageInterp Members with python interfaces:
    """
    def __init__(self, name="no_name", image=None, hwidth=1, err=None) :
        super(ImageInterp, self).__init__()
        Obit.CreateImageInterp(self.this, image.me, hwidth, err.me )
        self.myClass = myClass
    def __del__(self, DeleteImageInterp=_Obit.DeleteImageInterp):
        if _Obit!=None:
            DeleteImageInterp(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            if self.this!=None:
                Obit.ImageInterpUnref(Obit.ImageInterp_Get_me(self.this))
            # In with the new
            Obit.ImageInterp_Set_me(self.this,value)
            return
        # members
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, ImageInterp):
            return "Bogus Dude"+str(self.__class__)
        if name == "me" : 
            return Obit.ImageInterp_Get_me(self.this)
        raise AttributeError(name)
    def __repr__(self):
        if not isinstance(self, ImageInterp):
            return "Bogus Dude"+str(self.__class__)
        return "<C ImageInterp instance> " + Obit.ImageInterpGetName(self.me)
    def Value (self, ra, dec, err, rotate=0.0, plane=0 ):
        """
        Returns Interpolated pixel value

        * self     = the ImageInterp object
        * ra       = RA (or "X") (deg) coordinate
          Use ImageDesc.PMS2RA for convert from human form
        * dec      = Dec (or "Y") (deg) coordinate
          Use ImageDesc.PDMS2Dec for convert from human form
        * rotate   = (Parallactic) angle to rotate image by (deg)
        * plane    = plane in ImageInterp (from FindPlane)
        * err      = Obit error/message stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit ImageInterp")
        return Obit.ImageInterpValue(self.me, ra, dec, rotate, plane, err.me)
    # end Value
    
    def FindPlane(self, freq):
        """
        Returns nearest plane to a given frequency

        * self     = the ImageInterp object
        * freq     = frequency (Hz)
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit ImageInterp")
        return Obit.ImageInterpFindPlane(self.me, freq)
    # end FindPlane
    
    def ImageInterpIsA (self):
        """
        Tells if input really a Python Obit ImageInterp
        
        return True, False

        * self   = Python ImageInterp object
        """
        ################################################################
        return Obit.ImageInterpIsA(self.me)!=0
        # end ImageInterpIsA 
    # end class ImageInterp
    
def PIsA (inImageInterp):
    """
    Tells if input really a Python Obit ImageInterp
    
    return True, False

    * inImageInterp   = Python ImageInterp object
    """
    ################################################################
    if not isinstance(inImageInterp, ImageInterp):
        print("Actually",inImageInterp.__class__)
        return False
    return Obit.ImageInterpIsA(inImageInterp.me)!=0
    # end PIsA

def PCreate (name, image, err, hwidth=1):
    """
    Create the underlying structures of a ImageInterp
    
    return object created.

    * name      = Name to be given to object
    * image     = Python Image for which beam shape is desired
    * hwidth    = half width of interpolation (usually 1 or 2)
    * err       = Obit error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(image):
        raise TypeError("image MUST be a Python Obit Image")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    #
    out = ImageInterp("None",image=image,hwidth=hwidth,err=err);
    out.me = Obit.ImageInterpCreate(name, image.me, hwidth, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating beam object")
    return out;
    # end PCreate

