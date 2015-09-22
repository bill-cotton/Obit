# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2015
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

# Python shadow class to ObitFInterpolate class
import Obit, FArray, Image, ImageDesc, InfoList, OErr

class FInterpolatePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.FInterpolate_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.FInterpolate_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C FInterpolate instance>"
class FInterpolate(FInterpolatePtr):
    """
    Lagrangian interpolation in an FArray
    """
    def __init__(self, name, array, desc, hwidth) :
        self.this = Obit.new_FInterpolate(name, array.me, desc.me, hwidth)
    def __del__(self):
        if Obit!=None:
            Obit.delete_FInterpolate(self.this)
    
def PCopy  (inFI, outFI, err):
    """ 
    Make a deep copy of input object.

    * inFI    = Python Obit input FInterpolate
    * outFI   = Python Obit FInterpolate
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    if not PIsA(outFI):
        raise TypeError,"outFI MUST be a Python Obit FInterpolate"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    Obit.FInterpolateCopy (inFI.me, outFI.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error copying FInterpolate")
    # end PCopy

def PClone (inFI, outFI):
    """
    Make a shallow copy of a object (no data copied)

    * inFI    = Python Obit input FInterpolate
    * outFI   = Python Obit FInterpolate
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    if not PIsA(outFI):
        raise TypeError,"outFI MUST be a Python Obit FInterpolate"
    #
    Obit.FInterpolateClone (inFI.me, outFI.me)
    # end PClone 

def PReplace (inFI, newArray):
    """
    Replace the ObitFArray member to be interpolated.

    * inFI    = Python Obit input FInterpolate
    * newArray= Python Obit FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    if not FArray.PIsA(newArray):
        raise TypeError,"newArray MUST be a Python Obit FArray"
    #
    Obit.FInterpolateReplace (inFI.me, newArray.me)
    # end PReplace


def PPixel (inFI, pixel, err):
    """
    Interpolate pixel value
    
    return the value of a specified pixel

    * inFI    = Python Obit input FInterpolate
    * pixel   = pixel (1-rel) as array of float
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    ret = Obit.FInterpolatePixel (inFI.me, pixel, err.me)
    if err.isErr:
        printErrMsg(err, "Error interpolating pixel value")
    return ret
    # end PPixel

def P1D (inFI, pixel, err):
    """
    Interpolate pixel value in 1D array.
    
    return the value of a specified pixel

    * inFI    = Python Obit input FInterpolate
    * pixel   = pixel number as float
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    #
    return Obit.FInterpolate1D (inFI.me, pixel)
    # end P1D

def PPosition (inFI, coord, err):
    """
    Interpolate value at requested coordinate in array.
    
    Interpolate value at requested coordinate in array.
    The object must have an image descriptor to allow determing
    pixel coordinates.
    Interpolation between planes is not supported.
    return the value at specified coordinate

    * inFI    = Python Obit input FInterpolate
    * coord   = coordinate as array of float
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return None
    #
    ret = Obit.FInterpolatePosition (inFI.me, coord, err.me)
    if err.isErr:
        printErrMsg(err, "Error interpolating pixel value")
    return ret
    # end PPosition

def PGetList (inFI):
    """
    Get InfoList
    
    return InfoList reference

    * inFI    = Python Obit input FInterpolate
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    #
    out = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.FInterpolateGetList(inFI.me)
    return out
    # end PGetList


def PGetFArray (inFI):
    """
    Get Associated FArray reference
    
    return FArray (data being interpolated) reference 

    * inFI    = Python Obit input FInterpolate
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    #
    out    = FArray.FArray("None")
    out.me = Obit.FInterpolateGetFArray (inFI.me)
    return out
    # end PGetFArray 


def PGetDesc (inFI):
    """
    Get Image descriptor
    
    return ImageDesc reference

    * inFI    = Python Obit input FInterpolate
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
       raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    #
    out    = ImageDesc.ImageDesc("None")
    out.me = Obit.FInterpolateGetImageDesc (inFI.me)
    return out
    # end PGetDesc


def PSetDesc (inFI, desc):
    """
    Replace Image descriptor

    * inFI    = Python Obit input FInterpolate
    * desc    = Python Obit ImageDesc to use
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    if not ImageDesc.PIsA(desc):
        raise TypeError,"desc MUST be a Python Obit ImageDesc"
    #
    Obit.FInterpolateSetDesc (inFI.me, desc.me)
    # end PSetDesc


def PGetHwidth (inFI):
    """
    return Half width of interpolation kernal
    
    return hwidth member value

    * inFI    = Python Obit input FInterpolate
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    #
    return Obit.FInterpolateGetHwidth(inFI.me)
    # end PGetHwidth


def PSetHwidth (inFI, hwidth):
    """
    Set Half width of interpolation kernal

    * inFI    = Python Obit input FInterpolate
    * hwidth  = new half width of interpolation kernal
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit FInterpolate"
    #
    Obit.FInterpolateSetHwidth(inFI.me, hwidth)
    # end PSetHwidth


def PIsA (inFI):
    """
    Tells if object thinks it's a Python ObitFInterpolate
    
    return true, false (1,0)

    * inFI    = Python Obit input FInterpolate to test
    """
    ################################################################
    # Checks
    if inFI.__class__ !=  FInterpolate:
        return 0
    return Obit.FInterpolateIsA(inFI.me);

    # end PIsA
