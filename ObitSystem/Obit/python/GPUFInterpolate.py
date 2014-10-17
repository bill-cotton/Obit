# $Id: $
#-----------------------------------------------------------------------
#  Copyright (C) 2014
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

# Python shadow class to ObitGPUFInterpolate class
import Obit, FArray, ImageDesc, InfoList, OErr

class GPUFInterpolatePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.GPUFInterpolate_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.GPUFInterpolate_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C GPUFInterpolate instance>"
class GPUFInterpolate(GPUFInterpolatePtr):
    """
    Lagrangian interpolation in an GPUFArray
    """
    def __init__(self, name, inArray, xArray, yArray, hwidth, err) :
        self.this = Obit.new_GPUFInterpolate(name, inArray.me, 
                                          xArray.me, yArray.me, hwidth, err.me)
    def __del__(self):
        if Obit!=None:
            Obit.delete_GPUFInterpolate(self.this)
    

def PCreate (name, inArray, xArray, yArray, hwidth, err):
    """
    Create GPU 2D float array interpolator
    
    Returns GPUFInterpolatoe
    * name     = name for object
    * inArray  = 2D pixel array to be interpolated
    * xArray   = 2D array of input x pixel coordinates 
                 corresponding to output pixels
                 Use ImageUtil.PGetXYPixels to create
    * yArray   = 2D array of input x pixel coordinates 
                 corresponding to output pixels
    * hwidth   = half width of interpolation kernel [1,4]
    * err      = Obit Error stack
    """
    out = GPUFInterpolate(name, inArray, xArray, yArray, hwidth, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating GPU interpolator")
    return out
# end PCreate 

def PInterpolateImage (inFI, inArray, outArray, err):
    """
    Interpolate the pixels in inArray to outArray
    
    * inFI     = Python Obit input GPUFInterpolate
    * inArray  = 2D pixel array to be interpolated
    * outArray = 2D array for output, must match pixel arrays given 
                 to constructor
    * err      = Obit Error stack
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit GPUFInterpolate"
    #
    Obit.GPUFInterpolateImage(inFA.me, inArray.me, outArray.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error interpolating image")
   # end PInterpolateImage

def PCopy  (inFI, outFI, err):
    """ 
    Make a deep copy of input object.

    * inFI    = Python Obit input GPUFInterpolate
    * outFI   = Python Obit GPUFInterpolate
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit GPUFInterpolate"
    if not PIsA(outFI):
        raise TypeError,"outFI MUST be a Python Obit GPUFInterpolate"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    Obit.GPUFInterpolateCopy (inFI.me, outFI.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error copying GPUFInterpolate")
    # end PCopy

def PClone (inFI, outFI, err):
    """
    Make a shallow copy of a object (no data copied)

    * inFI    = Python Obit input GPUFInterpolate
    * outFI   = Python Obit GPUFInterpolate
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFI):
        raise TypeError,"inFI MUST be a Python Obit GPUFInterpolate"
    if not PIsA(outFI):
        raise TypeError,"outFI MUST be a Python Obit GPUFInterpolate"
    #
    Obit.GPUFInterpolateClone (inFI.me, outFI.me, err.me)
    # end PClone 

def PIsA (inFI):
    """
    Tells if object thinks it's a Python ObitGPUFInterpolate
    
    return true, false (1,0)
    * inFI    = Python Obit input GPUFInterpolate to test
    """
    ################################################################
    # Checks
    if inFI.__class__ !=  GPUFInterpolate:
        return 0
    return Obit.GPUFInterpolateIsA(inFI.me);
    # end PIsA
