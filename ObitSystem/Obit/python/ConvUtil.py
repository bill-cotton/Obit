# $Id: ConvUtil.py,v 1.1 2006/04/18 20:13:23 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2006
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

# Python utility package for convolving images
import Obit, Image, FArray, OErr

def PConv(inImage, convFn, doDivide, rescale, outImage, err):
    """ (de)Convolve an Image with an FArray and write outImage

    This routine convolves all selected planes in inImage with convFn if  
    doDivide is FALSE, else it does a linear deconvolution 
    Operations are performed using FFTs 
    inImage  = Obit Image to be convolved
    convFn   = Obit/FArray Convolving Function 
    doDovide = If true divide FT of convFn into FT of inImage, else multiply.
    rescale  = Multiplication factor to scale output to correct units
    outImage = Output ObitImage must be a clone of inImage
               Actual convolution size must be set externally 
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not FArray.PIsA(convFn):
        print "Actually ",convFn.__class__
        raise TypeError,"convFn MUST be a Python Obit FArray"
    if not Image.PIsA(outImage):
        print "Actually ",outImage.__class__
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ConvUtilConv (inImage.me, convFn.me, doDivide, rescale, outImage.me, err.me)
    # end PConv

def PGaus(inImage, Beam):
    """ Create an ObitFArray containing a unit area Gaussian in the center

    returns  Python FArray object with normalized Gaussian
    inImage  = Obit Image with geometry
    Beam     = [maj, min, PA] defining eliptical Gaussian
               size in image cell units or pixels if none given
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        print "Actually ",inImage.__class__
        raise TypeError,"inImage MUST be a Python Obit Image"
    if len(Beam) < 3:
        raise TypeError,"Beam MUST have 3 elements"
    #
    outFA = FArray("None")
    outFA.me = Obit.ConvUtilGaus (inImage.me, Beam[0], Beam[1], Beam[2])
    return outFA
    # end PGaus
