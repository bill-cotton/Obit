# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2005
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

# Python utility package for FArrays
import Obit, FArray, OErr

def PFitCGauss(inFA, FWHM, center, peak, err):
    """ Fit Circular Gaussian around peak in FArray

    returns  RMS residual to fit
    inFA    = Array to be fitted
    FWHM    = [in] as list, half width of box around peak to use in fitting
              0 = all.
              [out] as list, Full width Half Max (pixels) of fitted Gaussian
    center  = [out] [x,y] pixel (0-rel) coordinates of peak
    peak    = [out] as list, peak value in fitted Gaussian
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FArray.PIsA(inFA):
        raise TypeError,"inFA MUST be a Python Obit FArray"
    # results in retVal
    retVal = Obit.FArrayUtilFitCGauss(inFA.me, FWHM[0], center, peak[0], err.me)
    FWHM[0] = retVal[0];
    peak[0] = retVal[1];
    center[0] = retVal[2];
    center[1] = retVal[3];
    return retVal[4]
    # end PFitCGauss

def PConvolve(inFA1, inFA2, err):
    """ Return the convolution of two 2-D arrays

    inFA1 and inFA2 must be compatible size and should not contain magic value blanks.
    Convolved array returned.
    inFA1   = First array
    inFA2   = Second array
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FArray.PIsA(inFA1):
        raise TypeError,"inFA1 MUST be a Python Obit FArray"
    if not FArray.PIsA(inFA2):
        raise TypeError,"inFA2 MUST be a Python Obit FArray"
    # 
    outFA = FArray.FArray("None")
    outFA.me = Obit.FArrayUtilConvolve (inFA1.me, inFA2.me, err.me)
    return outFA
    # end PConvolve

