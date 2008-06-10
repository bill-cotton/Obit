"""
FFT Class for use with Obit in python.

Before performing an FFT, an FFT object must first be created.
This is of the form
>>> fft = FFT.FFT(name, dir, type, rank, dim)
where name is an arbitrary name for the object
      dir is the direction 1 => OBIT_FFT_Forward(R2C), 2 => OBIT_FFT_Reverse(C2R)
      type is the FFT type: 1 => OBIT_FFT_FullComplex,  2 => OBIT_FFT_HalfComplex
      rank is the dimensionality of the array being FFTed
      dim is an array giving the dimensionality of each of the rank axes.

Data passed to and from Obit FFT routines are as an FArray or CArray
which are stored in column major (Fortran) order.  For a half plane complex array,
the first axis is half plane axis, dimension n/2+1 where n is the corresponding
dimensionality in the real FFTed array and the first pixel on this axis is the zero
pixel.

Data are passed and returned in "center-at-the-edge" (i.e. unnatural) order
and there is NO transpose of the array axes.  Routine CArray.P2DCenter does
the shuffling for a half plane complex array and FArray.PCenter2D works for a
real array.

Typical usage for a half plane complex to real (nx x ny) FFT is:
>>> # Create FFT object
>>> FFTobj = FFT.FFT("FFT", 2, 2, 2, [nx,ny])
>>> # Shuffle center to corners of half plane complex input
>>> CArray.P2DCenter (aperGrid)
>>> # Do FFT
>>> FFT.PC2R(FFTobj, aperGrid, beamGrid)
>>> # Unshuffle real array beamGrid
>>> FArray.PCenter2D (beamGrid)


"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004,2006
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

# Python shadow class to ObitFFT class
import Obit

class FFTPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.FFTUnref(Obit.FFT_me_get(self.this))
            # In with the new
            Obit.FFT_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != FFT:
            return
        if name == "me" : 
            return Obit.FFT_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != FFT:
            return
        return "<C FFT instance> " + Obit.FFTGetName(self.me)
class FFT(FFTPtr):
    """ Python Obit FFT class
    
    This class is for performing FFT on memory resident data.
    """
    def __init__(self, name, dir, type, rank, dim):
        self.this = Obit.new_FFT(name, dir, type, rank, dim)
    def __del__(self):
        if Obit!=None:
            Obit.delete_FFT(self.this)

def PSuggestSize(length):
    """ Suggest efficient length of FFT equal or larger than length

    returns suggested FFT length (1-D)
    length = length of data array to be transformed.
    """
    ################################################################
    # Checks
    return Obit.FFTSuggestSize (length)
    #  end PSuggestSize


def PR2C(inFFT, inFA, outFA):
    """  Real to half plane complex FFT of multidimensional array

    inFFT  = input Python FFT
    inFA   = input ObitFArray
    outFA  = output ObitCArray
    """
    ################################################################
    # Checks
    if not PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    Obit.FFTR2C(inFFT.me, inFA.me, outFA.me)
    # end PR2C

def PC2R(inFFT, inCA, outFA):
    """  Half plane complex to Real FFT of multidimensional array

    inFFT  = input Python FFT
    inCA   = input ObitCArray
    outFA  = output ObitFArray
    """
    ################################################################
    # Checks
    if not PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    Obit.FFTC2R(inFFT.me, inCA.me, outFA.me)
    # end PC2R

def PC2C(inFFT, inCA, outCA):
    """  Full plane complex to complex FFT of multidimensional array

    inFFT  = input Python FFT
    inCA   = input ObitCArray
    outCA  = output ObitCArray
    """
    ################################################################
    # Checks
    if not PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    Obit.FFTC2C(inFFT.me, inCA.me, outCA.me)
    # end PC2C

def PIsA (inFFT):
    """ Tells if the input is a Python ObitFFT

    returns true or false (1,0)
    inFFT = Python Obit FFT to test
    """
    ################################################################
      # Checks
    if inFFT.__class__ != FFT:
        return 0
    return Obit.FFTIsA(inFFT.me)
    # end PIsA

def PGetName (inFFT):
    """ Get Name of an FFT onject

    returns object name
    inFFT = Python Obit FFT
    """
    ################################################################
      # Checks
    if not PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    return Obit.FFTGetName(inFFT.me)
    # end PGetName

def PGetRank (inFFT):
    """ Get rank of an FFT

    returns object rank
    inFFT = Python Obit FFT
    """
    ################################################################
      # Checks
    if not PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    return Obit.FFTGetRank(inFFT.me)
    # end PGetName

def PGetDim (inFFT):
    """ Get dimension of an FFT

    returns array of 7 elements
    inFFT = Python Obit FFT
    """
    ################################################################
    # Checks
    if not PIsA(inFFT):
        raise TypeError,"inFFT MUST be a Python Obit FFT"
    return Obit.FFTGetDim(inFFT.me)
    # end PGetDim

