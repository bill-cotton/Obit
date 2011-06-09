# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004,2005
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

# Python shadow class to ObitCArray class
import Obit, FArray

class CArrayPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.CArrayUnref(Obit.CArray_me_get(self.this))
            # In with the new
            Obit.CArray_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != CArray:
            return
        if name == "me" : 
            return Obit.CArray_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != CArray:
            return
        return "<C CArray instance> " + Obit.CArrayGetName(self.me)
class CArray(CArrayPtr):
    """
    Python Obit multidimensional array of complex class
    
    This class is for creating and manipulating a Array as a memory resident 
    multidimensional rectangular array of floats (x2).
    Elements are stored in order of the increasing axis order (the reverse of the
    usual c definition).
    Except as noted, magic value blanking is supported (OBIT_MAGIC).
    """
    def __init__(self, name, naxis=[1]):
        ndim = len(naxis)
        self.this = Obit.new_CArray(name, ndim, naxis)
    def __del__(self):
        if Obit!=None:
            Obit.delete_CArray(self.this)
    
def PGetVal(inCA, pos):
    """ 
    Return value of a cell in an CArray
    
    Returns cell contents. Returns a pair of Reals with the Real and Imaginary
    parts.

    * inCA  = input Python CArray
    * pos   = 0-rel cell number as an array, e.g. [10,24]
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    val = [0.0, 0.0] # create array
    Obit.CArrayGetVal (inCA.me, pos, val)
    return val;
    #  end PGetVal

def PSetVal(inCA, pos, val):
    """
    Set value of a cell in an CArray
    
    * inCA  = input Python CArray
    * pos   = 0-rel cell number as an array, e.g. [10,24]
    * val  = new value for cell as two reals
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArraySetVal(inCA.me, pos, val)
    # end PSetVal

def PCopy (inCA, err):
    """
    Make copy an CArray
    
    returns copy

    * inCA = input Python CArray
    * err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    outCA = CArray("None")
    outCA.me = Obit.CArrayCopy (inCA.me, outFA.me, err.me);
    if err.isErr:
        printErrMsg(err, "Error copying CArray")
    return outCA
    # end PCopy 

def PIsCompatable  (in1, in2):
    """
    Tells if two CArrays have compatable geometry
    
    returns true or false (1, 0)

    * in1 = first input Python CArray
    * in2 = second input Python CArray
    """
    ################################################################
     # Checks
    if not PIsA(in1):
        print "Actually ",in1.__class__
        raise TypeError,"in1 MUST be a Python Obit CArray"
    if not PIsA(in2):
        print "Actually ",in2.__class__
        raise TypeError,"in2 MUST be a Python Obit CArray"
    return Obit.CArrayIsCompatable(in1.me, in2.me)
    # end PIsCompatable


def PRealloc (inCA, naxis):
    """
    Change the geometry of an CArray
    
    Contents will be zeroed

    * inCA  = input Python CArray
    * naxis = array giving desired dimension, e.g. [20,30]
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    ndim = len(naxis)
    Obit.CArrayRealloc(inCA.me, ndim, naxis)
    # end PRealloc 


def PMaxAbs (inCA, pos) :
    """
    Find maximum abs value pixel value
    
    returns  maximum value (may have trouble w/ > 2 dim)

    * inCA  = first input Python CArray
    * pos   = [output] 0-rel position as array
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    lpos = [0,0]                          # Dummy
    ret = Obit.CArrayMax(inCA.me, lpos)   # Results in list ret
    out = ret[0]
    pos = ret[1]
    return out
    # end PMaxAbs

def PNeg (inCA):
    """
    Negate each element of the array.
    
    * inCA = input Python CArray
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArrayNeg(inCA.me)
    # end  PNeg


def PConjg (inCA):
    """
    Conjugate each element of the array.
    
    returns conjugate

    * inCA = input Python CArray
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArrayConjg(inCA.me)
    # end PConjg


def PFill (inCA, cmpx):
    """
    Fill the elements of an array with a complex scalar
    
    * in = cmpx
    * inCA   = input Python CArray
    * cmpx   = Value to fill as  (real,imaginary)
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArrayFill(inCA.me, cmpx)
    # end PFill


def PSAdd (inCA, scalar):
    """
    Add a scalar to each element of the array.

    in = in + scalar

    * inCA   = input Python CArray
    * scalar = Value to add
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArraySAdd(inCA.me, scalar)
    # end PSAdd


def PSMul (inCA, scalar):
    """
    Multiply each element of the array by a scalar
    
    in = in * scalar

    * inCA   = input Python CArray
    * scalar = Value to multiply
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArraySMul(inCA.me, scalar)
    # end PSMul

def PCSAdd (inCA, scalar):
    """
    Add a complex scalar to each element of the array.

    in = in + scalar

    * inCA   = input Python CArray
    * scalar = Value to add float[2]
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArraycSAdd(inCA.me, scalar)
    # end PCSAdd

def PCSMul (inCA, scalar):
    """
    Multiply each element of the array by a complex scalar
    
    in = in * scalar

    * inCA   = input Python CArray
    * scalar = Value to multiply float[2]
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArrayCSMul(inCA.me, scalar)
    # end PCSMul

def PAdd (in1, in2, out):
    """
    Add corresponding elements of two arrays.
    
    out = in1 + in2

    * in1 = first input Python CArray
    * in2 = second input Python CArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.CArrayAdd (in1.me, in2.me, out.me)
    # end  PAdd


def PSub (in1, in2, out):
    """
    Subtract corresponding elements of the arrays.
    
    out = in1 - in2

    * in1 = first input Python CArray
    * in2 = second input Python CArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.CArraySub (in1.me, in2.me, out.me)
    # end PSub

def PMul (in1, in2, out):
    """
    Multiply corresponding elements of the arrays.
    
    out = in1 * in2

    * in1 = first input Python CArray
    * in2 = second input Python CArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.CArrayMul (in1.me, in2.me, out.me)
    # end PMul


def PDiv (in1, in2, out):
    """
    Divide corresponding elements of the arrays.
    
    out = in1 / in2

    * in1 = first input Python CArray
    * in2 = second input Python CArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.CArrayDiv (in1.me, in2.me, out.me)
    # end PDiv


def PMakeF (Cin):
    """
    Make an FArray with the same geometry as a CArray
    
    return FArray size of CArray

    * Cin = input Python CArray
    """
    ################################################################
    # Checks
    if not PIsA(Cin):
        print "Actually ",iCn.__class__
        raise TypeError,"Cin MUST be a Python Obit CArray"
    outFA = FArray.FArray("None")
    outFA.me = Obit.CArrayMakeF(Cin.me)
    return outFA
    # end PMakeF


def PMakeC (Fin):
    """
    Make a CArray with the same geometry as a FArray
    
    return CArray size of Fin

    * Fin = input Python FArray
    """
    ################################################################
    # Checks
    if not FArray.PIsA(Fin):
        print "Actually ",Fin.__class__
        raise TypeError,"Fin MUST be a Python Obit FArray"
    outCA = CArray("None")
    outCA.me = Obit.CArrayMakeC(Fin.me)
    return outCA
    # end PMakeC


def PIsFCompatable  (Cin, Fin):
    """
    Tells if an FArray and a CArray are compatible
    
    returns true or false (1, 0)

    * Cin = input Python CArray
    * Fin = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(Cin):
        print "Actually ",Cin.__class__
        raise TypeError,"Cin MUST be a Python Obit CArray"
    if not FArray.PIsA(Fin):
        print "Actually ",Fin.__class__
        raise TypeError,"Fin MUST be a Python Obit FArray"
    return Obit.CArrayIsFCompatable(Cin.me, Fin.me)
    # end PIsFCompatable


def PFMul (Cin, Fin, out):
    """
    Multiply corresponding elements of a CArray and an FArray.

    * out = Cin * Fin
    * Cin = input Python CArray
    * Fin = input Python FArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (Cin, Fin):
        raise RuntimeError,"Cin and Fin have incompatable geometry"
    if not PIsCompatable  (Cin, out):
        raise RuntimeError,"Cin and out have incompatable geometry"
    Obit.CArrayFMul (Cin.me, Fin.me, out.me)
    # end PFMul

def PFDiv (Cin, Fin, out):
    """
    Divide corresponding elements of a CArray and an FArray.

    
    * out = Cin / Fin
    * Cin = input Python CArray
    * Fin = input Python FArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (Cin, Fin):
        raise RuntimeError,"Cin and Fin have incompatable geometry"
    if not PIsCompatable  (Cin, out):
        raise RuntimeError,"Cin and out have incompatable geometry"
    Obit.CArrayFDiv (Cin.me, Fin.me, out.me)
    # end PFDiv

def PFAdd (Cin, Fin, out):
    """
    Add corresponding elements of a CArray and an FArray.
    
    out = Cin + Fin

    * Cin = input Python CArray
    * Fin = input Python FArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (Cin, Fin):
        raise RuntimeError,"Cin and Fin have incompatable geometry"
    if not PIsCompatable  (Cin, out):
        raise RuntimeError,"Cin and out have incompatable geometry"
    Obit.CArrayFAdd (Cin.me, Fin.me, out.me)
    # end PFAdd


def PFRot (Cin, Fin, out):
    """
    Rotate phase of the elements of a CArray by the elements of an FArray
    
    out = Cin * exp(i*Fin)

    * Cin = input Python CArray
    * Fin = input Python FArray
    * out = output Python CArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (Cin, Fin):
        raise RuntimeError,"Cin and Fin have incompatable geometry"
    if not PIsCompatable  (Cin, out):
        raise RuntimeError,"Cin and out have incompatable geometry"
    Obit.CArrayFRot (Cin.me, Fin.me, out.me)
    # end PFRot


def PComplex (Fin1, Fin2, out):
    """
    Combine two FArrays into a CArray
    
    Fill A complex CArray from two FArrays

    * Fin1  = input Python FArray for Real part
    * Fin2  = input Python FArray for Imaginary part
    * out   = output Python CArray
    """
    ################################################################
    # Checks
    if not FArray.PIsCompatable  (Fin1, Fin2):
        raise RuntimeError,"Fin1 and Fin2 have incompatable geometry"
    if not PIsFCompatable  (out, Fin1):
        raise RuntimeError,"Fin1 and out have incompatable geometry"
    Obit.CArrayComplex (Fin1.me, Fin2.me, out.me)
    # end PComplex 

def PReal (inCA, outFA):
    """
    Copy the Real part of a CArray to an FArray
    
    out = real(in)

    * inCA  = input Python CArray
    * outFA = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (inCA, outFA):
        raise RuntimeError,"Fin1 and out have incompatable geometry"
    Obit.CArrayReal (inCA.me, outFA.me)
    # end PReal

def PImag (inCA, outFA):
    """
    Copy the Imaginary part of a CArray to an FArray
    
    out = imag(in)

    * inCA  = input Python CArray
    * outFA = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (inCA, outFA):
        raise RuntimeError,"Fin1 and out have incompatable geometry"
    Obit.CArrayImag (inCA.me, outFA.me)
    # end PImag

def PAmp (inCA, outFA):
    """
    Copy  amplitude of elements of a CArray to an FArray
    
    outFA = sqrt( real(inCA)^2 + imag(inCA)^2))

    * inCA  = input Python CArray
    * outFA = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (inCA, outFA):
        raise RuntimeError,"Fin1 and out have incompatable geometry"
    Obit.CArrayAmp (inCA.me, outFA.me)
    # end PAmp

def PPhase (inCA, outFA):
    """
    Copy phases of elements of a CArray to an FArray
    
    outFA = atan2(imag(inCA), real(inCA))

    * inCA  = input Python CArray
    * outFA = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsFCompatable  (inCA, outFA):
        raise RuntimeError,"Fin1 and out have incompatable geometry"
    Obit.CArrayPhase (inCA.me, outFA.me)
    # end PPhase


def P2DCenter (inCA):
    """
    Reorder half plane complex to center at the center

    * inCA  = input Python CArray to reorder
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    Obit.CArray2DCenter(inCA.me)
    # end P2DCenter


def PAddConjg (inCA, numConjRow):
    """
    Add conjugate rows to half plane complex image, copy
    from other half plane.
    
    out = in with conjugate rows added.

    * inCA  = input Python CArray
    * numConjRow = number of conjugate rows to copy
    """
    ################################################################
    # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    outCA = CArray("None")
    outCA.me = Obit.CArrayAddConjg (inCA.me, numConjRow)
    return outCA
    # end PAddConjg


def PIsA (inCA):
    """
    Tells if the input is a Python ObitCArray
    
    returns true or false (1,0)

    * inCA = Python Obit CArray to test
    """
    ################################################################
      # Checks
    if inCA.__class__ != CArray:
        return 0
    return Obit.CArrayIsA(inCA.me)
    # end PIsA

def PGetNdim (inCA):
    """
    Returns the number of dimensions in array
    
    * inCA = Python Obit CArray
    """
    ################################################################
      # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    return Obit.CArrayGetNdim(inCA.me)
    # end PGetNdim

def PGetNaxis (inCA):
    """
    Returns array of 7 elements with dimensions in array
    
    * inCA = Python Obit CArray
    """
    ################################################################
      # Checks
    if not PIsA(inCA):
        print "Actually ",inCA.__class__
        raise TypeError,"inCA MUST be a Python Obit CArray"
    return Obit.CArrayGetNaxis(inCA.me)
    # end PGetNaxis
