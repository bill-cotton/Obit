# Python interface to Obit float array class.
# $Id$
""" Python Obit multidimensional array of float class

This class is for creating and manipulating a Array as a memory resident 
multidimensional rectangular array of floats.
Elements are stored in order of the increasing axis order (the reverse of the
usual c definition).
Except as noted, magic value blanking is supported (OBIT_MAGIC) (returned to Python as NaN)

Virtual (read only) members (accessed as e.g. array.RMS
RMS     RMS of valid members (from histogram analysis)
RawRMS  RMS of valid members (from RMS about mean)
Mode    Mode of distribution of valid members 
Mean    Mode of distribution of valid members 
Sum     Sum of valid members 
Count   Count of valid members
Ndim    Number of axes in array
Naxis   list of axis dimensions (by storage order)
"""
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2007
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

# Python shadow class to ObitFArray class
import Obit, OErr

class FArrayPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name, value):
        if name == "me" :
            # Out with the old
            Obit.FArrayUnref(Obit.FArray_me_get(self.this))
            # In with the new
            Obit.FArray_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != FArray:
            return
        if name == "me" : 
            return Obit.FArray_me_get(self.this)
        # Virtual members
        if name=="RMS":
            return PRMS(self)
        if name=="RawRMS":
            return PRawRMS(self)
        if name=="Mode":
            return PMode(self)
        if name=="Mean":
            return PMean(self)
        if name=="Sum":
            return PSum(self)
        if name=="Count":
            return PCount(self)
        if name=="Ndim":
            return PGetNdim(self)
        if name=="Naxis":
            return PGetNaxis(self)
        if name=="Buf":
            return PGetBuf(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != FArray:
            return
        return "<C FArray instance> " + Obit.FArrayGetName(self.me)
class FArray(FArrayPtr):
    """ Python Obit multidimensional array of float class
    
    This class is for creating and manipulating a Array as a memory resident 
    multidimensional rectangular array of floats.
    Elements are stored in order of the increasing axis order (the reverse of the
    usual c definition).
    Except as noted, magic value blanking is supported (OBIT_MAGIC) (returned to Python as NaN)

    Virtual (read only) members (accessed as e.g. array.RMS
        RMS     RMS of valid members (from histogram analysis)
        RawRMS  RMS of valid members (from RMS about mean)
        Mode    Mode of distribution of valid members 
        Mean    Mode of distribution of valid members 
        Sum     Sum of valid members 
        Count   Count of valid members
        Ndim    Number of axes in array
        Naxis   list of axis dimensions (by storage order)
    """
    def __init__(self, name, naxis=[1]):
        ndim = len(naxis)
        self.this = Obit.new_FArray(name, ndim, naxis)
    def __del__(self):
        if Obit!=None:
            Obit.delete_FArray(self.this)
            
    def set (self, value, i1, i2=0, i3=0, i4=0, i5=0, i6=0):
        """ Set Array value [i1, i2, i3...] (0-rel)
        
        self      = Python FArray object
        value     = value for pixel (None = blanked)
        i1        = first axis index (0-rel)
        in        = nth axis index
        """
        # value, possible blanked
        if value==None:
            v = fblank
        else:
            v = value
        # Set value
        pos = [i1,i2,i3,i4,i5,i6]
        PSetVal(self, pos, v)
        # end set

    def get (self, i1, i2=0, i3=0, i4=0, i5=0, i6=0):
        """ Get Array value [i1, i2, i3...] (0-rel)
        
        Return value at pixel [i1,...in], None if blanked
        self      = Python FArray object
        i1        = first axis index (0-rel)
        in        = nth axis index
        """
        # Get value
        pos = [i1,i2,i3,i4,i5,i6]
        v = PGetVal(self, pos)
        # value, possible blanked
        if v==fblank:
            value = None
        else:
            value = v
        return value
        # end get
    # End Class FArray

def PGetBlank():
    """ Return Magic blanking value
    """
    ################################################################
    return Obit.FArrayGetBlank ()


# Module constants
fblank = PGetBlank()  # Blanked value


def PGetVal(inFA, pos):
    """ Return value of a cell in an FArray

    returns cell contents
    inFA  = input Python FArray
    pos   = 0-rel cell number as an array, e.g. [10,24]
    """
    ################################################################
    return Obit.FArrayGetVal (inFA.me, pos)


def PSetVal(inFA, pos, val):
    """  Set value of a cell in an FArray

    inFA  = input Python FArray
    pos   = 0-rel cell number as an array, e.g. [10,24]
    value = new value for cell
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArraySetVal(inFA.me, pos, val)


def PGetBuf(inFA):
    """  Get python memory buffer for data array

    returns python memory buffer
    inFA  = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayGetBuf(inFA.me)
    # end PGetBuf(


def PCopy (inFA, err):
    """  Make copy an FArray

    returns copy
    inFA = input Python FArray
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    outFA = FArray("None")
    outFA.me = Obit.FArrayCopy (inFA.me, outFA.me, err.me);
    if err.isErr:
        OErr.printErrMsg(err, "Error copying FArray")
    return outFA
    # end PCopy

def PClone (inFA, err):
    """  Make copy the structure of an FArray

    Zero fill and return FArray with same structure as in
    inFA  = input Python FArray
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    outFA = FArray("Clone")
    Obit.FArrayClone (inFA.me, outFA.me, err.me);
    if err.isErr:
        OErr.printErrMsg(err, "Error zero cloning FArray")
    return outFA
    # end PClone


def PSubArr  (inFA, blc, trc, err):
    """  Return a slice of an FArray

    returns Slice in FArray
    inFA = input Python FArray
    blc  = array giving (1-rel) lower index of first cell to copy, e.g. [1,1]
    trc  = array giving (1-rel) highest index of first cell to copy
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    outFA = FArray("None")
    outFA.me = Obit.FArraySubArr (inFA.me, blc, trc, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error slicing FArray")
    return outFA
    # emd PSubArr

def PTranspose  (inFA, order, err):
    """  Transpose an FArray

    returns Transposed FArray
    inFA  = input Python FArray
    order = output 1-rel order of the transposed axes, in storage order
            negative value = reverse order, 
            e,g, [2,1] = transpose 2D array
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    outFA = FArray("None")
    outFA.me = Obit.FArrayTranspose (inFA.me, order, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error transposing FArray")
    return outFA
    # end PTranspose

def PIsCompatable  (in1, in2):
    """  Tells if two FArrays have compatable geometry

    returns true or false (1, 0)
    in1 = first input Python FArray
    in2 = second input Python FArray
    """
    ################################################################
     # Checks
    if not PIsA(in1):
        print "Actually ",in1.__class__
        raise TypeError,"in1 MUST be a Python Obit FArray"
    if not PIsA(in2):
        print "Actually ",in2.__class__
        raise TypeError,"in2 MUST be a Python Obit FArray"
    return Obit.FArrayIsCompatable(in1.me, in2.me)


def PRealloc (inFA, naxis):
    """ Change the geometry of an FArray

    Contents will be zeroed
    inFA  = input Python FArray
    naxis = array giving desired dimension, e.g. [20,30]
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    ndim = len(naxis)
    Obit.FArrayRealloc(inFA.me, ndim, naxis)


def PMax (inFA, pos) :
    """ Find maximum pixel value

    returns  maximum value (may have trouble w/ > 2 dim)
    inFA  = first input Python FArray
    pos   = [output] 0-rel position as array
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    lpos = [0,0]                          # Dummy
    ret = Obit.FArrayMax(inFA.me, lpos)   # Results in list ret
    out = ret[0]
    pos[0]=ret[1]; pos[1]=ret[2]
    return out

def PMaxAbs (inFA, pos):
    """ Find maximum absolute pixel value

    returns  maximum abs value (may have trouble w/ > 2 dim)
    inFA  = first input Python FArray
    pos   = [output] 0-rel position as array
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    lpos = [0,0]                             # Dummy
    ret = Obit.FArrayMaxAbs(inFA.me, lpos)   # Results in list ret
    out = ret[0]
    pos[0]=ret[1]; pos[1]=ret[2]
    return out


def PMin (inFA, pos) :
    """ Find minimum pixel value

    returns  minimum value (may have trouble w/ > 2 dim)
    inFA  = first input Python FArray
    pos   = [output] 0-rel position as array
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    lpos = [0,0]                          # Dummy
    ret = Obit.FArrayMin(inFA.me, lpos)   # Results in list ret
    out = ret[0]
    pos[0]=ret[1]; pos[1]=ret[2]
    return out

def PDeblank (inFA, scalar):
    """ Replace any magic value blanks with scalar
    
    inFA   = input Python FArray
    scalar = value to replace magic blanks
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayDeblank (inFA.me, scalar)


def PRMS (inFA):
    """ Return RMS of pixel unblanked pixel values

    returns RMS value derived from a histogram
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayRMS(inFA.me)


def PRawRMS (inFA):
    """ Return RMS of pixel unblanked pixel values

    returns simple RMS about mean
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayRawRMS(inFA.me)


def PRMS0 (inFA):
    """ Return RMS of pixel unblanked pixel values about zero

    returns simple RMS about zero
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayRMS0(inFA.me)


def PMode (inFA):
    """ Return Mode of pixel unblanked pixel values
    
    returns Mode of values
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayMode(inFA.me)

def PMean (inFA):
    """ Return mean of pixel unblanked pixel values

    returns mean of values
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayMean(inFA.me)


def PFill (inFA, scalar):
    """ Fill all cells of an FArray with a scalar

    inFA   = input Python FArray
    scalar = Value to fill
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayFill(inFA.me, scalar)


def PNeg (inFA):
    """ Negate each element of the array.

    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayNeg(inFA.me)
    # end PNeg

def PSin (inFA):
    """ Sine of each element of the array.

    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArraySin(inFA.me)
    # end PSin


def PCos (inFA):
    """ Cosine of each element of the array.

    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayCos(inFA.me)
    # end PCos

def PSqrt (inFA):
    """ Square root of MAX (1.0e-20, each element of the array).

    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArraySqrt(inFA.me)
    # end PSqrt

def PSum (inFA):
    """ Sum each element of the array.

    returns sum
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArraySum(inFA.me)


def PCount (inFA):
    """ Give number of valid elements in the array.

    returns count
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayCount(inFA.me)

def PSAdd (inFA, scalar):
    """ Add a scalar to each element of the array.

    in = in + scalar
    inFA   = input Python FArray
    scalar = Value to add
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArraySAdd(inFA.me, scalar)


def PSMul (inFA, scalar):
    """ Multiply each element of the array by a scalar

    in = in * scalar
    inFA   = input Python FArray
    scalar = Value to multiply
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArraySMul(inFA.me, scalar)


def PSDiv (inFA, scalar):
    """ Divide each element of the array into a scalar.

    in = scalar / in
    No check for zeroes is made 
    inFA   = input Python FArray
    scalar = scalar
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArraySDiv(inFA.me, scalar)


def PClip (inFA, minVal, maxVal, newVal):
    """ Replace values outside of a given range with a new value

    in = newVal where in < minVal or in > maxVal
    inFA   = input Python FArray
    minVal = Minimum allowed value
    maxVal = Maximum allowed value
    newVal = Value to use if out of range
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayClip(inFA.me, minVal, maxVal, newVal)


def PInClip (inFA, minVal, maxVal, newVal):
    """ Replace values inside of a given range with a new value

    in = newVal where in >= minVal or in <= maxVal
    inFA   = input Python FArray
    minVal = Minimum allowed value
    maxVal = Maximum allowed value
    newVal = Value to use if in range
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayInClip(inFA.me, minVal, maxVal, newVal)
    # end PInClip


def PDivClip (inFA1, inFA2, minVal, outFA):
    """ Divide corresponding elements of the arrays with clipping.

    out =  in1 / in2 where in2>minVal, else blanked
    inFA1  = first input Python FArray
    inFA2  = second input Python FArray
    minVal = Minimum allowed value
    outFA  = Maximum allowed value
    """
    ################################################################
    # Checks
    if not PIsA(inFA1):
        print "Actually ",inFA1.__class__
        raise TypeError,"inFA1 MUST be a Python Obit FArray"
    if not PIsA(inFA2):
        print "Actually ",inFA2.__class__
        raise TypeError,"inFA2 MUST be a Python Obit FArray"
    if not PIsA(outFA):
        print "Actually ",outFA.__class__
        raise TypeError,"outFA MUST be a Python Obit FArray"
    Obit.FArrayDivClip(inFA1.me, inFA2.me, minVal, outFA.me)


def PClipBlank (inFA, minVal, maxVal):
    """ Replace values outside of a given range with blank value

    in = blank where in < minVal or in > maxVal
    inFA   = input Python FArray
    minVal = Minimum allowed value
    maxVal = Maximum allowed value
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayClipBlank(inFA.me, minVal, maxVal)


def PBlank (in1, in2, out):
    """  Blank elements of array in1 where array in2 is blanked

    out = in1 or blank where in2 is blank
    in1 = first input Python FArray
    in2 = second input Python FArray with blanking
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArrayBlank (in1.me, in2.me, out.me)


def PSumArr (in1, in2, out):
    """  SSum nonblanked elements of two arrays

    out = (in1 + in2) or whichever is not blanked
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArraySumArr (in1.me, in2.me, out.me)
    # end PSumArr


def PAvgArr (in1, in2, out):
    """  Average nonblanked elements of two arrays.

    out = (in1 + in2)/2 or whichever is not blanked
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArrayAvgArr (in1.me, in2.me, out.me)
    # end PAvgArr

def PMaxArr (in1, in2, out):
    """  Pick the larger nonblanked elements of two arrays.

    out = MAX (in1, in2) or whichever is not blanked
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArrayMaxArr (in1.me, in2.me, out.me)
    # end PMaxArr

def PMinArr (in1, in2, out):
    """  Pick the lesser nonblanked elements of two arrays.

    out = MIN (in1, in2) or whichever is not blanked
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArrayMinArr (in1.me, in2.me, out.me)
    # end PMinArr

def PAdd (in1, in2, out):
    """ Add corresponding elements of two arrays.

    out = in1 + in2
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArrayAdd (in1.me, in2.me, out.me)


def PSub (in1, in2, out):
    """ Subtract corresponding elements of the arrays.

    out = in1 - in2
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArraySub (in1.me, in2.me, out.me)

def PMul (in1, in2, out):
    """ Multiply corresponding elements of the arrays.

    out = in1 * in2
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArrayMul (in1.me, in2.me, out.me)


def PDiv (in1, in2, out):
    """ Divide corresponding elements of the arrays.

    out = in1 / in2
    in1 = first input Python FArray
    in2 = second input Python FArray
    out = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    Obit.FArrayDiv (in1.me, in2.me, out.me)


def PDot (in1, in2):
    """ Sum the products of the elements of two arrays
    
    return Sum (in1 x in2)
    in1 = first input Python FArray
    in2 = second input Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, in2):
        raise RuntimeError,"in1 and in2 have incompatable geometry"
    return Obit.FArrayDot(in1.me, in2.me)


def PMulColRow (inFA, row, col, out):
    """ Multiply elements of 2D array by row times column
    
    Multiply the elements of a 2D array by the corresponding elements
    of a row and column vector.returns cell contents
    out[i,j] = in[i,j] * row[j] * col[i]
    inFA  = input Python 2D FArray
    row   = "row" 1D Python FArray
    col   = "column" 1D Python FArray
    out   = output Python 2D FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayMulColRow (inFA.me, row.me, col.me, out.me)


def PCenter2D (inFA):
    """ Rearrange array for FFT

    In-place rearrangement of a center-at-the edges array to 
    center at the center, or the other way around.
    This is needed for the peculiar order of FFTs.
    FFTs don't like blanked values.
    inFA = input Python FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArray2DCenter (inFA.me)


def PSymInv2D (inFA):
    """ In-place inversion of a symmetric 2-D matrix

    return code, 0=>OK, else could not invert.
    Magic blanking not supported
    inFA   = input Python FArray with symmetric 2-D matrix
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArray2DSymInv (inFA.me)


def PCGauss2D (inFA, Cen, FWHM):
    """ Make 2-D Circular Gaussian

    Peak normalized to 1.0
    in   = Python FArray to be modified
    Cen  = 0-rel pixel center as an array, e,g, [25,26]
    FWMH = FWHM of Gaussian in pixels.
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayCGauss2D (inFA.me, Cen, FWHM)


def PEGauss2D (inFA, amp, Cen, GauMod):
    """ Make 2-D Eliptical Gaussian in FAArray

    Peak normalized to 1.0, model is added to previous contents.
    in     = Python FArray to be modified
    amp    = peak value of Gaussian
    Cen    = 0-rel pixel center as an array, e,g, [25.0,26.0]
    GauMod = Gaussian parameters, Major axis, FWHM, minor axis 
            FWHM (both in pixels) and rotation angle wrt "Y" axis (deg).
            e.g. [3.0,3.0,0.0]
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    Obit.FArrayEGauss2D (inFA.me, amp, Cen, GauMod)


def PShiftAdd (in1, pos1, in2, pos2, scalar, out):
    """ Shift and Add scaled arrays

    Two FArrays are aligned at specified pixels and the corresponding
    pixels are added with a scalar multiplied times the second.
    Only handles to 3 dimensions.
    If in1/out are 3D and in2 is 2D then the same plane in in2 is used
    for all planes in in1/out. 
    out = in1 + scalar * in2 in overlap, else in1
    in1    = first input Python FArray
    pos1   = allignment pixel (0-rel) in in1 as array
    in2    = second input Python FArray
    pos2   = allignment pixel (0-rel) in in2 as array
    scalar = scalar multiplier
    out    = output Python FArray
    """
    ################################################################
    # Checks
    if not PIsCompatable  (in1, out):
        raise RuntimeError,"in1 and out have incompatable geometry"
    if not PIsA(in2):
        print "Actually ",in2.__class__
        raise TypeError,"inn2 MUST be a Python Obit FArray"
    Obit.FArrayShiftAdd (in1.me, pos1, in2.me, pos2, scalar, out.me)
    # end PShiftAdd


def PPad (inFA, outFA, factor):
    """ Zero pad inFA into outFA and multiply by factor, deblank

    outFA is zero filled and the values in inFA are inserted, centered
    and multiplied by factor.  Blanks in the output are replaced by zeroes.

    inFA   = input Python FArray to be centered in outFA
    outFA  = inFA where given and zero elsewhere
    factor = scaling factor for inFA
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    if not PIsA(outFA):
        print "Actually ",outFA.__class__
        raise TypeError,"outFA MUST be a Python Obit FArray"
    Obit.FArrayPad (inFA.me, outFA.me, factor)
    # end 


def PSelInc  (inFA, outFA, blc, trc, inc, err):
    """  Select elements in an FArray by increment

    inFA = input Python FArray
    outFA= output Python FArray
    blc  = array giving (1-rel) lower index of first cell to copy, e.g. [1,1]
    trc  = array giving (1-rel) highest index of first cell to copy
    inc  = increment on each axis
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    if not PIsA(outFA):
        print "Actually ",outFA.__class__
        raise TypeError,"outFA MUST be a Python Obit FArray"
    Obit.FArraySelInc (inFA.me, outFA.me, blc, trc, inc, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error selecting FArray")
    return
    # emd PSelInc

def PIsA (inFA):
    """ Tells if the input is a Python ObitFArray

    returns true or false (1,0)
    inFA = Python Obit FArray to test
    """
    ################################################################
      # Checks
    if inFA.__class__ != FArray:
        return 0
    return Obit.FArrayIsA(inFA.me)
    # end PIsA

def PUnref (inFA):
    """ Decrement reference count

    Decrement reference count which will destroy object if it goes to zero
    Python object stays defined.
    inFA   = Python FArray object
    """
    ################################################################
     # Checks
    if not PIsA(inFA):
        raise TypeError,"inFA MUST be a Python Obit FArray"

    inFA.me = Obit.FArrayUnref(inFA.me)
    # end PUnref

def PGetNdim (inFA):
    """ Returns the number of dimensions in array

    inFA = Python Obit FArray
    """
    ################################################################
      # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayGetNdim(inFA.me)
    # end PGetNdim

def PGetNaxis (inFA):
    """ Returns array of 7 elements with dimensions in array

    inFA = Python Obit FArray
    """
    ################################################################
    # Checks
    if not PIsA(inFA):
        print "Actually ",inFA.__class__
        raise TypeError,"inFA MUST be a Python Obit FArray"
    return Obit.FArrayGetNaxis(inFA.me)
    # end PGetNaxis

