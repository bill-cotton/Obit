# This file was created automatically by SWIG 1.3.29.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Obit
import new
new_instancemethod = new.instancemethod
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class swig_globalvar(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, swig_globalvar, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, swig_globalvar, name)
    __repr__ = _swig_repr
    __swig_setmethods__["name"] = _Obit.swig_globalvar_name_set
    __swig_getmethods__["name"] = _Obit.swig_globalvar_name_get
    if _newclass:name = property(_Obit.swig_globalvar_name_get, _Obit.swig_globalvar_name_set)
    __swig_setmethods__["get_attr"] = _Obit.swig_globalvar_get_attr_set
    __swig_getmethods__["get_attr"] = _Obit.swig_globalvar_get_attr_get
    if _newclass:get_attr = property(_Obit.swig_globalvar_get_attr_get, _Obit.swig_globalvar_get_attr_set)
    __swig_setmethods__["set_attr"] = _Obit.swig_globalvar_set_attr_set
    __swig_getmethods__["set_attr"] = _Obit.swig_globalvar_set_attr_get
    if _newclass:set_attr = property(_Obit.swig_globalvar_set_attr_get, _Obit.swig_globalvar_set_attr_set)
    def __init__(self, *args): 
        this = _Obit.new_swig_globalvar(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_swig_globalvar
    __del__ = lambda self : None;
swig_globalvar_swigregister = _Obit.swig_globalvar_swigregister
swig_globalvar_swigregister(swig_globalvar)

class swig_varlinkobject(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, swig_varlinkobject, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, swig_varlinkobject, name)
    __repr__ = _swig_repr
    __swig_setmethods__["maxvars"] = _Obit.swig_varlinkobject_maxvars_set
    __swig_getmethods__["maxvars"] = _Obit.swig_varlinkobject_maxvars_get
    if _newclass:maxvars = property(_Obit.swig_varlinkobject_maxvars_get, _Obit.swig_varlinkobject_maxvars_set)
    def __init__(self, *args): 
        this = _Obit.new_swig_varlinkobject(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_swig_varlinkobject
    __del__ = lambda self : None;
swig_varlinkobject_swigregister = _Obit.swig_varlinkobject_swigregister
swig_varlinkobject_swigregister(swig_varlinkobject)

swig_varlink_repr = _Obit.swig_varlink_repr
swig_varlink_print = _Obit.swig_varlink_print
swig_varlink_getattr = _Obit.swig_varlink_getattr
swig_varlink_setattr = _Obit.swig_varlink_setattr
SWIG_newvarlink = _Obit.SWIG_newvarlink
SWIG_addvarlink = _Obit.SWIG_addvarlink
class SwigPtrType(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPtrType, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPtrType, name)
    __repr__ = _swig_repr
    __swig_setmethods__["name"] = _Obit.SwigPtrType_name_set
    __swig_getmethods__["name"] = _Obit.SwigPtrType_name_get
    if _newclass:name = property(_Obit.SwigPtrType_name_get, _Obit.SwigPtrType_name_set)
    __swig_setmethods__["len"] = _Obit.SwigPtrType_len_set
    __swig_getmethods__["len"] = _Obit.SwigPtrType_len_get
    if _newclass:len = property(_Obit.SwigPtrType_len_get, _Obit.SwigPtrType_len_set)
    __swig_setmethods__["cast"] = _Obit.SwigPtrType_cast_set
    __swig_getmethods__["cast"] = _Obit.SwigPtrType_cast_get
    if _newclass:cast = property(_Obit.SwigPtrType_cast_get, _Obit.SwigPtrType_cast_set)
    __swig_setmethods__["next"] = _Obit.SwigPtrType_next_set
    __swig_getmethods__["next"] = _Obit.SwigPtrType_next_get
    if _newclass:next = property(_Obit.SwigPtrType_next_get, _Obit.SwigPtrType_next_set)
    def __init__(self, *args): 
        this = _Obit.new_SwigPtrType(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_SwigPtrType
    __del__ = lambda self : None;
SwigPtrType_swigregister = _Obit.SwigPtrType_swigregister
SwigPtrType_swigregister(SwigPtrType)

class SwigCacheType(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigCacheType, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigCacheType, name)
    __repr__ = _swig_repr
    __swig_setmethods__["stat"] = _Obit.SwigCacheType_stat_set
    __swig_getmethods__["stat"] = _Obit.SwigCacheType_stat_get
    if _newclass:stat = property(_Obit.SwigCacheType_stat_get, _Obit.SwigCacheType_stat_set)
    __swig_setmethods__["tp"] = _Obit.SwigCacheType_tp_set
    __swig_getmethods__["tp"] = _Obit.SwigCacheType_tp_get
    if _newclass:tp = property(_Obit.SwigCacheType_tp_get, _Obit.SwigCacheType_tp_set)
    __swig_setmethods__["name"] = _Obit.SwigCacheType_name_set
    __swig_getmethods__["name"] = _Obit.SwigCacheType_name_get
    if _newclass:name = property(_Obit.SwigCacheType_name_get, _Obit.SwigCacheType_name_set)
    __swig_setmethods__["mapped"] = _Obit.SwigCacheType_mapped_set
    __swig_getmethods__["mapped"] = _Obit.SwigCacheType_mapped_get
    if _newclass:mapped = property(_Obit.SwigCacheType_mapped_get, _Obit.SwigCacheType_mapped_set)
    def __init__(self, *args): 
        this = _Obit.new_SwigCacheType(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_SwigCacheType
    __del__ = lambda self : None;
SwigCacheType_swigregister = _Obit.SwigCacheType_swigregister
SwigCacheType_swigregister(SwigCacheType)

swigsort = _Obit.swigsort
SWIG_RegisterMapping = _Obit.SWIG_RegisterMapping
SWIG_MakePtr = _Obit.SWIG_MakePtr
SWIG_GetPtr = _Obit.SWIG_GetPtr
SWIG_GetPtrObj = _Obit.SWIG_GetPtrObj
AIPSDirFindCNO = _Obit.AIPSDirFindCNO
AIPSDirHiSeq = _Obit.AIPSDirHiSeq
AIPSDirAlloc = _Obit.AIPSDirAlloc
AIPSDirRemoveEntry = _Obit.AIPSDirRemoveEntry
AIPSDirNumber = _Obit.AIPSDirNumber
AIPSDirInfo = _Obit.AIPSDirInfo
AIPSDirStatus = _Obit.AIPSDirStatus
AIPSSetDirname = _Obit.AIPSSetDirname
newBeamShape = _Obit.newBeamShape
BeamShapeCopy = _Obit.BeamShapeCopy
BeamShapeUnref = _Obit.BeamShapeUnref
BeamShapeRef = _Obit.BeamShapeRef
BeamShapeCreate = _Obit.BeamShapeCreate
BeamShapeGain = _Obit.BeamShapeGain
BeamShapeGainSym = _Obit.BeamShapeGainSym
BeamShapeGetName = _Obit.BeamShapeGetName
BeamShapeIsA = _Obit.BeamShapeIsA
class BeamShape(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, BeamShape, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, BeamShape, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.BeamShape_me_set
    __swig_getmethods__["me"] = _Obit.BeamShape_me_get
    if _newclass:me = property(_Obit.BeamShape_me_get, _Obit.BeamShape_me_set)
    def __init__(self, *args): 
        this = _Obit.new_BeamShape(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_BeamShape
    __del__ = lambda self : None;
BeamShape_swigregister = _Obit.BeamShape_swigregister
BeamShape_swigregister(BeamShape)
cvar = _Obit.cvar

CArrayCreate = _Obit.CArrayCreate
CArrayCopy = _Obit.CArrayCopy
CArrayIsCompatable = _Obit.CArrayIsCompatable
CArrayRealloc = _Obit.CArrayRealloc
CArrayGetVal = _Obit.CArrayGetVal
CArraySetVal = _Obit.CArraySetVal
CArrayMaxAbs = _Obit.CArrayMaxAbs
CArrayNeg = _Obit.CArrayNeg
CArrayConjg = _Obit.CArrayConjg
CArrayFill = _Obit.CArrayFill
CArraySAdd = _Obit.CArraySAdd
CArraySMul = _Obit.CArraySMul
CArrayCSAdd = _Obit.CArrayCSAdd
CArrayCSMul = _Obit.CArrayCSMul
CArrayAdd = _Obit.CArrayAdd
CArraySub = _Obit.CArraySub
CArrayMul = _Obit.CArrayMul
CArrayDiv = _Obit.CArrayDiv
CArrayMakeF = _Obit.CArrayMakeF
CArrayMakeC = _Obit.CArrayMakeC
CArrayIsFCompatable = _Obit.CArrayIsFCompatable
CArrayFMul = _Obit.CArrayFMul
CArrayFDiv = _Obit.CArrayFDiv
CArrayFAdd = _Obit.CArrayFAdd
CArrayFRot = _Obit.CArrayFRot
CArrayComplex = _Obit.CArrayComplex
CArrayReal = _Obit.CArrayReal
CArrayImag = _Obit.CArrayImag
CArrayAmp = _Obit.CArrayAmp
CArrayPhase = _Obit.CArrayPhase
CArray2DCenter = _Obit.CArray2DCenter
CArrayAddConjg = _Obit.CArrayAddConjg
CArrayGetName = _Obit.CArrayGetName
CArrayGetNdim = _Obit.CArrayGetNdim
CArrayGetNaxis = _Obit.CArrayGetNaxis
CArrayIsA = _Obit.CArrayIsA
CArrayRef = _Obit.CArrayRef
CArrayUnref = _Obit.CArrayUnref
class CArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CArray, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.CArray_me_set
    __swig_getmethods__["me"] = _Obit.CArray_me_get
    if _newclass:me = property(_Obit.CArray_me_get, _Obit.CArray_me_set)
    def __init__(self, *args): 
        this = _Obit.new_CArray(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_CArray
    __del__ = lambda self : None;
CArray_swigregister = _Obit.CArray_swigregister
CArray_swigregister(CArray)

TableMF2VL = _Obit.TableMF2VL
TableMFPrint = _Obit.TableMFPrint
TableVLAppend = _Obit.TableVLAppend
TableVLIndex = _Obit.TableVLIndex
TableVLMerge = _Obit.TableVLMerge
TableVLSelect = _Obit.TableVLSelect
TableVLPurge = _Obit.TableVLPurge
TableVLRedun = _Obit.TableVLRedun
TableVLPrint = _Obit.TableVLPrint
TableVL2VZ = _Obit.TableVL2VZ
TableVZSel = _Obit.TableVZSel
newCleanImage = _Obit.newCleanImage
CleanImageCopy = _Obit.CleanImageCopy
CleanImageUnref = _Obit.CleanImageUnref
CleanImageRef = _Obit.CleanImageRef
CleanImageGetList = _Obit.CleanImageGetList
CleanImageGetImageMosaic = _Obit.CleanImageGetImageMosaic
CleanImageSetImageMosaic = _Obit.CleanImageSetImageMosaic
CleanImageGetWindow = _Obit.CleanImageGetWindow
CleanImageSetWindow = _Obit.CleanImageSetWindow
CleanImageAddWindow = _Obit.CleanImageAddWindow
CleanImageCreate = _Obit.CleanImageCreate
CleanImageDeconvolve = _Obit.CleanImageDeconvolve
CleanImageRestore = _Obit.CleanImageRestore
CleanImageDefWindow = _Obit.CleanImageDefWindow
CleanImageGetName = _Obit.CleanImageGetName
CleanImageIsA = _Obit.CleanImageIsA
class CleanImage(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CleanImage, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CleanImage, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.CleanImage_me_set
    __swig_getmethods__["me"] = _Obit.CleanImage_me_get
    if _newclass:me = property(_Obit.CleanImage_me_get, _Obit.CleanImage_me_set)
    def __init__(self, *args): 
        this = _Obit.new_CleanImage(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_CleanImage
    __del__ = lambda self : None;
CleanImage_swigregister = _Obit.CleanImage_swigregister
CleanImage_swigregister(CleanImage)

newCleanVis = _Obit.newCleanVis
CleanVisCopy = _Obit.CleanVisCopy
CleanVisUnref = _Obit.CleanVisUnref
CleanVisRef = _Obit.CleanVisRef
CleanVisGetList = _Obit.CleanVisGetList
CleanVisGetImageMosaic = _Obit.CleanVisGetImageMosaic
CleanVisSetImageMosaic = _Obit.CleanVisSetImageMosaic
CleanVisGetSkyModel = _Obit.CleanVisGetSkyModel
CleanVisSetSkyModel = _Obit.CleanVisSetSkyModel
CleanVisGetWindow = _Obit.CleanVisGetWindow
CleanVisSetWindow = _Obit.CleanVisSetWindow
CleanVisAddWindow = _Obit.CleanVisAddWindow
CleanVisCreate = _Obit.CleanVisCreate
CleanVisDeconvolve = _Obit.CleanVisDeconvolve
CleanVisRestore = _Obit.CleanVisRestore
CleanVisDefWindow = _Obit.CleanVisDefWindow
CleanVisReimage = _Obit.CleanVisReimage
CleanVisGetName = _Obit.CleanVisGetName
CleanVisIsA = _Obit.CleanVisIsA
class CleanVis(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CleanVis, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CleanVis, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.CleanVis_me_set
    __swig_getmethods__["me"] = _Obit.CleanVis_me_get
    if _newclass:me = property(_Obit.CleanVis_me_get, _Obit.CleanVis_me_set)
    def __init__(self, *args): 
        this = _Obit.new_CleanVis(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_CleanVis
    __del__ = lambda self : None;
CleanVis_swigregister = _Obit.CleanVis_swigregister
CleanVis_swigregister(CleanVis)

ConvUtilConv = _Obit.ConvUtilConv
ConvUtilConvGauss = _Obit.ConvUtilConvGauss
ConvUtilGaus = _Obit.ConvUtilGaus
ConvUtilDeconv = _Obit.ConvUtilDeconv
FArrayCreate = _Obit.FArrayCreate
FArrayGetVal = _Obit.FArrayGetVal
FArrayGetBlank = _Obit.FArrayGetBlank
FArraySetVal = _Obit.FArraySetVal
FArrayGetBuf = _Obit.FArrayGetBuf
FArrayCopy = _Obit.FArrayCopy
FArrayClone = _Obit.FArrayClone
FArraySubArr = _Obit.FArraySubArr
FArrayTranspose = _Obit.FArrayTranspose
FArrayIsCompatable = _Obit.FArrayIsCompatable
FArrayRealloc = _Obit.FArrayRealloc
FArrayMax = _Obit.FArrayMax
FArrayMaxAbs = _Obit.FArrayMaxAbs
FArrayMin = _Obit.FArrayMin
FArrayDeblank = _Obit.FArrayDeblank
FArrayRMS = _Obit.FArrayRMS
FArrayRMS0 = _Obit.FArrayRMS0
FArrayRawRMS = _Obit.FArrayRawRMS
FArrayMode = _Obit.FArrayMode
FArrayMean = _Obit.FArrayMean
FArrayFill = _Obit.FArrayFill
FArrayNeg = _Obit.FArrayNeg
FArraySin = _Obit.FArraySin
FArrayCos = _Obit.FArrayCos
FArraySqrt = _Obit.FArraySqrt
FArraySum = _Obit.FArraySum
FArrayCount = _Obit.FArrayCount
FArraySAdd = _Obit.FArraySAdd
FArraySMul = _Obit.FArraySMul
FArraySDiv = _Obit.FArraySDiv
FArrayClip = _Obit.FArrayClip
FArrayInClip = _Obit.FArrayInClip
FArrayDivClip = _Obit.FArrayDivClip
FArrayClipBlank = _Obit.FArrayClipBlank
FArrayBlank = _Obit.FArrayBlank
FArraySumArr = _Obit.FArraySumArr
FArrayAvgArr = _Obit.FArrayAvgArr
FArrayMaxArr = _Obit.FArrayMaxArr
FArrayMinArr = _Obit.FArrayMinArr
FArrayAdd = _Obit.FArrayAdd
FArraySub = _Obit.FArraySub
FArrayMul = _Obit.FArrayMul
FArrayDiv = _Obit.FArrayDiv
FArrayDot = _Obit.FArrayDot
FArrayMulColRow = _Obit.FArrayMulColRow
FArray2DCenter = _Obit.FArray2DCenter
FArraySymInv2D = _Obit.FArraySymInv2D
FArrayCGauss2D = _Obit.FArrayCGauss2D
FArrayEGauss2D = _Obit.FArrayEGauss2D
FArrayShiftAdd = _Obit.FArrayShiftAdd
FArrayPad = _Obit.FArrayPad
FArraySelInc = _Obit.FArraySelInc
FArrayGetName = _Obit.FArrayGetName
FArrayGetNdim = _Obit.FArrayGetNdim
FArrayGetNaxis = _Obit.FArrayGetNaxis
FArrayIsA = _Obit.FArrayIsA
FArrayRef = _Obit.FArrayRef
FArrayUnref = _Obit.FArrayUnref
class FArray(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FArray, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FArray, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.FArray_me_set
    __swig_getmethods__["me"] = _Obit.FArray_me_get
    if _newclass:me = property(_Obit.FArray_me_get, _Obit.FArray_me_set)
    def __init__(self, *args): 
        this = _Obit.new_FArray(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_FArray
    __del__ = lambda self : None;
FArray_swigregister = _Obit.FArray_swigregister
FArray_swigregister(FArray)

FArrayUtilFitCGauss = _Obit.FArrayUtilFitCGauss
FArrayUtilConvolve = _Obit.FArrayUtilConvolve
FArrayUtilUVGaus = _Obit.FArrayUtilUVGaus
FFTCreate = _Obit.FFTCreate
FFTSuggestSize = _Obit.FFTSuggestSize
FFTR2C = _Obit.FFTR2C
FFTC2R = _Obit.FFTC2R
FFTC2C = _Obit.FFTC2C
FFTGetName = _Obit.FFTGetName
FFTGetRank = _Obit.FFTGetRank
FFTGetDim = _Obit.FFTGetDim
FFTIsA = _Obit.FFTIsA
FFTRef = _Obit.FFTRef
FFTUnref = _Obit.FFTUnref
class FFT(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FFT, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FFT, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.FFT_me_set
    __swig_getmethods__["me"] = _Obit.FFT_me_get
    if _newclass:me = property(_Obit.FFT_me_get, _Obit.FFT_me_set)
    def __init__(self, *args): 
        this = _Obit.new_FFT(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_FFT
    __del__ = lambda self : None;
FFT_swigregister = _Obit.FFT_swigregister
FFT_swigregister(FFT)

FInterpolateCreate = _Obit.FInterpolateCreate
FInterpolateCopy = _Obit.FInterpolateCopy
FInterpolateClone = _Obit.FInterpolateClone
FInterpolateReplace = _Obit.FInterpolateReplace
FInterpolatePixel = _Obit.FInterpolatePixel
FInterpolate1D = _Obit.FInterpolate1D
FInterpolatePosition = _Obit.FInterpolatePosition
FInterpolateGetList = _Obit.FInterpolateGetList
FInterpolateGetFArray = _Obit.FInterpolateGetFArray
FInterpolateGetDesc = _Obit.FInterpolateGetDesc
FInterpolateSetDesc = _Obit.FInterpolateSetDesc
FInterpolateGetHwidth = _Obit.FInterpolateGetHwidth
FInterpolateSetHwidth = _Obit.FInterpolateSetHwidth
FInterpolateRef = _Obit.FInterpolateRef
FInterpolateUnref = _Obit.FInterpolateUnref
FInterpolateIsA = _Obit.FInterpolateIsA
class FInterpolate(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FInterpolate, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FInterpolate, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.FInterpolate_me_set
    __swig_getmethods__["me"] = _Obit.FInterpolate_me_get
    if _newclass:me = property(_Obit.FInterpolate_me_get, _Obit.FInterpolate_me_set)
    def __init__(self, *args): 
        this = _Obit.new_FInterpolate(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_FInterpolate
    __del__ = lambda self : None;
FInterpolate_swigregister = _Obit.FInterpolate_swigregister
FInterpolate_swigregister(FInterpolate)

newFitModel = _Obit.newFitModel
FitModelCopy = _Obit.FitModelCopy
FitModelUnref = _Obit.FitModelUnref
FitModelRef = _Obit.FitModelRef
FitModelCreate = _Obit.FitModelCreate
FitModelGetType = _Obit.FitModelGetType
FitModelGetPeak = _Obit.FitModelGetPeak
FitModelGetDeltaX = _Obit.FitModelGetDeltaX
FitModelGetDeltaY = _Obit.FitModelGetDeltaY
FitModelGetNparm = _Obit.FitModelGetNparm
FitModelGetParms = _Obit.FitModelGetParms
FitModelGetePeak = _Obit.FitModelGetePeak
FitModelGeteDeltaX = _Obit.FitModelGeteDeltaX
FitModelGeteDeltaY = _Obit.FitModelGeteDeltaY
FitModelGeteParms = _Obit.FitModelGeteParms
FitModelSetType = _Obit.FitModelSetType
FitModelSetPeak = _Obit.FitModelSetPeak
FitModelSetDeltaX = _Obit.FitModelSetDeltaX
FitModelSetDeltaY = _Obit.FitModelSetDeltaY
FitModelSetNparm = _Obit.FitModelSetNparm
FitModelSetParms = _Obit.FitModelSetParms
FitModelSetePeak = _Obit.FitModelSetePeak
FitModelSeteDeltaX = _Obit.FitModelSeteDeltaX
FitModelSeteDeltaY = _Obit.FitModelSeteDeltaY
FitModelSeteParms = _Obit.FitModelSeteParms
DeconGau = _Obit.DeconGau
FitModelGaussErr = _Obit.FitModelGaussErr
FitModelGetName = _Obit.FitModelGetName
FitModelIsA = _Obit.FitModelIsA
class FitModel(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FitModel, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FitModel, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.FitModel_me_set
    __swig_getmethods__["me"] = _Obit.FitModel_me_get
    if _newclass:me = property(_Obit.FitModel_me_get, _Obit.FitModel_me_set)
    def __init__(self, *args): 
        this = _Obit.new_FitModel(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_FitModel
    __del__ = lambda self : None;
FitModel_swigregister = _Obit.FitModel_swigregister
FitModel_swigregister(FitModel)

newFitRegion = _Obit.newFitRegion
FitRegionCopy = _Obit.FitRegionCopy
FitRegionUnref = _Obit.FitRegionUnref
FitRegionRef = _Obit.FitRegionRef
FitRegionCreate = _Obit.FitRegionCreate
FitRegionSetCorner = _Obit.FitRegionSetCorner
FitRegionSetDim = _Obit.FitRegionSetDim
FitRegionSetPeak = _Obit.FitRegionSetPeak
FitRegionSetPeakResid = _Obit.FitRegionSetPeakResid
FitRegionSetRMSResid = _Obit.FitRegionSetRMSResid
FitRegionSetFluxResid = _Obit.FitRegionSetFluxResid
FitRegionSetNmodel = _Obit.FitRegionSetNmodel
FitRegionSetModels = _Obit.FitRegionSetModels
FitRegionGetCorner = _Obit.FitRegionGetCorner
FitRegionGetDim = _Obit.FitRegionGetDim
FitRegionGetPeak = _Obit.FitRegionGetPeak
FitRegionGetPeakResid = _Obit.FitRegionGetPeakResid
FitRegionGetRMSResid = _Obit.FitRegionGetRMSResid
FitRegionGetFluxResid = _Obit.FitRegionGetFluxResid
FitRegionGetNmodel = _Obit.FitRegionGetNmodel
FitRegionGetModels = _Obit.FitRegionGetModels
FitRegionGetName = _Obit.FitRegionGetName
FitRegionIsA = _Obit.FitRegionIsA
class FitRegion(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FitRegion, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FitRegion, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.FitRegion_me_set
    __swig_getmethods__["me"] = _Obit.FitRegion_me_get
    if _newclass:me = property(_Obit.FitRegion_me_get, _Obit.FitRegion_me_set)
    def __init__(self, *args): 
        this = _Obit.new_FitRegion(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_FitRegion
    __del__ = lambda self : None;
FitRegion_swigregister = _Obit.FitRegion_swigregister
FitRegion_swigregister(FitRegion)

FITSFileExist = _Obit.FITSFileExist
FITSAddDir = _Obit.FITSAddDir
FITSSetDir = _Obit.FITSSetDir
newFullBeam = _Obit.newFullBeam
FullBeamCopy = _Obit.FullBeamCopy
FullBeamUnref = _Obit.FullBeamUnref
FullBeamRef = _Obit.FullBeamRef
FullBeamCreate = _Obit.FullBeamCreate
FullBeamValue = _Obit.FullBeamValue
FullBeamFindPlane = _Obit.FullBeamFindPlane
FullBeamGetName = _Obit.FullBeamGetName
FullBeamIsA = _Obit.FullBeamIsA
class FullBeam(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, FullBeam, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, FullBeam, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.FullBeam_me_set
    __swig_getmethods__["me"] = _Obit.FullBeam_me_get
    if _newclass:me = property(_Obit.FullBeam_me_get, _Obit.FullBeam_me_set)
    def __init__(self, *args): 
        this = _Obit.new_FullBeam(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_FullBeam
    __del__ = lambda self : None;
FullBeam_swigregister = _Obit.FullBeam_swigregister
FullBeam_swigregister(FullBeam)

HistoryCreate = _Obit.HistoryCreate
HistoryZap = _Obit.HistoryZap
HistoryCopy = _Obit.HistoryCopy
HistoryCopyHeader = _Obit.HistoryCopyHeader
HistoryCopy2Header = _Obit.HistoryCopy2Header
HistoryHeader2Header = _Obit.HistoryHeader2Header
HistoryOpen = _Obit.HistoryOpen
HistoryClose = _Obit.HistoryClose
HistoryReadRec = _Obit.HistoryReadRec
HistoryWriteRec = _Obit.HistoryWriteRec
HistoryTimeStamp = _Obit.HistoryTimeStamp
HistoryUnref = _Obit.HistoryUnref
HistoryRef = _Obit.HistoryRef
HistoryGetList = _Obit.HistoryGetList
HistoryIsA = _Obit.HistoryIsA
HistoryGetName = _Obit.HistoryGetName
class History(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, History, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, History, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.History_me_set
    __swig_getmethods__["me"] = _Obit.History_me_get
    if _newclass:me = property(_Obit.History_me_get, _Obit.History_me_set)
    def __init__(self, *args): 
        this = _Obit.new_History(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_History
    __del__ = lambda self : None;
History_swigregister = _Obit.History_swigregister
History_swigregister(History)

ImageDescCreate = _Obit.ImageDescCreate
ImageDescCopy = _Obit.ImageDescCopy
ImageDescCopyDesc = _Obit.ImageDescCopyDesc
ImageDescDefault = _Obit.ImageDescDefault
ImageDescIndex = _Obit.ImageDescIndex
ImageDescOverlap = _Obit.ImageDescOverlap
ImageDescCvtPixel = _Obit.ImageDescCvtPixel
ImageDescGetPos = _Obit.ImageDescGetPos
ImageDescGetPixel = _Obit.ImageDescGetPixel
ImageDescGetList = _Obit.ImageDescGetList
ImageDescGetDict = _Obit.ImageDescGetDict
ImageDescSetDict = _Obit.ImageDescSetDict
ImageDescRef = _Obit.ImageDescRef
ImageDescUnref = _Obit.ImageDescUnref
ImageDescIsA = _Obit.ImageDescIsA
class ImageDesc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ImageDesc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ImageDesc, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.ImageDesc_me_set
    __swig_getmethods__["me"] = _Obit.ImageDesc_me_get
    if _newclass:me = property(_Obit.ImageDesc_me_get, _Obit.ImageDesc_me_set)
    def __init__(self, *args): 
        this = _Obit.new_ImageDesc(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_ImageDesc
    __del__ = lambda self : None;
ImageDesc_swigregister = _Obit.ImageDesc_swigregister
ImageDesc_swigregister(ImageDesc)

newImageFit = _Obit.newImageFit
ImageFitCopy = _Obit.ImageFitCopy
ImageFitUnref = _Obit.ImageFitUnref
ImageFitRef = _Obit.ImageFitRef
ImageFitGetList = _Obit.ImageFitGetList
ImageFitFit = _Obit.ImageFitFit
ImageFitGetName = _Obit.ImageFitGetName
ImageFitIsA = _Obit.ImageFitIsA
class ImageFit(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ImageFit, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ImageFit, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.ImageFit_me_set
    __swig_getmethods__["me"] = _Obit.ImageFit_me_get
    if _newclass:me = property(_Obit.ImageFit_me_get, _Obit.ImageFit_me_set)
    def __init__(self, *args): 
        this = _Obit.new_ImageFit(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_ImageFit
    __del__ = lambda self : None;
ImageFit_swigregister = _Obit.ImageFit_swigregister
ImageFit_swigregister(ImageFit)

ImageSetFITS = _Obit.ImageSetFITS
ImageSetAIPS = _Obit.ImageSetAIPS
ImageCastData = _Obit.ImageCastData
ImageCreate = _Obit.ImageCreate
ImageScratch = _Obit.ImageScratch
ImageInfo = _Obit.ImageInfo
ImageZap = _Obit.ImageZap
ImageRename = _Obit.ImageRename
ImageCopy = _Obit.ImageCopy
ImageClone = _Obit.ImageClone
ImageClone2 = _Obit.ImageClone2
ImageCloneMem = _Obit.ImageCloneMem
ImageGetTable = _Obit.ImageGetTable
ImageOpen = _Obit.ImageOpen
ImageDirty = _Obit.ImageDirty
ImageClose = _Obit.ImageClose
ImageRead = _Obit.ImageRead
ImageWrite = _Obit.ImageWrite
ImageReadFA = _Obit.ImageReadFA
ImageWriteFA = _Obit.ImageWriteFA
ImageGetPlane = _Obit.ImageGetPlane
ImagePutPlane = _Obit.ImagePutPlane
ImageZapTable = _Obit.ImageZapTable
ImageCopyTables = _Obit.ImageCopyTables
ImageUpdateTables = _Obit.ImageUpdateTables
ImagefullInstantiate = _Obit.ImagefullInstantiate
ImageUnref = _Obit.ImageUnref
ImageRef = _Obit.ImageRef
ImageGetList = _Obit.ImageGetList
ImageGetTableList = _Obit.ImageGetTableList
ImageGetDesc = _Obit.ImageGetDesc
ImageSetDesc = _Obit.ImageSetDesc
ImageGetFArray = _Obit.ImageGetFArray
ImageSetFArray = _Obit.ImageSetFArray
ImageGetBeam = _Obit.ImageGetBeam
ImageSetBeam = _Obit.ImageSetBeam
ImageGetHighVer = _Obit.ImageGetHighVer
ImageisScratch = _Obit.ImageisScratch
ImageIsA = _Obit.ImageIsA
ImageGetName = _Obit.ImageGetName
class Image(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Image, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Image, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.Image_me_set
    __swig_getmethods__["me"] = _Obit.Image_me_get
    if _newclass:me = property(_Obit.Image_me_get, _Obit.Image_me_set)
    def __init__(self, *args): 
        this = _Obit.new_Image(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_Image
    __del__ = lambda self : None;
Image_swigregister = _Obit.Image_swigregister
Image_swigregister(Image)

newImageInterp = _Obit.newImageInterp
ImageInterpCopy = _Obit.ImageInterpCopy
ImageInterpUnref = _Obit.ImageInterpUnref
ImageInterpRef = _Obit.ImageInterpRef
ImageInterpCreate = _Obit.ImageInterpCreate
ImageInterpValue = _Obit.ImageInterpValue
ImageInterpFindPlane = _Obit.ImageInterpFindPlane
ImageInterpGetName = _Obit.ImageInterpGetName
ImageInterpIsA = _Obit.ImageInterpIsA
class ImageInterp(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ImageInterp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ImageInterp, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.ImageInterp_me_set
    __swig_getmethods__["me"] = _Obit.ImageInterp_me_get
    if _newclass:me = property(_Obit.ImageInterp_me_get, _Obit.ImageInterp_me_set)
    def __init__(self, *args): 
        this = _Obit.new_ImageInterp(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_ImageInterp
    __del__ = lambda self : None;
ImageInterp_swigregister = _Obit.ImageInterp_swigregister
ImageInterp_swigregister(ImageInterp)

newImageMosaic = _Obit.newImageMosaic
ImageMosaicCopy = _Obit.ImageMosaicCopy
ImageMosaicUnref = _Obit.ImageMosaicUnref
ImageMosaicRef = _Obit.ImageMosaicRef
ImageMosaicZapImage = _Obit.ImageMosaicZapImage
ImageMosaicGetList = _Obit.ImageMosaicGetList
ImageMosaicGetFullImage = _Obit.ImageMosaicGetFullImage
ImageMosaicSetFullImage = _Obit.ImageMosaicSetFullImage
ImageMosaicGetImage = _Obit.ImageMosaicGetImage
ImageMosaicSetImage = _Obit.ImageMosaicSetImage
ImageMosaicCreate = _Obit.ImageMosaicCreate
ImageMosaicDefine = _Obit.ImageMosaicDefine
ImageMosaicFlatten = _Obit.ImageMosaicFlatten
ImageMosaicGetName = _Obit.ImageMosaicGetName
ImageMosaicGetNumber = _Obit.ImageMosaicGetNumber
ImageMosaicIsA = _Obit.ImageMosaicIsA
class ImageMosaic(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ImageMosaic, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ImageMosaic, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.ImageMosaic_me_set
    __swig_getmethods__["me"] = _Obit.ImageMosaic_me_get
    if _newclass:me = property(_Obit.ImageMosaic_me_get, _Obit.ImageMosaic_me_set)
    def __init__(self, *args): 
        this = _Obit.new_ImageMosaic(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_ImageMosaic
    __del__ = lambda self : None;
ImageMosaic_swigregister = _Obit.ImageMosaic_swigregister
ImageMosaic_swigregister(ImageMosaic)

ImageUtilCreateImage = _Obit.ImageUtilCreateImage
ImageUtilMakeImage = _Obit.ImageUtilMakeImage
ImageUtilInterpolateImage = _Obit.ImageUtilInterpolateImage
ImageUtilPBApply = _Obit.ImageUtilPBApply
ImageUtilPBImage = _Obit.ImageUtilPBImage
ImageUtilPBCorr = _Obit.ImageUtilPBCorr
ImageUtilQuanFITS = _Obit.ImageUtilQuanFITS
ImageUtilCCScale = _Obit.ImageUtilCCScale
ImageUtilUVFilter = _Obit.ImageUtilUVFilter
class InfoListBlob(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, InfoListBlob, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, InfoListBlob, name)
    __repr__ = _swig_repr
    __swig_setmethods__["name"] = _Obit.InfoListBlob_name_set
    __swig_getmethods__["name"] = _Obit.InfoListBlob_name_get
    if _newclass:name = property(_Obit.InfoListBlob_name_get, _Obit.InfoListBlob_name_set)
    __swig_setmethods__["type"] = _Obit.InfoListBlob_type_set
    __swig_getmethods__["type"] = _Obit.InfoListBlob_type_get
    if _newclass:type = property(_Obit.InfoListBlob_type_get, _Obit.InfoListBlob_type_set)
    __swig_setmethods__["dim"] = _Obit.InfoListBlob_dim_set
    __swig_getmethods__["dim"] = _Obit.InfoListBlob_dim_get
    if _newclass:dim = property(_Obit.InfoListBlob_dim_get, _Obit.InfoListBlob_dim_set)
    __swig_setmethods__["data"] = _Obit.InfoListBlob_data_set
    __swig_getmethods__["data"] = _Obit.InfoListBlob_data_get
    if _newclass:data = property(_Obit.InfoListBlob_data_get, _Obit.InfoListBlob_data_set)
    def __init__(self, *args): 
        this = _Obit.new_InfoListBlob(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_InfoListBlob
    __del__ = lambda self : None;
InfoListBlob_swigregister = _Obit.InfoListBlob_swigregister
InfoListBlob_swigregister(InfoListBlob)

makeInfoListBlob = _Obit.makeInfoListBlob
freeInfoListBlob = _Obit.freeInfoListBlob
InfoListCreate = _Obit.InfoListCreate
freeInfoList = _Obit.freeInfoList
InfoListCopy = _Obit.InfoListCopy
InfoListRef = _Obit.InfoListRef
InfoListUnref = _Obit.InfoListUnref
InfoListCopyData = _Obit.InfoListCopyData
InfoListRemove = _Obit.InfoListRemove
InfoListItemResize = _Obit.InfoListItemResize
InfoListIsA = _Obit.InfoListIsA
InfoListPutInt = _Obit.InfoListPutInt
InfoListAlwaysPutInt = _Obit.InfoListAlwaysPutInt
InfoListPutLong = _Obit.InfoListPutLong
InfoListAlwaysPutLong = _Obit.InfoListAlwaysPutLong
InfoListPutFloat = _Obit.InfoListPutFloat
InfoListAlwaysPutFloat = _Obit.InfoListAlwaysPutFloat
InfoListPutDouble = _Obit.InfoListPutDouble
InfoListAlwaysPutDouble = _Obit.InfoListAlwaysPutDouble
InfoListPutBoolean = _Obit.InfoListPutBoolean
InfoListAlwaysPutBoolean = _Obit.InfoListAlwaysPutBoolean
InfoListPutString = _Obit.InfoListPutString
InfoListAlwaysPutString = _Obit.InfoListAlwaysPutString
InfoListPutSString = _Obit.InfoListPutSString
InfoListAlwaysPutSString = _Obit.InfoListAlwaysPutSString
InfoListGetHelper = _Obit.InfoListGetHelper
Info2List = _Obit.Info2List
InfoListGet = _Obit.InfoListGet
InfoListGetNumber = _Obit.InfoListGetNumber
InfoListGetDict = _Obit.InfoListGetDict
class InfoList(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, InfoList, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, InfoList, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.InfoList_me_set
    __swig_getmethods__["me"] = _Obit.InfoList_me_get
    if _newclass:me = property(_Obit.InfoList_me_get, _Obit.InfoList_me_set)
    def __init__(self, *args): 
        this = _Obit.new_InfoList(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_InfoList
    __del__ = lambda self : None;
InfoList_swigregister = _Obit.InfoList_swigregister
InfoList_swigregister(InfoList)

IoN2SolNTableConvert = _Obit.IoN2SolNTableConvert
ObitErrUnref = _Obit.ObitErrUnref
ObitErrRef = _Obit.ObitErrRef
ObitErrLog = _Obit.ObitErrLog
ObitErrClear = _Obit.ObitErrClear
ObitErrIsA = _Obit.ObitErrIsA
ObitErrCreate = _Obit.ObitErrCreate
isError = _Obit.isError
SetError = _Obit.SetError
ErrorInit = _Obit.ErrorInit
LogError = _Obit.LogError
OErrMsg = _Obit.OErrMsg
Bomb = _Obit.Bomb
class OErr(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OErr, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OErr, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.OErr_me_set
    __swig_getmethods__["me"] = _Obit.OErr_me_get
    if _newclass:me = property(_Obit.OErr_me_get, _Obit.OErr_me_set)
    def __init__(self, *args): 
        this = _Obit.new_OErr(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_OErr
    __del__ = lambda self : None;
OErr_swigregister = _Obit.OErr_swigregister
OErr_swigregister(OErr)

Startup = _Obit.Startup
Shutdown = _Obit.Shutdown
SystemIsInit = _Obit.SystemIsInit
SystemToday = _Obit.SystemToday
SystemGetPgmName = _Obit.SystemGetPgmName
SystemSetPgmName = _Obit.SystemSetPgmName
SystemGetPgmNumber = _Obit.SystemGetPgmNumber
SystemSetPgmNumber = _Obit.SystemSetPgmNumber
SystemGetAIPSuser = _Obit.SystemGetAIPSuser
SystemSetAIPSuser = _Obit.SystemSetAIPSuser
MemPrint = _Obit.MemPrint
SystemAllowThreads = _Obit.SystemAllowThreads
SystemGetNoThreads = _Obit.SystemGetNoThreads
class OSystem(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OSystem, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OSystem, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.OSystem_me_set
    __swig_getmethods__["me"] = _Obit.OSystem_me_get
    if _newclass:me = property(_Obit.OSystem_me_get, _Obit.OSystem_me_set)
    def __init__(self, *args): 
        this = _Obit.new_OSystem(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_OSystem
    __del__ = lambda self : None;
OSystem_swigregister = _Obit.OSystem_swigregister
OSystem_swigregister(OSystem)

ODataSetFITS = _Obit.ODataSetFITS
ODataSetAIPS = _Obit.ODataSetAIPS
ODataScratch = _Obit.ODataScratch
ODataZap = _Obit.ODataZap
ODataRename = _Obit.ODataRename
ODataCopy = _Obit.ODataCopy
ODataClone = _Obit.ODataClone
newODataTable = _Obit.newODataTable
newODataHistory = _Obit.newODataHistory
ODataZapTable = _Obit.ODataZapTable
ODataCopyTables = _Obit.ODataCopyTables
ODataUpdateTables = _Obit.ODataUpdateTables
ODataFullInstantiate = _Obit.ODataFullInstantiate
ODataOpen = _Obit.ODataOpen
ODataDirty = _Obit.ODataDirty
ODataClose = _Obit.ODataClose
ODataUnref = _Obit.ODataUnref
ODataRef = _Obit.ODataRef
ODataGetList = _Obit.ODataGetList
ODataGetTableList = _Obit.ODataGetTableList
ODataGetHighVer = _Obit.ODataGetHighVer
ODataisScratch = _Obit.ODataisScratch
ODataIsA = _Obit.ODataIsA
ODataGetName = _Obit.ODataGetName
class OData(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OData, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OData, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.OData_me_set
    __swig_getmethods__["me"] = _Obit.OData_me_get
    if _newclass:me = property(_Obit.OData_me_get, _Obit.OData_me_set)
    def __init__(self, *args): 
        this = _Obit.new_OData(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_OData
    __del__ = lambda self : None;
OData_swigregister = _Obit.OData_swigregister
OData_swigregister(OData)

ODisplayCreate = _Obit.ODisplayCreate
ODisplayImage = _Obit.ODisplayImage
ODisplayMosaic = _Obit.ODisplayMosaic
ODisplayImageEdit = _Obit.ODisplayImageEdit
ODisplayMosaicEdit = _Obit.ODisplayMosaicEdit
ODisplayIsA = _Obit.ODisplayIsA
ODisplayRef = _Obit.ODisplayRef
ODisplayUnref = _Obit.ODisplayUnref
ODisplayGetName = _Obit.ODisplayGetName
class ODisplay(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ODisplay, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ODisplay, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.ODisplay_me_set
    __swig_getmethods__["me"] = _Obit.ODisplay_me_get
    if _newclass:me = property(_Obit.ODisplay_me_get, _Obit.ODisplay_me_set)
    def __init__(self, *args): 
        this = _Obit.new_ODisplay(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_ODisplay
    __del__ = lambda self : None;
ODisplay_swigregister = _Obit.ODisplay_swigregister
ODisplay_swigregister(ODisplay)

OPlotCreate = _Obit.OPlotCreate
PlotInitPlot = _Obit.PlotInitPlot
PlotFinishPlot = _Obit.PlotFinishPlot
PlotCopy = _Obit.PlotCopy
PlotXYPlot = _Obit.PlotXYPlot
PlotXYOver = _Obit.PlotXYOver
PlotXYErr = _Obit.PlotXYErr
PlotContour = _Obit.PlotContour
PlotGrayScale = _Obit.PlotGrayScale
PlotMarkCross = _Obit.PlotMarkCross
PlotSetPlot = _Obit.PlotSetPlot
PlotLabel = _Obit.PlotLabel
PlotDrawAxes = _Obit.PlotDrawAxes
PlotSetCharSize = _Obit.PlotSetCharSize
PlotSetLineWidth = _Obit.PlotSetLineWidth
PlotSetLineStyle = _Obit.PlotSetLineStyle
PlotSetColor = _Obit.PlotSetColor
PlotSetPage = _Obit.PlotSetPage
PlotText = _Obit.PlotText
PlotRelText = _Obit.PlotRelText
PlotDrawLine = _Obit.PlotDrawLine
PlotDrawCurve = _Obit.PlotDrawCurve
PlotDrawSymbol = _Obit.PlotDrawSymbol
OPlotIsA = _Obit.OPlotIsA
OPlotRef = _Obit.OPlotRef
OPlotUnref = _Obit.OPlotUnref
OPlotGetName = _Obit.OPlotGetName
PlotGetList = _Obit.PlotGetList
class OPlot(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OPlot, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OPlot, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.OPlot_me_set
    __swig_getmethods__["me"] = _Obit.OPlot_me_get
    if _newclass:me = property(_Obit.OPlot_me_get, _Obit.OPlot_me_set)
    def __init__(self, *args): 
        this = _Obit.new_OPlot(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_OPlot
    __del__ = lambda self : None;
OPlot_swigregister = _Obit.OPlot_swigregister
OPlot_swigregister(OPlot)

OWindowCreate = _Obit.OWindowCreate
OWindowCreate1 = _Obit.OWindowCreate1
OWindowCopy = _Obit.OWindowCopy
OWindowSetList = _Obit.OWindowSetList
OWindowGetList = _Obit.OWindowGetList
OWindowAdd = _Obit.OWindowAdd
OWindowUpdate = _Obit.OWindowUpdate
OWindowDel = _Obit.OWindowDel
OWindowGetMaxID = _Obit.OWindowGetMaxID
OWindowRef = _Obit.OWindowRef
OWindowUnref = _Obit.OWindowUnref
OWindowIsA = _Obit.OWindowIsA
class OWindow(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OWindow, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OWindow, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.OWindow_me_set
    __swig_getmethods__["me"] = _Obit.OWindow_me_get
    if _newclass:me = property(_Obit.OWindow_me_get, _Obit.OWindow_me_set)
    def __init__(self, *args): 
        this = _Obit.new_OWindow(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_OWindow
    __del__ = lambda self : None;
OWindow_swigregister = _Obit.OWindow_swigregister
OWindow_swigregister(OWindow)

Parse = _Obit.Parse
Dump = _Obit.Dump
SkyGeomShiftXY = _Obit.SkyGeomShiftXY
SkyGeomXYShift = _Obit.SkyGeomXYShift
SkyGeomNewPos = _Obit.SkyGeomNewPos
SkyGeomWorldPos = _Obit.SkyGeomWorldPos
SkyGeomCDpos = _Obit.SkyGeomCDpos
SkyGeomXYpix = _Obit.SkyGeomXYpix
SkyGeomCDpix = _Obit.SkyGeomCDpix
SkyGeomWorldPosLM = _Obit.SkyGeomWorldPosLM
SkyGeomXYPixLM = _Obit.SkyGeomXYPixLM
SkyGeomBtoJ = _Obit.SkyGeomBtoJ
SkyGeomJtoB = _Obit.SkyGeomJtoB
SkyGeomEq2Gal = _Obit.SkyGeomEq2Gal
SkyGeomGal2Eq = _Obit.SkyGeomGal2Eq
SkyGeomEq2Ec = _Obit.SkyGeomEq2Ec
SkyGeomEc2Eq = _Obit.SkyGeomEc2Eq
SkyGeomRADec2Zern = _Obit.SkyGeomRADec2Zern
PBUtilPoly = _Obit.PBUtilPoly
PBUtilJinc = _Obit.PBUtilJinc
PBUtilRelPB = _Obit.PBUtilRelPB
newSkyModel = _Obit.newSkyModel
SkyModelCopy = _Obit.SkyModelCopy
SkyModelUnref = _Obit.SkyModelUnref
SkyModelRef = _Obit.SkyModelRef
SkyModelGetList = _Obit.SkyModelGetList
SkyModelGetImageMosaic = _Obit.SkyModelGetImageMosaic
SkyModelSetImageMosaic = _Obit.SkyModelSetImageMosaic
SkyModelCreate = _Obit.SkyModelCreate
SkyModelSubUV = _Obit.SkyModelSubUV
SkyModelDivUV = _Obit.SkyModelDivUV
SkyModelCompressCC = _Obit.SkyModelCompressCC
SkyModelGetName = _Obit.SkyModelGetName
SkyModelIsA = _Obit.SkyModelIsA
class SkyModel(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SkyModel, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SkyModel, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.SkyModel_me_set
    __swig_getmethods__["me"] = _Obit.SkyModel_me_get
    if _newclass:me = property(_Obit.SkyModel_me_get, _Obit.SkyModel_me_set)
    def __init__(self, *args): 
        this = _Obit.new_SkyModel(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_SkyModel
    __del__ = lambda self : None;
SkyModel_swigregister = _Obit.SkyModel_swigregister
SkyModel_swigregister(SkyModel)

newSkyModelVMBeam = _Obit.newSkyModelVMBeam
SkyModelVMBeamCopy = _Obit.SkyModelVMBeamCopy
SkyModelVMBeamUnref = _Obit.SkyModelVMBeamUnref
SkyModelVMBeamRef = _Obit.SkyModelVMBeamRef
SkyModelVMBeamGetList = _Obit.SkyModelVMBeamGetList
SkyModelVMBeamCreate = _Obit.SkyModelVMBeamCreate
SkyModelVMBeamGetName = _Obit.SkyModelVMBeamGetName
SkyModelVMBeamIsA = _Obit.SkyModelVMBeamIsA
class SkyModelVMBeam(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SkyModelVMBeam, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SkyModelVMBeam, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.SkyModelVMBeam_me_set
    __swig_getmethods__["me"] = _Obit.SkyModelVMBeam_me_get
    if _newclass:me = property(_Obit.SkyModelVMBeam_me_get, _Obit.SkyModelVMBeam_me_set)
    def __init__(self, *args): 
        this = _Obit.new_SkyModelVMBeam(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_SkyModelVMBeam
    __del__ = lambda self : None;
SkyModelVMBeam_swigregister = _Obit.SkyModelVMBeam_swigregister
SkyModelVMBeam_swigregister(SkyModelVMBeam)

newSkyModelVMIon = _Obit.newSkyModelVMIon
SkyModelVMIonCopy = _Obit.SkyModelVMIonCopy
SkyModelVMIonUnref = _Obit.SkyModelVMIonUnref
SkyModelVMIonRef = _Obit.SkyModelVMIonRef
SkyModelVMIonGetList = _Obit.SkyModelVMIonGetList
SkyModelVMIonCreate = _Obit.SkyModelVMIonCreate
SkyModelVMIonGetName = _Obit.SkyModelVMIonGetName
SkyModelVMIonIsA = _Obit.SkyModelVMIonIsA
class SkyModelVMIon(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SkyModelVMIon, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SkyModelVMIon, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.SkyModelVMIon_me_set
    __swig_getmethods__["me"] = _Obit.SkyModelVMIon_me_get
    if _newclass:me = property(_Obit.SkyModelVMIon_me_get, _Obit.SkyModelVMIon_me_set)
    def __init__(self, *args): 
        this = _Obit.new_SkyModelVMIon(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_SkyModelVMIon
    __del__ = lambda self : None;
SkyModelVMIon_swigregister = _Obit.SkyModelVMIon_swigregister
SkyModelVMIon_swigregister(SkyModelVMIon)

newSpectrumFit = _Obit.newSpectrumFit
SpectrumFitCopy = _Obit.SpectrumFitCopy
SpectrumFitUnref = _Obit.SpectrumFitUnref
SpectrumFitRef = _Obit.SpectrumFitRef
SpectrumFitCreate = _Obit.SpectrumFitCreate
SpectrumFitCube = _Obit.SpectrumFitCube
SpectrumFitImArr = _Obit.SpectrumFitImArr
SpectrumFitEval = _Obit.SpectrumFitEval
SpectrumFitSingle = _Obit.SpectrumFitSingle
SpectrumFitGetList = _Obit.SpectrumFitGetList
SpectrumFitGetName = _Obit.SpectrumFitGetName
SpectrumFitIsA = _Obit.SpectrumFitIsA
class SpectrumFit(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SpectrumFit, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SpectrumFit, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.SpectrumFit_me_set
    __swig_getmethods__["me"] = _Obit.SpectrumFit_me_get
    if _newclass:me = property(_Obit.SpectrumFit_me_get, _Obit.SpectrumFit_me_set)
    def __init__(self, *args): 
        this = _Obit.new_SpectrumFit(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_SpectrumFit
    __del__ = lambda self : None;
SpectrumFit_swigregister = _Obit.SpectrumFit_swigregister
SpectrumFit_swigregister(SpectrumFit)

TableAN = _Obit.TableAN
TableANGetHeadKeys = _Obit.TableANGetHeadKeys
TableANSetHeadKeys = _Obit.TableANSetHeadKeys
TableAT = _Obit.TableAT
TableATGetHeadKeys = _Obit.TableATGetHeadKeys
TableATSetHeadKeys = _Obit.TableATSetHeadKeys
TableBL = _Obit.TableBL
TableBLGetHeadKeys = _Obit.TableBLGetHeadKeys
TableBLSetHeadKeys = _Obit.TableBLSetHeadKeys
TableBP = _Obit.TableBP
TableBPGetHeadKeys = _Obit.TableBPGetHeadKeys
TableBPSetHeadKeys = _Obit.TableBPSetHeadKeys
TableCC = _Obit.TableCC
TableCCGetHeadKeys = _Obit.TableCCGetHeadKeys
TableCCSetHeadKeys = _Obit.TableCCSetHeadKeys
TableCD = _Obit.TableCD
TableCDGetHeadKeys = _Obit.TableCDGetHeadKeys
TableCDSetHeadKeys = _Obit.TableCDSetHeadKeys
TableCL = _Obit.TableCL
TableCLGetHeadKeys = _Obit.TableCLGetHeadKeys
TableCLSetHeadKeys = _Obit.TableCLSetHeadKeys
TableCQ = _Obit.TableCQ
TableCQGetHeadKeys = _Obit.TableCQGetHeadKeys
TableCQSetHeadKeys = _Obit.TableCQSetHeadKeys
TableCT = _Obit.TableCT
TableCTGetHeadKeys = _Obit.TableCTGetHeadKeys
TableCTSetHeadKeys = _Obit.TableCTSetHeadKeys
TableDescDeeBlank = _Obit.TableDescDeeBlank
TableDescCreate = _Obit.TableDescCreate
TableDescCopy = _Obit.TableDescCopy
TableDescCopyDesc = _Obit.TableDescCopyDesc
TableDescIndex = _Obit.TableDescIndex
TableDescGetList = _Obit.TableDescGetList
TableDescGetDict = _Obit.TableDescGetDict
TableDescSetDict = _Obit.TableDescSetDict
TableDescDef = _Obit.TableDescDef
TableDescRef = _Obit.TableDescRef
TableDescUnref = _Obit.TableDescUnref
TableDescIsA = _Obit.TableDescIsA
class TableDesc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TableDesc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TableDesc, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.TableDesc_me_set
    __swig_getmethods__["me"] = _Obit.TableDesc_me_get
    if _newclass:me = property(_Obit.TableDesc_me_get, _Obit.TableDesc_me_set)
    def __init__(self, *args): 
        this = _Obit.new_TableDesc(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_TableDesc
    __del__ = lambda self : None;
TableDesc_swigregister = _Obit.TableDesc_swigregister
TableDesc_swigregister(TableDesc)

TableFG = _Obit.TableFG
TableFGGetHeadKeys = _Obit.TableFGGetHeadKeys
TableFGSetHeadKeys = _Obit.TableFGSetHeadKeys
TableFQ = _Obit.TableFQ
TableFQGetHeadKeys = _Obit.TableFQGetHeadKeys
TableFQSetHeadKeys = _Obit.TableFQSetHeadKeys
TableGC = _Obit.TableGC
TableGCGetHeadKeys = _Obit.TableGCGetHeadKeys
TableGCSetHeadKeys = _Obit.TableGCSetHeadKeys
TableHistory = _Obit.TableHistory
TableHistoryGetHeadKeys = _Obit.TableHistoryGetHeadKeys
TableHistorySetHeadKeys = _Obit.TableHistorySetHeadKeys
TableIDI_ANTENNA = _Obit.TableIDI_ANTENNA
TableIDI_ANTENNAGetHeadKeys = _Obit.TableIDI_ANTENNAGetHeadKeys
TableIDI_ANTENNASetHeadKeys = _Obit.TableIDI_ANTENNASetHeadKeys
TableIDI_ARRAY_GEOMETRY = _Obit.TableIDI_ARRAY_GEOMETRY
TableIDI_ARRAY_GEOMETRYGetHeadKeys = _Obit.TableIDI_ARRAY_GEOMETRYGetHeadKeys
TableIDI_ARRAY_GEOMETRYSetHeadKeys = _Obit.TableIDI_ARRAY_GEOMETRYSetHeadKeys
TableIDI_BANDPASS = _Obit.TableIDI_BANDPASS
TableIDI_BANDPASSGetHeadKeys = _Obit.TableIDI_BANDPASSGetHeadKeys
TableIDI_BANDPASSSetHeadKeys = _Obit.TableIDI_BANDPASSSetHeadKeys
TableIDI_CALIBRATION = _Obit.TableIDI_CALIBRATION
TableIDI_CALIBRATIONGetHeadKeys = _Obit.TableIDI_CALIBRATIONGetHeadKeys
TableIDI_CALIBRATIONSetHeadKeys = _Obit.TableIDI_CALIBRATIONSetHeadKeys
TableIDI_FLAG = _Obit.TableIDI_FLAG
TableIDI_FLAGGetHeadKeys = _Obit.TableIDI_FLAGGetHeadKeys
TableIDI_FLAGSetHeadKeys = _Obit.TableIDI_FLAGSetHeadKeys
TableIDI_FREQUENCY = _Obit.TableIDI_FREQUENCY
TableIDI_FREQUENCYGetHeadKeys = _Obit.TableIDI_FREQUENCYGetHeadKeys
TableIDI_FREQUENCYSetHeadKeys = _Obit.TableIDI_FREQUENCYSetHeadKeys
TableIDI_GAIN_CURVE = _Obit.TableIDI_GAIN_CURVE
TableIDI_GAIN_CURVEGetHeadKeys = _Obit.TableIDI_GAIN_CURVEGetHeadKeys
TableIDI_GAIN_CURVESetHeadKeys = _Obit.TableIDI_GAIN_CURVESetHeadKeys
TableIDI_INTERFEROMETER_MODEL = _Obit.TableIDI_INTERFEROMETER_MODEL
TableIDI_INTERFEROMETER_MODELGetHeadKeys = _Obit.TableIDI_INTERFEROMETER_MODELGetHeadKeys
TableIDI_INTERFEROMETER_MODELSetHeadKeys = _Obit.TableIDI_INTERFEROMETER_MODELSetHeadKeys
TableIDI_PHASE_CAL = _Obit.TableIDI_PHASE_CAL
TableIDI_PHASE_CALGetHeadKeys = _Obit.TableIDI_PHASE_CALGetHeadKeys
TableIDI_PHASE_CALSetHeadKeys = _Obit.TableIDI_PHASE_CALSetHeadKeys
TableIDI_SOURCE = _Obit.TableIDI_SOURCE
TableIDI_SOURCEGetHeadKeys = _Obit.TableIDI_SOURCEGetHeadKeys
TableIDI_SOURCESetHeadKeys = _Obit.TableIDI_SOURCESetHeadKeys
TableIDI_SYSTEM_TEMPERATURE = _Obit.TableIDI_SYSTEM_TEMPERATURE
TableIDI_SYSTEM_TEMPERATUREGetHeadKeys = _Obit.TableIDI_SYSTEM_TEMPERATUREGetHeadKeys
TableIDI_SYSTEM_TEMPERATURESetHeadKeys = _Obit.TableIDI_SYSTEM_TEMPERATURESetHeadKeys
TableIDI_UV_DATA = _Obit.TableIDI_UV_DATA
TableIDI_UV_DATAGetHeadKeys = _Obit.TableIDI_UV_DATAGetHeadKeys
TableIDI_UV_DATASetHeadKeys = _Obit.TableIDI_UV_DATASetHeadKeys
TableIDI_WEATHER = _Obit.TableIDI_WEATHER
TableIDI_WEATHERGetHeadKeys = _Obit.TableIDI_WEATHERGetHeadKeys
TableIDI_WEATHERSetHeadKeys = _Obit.TableIDI_WEATHERSetHeadKeys
TableIM = _Obit.TableIM
TableIMGetHeadKeys = _Obit.TableIMGetHeadKeys
TableIMSetHeadKeys = _Obit.TableIMSetHeadKeys
TableDeeBlank = _Obit.TableDeeBlank
TableCreate = _Obit.TableCreate
TableZap = _Obit.TableZap
TableCopy = _Obit.TableCopy
TableClone = _Obit.TableClone
TableConcat = _Obit.TableConcat
TablefullInstantiate = _Obit.TablefullInstantiate
TableOpen = _Obit.TableOpen
TableDirty = _Obit.TableDirty
TableClose = _Obit.TableClose
TableReadRow = _Obit.TableReadRow
TableWriteRow = _Obit.TableWriteRow
TableUnref = _Obit.TableUnref
TableRef = _Obit.TableRef
TableGetList = _Obit.TableGetList
TableGetIOList = _Obit.TableGetIOList
TableGetDesc = _Obit.TableGetDesc
TableGetIODesc = _Obit.TableGetIODesc
TableSetDesc = _Obit.TableSetDesc
TableGetVer = _Obit.TableGetVer
TableIsA = _Obit.TableIsA
TableGetName = _Obit.TableGetName
TableUtilSort = _Obit.TableUtilSort
class Table(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Table, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Table, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.Table_me_set
    __swig_getmethods__["me"] = _Obit.Table_me_get
    if _newclass:me = property(_Obit.Table_me_get, _Obit.Table_me_set)
    def __init__(self, *args): 
        this = _Obit.new_Table(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_Table
    __del__ = lambda self : None;
Table_swigregister = _Obit.Table_swigregister
Table_swigregister(Table)

TableListCreate = _Obit.TableListCreate
TableListCopy = _Obit.TableListCopy
TableListGetList = _Obit.TableListGetList
TableListGetHigh = _Obit.TableListGetHigh
TableListPutHi = _Obit.TableListPutHi
TableListRef = _Obit.TableListRef
TableListUnref = _Obit.TableListUnref
TableListIsA = _Obit.TableListIsA
class TableList(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TableList, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TableList, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.TableList_me_set
    __swig_getmethods__["me"] = _Obit.TableList_me_get
    if _newclass:me = property(_Obit.TableList_me_get, _Obit.TableList_me_set)
    def __init__(self, *args): 
        this = _Obit.new_TableList(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_TableList
    __del__ = lambda self : None;
TableList_swigregister = _Obit.TableList_swigregister
TableList_swigregister(TableList)

TableMC = _Obit.TableMC
TableMCGetHeadKeys = _Obit.TableMCGetHeadKeys
TableMCSetHeadKeys = _Obit.TableMCSetHeadKeys
TableMF = _Obit.TableMF
TableMFGetHeadKeys = _Obit.TableMFGetHeadKeys
TableMFSetHeadKeys = _Obit.TableMFSetHeadKeys
TableNI = _Obit.TableNI
TableNIGetHeadKeys = _Obit.TableNIGetHeadKeys
TableNISetHeadKeys = _Obit.TableNISetHeadKeys
TableNX = _Obit.TableNX
TableNXGetHeadKeys = _Obit.TableNXGetHeadKeys
TableNXSetHeadKeys = _Obit.TableNXSetHeadKeys
TableOB = _Obit.TableOB
TableOBGetHeadKeys = _Obit.TableOBGetHeadKeys
TableOBSetHeadKeys = _Obit.TableOBSetHeadKeys
TableOF = _Obit.TableOF
TableOFGetHeadKeys = _Obit.TableOFGetHeadKeys
TableOFSetHeadKeys = _Obit.TableOFSetHeadKeys
TableOT = _Obit.TableOT
TableOTGetHeadKeys = _Obit.TableOTGetHeadKeys
TableOTSetHeadKeys = _Obit.TableOTSetHeadKeys
TablePC = _Obit.TablePC
TablePCGetHeadKeys = _Obit.TablePCGetHeadKeys
TablePCSetHeadKeys = _Obit.TablePCSetHeadKeys
TablePS = _Obit.TablePS
TablePSGetHeadKeys = _Obit.TablePSGetHeadKeys
TablePSSetHeadKeys = _Obit.TablePSSetHeadKeys
TableSN = _Obit.TableSN
TableSNGetHeadKeys = _Obit.TableSNGetHeadKeys
TableSNSetHeadKeys = _Obit.TableSNSetHeadKeys
TableSU = _Obit.TableSU
TableSUGetHeadKeys = _Obit.TableSUGetHeadKeys
TableSUSetHeadKeys = _Obit.TableSUSetHeadKeys
TableSY = _Obit.TableSY
TableSYGetHeadKeys = _Obit.TableSYGetHeadKeys
TableSYSetHeadKeys = _Obit.TableSYSetHeadKeys
TableTY = _Obit.TableTY
TableTYGetHeadKeys = _Obit.TableTYGetHeadKeys
TableTYSetHeadKeys = _Obit.TableTYSetHeadKeys
TableCCUtilMerge = _Obit.TableCCUtilMerge
TableVL = _Obit.TableVL
TableVLGetHeadKeys = _Obit.TableVLGetHeadKeys
TableVLSetHeadKeys = _Obit.TableVLSetHeadKeys
TableVZ = _Obit.TableVZ
TableVZGetHeadKeys = _Obit.TableVZGetHeadKeys
TableVZSetHeadKeys = _Obit.TableVZSetHeadKeys
TableWX = _Obit.TableWX
TableWXGetHeadKeys = _Obit.TableWXGetHeadKeys
TableWXSetHeadKeys = _Obit.TableWXSetHeadKeys
TimeFilterCreate = _Obit.TimeFilterCreate
TimeFilterResize = _Obit.TimeFilterResize
TimeFilterGridTime = _Obit.TimeFilterGridTime
TimeFilterUngridTime = _Obit.TimeFilterUngridTime
TimeFilter2Freq = _Obit.TimeFilter2Freq
TimeFilter2Time = _Obit.TimeFilter2Time
TimeFilterFilter = _Obit.TimeFilterFilter
TimeFilterDoFilter = _Obit.TimeFilterDoFilter
TimeFilterPlotPower = _Obit.TimeFilterPlotPower
TimeFilterPlotTime = _Obit.TimeFilterPlotTime
TimeFilterIsA = _Obit.TimeFilterIsA
TimeFilterRef = _Obit.TimeFilterRef
TimeFilterUnref = _Obit.TimeFilterUnref
TimeFilterGetName = _Obit.TimeFilterGetName
TimeFilterGetTime = _Obit.TimeFilterGetTime
TimeFilterGetFreq = _Obit.TimeFilterGetFreq
TimeFilterGetPower = _Obit.TimeFilterGetPower
TimeFilterSetTime = _Obit.TimeFilterSetTime
TimeFilterSetFreq = _Obit.TimeFilterSetFreq
class TimeFilter(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TimeFilter, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TimeFilter, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.TimeFilter_me_set
    __swig_getmethods__["me"] = _Obit.TimeFilter_me_get
    if _newclass:me = property(_Obit.TimeFilter_me_get, _Obit.TimeFilter_me_set)
    def __init__(self, *args): 
        this = _Obit.new_TimeFilter(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_TimeFilter
    __del__ = lambda self : None;
TimeFilter_swigregister = _Obit.TimeFilter_swigregister
TimeFilter_swigregister(TimeFilter)

UVDescCreate = _Obit.UVDescCreate
UVDescCopy = _Obit.UVDescCopy
UVDescCopyDesc = _Obit.UVDescCopyDesc
UVDescIndex = _Obit.UVDescIndex
UVDescGetList = _Obit.UVDescGetList
UVDescGetDict = _Obit.UVDescGetDict
UVDescSetDict = _Obit.UVDescSetDict
UVDescRef = _Obit.UVDescRef
UVDescUnref = _Obit.UVDescUnref
UVDescIsA = _Obit.UVDescIsA
class UVDesc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UVDesc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UVDesc, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.UVDesc_me_set
    __swig_getmethods__["me"] = _Obit.UVDesc_me_get
    if _newclass:me = property(_Obit.UVDesc_me_get, _Obit.UVDesc_me_set)
    def __init__(self, *args): 
        this = _Obit.new_UVDesc(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_UVDesc
    __del__ = lambda self : None;
UVDesc_swigregister = _Obit.UVDesc_swigregister
UVDesc_swigregister(UVDesc)

newUVGSolve = _Obit.newUVGSolve
UVGSolveCreate = _Obit.UVGSolveCreate
UVGSolveCopy = _Obit.UVGSolveCopy
UVGSolveUnref = _Obit.UVGSolveUnref
UVGSolveRef = _Obit.UVGSolveRef
UVGSolveGetList = _Obit.UVGSolveGetList
UVGSolveIsA = _Obit.UVGSolveIsA
UVGSolveCal = _Obit.UVGSolveCal
class UVGSolve(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UVGSolve, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UVGSolve, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.UVGSolve_me_set
    __swig_getmethods__["me"] = _Obit.UVGSolve_me_get
    if _newclass:me = property(_Obit.UVGSolve_me_get, _Obit.UVGSolve_me_set)
    def __init__(self, *args): 
        this = _Obit.new_UVGSolve(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_UVGSolve
    __del__ = lambda self : None;
UVGSolve_swigregister = _Obit.UVGSolve_swigregister
UVGSolve_swigregister(UVGSolve)

newUVImager = _Obit.newUVImager
UVImagerCopy = _Obit.UVImagerCopy
UVImagerUnref = _Obit.UVImagerUnref
UVImagerRef = _Obit.UVImagerRef
UVImagerGetUV = _Obit.UVImagerGetUV
UVImagerGetMosaic = _Obit.UVImagerGetMosaic
UVImagerCreate = _Obit.UVImagerCreate
UVImagerWeight = _Obit.UVImagerWeight
UVImagerImage = _Obit.UVImagerImage
UVImagerFlatten = _Obit.UVImagerFlatten
UVImagerGetName = _Obit.UVImagerGetName
UVImagerIsA = _Obit.UVImagerIsA
class UVImager(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UVImager, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UVImager, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.UVImager_me_set
    __swig_getmethods__["me"] = _Obit.UVImager_me_get
    if _newclass:me = property(_Obit.UVImager_me_get, _Obit.UVImager_me_set)
    def __init__(self, *args): 
        this = _Obit.new_UVImager(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_UVImager
    __del__ = lambda self : None;
UVImager_swigregister = _Obit.UVImager_swigregister
UVImager_swigregister(UVImager)

UVSetFITS = _Obit.UVSetFITS
UVSetAIPS = _Obit.UVSetAIPS
UVCastData = _Obit.UVCastData
UVCreate = _Obit.UVCreate
UVScratch = _Obit.UVScratch
UVInfo = _Obit.UVInfo
UVZap = _Obit.UVZap
UVRename = _Obit.UVRename
UVCopy = _Obit.UVCopy
UVClone = _Obit.UVClone
newUVTable = _Obit.newUVTable
UVZapTable = _Obit.UVZapTable
UVCopyTables = _Obit.UVCopyTables
UVUpdateTables = _Obit.UVUpdateTables
UVfullInstantiate = _Obit.UVfullInstantiate
UVOpen = _Obit.UVOpen
UVRead = _Obit.UVRead
UVWrite = _Obit.UVWrite
UVRewrite = _Obit.UVRewrite
UVGetVisBuf = _Obit.UVGetVisBuf
UVDirty = _Obit.UVDirty
UVClose = _Obit.UVClose
UVWeightData = _Obit.UVWeightData
UVGetFreq = _Obit.UVGetFreq
UVGetSubA = _Obit.UVGetSubA
UVUnref = _Obit.UVUnref
UVRef = _Obit.UVRef
UVGetList = _Obit.UVGetList
UVGetTableList = _Obit.UVGetTableList
UVGetHighVer = _Obit.UVGetHighVer
UVGetDesc = _Obit.UVGetDesc
UVGetIODesc = _Obit.UVGetIODesc
UVSetDesc = _Obit.UVSetDesc
UVisScratch = _Obit.UVisScratch
UVIsA = _Obit.UVIsA
UVGetName = _Obit.UVGetName
UVUtilUVWExtrema = _Obit.UVUtilUVWExtrema
UVUtilCopyZero = _Obit.UVUtilCopyZero
UVUtilVisDivide = _Obit.UVUtilVisDivide
UVUtilVisSub = _Obit.UVUtilVisSub
UVUtilVisCompare = _Obit.UVUtilVisCompare
UVUtilIndex = _Obit.UVUtilIndex
UVUtilQuack = _Obit.UVUtilQuack
UVUtilHann = _Obit.UVUtilHann
UVUtilAvgF = _Obit.UVUtilAvgF
UVUtilAvgT = _Obit.UVUtilAvgT
UVUtilAvgTF = _Obit.UVUtilAvgTF
UVUtilCount = _Obit.UVUtilCount
UVUtilSplitCh = _Obit.UVUtilSplitCh
UVUtilNoise = _Obit.UVUtilNoise
UVUtilAppend = _Obit.UVUtilAppend
UVEditTD = _Obit.UVEditTD
UVEditFD = _Obit.UVEditFD
UVEditStokes = _Obit.UVEditStokes
UVEditClip = _Obit.UVEditClip
UVEditClipStokes = _Obit.UVEditClipStokes
UVUtilFlag = _Obit.UVUtilFlag
TableCLGetDummy = _Obit.TableCLGetDummy
TableSNGetZeroFR = _Obit.TableSNGetZeroFR
class UV(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UV, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UV, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.UV_me_set
    __swig_getmethods__["me"] = _Obit.UV_me_get
    if _newclass:me = property(_Obit.UV_me_get, _Obit.UV_me_set)
    def __init__(self, *args): 
        this = _Obit.new_UV(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_UV
    __del__ = lambda self : None;
UV_swigregister = _Obit.UV_swigregister
UV_swigregister(UV)

newUVRFIXize = _Obit.newUVRFIXize
UVRFIXizeCreate = _Obit.UVRFIXizeCreate
UVRFIXizeUnref = _Obit.UVRFIXizeUnref
UVRFIXizeRef = _Obit.UVRFIXizeRef
UVRFIXizeCounterRot = _Obit.UVRFIXizeCounterRot
UVRFIXizeFilter = _Obit.UVRFIXizeFilter
UVRFIXizeCorrect = _Obit.UVRFIXizeCorrect
UVRFIXizeGetList = _Obit.UVRFIXizeGetList
UVRFIXizeGetRFI = _Obit.UVRFIXizeGetRFI
UVRFIXizeIsA = _Obit.UVRFIXizeIsA
class UVRFIXize(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UVRFIXize, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UVRFIXize, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.UVRFIXize_me_set
    __swig_getmethods__["me"] = _Obit.UVRFIXize_me_get
    if _newclass:me = property(_Obit.UVRFIXize_me_get, _Obit.UVRFIXize_me_set)
    def __init__(self, *args): 
        this = _Obit.new_UVRFIXize(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_UVRFIXize
    __del__ = lambda self : None;
UVRFIXize_swigregister = _Obit.UVRFIXize_swigregister
UVRFIXize_swigregister(UVRFIXize)

newUVSelfCal = _Obit.newUVSelfCal
UVSelfCalCreate = _Obit.UVSelfCalCreate
UVSelfCalCopy = _Obit.UVSelfCalCopy
UVSelfCalUnref = _Obit.UVSelfCalUnref
UVSelfCalRef = _Obit.UVSelfCalRef
UVSelfCalGetList = _Obit.UVSelfCalGetList
UVSelfCalGetSkyModel = _Obit.UVSelfCalGetSkyModel
UVSelfCalSetSkyModel = _Obit.UVSelfCalSetSkyModel
UVSelfCalIsA = _Obit.UVSelfCalIsA
class UVSelfCal(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UVSelfCal, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UVSelfCal, name)
    __repr__ = _swig_repr
    __swig_setmethods__["me"] = _Obit.UVSelfCal_me_set
    __swig_getmethods__["me"] = _Obit.UVSelfCal_me_get
    if _newclass:me = property(_Obit.UVSelfCal_me_get, _Obit.UVSelfCal_me_set)
    def __init__(self, *args): 
        this = _Obit.new_UVSelfCal(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Obit.delete_UVSelfCal
    __del__ = lambda self : None;
UVSelfCal_swigregister = _Obit.UVSelfCal_swigregister
UVSelfCal_swigregister(UVSelfCal)

UVSoln2Cal = _Obit.UVSoln2Cal
UVSolnRefAnt = _Obit.UVSolnRefAnt
UVSolnSNSmo = _Obit.UVSolnSNSmo
UVSolnDeselSN = _Obit.UVSolnDeselSN
UVSolnDeselCL = _Obit.UVSolnDeselCL
SNInvert = _Obit.SNInvert
UVVisGet = _Obit.UVVisGet
UVVisSet = _Obit.UVVisSet
getVersion = _Obit.getVersion
Zernike = _Obit.Zernike
ZernikeGradX = _Obit.ZernikeGradX
ZernikeGradY = _Obit.ZernikeGradY
ZernikePolar = _Obit.ZernikePolar

version = cvar.version

