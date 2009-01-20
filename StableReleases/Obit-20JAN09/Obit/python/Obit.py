# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _Obit

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types



AIPSDirFindCNO = _Obit.AIPSDirFindCNO

AIPSDirAlloc = _Obit.AIPSDirAlloc

AIPSDirRemoveEntry = _Obit.AIPSDirRemoveEntry

AIPSDirNumber = _Obit.AIPSDirNumber

AIPSDirInfo = _Obit.AIPSDirInfo

AIPSDirInfoClean = _Obit.AIPSDirInfoClean

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
    def __repr__(self):
        return "<C CArray instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.CArray_me_set
    __swig_getmethods__["me"] = _Obit.CArray_me_get
    if _newclass:me = property(_Obit.CArray_me_get, _Obit.CArray_me_set)
    def __init__(self, *args):
        _swig_setattr(self, CArray, 'this', _Obit.new_CArray(*args))
        _swig_setattr(self, CArray, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_CArray):
        try:
            if self.thisown: destroy(self)
        except: pass

class CArrayPtr(CArray):
    def __init__(self, this):
        _swig_setattr(self, CArray, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, CArray, 'thisown', 0)
        _swig_setattr(self, CArray,self.__class__,CArray)
_Obit.CArray_swigregister(CArrayPtr)


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

CleanImageDefWindow = _Obit.CleanImageDefWindow

CleanImageGetName = _Obit.CleanImageGetName

CleanImageIsA = _Obit.CleanImageIsA
class CleanImage(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CleanImage, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CleanImage, name)
    def __repr__(self):
        return "<C CleanImage instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.CleanImage_me_set
    __swig_getmethods__["me"] = _Obit.CleanImage_me_get
    if _newclass:me = property(_Obit.CleanImage_me_get, _Obit.CleanImage_me_set)
    def __init__(self, *args):
        _swig_setattr(self, CleanImage, 'this', _Obit.new_CleanImage(*args))
        _swig_setattr(self, CleanImage, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_CleanImage):
        try:
            if self.thisown: destroy(self)
        except: pass

class CleanImagePtr(CleanImage):
    def __init__(self, this):
        _swig_setattr(self, CleanImage, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, CleanImage, 'thisown', 0)
        _swig_setattr(self, CleanImage,self.__class__,CleanImage)
_Obit.CleanImage_swigregister(CleanImagePtr)


newCleanVis = _Obit.newCleanVis

CleanVisCopy = _Obit.CleanVisCopy

CleanVisUnref = _Obit.CleanVisUnref

CleanVisRef = _Obit.CleanVisRef

CleanVisGetList = _Obit.CleanVisGetList

CleanVisGetImageMosaic = _Obit.CleanVisGetImageMosaic

CleanVisSetImageMosaic = _Obit.CleanVisSetImageMosaic

CleanVisGetWindow = _Obit.CleanVisGetWindow

CleanVisSetWindow = _Obit.CleanVisSetWindow

CleanVisAddWindow = _Obit.CleanVisAddWindow

CleanVisCreate = _Obit.CleanVisCreate

CleanVisDeconvolve = _Obit.CleanVisDeconvolve

CleanVisDefWindow = _Obit.CleanVisDefWindow

CleanVisGetName = _Obit.CleanVisGetName

CleanVisIsA = _Obit.CleanVisIsA
class CleanVis(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CleanVis, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CleanVis, name)
    def __repr__(self):
        return "<C CleanVis instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.CleanVis_me_set
    __swig_getmethods__["me"] = _Obit.CleanVis_me_get
    if _newclass:me = property(_Obit.CleanVis_me_get, _Obit.CleanVis_me_set)
    def __init__(self, *args):
        _swig_setattr(self, CleanVis, 'this', _Obit.new_CleanVis(*args))
        _swig_setattr(self, CleanVis, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_CleanVis):
        try:
            if self.thisown: destroy(self)
        except: pass

class CleanVisPtr(CleanVis):
    def __init__(self, this):
        _swig_setattr(self, CleanVis, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, CleanVis, 'thisown', 0)
        _swig_setattr(self, CleanVis,self.__class__,CleanVis)
_Obit.CleanVis_swigregister(CleanVisPtr)


FArrayCreate = _Obit.FArrayCreate

FArrayGetVal = _Obit.FArrayGetVal

FArrayGetBlank = _Obit.FArrayGetBlank

FArraySetVal = _Obit.FArraySetVal

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

FArrayMode = _Obit.FArrayMode

FArrayMean = _Obit.FArrayMean

FArrayFill = _Obit.FArrayFill

FArrayNeg = _Obit.FArrayNeg

FArraySin = _Obit.FArraySin

FArrayCos = _Obit.FArrayCos

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
    def __repr__(self):
        return "<C FArray instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.FArray_me_set
    __swig_getmethods__["me"] = _Obit.FArray_me_get
    if _newclass:me = property(_Obit.FArray_me_get, _Obit.FArray_me_set)
    def __init__(self, *args):
        _swig_setattr(self, FArray, 'this', _Obit.new_FArray(*args))
        _swig_setattr(self, FArray, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_FArray):
        try:
            if self.thisown: destroy(self)
        except: pass

class FArrayPtr(FArray):
    def __init__(self, this):
        _swig_setattr(self, FArray, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FArray, 'thisown', 0)
        _swig_setattr(self, FArray,self.__class__,FArray)
_Obit.FArray_swigregister(FArrayPtr)


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
    def __repr__(self):
        return "<C FFT instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.FFT_me_set
    __swig_getmethods__["me"] = _Obit.FFT_me_get
    if _newclass:me = property(_Obit.FFT_me_get, _Obit.FFT_me_set)
    def __init__(self, *args):
        _swig_setattr(self, FFT, 'this', _Obit.new_FFT(*args))
        _swig_setattr(self, FFT, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_FFT):
        try:
            if self.thisown: destroy(self)
        except: pass

class FFTPtr(FFT):
    def __init__(self, this):
        _swig_setattr(self, FFT, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FFT, 'thisown', 0)
        _swig_setattr(self, FFT,self.__class__,FFT)
_Obit.FFT_swigregister(FFTPtr)


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
    def __repr__(self):
        return "<C FInterpolate instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.FInterpolate_me_set
    __swig_getmethods__["me"] = _Obit.FInterpolate_me_get
    if _newclass:me = property(_Obit.FInterpolate_me_get, _Obit.FInterpolate_me_set)
    def __init__(self, *args):
        _swig_setattr(self, FInterpolate, 'this', _Obit.new_FInterpolate(*args))
        _swig_setattr(self, FInterpolate, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_FInterpolate):
        try:
            if self.thisown: destroy(self)
        except: pass

class FInterpolatePtr(FInterpolate):
    def __init__(self, this):
        _swig_setattr(self, FInterpolate, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, FInterpolate, 'thisown', 0)
        _swig_setattr(self, FInterpolate,self.__class__,FInterpolate)
_Obit.FInterpolate_swigregister(FInterpolatePtr)


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
    def __repr__(self):
        return "<C History instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.History_me_set
    __swig_getmethods__["me"] = _Obit.History_me_get
    if _newclass:me = property(_Obit.History_me_get, _Obit.History_me_set)
    def __init__(self, *args):
        _swig_setattr(self, History, 'this', _Obit.new_History(*args))
        _swig_setattr(self, History, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_History):
        try:
            if self.thisown: destroy(self)
        except: pass

class HistoryPtr(History):
    def __init__(self, this):
        _swig_setattr(self, History, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, History, 'thisown', 0)
        _swig_setattr(self, History,self.__class__,History)
_Obit.History_swigregister(HistoryPtr)


ImageDescCreate = _Obit.ImageDescCreate

ImageDescCopy = _Obit.ImageDescCopy

ImageDescCopyDesc = _Obit.ImageDescCopyDesc

ImageDescDefault = _Obit.ImageDescDefault

ImageDescIndex = _Obit.ImageDescIndex

ImageDescCvtPixel = _Obit.ImageDescCvtPixel

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
    def __repr__(self):
        return "<C ImageDesc instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.ImageDesc_me_set
    __swig_getmethods__["me"] = _Obit.ImageDesc_me_get
    if _newclass:me = property(_Obit.ImageDesc_me_get, _Obit.ImageDesc_me_set)
    def __init__(self, *args):
        _swig_setattr(self, ImageDesc, 'this', _Obit.new_ImageDesc(*args))
        _swig_setattr(self, ImageDesc, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_ImageDesc):
        try:
            if self.thisown: destroy(self)
        except: pass

class ImageDescPtr(ImageDesc):
    def __init__(self, this):
        _swig_setattr(self, ImageDesc, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ImageDesc, 'thisown', 0)
        _swig_setattr(self, ImageDesc,self.__class__,ImageDesc)
_Obit.ImageDesc_swigregister(ImageDescPtr)


ImageSetFITS = _Obit.ImageSetFITS

ImageSetAIPS = _Obit.ImageSetAIPS

ImageCastData = _Obit.ImageCastData

ImageCreate = _Obit.ImageCreate

ImageScratch = _Obit.ImageScratch

ImageZap = _Obit.ImageZap

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
    def __repr__(self):
        return "<C Image instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.Image_me_set
    __swig_getmethods__["me"] = _Obit.Image_me_get
    if _newclass:me = property(_Obit.Image_me_get, _Obit.Image_me_set)
    def __init__(self, *args):
        _swig_setattr(self, Image, 'this', _Obit.new_Image(*args))
        _swig_setattr(self, Image, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_Image):
        try:
            if self.thisown: destroy(self)
        except: pass

class ImagePtr(Image):
    def __init__(self, this):
        _swig_setattr(self, Image, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Image, 'thisown', 0)
        _swig_setattr(self, Image,self.__class__,Image)
_Obit.Image_swigregister(ImagePtr)


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
    def __repr__(self):
        return "<C ImageMosaic instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.ImageMosaic_me_set
    __swig_getmethods__["me"] = _Obit.ImageMosaic_me_get
    if _newclass:me = property(_Obit.ImageMosaic_me_get, _Obit.ImageMosaic_me_set)
    def __init__(self, *args):
        _swig_setattr(self, ImageMosaic, 'this', _Obit.new_ImageMosaic(*args))
        _swig_setattr(self, ImageMosaic, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_ImageMosaic):
        try:
            if self.thisown: destroy(self)
        except: pass

class ImageMosaicPtr(ImageMosaic):
    def __init__(self, this):
        _swig_setattr(self, ImageMosaic, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ImageMosaic, 'thisown', 0)
        _swig_setattr(self, ImageMosaic,self.__class__,ImageMosaic)
_Obit.ImageMosaic_swigregister(ImageMosaicPtr)


ImageUtilCreateImage = _Obit.ImageUtilCreateImage

ImageUtilMakeImage = _Obit.ImageUtilMakeImage

ImageUtilInterpolateImage = _Obit.ImageUtilInterpolateImage

ImageUtilPBApply = _Obit.ImageUtilPBApply

ImageUtilPBImage = _Obit.ImageUtilPBImage

ImageUtilPBCorr = _Obit.ImageUtilPBCorr
class InfoListBlob(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, InfoListBlob, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, InfoListBlob, name)
    def __repr__(self):
        return "<C InfoListBlob instance at %s>" % (self.this,)
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
        _swig_setattr(self, InfoListBlob, 'this', _Obit.new_InfoListBlob(*args))
        _swig_setattr(self, InfoListBlob, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_InfoListBlob):
        try:
            if self.thisown: destroy(self)
        except: pass

class InfoListBlobPtr(InfoListBlob):
    def __init__(self, this):
        _swig_setattr(self, InfoListBlob, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, InfoListBlob, 'thisown', 0)
        _swig_setattr(self, InfoListBlob,self.__class__,InfoListBlob)
_Obit.InfoListBlob_swigregister(InfoListBlobPtr)


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

InfoListGet = _Obit.InfoListGet
class InfoList(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, InfoList, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, InfoList, name)
    def __repr__(self):
        return "<C InfoList instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.InfoList_me_set
    __swig_getmethods__["me"] = _Obit.InfoList_me_get
    if _newclass:me = property(_Obit.InfoList_me_get, _Obit.InfoList_me_set)
    def __init__(self, *args):
        _swig_setattr(self, InfoList, 'this', _Obit.new_InfoList(*args))
        _swig_setattr(self, InfoList, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_InfoList):
        try:
            if self.thisown: destroy(self)
        except: pass

class InfoListPtr(InfoList):
    def __init__(self, this):
        _swig_setattr(self, InfoList, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, InfoList, 'thisown', 0)
        _swig_setattr(self, InfoList,self.__class__,InfoList)
_Obit.InfoList_swigregister(InfoListPtr)


ObitErrUnref = _Obit.ObitErrUnref

ObitErrRef = _Obit.ObitErrRef

ObitErrLog = _Obit.ObitErrLog

ObitErrClear = _Obit.ObitErrClear

ObitErrIsA = _Obit.ObitErrIsA

ObitErrCreate = _Obit.ObitErrCreate

isError = _Obit.isError

OErrMsg = _Obit.OErrMsg

Bomb = _Obit.Bomb
class OErr(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OErr, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OErr, name)
    def __repr__(self):
        return "<C OErr instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.OErr_me_set
    __swig_getmethods__["me"] = _Obit.OErr_me_get
    if _newclass:me = property(_Obit.OErr_me_get, _Obit.OErr_me_set)
    def __init__(self, *args):
        _swig_setattr(self, OErr, 'this', _Obit.new_OErr(*args))
        _swig_setattr(self, OErr, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_OErr):
        try:
            if self.thisown: destroy(self)
        except: pass

class OErrPtr(OErr):
    def __init__(self, this):
        _swig_setattr(self, OErr, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OErr, 'thisown', 0)
        _swig_setattr(self, OErr,self.__class__,OErr)
_Obit.OErr_swigregister(OErrPtr)


Startup = _Obit.Startup

Shutdown = _Obit.Shutdown

SystemGetPgmName = _Obit.SystemGetPgmName

SystemSetPgmName = _Obit.SystemSetPgmName

SystemGetPgmNumber = _Obit.SystemGetPgmNumber

SystemSetPgmNumber = _Obit.SystemSetPgmNumber

SystemGetAIPSuser = _Obit.SystemGetAIPSuser

SystemSetAIPSuser = _Obit.SystemSetAIPSuser

MemPrint = _Obit.MemPrint
class OSystem(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OSystem, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OSystem, name)
    def __repr__(self):
        return "<C OSystem instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.OSystem_me_set
    __swig_getmethods__["me"] = _Obit.OSystem_me_get
    if _newclass:me = property(_Obit.OSystem_me_get, _Obit.OSystem_me_set)
    def __init__(self, *args):
        _swig_setattr(self, OSystem, 'this', _Obit.new_OSystem(*args))
        _swig_setattr(self, OSystem, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_OSystem):
        try:
            if self.thisown: destroy(self)
        except: pass

class OSystemPtr(OSystem):
    def __init__(self, this):
        _swig_setattr(self, OSystem, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OSystem, 'thisown', 0)
        _swig_setattr(self, OSystem,self.__class__,OSystem)
_Obit.OSystem_swigregister(OSystemPtr)


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
    def __repr__(self):
        return "<C ODisplay instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.ODisplay_me_set
    __swig_getmethods__["me"] = _Obit.ODisplay_me_get
    if _newclass:me = property(_Obit.ODisplay_me_get, _Obit.ODisplay_me_set)
    def __init__(self, *args):
        _swig_setattr(self, ODisplay, 'this', _Obit.new_ODisplay(*args))
        _swig_setattr(self, ODisplay, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_ODisplay):
        try:
            if self.thisown: destroy(self)
        except: pass

class ODisplayPtr(ODisplay):
    def __init__(self, this):
        _swig_setattr(self, ODisplay, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ODisplay, 'thisown', 0)
        _swig_setattr(self, ODisplay,self.__class__,ODisplay)
_Obit.ODisplay_swigregister(ODisplayPtr)


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
    def __repr__(self):
        return "<C OWindow instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.OWindow_me_set
    __swig_getmethods__["me"] = _Obit.OWindow_me_get
    if _newclass:me = property(_Obit.OWindow_me_get, _Obit.OWindow_me_set)
    def __init__(self, *args):
        _swig_setattr(self, OWindow, 'this', _Obit.new_OWindow(*args))
        _swig_setattr(self, OWindow, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_OWindow):
        try:
            if self.thisown: destroy(self)
        except: pass

class OWindowPtr(OWindow):
    def __init__(self, this):
        _swig_setattr(self, OWindow, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, OWindow, 'thisown', 0)
        _swig_setattr(self, OWindow,self.__class__,OWindow)
_Obit.OWindow_swigregister(OWindowPtr)


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

SkyModelGetName = _Obit.SkyModelGetName

SkyModelIsA = _Obit.SkyModelIsA
class SkyModel(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SkyModel, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SkyModel, name)
    def __repr__(self):
        return "<C SkyModel instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.SkyModel_me_set
    __swig_getmethods__["me"] = _Obit.SkyModel_me_get
    if _newclass:me = property(_Obit.SkyModel_me_get, _Obit.SkyModel_me_set)
    def __init__(self, *args):
        _swig_setattr(self, SkyModel, 'this', _Obit.new_SkyModel(*args))
        _swig_setattr(self, SkyModel, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_SkyModel):
        try:
            if self.thisown: destroy(self)
        except: pass

class SkyModelPtr(SkyModel):
    def __init__(self, this):
        _swig_setattr(self, SkyModel, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, SkyModel, 'thisown', 0)
        _swig_setattr(self, SkyModel,self.__class__,SkyModel)
_Obit.SkyModel_swigregister(SkyModelPtr)


TableAN = _Obit.TableAN

TableBL = _Obit.TableBL

TableBP = _Obit.TableBP

TableCC = _Obit.TableCC

TableCL = _Obit.TableCL

TableCQ = _Obit.TableCQ

TableDescDeeBlank = _Obit.TableDescDeeBlank

TableDescCreate = _Obit.TableDescCreate

TableDescCopy = _Obit.TableDescCopy

TableDescCopyDesc = _Obit.TableDescCopyDesc

TableDescIndex = _Obit.TableDescIndex

TableDescGetList = _Obit.TableDescGetList

TableDescGetDict = _Obit.TableDescGetDict

TableDescSetDict = _Obit.TableDescSetDict

TableDescRef = _Obit.TableDescRef

TableDescUnref = _Obit.TableDescUnref

TableDescIsA = _Obit.TableDescIsA
class TableDesc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TableDesc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TableDesc, name)
    def __repr__(self):
        return "<C TableDesc instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.TableDesc_me_set
    __swig_getmethods__["me"] = _Obit.TableDesc_me_get
    if _newclass:me = property(_Obit.TableDesc_me_get, _Obit.TableDesc_me_set)
    def __init__(self, *args):
        _swig_setattr(self, TableDesc, 'this', _Obit.new_TableDesc(*args))
        _swig_setattr(self, TableDesc, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_TableDesc):
        try:
            if self.thisown: destroy(self)
        except: pass

class TableDescPtr(TableDesc):
    def __init__(self, this):
        _swig_setattr(self, TableDesc, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, TableDesc, 'thisown', 0)
        _swig_setattr(self, TableDesc,self.__class__,TableDesc)
_Obit.TableDesc_swigregister(TableDescPtr)


TableFG = _Obit.TableFG

TableFQ = _Obit.TableFQ

TableHistory = _Obit.TableHistory

TableDeeBlank = _Obit.TableDeeBlank

TableCreate = _Obit.TableCreate

TableZap = _Obit.TableZap

TableCopy = _Obit.TableCopy

TableClone = _Obit.TableClone

TableConcat = _Obit.TableConcat

TablefullInstantiate = _Obit.TablefullInstantiate

TableOpen = _Obit.TableOpen

TableClose = _Obit.TableClose

TableReadRow = _Obit.TableReadRow

TableWriteRow = _Obit.TableWriteRow

TableUnref = _Obit.TableUnref

TableRef = _Obit.TableRef

TableGetList = _Obit.TableGetList

TableGetIOList = _Obit.TableGetIOList

TableGetDesc = _Obit.TableGetDesc

TableGetIODesc = _Obit.TableGetIODesc

TableGetVer = _Obit.TableGetVer

TableIsA = _Obit.TableIsA

TableGetName = _Obit.TableGetName

TableUtilSort = _Obit.TableUtilSort
class Table(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Table, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Table, name)
    def __repr__(self):
        return "<C Table instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.Table_me_set
    __swig_getmethods__["me"] = _Obit.Table_me_get
    if _newclass:me = property(_Obit.Table_me_get, _Obit.Table_me_set)
    def __init__(self, *args):
        _swig_setattr(self, Table, 'this', _Obit.new_Table(*args))
        _swig_setattr(self, Table, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_Table):
        try:
            if self.thisown: destroy(self)
        except: pass

class TablePtr(Table):
    def __init__(self, this):
        _swig_setattr(self, Table, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Table, 'thisown', 0)
        _swig_setattr(self, Table,self.__class__,Table)
_Obit.Table_swigregister(TablePtr)


TableListCreate = _Obit.TableListCreate

TableListCopy = _Obit.TableListCopy

TableListGetList = _Obit.TableListGetList

TableListGetHigh = _Obit.TableListGetHigh

TableListRef = _Obit.TableListRef

TableListUnref = _Obit.TableListUnref

TableListIsA = _Obit.TableListIsA
class TableList(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, TableList, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, TableList, name)
    def __repr__(self):
        return "<C TableList instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.TableList_me_set
    __swig_getmethods__["me"] = _Obit.TableList_me_get
    if _newclass:me = property(_Obit.TableList_me_get, _Obit.TableList_me_set)
    def __init__(self, *args):
        _swig_setattr(self, TableList, 'this', _Obit.new_TableList(*args))
        _swig_setattr(self, TableList, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_TableList):
        try:
            if self.thisown: destroy(self)
        except: pass

class TableListPtr(TableList):
    def __init__(self, this):
        _swig_setattr(self, TableList, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, TableList, 'thisown', 0)
        _swig_setattr(self, TableList,self.__class__,TableList)
_Obit.TableList_swigregister(TableListPtr)


TableMF = _Obit.TableMF

TableNI = _Obit.TableNI

TableNX = _Obit.TableNX

TableSN = _Obit.TableSN

TableSU = _Obit.TableSU

TableVL = _Obit.TableVL

TableVZ = _Obit.TableVZ

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
    def __repr__(self):
        return "<C UVDesc instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.UVDesc_me_set
    __swig_getmethods__["me"] = _Obit.UVDesc_me_get
    if _newclass:me = property(_Obit.UVDesc_me_get, _Obit.UVDesc_me_set)
    def __init__(self, *args):
        _swig_setattr(self, UVDesc, 'this', _Obit.new_UVDesc(*args))
        _swig_setattr(self, UVDesc, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_UVDesc):
        try:
            if self.thisown: destroy(self)
        except: pass

class UVDescPtr(UVDesc):
    def __init__(self, this):
        _swig_setattr(self, UVDesc, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, UVDesc, 'thisown', 0)
        _swig_setattr(self, UVDesc,self.__class__,UVDesc)
_Obit.UVDesc_swigregister(UVDescPtr)


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
    def __repr__(self):
        return "<C UVImager instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.UVImager_me_set
    __swig_getmethods__["me"] = _Obit.UVImager_me_get
    if _newclass:me = property(_Obit.UVImager_me_get, _Obit.UVImager_me_set)
    def __init__(self, *args):
        _swig_setattr(self, UVImager, 'this', _Obit.new_UVImager(*args))
        _swig_setattr(self, UVImager, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_UVImager):
        try:
            if self.thisown: destroy(self)
        except: pass

class UVImagerPtr(UVImager):
    def __init__(self, this):
        _swig_setattr(self, UVImager, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, UVImager, 'thisown', 0)
        _swig_setattr(self, UVImager,self.__class__,UVImager)
_Obit.UVImager_swigregister(UVImagerPtr)


UVSetFITS = _Obit.UVSetFITS

UVSetAIPS = _Obit.UVSetAIPS

UVCastData = _Obit.UVCastData

UVCreate = _Obit.UVCreate

UVScratch = _Obit.UVScratch

UVZap = _Obit.UVZap

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

UVSetDesc = _Obit.UVSetDesc

UVisScratch = _Obit.UVisScratch

UVIsA = _Obit.UVIsA

UVGetName = _Obit.UVGetName

UVUtilUVWExtrema = _Obit.UVUtilUVWExtrema

UVUtilCopyZero = _Obit.UVUtilCopyZero

UVUtilVisDivide = _Obit.UVUtilVisDivide

UVUtilVisCompare = _Obit.UVUtilVisCompare

UVEditTD = _Obit.UVEditTD

UVEditStokes = _Obit.UVEditStokes

UVEditClip = _Obit.UVEditClip

UVEditClipStokes = _Obit.UVEditClipStokes
class UV(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UV, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UV, name)
    def __repr__(self):
        return "<C UV instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.UV_me_set
    __swig_getmethods__["me"] = _Obit.UV_me_get
    if _newclass:me = property(_Obit.UV_me_get, _Obit.UV_me_set)
    def __init__(self, *args):
        _swig_setattr(self, UV, 'this', _Obit.new_UV(*args))
        _swig_setattr(self, UV, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_UV):
        try:
            if self.thisown: destroy(self)
        except: pass

class UVPtr(UV):
    def __init__(self, this):
        _swig_setattr(self, UV, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, UV, 'thisown', 0)
        _swig_setattr(self, UV,self.__class__,UV)
_Obit.UV_swigregister(UVPtr)


newUVSelfCal = _Obit.newUVSelfCal

UVSelfCalCreate = _Obit.UVSelfCalCreate

UVSelfCalCopy = _Obit.UVSelfCalCopy

UVSelfCalUnref = _Obit.UVSelfCalUnref

UVSelfCalRef = _Obit.UVSelfCalRef

UVSelfCalGetList = _Obit.UVSelfCalGetList

UVSelfCalGetSkyModel = _Obit.UVSelfCalGetSkyModel

UVSelfCalSetSkyModel = _Obit.UVSelfCalSetSkyModel

UVSelfCalIsA = _Obit.UVSelfCalIsA

UVSelfCalCal = _Obit.UVSelfCalCal

UVSelfCalRefAnt = _Obit.UVSelfCalRefAnt

UVSelfSNSmo = _Obit.UVSelfSNSmo

UVSelfCalDeselSN = _Obit.UVSelfCalDeselSN
class UVSelfCal(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, UVSelfCal, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, UVSelfCal, name)
    def __repr__(self):
        return "<C UVSelfCal instance at %s>" % (self.this,)
    __swig_setmethods__["me"] = _Obit.UVSelfCal_me_set
    __swig_getmethods__["me"] = _Obit.UVSelfCal_me_get
    if _newclass:me = property(_Obit.UVSelfCal_me_get, _Obit.UVSelfCal_me_set)
    def __init__(self, *args):
        _swig_setattr(self, UVSelfCal, 'this', _Obit.new_UVSelfCal(*args))
        _swig_setattr(self, UVSelfCal, 'thisown', 1)
    def __del__(self, destroy=_Obit.delete_UVSelfCal):
        try:
            if self.thisown: destroy(self)
        except: pass

class UVSelfCalPtr(UVSelfCal):
    def __init__(self, this):
        _swig_setattr(self, UVSelfCal, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, UVSelfCal, 'thisown', 0)
        _swig_setattr(self, UVSelfCal,self.__class__,UVSelfCal)
_Obit.UVSelfCal_swigregister(UVSelfCalPtr)


