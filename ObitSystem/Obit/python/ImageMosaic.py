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

# Obit ImageMosaic
import Obit, OErr, Image, InfoList, UV

# Python shadow class to ObitImageMosaic class
 
class ImageMosaicPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.ImageMosaicUnref(Obit.ImageMosaic_me_get(self.this))
            # In with the new
            Obit.ImageMosaic_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != ImageMosaic:
            return
        if name == "me" : 
            return Obit.ImageMosaic_me_get(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        if name=="Number":
            return PGetNumber(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != ImageMosaic:
            return
        return "<C ImageMosaic instance> " + Obit.ImageMosaicGetName(self.me)
#
class ImageMosaic(ImageMosaicPtr):
    """ Python Obit Image class

    This class contains an array of astronomical images and allows access.
    Both FITS and AIPS cataloged images are supported.

    ImageMosaic Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List (readonly)
    mosaic    - use PGetMosaic (readonly)
    uvdata    - use PGetUV (readonly)
    """
    def __init__(self, name, number) :
        self.this = Obit.new_ImageMosaic(name, number)
    def __del__(self):
        if Obit!=None:
            Obit.delete_ImageMosaic(self.this)

# Commonly used, dangerous variables
dim=[1,1,1,1,1]
blc=[1,1,1,1,1,1,1]
trc=[0,0,0,0,0,0,0]
err=OErr.OErr()

def newObit(name, number, err):
    """ Create and initialize an ImageMosaic structure

    Create array of images
    Returns the Python ImageMosaic object
    name     = name desired for object (labeling purposes)
    number   = Number of images in Mosaic
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = ImageMosaic (name, number)
    return out      # seems OK
    # end newObit

    
def PZapImage (inImageM, number, err):
    """ Zap (delete with underlying structures) selected image member(s).

    inImageM  = Python ImageMosaic object
    number    =  The 0-rel image number, -1=> all
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.ImageMosaicZapImage (inImageM.me, number, err.me)
    # end PZapImage

def PCopy (inImageM, outImageM, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inImage, copies pointers
    inImageM  = Python ImageMosaic object to copy
    outImageM = Output Python ImageMosaic object, must be defined
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    if not PIsA(outImageM):
        raise TypeError,"outImageM MUST be a Python Obit ImageMosaic"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.ImageMosaicCopy (inImageM.me, outImageM.me, err.me)
    # end PCopy


def PGetList (inImageM):
    """ Return the member InfoList

    returns InfoList
    inImageM  = Python ImageMosaic object
    """
    ################################################################
     # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.ImageMosaicGetList(inImageM.me)
    return out
    # end PGetList


def PGetImage (inImageM, number, err):
    """ Return the member Image

    returns Image
    inImageM  = Python ImageMosaic object
    number    = 0-rel image index
    err       = Python Obit Error/message stack
    """
    ################################################################
     # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    if err.isErr: # existing error?
        return
    #
    out    = Image.Image("None")
    out.me = Obit.ImageMosaicGetImage(inImageM.me, number, err.me)
    return out
    # end PGetImage

def PSetImage (inImageM, number, image):
    """ Replace an Image in the Mosaic

    inImageM  = Python ImageMosaic object
    number    = 0-rel image index
    image     = Python Image to attach
    """
    ################################################################
    # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    if not Image.PIsA(image):
        raise TypeError,"array MUST be a Python Obit Image"
    #
    Obit.ImageMosaicSetImage(inImageM.me, number, image.me, err.me)
    # end PSetImage

def PGetFullImage (inImageM, err):
    """ Return the full field member Image

    returns Image
    inImageM  = Python ImageMosaic object
    err       = Python Obit Error/message stack
    """
    ################################################################
     # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    #
    out    = Image.Image("None")
    if err.isErr: # existing error?
        return out
    out.me = Obit.ImageMosaicGetFullImage(inImageM.me, err.me)
    return out
    # end PGetFullImage

def PSetFullImage (inImage, image):
    """ Replace the full field member Image in the Mosaic

    inImageM  = Python ImageMosaic object
    image     = Python Image to attach
    """
    ################################################################
    # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    if not Image.PIsA(image):
        raise TypeError,"array MUST be a Python Obit Image"
    #
    Obit.ImageMosaicSetFullImage(inImageM.me, image.me, err.me)
    # end PSetFullImage

def PCreate (name, uvData, err):
    """ Create the parameters and underlying structures of a set of images.

    name      = Name to be given to object
    uvData    = Python uv data from which the image mosaic will be derived
                Most control parameters are in InfoList member
    err       = Python Obit err stack.
    """
    ################################################################
    # Checks
    if not UV.PIsA(uvData):
        raise TypeError,"uvData MUST be a Python Obit UV"
    #
    out = ImageMosaic("None", 1);
    if err.isErr: # existing error?
        return out
    out.me = Obit.ImageMosaicCreate(name, uvData.me, err.me)
    return out;
    # end PCreate

def PDefine (inImageM, uvData, doBeam, err):
    """ Define the parameters and underlying structures of a set of images.

    inImageM  = Python ImageMosaic object
    uvData    = Python uv data from which the image mosaic will be derived
    doBeam    = if True then make dirty beams
    err       = Python Obit err stack.
    """
    ################################################################
    # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    if not UV.PIsA(uvData):
        raise TypeError,"uvData MUST be a Python Obit UV"
    if err.isErr: # existing error?
        return
    #
    Obit.ImageMosaicDefine(inImageM.me, uvData.me, doBeam, err.me)
    # end PDefine

def PFlatten (inImageM, err):
    """ Project the tiles of a Mosaic to the full field flattened image.

    inImageM  = Python ImageMosaic object
    err       = Python Obit err stack.
    """
    ################################################################
    # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    if err.isErr: # existing error?
        return
    #
    Obit.ImageMosaicFlatten(inImageM.me, err.me)
    # end PFlatten

def PGetName (inImageM):
    """ Tells Image object name (label)

    returns name as character string
    inImageM  = Python ImageMosaic object
    """
    ################################################################
     # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    #
    return Obit.ImageMosaicGetName(inImageM.me)
    # end PGetName

def PGetNumber (inImageM):
    """ Tells number of images in Mosais

    returns number of images
    inImageM  = Python ImageMosaic object
    """
    ################################################################
     # Checks
    if not PIsA(inImageM):
        raise TypeError,"inImageM MUST be a Python Obit ImageMosaic"
    #
    return Obit.ImageMosaicGetNumber(inImageM.me)
    # end PGetName

def PIsA (inImageM):
    """ Tells if input really a Python Obit ImageMosaic

    return true, false (1,0)
    inImageM   = Python ImageMosaic object
    """
    ################################################################
    # Checks
    if inImageM.__class__ != ImageMosaic:
        return 0
    return Obit.ImageMosaicIsA(inImageM.me)
    # end PIsA

