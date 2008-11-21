""" This class is for performing CLEAN on images.

This implements an OTF image plane CLEAN
It is mostly appropriate for single dish images where the support of
the dirty beam is extremely limited.
There are no restrictions on the relative sizes of the dirty image and beam.

Arguments to the constructor:
name   - Name of the CLEAN object (a label)
dirty  - Python Obit Image dirty image object
beam   - Python Obit Image dirty beam object
clean  - Extant  Python Obit Image to receive, should be cloned from dirty
err    - Python Obit Error/message stack
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2005,2008
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

# Python shadow class to ObitDConCleanOTF class
import Obit, OErr, Image, FArray, Table, InfoList, OWindow, ODisplay

class CleanOTFPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
    #def __del__(self):
    #    if self.thisown == 1 :
    #        # If Obit has been unloaded don't bother
    #        if Obit.__class__ == Obit:
    #            Obit.delete_CleanOTF(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.CleanOTF_me_set(self.this,value)
            return
        if name=="Dirty":
            PSetDirty(self, value)
            return 
        if name=="Beam":
            PSetBeam(self, value)
            return 
        if name=="Clean":
            PSetClean(self, value)
            return 
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.CleanOTF_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="Dirty":
            return PGetDirty(self)
        if name=="Beam":
            return PGetBeam(self)
        if name=="Clean":
            return PGetClean(self)
        if name=="Size":
            return PGetCleanSize(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C CleanOTF instance>"
class CleanOTF(CleanOTFPtr):
    """ This class is for performing CLEAN on images.

    This implements an OTF image plane CLEAN
    It is mostly appropriate for single dish images where the support of
    the dirty beam is extremely limited.
    There are no restrictions on the relative sizes of the dirty image and beam.

    Arguments to the constructor:
    name   - Name of the CLEAN object (a label)
    dirty  - Python Obit Image dirty image object
    beam   - Python Obit Image dirty beam object
    clean  - Extant  Python Obit Image to receive, should be cloned from dirty
    err    - Python Obit Error/message stack
    """
    def __init__(self, name, dirty, beam, clean, err) :
        self.this = Obit.new_CleanOTF(name, dirty, beam, clean, err)
        #self.me   = Obit.new_CleanOTF(name, dirty, beam, clean, err)
        self.thisown = 1
    def __del__(self):
        if Obit!=None:
            Obit.delete_CleanOTF(self.this)

    def DefWindow(self, err):
        """ Set default window (all image)
        
        self   = Python OTF object
        err       = Python Obit Error/message stack
        """
        PDefWindow(self, err)
        # end DefWindow

    def AddWindow(self, window, err):
        """ Add a  window
        
        self   = Python OTF object
        window = set of 4 integers:
                 if window[0]<0 box is round and
                 window[1]=radius, [2,3] = center
                 else rectangular and
                 blc=(window[0],window[1]), trc= blc=(window[2],window[3])
        err    = Python Obit Error/message stack
        """
        PAddWindow(self, window, err)
        # end AddWindow

def PCreate (name, dirty, beam, clean, err):
    """ Create CleanOTF  Object

    returns CleanOTF object
    name    = Name for clean
    dirty   = Python Obit dirty image
    beam    = Python Obit dirty beam
              Must have same cell spacing are dirty but need not be same size
              if None, use Gaussian the size of beamMaj in dirty
    clean   = Python Obit CLEAN image
              Should be defined but need not be instantiated.
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks 
    if not Image.PIsA(dirty):
        raise TypeError,"dirty MUST be a Python Obit Image"
    if not Image.PIsA(clean):
        raise TypeError,"clean MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return None
    #
    if beam==None:
        lbeam = Image.Image("NoBeam")
    else:
        if not Image.PIsA(beam):
            raise TypeError,"beam MUST be a Python Obit Image"
        lbeam = beam
    out = CleanOTF(name, dirty.me, lbeam.me, clean.me, err.me)
    if beam:
        dirty.Beam = beam
  
    return out
    # end PCreate


def PGetWindow (inCleanOTF):
    """ Return the member OWindow

    returns OWindow
    inCleanOTF  = Python CleanOTF object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanOTF):
        raise TypeError,"inCleanOTF MUST be a Python Obit CleanOTF"
    #
    out    = OWindow.OWindow()
    out.me = Obit.CleanOTFGetWindow(inCleanOTF.me)
    return out
    # end PGetWindow

def PSetWindow (inCleanOTF, window):
    """ Replace OWindow in the CleanOTF

    inCleanOTF  = Python CleanOTF object
    window      = Python OWindow to attach
    """
    ################################################################
    # Checks
    if not PIsA(inCleanOTF):
        raise TypeError,"inCleanOTF MUST be a Python ObitCleanOTF"
    if not OWindow.PIsA(window):
        raise TypeError,"array MUST be a Python Obit OWindow"
    #
    Obit.CleanOTFSetWindow(inCleanOTF.me, window.me)
    # end PSetWindow

def PDefWindow (clean, err):
    """ Set default windows on image mosaic member.

    If mosaic member Radius>0 then make round boxes on Fly's eye field
    with this radius, else use rectangular box including all but outer 5 pixels
    On outlier fields, use rectangular box of width OutlierSize.
    Assumes all images in mosaic have descriptors defined.
    clean     = Clean object containing mosaic
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(clean):
        raise TypeError,"mosaic MUST be a Python Obit CleanOTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.CleanOTFDefWindow(clean.me,  err.me)
    # end PDefWindow

def PAddWindow (inCleanOTF, window, err):
    """ Add a window to be CLEANed

    inCleanOTF = Python CleanOTF object
    window     = set of 4 integers:
                 if window[0]<0 box is round and
                 window[1]=radius, [2,3] = center
                 else rectangular and
                 blc=(window[0],window[1]), trc= blc=(window[2],window[3])
    err        = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCleanOTF):
        raise TypeError,"inCleanOTF MUST be a Python ObitCleanOTF"
    if err.isErr: # existing error?
        return 
    #
    Obit.CleanOTFAddWindow(inCleanOTF.me, window, err.me)
    # end PAddWindow

# Perform Clean
CleanInput={'structure':['Clean',[('CleanOTF','CleanOTF Object'),
                                  ('disp','Image display to edit window'),
                                  ('Niter','Maximum number of CLEAN iterations'),
                                  ('Patch','Beam patch in pixels [def 100]'),
                                  ('BeamSize','Restoring beam FWHM (deg)'),
                                  ('Gain','CLEAN loop gain'),
                                  ('minFlux','Minimun flux density (Jy)'),
                                  ('noResid','If True do not include residuals in restored image'),
                                  ('doRestore','If True restore components'),
                                  ('doScale','If True scale residuals in restored image by beam areas'),
                                  ('Factor','CLEAN depth factor'),
                                  ('Plane','Plane being processed, 1-rel indices of axes 3-?'),
                                  ('autoWindow','Automatically set Windows?'),
                                  ('CCVer','CC table version number [0 => highest]'),
                                  ('scale','if True, scale CCs to units of restored CLEAN image')]],
            # defaults
            'CleanOTF':None,
            'disp':None,
            'Niter':100,
            'Patch':100,
            'BeamSize':0.0,
            'Gain':0.1,
            'minFlux':0.0,
            'noResid':False,
            'doRestore':True,
            'doScale':True,
            'Factor':0.0,
            'Plane':[1,1,1,1,1],
            'autoWindow':False,
            'CCVer':0,
            'scale':True}
def PClean (err, input=CleanInput):
    """ Performs image based CLEAN

    The peak in the image is iteratively found and then the beam
    times a fraction of the peak is subtracted and the process is iterated.  
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
    Input dictionary entries:
    CleanOTF    = Input CleanOTF,
    disp        = Image display to edit window
    Niter       = Maximum number of CLEAN iterations
    Patch       = Beam patch in pixels [def 100]
    maxPixel    = Maximum number of residuals [def 20000]
    BeamSize    = Restoring beam (deg)
    Gain        = CLEAN loop gain
    minFlux     = Minimun flux density (Jy)
    noResid     = If True do not include residuals in restored image
    doRestore   = If True restore components
    doScale     = If True scale residuals in restored image by beam areas
    Factor      = CLEAN depth factor
    Plane       = Plane being processed, 1-rel indices of axes 3-?
    autoWindow  = True if autoWindow feature wanted.
    CCVer       = CC table version number
    scale       = If True, scale CCs to units of restored CLEAN image
    """
    ################################################################
    # Get input parameters
    inCleanOTF  = input["CleanOTF"]
    # Checks
    if not PIsA(inCleanOTF):
        print "Really is",inCleanOTF.__class__
        raise TypeError,"inCleanOTF MUST be a Python Obit CleanOTF"
    if err.isErr: # existing error?
        return 
    #
    dim = [1,1,1,1,1]
    #
    # Set control values on CleanOTF
    dim[0] = 1;
    inInfo = PGetList(inCleanOTF)    # 
    InfoList.PAlwaysPutInt   (inInfo, "Niter",    dim, [input["Niter"]])
    InfoList.PAlwaysPutInt   (inInfo, "Patch",    dim, [input["Patch"]])
    InfoList.PAlwaysPutInt   (inInfo, "CCVer",    dim, [input["CCVer"]])
    InfoList.PAlwaysPutFloat (inInfo, "BeamSize", dim, [input["BeamSize"]])
    InfoList.PAlwaysPutFloat (inInfo, "Gain",     dim, [input["Gain"]])
    InfoList.PAlwaysPutFloat (inInfo, "minFlux",  dim, [input["minFlux"]])
    InfoList.PAlwaysPutFloat (inInfo, "Factor",   dim, [input["Factor"]])
    InfoList.PAlwaysPutBoolean (inInfo, "noResid",    dim, [input["noResid"]])
    InfoList.PAlwaysPutBoolean (inInfo, "doRestore",  dim, [input["doRestore"]])
    InfoList.PAlwaysPutBoolean (inInfo, "doScale",    dim, [input["doScale"]])
    InfoList.PAlwaysPutBoolean (inInfo, "doScaleCC",  dim, [input["scale"]])
    InfoList.PAlwaysPutBoolean (inInfo, "autoWindow", dim, [input["autoWindow"]])
    dim[0] = len(input["Plane"])
    InfoList.PAlwaysPutInt   (inInfo, "Plane",    dim, input["Plane"])
    #
    # show any errors 
    OErr.printErrMsg(err, "Clean: Error setting parameters")
    # Edit CLEAN window?
    disp = input["disp"]
    if disp!=None:
        window = PGetWindow(inCleanOTF)
        print "Display Dirty image for editing"
        ODisplay.PImage(disp, inCleanOTF.Dirty, err, window=window)
        OErr.printErrMsg(err, "Error editing CLEAN boxes")
    #
    # if Beam Given set on dirty image
    if inCleanOTF.Beam:
        dirty.Beam = inCleanOTF.Beam
    # Do operation
    Obit.CleanOTFClean(inCleanOTF.me, err.me)
    # end PClean

def PRestore (inCln, err):
    """ Restores components

    This is done automatically unless the restoring beam size is negative
    inCln   = Python Obit input OTFClean
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.CleanOTFRestore (inCln.me, err.me)
    # end PRestore


def PGetDirty (inCln):
    """ Get Dirty image

    returns Dirty image as Python Obit Image
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = Image.Image("None")
    out.me = Obit.CleanOTFGetDirty (inCln.me)
    return out
    # end PGetDirty

def PSetDirty (inCln, image):
    """ Set Dirty image

    inCln   = Python Obit input OTFClean
    image   = Python Obit Image for dirty image
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not Image.PIsA(image):
        raise TypeError,"Image MUST be a Python Obit Image"
    #
    Obit.CleanOTFSetDirty (inCln.me, image.me)
    # end PSetDirty


def PGetBeam (inCln):
    """ Get Beam image

    returns Beam image as Python Obit Image
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = Image.Image("None")
    out.me = Obit.CleanOTFGetBeam (inCln.me)
    return out
    # end PGetBeam


def PSetBeam (inCln, image):
    """ Set Beam image

    inCln   = Python Obit input OTFClean
    image   = Python Obit Image for dirty image
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not Image.PIsA(image):
        raise TypeError,"Image MUST be a Python Obit Image"
    #
    Obit.CleanOTFSetBeam (inCln.me, image.me)
    # end PSetDirty

def PGetClean (inCln):
    """ Get Clean image

    returns Clean image as Python Obit Image
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = Image.Image("None")
    out.me = Obit.CleanOTFGetClean (inCln.me)
    return out
    # end PGetClean


def PSetClean (inCln, image):
    """ Set Clean image

    inCln   = Python Obit input OTFClean
    image   = Python Obit Image for clean image
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    if not Image.PIsA(image):
        raise TypeError,"Image MUST be a Python Obit Image"
    #
    Obit.CleanOTFSetClean (inCln.me, image.me)
    # end PSetClean


def PGetList (inCln):
    """ Get InfoList

    return InfoList
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.CleanOTFGetList (inCln.me)
    return out
    # end PGetList


def PGetNiter (inCln):
    """ Get maximum number of CLEAN iterations

    This is only set after the CLEAN has run.
    returns maximum number of CLEAN iterations (int)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFGetNiter (inCln.me)
    # end PGetNiter

def PGetGain (inCln):
    """ Get CLEAN loop gain
 
    This is only set after the CLEAN has run.
    returns CLEAN loop gain (float)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFGetGain (inCln.me)
    # end PGetGain

def PGetFlux (inCln):
    """  Get min abs flux density for CLEAN

    This is only set after the CLEAN has run.
    returns min abs flux density for CLEAN (float)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFGetFlux (inCln.me)
    # end PGetFlux 

def PGetCleanSize (inCln):
    """ Get CLEAN restoring beam size

    This is only set after the CLEAN has run.
    returns CLEAN restoring beam size in pixels (float)
    inCln  = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if not PIsA(inCln):
        raise TypeError,"inCln MUST be a Python Obit OTFClean"
    #
    return Obit.CleanOTFGetCleanSize (inCln.me)
    # end PGetCleanSize


def PIsA (inCln):
    """ Tells if object thinks it's a Python ObitOTFClean

    return true, false (1,0)
    inCln   = Python Obit input OTFClean
    """
    ################################################################
    # Checks
    if inCln.__class__ != CleanOTF:
        return 0
    return Obit.CleanOTFIsA(inCln.me)
    # end PIsA

