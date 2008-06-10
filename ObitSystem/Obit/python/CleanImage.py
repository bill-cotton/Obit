""" Python Obit Image Clean class

This class does BGC-like image based cleans

CleanImage Members with python interfaces:
InfoList  - used to pass instructions to processing
Member List (readonly)
mosaic    - ImageMosaic, use PGetMosaic, PSetMosaic
"""
# $Id: CleanImage.py,v 1.12 2007/02/16 16:57:03 bcotton Exp $
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

# Obit CleanImage
import Obit, OErr, ImageMosaic, InfoList, UV, OWindow

# Python shadow class to ObitDConCleanImage class
 
class CleanImagePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.CleanImageUnref(Obit.CleanImage_me_get(self.this))
            # In with the new
            Obit.CleanImage_me_set(self.this,value)
            return
        if name=="Mosaic":
            PSetMosaic(self, value)
            return 
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != CleanImage:
            return
        if name == "me" : 
            return Obit.CleanImage_me_get(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        if name=="Number":
            return PGetNumber(self)
        if name=="Mosaic":
            return PGetMosaic(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != CleanImage:
            return
        return "<C CleanImage instance> " + Obit.CleanImageGetName(self.me)
#
class CleanImage(CleanImagePtr):
    """ Python Obit Image Clean class

    This class does BGC-like image based cleans

    CleanImage Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List (readonly)
    mosaic    - ImageMosaic, use PGetMosaic, PSetMosaic
    """
    def __init__(self, name) :
        self.this = Obit.new_CleanImage(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_CleanImage(self.this)

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

# Commonly used, dangerous variables
dim=[1,1,1,1,1]
blc=[1,1,1,1,1,1,1]
trc=[0,0,0,0,0,0,0]
err=OErr.OErr()

def input(inputDict):
    """ Print the contents of an input Dictionary

    inputDict = Python Dictionary containing the parameters for a routine
    There should be a member of the dictionary ('structure') with a value
    being a list containing:
    1) The name for which the input is intended (string)
    2) a list of tuples consisting of (parameter name, doc string)
       with an entry for each parameter in the dictionary.
       The display of the the inputs dictionary will be in the order of
       the tuples and display the doc string after the value.
       An example:
       Soln2CalInput={'structure':['Soln2Cal',[('InData','Input OTF'),
                                               ('soln','input soln table version'),
                                               ('oldCal','input cal table version, -1=none'),
                                               ('newCal','output cal table')]],
                      'InData':None, 'soln':0, 'oldCal':-1, 'newCal':0}
    """
    ################################################################
    structure = inputDict['structure']  # Structure information
    print 'Inputs for ',structure[0]
    for k,v in structure[1]:
        print '  ',k,' = ',inputDict[k],' : ',v
        
    # end input

def newObit(name, err):
    """ Create and initialize an CleanImage structure

    Create sky model object
    Returns the Python CleanImage object
    name     = name desired for object (labeling purposes)
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = CleanImage (name)
    return out      # seems OK
    # end newObit

def PCopy (inCleanImage, outCleanImage, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inCleanImage, copies pointers
    inCleanImage  = Python CleanImage object to copy
    outCleanImage = Output Python CleanImage object, must be defined
    err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python Obit CleanImage"
    if not PIsA(outCleanImage):
        raise TypeError,"outCleanImage MUST be a Python Obit CleanImage"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.CleanImageCopy (inCleanImage.me, outCleanImage.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error copyint CleanImage")
    # end PCopy

def PGetList (inCleanImage):
    """ Return the member InfoList

    returns InfoList
    inCleanImage  = Python CleanImage object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python Obit CleanImage"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.CleanImageGetList(inCleanImage.me)
    return out
    # end PGetList

def PGetMosaic (inCleanImage):
    """ Return the member mosaic

    returns ImageMosaic
    inCleanImage  = Python CleanImage object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python Obit CleanImage"
    #
    out    = ImageMosaic.ImageMosaic("None", 1)
    out.me = Obit.CleanImageGetImageMosaic(inCleanImage.me)
    return out
    # end PGetMosaic

def PSetMosaic (inCleanImage, mosaic):
    """ Replace an ImageMosaic in the CleanImage

    inCleanImage  = Python CleanImage object
    mosaic      = Python ImageMosaic to attach
    """
    ################################################################
    # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python ObitCleanImage"
    if not ImageMosaic.PIsA(mosaic):
        raise TypeError,"array MUST be a Python Obit ImageMosaic"
    #
    Obit.CleanImageSetImageMosaic(inCleanImage.me, mosaic.me)
    # end PSetMosaic

def PGetWindow (inCleanImage):
    """ Return the member OWindow

    returns OWindow
    inCleanImage  = Python CleanImage object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python Obit CleanImage"
    #
    out    = OWindow.OWindow()
    out.me = Obit.CleanImageGetWindow(inCleanImage.me)
    return out
    # end PGetWindow

def PSetWindow (inCleanImage, window):
    """ Replace OWindow in the CleanImage

    inCleanImage  = Python CleanImage object
    window        = Python OWindow to attach
    """
    ################################################################
    # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python ObitCleanImage"
    if not OWindow.PIsA(window):
        raise TypeError,"array MUST be a Python Obit OWindow"
    #
    Obit.CleanImageSetWindow(inCleanImage.me, window.me)
    # end PSetWindow

def PAddWindow (inCleanImage, field, window, err):
    """ Add a window to a field to be CLEANed

    inCleanImage  = Python CleanImage object
    field         = Which field is the window in?
    window        = set of 4 integers:
                    if window[0]<0 box is round and
                    window[1]=radius, [2,3] = center
                    else rectangular and
                    blc=(window[0],window[1]), trc= blc=(window[2],window[3])
    """
    ################################################################
    # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python ObitCleanImage"
    #
    Obit.CleanImageAddWindow(inCleanImage.me, field, window, err.me)
    if err.isErr:
        printErrMsg(err, "Error adding window")
    # end PAddWindow

def PCreate (name, mosaic, err):
    """ Create the parameters and underlying structures of a CleanImage.

    Note: The dirty image will be replaced by the CLEAN image by
    the Clean.
    name      = Name to be given to object
                Most control parameters are in InfoList member
    mosaic    = Python ImageMosaic to attach
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not ImageMosaic.PIsA(mosaic):
        raise TypeError,"mosaic MUST be a Python Obit ImageMosaic"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
   #
    out = CleanImage(name);
    out.me = Obit.CleanImageCreate(name, mosaic.me,  err.me)
    if err.isErr:
        printErrMsg(err, "Error creating CleanImage")
    return out
    # end PCreate

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
        raise TypeError,"mosaic MUST be a Python Obit CleanImage"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.CleanImageDefWindow(clean.me,  err.me)
    if err.isErr:
        printErrMsg(err, "Error defining window")
    # end PDefWindow

# Perform Clean
CleanInput={'structure':['Clean',[('CleanImage','CleanImage Object'),
                                  ('Niter','Maximum number of CLEAN iterations'),
                                  ('minPatch','Minimum beam patch in pixels [def 100]'),
                                  ('maxPixel','Maximum number of residuals [def 20000]'),
                                  ('BMAJ','Restoring beam major axis (deg)'),
                                  ('BMIN','Restoring beam minor axis (deg)'),
                                  ('BPA','Restoring beam position angle (deg)'),
                                  ('Gain','CLEAN loop gain'),
                                  ('minFlux','Minimun flux density (Jy)'),
                                  ('Factor','CLEAN depth factor'),
                                  ('Plane','Plane being processed, 1-rel indices of axes 3-?'),
                                  ('CCVer','CC table version number [0 => highest]')]],
            # defaults
            'CleanImage':None,
            'Niter':100,
            'minPatch':100,
            'maxPixel':20000,
            'BMAJ':0.0,
            'BMIN':0.0,
            'BPA':0.0,
            'Gain':0.1,
            'minFlux':0.0,
            'Factor':0.0,
            'Plane':[1,1,1,1,1],
            'CCVer':0}
def PClean (err, input=CleanInput):
    """ Performs image based CLEAN

    The peak in the image is iteratively found and then a limited
    region of the beam times a fraction of the peak is subtracted
    and the process is iterated.  Occasionally a proper residual
    image is formed by gridding the components found, FFTing them,
    multiplying by the FFT of the dirty beam, FFTing back and
    subtracting from the previous residual image.
    This process only can deconvolve a single field and a quarter of
    that.  However, assymetric (but stationary) dirty beams are
    allowed.  The dirty bean must be attached to its dirty image.
    The dirty image passed is replaced by the CLEAN image.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary
    
    Input dictionary entries:
    CleanImage  = Input CleanImage,
    Niter       = Maximum number of CLEAN iterations
    minPatch    = Minimum beam patch in pixels [def 100]
    maxPixel    = Maximum number of residuals [def 20000]
    BMAJ        = Restoring beam major axis (deg)
    BMIN        = Restoring beam minor axis (deg)
    BPA         = Restoring beam position angle (deg)
    Gain        = CLEAN loop gain
    minFlux     = Minimun flux density (Jy)
    Factor      = CLEAN depth factor
    Plane       = Plane being processed, 1-rel indices of axes 3-?
    CCVer       = CC table version number

    """
    ################################################################
    # Get input parameters
    inCleanImage  = input["CleanImage"]
    # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python Obit CleanImage"
    #
    dim = [1,1,1,1,1]
    #
    # Set control values on CleanImage
    dim[0] = 1;
    inInfo = PGetList(inCleanImage)    # 
    InfoList.PPutInt   (inInfo, "Niter",    dim, [input["Niter"]],    err)
    InfoList.PPutInt   (inInfo, "minPatch", dim, [input["minPatch"]], err)
    InfoList.PPutInt   (inInfo, "maxPixel", dim, [input["maxPixel"]], err)
    InfoList.PPutInt   (inInfo, "CCVer",    dim, [input["CCVer"]],    err)
    InfoList.PPutFloat (inInfo, "BMAJ",     dim, [input["BMAJ"]],     err)
    InfoList.PPutFloat (inInfo, "BMIN",     dim, [input["BMIN"]],     err)
    InfoList.PPutFloat (inInfo, "BPA",      dim, [input["BPA"]],      err)
    InfoList.PPutFloat (inInfo, "Gain",     dim, [input["Gain"]],     err)
    InfoList.PPutFloat (inInfo, "minFlux",  dim, [input["minFlux"]],  err)
    InfoList.PPutFloat (inInfo, "Factor",   dim, [input["Factor"]],   err)
    dim[0] = len(input["Plane"])
    InfoList.PPutInt   (inInfo, "Plane",    dim, input["Plane"],      err)
    #
    # show any errors 
    #OErr.printErrMsg(err, "Clean: Error setting parameters")
    #
    # Do operation
    Obit.CleanImageDeconvolve(inCleanImage.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error deconvolving")
    # end PClean

def PGetName (inCleanImage):
    """ Tells Image object name (label)

    returns name as character string
    inCleanImage  = Python CleanImage object
    """
    ################################################################
     # Checks
    if not PIsA(inCleanImage):
        raise TypeError,"inCleanImage MUST be a Python Obit CleanImage"
    #
    return Obit.CleanImageGetName(inCleanImage.me)
    # end PGetName

def PIsA (inCleanImage):
    """ Tells if input really a Python Obit CleanImage

    return true, false (1,0)
    inCleanImage   = Python CleanImage object
    """
    ################################################################
    # Checks
    if inCleanImage.__class__ != CleanImage:
        return 0
    return Obit.CleanImageIsA(inCleanImage.me)
    # end PIsA
