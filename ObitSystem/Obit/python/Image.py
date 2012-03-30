""" Python Obit Image class

This class contains an astronomical image and allows access.
An ObitImage is the front end to a persistent disk resident structure.
Magic value blanking is supported, blanked pixels have the value
OBIT_MAGIC (ObitImageDesc.h).
Pixel data are kept in an FArray structure which is how Python acceses the data.
There may be associated tables (e.g. "AIPS CC" tables).
Both FITS and AIPS cataloged images are supported.

Image Members with python interfaces:

==========  ==========================================================
exist       True if object previously existed prior to object creation
InfoList    used to pass instructions to processing
ImageDesc   Astronomical labeling of the image Member Desc 
FArray      Container used for pixel data Member FArray
PixBuf      memory pointer into I/O Buffer
Additional  Functions are available in ImageUtil.
==========  ==========================================================
"""
# Python/Obit Astronomical Image class
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2009
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

# Obit Image manupulation
import Obit, Table, FArray, OErr, ImageDesc, InfoList, History, AIPSDir, OSystem
import TableList, AIPSDir, FITSDir, ImageFit, FitRegion, FitModel
import OData
import SkyGeom, math
#import AIPSData

# Python shadow class to ObitImage class
 
# class name in C
myClass = "ObitImage"

class Image(OData.OData):
    """
    Python Obit Image class
    
    Additional Functions are available in ImageUtil.
    """
    def __init__(self, name) :
        self.myClass = myClass
        self.this = Obit.new_Image(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_Image(self.this)

    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.ImageUnref(Obit.Image_me_get(self.this))
            # In with the new
            Obit.Image_me_set(self.this,value)
            return
        if name=="FArray":
            PSetFArray(self, value)
        if name=="Beam":
            return PSetBeam(self, value)
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != Image:
            return
        if name == "me" : 
            return Obit.Image_me_get(self.this)
        # Functions to return members
        if name=="List":
            if not self.ImageIsA():
                raise TypeError,"input MUST be a Python Obit Image"
            out    = InfoList.InfoList()
            out.me = Obit.InfoListUnref(out.me)
            out.me = Obit.ImageGetList(self.cast(myClass))
            return out
        if name=="TableList":
            if not self.ImageIsA():
                raise TypeError,"input MUST be a Python Obit Image"
            out    = TableList.TableList("TL")
            out.me = Obit.TableListUnref(out.me)
            out.me = Obit.ImageGetTableList(self.cast(myClass))
            return out
        if name=="Desc":
            if not self.ImageIsA():
                raise TypeError,"input MUST be a Python Obit Image"
            out    = ImageDesc.ImageDesc("None")
            out.me = Obit.ImageGetDesc(self.cast(myClass))
            return out
        if name=="FArray":
            return PGetFArray(self)
        if name=="Beam":
            return PGetBeam(self)
        if name=="PixBuf":
            if not self.ImageIsA():
                raise TypeError,"input MUST be a Python Obit Image"
            fa = self.FArray
            return Obit.FArrayGetBuf(fa.me)
        raise AttributeError,str(name)  # Unknown
    def __repr__(self):
        if self.__class__ != Image:
            return
        return "<C Image instance> " + Obit.ImageGetName(self.cast(myClass))
    
    def cast(self, toClass):
        """
        Casts object pointer to specified class

        * self     = object whose cast pointer is desired
        * toClass  = Class string to cast to ("ObitImage")
        """
        # Get pointer with type of this class
        out = self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
            
    def Open (self, access, err, blc=None, trc=None):
        """
        Open an image persistent (disk) form

        * self   = Python Image object
        * access    = access READONLY (1), WRITEONLY (2), READWRITE(3)
        * err       = Python Obit Error/message stack
        * blc       = if given and a list of integers (min 2) giving
          bottom left corner (1-rel) of subimage
        * trc       = if given and a list of integers (min 2) giving
          top right corner (1-rel) of subimage
        """
        POpen(self, access, err, blc=blc, trc=trc)
        # end Open

    def Close (self, err):
        """
        Close an image  persistent (disk) form

        * self      = Python Image object
        * err       = Python Obit Error/message stack
        """
        PClose (self, err)
        # end Close

    def Read (self, err):
        """
        Read an image persistent (disk) form
        
        The data to be read is specified in the InfoList mamber
        Uses FArray member as buffer.

        * self      = Python Image object
        * err       = Python Obit Error/message stack
        """
        PRead (self, err)
        # end Read

    def Write (self, err):
        """
        Write an image  persistent (disk) form
        
        The data to be written is specified in the InfoList member
        Uses FArray member as buffer.

        * self      = Python Image object
        * err       = Python Obit Error/message stack
        """
        PWrite (self, err)
        # end Write
    
    def ReadFA (self, array, err):
        """
        Read an image  persistent (disk) form to a specified FArray
        
        The data to be read is specified in the InfoList member

        * self   = Python Image object
        * array  = Python FArray to accept data
        * err    = Python Obit Error/message stack
        """
        PReadFA (self, array, err)
        # end ReadFA

    def WriteFA (self, array, err):
        """
        Write an image  persistent (disk) form from a specified FArray
        
        The data to be written is specified in the InfoList member

        * self      = Python Image object
        * array     = Python FArray to write
        * err       = Python Obit Error/message stack
        """
        PWriteFA (self, array, err)
        # end WriteFA

    def ReadPlane (self, err, blc=None, trc=None):
        """
        Read an image plane into the FArray 
        
        Reads the plane specified by blc, trc
        into the FArray associated with the image

        * self     = Python Image object
        * err      = Python Obit Error/message stack
        * blc      = if given and a list of integers (min 2) giving
          bottom left corner (1-rel) of subimage
        * trc      = if given and a list of integers (min 2) giving
          top right corner (1-rel) of subimage

        returns Python  FArray from Image with data read
        """
        return PReadPlane (self, err, blc, trc)
        # end PReadPlane
   
    def WritePlane (self, imageData, err):
        """
        Write an image plane.
        
        Writes the plane specified by blc, trc on image infoList
        Checks if the current FArray on Image is compatable with
        imageData.

        * self      = Python Image object
        * imageData = Python FArray with data to write
        * err       = Python Obit Error/message stack
        """
        PWritePlane (self, imageData, err)
        # end PWritePlane

    def GetPlane (self, array, plane, err):
        """
        Read an image persistent (disk) form to an (optional) specified FArray
        
        The data to be read is specified in the InfoList member as modified by plane

        * self   = Python Image object
        * array  = Python FArray to accept data, if None use inImage buffer
        * plane  = array of 5 integers giving (1-rel) pixel numbers
        * err    = Python Obit Error/message stack
        """
        PGetPlane (self, array, plane, err)
        # end PGetPlane

    def PutPlane (self, array, plane, err):
        """
        Write an image persistent (disk) form from an (optional) specified FArray
        
        The data to be written is specified in the InfoList member as modified by plane

        * self      = Python Image object
        * array     = Python FArray to provide data, if None use inImage buffer
        * plane     = array of 5 integers giving (1-rel) pixel numbers
        * err       = Python Obit Error/message stack
        """
        PPutPlane (self, array, plane, err)
        # end PutPlane
        
    def GetPixel (self, pixel, err):
        """
        Return the specified pixel value
        
        Return pixel value

        * self     = Image
        * pixel    = pixel coordinate (1-rel int) as [1,1,1,1,1,1,1]
        * err      = Obit error stack
        """
        ################################################################
        # Checks
        if not self.ImageIsA():
            raise TypeError,"self MUST be a Python Obit Image"
        #
        self.Open(READONLY,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error opening image")
        plane = pixel[2:]
        self.GetPlane(None, plane, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error reading image")
        out = self.FArray.get(pixel[0]-1, pixel[1]-1)
        self.Close(err)
        if err.isErr:
            OErr.printErrMsg(err, "Error closing image")
        return out
    # end GetPixel

    def Copy (self, outImage, err):
        """
        Make a deep copy of input object.
        
        Makes structure the same as self, copies data, tables

        * self      = Python Image object to copy
        * outImage  = Output Python Image object, must be defined
        * err       = Python Obit Error/message stack
        """
        PCopy (self, outImage, err)
    # end Copy

    def Clone (self, outImage, err):
        """
        Make a copy of a object but do not copy the actual data
        
        This is useful to create an Image similar to the input one.

        * self      = Python Image object
        * outImage  = Output Python Image object, must be defined
        * err       = Python Obit Error/message stack
        """
        PClone (self, outImage, err)
        # end Clone

    def Scratch (self, err):
        """
        Create a scratch file suitable for accepting the data to be read from self
        
        A scratch Image is more or less the same as a normal Image except that it is
        automatically deleted on the final unreference.

        * self      = Python Image object
        * err       = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not self.ImageIsA():
            raise TypeError,"self MUST be a Python Obit Image"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        outImage    = Image("None")
        outImage.me = Obit.ImageScratch (self.cast(myClass), err.me);
        outImage.Info(err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating scratch file")

        return outImage
    # end Scratch

    def Header (self, err):
        """
        Write image header on output

        * self   = Python Obit Image object
        * err    = Python Obit Error/message stack
        """
        PHeader (self, err)
        # end Header

    def Info (self, err):
        """
        Get underlying data file info

        * self   = Python Obit Image object
        * err    = Python Obit Error/message stack
        """
        PImageInfo(self, err)
        # end Info
        
    def UpdateDesc (self, err, Desc=None):
        """
        Update any disk resident structures about descriptor

        * self      = Python Image object
        * err       = Python Obit Error/message stack
        * Desc      = Descriptor, if None then use current descriptor
          Contents can be accessed throuth the Dict member
        """
        PUpdateDesc (self, err, Desc=Desc)
        # end UpdateDesc
        
    def ImageIsA (self):
        """
        Tells if input really a Python Obit Image
        
        return true, false (1,0)

        * self   = Python UV object
        """
        ################################################################
        # Allow derived types
        return Obit.ImageIsA(self.cast(myClass))

    def TVFit (self, disp, err, file=None):
        """
        Fit Gaussian models. setting initial model from the TV display
        
        Returns FitModel object after fitting

        * self      = Python Image object
        * disp      = Display to use to interactively set initial model
        * err       = Python Obit Error/message stack
        * file      = If given, the file to write the results to, else terminal
        """
        # Interactively set initial model
        fr = FitRegion.PSetup(self,disp,err)
        # Setup for fitting - define parameters
        imf = ImageFit.ImageFit("ImageFit")
        input = ImageFit.FitInput
        input["fitImage"]  = self
        input["fitRegion"] = fr
        input["MaxIter"]   = 100
        input["prtLv"]     = 0
        input["PosGuard"]  = 1.
        # Fit
        imf.Fit(err,input)
        # Show results
        fr.Print(self.Desc,file=file)
        return fr
        # end TVFit

    def GaussFit (self, err, \
                  plane=1, cen=None, dim=[20,20], x=[0.0], y=[0.0], flux=[100.0], \
                  gparm=[[3.,3.,0.]], file=None):
        """
        Fit Gaussian models, setting initial model from parameters
        
        The initial model is defined by x,y,flux, gparm, all lists of the same dimension
        giving the location, fluxes and sizes of the initial models.
        Defaults OK for single source at reference pixel in image.
        Only uses first plane
        Returns FitModel object after fitting

        * self      = Python Image object
        * err       = Python Obit Error/message stack
        * plane     = 1-rel plane
        * cen       = If given the 1-rel center pixel of the region to be fit
          If not given, the reference pixel of the image is used
        * dim       = dimension in pixels of the region to be fit
        * x         = offset in x (pixels) of initial Gaussians from cen
        * y         = offset in y (pixels) of initial Gaussians from cen
        * flux      = fluxes of initial Gaussians
        * gparm     = Initial Gaussian size [major, minor, PA] (pixel,pixel,deg)
        * file      = If given, the file to write the results to, else terminal
        """
        self.List.set("BLC",[1,1,plane])
        self.List.set("TRC",[0,0,plane])
        self.Open(1,err)
        self.Read(err)
        self.Close(err)
        # Fitting region
        d = self.Desc.Dict
        if cen==None:
            cen = [d["crpix"][0], d["crpix"][1]]
        corner = [int(cen[0]-dim[0]/2), int(cen[1]-dim[1]/2)]
        # set initial model
        i=0
        fm = []
        for i in range(0,len(x)):
            flx = min(flux[i],d["maxval"])
            fm.append(FitModel.FitModel(type=FitModel.GaussMod, Peak=flx, \
                                        parms=gparm[i], \
                                        DeltaX=x[i]+dim[0]/2, DeltaY=y[i]+dim[1]/2));
        fr = FitRegion.FitRegion(corner=corner,dim=dim,models=fm)
        # Setup for fitting - define parameters
        imf = ImageFit.ImageFit("ImageFit")
        input = ImageFit.FitInput
        input["fitImage"]  = self
        input["fitRegion"] = fr
        input["MaxIter"]   = 100
        input["prtLv"]     = 0
        input["PosGuard"]  = 1.
        # Fit
        imf.Fit(err,input)
        # Show results
        fr.Print(self.Desc,file=file)
        return fr
        # end GaussFit

    # End of class member functions (i.e. invoked by x.func()(


# Commonly used, dangerous variables
dim = [1,1,1,1,1]
blc = [1,1,1,1,1,1,1]
trc = [0,0,0,0,0,0,0]
err = OErr.OErr()

# Symbolic names for access codes
READONLY  = OData.READONLY  # 1
WRITEONLY = OData.WRITEONLY # 2
READWRITE = OData.READWRITE # 3

def ObitName(ObitObject):
    """
    Return name of an Obit object or input if not an Obit Object
    """
    ################################################################
    out = ObitObject    # in case
    if type(out) == types.StringType:
        return out
    if ObitObject.me.find("_ObitImage_p") >= 0:
        return Obit.ImageGetName(ObitObject.me)
    if ObitObject.me.find("_ObitOTF_p") >= 0:
        return Obit.OTFGetName(ObitObject.me)
    if ObitObject.me.find("_ObitTable_p") >= 0:
        return Obit.TableGetName(ObitObject.me)
    if ObitObject.me.find("_Obit_p") >= 0:
        return Obit.GetName(ObitObject.me)
    return out
    # end ObitName
        

def input(inputDict):
    """
    Print the contents of an input Dictionary

    * inputDict = Python Dictionary containing the parameters for a routine
    """
    ################################################################
    print 'Current values of entries'
    myList = inputDict.items()
    myList.sort()
    for k,v in myList:
        # Print name of Obit objects (or string)
        #if (type(v)==types.StringType):
        #    print '  ',k,' = ',ObitName(v)
        #else:
        print '  ',k,' = ',v
        
    # end input

def newObit(name, filename, disk, exists, err):
    """
    Create and initialize an Image structure
    
    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    Returns the Python Image object

    * name     = name desired for object (labeling purposes)
    * filename = name of FITS file
    * disk     = FITS directory number
    * exists   = if true then the file is opened and closed to verify
    * err      = Python Obit Error/message stack
    """
    ################################################################
    out = Image (name)
    Obit.ImageSetFITS(out.me, 2, disk, filename, blc, trc, err.me)
    if exists:
        Obit.ImagefullInstantiate (out.me, 1, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating Image object")
    # show any errors 
    #OErr.printErrMsg(err, "newObit: Error verifying file")
    out.FileType = 'FITS'
    out.FileName = filename
    out.Disk = disk
    return out      # seems OK
    # end newObit

    
def newPFImage(name, filename, disk, exists, err, verbose=True):
    """
    Create and initialize an FITS based Image structure
    
    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    isOK member set to indicate success
    Returns the Python Image object

    * name     = name desired for object (labeling purposes)
    * filename = name of FITS file
    * disk     = FITS directory number
    * exists   = if true then the file is opened and closed to verify
    * err      = Python Obit Error/message stack
    * verbose  = If true any give error messages, else suppress
    """
    ################################################################
    out = Image (name)
    out.isOK = True  # until proven otherwise
    # Does it really previously exist?
    out.exist = FITSDir.PExist(filename, disk, err)
    Obit.ImageSetFITS(out.me, 2, disk, filename, blc, trc, err.me)
    if exists:
        Obit.ImagefullInstantiate (out.me, 1, err.me)
    # show any errors if wanted
    if  verbose and err.isErr:
        out.isOK = False
        OErr.printErrMsg(err, "Error creating FITS image object")
    elif err.isErr:
        out.isOK = False
        OErr.PClear(err)  # Clear unwanted messages
    # Check if really uvtab data and not an image
    outd = out.Desc.Dict
    if outd["inaxes"][0]==777777701:
        out.isOK = False
        raise TypeError,"Error: Object probably uvtab (UV) data"

    # It work?
    if not out.isOK:
        return out
    
    out.FileType = 'FITS'
    out.FileName = filename
    out.Fname    = filename
    out.Disk = disk
    out.Otype  = "Image"
    return out      # seems OK
    # end newPFImage

    
def newPAImage(name, Aname, Aclass, disk, seq, exists, err, verbose=False):
    """
    Create and initialize an AIPS based Image structure
    
    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    Returns the Python Image object
    isOK member set to indicate success

    * name     = name desired for object (labeling purposes)
    * Aname    = AIPS name of file
    * Aclass   = AIPS class of file
    * seq      = AIPS sequence number of file
    * disk     = FITS directory number
    * exists   = if true then the file is opened and closed to verify
    * err      = Python Obit Error/message stack
    * verbose  = If true any give error messages, else suppress
    """
    ################################################################
    out = Image (name)
    out.isOK = True  # until proven otherwise
    cno = -1
    user = OSystem.PGetAIPSuser()
    # print "disk, aseq", disk, seq
    # Does it really previously exist?
    test = AIPSDir.PTestCNO(disk, user, Aname, Aclass, "MA", seq, err)
    out.exist = test>0
    if exists: # If user thinks file exists...
        if out.exist: # If file is defined in catalog -> verify that file exists
            OErr.PLog(err, OErr.Info, Aname + " image found. Now verifying...")
            if verbose: OErr.printErr(err)
            cno = AIPSDir.PFindCNO(disk, user, Aname, Aclass, "MA", seq, err)
            Obit.ImageSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)
            Obit.ImagefullInstantiate (out.me, 1, err.me)
            #print "found",Aname,Aclass,"as",cno
        else: # If file not defined in catalog -> error
            OErr.PLog(err, OErr.Error, Aname + " image does not exist")
            out.isOK = False
    else: # exists=False
        # Create new image entry in catalog; if image already defined, this
        # has no effect
        OErr.PLog(err, OErr.Info, "Creating new image: "+Aname+", "+Aclass)
        if verbose: OErr.printErr(err)
        cno = AIPSDir.PAlloc(disk, user, Aname, Aclass, "MA", seq, err)
        Obit.ImageSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)
        #print "assigned",Aname,Aclass,"to",cno

    # show any errors if wanted
    if verbose and err.isErr:
        out.isOK = False
        OErr.printErrMsg(err, "Error creating AIPS Image object")
    elif err.isErr:
        out.isOK = False
        OErr.PClear(err)  # Clear unwanted messages
    else: OErr.PClear(err) # Clear non-error messages

    # It work?
    if not out.isOK:
        return out
    
    # Add File info
    out.FileType = 'AIPS'
    out.Disk   = disk
    out.Aname  = Aname
    out.Aclass = Aclass
    out.Aseq   = seq 
    out.Otype  = "Image"
    out.Acno   = cno
    return out      # seems OK
    # end newPAImage

    
def newPACNO(disk, cno, exists, err, verbose=True):
    """
    Create and initialize an AIPS based Image structure
    
    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    Returns the Python Image object
    isOK member set to indicate success

    * disk     = AIPS directory number
    * cno      = AIPS catalog number
    * exists   = if true then the file is opened and closed to verify
    * err      = Python Obit Error/message stack
    * verbose  = If true any give error messages, else suppress
    """
    ################################################################
    out = Image ("AIPS Image")
    user = OSystem.PGetAIPSuser()
    out.isOK = True  # until proven otherwise
    # print "disk, aseq", disk, seq
    # Does it really previously exist?
    test = AIPSDir.PInfo(disk, user, cno, err)
    out.exist = test != None
    
    if exists:
        Obit.ImageSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)
        Obit.ImagefullInstantiate (out.me, 1, err.me)
    else:
        Obit.ImageSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)

    # show any errors if wanted
    if verbose and err.isErr:
        out.isOK = False
        OErr.printErrMsg(err, "Error finding AIPS catalog entry")
    elif err.isErr:
        out.isOK = False
        OErr.PClear(err)  # Clear unwanted messages

    # It work?
    if not out.isOK:
        return out
    
    # Add File info
    out.FileType = 'AIPS'
    out.Disk   = disk
    out.Acno   = cno
    # Lookup name etc
    s = AIPSDir.PInfo (disk, user, cno, err)
    # parse returned string
    Aname = s[0:12]
    Aclass = s[13:19]
    Aseq = int(s[20:25])
    Atype = s[26:28]
    out.Aname  = Aname
    out.Aclass = Aclass
    out.Aseq   = Aseq 
    out.Otype  = "Image"
    return out      # seems OK
    # end newPACNO

    
# Image utilities
def PReadPlane (inImage, err, blc=None, trc=None):
    """
    Read an image plane into the FArray 
    
    Reads the plane specified by blc, trc
    into the FArray associated with the image

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    * blc       = if given and a list of integers (min 2) giving
      bottom left corner (1-rel) of subimage
    * trc       = if given and a list of integers (min 2) giving
      top right corner (1-rel) of subimage

    returns Python  FArray from Image with data read
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,'inImage MUST be a Python Obit Image'
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Read image
    inImage.Open(READONLY, err, blc, trc)
    inImage.Read(err)
    #OErr.printErrMsg(err, "Error reading image")
    imageData = inImage.FArray         # Image FArray (data)
    inImage.Close(err)
    return imageData
    # end PReadPlane
   
def PWritePlane (Image, imageData, err):
    """
    Write an image plane.
    
    Writes the plane specified by blc, trc on image infoList
    Checks if the current FArray on Image is compatable with
    imageData.

    * Image     = Python Image object
    * imageData = Python FArray with data to write
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.ImageIsA():
        raise TypeError,"Image MUST be a Python Obit Image"
    if not FArray.PIsA(imageData):
        raise TypeError,"imageData MUST be a Python Obit FArray"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Write image
    Image.Open(READWRITE, err)
    # Check that FArrays are compatible
    idtest = PGetFArray(Image)
    if not FArray.PIsCompatable (imageData, idtest):
        raise RuntimeError,'Images incompatable'
    Image.WriteFA(imageData, err)
    Image.Close(err)
    #OErr.printErrMsg(err, "Error writing images")
    # end PWritePlane

def PZap (inImage, err):
    """
    Delete underlying files and the basic object.

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    """
    ################################################################
    inImage.Zap(err)
    # end PZap

def PScratch (inImage, err):
    """
    Create a scratch file suitable for accepting the data to be read from inImage
    
    A scratch Image is more or less the same as a normal Image except that it is
    automatically deleted on the final unreference.

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    """
    ################################################################
    return inImage.Scratch(err)
    # end PScratch

def PCopy (inImage, outImage, err):
    """
    Make a deep copy of input object.
    
    Makes structure the same as inImage, copies data, tables

    * inImage   = Python Image object to copy
    * outImage  = Output Python Image object, must be defined
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not outImage.ImageIsA():
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageCopy (inImage.cast(myClass), outImage.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying Image")
    # end PCopy

def PClone (inImage, outImage, err):
    """
    Make a copy of a object but do not copy the actual data
    
    This is useful to create an Image similar to the input one.

    * inImage   = Python Image object
    * outImage  = Output Python Image object, must be defined
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not outImage.ImageIsA():
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageClone (inImage.cast(myClass), outImage.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error cloning Image")
    # end PClone

def PClone2 (inImage1, inImage2, outImage, err):
    """
    Make a copy of a object but do not copy the actual data

    * inImage1  = Python Image object to clone
    * inImage2  = Python Image object whose geometry is to be used
    * outImage  = Output Python Image object, must be defined,
      will be defined as Memory only
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage1.ImageIsA():
        raise TypeError,"inImage1 MUST be a Python Obit Image"
    if not inImage2.ImageIsA():
        raise TypeError,"inImage2 MUST be a Python Obit Image"
    if not outImage.ImageIsA():
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageClone2 (inImage1.cast(myClass), inImage2.cast(myClass), \
                      outImage.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error cloning Image")
    # end PClone2

def PCloneMem (inImage, outImage, err):
    """
    Make a Memory only clone of an Image structure
    
    This is useful for temporary structures

    * inImage   = Python Image object
    * outImage  = Output Python Image object, must be defined
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not outImage.ImageIsA():
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageCloneMem (inImage.cast(myClass), outImage.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error cloning Image")
    # end PCloneMem

def PCopyQuantizeFITS (inImage, outImage, err, fract=0.25, quant=None, \
                       inHistory=None):
    """
    Make a copy of an image quantizing to a 16 or 32 bit integer
        FITS image

    * inImage   = Python Image object
    * outImage  = Output Python Image object, must be defined
      but not fully created
    * err       = Python Obit Error/message stack
    * fract     = quantization level as a fraction of the plane min. RMS
    * quant     = quantization level in image units, has precedence over fract
      None or <= 0 => use fract.
    * inHistory = if given a History object to copy to the output FITS header
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not outImage.ImageIsA():
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Add parameters
    inInfo = PGetList(inImage)    #
    dim = [1,1,1,1,1]
    if fract:
        InfoList.PAlwaysPutFloat  (inInfo, "factor",  dim, [fract])
        InfoList.PRemove          (inInfo, "quant")
    if quant:
        InfoList.PAlwaysPutFloat  (inInfo, "quant",  dim, [quant])
        InfoList.PRemove          (inInfo, "factor")

    # Copy data
    outImage.me = Obit.ImageUtilQuanFITS (inImage.cast(myClass), outImage.FileName, \
                                          outImage.Disk, err.me)

    # Input image info
    inImage.Open(READONLY, err)
    #OErr.printErrMsg(err, "Error opening input")
    
    # Open output
    outImage.Open(WRITEONLY, err)                # Open
    #OErr.printErrMsg(err, "Error opening output")

    # Copy history if requested
    if inHistory!=None:
        outHistory  = History.History("Output history", outImage.List, err)
        History.PCopy2Header (inHistory, outHistory, err)
        OErr.printErrMsg(err, "Error with history")
  
    # Close images
    inImage.Close(err)
    outImage.Close(err)
# end PCopyQuantizeFITS

def PCompare (in1Image, in2Image, err, plane=[1,1,1,1,1]):
    """
    Compare a plane of two images
    
    returns list [max. abs in1Image, max abs difference, RMS difference]

    * in1Image  = Python Image object
    * in2Image  = Python Image object, on output, the FArray contains the difference.
    * err       = Python Obit Error/message stack
    * plane     = plane to compare
    """
    ################################################################
    # Checks
    if not in1Image.ImageIsA():
        raise TypeError,"in1Image MUST be a Python Obit Image"
    if not in2Image.ImageIsA():
        raise TypeError,"in2Image MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Input images, read plane 
    in1Image.Open(READONLY, err)
    #OErr.printErrMsg(err, "Error opening input")
    in1Image.GetPlane(None, plane, err)
    I1Data = in1Image.FArray 
    in2Image.Open(READONLY, err)
    #OErr.printErrMsg(err, "Error opening input")
    in2Image.GetPlane(None, plane, err)
    I2Data = in2Image.FArray
   
    # Close images
    in1Image.Close(err)
    in2Image.Close(err)
    #OErr.printErrMsg(err, "Error closing files")

    # Get difference
    FArray.PSub(I1Data, I2Data, I2Data)
    pos=[1,1]
    result = [FArray.PMaxAbs(I1Data,pos), FArray.PMaxAbs(I2Data,pos), FArray.PRMS(I2Data)]
    return result
# end PCompare

def PImageGetTable (inImage, access, tabType, tabVer, err):
    """
    Obsolete use PGetTable
    """
    return  PGetTable (inImage, access, tabType, tabVer, err)
# end  PImageGetTable

def PGetTable (inImage, access, tabType, tabVer, err,\
               noParms=0):
    """
    Return (create) the specified associated table
    
    Specific table types are recognized and the appropriate constructor
    called, these may have additional parameters.  This allows creating
    new tables of the appropriate type.
    returns Python Obit Table

    * inImage   = Python Image object
    * access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
    * tabType   = Table type, e.g. "AIPS AN", or "OTFSoln"
    * tabVer    = table version, if > 0 on input that table returned,
      if 0 on input, the highest version is used.
    * err       = Python Obit Error/message stack
    * noParms   = Number of parameters in CC table model
    """
    ################################################################
    return inImage.NewTable(access, tabType, tabVer, err, noParms=noParms)
    # end PGetTable
    

def PHeader (inImage, err):
    """
    Print image descriptor

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # ObitTalk or AIPSImage data?
    if inImage.myClass == 'AIPSImage':
        # header
        dict = inImage.header()
        ImageDesc.PHeaderDict(dict)
        # tablelist
        Tlist = inImage.tables()
        Tdict = {}
        # Once to find everything
        for item in Tlist:
            Tdict[item[1]] = item[0]
        # Again to get Max
        for item in Tlist:
            count = max (Tdict[item[1]], item[0])
            Tdict[item[1]] = count
        for item,count in Tdict.items():
            print "Maximum version number of %s tables is %d " % \
                  (item, count)
        return
        # End AIPSImage
    elif inImage.myClass == 'FITSImage':
        # header
        dict = inImage.header()
        ImageDesc.PHeaderDict(dict)
        # tablelist
        Tlist = inImage.tables()
        Tdict = {}
        # Once to find everything
        for item in Tlist:
            Tdict[item[1]] = item[0]
        # Again to get Max
        for item in Tlist:
            count = max (Tdict[item[1]], item[0])
            Tdict[item[1]] = count
        for item,count in Tdict.items():
            print "Maximum version number of %s tables is %d " % \
                  (item, count)
        return
        # End FITSImage
    
    # ObitTalk Image Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    #
    # Fully instantiate
    PFullInstantiate (inImage, READONLY, err)
    # File info
    if inImage.FileType=="AIPS":
        print "AIPS Image Name: %12s Class: %6s seq: %8d disk: %4d" % \
              (inImage.Aname, inImage.Aclass, inImage.Aseq, inImage.Disk)
    elif inImage.FileType=="FITS":
        print "FITS Image Disk: %5d File Name: %s " % \
              (inImage.Disk, inImage.FileName)
    # print in ImageDesc
    ImageDesc.PHeader(inImage.Desc)
    # Tables
    TL = inImage.TableList
    Tlist = TableList.PGetList(TL, err)
    Tdict = {}
    # Once to find everything
    for item in Tlist:
        Tdict[item[1]] = item[0]
    # Again to get Max
    for item in Tlist:
        count = max (Tdict[item[1]], item[0])
        Tdict[item[1]] = count
    for item,count in Tdict.items():
        print "Maximum version number of %s tables is %d " % \
              (item, count)
    # end PHeader
    

def POpen (inImage, access, err, blc=None, trc=None):
    """
    Open an image persistent (disk) form

    * inImage   = Python Image object
    * access    = access READONLY (1), WRITEONLY (2), READWRITE(3)
    * err       = Python Obit Error/message stack
    * blc       = if given and a list of integers (min 2) giving
      bottom left corner (1-rel) of subimage
    * trc       = if given and a list of integers (min 2) giving
      top right corner (1-rel) of subimage
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Set subimage if given
    if (blc.__class__==list) | (trc.__class__==list):
        inInfo = inImage.List
        if (blc.__class__==list) & (len(blc)>1) & (blc[0].__class__==int):
            dim = [len(blc),1,1,1,1]
            InfoList.PAlwaysPutInt  (inInfo, "BLC", dim, blc)
        if (trc.__class__==list) & (len(trc)>1) & (trc[0].__class__==int):
            dim = [len(trc),1,1,1,1]
            InfoList.PAlwaysPutInt  (inInfo, "TRC", dim, trc)
    Obit.ImageOpen(inImage.cast(myClass), access, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening Image")
    # end POpen

def PClose (inImage, err):
    """
    Close an image  persistent (disk) form

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageClose (inImage.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing Image")
    # end PClose

def PDirty (inImage):
    """
    Mark Image as needing a header update to disk file

    * inImage     = Python Image object
    """
    ################################################################
    inImage.Dirty()
    # end PDirty

def PRead (inImage, err):
    """
    Read an image  persistent (disk) form
    
    The data to be read is specified in the InfoList mamber
    Uses FArray member as buffer.

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageRead (inImage.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading Image")
    # end PRead

def PWrite (inImage, err):
    """
    Write an image  persistent (disk) form
    
    The data to be written is specified in the InfoList member
    Uses FArray member as buffer.

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageWrite (inImage.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing Image")
    # end PWrite
    
def PReadFA (inImage, array, err):
    """
    Read an image  persistent (disk) form to a specified FArray
    
    The data to be read is specified in the InfoList member

    * inImage   = Python Image object
    * array     = Python FArray to accept data
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not FArray.PIsA(array):
        raise TypeError,"array MUST be a Python Obit FArray"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageReadFA (inImage.cast(myClass), array.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading Image")
    # end PReadFA

def PWriteFA (inImage, array, err):
    """
    Write an image  persistent (disk) form from a specified FArray
    
    The data to be written is specified in the InfoList member

    * inImage   = Python Image object
    * array     = Python FArray to write
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not FArray.PIsA(array):
        raise TypeError,"array MUST be a Python Obit FArray"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageWriteFA (inImage.cast(myClass), array.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing Image")
    # end PWriteFA

def PGetPlane (inImage, array, plane, err):
    """
    Read an image  persistent (disk) form to an (optional) specified FArray
    
    The data to be read is specified in the InfoList member as modified by plane

    * inImage   = Python Image object
    * array     = Python FArray to accept data, if None use inImage buffer
    * plane     = array of 5 integers giving (1-rel) pixel numbers
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not ((array == None) or FArray.PIsA(array)):
        raise TypeError,"array MUST be a Python Obit FArray or None"
    if len(plane) != 5:
        raise TypeError,"plane must have 5 integer elements"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    if array == None:
        tarray = inImage.FArray   # use image objects FArray
        larray = tarray.me
    else:
        larray = array.me
    Obit.ImageGetPlane (inImage.cast(myClass), larray, plane, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading Image plane")
    # end PGetPlane

def PPutPlane (inImage, array, plane, err):
    """
    Write an image persistent (disk) form from an (optional) specified FArray
    
    The data to be written is specified in the InfoList member as modified by plane

    * inImage   = Python Image object
    * array     = Python FArray to provide data, if None use inImage buffer
    * plane     = array of 5 integers giving (1-rel) pixel numbers
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not ((array == None) or FArray.PIsA(array)):
        raise TypeError,"array MUST be a Python Obit FArray or None"
    if len(plane) != 5:
        raise TypeError,"plane must have 5 integer elements"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    if array == None:
        tarray =  inImage.FArray  # use image objects FArray
        larray = tarray.me
    else:
        larray = array.me
    Obit.ImagePutPlane (inImage.cast(myClass), larray, plane, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing Image plane")
    # end PPutPlane

def PZapTable (inImage, tabType, tabVer, err):
    """
    Destroy specified table

    * inImage   = Python Image object
    * tabType   = Table type, e.g. "AIPS CC"
    * tabVer    = table version, integer
    * err       = Python Obit Error/message stack
    """
    ################################################################
    inImage.ZapTable (tabType, tabVer, err)
    # end PZapTable

def PCopyTables (inImage, outImage, exclude, include, err):
    """
    Copy Tabeles from one image to another

    * inImage   = Python Image object
    * outImage  = Output Python Image object, must be defined
    * exclude   = list of table types to exclude (list of strings)
      has priority
    * include   = list of table types to include (list of strings)
    * err       = Python Obit Error/message stack
    """
    ################################################################
    inImage.CopyTables (outImage, exclude, include, err)
    # end PCopyTables

def PUpdateTables (inImage, err):
    """
    Update any disk resident structures about the current tables

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    """
    ################################################################
    inImage.UpdateTables (err)
    # end PUpdateTables

def PUpdateDesc (inImage, err, Desc=None):
    """
    Update external representation of descriptor

    * inImage   = Python Image object
    * err       = Python Obit Error/message stack
    * Desc      = Descriptor, if None then use current descriptor
      Contents can be accessed throuth the Dict member
    """
    ################################################################
     # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    #
    # if Desc=None make copy of current contents
    if Desc == None:
        d = inImage.Desc.Dict
    else:
        d = Desc.Dict
    # Open for write
    inImage.Open(READWRITE,err)   # Open
    inImage.Desc.Dict = d         # Update header
    Obit.ImageDirty(inImage.cast(myClass))   # force update
    inImage.Close(err)            # Close to update
    # end PUpdateDesc

def PImageInfo (inImage, err):
    """
    Get file info for extant uv data object
    
    Fills in information on object, useful for scratch files

    * inImage = Python Image object
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # file info
    info = Obit.ImageInfo (inImage.cast(myClass), err.me);
    if err.isErr:
        OErr.printErrMsg(err, "Error creating scratch file")
    if info["type"]=="AIPS":
        inImage.FileType = 'AIPS'
        inImage.Disk   = info["disk"]
        inImage.Otype  = "Image"
        inImage.Acno   = info["CNO"]
        # Lookup name etc
        s = AIPSDir.PInfo (inImage.Disk, info["user"], inImage.Acno, err)
        # parse returned string
        inImage.Aname  = s[0:12]
        inImage.Aclass = s[13:19]
        inImage.Aseq   = 0
    
    if info["type"]=="FITS":
        inImage.FileType = 'FITS'
        inImage.Disk   = info["disk"]
        inImage.Otype  = "Image"
        inImage.FileName = info["filename"]
        inImage.Fname    = info["filename"]
    # end PImageInfo

def PFullInstantiate (inImage, access, err):
    """
    Fully instantiate an Image by opening and closing
    
    return 0 on success, else failure

    * inImage   = Python Image object
    * access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.ImagefullInstantiate (inImage.cast(myClass), access, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error verifying Images")
    return ret
    # end PFullInstantiate

def PFArray2Image (inArray, outImage, err):
    """
    Attach an FArray to an image and write it
    
    Very rudimentary header attached

    * inArray   = Python Image object
    * outImage  = Python Image to write
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FArray.PIsA(inArray):
        raise TypeError,"inArray MUST be a Python Obit FArray"
    if not outImage.ImageIsA():
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Get array info
    naxis = inArray.Naxis[0:2]
    naxis = [naxis[0], naxis[1], 1, 1, 1, 1, 1]
    #    
    # Set size on image descriptor
    desc = ImageDesc.PDefault("temp")
    descDict = desc.Dict                      # Input Python dict object
    descDict["inaxes"] = naxis                # Update size
    descDict["crpix"]  = [1.0+naxis[0]*0.5, 1.0+naxis[1]*0.5, 1.0, 1.0, 1.0, 1.0, 1.0]
    descDict["bitpix"] = -32  # output floating
    descDict["object"] = "Temp"
    outImage.Desc.Dict = descDict
    #
    # Write output image
    outImage.Open(WRITEONLY, err)
    outImage.WriteFA(inArray, err)
    outImage.Close(err)
    #OErr.printErrMsg(err, "Error writing padded image for "+PGetName(outImage))
    # end PFArray2Image

def PFArray2FITS (inArray, outFile, err, outDisk=1, oDesc=None ):
    """
    Write an FArray to a FITS image
    
    Very rudimentary header attached
    Returns image object

    * inArray   = Python FArray object
    * outFile   = Name of FITS file
    * outDisk   = FITS disk number
    * oDesc     = None or ImageDescriptor to be written
    * err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not FArray.PIsA(inArray):
        raise TypeError,"inArray MUST be a Python Obit FArray"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create image
    outImage  = newObit(outFile, outFile, outDisk, 0, err)
    #OErr.printErrMsg(err, "Error creating FITS image "+outFile)
    #
    # Get array info
    naxis = inArray.Naxis[0:2]
    naxis = [naxis[0], naxis[1], 1, 1, 1, 1, 1]
    #
    # Set size on image descriptor
    if oDesc!=None:
        desc = oDesc
        descDict = desc.Dict                      # Input Python dict object
    else:
        desc = ImageDesc.PDefault("temp")
        descDict = desc.Dict                      # Input Python dict object
        descDict["inaxes"] = naxis                # Update size
        descDict["crpix"]  = [1.0+naxis[0]*0.5, 1.0+naxis[1]*0.5, 1.0, 1.0, 1.0, 1.0, 1.0]
        descDict["bitpix"] = -32  # output floating
        descDict["object"] = "Temp"
        #desc.Dict = descDict
    pos=[0,0]# Update descriptor
    descDict["minval"] = FArray.PMin(inArray, pos)
    descDict["maxval"] = FArray.PMax(inArray, pos)
    outImage.Desc.Dict = descDict                 # Update descriptor on image
    #
    # Write output image
    outImage.Open(WRITEONLY, err)
    outImage.WriteFA(inArray, err)
    outImage.Close(err)
    return outImage
    # end PFArray2FITS

def PSwapAxis (inImage, err, ax1=3, ax2=4):
    """
    Swap axes on an image
    
    The order of two adjacent axes may be swapped if the dimensionality
    of at least one of them is 1

    * inImage  = Image whose axes are to be swapped
    * err      = Python Obit Error/message stack
    * ax1      = first (1-rel) axis number
    * ax2      = second (1-rel) axis number
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,'inImage MUST be a Python Obit Image'
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if not (abs(ax1-ax2)==1):
        raise TypeError,"Axes "+str(ax1)+" and "+str(ax2)+" not adjacent"
    PFullInstantiate (inImage, READWRITE, err)
    # Get descriptor as dict
    d = inImage.Desc.Dict
    naxes = d['inaxes']
    if not ((naxes[ax1-1]==1) or (naxes[ax2-1]==1)):
        raise TypeError,"Both axes n>1: "+str(naxes[ax1-1])+","+str(naxes[ax2-1])
    #
    items = ['inaxes', 'ctype', 'crval', 'crpix', 'cdelt', 'crota']
    items = ['inaxes']
    for item in items:
        v = d[item][ax1-1]
        d[item][ax1-1] = d[item][ax2-1]
        d[item][ax2-1] = v
        
    inImage.Open(READWRITE, err)
    inImage.Desc.Dict = d
    PUpdateDesc(inImage, err)
    inImage.Close(err)
    # End  PSwapAxis
    
def PRelabelGal (inImage, err):
    """
    Relabel an image in galactic coordinates

    Change descriptor from Equatorial (J2000) to Galactic coordinates
    Returns modified image
    * inImage  = Image to be relabeled
    * err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,'inImage MUST be a Python Obit Image'
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    d = inImage.Desc.Dict
    ra2000  = d["crval"][0]
    dec2000 = d["crval"][1]
    # to B1950
    (ra1950,dec1950) = SkyGeom.PJtoB (ra2000, dec2000)
    # to galactic
    (glon,glat) = SkyGeom.PEq2Gal (ra1950, dec1950)
    # Need rotation
    (glonr,glatr) = SkyGeom.PEq2Gal (ra1950, dec1950+10.0/3600)
    rot = -57.296*math.atan2((glonr-glon), (glatr-glat))
    # rot=-58.3 # V. close
    #print "rot",rot
    d["ctype"][0] = "GLON-SIN"
    d["ctype"][1] = "GLAT-SIN"
    d["crval"][0] = glon
    d["crval"][1] = glat
    d["crota"][1] = rot
    d["obsra"]    = glon
    d["obsdec"]   = glat
    inImage.Desc.Dict = d
    inImage.UpdateDesc(err)
    return inImage
    # End  PRelabelGal
    
def PGetList (inImage):
    """
    Return the member InfoList
    
    returns InfoList

    * inImage   = Python Image object
    """
    ################################################################
    return inImage.List
    # end PGetList

def PGetTableList (inImage):
    """
    Return the member tableList
    
    returns tableList

    * inImage   = Python Image object
    """
    ################################################################ 
    return inImage.TableList
    # end PGetTableList


def PGetDesc (inImage):
    """
    Return the member ImageDesc
    
    returns ImageDesc as a Python Dictionary

    * inImage   = Python Image object
    """
    ################################################################
    return inImage.Desc
    # end PGetDesc

def PUpdateDesc (inImage, err, Desc=None):
    """
    Update external representation of descriptor

    * inImage = Python Image object
    * err     = Python Obit Error/message stack
    * Desc    = Image descriptor, if None then use current descriptor
    """
    ################################################################
     # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    #
    # if Desc=None make copy of current contents
    if Desc == None:
        d = inImage.Desc.Dict
        l = inImage.Desc.List
    else:
        d = Desc.Dict
        l = Desc.List
    # Open for write
    inImage.Open(READWRITE,err)    # Open
    inImage.Desc.Dict = d          # Update header
    InfoList.PCopyData(l, inImage.Desc.List)  # Copy list
    Obit.ImageDirty(inImage.cast(myClass))    # force update
    inImage.Close(err)             # Close to update
    # end PUpdateDesc

def PGetFArray (inImage):
    """
    Return FArray used to buffer Image data
    
    returns FArray with image pixel data

    * inImage   = Python Image object
    """
    ################################################################
     # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    #
    out    = FArray.FArray("None")
    out.me = Obit.ImageGetFArray(inImage.cast(myClass))
    return out
    # end PGetFArray

def PSetFArray (inImage, array):
    """
    Replace the FArray on an Image

    * inImage   = Python Image object
    * array     = Python FArray to attach
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not FArray.PIsA(array):
        raise TypeError,"array MUST be a Python Obit FArray"
    #
    Obit.ImageSetFArray(inImage.cast(myClass), array.me)
    # end PSetFArray

def PGetPixBuf (inImage):
    """
    Return python memory buffer for pixel array in memory

    * inImage   = Python Image object
    """
    ################################################################
    return inImage.PixBuf
    # end PGetPixBuf

def PGetBeam (inImage):
    """
    Return Beam attached to Image
    
    returns Beam with image pixel data

    * inImage   = Python Image object
    """
    ################################################################
     # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    #
    out    = Image("None")
    out.me = Obit.ImageGetBeam(inImage.cast(myClass))
    return out
    # end PGetBeam

def PSetBeam (inImage, beam):
    """
    Replace the Beam attached to an Image

    * inImage   = Python Image object
    * beam      = Python Beam Image to attach
    """
    ################################################################
    # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not beam.ImageIsA():
        raise TypeError,"array MUST be a Python Obit Image"
    #
    Obit.ImageSetBeam(inImage.cast(myClass), beam.me)
    # end PSetBeam

def PGetHighVer (inImage, tabType):
    """
    Get highest version number of a specified Table
    
    returns highest tabType version number, 0 if none.

    * inImage   = Python Image object
    * tabType   = Table type, e.g. "OTFSoln"
    """
    ################################################################
    return inImage.GetHighVer(tabType)
    # end PGetHighVer

def PIsScratch (inImage):
    """
    Tells if Image is a scratch object
    
    return true, false (1,0)

    * inImage   = Python Image object
    """
    ################################################################
    return inImage.IsScratch ()
    # end PIsScratch

def PIsA (inImage):
    """
    Tells if input really a Python Obit Image
    
    return True, False (1,0)

    * inImage   = Python Image object
    """
    ################################################################
    try:
        return inImage.ImageIsA()
    except:
        return False
    # end PIsA

def PUnref (inImage):
    """
    Decrement reference count
    
    Decrement reference count which will destroy object if it goes to zero
    Python object stays defined.

    * inImage   = Python Image object
    """
    ################################################################
     # Checks
    if not inImage.ImageIsA():
        raise TypeError,"inImage MUST be a Python Obit Image"

    inImage.me = Obit.ImageUnref(inImage.cast(myClass))
    # end PUnref

def PGetName (inImage):
    """
    Tells Image object name (label)
    
    returns name as character string

    * inImage   = Python Image object
    """
    ################################################################
    return inImage.GetName()
    # end PGetName
