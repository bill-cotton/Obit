""" Python Obit ImageMF class

This class contains an astronomical image and allows access.
An ObitImageMF is the front end to a persistent disk resident structure.
Magic value blanking is supported, blanked pixels have the value
OBIT_MAGIC (ObitImageDesc.h).
Pixel data are kept in an FArray structure which is how Python acceses the data.
There may be associated tables (e.g. "AIPS CC" tables).
Both FITS and AIPS cataloged images are supported.

ImageMF Members with python interfaces:
exist     - True if object previously existed prior to object creation
InfoList  - used to pass instructions to processing
ImageDesc - Astronomical labeling of the image Member Desc 
FArray    - Container used for pixel data Member FArray
PixBuf    - memory pointer into I/O Buffer
Additional Functions are available in ImageUtil.
"""
# Python/Obit Astronomical ImageMF class
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2010
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

# Obit ImageMF manupulation
import Obit, Table, FArray, OErr, Image, ImageDesc, InfoList, History, AIPSDir, OSystem
import TableList, AIPSDir, FITSDir, ImageFit, FitRegion, FitModel
import OData
#import AIPSData

# Python shadow class to ObitImageMF class
 
# class name in C
myClass = "ObitImageMF"

class ImageMF(Image.Image):
    """ Python Obit ImageMF class

    Additional Functions are available in ImageUtil.
    """
    def __init__(self, name) :
        self.myClass = myClass
        self.this = Obit.new_ImageMF(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_ImageMF(self.this)

    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.ImageMFUnref(Obit.ImageMF_me_get(self.this))
            # In with the new
            Obit.ImageMF_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != ImageMF:
            return
        if name == "me" : 
            return Obit.ImageMF_me_get(self.this)
        # Functions to return members
        if name=="List":
            if not self.ImageMFIsA():
                raise TypeError,"input MUST be a Python Obit ImageMF"
            out    = InfoList.InfoList()
            out.me = Obit.InfoListUnref(out.me)
            out.me = Obit.ImageGetList(self.cast("ObitImage"))
            return out
        if name=="TableList":
            if not self.ImageMFIsA():
                raise TypeError,"input MUST be a Python Obit ImageMF"
            out    = TableList.TableList("TL")
            out.me = Obit.TableListUnref(out.me)
            out.me = Obit.ImageGetTableList(self.cast("ObitImage"))
            return out
        if name=="Desc":
            if not self.ImageMFIsA():
                raise TypeError,"input MUST be a Python Obit ImageMF"
            out    = ImageDesc.ImageDesc("None")
            out.me = Obit.ImageGetDesc(self.cast("ObitImage"))
            return out
        if name=="FArray":
            return PGetFArray(self.cast("ObitImage"))
        if name=="Beam":
            return PGetBeam(self.cast("ObitImage"))
        if name=="PixBuf":
            if not self.ImageMFIsA():
                raise TypeError,"input MUST be a Python Obit ImageMF"
            fa = self.FArray
            return Obit.FArrayGetBuf(fa.me)
        raise AttributeError,str(name)  # Unknown
    def __repr__(self):
        if self.__class__ != ImageMF:
            return
        return "<C ImageMF instance> " + Obit.ImageGetName(self.cast("ObitImage"))
    
    def cast(self, toClass):
        """ Casts object pointer to specified class
        
        self     = object whose cast pointer is desired
        toClass  = Class string to cast to ("ObitImageMF")
        """
        # Get pointer with type of this class
        out = self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
            
    def ImageMFIsA (self):
        """ Tells if input really a Python Obit ImageMF
        
        return true, false (1,0)
        self   = Python ImageMF object
        """
        ################################################################
        # Allow derived types
        return Obit.ImageMFIsA(self.cast(myClass))

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
    """Return name of an Obit object or input if not an Obit Object
    """
    ################################################################
    out = ObitObject    # in case
    if type(out) == types.StringType:
        return out
    if ObitObject.me.find("_ObitImageMF_p") >= 0:
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
    """ Print the contents of an input Dictionary

    inputDict = Python Dictionary containing the parameters for a routine
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

def newPFImageMF(name, filename, disk, exists, err, verbose=True):
    """ Create and initialize an FITS based ImageMF structure

    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    isOK member set to indicate success
    Returns the Python ImageMF object
    name     = name desired for object (labeling purposes)
    filename = name of FITS file
    disk     = FITS directory number
    exists   = if true then the file is opened and closed to verify
    err      = Python Obit Error/message stack
    verbose  = If true any give error messages, else suppress
    """
    ################################################################
    out = ImageMF (name)
    out.isOK = True  # until proven otherwise
    # Does it really previously exist?
    out.exist = FITSDir.PExist(filename, disk, err)
    Obit.ImageMFSetFITS(out.me, 2, disk, filename, blc, trc, err.me)
    if exists:
        Obit.ImagefullInstantiate (out.cast("ObitImage"), 1, err.me)
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
    # end newPFImageMF

    
def newPAImage(name, Aname, Aclass, disk, seq, exists, err, verbose=False):
    """ Create and initialize an AIPS based Image structure

    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    Returns the Python Image object
    isOK member set to indicate success
    name     = name desired for object (labeling purposes)
    Aname    = AIPS name of file
    Aclass   = AIPS class of file
    seq      = AIPS sequence number of file
    disk     = FITS directory number
    exists   = if true then the file is opened and closed to verify
    err      = Python Obit Error/message stack
    verbose  = If true any give error messages, else suppress
    """
    ################################################################
    out = ImageMF (name)
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
            Obit.ImageMFSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)
            Obit.ImagefullInstantiate (out.cast("ObitImage"), 1, err.me)
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
        Obit.ImageMFSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)
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
    """ Create and initialize an AIPS based Image structure

    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    Returns the Python Image object
    isOK member set to indicate success
    disk     = AIPS directory number
    cno      = AIPS catalog number
    exists   = if true then the file is opened and closed to verify
    err      = Python Obit Error/message stack
    verbose  = If true any give error messages, else suppress
    """
    ################################################################
    out = ImageMF ("AIPS Image")
    user = OSystem.PGetAIPSuser()
    out.isOK = True  # until proven otherwise
    # print "disk, aseq", disk, seq
    # Does it really previously exist?
    test = AIPSDir.PInfo(disk, user, cno, err)
    out.exist = test != None
    
    if exists:
        Obit.ImageMFSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)
        Obit.ImagefullInstantiate (out.cast("ObitImage"), 1, err.me)
    else:
        Obit.ImageMFSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)

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

def PHeader (inImage, err):
    """ Print image descriptor

    inImage   = Python Image object
    err       = Python Obit Error/message stack
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
    

def PFitSpec (inImage, err, antSize=0.0, nOrder=1):
    """ Fit spectrum to each pixel of an ImageMF

    inImage   = Python Image object
    antSize   = If > 0 make primary beam corrections assuming antenna
                diameter (m) antSize
    nOrder    = Order of fit, 0=intensity, 1=spectral index,
                2=also curvature
    err       = Python Obit Error/message stack
    """
    ################################################################
    inImage.List.set("nOrder",nOrder)
    Obit.ImageMFFitSpec(inImage.me, antSize, err.me)
    # end PFitSpec

def PIsA (inImage):
    """ Tells if input really a Python Obit ImageMF

    return True, False (1,0)
    inImage   = Python Image object
    """
    ################################################################
    try:
        return inImage.ImageMFIsA()
    except:
        return False
    # end PIsA

