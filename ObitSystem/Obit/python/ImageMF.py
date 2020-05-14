""" Python Obit ImageMF class

This class contains an astronomical image and allows access.
An ObitImageMF is the front end to a persistent disk resident structure.
Magic value blanking is supported, blanked pixels have the value
OBIT_MAGIC (ObitImageDesc.h).
Pixel data are kept in an FArray structure which is how Python acceses the data.
There may be associated tables (e.g. "AIPS CC" tables).
Both FITS and AIPS cataloged images are supported.

ImageMF Members with python interfaces:

==========  ==========================================================
exist       True if object previously existed prior to object creation
InfoList    used to pass instructions to processing
ImageDesc   Astronomical labeling of the image Member Desc 
FArray      Container used for pixel data Member FArray
PixBuf      memory pointer into I/O Buffer
Additional  Functions are available in ImageUtil.
==========  ==========================================================
"""
# Python/Obit Astronomical ImageMF class
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2010-2020
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
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, Table, FArray, OErr, Image, ImageDesc, InfoList, History, AIPSDir, OSystem
import TableList, AIPSDir, FITSDir, ImageFit, FitRegion, FitModel
import OData
#import AIPSData

# Python shadow class to ObitImageMF class
 
# class name in C
myClass = "ObitImageMF"

class ImageMF(Obit.ImageMF, Image.Image):
    """
    Python Obit ImageMF class
    
    Additional Functions are available in ImageUtil.
    """
    def __init__(self, name) :
        super(ImageMF, self).__init__()
        Obit.CreateImageMF(self.this, name)
        self.myClass = myClass
    def __del__(self, DeleteImageMF=_Obit.DeleteImageMF):
        if _Obit!=None:
            DeleteImageMF(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            if value==None:
                raise TypeError("None given for ImageMF object")
            # Out with the old
            if self.this!=None:
                Obit.ImageMFUnref(Obit.ImageMF_Get_me(self.this))
            # In with the new
            Obit.ImageMF_Set_me(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, ImageMF):
            return "Bogus dude "+str(self.__class__)
        if name == "me" : 
            return Obit.ImageMF_Get_me(self.this)
        # Functions to return members
        if name=="List":
            if not self.ImageMFIsA():
                raise TypeError("input MUST be a Python Obit ImageMF")
            out    = InfoList.InfoList()
            out.me = Obit.ImageMFGetList(self.me)
            return out
        if name=="TableList":
            if not self.ImageMFIsA():
                raise TypeError("input MUST be a Python Obit ImageMF")
            out    = TableList.TableList("TL")
            out.me = Obit.ImageMFGetTableList(self.me)
            return out
        if name=="Desc":
            if not self.ImageMFIsA():
                raise TypeError("input MUST be a Python Obit ImageMF")
            out    = ImageDesc.ImageDesc("None")
            out.me = Obit.ImageMFGetDesc(self.me)
            return out
        if name=="FArray":
            return PGetFArray(self.me)
        if name=="Beam":
            return PGetBeam(self.me)
        if name=="PixBuf":
            if not self.ImageMFIsA():
                raise TypeError("input MUST be a Python Obit ImageMF")
            fa = self.FArray
            return Obit.FArrayGetBuf(fa.me)
        raise AttributeError(str(name))  # Unknown
    def __repr__(self):
        if not isinstance(self, ImageMF):
            return "Bogus dude "+str(self.__class__)
        return "<C ImageMF instance> " + Obit.ImageMFGetName(self.me)
    
    def ImageMFIsA (self):
        """
        Tells if input really a Python Obit ImageMF
        
        return True, False 

        * self   = Python ImageMF object
        """
        ################################################################
        return Obit.ImageMFIsA(self.me)

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
    if type(out) == bytes:
        return out
    if ObitObject.me.find("ObitImageMF_p") >= 0:
        return Obit.ImageGetName(ObitObject.me)
    if ObitObject.me.find("ObitOTF_p") >= 0:
        return Obit.OTFGetName(ObitObject.me)
    if ObitObject.me.find("ObitTable_p") >= 0:
        return Obit.TableGetName(ObitObject.me)
    if ObitObject.me.find("Obit_p") >= 0:
        return Obit.GetName(ObitObject.me)
    return out
    # end ObitName
        

def input(inputDict):
    """
    Print the contents of an input Dictionary

    * inputDict = Python Dictionary containing the parameters for a routine
    """
    ################################################################
    print('Current values of entries')
    myList = list(inputDict.items())
    myList.sort()
    for k,v in myList:
        # Print name of Obit objects (or string)
        #if (type(v)==types.StringType):
        #    print '  ',k,' = ',ObitName(v)
        #else:
        print('  ',k,' = ',v)
        
    # end input

def newPFImageMF(name, filename, disk, exists, err, verbose=True):
    """
    Create and initialize an FITS based ImageMF structure
    
    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    isOK member set to indicate success
    Returns the Python ImageMF object

    * name     = name desired for object (labeling purposes)
    * filename = name of FITS file
    * disk     = FITS directory number
    * exists   = if true then the file is opened and closed to verify
    * err      = Python Obit Error/message stack
    * verbose  = If true any give error messages, else suppress
    """
    ################################################################
    out = ImageMF (name)
    out.isOK = True  # until proven otherwise
    # Does it really previously exist?
    out.exist = FITSDir.PExist(filename, disk, err)
    Obit.ImageMFSetFITS(out.me, 2, disk, filename, blc, trc, err.me)
    if exists:
        Obit.ImageMFfullInstantiate (out.me, 1, err.me)
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
        raise TypeError("Error: Object probably uvtab (UV) data")

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
            Obit.ImageMFfullInstantiate (out.me, 1, err.me)
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
    out = ImageMF ("AIPS Image")
    user = OSystem.PGetAIPSuser()
    out.isOK = True  # until proven otherwise
    # print "disk, aseq", disk, seq
    # Does it really previously exist?
    test = AIPSDir.PInfo(disk, user, cno, err)
    out.exist = test != None
    
    if exists:
        Obit.ImageMFSetAIPS(out.me, 2, disk, cno, user, blc, trc, err.me)
        Obit.ImageMFfullInstantiate (out.me, 1, err.me)
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
            print("Maximum version number of %s tables is %d " % \
                  (item, count))
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
            print("Maximum version number of %s tables is %d " % \
                  (item, count))
        return
        # End FITSImage
    
    # ObitTalk Image Checks
    if not inImage.ImageMFIsA():
        raise TypeError("inImage MUST be a Python Obit ImageMF")
    #
    # Fully instantiate
    #DAMN PFullInstantiate (inImage, READONLY, err)
    # File info
    if inImage.FileType=="AIPS":
        print("AIPS Image Name: %12s Class: %6s seq: %8d disk: %4d" % \
              (inImage.Aname, inImage.Aclass, inImage.Aseq, inImage.Disk))
    elif inImage.FileType=="FITS":
        print("FITS Image Disk: %5d File Name: %s " % \
              (inImage.Disk, inImage.FileName))
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
        print("Maximum version number of %s tables is %d " % \
              (item, count))
    # end PHeader
    

def PFitSpec (inImage, err, antSize=0.0, nOrder=1, corAlpha=0.0):
    """
    Fit spectrum to each pixel of an ImageMF

    * inImage   = Python Image object
    * antSize   = If > 0 make primary beam corrections assuming antenna
                  diameter (m) antSize
    * nOrder    = Order of fit, 0=intensity, 1=spectral index,
                  2=also curvature
    * corAlpha  = Spectral index correction to apply before fitting
    * err       = Python Obit Error/message stack
    """
    ################################################################
    inImage.List.set("nOrder",nOrder)
    inImage.List.set("corAlpha",corAlpha, ttype='float')
    Obit.ImageMFFitSpec(inImage.me, antSize, err.me)
    # end PFitSpec

def PFitSpec2 (inImage, outImage, err, nterm=2, \
               refFreq=None, maxChi2=2.0, doError=False, doBrokePow=False, \
               calFract=None, doPBCor=False, PBmin=None, antSize=None, \
               corAlpha=0.0, minWt=0.5, doTab=False):
    """
    Fit spectrum to each pixel of an ImageMF writing a new cube

    Fitted spectral polynomials returned in outImage
    Can run with multiple threads if enabled:
    OSystem.PAllowThreads(2)  # 2 threads
     * inImage   = Python Image object; 
     * outImage  = Image cube with fitted spectra.
                   Should be defined but not created.
                   Planes 1->nterm are coefficients per pixel
                   Planes nterm+1->2*nterm are uncertainties in coefficients
                   Plane 2*nterm+1 = Chi squared of fit
    * err       = Python Obit Error/message stack
    * nterm     = Number of terms in spectral fit, 2=SI, 3=curvature
    * refFreq   = Reference frequency for fit [def ref for inImage]
    * maxChi2   = Max. Chi Sq for accepting a partial spectrum [def 2.0]
    * doError   = If true do error analysis [def False]
    * doBrokePow= If true do broken power law (3 terms). [def False]
    * calFract  = Calibration error as fraction of flux
                  One per frequency or one for all, def 0.05
    * doPBCor   = If true do primary beam correction. [def False]
    * PBmin     = Minimum beam gain correction
                  One per frequency or one for all, def 0.05,
                  1.0 => no gain corrections
    * antSize   = Antenna diameter (m) for PB gain corr, 
                  One per frequency or one for all, def 25.0
    * corAlpha  = Spectral index correction to apply before fitting
    * minWt     = min. fract of possible weight per pixel
    * doTab     = Use tabulated beam if available
    """
    ################################################################
    inImage.List.set("nterm",      nterm,      ttype='long')
    inImage.List.set("maxChi2",    maxChi2,    ttype='float')
    inImage.List.set("doError",    doError,    ttype='boolean')
    inImage.List.set("doBrokePow", doBrokePow, ttype='boolean')
    inImage.List.set("doPBCor",    doPBCor,    ttype='boolean')
    inImage.List.set("corAlpha" ,  corAlpha,   ttype='float')
    inImage.List.set("minWt",      minWt,      ttype='float')
    inImage.List.set("doTab"  ,    doTab,      ttype='boolean')
    if refFreq:
        inImage.List.set("refFreq",  refFreq,  ttype='double')
    if calFract:
        inImage.List.set("calFract", calFract, ttype='float')
    if PBmin:
        inImage.List.set("PBmin",    PBmin,    ttype='float')
    if antSize:
        inImage.List.set("antSize",   antSize, ttype='float')

    Obit.ImageMFFitSpec2(inImage.me, outImage.me, err.me)
    # end PFitSpec2

def PEffFqCorr (fitImage, rawImage, err, corAlpha=-0.7, \
                doPBCor=False, PBmin=None, antSize=None, doTab=False):
    """
    Correct flux density plane for effective frequency

    Calculate the effective frequency in each pixel using  (optional) 
    primary beam correction and RMS weighting and correct the flux 
    density in each pixel to the reference frequency using alpha=corAlpha.
    Multiply flux density plane by 
    exp(alpha*ln(nu/nu_0) where nu is determined for each pixel by a weighted average
    weighted by (1/RMS^2)*(PB^2)
    Can run with multiple threads if enabled:
    OSystem.PAllowThreads(2)  # 2 threads
    * fitImage  = Python Image with fitted spectrum to be modified
    * rawImage  = Raw ImageMF
    * err       = Python Obit Error/message stack
    * corAlpha  = Spectral index correction to apply before fitting
    * doPBCor   = If true do primary beam correction. [def False]
    * PBmin     = Minimum beam gain correction
                  One per frequency or one for all, def 0.05,
                  1.0 => no gain corrections
    * antSize   = Antenna diameter (m) for PB gain corr, 
                  One per frequency or one for all, def 25.0
    * doTab     = Use tabulated beam if available
     """
    ################################################################
    fitImage.List.set("doPBCor",    doPBCor,    ttype='boolean')
    fitImage.List.set("corAlpha" ,  corAlpha,   ttype='float')
    fitImage.List.set("doTab"  ,    doTab,    ttype='boolean')
    if PBmin:
        fitImage.List.set("PBmin",    PBmin,    ttype='float')
    if antSize:
        fitImage.List.set("antSize",   antSize, ttype='float')

    Obit.ImageMFEffFqCorr(fitImage.me, rawImage.me, err.me)
    # end PEffFqCorr

def PMFPBCor (inIm, err, antSize=25, minGain=0.05):
    """
    Apply primary beam corrections to an ImageMF

    WARNING: This routine modifies the input image;
    ONLY RUN IT ONCE.
    * inIm    = Image to be modified
    * err     = Python Obit Error/message stack
    * antSize = Antenna diameter in m.
    * minGain = Minimum antenna gain
    """
    # Image info
    nterm = inIm.Desc.List.Dict['NTERM'][2][0]
    nspec = inIm.Desc.List.Dict['NSPEC'][2][0]
    freqs = []
    for i in range(1,nspec+1):
        key = 'FREQ%4.4d'%i
        freqs.append(inIm.Desc.List.Dict[key][2][0])
    
    # end loop
    # Make scratch image for beam
    beam = Image.Image("PBeam")
    Image.PCloneMem(inIm, beam, err)
    OErr.printErrMsg(err, "Error with scratch beam image")
    # Debug
    xf = Image.PFArray2FITS(beam.FArray, "Beam.fits", err, 0, oDesc=beam.Desc)
    # Loop over planes
    for i in range(1,nspec+1):
        # Read plane
        plane=[i+nterm,1,1,1,1]
        Image.PGetPlane (inIm, None, plane, err)
        OErr.printErrMsg(err, "Error reading image")
        # Set frequency for PBCor
        d = inIm.Desc.Dict; d['crval'][2] = freqs[i-1]; inIm.Desc.Dict = d
        # Make PB Image
        ImageUtil.PPBImage(beam, beam, err, minGain=minGain, antSize=antSize, outPlane=plane)
        OErr.printErrMsg(err, "Error making PB image")
        # Debug
        oldVal = inIm.FArray.get(5046,1533)
        gain = beam.FArray.get(5046,1533)
        # Divide
        FArray.PDivClip (inIm.FArray, beam.FArray, minGain, inIm.FArray)
        newVal = inIm.FArray.get(5046,1533)
        # Rewrite plane
        Image.PPutPlane (inIm, None, plane, err)
        # Debug
        print (i,  plane[0], "nu=",freqs[i-1],'g=', gain, oldVal, newVal)
        Image.PPutPlane (xf, beam.FArray, plane, err)
        OErr.printErrMsg(err, "Error writing image")
   
    # end loop

# End PMFPBCor
def PGetFArray (inImageMF):
    """
    Return FArray used to buffer Image data
    
    returns FArray with image pixel data
    * inImageMF   = Python Image object
    """
    ################################################################
    if inImageMF.myClass=='AIPSImage':
        raise TypeError("Function unavailable for "+inImage.myClass)
    # Checks
    if not inImageMF.ImageMFIsA():
        raise TypeError("inImageMF MUST be a Python Obit Image")
    #
    out    = FArray.FArray("None")
    out.me = Obit.ImageMFGetFArray(inImageMF.me)
    return out
    # end PGetFArray

def PGetBeam (inImage):
    """
    Return Beam attached to Image
    
    returns Beam with image pixel data
    * inImage   = Python Image object
    """
    ################################################################
    if inImage.myClass=='AIPSImage':
        raise TypeError("Function unavailable for "+inImage.myClass)
    # Checks
    if not inImage.ImageMFIsA():
        raise TypeError("inImage MUST be a Python Obit ImageMF")
    #
    out    = Image("None")
    out.me = Obit.ImageMFGetBeam(inImage.me)
    return out
    # end PGetBeam

def PIsA (inImage):
    """
    Tells if input really a Python Obit ImageMF
    
    return True, False
    * inImage   = Python Image object
    """
    ################################################################
    try:
        return inImage.ImageMFIsA()!=0
    except:
        return False
    # end PIsA

