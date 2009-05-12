"""
Utility package for use in ObitTalk

   ObitTalk is derived from the ParselTongue project at JIVE and
provides a scripting and interactive command line interface to
astronomical data and processing software.  In particular, AIPS and
FITS data structures as used in the AIPS and Obit software packages
are supported as well as AIPS tasks and Obit tasks and other python
enabled software.

This utility  package facilitates the access to data and images from
python as well as various interactive features.  The details of the
functions in this package are given later in this help.  Many of
these functions have equivalents in POPS although adapted to python.

   AIPS tasks will use the AIPS XAS TV which must be started separately
(Haven't checked that this works).  Obit tasks and ObitTalk use the
ObitView image display which must also be started independently.

   ObitTalk can start tasks either locally or on a remote machine
which has an ObitTalkServer process running.  Some data access is
supported through the AIPSUVData, AIPSImage, FITSUVData and FITSImage
classes.  Currently other python functions only work locally.

   Both Tasks and more detailed access to and manipulation of data
are available.  These are described breifly below and methods of
obtaining more detailed descriptions are described.

Obit Basics
-----------
   Obit consists of class libraries and a number of prepackaged
tasks similar to AIPS tasks.  The classes are implemented in c but
there are python bindings to much of the high-level functionality
allowing python scripts a high degree of flexibility in accessing
and manipulating data.
   Obit can support multiple physical data formats as long as they are
uniquely mapable to a common data model.  Above a data access level,
the underlying physical data representation is (mostly) hidden.
Currently, AIPS and FITS (as practiced by AIPS) are supported.
Only FITS format OTF data is supported.
   Control parameters to Obit routines are largely passed in an
InfoList structure (a type of associative array) but many of the
python interface routines take care of this detail and their
parameters are passed through a python dictionary.

Tasks:
-----------

AIPS tasks:
All AIPS tasks

ObitTasks:
AutoFlag  Radio interferometry data editing software
Calib     Calibrate visibility data (amp & phase)
CLCal     Apply gain solutions to a CL table
Convol    Convolve images
CubeClip  Blank insignificant pixels in an image
CubeVel   Make velocity image from spectral cube
FndSou    Make a source catalog from an image.
Feather   Task to Feather together images
GetJy     Determine calibrator flux densities
HGeom     Task to make an image consistent with another image
Imager    Radio interferometry imaging task
IonCal    Low frequency Field based calibration
IonImage  Low frequency Field based calibration and imaging
SCMap     Interferometry self calibration imaging
SetJy     Modify Source (SU) table
SNCor     Modify Visibility gain (SN) table
SNFilt    Fits for instrumental phases in SN table.
SNSmo     Smooth Visibility gain (SN) table
Squint    Imaging with beam squint correction
Squish    Compress image cube along third axis
SubImage  Task to copy a sub region of an image
TabCopy   Task to copy tables
Template  Task to print the mean, rms and extrema in an image
UVSub     Task to subtract a clean model from a uv data base
VL2VZ     Convert VL format source catalogs to VZ format
VLSSFix   Corrects residual geometry in low frequency images

Obit SD Tasks:
CCBCalib  Calibrate GBT CCB OTF format data
OTFImage  Image OTF format data

   To see task documentation either a python task object may first
be created is its documentation viewed or more directly:
AIPSHelp("AIPS_task_name")
or
ObitHelp("Obit_task_name")

   To create a task object:
>>> im=AIPSTask("IMEAN")
to create an AIPS task object for task IMEAN, or

>>> fe=ObitTask("Feather")
to create an Obit Task object for Feather
Note the names of the objects are arbitrary.

Task parameters can be set using the form object.parameter=value:
>>> im.inname="MY FILE"
where the parameter names are subject to minimum match.
Array values are given in square brackest "[  ]", the usual form
for a python list.  AIPS array values are indexed 1-relative and Obit
arrays 0-relative but this is largely transparent.
Note: unlike POPS, ALL strings are case sensitive.

Task parameters can be reviewed using the inputs() function:
>>> im.inputs()
or
>>> inputs(im)
Note: there is NO minimum match on functions and you must give
the parentheses.

POPS style help can be viewed:
>>> im.help()
or
>>> help(im)

or EXPLAIN (if available) by:
>>> im.explain()
or
>>> explain(im)

Tasks can be run using the go function:
>>> log=im.go()
The go function currently runs synchronously and does not return
until the task finishes.  The go function returns a list of the
task messages which in the above example are stored in the variable
log as a list of strings. This log may be saved to a disk file using
the saveLog function.

After a task is run which generates output values, these can be
viewed using the outputs function:
>> im.outputs()
and the values can be accessed through the task parameter.

The task functions work for both AIPS and Obit tasks.

ObitTalk DATA Classes
---------------------
   The ObitTalk classes  AIPSUVData, AIPSImage, FITSUVData and
FITSImage allow local or remote access to AIPS and FITS Images
and UV data.  Details of these classes interfaces can be viewed
using:
>>> help(AIPSUVData)
>>> help(AIPSImage)
>>> help(FITSUVData)
>>> help(FITSImage)

Obit classes and utility packages with python interfaces:
---------------------------------------------------------
   There are a number of Obit functions with high level python
interfaces.  To see more details import and view the help for each:

>>> import History
>>> help(History)
   

Obit/AIPS/Radio Interferometry/Image classes and utilities
AIPSDir        AIPS directory class
CArray         Complex array class
Catalog        Source Catalog package
CleanImage     Image CLEAN
CleanVis       Visibility based CLEAN
ConvUtil       Image convolution utilities
FArray         float array class
FArrayUtil     FArray utilities
FeatherUtil    Image feathering utilities
FFT            Fast Fourier Transform class
FInterpolate   Float array interpolator
FITSDir        FITS directory routines
History        History class
ImageDesc      Image Descriptor (header)
ImageMosaic    Image Mosaic class
Image          Image class
ImageUtil      Image utilities
InfoList       Obit associative array for control info
IonCal         Ionospheric calibration
MosaicUtil     Image mosaicing utilities
ODisplay       Interface to ObitView display
OPlot          pgplot interface
OErr           Obit message/error class
OSystem        Obit System class
OWindow        (CLEAN) image window class
SkyGeom        Celestial geometry
SkyModel       Sky model class
TableDesc      Table descriptor (header) class
TableList      Table list for data object (Image, UVData, OTF)
Table          Table class
TableUtil      Table utilities
TaskWindow     Task message window class
UVDesc         UV data descriptor (header)
UVGSolve       UV gain solutions
UVImager       UV data imager class
UV             UV data class
UVSelfCal      UV Self calibration class
UVSoln2Cal     UV SN to CL table routines.
ZernikeUtil    Zernike polynomial utilities

               Single dish/OTF imaging classes and utilities
CCBUtil        GBT CCB utility package
CleanOTF       Single dish (Hogbom) CLEAN
GBTDCROTF      Convert GBD DCR data to OTF format
OTFDesc        OTF Descriptor
OTFGetAtmCor   OTF Atmospheric correction utilities
OTFGetSoln     OTF calibration solution utilities
OTF            OTF ("On the Fly") data
OTFSoln2Cal    Utilities to convert OTF solutiion to calibration tables
OTFUtil        OTF Utilities

Obit Python Example

   Accessing and manipulating data through Obit requires an initialization
of the Obit System which is done when this (OTObit) module is loaded into
python so is not described here.
   The following example reads an AIPS image and writes a integerized FITS
image with the pixel values truncated at a set fraction of the RMS "noise"
in the image.  This operation creates an image which is more compressible
but with a controlled loss of precision.

# Specify input and output
inDisk   = 1
Aname    = "INPUT IMAGE"
Aclass   = "CLASS"
Aseq     = 1
outDisk  = 1
outFile  = "Quantized.fits"

# Create Images
inImage   = Image.newPAImage("Input image", Aname, Aclass, inDisk, Aseq, True, err)
# Note: inImage can also be created using getname(cno,disk)
outImage  = Image.newPFImage("Output image",   outFile,  outDisk,  False, err)
Image.PClone(inImage, outImage, err)   # Same structure etc.
OErr.printErrMsg(err, "Error initializing")

# Fraction of RMS if given
fract = 0.25

# Copy to quantized integer image with history
inHistory  = History.History("history", inImage.List, err)
Image.PCopyQuantizeFITS (inImage, outImage, err, fract=fract, inHistory=inHistory)
OErr.printErrMsg(err, "Writing to FITS")


Obit Python Table access example

Obit Table objects can be created as shown in the following:
inUV=UV.newPAUV("UV", "20050415", "LINE", 1, 1, True,err)
tabType="AIPS SU"
tabVer=1
access=UV.READONLY
su = inUV.NewTable(access,tabType,tabVer,err)

The table header (descriptor) can be obtained as a python Dict:
h= su.Desc.Dict

Data from a row in the table can be obtained as a python Dict:
su.Open(access,err)
row1=su.ReadRow(access,err)
OErr.printErrMsg(err, "Error reading")
su.Close(err)
print "row1",row1

Note: these dict are independent of the underlying data structures.

"""
# Interactive routines to Obit use from ObitTalk
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2005-2009
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
import Obit, Table, FArray, OErr, InfoList, History, OSystem
#import FITSDir, AIPSDir
import Image, ImageDesc, ImageUtil, TableList, ODisplay, UV, OWindow
import re
import TaskMsgBuffer
try:
    import TaskWindow
except:
    pass
else:
    pass
import os, AIPS, pickle, string, pydoc
import FITS
# ObitTalk classes
import AIPSData, FITSData, AIPSTask, ObitTask
from AIPSTask import AIPSTask
from ObitTask import ObitTask
# from OTObit import *
# ObitStart()

global Adisk, Fdisk, FITSdisks, nFITS, AIPSdisks, nAIPS
from FITSDir import FITSdisks, nFITS
from AIPSDir import AIPSdisks, nAIPS

# Init Obit
userno =  AIPS.AIPS.userno
popsno = 1
from ObitInit import *

# Display connection
disp = ODisplay.ODisplay("ObitView", "ObitView", err)

def ShowErr(err=err):
    """ Print any errors and clear stack
    
    err  = Python Obit Error/message stack, default of OTObit version
    """
    ################################################################
    OErr.printErrMsg(err, "Error")
    # end ShowErr

def ClearErr(err=err):
    """ Print any errors and clear stack
    
    err  = Python Obit Error/message stack, default is OTObit version
    """
    ################################################################
    OErr.printErrMsg(err, "Error")
    # end ClearErr

def Acat(disk=None, first=1, last=1000, Aname=None, Aclass=None, Aseq=0,
         giveList=False):
    """ Catalog listing of AIPS files on disk disk

    The class remembers the last disk accessed
    Strings use AIPS wild cards:
        blank => any
        '?'   => one of any character
        "*"   => arbitrary string
    If giveList then return list of CNOs
    disk      = AIPS disk number to list
    first     = lowest slot number to list
    last      = highest slot number to list
    Aname     = desired AIPS name, using AIPS wildcards, None -> don't check
    Aclass    = desired AIPS class, using AIPS wildcards, None -> don't check
    Aseq      = desired AIPS sequence, 0=> any
    giveList = If true, return list of CNOs matching
    """
    ################################################################
    global Adisk
    if disk==None:
        disk = Adisk
    else:
        Adisk = disk
    # Get catalog
    cat = AIPSData.AIPSCat(disk)    
    olist=AIPSDir.PListCat(cat.catalog, disk, first=first, last=last,
                           Aname=Aname, Aclass=Aclass, Aseq=Aseq,
                           giveList=giveList)
    OErr.printErrMsg(err, "Error with AIPS catalog")
    return olist
    # end Acat

def Fdir(disk=None, dir=None):
    """ Catalog listing of FITS files on disk disk

    The class remembers the last disk accessed
    disk      = AIPS disk number to list
    dir       = relative or abs. path of directory, def. = cwd
                Only used if disk == 0
    """
    ################################################################
    global Fdisk
    if disk==None:
        disk = Fdisk
    else:
        Fdisk = disk
    if dir==None:
        dir = "./"
    # Get catalog
    cat = FITSData.FITSCat(disk, dir=dir)
    if disk==0:
        catList = "Catalog for FITS disk "+dir+"\n"
    else:
        catList = "Catalog for FITS disk "+str(disk)+"\n"
    for entry in cat.catalog:
        catList = catList+("%d %s" % (entry[0],entry[1]))+"\n"
    pydoc.ttypager(catList)
    # end Fdir

def AMcat(disk=None, first=1, last=1000, Aname=None, Aclass=None, Aseq=0,
          giveList=False):
    """ Catalog listing of AIPS Image files on disk disk

    Strings use AIPS wild cards:
        blank => any
        '?'   => one of any character
        "*"   => arbitrary string
    If giveList then return list of CNOs
    disk      = AIPS disk number to list
    first     = lowest slot number to list
    last      = highest slot number to list
    Aname     = desired AIPS name, using AIPS wildcards, None -> don't check
    Aclass    = desired AIPS class, using AIPS wildcards, None -> don't check
    Aseq      = desired AIPS sequence, 0=> any
    giveList  = If true, return list of CNOs matching
    """
    ################################################################
    global Adisk
    if disk==None:
        disk = Adisk
    else:
        Adisk = disk
    # Get catalog
    cat = AIPSData.AIPSCat(disk)    
    olist=AIPSDir.PListCat(cat.catalog, disk, type="MA", first=first, last=last,
                           Aname=Aname, Aclass=Aclass, Aseq=Aseq,
                           giveList=giveList)
    OErr.printErrMsg(err, "Error with AIPS catalog")
    return olist
    # end AMcat

def AUcat(disk=None, first=1, last=1000, Aname=None, Aclass=None, Aseq=0,
          giveList=False):
    """ Catalog listing of AIPS UV data files on disk disk

    Strings use AIPS wild cards:
        blank => any
        '?'   => one of any character
        "*"   => arbitrary string
    If giveList then return list of CNOs
    disk      = AIPS disk number to list
    first     = lowest slot number to list
    last      = highest slot number to list
    Aname     = desired AIPS name, using AIPS wildcards, None -> don't check
    Aclass    = desired AIPS class, using AIPS wildcards, None -> don't check
    Aseq      = desired AIPS sequence, 0=> any
    giveList  = If true, return list of CNOs matching
    """
    ################################################################
    global Adisk
    if disk==None:
        disk = Adisk
    else:
        Adisk = disk
    # Get catalog
    cat = AIPSData.AIPSCat(disk)    
    olist = AIPSDir.PListCat(cat.catalog, disk, type="UV", first=first, last=last,
                             Aname=Aname, Aclass=Aclass, Aseq=Aseq,
                             giveList=giveList)
    OErr.printErrMsg(err, "Error with AIPS catalog")
    return olist
    # end AUcat

def AllDest (disk=None, Atype="  ", Aname="            ", Aclass="      ", Aseq=0):
    """ Delete AIPS files matching a pattern

    Strings use AIPS wild cards:
        blank => any
        '?'   => one of any character
        "*"   => arbitrary string
    disk      = AIPS disk number, 0=>all
    Atype     = AIPS entry type, 'MA' or 'UV'; '  => all
    Aname     = desired AIPS name, using AIPS wildcards, None -> don't check
    Aclass    = desired AIPS class, using AIPS wildcards, None -> don't check
    Aseq      = desired AIPS sequence, 0=> any
    """
    ################################################################
    global Adisk
    if disk==None:
        disk = Adisk
    else:
        if disk>0:
            Adisk = disk
    # Confirm
    prompt = "Do you really want to delete all AIPS "+Atype+" files on disk(s) "\
             +str(disk)+"\nwith names matching "+Aname+"."+Aclass+"."+str(Aseq)+ \
             " y/n "
    ans = raw_input(prompt)
    if ans.startswith('y'):
        AIPSDir.PAllDest (disk, err, Atype=Atype, Aname=Aname, Aclass=Aclass,
                          Aseq=Aseq)
    else:
        print "Not confirmed"
    OErr.printErrMsg(err, "Error with destroying AIPS enreies")
    # end AllDest

def getname(cno, disk=None):
    """ Return Obit object for AIPS file in cno on disk

    cno       = AIPS catalog slot number 
    disk      = AIPS disk number
    """
    ################################################################
    global Adisk
    if disk==None:
        disk = Adisk
    else:
        Adisk = disk
    user =  AIPS.AIPS.userno
    #s = AIPSDir.PInfo(disk, user, cno, err)
    #OErr.printErrMsg(err, "Error with AIPS catalog")
    cat = AIPSData.AIPSCat(disk)
    s = cat.info(cno)
    if not s:   # Not found
        raise RuntimeError,"Cannot find slot "+str(cno)+",disk "+str(disk)+",user "+str(user)
    # parse returned string
    Aname = s[0:12]
    Aclass = s[13:19]
    Aseq = int(s[20:25])
    Atype = s[26:28]
    if Atype == 'MA':
        # Create AIPSData or ObitTalk if remote/local
        if AIPS.AIPS.disks[disk].url:
            out = AIPSData.AIPSImage(Aname, Aclass, disk, Aseq)
            out.Otype  = 'Image'
            out.isOK   = True
        else:
            out = Image.newPAImage("AIPS image", Aname, Aclass, disk, Aseq, True, err, \
                                   verbose=False)
        print "AIPS Image",Aname, Aclass, disk, Aseq
    elif Atype == 'UV':
         # Create AIPSData or ObitTalk if remote/local
        if AIPS.AIPS.disks[disk].url:
            out = AIPSData.AIPSUVData(Aname, Aclass, disk, Aseq)
            out.Otype  = 'UV'
            out.isOK   = True
        else:
            out = UV.newPAUV("AIPS UV data", Aname, Aclass, disk, Aseq, True, err, \
                             verbose=False)
        print "AIPS UV",Aname, Aclass, disk, Aseq
    out.Aname  = Aname
    out.Aclass = Aclass
    out.Aseq   = Aseq 
    out.Atype  = Atype
    out.Disk   = disk
    out.Acno   = cno
    return out
    # end getname

def getFITS(file, disk=None, Ftype='Image'):
    """ Return Obit object for FITS file in file on disk

    file      = FITS file name
    disk      = FITS disk number
    Ftype     = FITS data type: 'Image', 'UV'
    """
    ################################################################
    global Fdisk
    if disk==None:
        disk = Fdisk
    if Ftype == 'Image':
        # ObitTalk (local) or AIPSData (remote)?
        if (disk>0) and FITS.FITS.disks[disk].url:
            out = FITSData.FITSImage(file, disk)
        else:
            out = Image.newPFImage("FITS image", file, disk, True, err)
    elif Ftype == 'UV':
        if (disk>0) and FITS.FITS.disks[disk].url:
            out = FITSData.FITSUVData(file, disk)
        else:
            out = UV.newPFUV("FITS UV data", file, disk, True, err)
    out.Fname  = file
    out.Disk   = disk 
    out.Otype  = Ftype
    return out
    # end getFITS

def tvlod(image, window=None):
    """ display image

    image  = Obit Image or AIPSImage, created with getname, getFITS
    window = Optional window for image to edit
    """
    ################################################################
    if image.__class__==AIPSData.AIPSImage:
        # AIPS Image
        image.display(disp.url)
        return
    # Otherwise must be Obit Image
    if Image.PIsA(image):
        # Obit/Image
        ODisplay.PImage(disp, image, err, window=window)
    elif image.__class__==AIPSData.AIPSImage:
        # AIPS Image
        image.display(disp.url)
    elif image.__class__==FITSData.FITSImage:
        # FITS Image - only works locally
        tmp = Image.newPFImage("FITS Image",image.filename, image.disk, True, err)
        ODisplay.PImage(disp, tmp, err, window=window)
        del tmp
    OErr.printErr(err)  # In case of failure
    # end tvlod

def newDisplay(port=8765, URL=None):
    """ Recreate display to another display server

    port   = port number on local machine
    URL    = Full URL (e.g. http://localhost:8765/RPC2)
    """
    ################################################################
    global disp
    del disp
    if URL:
        disp = ODisplay.ODisplay("ObitView", URL, err)
    else:
        url = "http://localhost:"+str(port)+"/RPC2"
        disp = ODisplay.ODisplay("ObitView", url, err)
    # end newDisplay

def window (image):
    """ Make a window object for an image

    Returns OWindow object
    image  = Obit image object
    """
    ################################################################
    if Image.PIsA(image):
        # Obit/Image
        naxis = image.Desc.Dict["inaxes"][0:2]
    elif image.__class__==AIPSData.AIPSImage:
        # AIPS Image
        tmp = Image.newPAImage("AIPS Image",image.name, image.klass, image.disk, \
                               image.seq, True, err)
        naxis = tmp.Desc.Dict["inaxes"][0:2]
        del tmp
    elif image.__class__==FITSData.FITSImage:
        # FITS Image
        tmp = Image.newPFImage("FITS Image",image.filename, image.disk, True, err)
        naxis = tmp.Desc.Dict["inaxes"][0:2]
        del tmp
    return OWindow.PCreate1("Window", naxis, err)
    # end window

def go (TaskObj, MsgBuf=False, URL="http://localhost:8777/RPC2"):
    """ Execute task or script

    Returns TaskWindow object if run asynchronously (doWait=True)
    or the task message log if run synchronously (doWait=False)
    The wait() function on the TaskWindow or TaskMsgBuffer will hang
    until the task finishes.
    Returns either a TaskWindow or TaskMsgBuffer depending on MsgBuf
    TaskObj    = Task object to execute
                 If doWait member is true run synchronously,
                 else run with messages in a separate Message window
    MsgBuf     = if true and  TaskObj.doWait=False run asynchronously
                 using a TaskMsgBuffer
    URL        = URL of ObitMess message server if MsgBuf=False
    """
    ################################################################
    import OTWindow

    if TaskObj.doWait:
        return TaskObj.go()
    else:
        if MsgBuf:
            # use TaskMsgBuffer
            tb = TaskMsgBuffer.TaskMsgBuffer(TaskObj)
            tb.start()
            return tb
        else:
            # use TaskWindow
            tw = TaskWindow.TaskWindow(TaskObj,URL=URL)
            #OTWindow.newMsgWin(tw)
            tw.start()
            return tw
    # end go
   
def inputs (TaskObj):
    """ List task inputs

    TaskObj    = Task object whose inputs to list
    """
    ################################################################
    TaskObj.inputs()
    # end inputs
   
def explain (TaskObj):
    """ Give explanation for a task if available

    TaskObj    = Task object whose inputs to list
    """
    ################################################################
    TaskObj.explain()
    # end explain
   
def AIPSHelp (Task):
    """ Give Help for AIPS task Task

    Task    = AIPSTask name to give (e.g. "imean")
    """
    ################################################################
    t=AIPSTask(Task)
    t.help()
    # end  AIPSHelp
   
def ObitHelp (Task):
    """ Give Help for OBIT task Task

    Task    = ObitTask name to give (e.g. "Feather")
    """
    ################################################################
    t=ObitTask(Task)
    t.help()
    # end  ObitHelp
   
def copyInputs (inTask, outTask):
    """ Copy values from one task object to another

    Copies parameter values from inTask to outTask which are in both the
    inTask and outTask _input_list.
    Need not be the same task.
    inTask    = Task object to copy from
    outTask   = Task object to copy to
    """
    ################################################################
    for x in inTask._input_list:
        if x in outTask._input_list:
            y = inTask.__dict__[x]
            outTask.__dict__[x] = y
    
    # end copyInputs

def imhead (ObitObj):
    """ List header

    ObitObj    = Obit or ObitTalk data object
    """
    ################################################################
    if ObitObj.__class__==AIPSData.AIPSImage:
        # AIPS Image
        Image.PHeader(ObitObj, err)
    elif ObitObj.__class__==FITSData.FITSImage:
        # FITS Image
        Image.PHeader(ObitObj, err)
    elif ObitObj.__class__==AIPSData.AIPSUVData:
        # AIPS UVData
        UV.PHeader(ObitObj, err)
    elif ObitObj.__class__==FITSData.FITSUVData:
        # FITS UVData
        UV.PHeader(ObitObj, err)
    else:
        # Presume it's an Obit object
        ObitObj.Header(err)
    # end imhead
   
def setname (inn, out):
    """ Copy file definition from inn to out as in...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    inn  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask:
        out.DataType = inn.FileType
        out.inDisk   = int(inn.Disk)
        if inn.FileType == 'FITS':
            out.inFile = inn.Fname
        else:   # AIPS
            out.inName  = inn.Aname
            out.inClass = inn.Aclass
            out.inSeq   = int(inn.Aseq)
    else:  # AIPS
        out.inname  = inn.Aname
        out.inclass = inn.Aclass
        out.inseq   = inn.Aseq
        out.indisk  = inn.Disk
    # end setname
   
def set2name (in2, out):
    """ Copy file definition from in2 to out as in2...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    in2  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask:
        out.DataType  = in2.FileType
        out.in2Disk   = int(in2.Disk)
        if in2.FileType == 'FITS':
            out.in2File = in2.Fname
        else:   # AIPS
            out.in2Name  = in2.Aname
            out.in2Class = in2.Aclass
            out.in2Seq   = int(in2.Aseq)
    else: # AIPS
        out.in2name  = in2.Aname
        out.in2class = in2.Aclass
        out.in2seq   = in2.Aseq
        out.in2disk  = in2.Disk
    # end set2name
   
def set3name (in3, out):
    """ Copy file definition from in3 to out as in3...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    in3  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask:
        out.DataType  = in3.FileType
        out.in3Disk   = int(in3.Disk)
        if in3.FileType == 'FITS':
            out.in3File = in3.Fname
        else:   # AIPS
            out.in3Name  = in3.Aname
            out.in3Class = in3.Aclass
            out.in3Seq   = int(in3.Aseq)
    else: # AIPS
        out.in3name  = in3.Aname
        out.in3class = in3.Aclass
        out.in3seq   = in3.Aseq
        out.in3disk  = in3.Disk
    # end set3name
   
def set4name (in4, out):
    """ Copy file definition from in4 to out as in4...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    in4  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask:
        out.DataType  = in4.FileType
        out.in4Disk   = int(in4.Disk)
        if in4.FileType == 'FITS':
            out.in4File = in4.Fname
        else:   # AIPS
            out.in4Name  = in4.Aname
            out.in4Class = in4.Aclass
            out.in4Seq   = int(in4.Aseq)
    else: # AIPS
        out.in4name  = in4.Aname
        out.in4class = in4.Aclass
        out.in4seq   = in4.Aseq
        out.in4disk  = in4.Disk
    # end set4name
   
def setoname (inn, out):
    """ Copy file definition from inn to out as outdisk...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    inn  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    # AIPS or Obit?
    if out.__class__ == ObitTask:
        out.DataType  = inn.FileType
        out.outDisk   = int(inn.Disk)
        if inn.FileType == 'FITS':
            out.outFile = inn.Fname
        else:   # AIPS
            out.outName  = inn.Aname
            out.outClass = inn.Aclass
            out.outSeq   = int(inn.Aseq)
    else:  # AIPS
        out.outname  = inn.Aname
        out.outclass = inn.Aclass
        out.outseq   = inn.Aseq
        out.outdisk  = inn.Disk
    # end setoname
   
def setwindow (w, out):
    """ Set BLC and TRC members on out from OWindow w

    Uses first window in first field on w which must be a rectangle
    This may be set interactively using tvlod
    w    = OWindow object
    out  = ObitTask object, BLC and TRC members [0] and [1] are modified
    """
    ################################################################
    # Must be rectangle
    l =  OWindow.PGetList(w, 1, err)
    if l[0][1] !=0:
        raise TypeError,"Window MUST be a rectangle"
    # AIPS or Obit?
    if out.__class__ == ObitTask:
        out.BLC[0] = l[0][2]+1  # make 1-rel
        out.BLC[1] = l[0][3]+1  
        out.TRC[0] = l[0][4]+1  
        out.TRC[1] = l[0][5]+1
    else:  # AIPS
        out.blc[1] = l[0][2]+1  # make 1-rel
        out.blc[2] = l[0][3]+1  
        out.trc[1] = l[0][4]+1  
        out.trc[2] = l[0][5]+1
       
    # end setwindow 
   
def zap (o):
    """ Zap object o

    Removes all external components (files)
    o    = Obit Data object to delete
    """
    ################################################################
    if o.__class__==AIPSData.AIPSImage:
        # AIPS Image
        o.zap()
    elif o.__class__==AIPSData.AIPSUVData:
        # AIPS UVData
        o.zap()
    else:
        # Presume it's an Obit object
        o.Zap(err)
        ShowErr()
    # end zap
   
def clearstat (o, code=4):
    """ Clear status of AIPS catalog entry

    Clears AIPS status of object o,
    Optionally sets status using code parameter
    o    = Obit AIPS Data object
    code = status code:
        0 = Add write status
        1 = Clear write status
        2 = Increment Read Status
        3 = Decrement Read Status
        4 = Clear All Status
        
    """
    ################################################################
    # Get APS info
    adisk=None
    user =  AIPS.AIPS.userno
    if o.__class__==AIPSData.AIPSImage:
        # AIPS Image
        o.clearstat(code)
    elif o.__class__==AIPSData.AIPSUVData:
        # AIPS UVData
        o.clearstat(code)
    else:
        # Presume it's an Obit object
        if o.FileType == "AIPS":
            cno = AIPSDir.PTestCNO (o.Disk, user, o.Aname, o.Aclass, o.Atype, o.Aseq, err)
            if cno>0:
                adisk = o.Disk
    
    # Clear status - ignore if not AIPS data (adisk==None)
    if adisk:
        e = Obit.AIPSDirStatus(adisk, user, cno, code, err.me)
    ShowErr()
    # end clearstat
   
def alldest(Aname=".*", Aclass=".*", Atype=".?", Adisk=0, Aseq=0, test=False):
    """ Delete AIPS files matching a pattern

    Uses regular expression matching for strings
    Note: "+" values are escaped
    Clears any status before deleting
    Aname    = AIPS file name , " " => any
    Aclass   = AIPS class name,  " " => any
    Atype    = 'MA', 'UV' or any
    Adisk    = AIPS disk number, 0=> any
    Aseq     = AIPS sequence number; 0=> any
    test    = if true only list and not delete
    """
    ################################################################
    # confirm
    prompt = "Do you really want to delets all AIPS files matching\n" +\
             "name:"+Aname+" class:"+Aclass+" disk:"+str(Adisk)+ \
             " seq:"+str(Aseq)+ " ?\n" + \
             "yes/no: "
    ans = raw_input(prompt)
    if ans != "yes":
        print "confirmination negative"
        return
    # How many disks
    # "+" is a special symbol in regular expressions, escape any
    Aname = re.sub("\+", "\\+", Aname)
    Aclass = re.sub("\+", "\\+", Aclass)
    ndisk = len(AIPSdisks)
    if Adisk>0:
        disk1 = Adisk; disk2 = Adisk
    else:
        disk1 = 1; disk2 = ndisk
    user =  AIPS.AIPS.userno  # AIPS user no
    for idisk  in range (disk1, disk2+1):
        # loop over directory 
        maxcno = AIPSDir.PNumber (idisk, user, err)
        for cno in range(1, maxcno+1):
            line=AIPSDir.PInfo(idisk, user, cno, err);
            if err.isErr:
                break
            obj = None
            if (line!=None):
                tname  = line[0:12]
                tclass = line[13:19]
                tseq   = int(line[20:25])
                ttype  = line[26:28]
                # Check type
                mat = re.match(Atype, ttype) or \
                      ((Atype!="MA") and (Atype!="UV"))
                # Check seq
                mat = mat and ((Aseq<=0) or (tseq==Aseq))
                # Check Name
                mat = mat and re.match(Aname,tname)
                # Check Class
                mat = mat and re.match(Aclass,tclass)
                if mat:   # Found a match?
                    print cno,tname, tclass, tseq, ttype # debug
                    #obj = getname(cno, disk=idisk)
                    obj = Image.newPAImage("zap", tname, tclass, idisk, tseq,
                                           True, err, verbose=False)
                    if not test:
                        obj.Atype = ttype
                        obj.Aseq  = tseq
                        clearstat(obj)
                        obj.Zap(err)
        # end loop over cno
        if err.isErr:
            OErr.printErrMsg(err, "Error getting AIPS directory listing")
    # end loop over disk
    # end alldest

def tput (to, file=None):
    """ save task object

    save values in task object
    to    = task object to save
    file  = optional file name, the default is <task_name>.pickle
            in the current working directory
    """
    ################################################################
    saveit = {}
    for adverb in to._input_list:
        value = to._retype(to.__dict__[adverb])
        saveit[adverb] = value

    # Save type
    saveit["TaskType"] = to.__class__
    
    tfile = file
    if file==None:
        tfile = to._name+".pickle"
    fd = open(tfile, "w")
    pickle.dump(saveit, fd)
    fd.close()
    # end tput
   
def tget (inn, file=None):
    """ Restore task object from disk

    Restore values in task object
    inn   = task name, or a task object of the desired type
            in the latter case, the input object will NOT be modified
    file  = optional file name, the default is <task_name>.pickle
            in the current working directory
    """
    ################################################################
    # Get name
    if type(inn)==str:
        name = inn
    else:
        # it better be a task object
        name = inn._name

    # unpickle file
    tfile = file
    if file==None:
        tfile = name+".pickle"
    fd = open(tfile, "r")
    saveit = pickle.load(fd)
    fd.close()

    myType = saveit["TaskType"]
    if myType == AIPSTask:
        to = AIPSTask(name)
    else:
        to = ObitTask(name)
    for adverb in saveit:
        if adverb == "TaskType":
            pass
        else:
            to.__dict__[adverb] = saveit[adverb]

    return to
    # end tget
   
def saveLog (log, file=None):
    """ Save AIPS or Obit task log to file

    Save log messages
    log   = log string array returned from the execution of an
            Obit or AIPS task.
    file  = optional file name, the default is "task.log"
            in the current working directory.
            New entries are appended.
    """
    ################################################################
    # log must be list
    if type(log) != list:
        raise TypeError,"log MUST be a list"
    # save
    tfile = file
    if file==None:
        tfile = "task.log"
    fd = open(tfile, "a")
    for msg in log:
        if type(msg)==str:
            fd.write(msg+"\n")
        else:
            fd.write(msg[1]+"\n")
    fd.close()
    # end saveLog

def PrintHistory (ObitObj, hiStart=1, hiEnd=1000000, file=None):
    """ Display history log or write to file

    Reads selected history records and displays with "more"
    ObitObj   = Python Obit object with history
    err       = Python Obit Error/message stack
    hiStart   = if given the first (1-rel) history record
    hiEnd     = if given the highest (1-rel) history record
    file      = if present, the name of a file into which to write
                the history rather than displaying it on the screen
    """
    ################################################################
    # Figure out label for data
    label = "unknown"
    if ObitObj.__class__==AIPSData.AIPSImage:
        # AIPS Image
        label = "AIPS Image:"+ObitObj.name,+"."+ObitObj.klass+"."+\
                str(ObitObj.seq)+"."+str(ObitObj.disk)
        tmp = Image.newPAImage("AIPS Image",ObitObj.name, ObitObj.klass, ObitObj.disk, \
              ObitObj.seq, True, err)
    elif ObitObj.__class__==FITSData.FITSImage:
        # FITS Image
        label = "FITS Image:"+ObitObj.filename,+"."+str(ObitObj.disk)
        tmp = Image.newPFImage("FITS Image",ObitObj.filename, ObitObj.disk, True, err)
    elif ObitObj.__class__==AIPSData.AIPSUVData:
        # AIPS UVData
        label = "AIPS UV:"+ObitObj.name,+"."+ObitObj.klass+"."+\
                str(ObitObj.seq)+"."+str(ObitObj.disk)
        tmp = UV.newPAImage("AIPS UVData",ObitObj.name, ObitObj.klass, ObitObj.disk, \
              ObitObj.seq, True, err)
    elif ObitObj.__class__==FITSData.FITSUVData:
        # FITS UVData
        label = "FITS UV:"+ObitObj.filename,+"."+str(ObitObj.disk)
        tmp = UV.newPFImage("FITS UVData",ObitObj.filename, ObitObj.disk, True, err)
    else:
        # Presume it's an Obit object
        tmp = None
        if ObitObj.FileType == "AIPS":
            label = ObitObj.FileType+":"+ObitObj.Otype+":"+ObitObj.Aname+"."+ObitObj.Aclass+ \
                    "."+str(ObitObj.Aseq)+"."+str(ObitObj.Disk)
        if (ObitObj.FileType=="FITS"):
            label = ObitObj.FileType+":"+ObitObj.Otype+":"+ObitObj.FileName+"."+str(ObitObj.Disk)
    
            
    # Create history
    if tmp:
        hist = History.History("history", tmp.List, err)
    else:
        hist = History.History("history", ObitObj.List, err)
        
    # read and swallow history
    r=hist.Open(History.READONLY,err)
    recno=hiStart-1
    hilist = "History for"+label+"\n"
    x = hist.ReadRec(recno,err)
    while (len(x)>0) and (recno<hiEnd):
        hilist = hilist+string.rjust(str(recno+1),6)+" "+x+"\n"
        recno = recno+1
        x = hist.ReadRec(recno,err)
    #print x
    
    r=hist.Close(err)

    # Display or log
    if file:
        fd = open(file, "a")
        fd.write(hilist)
        fd.close()
    else:
        pydoc.ttypager(hilist)
    del hilist
    # end PrintHistory

def tvstat (inImage):
    """ Set region in an image using the display and tell mean, rms

    Returns dictionary with statistics of selected region with entries:
        Mean    = Mean value
        RMSHist = RMS value from a histogram analysis
        RMS     = Simple RMS value
        Max     = maximum value
        MaxPos  = pixel of maximum value
        Min     = minimum value
        MinPos  = pixel of minimum value
        Flux    = Flux density if CLEAN beam given, else -1
        BeamArea= CLEAN Beam area in pixels
    inImage   = Python Image object, created with getname, getFITS
    """
    ################################################################
    # Get window
    w = window(inImage)
    tvlod(inImage,w)
    ShowErr()
    # Must be rectangle
    l =  OWindow.PGetList(w, 1, err)
    if l[0][1] !=0:
        raise TypeError,"Window MUST be a single rectangle"
    blc = [min(l[0][2]+1, l[0][4]+1), min(l[0][3]+1, l[0][5]+1)]
    trc = [max(l[0][2]+1, l[0][4]+1), max(l[0][3]+1, l[0][5]+1)]
    
    # Read plane
    p    = Image.PReadPlane(inImage,err,blc=blc,trc=trc)
    ShowErr()
    head = inImage.Desc.Dict  # Header

    # Get statistics
    Mean = p.Mean
    RMS  = p.RMS
    RawRMS  = p.RawRMS
    MaxPos=[0,0]
    Max = FArray.PMax(p, MaxPos)
    MaxPos[0] = MaxPos[0]+blc[0]
    MaxPos[1] = MaxPos[1]+blc[1]
    MinPos=[0,0]
    Min = FArray.PMin(p, MinPos)
    MinPos[0] = MinPos[0]+blc[0]
    MinPos[1] = MinPos[1]+blc[1]
    # Integrated flux density
    Flux = -1.0
    if (head["beamMaj"]>0.0) :
        beamarea = 1.1331*(head["beamMaj"]/abs(head["cdelt"][0])) * \
                   (head["beamMin"]/abs(head["cdelt"][1]))
        Flux = p.Sum/beamarea
    print "Region Mean %g, RMSHist %g RMS %g" % (Mean, RMS, RawRMS)
    print "  Max %g @ pixel " % Max, MaxPos
    print "  Min %g @ pixel " % Min, MinPos
    print "  Integrated Flux density %g, beam area = %7.1f pixels" % (Flux, beamarea)
    
    # Reset BLC, TRC
    blc = [1,1,1,1,1]
    trc = [0,0,0,0,0]
    Image.POpen(inImage, Image.READONLY, err, blc=blc, trc=trc)
    Image.PClose (inImage, err)
    ShowErr()
    
    del w, p, blc, trc
    return {"Mean":Mean,"RMSHist":RMS,"RMS":RawRMS,"Max":Max, \
            "MaxPos":MaxPos,"Min":Min,"MinPos":MinPos,"Flux":Flux,
            "BeamArea":beamarea}
    # end tvstat
   

def imstat (inImage, blc=[1,1,1,1,1], trc=[0,0,0,0,0]):
    """ Get statistics in a specified region of an image plane

    Returns dictionary with statistics of selected region with entries:
        Mean    = Mean value
        RMSHist = RMS value from a histogram analysis
        RMS     = Simple RMS value
        Max     = maximum value
        MaxPos  = pixel of maximum value
        Min     = minimum value
        MinPos  = pixel of minimum value
        Flux    = Flux density if CLEAN beam given, else -1
        BeamArea= CLEAN Beam area in pixels
    inImage   = Python Image object, created with getname, getFITS
    """
    ################################################################
    # Read plane
    p    = Image.PReadPlane(inImage,err,blc=blc,trc=trc)
    ShowErr()
    head = inImage.Desc.Dict  # Header

    # Get statistics
    Mean = p.Mean
    RMS  = p.RMS
    RawRMS  = p.RawRMS
    MaxPos=[0,0]
    Max = FArray.PMax(p, MaxPos)
    MaxPos[0] = MaxPos[0]+blc[0]
    MaxPos[1] = MaxPos[1]+blc[1]
    MinPos=[0,0]
    Min = FArray.PMin(p, MinPos)
    MinPos[0] = MinPos[0]+blc[0]
    MinPos[1] = MinPos[1]+blc[1]
    # Integrated flux density
    Flux = -1.0
    if (head["beamMaj"]>0.0) :
        beamarea = 1.1331*(head["beamMaj"]/abs(head["cdelt"][0])) * \
                   (head["beamMin"]/abs(head["cdelt"][1]))
        Flux = p.Sum/beamarea
    print "Region Mean %g, RMSHist %g RMS %g" % (Mean, RMS, RawRMS)
    print "  Max %g @ pixel " % Max, MaxPos
    print "  Min %g @ pixel " % Min, MinPos
    print "  Integrated Flux density %g, beam area = %7.1f pixels" % (Flux, beamarea)
   
    # Reset BLC, TRC
    blc = [1,1,1,1,1]
    trc = [0,0,0,0,0]
    Image.POpen(inImage, Image.READONLY, err, blc=blc, trc=trc)
    Image.PClose (inImage, err)
    ShowErr()
    
    del p, blc, trc
    return {"Mean":Mean,"RMSHist":RMS,"RMS":RawRMS,"Max":Max, \
            "MaxPos":MaxPos,"Min":Min,"MinPos":MinPos,"Flux":Flux,
            "BeamArea":beamarea}
    # end imstat
   

def altswitch(inImage):
    """ Switch frequency and velocity

    Algorithm lifted from AIPS AU7.FOR
    inImage   = Python Image object, created with getname, getFITS
    """
    ################################################################
    # Get Header dictionary
    hd = inImage.Desc.Dict 
    velite = 2.997924562e8   # Speed of light
    
    i = hd["jlocf"]
    frqType = hd["ctype"][i]
    if frqType[0:4]=="FREQ":
        # Has Frequency -  convert to velocity
        if  hd["VelDef"]==1:
            vsign = -1.0
            ctype = "VELO"
        else:
            vsign = 1.0
            ctype = "FELO"

        if hd["VelReference"]==1:
            hd["ctype"][i] = ctype+"-LSR"
        elif hd["VelReference"]==2:
            hd["ctype"][i] = ctype+"-HEL"
        elif hd["VelReference"]==3:
            hd["ctype"][i] = ctype+"-OBS"
        else:
            hd["ctype"][i] = ctype+"-LSR"
            hd["VelReference"] = 1  # Fix?
        tCrpix = hd["crpix"][i]
        tCrval = hd["crval"][i]
        hd["crpix"][i] = hd["altCrpix"]
        hd["crval"][i] = hd["altRef"]
        delnu   = hd["cdelt"][i]
        refnu   = tCrval
        frline  = hd["altCrpix"]
        dvzero  = hd["altRef"]
        hd["cdelt"][i] = -delnu * (velite + vsign * dvzero) / \
                         (refnu + delnu * (frline - tCrpix))
        hd["altCrpix"] = tCrpix
        hd["altRef"]   = tCrval
        
    # Look for velocity
    if (frqType[0:4]=="VELO") or (frqType[0:4]=="FELO"):
        # Has Velocity convert to frequency
        if frqType[0:4]=="FELO":
            vsign = 1.0
        else:
            vsign = -1.0
        hd["ctype"][i] = "FREQ    "
        tCdelt = hd["cdelt"][i]
        tCrpix = hd["crpix"][i]
        tCrval = hd["crval"][i]
        delnu  = -tCdelt * hd["altRef"] / \
                (velite + vsign * tCrval + tCdelt * (tCrpix - hd["altCrpix"]))
        hd["crpix"][i] = hd["altCrpix"]
        hd["crval"][i] = hd["altRef"]
        hd["altCrpix"] = tCrpix
        hd["altRef"]   = tCrval
        hd["cdelt"][i] = delnu

    inImage.Desc.Dict = hd     # Update header
    inImage.UpdateDesc(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error updating Image descriptor")
    return
    # end altswitch

# Functions to copy between AIPS and FITS
def uvlod(filename, inDisk, Aname, Aclass, Adisk, Aseq, err):
    """ Load FITS UV data to AIPS

    Read a UVTAB FITS UV data file and write an AIPS data set
    filename   = name of FITS file
    inDisk     = FITS directory number
    Aname      = AIPS name of file
    Aclass     = AIPS class of file
    Aseq       = AIPS sequence number of file, 0=> create new
    Adisk      = FITS directory number
    err        = Python Obit Error/message stack
    returns AIPS UV data object
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Get input
    inUV = UV.newPFUV("FITS UV DATA", filename, inDisk, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with FITS data")
    # Get output, create new if seq=0
    if Aseq<1:
        OErr.printErr(err)   # Print any outstanding messages
        user = OSystem.PGetAIPSuser()
        Aseq=AIPSDir.PHiSeq(Adisk,user,Aname,Aclass,"MA",err)
        # If it already exists, increment seq
        if AIPSDir.PTestCNO(Adisk,user,Aname,Aclass,"MA",Aseq,err)>0:
            Aseq = Aseq+1
        OErr.PClear(err)     # Clear any message/error
    print "Creating AIPS UV file",Aname,".",Aclass,".",Aseq,"on disk",Adisk
    outUV = UV.newPAUV("AIPS UV DATA", Aname, Aclass, Adisk, Aseq, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")
    # Copy
    UV.PCopy (inUV, outUV, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying UV data to AIPS")
    # Copy History
    inHistory  = History.History("inhistory",  inUV.List, err)
    outHistory = History.History("outhistory", outUV.List, err)
    History.PCopyHeader(inHistory, outHistory, err)
    # Add history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit uvlod",err)
    outHistory.WriteRec(-1,"uvlod   / FITS file "+filename+" disk "+str(inDisk),err)
    outHistory.Close(err)
   #
    # Copy Tables
    exclude=["AIPS HI", "AIPS AN", "AIPS FQ", "AIPS SL", "AIPS PL", "History"]
    include=[]
    UV.PCopyTables (inUV, outUV, exclude, include, err)
    return outUV  # return new object
    # end uvlod

def uvtab(inUV, filename, outDisk, err, compress=False, \
          exclude=["AIPS HI", "AIPS AN", "AIPS FQ", "AIPS SL", "AIPS PL"], \
          include=[], headHi=False):
    """ Write UV data as FITS file

    Write a UV data set as a FITAB format file
    History written to header
    inUV       = UV data to copy
    filename   = name of FITS file
    inDisk     = FITS directory number
    err        = Python Obit Error/message stack
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy (FQ, AN always done )
                 Exclude has presidence over include
    headHi     = if True move history to header, else leave in History table
    returns FITS UV data object
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Set output
    outUV = UV.newPFUV("FITS UV DATA", filename, outDisk, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating FITS data")
    #Compressed?
    if compress:
        inInfo = UV.PGetList(outUV)    # 
        dim = [1,1,1,1,1]
        InfoList.PAlwaysPutBoolean (inInfo, "Compress", dim, [True])        
    # Copy
    UV.PCopy (inUV, outUV, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying UV data to FITS")
    # History
    inHistory  = History.History("inhistory",  outUV.List, err)
    outHistory = History.History("outhistory", outUV.List, err)
    # Add history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit uvtab",err)
    outHistory.WriteRec(-1,"uvtab   / FITS file "+filename+" disk "+str(outDisk),err)
    outHistory.Close(err)
    # History in header?
    if headHi:
        History.PCopy2Header (inHistory, outHistory, err)
        OErr.printErrMsg(err, "Error with history")
        # zap table
        outHistory.Zap(err)
    # Copy Tables
    UV.PCopyTables (inUV, outUV, exclude, include, err)
    return outUV  # return new object
    # end uvtab

def uvTabSave(inUV, filename, outDisk, err, \
          exclude=["AIPS HI", "AIPS_AN", "AIPS FQ","AIPS PL", "AIPS SL"], include=[]):
    """ Write UV data tables (but not data) to a FITS file

    Write tables associated with UV data set as a FITAB format file
    History written to header
    inUV       = UV data to copy
    filename   = name of FITS file
    inDisk     = FITS directory number
    err        = Python Obit Error/message stack
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy (FQ, AN always done )
    returns FITS UV data object
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Set output
    outUV = UV.newPFUV("FITS UV DATA", filename, outDisk, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating FITS data")
    # Clone
    UV.PClone (inUV, outUV, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error cloning UV data to FITS")
    # Copy back to header
    inHistory  = History.History("inhistory",  outUV.List, err)
    outHistory = History.History("outhistory", outUV.List, err)
    # Add history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit uvTabSave",err)
    outHistory.WriteRec(-1,"uvTabSave / FITS file "+filename+" disk "+str(outDisk),err)
    outHistory.Close(err)
    History.PCopy2Header (inHistory, outHistory, err)
    OErr.printErrMsg(err, "Error with history")
    # zap table
    outHistory.Zap(err)
    # Copy Tables
    UV.PCopyTables (inUV, outUV, exclude, include, err)
    return outUV  # return new object
    # end uvTabSave

def imlod(filename, inDisk, Aname, Aclass, Adisk, Aseq, err):
    """ Load FITS Image data to AIPS

    Read a ImageTAB FITS Image data file and write an AIPS data set
    filename   = name of FITS file
    inDisk     = FITS directory number
    Aname      = AIPS name of file
    Aclass     = AIPS class of file
    Adisk      = AIPS disk number
    Aseq       = AIPS sequence number of file, 0=> create new
    err        = Python Obit Error/message stack
    returns AIPS Image data object
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Get input
    inImage = Image.newPFImage("FITS Image DATA", filename, inDisk, True, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with FITS data")
    # Get output, create new if seq=0
    if Aseq<1:
        user = OSystem.PGetAIPSuser()
        OErr.printErr(err)   # Print any outstanding messages
        Aseq=AIPSDir.PHiSeq(Adisk,user,Aname,Aclass,"MA",err)
        OErr.PClear(err)     # Clear any message/error
        # If it already exists, increment seq
        if AIPSDir.PTestCNO(Adisk,user,Aname,Aclass,"MA",Aseq,err)>0:
            Aseq = Aseq+1
        OErr.PClear(err)     # Clear any message/error
    print "Creating AIPS MA file",Aname,".",Aclass,".",Aseq,"on disk",Adisk
    outImage = Image.newPAImage("AIPS Image DATA", Aname, Aclass, Adisk, Aseq, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating AIPS data")
    # Copy
    Image.PCopy (inImage, outImage, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying Image data to FITS")
    # Copy History
    inHistory  = History.History("inhistory",  inImage.List, err)
    outHistory = History.History("outhistory", outImage.List, err)
    History.PCopyHeader(inHistory, outHistory, err)
    # Add history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit uvlod",err)
    outHistory.WriteRec(-1,"imlod   / FITS file "+filename+" disk "+str(inDisk),err)
    outHistory.Close(err)
    # Copy Tables
    exclude=["AIPS HI", "AIPS SL", "AIPS PL", "History"]
    include=[]
    Image.PCopyTables (inImage, outImage, exclude, include, err)
    return outImage  # return new Object
    # end imlod

def imtab(inImage, filename, outDisk, err, fract=None, quant=None, \
          exclude=["AIPS HI","AIPS PL","AIPS SL"], include=["AIPS CC"],
          headHi=False):
    """ Write Image data as FITS file

    Write a Image data set as a integer FITAB format file
    History written to header
    inImage    = Image data to copy
    filename   = name of FITS file
    outDisk     = FITS directory number
    err        = Python Obit Error/message stack
    fract      = Fraction of RMS to quantize
    quant      = quantization level in image units, has precedence over fract
                 None or <= 0 => use fract.
    exclude    = List of table types NOT to copy
                 NB: "AIPS HI" isn't really a table and gets copied anyway
    include    = List of table types to copy
    headHi     = if True move history to header, else leave in History table
    returns FITS Image data object
    """
    ################################################################
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Set output
    outImage = Image.newPFImage("FITS Image DATA", filename, outDisk, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating FITS data")
    # Check for valid pixels
    if inImage.Desc.Dict["maxval"]<=inImage.Desc.Dict["minval"]:
        fract=None; quant=None
    # Copy
    if fract or quant:
        Image.PCopyQuantizeFITS (inImage, outImage, err, fract=fract, quant=quant)
    else:
        Image.PCopy (inImage, outImage, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying Image data to FITS")
    # Copy History
    inHistory  = History.History("inhistory",  inImage.List, err)
    outHistory = History.History("outhistory", outImage.List, err)
    History.PCopy(inHistory, outHistory, err)
    # Add this programs history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit imtab",err)
    if fract:
        outHistory.WriteRec(-1,"imtab   / Quantized at "+str(fract)+" RMS",err)
    outHistory.WriteRec(-1,"imtab   / FITS file "+filename+", disk "+str(outDisk),err)
    outHistory.Close(err)
    # History in header?
    if headHi:
        # Copy back to header
        inHistory  = History.History("inhistory",  outImage.List, err)
        History.PCopy2Header (inHistory, outHistory, err)
        # zap table
        outHistory.Zap(err)
    OErr.printErrMsg(err, "Error with history")
    # Copy Tables
    Image.PCopyTables (inImage, outImage, exclude, include, err)
    return outImage   # return new object
    # end imtab

def tabdest (ObitObj, tabType, tabVer):
    """ Delete a table

    Deletes associated tables
    ObitObj   = Python Obit object with tables
    tabType   = Table type,  NB AIPS tables names start with "AIPS "
                e.g. "AIPS CC"
    tabVer    = table version, 0=> highest, <0 => all
    """
    ################################################################
    ObitObj.ZapTable (tabType, tabVer, err)
    # end tabdest

def dhms2day(st):
    """ convert a time string in d/hh:mm:ss.s to days

    Returns time in days
    st        time string as "d/hh:mm:ss.s"
    """
    ################################################################
    stt = st
    if st.__contains__("/"):
        pp=stt.split("/")
        day = int(pp[0])
        stt = pp[1]
    else:
        day = 0
    pp=stt.split(":")
    if len(pp)>0:
        hour = int(pp[0])
    else:
        hour = 0
    if len(pp)>1:
        min = int(pp[1])
    else:
        min = 0
    if len(pp)>2:
        ssec = float(pp[2])
    else:
        ssec = 0.0
    tim = day + hour/24.0 + min/1440.0 + ssec/86400.0
    return tim
    # end dhms2day

def day2dhms(tim):
    """ convert a time in days to a string as d/hh:mm:ss.s

    Returns time as string:  "d/hh:mm:ss.s"
    tim       time in days
    """
    ################################################################
    day=int(tim)
    ttim = 24.0*(tim - day)
    thour = min (int(ttim), 23)
    ttim = 60.0*(ttim - thour)
    tmin = min (int(ttim), 59)
    ssec = 60.0*(ttim - tmin)
    return str(day)+"/"+str(thour).zfill(2)+":"+str(tmin).zfill(2)+\
           ":"+str(ssec)
    # end day2dhms

