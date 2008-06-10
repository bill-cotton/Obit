# Interactive routines to Obit use from ParselTongue
# $Id: PTObit.py,v 1.3 2005/09/21 14:09:21 bcotton Exp $
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
import Obit, Table, FArray, OErr, InfoList, History, AIPSDir, OSystem
import Image, ImageDesc, TableList, ODisplay, UV, OWindow
import os, AIPS
# ParselTongue classes
import AIPSData, FITSData
# from PTObit import *
# ObitStart()

#global Adisk

err=OErr.OErr()
ObitSys=None
Adisk = 1

# Display connection
disp = ODisplay.ODisplay("ObitView", "ObitView", err)

#Initialize Obit system
# Get list of FITS disks
FITSdisks = []
for dsk in ["FITS","FITS01","FITS02","FITS03","FITS04","FITS05","FITS06"]:
    dir = os.getenv(dsk)
    if dir:
        FITSdisks.append(dir)
        nFITS = len(FITSdisks)
        
# Get list of AIPS disks
AIPSdisks = []
for dsk in ["DA01","DA02","DA03","DA04","DA05","DA06","DA07","DA08","DA09","DA10"]:
    dir = os.getenv(dsk)
    if dir:
        AIPSdisks.append(dir)
        nAIPS = len(AIPSdisks)
        
# Init Obit
userno =  AIPS.AIPS.userno
popsno = 1
ObitSys=OSystem.OSystem ("Interactive", popsno, userno, nAIPS, AIPSdisks, \
                         nFITS, FITSdisks, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

def ShowErr(err=err):
    """ Print any errors and clear stack
    
    err  = Python Obit Error/message stack, default of PTObit version
    """
    ################################################################
    OErr.printErrMsg(err, "Error")
    # end ShowErr

def ClearErr(err=err):
    """ Print any errors and clear stack
    
    err  = Python Obit Error/message stack, default of PTObit version
    """
    ################################################################
    OErr.printErrMsg(err, "Error")
    # end ClearErr

def Acat(disk=Adisk, first=1, last=1000):
    """ Catalog listing of AIPS files on disk disk

    The class remembers the last disk accessed
    disk      = AIPS disk number to list
    first     = lowest slot number to list
    last      =highest slot number to list
    """
    ################################################################
    Adisk = disk
    AIPSDir.PListDir(disk, err, first=first, last=last)
    OErr.printErrMsg(err, "Error with AIPS catalog")
    # end Acat

def AMcat(disk=Adisk, first=1, last=1000):
    """ Catalog listing of AIPS Image files on disk disk

    disk      = AIPS disk number to list
    first     = lowest slot number to list
    last      =highest slot number to list
    """
    ################################################################
    Adisk = disk
    AIPSDir.PListDir(disk, err, type=AIPSDir.MAType, first=first, last=last)
    OErr.printErrMsg(err, "Error with AIPS catalog")
    # end AMcat

def AUcat(disk=Adisk, first=1, last=1000):
    """ Catalog listing of AIPS UV data files on disk disk

    disk      = AIPS disk number to list
    first     = lowest slot number to list
    last      = highest slot number to list
    """
    ################################################################
    Adisk = disk
    AIPSDir.PListDir(disk, err, type=AIPSDir.UVType, first=first, last=last)
    OErr.printErrMsg(err, "Error with AIPS catalog")
    # end AUcat

def getname(cno, disk=Adisk):
    """ Return Obit object for AIPS file in cno on disk

    cno       = AIPS catalog slot number 
    disk      = AIPS disk number
    """
    ################################################################
    Adisk = disk
    user =  AIPS.AIPS.userno
    s = AIPSDir.PInfo(disk, user, cno, err)
    OErr.printErrMsg(err, "Error with AIPS catalog")
    # parse returned string
    Aname = s[0:12]
    Aclass = s[13:19]
    Aseq = int(s[20:25])
    Atype = s[26:28]
    if Atype == 'MA':
        out = Image.newPAImage("AIPS image", Aname, Aclass, disk, Aseq, True, err)
        print "AIPS Image",Aname, Aclass, disk, Aseq
    elif Atype == 'UV':
        out = UV.newPAUV("AIPS UV data", Aname, Aclass, disk, Aseq, True, err)
        print "AIPS UV",Aname, Aclass, disk, Aseq
    out.Aname  = Aname
    out.Aclass = Aclass
    out.Aseq   = Aseq 
    out.Atype  = Atype
    out.Disk   = disk
    out.Acno   = cno
    return out
    # end getname

def getFITS(file, disk=Adisk, Ftype='Image'):
    """ Return Obit object for FITS file in file on disk

    file      = FITS file name
    disk      = FITS disk number
    Ftype     = FITS data type: 'Image', 'UV'
    """
    ################################################################
    if Ftype == 'Image':
        out = Image.newPFImage("FITS image", file, disk, True, err)
    elif Ftype == 'UV':
        out = UV.newPFUV("FITS UV data", file, disk, True, err)
    out.Fname  = file
    out.Disk   = disk 
    out.Otype  = Ftype
    return out
    # end getFITS

def tvlod(image, window=None):
    """ display image

    image  = Obit Image, created with getname, getFITS
    window = Optional window for image to edit
    """
    ################################################################
    if Image.PIsA(image):
        # Obit/Image
        ODisplay.PImage(disp, image, err, window=window)
    elif image.__class__==AIPSData.AIPSImage:
        # AIPS Image
        tmp = Image.newPAImage("AIPS Image",image.name, image.klass, image.disk, \
                               image.seq, True, err)
        ODisplay.PImage(disp, tmp, err, window=window)
        del tmp
    elif image.__class__==FITSData.FITSImage:
        # FITS Image
        tmp = Image.newPFImage("FITS Image",image.filename, image.disk, True, err)
        ODisplay.PImage(disp, tmp, err, window=window)
        del tmp
    # end tvlod

def window (image):
    """ Make a window object for an image

    Returns OWindow object
    image  = Obit image object
    """
    ################################################################
    if Image.IsA(image):
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

def imhead (ObitObj):
    """ List header

    ObitObj    = Obit or ParselTongue data object
    """
    ################################################################
    if ObitObj.__class__==AIPSData.AIPSImage:
        # AIPS Image
        tmp = Image.newPAImage("AIPS Image",ObitObj.name, ObitObj.klass, ObitObj.disk, \
                               ObitObj.seq, True, err)
        tmp.Header(err)
        del tmp
    elif ObitObj.__class__==FITSData.FITSImage:
        # FITS Image
        tmp = Image.newPFImage("FITS Image",ObitObj.filename, ObitObj.disk, True, err)
        tmp.Header(err)
        del tmp
    elif ObitObj.__class__==AIPSData.AIPSUVData:
        # AIPS UVData
        tmp = UV.newPAImage("AIPS UVData",ObitObj.name, ObitObj.klass, ObitObj.disk, \
                               ObitObj.seq, True, err)
        tmp.Header(err)
        del tmp
    elif ObitObj.__class__==FITSData.FITSUVData:
        # FITS UVData
        tmp = UV.newPFImage("FITS UVData",ObitObj.filename, ObitObj.disk, True, err)
        tmp.Header(err)
        del tmp
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
    out.DataType = inn.FileType
    out.inDisk   = inn.Disk
    if inn.FileType == 'FITS':
        out.inFile = inn.Fname
    else:   # AIPS
         out.inName  = inn.Aname
         out.inClass = inn.Aclass
         out.inSeq   = inn.Aseq 
    # end setname
   
def set2name (in2, out):
    """ Copy file definition from in2 to out as in2...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    in2  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    out.DataType  = in2.FileType
    out.in2Disk   = in2.Disk
    if in2.FileType == 'FITS':
        out.in2File = in2.Fname
    else:   # AIPS
         out.in2Name  = in2.Aname
         out.in2Class = in2.Aclass
         out.in2Seq   = in2.Aseq 
    # end set2name
   
def setoname (inn, out):
    """ Copy file definition from inn to out as outdisk...

    Supports both FITS and AIPS
    Copies Data type and file name, disk, class etc
    inn  = Obit data object, created with getname, getFITS
    out  = ObitTask object,
    """
    ################################################################
    out.DataType  = inn.FileType
    out.outDisk   = inn.Disk
    if inn.FileType == 'FITS':
        out.outFile = inn.Fname
    else:   # AIPS
         out.outName  = inn.Aname
         out.outClass = inn.Aclass
         out.outSeq   = inn.Aseq 
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
    out.BLC[0] = l[0][2]+1  # make 1-rel
    out.BLC[1] = l[0][3]+1  
    out.TRC[0] = l[0][4]+1  
    out.TRC[1] = l[0][5]+1  
    # end setwindow 
   
def zap (o):
    """ Zap object o

    Removes all external components (files)
    o    = Obit Data object to delete
    """
    ################################################################
    o.Zap(err)
    # end zap
   
