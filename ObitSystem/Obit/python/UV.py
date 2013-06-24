""" 
Python Obit inteferometer (UV) data class

This class contains interoferometric data and allows access.
An ObitUV is the front end to a persistent disk resident structure.
There maybe (usually are) associated tables which either describe
the data or contain calibration and/or editing information.
Both FITS (as Tables) and AIPS cataloged data are supported.
Most access to UV data is through functions as the volume of the data is
inappropriate to be processed directly in python.

UV Members with python interfaces:

=========  ==========================================================
exist      True if object previously existed prior to object creation
List       used to pass instructions to processing
Desc       Astronomical labeling of the data
TableList  List of tables attached
VisBuf     memory pointer into I/O Buffer
=========  ==========================================================

Data selection, calibration and editing parameters on List member:

============== ============== ==================================================
"doCalSelect"  bool (1,1,1)   Select/calibrate/edit data?
"Stokes"       string (4,1,1) Selected output Stokes parameters:
                              "    "=> no translation,"I   ","V   ","Q   ", "U
                              ", "IQU ", "IQUV",  "IV  ", "RR  ", "LL  ", "RL
                              ", "LR  ", "HALF" = RR,LL, "FULL"=RR,LL,RL,LR.
                              [default "    "] In the above 'F' can substitute
                              for "formal" 'I' (both RR+LL).
"BChan"        int (1,1,1)    First spectral channel selected. [def all]
"EChan"        int (1,1,1)    Highest spectral channel selected. [def all]
"chanInc"      int (1,1,1)    channel selected increment. [def 1]
"BIF"          int (1,1,1)    First "IF" selected. [def all]
"EIF"          int (1,1,1)    Highest "IF" selected. [def all]
"IFInc"        int (1,1,1)    "IF" selected increment. [def 1]
"doPol"        int (1,1,1)    >0 -> calibrate polarization.
"doCalib"      int (1,1,1)    >0 -> calibrate, 2=> also calibrate Weights
"gainUse"      int (1,1,1)    SN/CL table version number, 0-> use highest
"flagVer"      int (1,1,1)    Flag table version, 0-> use highest, <0-> none
"BLVer"        int (1,1,1)    BL table version, 0> use highest, <0-> none
"BPVer"        int (1,1,1)    Band pass (BP) table version, 0-> use highest
"Subarray"     int (1,1,1)    Selected subarray, <=0->all [default all]
"dropSubA"     bool (1,1,1)   Drop subarray info?
"FreqID"       int (1,1,1)    Selected Frequency ID, <=0->all [default all]
"timeRange"    float (2,1,1)  Selected timerange in days.
"UVRange"      float (2,1,1)  Selected UV range in kilowavelengths.
"InputAvgTime" float (1,1,1)  Input data averaging time (sec).
                              used for fringe rate decorrelation correction.
"Sources"      string (?,?,1) Source names selected unless any starts with
                              a '-' in which case all are deselected (with '-'
                              sgFtripped).
"souCode"      string (4,1,1) Source Cal code desired, '    ' => any code 
                              selected
                              '*   ' => any non blank code (calibrators only)
                              '-CAL' => blank codes only (no calibrators)
"Qual"         int (1,1,1)    Source qualifier, -1 [default] = any
"Antennas"     int (?,1,1)    a list of selected antenna numbers, if any is 
                              negative then the absolute values are used and
                              the specified antennas are deselected.
"corrtype"     int (1,1,1)    Correlation type, 0=cross corr only, 1=both, 
                              2=auto only.
"passAll"      bool (1,1,1)   If True, pass along all data when 
                              selecting/calibration even if it's all flagged,
                              data deselected by time, source, antenna etc. is
                              not passed.
"doBand"       int (1,1,1)    Band pass application type <0-> none
                              
                              (1) if = 1 then all the bandpass data for each 
                                  antenna will be averaged to form a composite
                                  bandpass spectrum, this will then be used to
                                  correct the data.
                              (2) if = 2 the bandpass spectra nearest in time 
                                  (in a weighted sense) to the uv data point
                                  will be used to correct the data.
                              (3) if = 3 the bandpass data will be interpolated 
                                  in time using the solution weights to form a
                                  composite bandpass spectrum, this
                                  interpolated spectrum will then be used to
                                  correct the data.
                              (4) if = 4 the bandpass spectra nearest in time 
                                  (neglecting weights) to the uv data point
                                  will be used to correct the data.
                              (5) if = 5 the bandpass data will be interpolated 
                                  in time ignoring weights to form a composite
                                  bandpass spectrum, this interpolated spectrum
                                  will then be used to correct the data.
"Smooth"       float (3,1,1)  specifies the type of spectral smoothing
                              
                              Smooth(1) = type of smoothing to apply:
                                0) => no smoothing
                                1) => Hanning
                                2) => Gaussian
                                3) => Boxcar
                                4) => Sinc (i.e. sin(x)/x)
                              Smooth(2) = the "diameter" of the function, i.e.
                                width between first nulls of Hanning triangle
                                and sinc function, FWHM of Gaussian, width of
                                Boxcar. Defaults (if < 0.1) are 4, 2, 2 and 3
                                channels for Smooth(1) = 1 - 4.
                              Smooth(3) = the diameter over which the convolving
                                function has value - in channels.
                                Defaults: 1, 3, 1, 4 times Smooth(2) used when
"Alpha"        float scalar   If != 0.0 then correct data by spectral index 
                              Alpha -0.7 is typical for synchrotron.
"SubScanTime"  float scalar   [Optional] if given, this is the 
                              desired time (days) of a sub scan.  This is used 
                              by the selector to suggest a value close to this
                              which will evenly divide the current scan.  0 =>
                              Use scan average.  This is only useful for
                              ReadSelect operations on indexed ObitUVs.
============== ============== ==================================================
"""
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
 
# Python shadow class to ObitUV class
import OErr, InfoList, Table, AIPSDir, FITSDir, OSystem
import Obit, TableList,  UVDesc, string, UVVis
import OData
#import AIPSData

# class name in C
myClass = "ObitUV"

class UV (OData.OData):
    """ 
Python Obit inteferometer (UV) data class.  UV Members with python interfaces:

=========  ==========================================================
List       used to pass instructions to processing
Desc       Astronomical labeling of the data
TableList  List of tables attached
VisBuf     memory pointer into I/O Buffer
=========  ==========================================================
    """
    def __init__(self, name) :
        self.this = Obit.new_UV(name)
        self.myClass = myClass
    def __del__(self):
        if Obit!=None:
            Obit.delete_UV(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.UVUnref(Obit.UV_me_get(self.this))
            # In with the new
            Obit.UV_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != UV:
            return
        if name == "me" : 
            return Obit.UV_me_get(self.this)
        # Functions to return members
        if name=="List":
            if not self.UVIsA():
                raise TypeError,"input MUST be a Python Obit UV"
            out    = InfoList.InfoList()
            out.me = Obit.InfoListUnref(out.me)
            out.me = Obit.UVGetList(self.cast(myClass))
            return out
        if name=="TableList":
            if not self.UVIsA():
                raise TypeError,"input MUST be a Python Obit UV"
            out    = TableList.TableList("TL")
            out.me = Obit.TableListUnref(out.me)
            out.me = Obit.UVGetTableList(self.cast(myClass))
            return out
        if name=="Desc":
            if not self.UVIsA():
                raise TypeError,"input MUST be a Python Obit UV"
            out    = UVDesc.UVDesc("None")
            out.me = Obit.UVGetDesc(self.cast(myClass))
            return out
        if name=="IODesc":
            if not self.UVIsA():
                raise TypeError,"input MUST be a Python Obit UV"
            out    = UVDesc.UVDesc("None")
            out.me = Obit.UVGetIODesc(self.cast(myClass))
            return out
        if name=="VisBuf":
            if not self.UVIsA():
                raise TypeError,"input MUST be a Python Obit UV"
            return Obit.UVGetVisBuf(self.cast(myClass))
        raise AttributeError,str(name)  # Unknown
    def __repr__(self):
        if self.__class__ != UV:
            return
        return "<C UV instance> " + Obit.UVGetName(self.cast(myClass))
    
    def cast(self, toClass):
        """ 
Casts object pointer to specified class

* self     = object whose cast pointer is desired
* toClass  = Class string to cast to ("ObitUV")
        """
        # Get pointer with type of this class
        out =  self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
            
    def Open (self, access, err):
        """ 
Open a UV data persistent (disk) form

* Returns: 0 on success, else failure
* self   = Python UV object
* access = access READONLY (1), WRITEONLY (2), READWRITE(3), READCAL(4)
* err    = Python Obit Error/message stack
        """
        inUV = self
        # Checks
        if not self.UVIsA():
            raise TypeError,"input MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        ret = Obit.UVOpen(inUV.cast(myClass), access, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Opening UV data")
        return ret
        # end Open
        
    def Read (self, err, firstVis=None):
        """ 
Read a UV  persistent (disk) form.  Reads into buffer attached to UV data, use 
VisBuf for access.

* Returns: 0 on success, else failure
* self     = Python UV object
* err      = Python Obit Error/message stack
* firstVis = If given the first 1-rel visibility in data set
        """
        inUV = self
        # Checks
        if not self.UVIsA():
            raise TypeError,"input MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        # Set first vis?
        if firstVis:
            d = self.IODesc.Dict
            d["firstVis"] = firstVis-1
            self.IODesc.Dict = d
        # Do I/O
        ret = Obit.UVRead (inUV.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Reading UV data")
        return ret
        # end Read
        
    def Write (self, err, firstVis=None):
        """ 
Write a UV  persistent (disk) form. Writes buffer attached to UV data, use 
VisBuf for access.

* returns: 0 on success, else failure
* self     = Python UV object
* err      = Python Obit Error/message stack
* firstVis = If given the first 1-rel visibility in data set
        """
        inUV = self
        # Checks
        if not self.UVIsA():
            raise TypeError,"input MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        # Set first vis?
        if firstVis:
            d = self.IODesc.Dict
            d["firstVis"] = firstVis-1
            self.IODesc.Dict = d
        # Do I/O
        ret = Obit.UVWrite (inUV.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Reading UV data")
        return ret
        # end Write
        
    def ReadVis (self, err, firstVis=None):
        """ 
Read a UV  persistent (disk) form.  Reads into UVVis structure.

* Returns: UVVis structure
* self     = Python UV object
* err      = Python Obit Error/message stack
* firstVis = If given the first 1-rel visibility in data set
        """
        # Checks
        if not self.UVIsA():
            raise TypeError,"input MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        # Set first vis?
        if firstVis:
            d = self.IODesc.Dict
            d["firstVis"] = firstVis-1
            self.IODesc.Dict = d
        # Do I/O
        ret = UVVis.PGet(self, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error Reading UV data")
        return ret
        # end ReadVis
        
    def WriteVis (self, outVis, err, firstVis=None):
        """ 
Write a UVVis to UV  persistent (disk) form. Writes visibility to UV data.

* self      = Python UV object
* outVis    = Vis structure to write
* err       = Python Obit Error/message stack
* firstVis  = If given the first 1-rel visibility in data set
        """
        # Checks
        if not self.UVIsA():
            raise TypeError,"input MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        # Set first vis?
        if firstVis:
            d = self.IODesc.Dict
            d["firstVis"] = firstVis-1
            self.IODesc.Dict = d
        # Do I/O
        UVVis.PSet(outVis, self, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error Writing UV data")
        return
        # end WriteVis
        
    def Close (self, err):
        """ 
Close a UV  persistent (disk) form.

* returns: 0 on success, else failure
* self      = Python UV object
* err       = Python Obit Error/message stack
        """
        inUV = self
        # Checks
        if not self.UVIsA():
            raise TypeError,"input MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        ret = Obit.UVClose (inUV.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Closing UV data")
        return ret
        # end Close
        
    def Copy (self, outUV, err):
        """ 
Make a deep copy of input object. Makes structure the same as self, copies 
data, tables

* self   = Python UV object to copy
* outUV  = Output Python UV object, must be defined
* err    = Python Obit Error/message stack
        """
        # Checks
        if not self.UVIsA():
            raise TypeError,"self MUST be a Python Obit UV"
        if not outUV.UVIsA():
            raise TypeError,"outUV MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        Obit.UVCopy (self.cast(myClass), outUV.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error copying UV data")
    # end Copy

    def Clone (self, outUV, err):
        """ 
Make a copy of a object but do not copy the actual data. This is useful to 
create an UV similar to the input one.

* self   = Python UV object
* outUV  = Output Python UV object, must be defined
* err    = Python Obit Error/message stack
        """
        # Checks
        if not self.UVIsA():
            raise TypeError,"self MUST be a Python Obit UV"
        if not outUV.UVIsA():
            raise TypeError,"outUV MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        Obit.UVClone (self.cast(myClass), outUV.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error copying UV data")
    # end Clone

    def Scratch (self, err):
        """ 
Create a scratch file suitable for accepting the data to be read from self.  A
scratch UV is more or less the same as a normal UV except that it is
automatically deleted on the final unreference.

* self      = Python UV object
* err       = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not self.UVIsA():
            raise TypeError,"self MUST be a Python Obit UV"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        outUV    = UV("None")
        outUV.me = Obit.UVScratch (self.cast(myClass), err.me);
        if err.isErr:
            OErr.printErrMsg(err, "Error creating scratch file")
            
        return outUV
    # end Scratch

    def Header (self, err):
        """ 
Write image header on output

* self   = Python Obit UV object
* err    = Python Obit Error/message stack
        """
        PHeader (self, err)
        # end Header
        
    def Info (self, err):
        """ 
Get underlying data file info

* self   = Python Obit UV object
* err    = Python Obit Error/message stack
        """
        PUVInfo(self, err)
        # end Info
        
    def UpdateDesc (self, err, Desc=None):
        """ Update any disk resident structures about descriptor
        
        self      = Python UV object
        err       = Python Obit Error/message stack
        Desc      = Descriptor, if None then use current descriptor
                    Contents can be accessed throuth the Dict member
        """
        # Checks
        inUV = self
        if not self.UVIsA():
            raise TypeError,"input MUST be a Python Obit UV"
        #
        # if Desc=None make copy of current contents
        if Desc == None:
            d = inUV.Desc.Dict
        else:
            d = Desc.Dict
        # Open for write
        inUV.Open(READWRITE,err)           # Open
        inUV.Desc.Dict = d                 # Update header
        Obit.UVDirty(inUV.cast(myClass))   # force update
        inUV.Close(err)                    # Close to update
        # end UpdateDesc
        
    def UVIsA (self):
        """ Tells if input really a Python Obit UV
        
        return true, false (1,0)
        self   = Python UV object
        """
        ################################################################
        # Allow derived types
        return Obit.UVIsA(self.cast(myClass))
    # end UVIsA

    # End of class member functions (i.e. invoked by x.func())

# Symbolic names for access codes
READONLY  = OData.READONLY  # 1
WRITEONLY = OData.WRITEONLY # 2
READWRITE = OData.READWRITE # 3
READCAL   = 4

def newPFUV(name, filename, disk, exists, err, verbose=True, nvis=1000):
    """ 
Create and initialize an FITS based UV structure.
Create, set initial access information (full image, plane at a time)
and if exists verifies the file. Sets buffer to hold nvis vis.

* Return: the Python UV object. isOK member set to indicate success.
* name     = name desired for object (labeling purposes)
* filename = name of FITS file
* disk     = FITS directory number
* exists   = if true then the file is opened and closed to verify
* err      = Python Obit Error/message stack
* verbose  = If true any give error messages, else suppress
* nvis     = Number of visibilities read/written per call
    """
    ################################################################
    out = UV (name)
    out.isOK = True  # until proven otherwise
    # Does it really previously exist?
    out.exist = FITSDir.PExist(filename, disk, err)
    Obit.UVSetFITS(out.me, nvis, disk, filename, err.me)
    if exists:
        Obit.UVfullInstantiate (out.me, 1, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error verifying UV data")
    # show any errors if wanted
    if  verbose and err.isErr:
        out.isOK = False
        OErr.printErrMsg(err, "newPFUV: Error verifying file")
    elif err.isErr:
        out.isOK = False
        OErr.PClear(err)  # Clear unwanted messages
    
    # It work?
    if not out.isOK:
        return out
    
    out.FileType = 'FITS'
    out.FileName = filename
    out.Fname    = filename
    out.Disk     = disk
    out.Otype    = "UV"
    return out
    # end newPFUV
    
    
def newPAUV(name, Aname, Aclass, disk, seq, exists, err, verbose=True, nvis=1000):
    """ Create and initialize an AIPS based UV structure
    
    Create, set initial access information (full image, plane at a time)
    and if exists verifies the file.
    Sets buffer to hold nvis vis.
    Returns the Python UV object
    isOK member set to indicate success
    name     = name desired for object (labeling purposes)
    Aname    = AIPS name of file
    Aclass   = AIPS class of file
    seq      = AIPS sequence number of file
    disk     = FITS directory number
    exists   = if true then the file is opened and closed to verify
    err      = Python Obit Error/message stack
    verbose  = If true any give error messages, else suppress
    nvis     = Number of visibilities read/written per call
    """
    ################################################################
    out = UV (name)
    user = OSystem.PGetAIPSuser()
    out.isOK = True  # until proven otherwise
    # Does it really previously exist?
    test = AIPSDir.PTestCNO(disk, user, Aname, Aclass, "UV", seq, err)
    out.exist = test>0
    if exists:
        cno = AIPSDir.PFindCNO(disk, user, Aname, Aclass, "UV", seq, err)
        Obit.UVSetAIPS(out.me, nvis, disk, cno, user, err.me)
        Obit.UVfullInstantiate (out.me, 1, err.me)
    else:
        cno = AIPSDir.PAlloc(disk, user, Aname, Aclass, "UV", seq, err)
        Obit.UVSetAIPS(out.me, nvis, disk, cno, user, err.me)

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
    
    out.FileType = 'AIPS'
    out.Disk   = disk
    out.Aname  = Aname
    out.Aclass = Aclass
    out.Aseq   = seq 
    out.Otype  = "UV"
    out.Atype  = "UV"
    out.Acno   = cno
    return out
    # end newPAUV

def newPACNO(disk, cno, exists, err, verbose=True, nvis=1000):
    """ Create and initialize an AIPS based UV structure

    Create, set initial access information 
    and if exists verifies the file.
    Sets buffer to hold nvis vis.
    Returns the Python UV object
    isOK member set to indicate success
    disk     = AIPS directory number
    cno      = AIPS catalog number
    exists   = if true then the file is opened and closed to verify
    err      = Python Obit Error/message stack
    verbose  = If true any give error messages, else suppress
    nvis     = Number of visibilities read/written per call
    """
    ################################################################
    out = UV ("AIPS UVdata")
    user = OSystem.PGetAIPSuser()
    out.isOK = True  # until proven otherwise
    # print "disk, aseq", disk, seq
    # Does it really previously exist?
    test = AIPSDir.PInfo(disk, user, cno, err)
    out.exist = test != None
    
    if exists:
        Obit.UVSetAIPS(out.me, nvis, disk, cno, user, err.me)
        Obit.UVfullInstantiate (out.me, 1, err.me)

    else:
        Obit.UVSetAIPS(out.me, nvis, disk, cno, user, err.me)

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
    out.Otype  = "UV"
    return out   # seems OK
    # end newPACNO

    
def PScratch (inUV, err):
    """ Create a scratch file suitable for accepting the data to be read from inUV
    
    A scratch UV is more or less the same as a normal UV except that it is
    automatically deleted on the final unreference.
    inUV      = Python UV object
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inUV.Scratch(err)
    # end PScratch

def PZap (inUV, err):
    """ Delete underlying files and the basic object.
    
    inUV      = Python UV object
    err       = Python Obit Error/message stack
    """
    ################################################################
    inUV.Zap(err)
    # end PZap

def PRename (inUV, err, newFITSName=None, \
             newAIPSName="            ", \
             newAIPSClass="      ", newAIPSSeq=0):
    """ Rename underlying files

    inUV   = Python UV object
    err       = Python Obit Error/message stack
    For FITS files:
    newFITSName = new name for FITS file

    For AIPS:
    newAIPSName  = New AIPS Name (max 12 char) Blank => don't change.
    newAIPSClass = New AIPS Class (max 6 char) Blank => don't change.
    newAIPSSeq   = New AIPS Sequence number, 0 => unique value
    """
    ################################################################
    inUV.Rename(err,newFITSName=newFITSName, \
                newAIPSName=newAIPSName, newAIPSClass=newAIPSClass,
                newAIPSSeq=newAIPSSeq)
    # end PRename

def PCopy (inUV, outUV, err):
    """ Make a deep copy of input object.

    Makes structure the same as inUV, copies data, tables
    inUV   = Python UV object to copy
    outUV  = Output Python UV object, must be defined
    err    = Python Obit Error/message stack
    """
    ################################################################
    inUV.Copy(outUV, err)
    # end PCopy

def PClone (inUV, outUV, err):
    """ Make a copy of a object but do not copy the actual data

    This is useful to create an UV similar to the input one.
    inUV   = Python UV object
    outUV  = Output Python UV object, must be defined
    err    = Python Obit Error/message stack
    """
    ################################################################
    inUV.Clone(outUV, err)
    # end PClone

def PNewUVTable (inUV, access, tabType, tabVer, err):
    """ Obsolete use PGetTable
    """
    return  PGetTable (inUV, access, tabType, tabVer, err)
# end  PNewUVTable

def PGetTable (inUV, access, tabType, tabVer, err, \
               numOrb=0, numPCal=3, numIF=1, numPol=1, \
               numTerm=0, numChan=1, numTones=1, numBand=1, \
               numTabs=1, npoly=1, numCoef=5,
               maxis1=2, maxis2=1, maxis3=1, maxis4=1, maxis5=1):
    """ Return (create)the specified associated table
    
    Specific table types are recognized and the appropriate constructor
    called, these may have additional parameters.  This allows creating
    new tables of the appropriate type.
    returns Python Obit Table
    inUV      = Python UV object
    access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
    tabType   = Table type, e.g. "AIPS AN"
    tabVer    = table version, if > 0 on input that table returned,
                if 0 on input, the highest version is used.
    err       = Python Obit Error/message stack
    Optional parameters, values only used if table created
    numOrb    = Number of orbital parameters (AN)
    numPCal   = Number of polarization parameters (AN)
    numIF     = Number of IFs (FQ, SN, CL, BP, BL, TY, CQ)
    numPol    = Number of Stokes' (SN, CL, BP, BL, PC, TY, GC, MC, IM)
    numTerm   = Number of terms in model polynomial (CL)
    numChan   = Number of spectral channels (BP)
    numTomes  = Number of Phase cal tones (PC)
    numTabs   = Number of ??? (GC)
    numCoef   = Number of polynomial coefficents (NI)
    numBand   = Number  Bands(?) (IM, GC)
    npoly     = number of polynomial terms (IM)
    maxis1-5  = Dimension of axes of IDI data matrix
    """
    ################################################################
    return inUV.NewTable(access, tabType, tabVer, err,\
                         numOrb=numOrb, numPCal=numPCal, numIF=numIF, \
                         numPol=numPol, numTerm=numTerm, numChan=numChan,\
                         numTones=numTones, numBand=numBand, \
                         numTabs=numTabs, npoly=npoly, numCoef=numCoef)
    # end PGetTable


def POpen (inUV, access, err):
    """ Open an image persistent (disk) form

    inUV   = Python UV object
    access    = access 1=READONLY, 2=WRITEONLY, 3=READWRITE
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inUV.Open(access, err)
    # end POpen

def PDirty (inUV):
    """ Mark UV as needing a header update to disk file

    inUV     = Python UV object
    """
    ################################################################
    inUV.Dirty()
    # end PDirty

def PClose (inUV, err):
    """ Close an image  persistent (disk) form

    inUV   = Python UV object
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inUV.Close( err)
    # end PClose

def PZapTable (inUV, tabType, tabVer, err):
    """ Destroy specified table

    Returns 0 on success
    inUV      = Python UV object
    tabType   = Table type, e.g. "AIPS AN"
    tabVer    = table version, integer
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inUV.ZapTable (tabType, tabVer, err)
    # end PZapTable

def PCopyTables (inUV, outUV, exclude, include, err):
    """ Copy Tabeles from one image to another

    inUV      = Python UV object
    outUV     = Output Python UV object, must be defined
    exclude   = list of table types to exclude (list of strings)
                has priority
    include   = list of table types to include (list of strings)
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inUV.CopyTables (outUV, exclude, include, err)
    # end PCopyTables

def PUpdateTables (inUV, err):
    """ Update any disk resident structures about the current tables

    inUV      = Python UV object
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inUV.UpdateTables (err)
    # end PUpdateTables

def PFullInstantiate (inUV, access, err):
    """ Fully instantiate an UV by opening and closing

    return 0 on success, else failure
    inUV   = Python UV object
    access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inUV.FullInstantiate (access, err)
    # end PfullInstantiate

def PGetList (inUV):
    """ Return the member InfoList

    returns InfoList
    inUV   = Python UV object
    """
    ################################################################
    return inUV.List
    # end PGetList

def PGetTableList (inUV):
    """ Return the member tableList

    returns tableList
    inUV   = Python UV object
    """
    ################################################################
    return inUV.TableList
    # end PGetTableList


def PHeader (inUV, err):
    """ Print data descriptor

    inUV      = Python Obit UV object
    err       = Python Obit Error/message stack
    """
    ################################################################
    # ObitTalk or AIPSUVData data?
    if inUV.myClass == 'AIPSUVData':
        dict = inUV.header()
        UVDesc.PHeaderDict(dict)
        # tablelist
        Tlist = inUV.tables()
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
        # end  AIPSUVData
    elif inUV.myClass == 'FITSUVData':
        dict = inUV.header()
        UVDesc.PHeaderDict(dict)
        # tablelist
        Tlist = inUV.tables()
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
        # end  FITSUVData
    # ObitTalk Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV data"
    #
    # Fully instantiate
    PFullInstantiate (inUV, READONLY, err)
    # File info
    if inUV.FileType=="AIPS":
        print "AIPS UV Data Name: %12s Class: %6s seq: %8d disk: %4d" % \
              (inUV.Aname, inUV.Aclass, inUV.Aseq, inUV.Disk)
    elif inUV.FileType=="FITS":
        print "FITS UV Data Disk: %5d File Name: %s " % \
              (inUV.Disk, inUV.FileName)
    # print in UVDesc
    UVDesc.PHeader(inUV.Desc)
    # Tables
    TL = inUV.TableList
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
    

def PGetDesc (inUV):
    """ Return the member UVDesc

    returns UVDesc as a Python Dictionary
    inUV   = Python UV object
    """
    ################################################################
    return inUV.Desc
    # end PGetDesc

def PGetIODesc (inUV):
    """ Return the member IO UVDesc

    returns IO UVDesc as a Python Dictionary
    inUV   = Python UV object
    """
    ################################################################
     # Checks
    if not PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    #
    out    = UVDesc.UVDesc("None")
    out.me = Obit.UVGetIODesc(inUV.me)
    return out
    # end PGetIODesc

def PUpdateDesc (inUV, err, Desc=None):
    """ Update external representation of descriptor
    inUV   = Python UV object
    err    = Python Obit Error/message stack
    Desc   = UV descriptor, if None then use current descriptor
             Contents can be accessed throuth the Dict member
   """
    ################################################################
    inUV.UpdateDesc (err, Desc=Desc)
    # end PUpdateDesc

def PUVInfo (inUV, err):
    """ Get file info for extant uv data object

    Fills in information on object, useful for scratch files
    inUV   = Python UV object
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # file info
    info = Obit.UVInfo (inUV.cast(myClass), err.me);
    if err.isErr:
        OErr.printErrMsg(err, "Error creating scratch file")
    if info["type"]=="AIPS":
        inUV.FileType = 'AIPS'
        inUV.Disk   = info["disk"]
        inUV.Otype  = "UV"
        inUV.Acno   = info["CNO"]
        # Lookup name etc
        s = AIPSDir.PInfo (inUV.Disk, info["user"], inUV.Acno, err)
        # parse returned string
        inUV.Aname  = s[0:12]
        inUV.Aclass = s[13:19]
        inUV.Aseq   = 0
    
    if info["type"]=="FITS":
        inUV.FileType = 'FITS'
        inUV.Disk   = info["disk"]
        inUV.Otype  = "UV"
        inUV.FileName = info["filename"]
        inUV.Fname    = info["filename"]
    # end PUVInfo

def PGetVisBuf (inUV):
    return inUV.VisBuf

def PGetHighVer (inUV, tabType):
    """ Get highest version number of a specified Table

    returns highest tabType version number, 0 if none.
    inUV   = Python UV object
    tabType   = Table type, e.g. "OTFSoln"
    """
    ################################################################
    return  inUV.GetHighVer (tabType)
    # end PGetHighVer

def PGetSubA (inUV, err):
    """ Get Subarray information

    returns 0 on success, else 1
    inUV   = Python UV object
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.UVGetSubA(inUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error getting subarray info")
    # end PGetSubA 

def PGetFreq (inUV, err):
    """ Get Frequency information

    inUV   = Python UV object
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.UVGetSubA(inUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error getting frequency info")
    return ret
    # end PGetFreq


def PIsScratch (inUV):
    """ Tells if UV is a scratch object

    return true, false (1,0)
    inUV   = Python UV object
    """
    ################################################################
    return inUV.IsScratch ()
    # end PIsScratch

def PIsA (inUV):
    """ Tells if input really a Python Obit UV

    return true, false (1,0)
    inUV   = Python UV object
    """
    ################################################################
    try:
        return inUV.UVIsA()
    except:
        return False
   # end PIsA

def PGetName (inUV):
    """ Tells UV object name (label)

    returns name as character string
    inUV   = Python UV object
    """
    ################################################################
    return inUV.GetName()
    # end PGetName

#----------------------  UVUtil routines    ------------------------

def PUtilUVWExtrema (inUV, err):
    """ Get UV coverage information

    returns array [0]=maximum baseline length (in U,V), [1] = maximum W
    inUV   = Python UV object
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    val = [0.0, 0.0] # create array
    Obit.UVUtil(inUV.cast(myClass), err.me, val)
    if err.isErr:
        OErr.printErrMsg(err, "Error getting  uv coverage info")
    return val
    # end PUTilUVWExtrema

def PUtilCopyZero (inUV, scratch, outUV, err):
    """ Copy a UV data set replacing data by zero, weight 1

    returns UV data object
    inUV   = Python UV object to copy
    scratch= True if this is to be a scratch file (same type as inUV)
    outUV  = Predefined UV data if scratch is False
             ignored if scratch True.
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output for scratch
    if scratch:
        outUV = UV("None")
        outUV.me = Obit.UVUtilCopyZero(inUV.cast(myClass), scratch, outUV.me, err.me)
    else:
        Obit.UVUtilCopyZero(inUV.cast(myClass), scratch, outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying UV data to zeroed UV data")
    
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PUTilCopyZero

def PUtilVisDivide (in1UV, in2UV, outUV, err):
    """ Divides the visibilites in in1UV by those in in2UV

    outUV = in1UV / in2UV
    in1UV   = Numerator  Python UV object, no calibration/selection
    in2UV   = Denominator Python UV object
    outUV   = Output python UV object
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not in1UV.UVIsA():
        raise TypeError,"in1UV MUST be a Python Obit UV"
    if not in2UV.UVIsA():
        raise TypeError,"in2UV MUST be a Python Obit UV"
    if not outUV.UVIsA():
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.UVUtilVisDivide(in1UV.cast(myClass), in2UV.cast(myClass), \
                         outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error dividing UV data")
    # end PUtilVisDivide

def PUtilVisSub (in1UV, in2UV, outUV, err):
    """ Subtracts the visibilites in in2UV from those in in1UV

    outUV = in1UV - in2UV
    in1UV   = First python UV object, no calibration/selection
    in2UV   = Second python UV object, calibration allowed
    outUV   = Output Python UV object, may be same as in1UV
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not in1UV.UVIsA():
        raise TypeError,"in1UV MUST be a Python Obit UV"
    if not in2UV.UVIsA():
        raise TypeError,"in2UV MUST be a Python Obit UV"
    if not outUV.UVIsA():
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.UVUtilVisSub(in1UV.cast(myClass), in2UV.cast(myClass), outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error subtracting UV data")
    # end PUtilVisSub

def PUtilVisCompare (in1UV, in2UV, err):
    """ 
Compares the visibilites in *in1UV* with those in *in2UV*.

* return:  RMS ( real, imaginary differences / amplitude of *inUV2* )
* in1UV = Numerator Python UV object. Possible infoList parameter printRat 
          scalar float if given and >0.0 then tell about entries
          with a real or imaginary difference ratio > printRat
* in2UV = Denominator Python UV object
* err   = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not in1UV.UVIsA():
        raise TypeError,"in1UV MUST be a Python Obit UV"
    if not in2UV.UVIsA():
        raise TypeError,"in2UV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.UVUtilVisCompare(in1UV.cast(myClass), in2UV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error comparing UV data")
    return ret
    # end PUtilVisCompare

def PUtilIndex (inUV, err, maxScan=None, maxGap=None):
    """ Indexes a uv data

    inUV    = Python UV object to index
    err     = Python Obit Error/message stack
    maxScan = max. scan length in min. [def. long]
    maxGap  = max. scan gap in min. [def. long]
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Max scan time specified?
    if maxScan !=None:
        inInfo = PGetList(inUV)    # 
        dim = [1,1,1,1,1]
        InfoList.PAlwaysPutFloat (inInfo, "maxScan", dim, [maxScan])
    # Max scan gap specified?
    if maxGap !=None:
        inInfo = PGetList(inUV)    # 
        dim = [1,1,1,1,1]
        InfoList.PAlwaysPutFloat (inInfo, "maxGap", dim, [maxGap])
    # Index
    Obit.UVUtilIndex(inUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error indexing UV data")
    # end PUtilIndex

def PQuack (inUV, err, begDrop=0.0, endDrop=0.0, Reason="    ", flagVer=1):
    """ Flags beginning and/or end of scans

    inUV    = Python UV object to index
              Selection of source, timerange, IF, channels etc. honored
              Scans are as defined in the iNdeX table
    err     = Python Obit Error/message stack
    begDrop = time in min to drop from the start of each scan
    endDrop = time in min to drop from the end of each scan
    Reason  = Reason string for FG tabls (max. 24 char)
    flagVer = AIPS FG table in which to put entries.
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.UVUtilQuack(inUV.cast(myClass), begDrop, endDrop, Reason,
                     flagVer, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error Quacking UV data")
    # end PQuack

def PUtilHann (inUV, outUV, err, scratch=False):
    """ Hanning smooth a UV data set

    returns smoothed UV data object
    inUV   = Python UV object to smooth
             Any selection editing and calibration applied before average.
    outUV  = Predefined UV data if scratch is False, ignored if
             scratch is True.
    err    = Python Obit Error/message stack
    scratch  = True if this is to be a scratch file (same type as inUV)
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVUtilHann(inUV.cast(myClass), scratch, outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error Hanning UV data")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PUtilHann

def PUtilAvgF (inUV, outUV, err, scratch=False, 
               NumChAvg=0, doAvgAll=False, ChanSel=None):
    """ Average A UV data set in Frequency

    returns Averaged UV data object
    inUV   = Python UV object to copy
             Any selection editing and calibration applied before average.
    outUV  = Predefined UV data if scratch is False, ignored if
             scratch is True.
    err    = Python Obit Error/message stack
    scratch  = True if this is to be a scratch file (same type as inUV)
    NumChAvg = Number of channels to average, [def.0  = all]
    doAvgAll = If TRUE then average all channels and IF.
    ChanSel  =  Groups of channels to consider (relative to channels &
                IFs selected by BChan, EChan, BIF, EIF)
                (start, end, increment, IF) as array of tuples
                where start and end at the beginning and ending
                channel numbers (1-rel) of the group to be included,
                increment is the increment between selected channels
                and IF is the IF number (1-rel)
                default increment is 1, IF=0 means all IF.
                Default is all channels in each IF.
                Example [(3,14,1,0),(25,30,1,0)] averages channels
                3 through 14 and 25 through 30 in each IF.
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Save parameters
    dim = [1,1,1,1,1]
    inInfo = PGetList(inUV)    # 
    InfoList.PAlwaysPutBoolean (inInfo, "doAvgAll",  dim, [doAvgAll])
    InfoList.PAlwaysPutInt     (inInfo, "NumChAvg",  dim, [NumChAvg])
    if ChanSel:
        temp=[]
        for x in ChanSel:
            temp.append(x[0]); temp.append(x[1]); 
            temp.append(x[2]); temp.append(x[3]);
        temp.append(0); temp.append(0); 
        temp.append(0); temp.append(0);
        dim[0] = 4; dim[1] = len(ChanSel)+1
        InfoList.PAlwaysPutInt  (inInfo, "ChanSel",  dim, temp)

    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVUtilAvgF(inUV.cast(myClass), scratch, outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error averaging UV data")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PUtilAvgF

def PUtilAvgT (inUV, outUV, err, scratch=False, timeAvg=1.0):
    """ Average A UV data set in Time

    returns Averaged UV data object
    inUV   = Python UV object to copy
             Any selection editing and calibration applied before average.
    outUV  = Predefined UV data if scratch is False, ignored if
             scratch is True.
    err    = Python Obit Error/message stack
    scratch  = True if this is to be a scratch file (same type as inUV)
    timeAvg  = Averaging time in min
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Save parameter
    dim = [1,1,1,1,1]
    inInfo = PGetList(inUV)    # 
    InfoList.PAlwaysPutFloat (inInfo, "timeAvg",  dim, [timeAvg])
    
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVUtilAvgT(inUV.cast(myClass), scratch, outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error averaging UV data")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PUtilAvgT

def PUtilAvgTF (inUV, outUV, err, scratch=False, \
                FOV=0.333, maxInt=1., maxFact=1.01, \
                NumChAvg=1, doAvgAll=False, ChanSel=[1,0,1,1]):
    """ Average A UV data set in Time and or frequency

    Time averaging is baseline dependent
    inUV   = Python UV object to copy
             Any selection editing and calibration applied before average.
    outUV  = Predefined UV data if scratch is False, ignored if
             scratch is True.
    err    = Python Obit Error/message stack
    scratch = True if this is to be a scratch file (same type as inUV)
    FOV     = Field of view (radius, deg)
    maxInt  = Maximum integration (min)
    maxFact = Maximum time smearing factor
    NumChAvg= Number of channels to average, [0 = all]
    doAvgAll= If TRUE then average all channels and IFs
    ChanSel = Groups of channels to consider (relative to
              channels & IFs selected by BChan, EChan, BIF, EIF)
              (start, end, increment, IF) where start and end at the 
              beginning and ending channel numbers (1-rel) of the group
              to be included, increment is the increment between
              selected channels and IF is the IF number (1-rel)
              default increment is 1, IF=0 means all IF.
              The list of groups is terminated by a start <=0
              """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Save parameters
    dim = [1,1,1,1,1]
    inInfo = PGetList(inUV) 
    InfoList.PAlwaysPutFloat   (inInfo, "FOV",      dim, [FOV])
    InfoList.PAlwaysPutFloat   (inInfo, "maxInt",   dim, [maxInt])
    InfoList.PAlwaysPutFloat   (inInfo, "maxFact",  dim, [maxFact])
    InfoList.PAlwaysPutFloat   (inInfo, "NumChAvg", dim, [NumChAvg])
    InfoList.PAlwaysPutBoolean (inInfo, "doAvgAll", dim, [doAvgAll])
    dim[0] = 4; dim[1] = len(ChanSel)/4
    InfoList.PAlwaysPutInt     (inInfo, "ChanSel",  dim, ChanSel)
    
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVUtilAvgTF(inUV.cast(myClass), scratch, outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error averaging UV data")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PUtilAvgTF

def PUtilCount (inUV, err, timeInt=1440.0):
    """ Count data values by interval in a UV dataset

    Each new source starts a new interval
    returns a dist with entries:
    numTime  = Number of time intervals
    numCorr  = Number of Correlations per vis
    Count    = Number of good correlation/visibilities
    Bad      = Number of flagged correlation/visibilities
    Source   = Source ID per interval (or 0 if no source ID)
    LST      = Average LST (days) per interval
    
    inUV   = Python UV object to copy
             Any selection editing and calibration applied before average.
    err    = Python Obit Error/message stack
    timeInt  = interval  in min (max 500 intervals)
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    #
    list    = InfoList.InfoList()
    # Get values in info List
    list.me = Obit.UVUtilCount (inUV.cast(myClass), timeInt, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error counting UV data")
    # Save values
    out = {}
    dlist = InfoList.PGetDict(list)  # List as dict
    for x in dlist:
        if dlist[x][1][0]==1:  # Scalar or array
            out[x] = dlist[x][2][0]
        else:
            out[x] = dlist[x][2]
    del list, dlist
    
    return out
    # end PUtilCount

def PEditTD (inUV, outUV, err):
    """ Time-domain editing of UV data - produces FG table

    Fill flagging table with clipping by RMS values of the real and imaginary
    parts.  All correlations are clipped on each baseline if the RMS is
    larger than  the maximum.  The clipping is done independently in
    each time interval defined by timeAvg. 
       The clipping level is given by MIN (A, MAX (B,C)) where:
    A = sqrt (maxRMS[0]**2 + (avg_amp * maxRMS[1])**2)
       and avg_amp is the average amplitude on each baseline.
    B = median RMS + 3 * sigma of the RMS distribution.
    C = level corresponding to 3% of the data.
       All data on a given baseline/correlator are flagged if the RMS
    exceeds the limit.  If a fraction of bad baselines on any correlator
    exceeds maxBad, then all data to that correlator is flagged.  In
    addition, if the offending correlator is a parallel hand correlator
    then any corresponding cross hand correlations are also flagged.
    Flagging entries are written into FG table flagTab.
    Control parameters on inUV info member
      "flagTab" int    (1,1,1) FG table version number [ def. 1]
      "timeAvg" float  (1,1,1) Time interval over which to determine 
                data to be flagged (days) [def = 1 min.]
                NB: this should be at least 2 integrations.
      "maxRMS"  float (2,1,1) Maximum RMS allowed, constant plus 
                amplitude coefficient. 
      "maxBad"  float (1,1,1) Fraction of allowed flagged baselines 
                [default 0.25]

    inUV   = Python UV object to clip/flag
    outUV  = UV data onto which the FG table is to be attached.
             May be the same as inUV.
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVEditTD(inUV.cast(myClass), outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error in Time domain edit")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PEditTD

def PEditFD (inUV, outUV, err):
    """ Frequency-domain editing of UV data - produces FG table

    Editing is done independently for each visibility channel.  
    First clipping is done on correlator and Vpol amplitudes.  
    Following this, an average and RMS is determined for each channel 
    in each timeAvg period and a spectral baseline is established
    for the average values, either using a  median window filter (FDwidMW>0) 
    or a linear baseline fit (FDwidMW<=0) to specified channels.  
    Channels with excessive RMSes or residual amplitudes are flagged.
    Flagging is done by entering the offending data in FG table flagTab
    on outUV.
    Control parameters on inUV info member
      "flagTab" int    (1,1,1) FG table version number [ def. 1]
      "timeAvg" float  (1,1,1) Time interval over which to average 
                data to be flagged (days) [def = 1 min.]
      "FDmaxAmp"  float (1,1,1) Maximum average amplitude allowed in the
                spectrum before fitting.  Any channel exceeding this is
                flagged in advance of the  baseline fitting or median
                filtering,. default = infinite 
      "FDmaxV"  float (1,1,1) Maximum average amplitude allowed in V
                polarization; any channel exceeding this is flagged in
                advance of the  baseline fitting or median filtering, 
                Calculates V from difference in amplitudes.   
                default = infinite 
      "FDwidMW" int (1,1,1) If > 0 the width of the median window in channels. 
                An odd number (5) is recommended,  default or 0 => linear baseline
      "FDmaxRMS" float (2,1,1) Flag all channels having RMS 
                values > maxRMS[0] of the channel median sigma.[default = 6.]
                plus maxRMS[1] (default 0.1) of the channel average in quadrature
      "FDmaxRes" float (1,1,1) Max. residual flux in sigma allowed for 
                channels outside the baseline fitting regions.  
                default = 6.
      "FDmaxResBL"  float (1,1,1) Max. residual flux in sigma allowed for 
                channels within the baseline fitting regions. 
                Default = FDmaxRes
      "FDbaseSel"  int (4,*,1) Channel selection to define spectral baseline 
                Used only for linear baseline fitting.
                Select groups of channels/IF(s) to fit as sets 
                of (Start,end,inc,IF), i.e., chanSel = 6,37,1,0, 
                92,123,1,0 for two regions applying to all IFs.  
                Channel and IF numbers 1 -rel
                The first group for which the end channel == 0 terminates the list
                Channel increments defaults to 1
                If the IF==0 then the group applies to all IF.
                Default is channels 2 => nchan-1 all IFs
    inUV   = Python UV object to flag
             Any prior selection and editing is applied.
    outUV  = UV data onto which the FG table is to be attached.
             May be the same as inUV.
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVEditFD(inUV.cast(myClass), outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error in Frequency domain edit")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PEditFD

def PEditStokes (inUV, outUV, err):
    """ Stokes editing of UV data, FG table out

       All data on a given baseline/correlator are flagged if the 
    amplitude of the datatype "FlagStok"  exceeds maxAmp.  
    If a fraction of bad baselines on any antenna/channel/IF exceeds 
    maxBad, then all data to that correlator is flagged.  
    Flagging entries are written into FG table flagTab.
    Results are unpredictable for uncalibrated data.
    Control parameters on info member of inUV:
      "flagStok" string (1,1,1) Stokes value to clip (I, Q, U, V, R, L)
                 default = "V"
      "flagTab" int    (1,1,1) FG table version number [ def. 1]
                NB: this should not also being used to flag the input data!
      "timeAvg" float  (1,1,1) Time interval over which to determine
                data to be flagged (days) [def = 1 min.]
      "maxAmp"  float (1,1,1) Maximum VPol allowed
      "maxBad"  float (1,1,1) Fraction of allowed flagged baselines 
                to an antenna above which all baselines are flagged.
                [default 0.25]
 
    inUV   = Python UV object to clip/flag
    outUV  = UV data onto which the FG table is to be attached.
             May be the same as inUV.
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVEditStokes(inUV.cast(myClass), outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error in clipping")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PEditStokes

def PEditClip (inUV, scratch, outUV, err):
    """ Clip raw visibilities

    control parameters on inUV info member
       "maxAmp" float  (1,1,1) Maximum allowed amplitude
       "oper"   string (4,1,1) operation type:
            "flag" flag data with amplitudes in excess of maxAmp
            "clip" clip amplitudes at maxAmp and preserve phase
               default is "flag"
    returns UV data object
    inUV   = Python UV object to clip/flag
    scratch= True if this is to be a scratch file (same type as inUV)
    outUV  = Predefined UV data if scratch is False, may be inUV
             ignored if scratch True.
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVEditClip(inUV.cast(myClass), scratch, outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error in clipping")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PEditClip

def PEditClipStokes (inUV, scratch, outUV, err):
    """ Flag visibilities by Stokes

    Clip a uv data set.  Data with amplitudes of the selected stokes
    in excess of maxAmp are flagged.  Optionally all correlations associated
    may be flagged.  Stokes conversion as needed for test.
    Control parameters are on the inUV info member:
      "clipStok" string (1,1,1) Stokes value to clip (I, Q, U, V, R, L)
                  default = "I"
       "flagAll"  Obit_bool   (1,1,1) if true, flag all associated correlations
                  default = True
       "maxAmp"   float  (1,1,1) Maximum allowed amplitude
    returns UV data object
    inUV   = Python UV object to clip/flag
    scratch= True if this is to be a scratch file (same type as inUV)
    outUV  = Predefined UV data if scratch is False, may be inUV
             ignored if scratch True.
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if ((not scratch) and (not outUV.UVIsA())):
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Create output for scratch
    if scratch:
        outUV = UV("None")
    outUV.me = Obit.UVEditClipStokes(inUV.cast(myClass), scratch, \
                                     outUV.cast(myClass), err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error in clipping")
    # Get scratch file info
    if scratch:
        PUVInfo (outUV, err)
    return outUV
    # end PEditClipStokes

def PNoise(inUV, outUV, scale, sigma, err):
    """ 
Scale and add Gaussian noise to data in a UV dataset.

out = in*scale + noise(sigma) for each real,imag

* inUV    = input Python Obit UV
* outUV   = output Python Obit UV, must be previously defined
* scale   = multiplicative term
* sigma   = Std. deviation of noise to be added
* err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not outUV.UVIsA():
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"

    if err.isErr: # existing error?
        return
    #
    Obit.UVUtilNoise(inUV.me, outUV.me, scale, sigma, err.me)
    # end PNoise

def PAppend(inUV, outUV, err):
    """ Append the contents of inUV to the end of outUV

    inUV    = input Python Obit UV
    outUV   = output Python Obit UV, must be previously defined
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not outUV.UVIsA():
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"

    if err.isErr: # existing error?
        return
    #
    Obit.UVUtilAppend(inUV.me, outUV.me, err.me)
    # end PAppend

def PFlag (inUV, err,
           flagVer=1, timeRange=[0.0,1.0e20], Ants=[0,0], Source="Any",
           Chans=[1,0], IFs=[1,0], freqID=0, subA=0, Stokes="1111", Reason=" "):
    """ Adds entry to flag table

    Adds flagging table entry.
    inUV      = Python Obit UV on which to write flags
    err       = Python Obit Error/message stack
    flagVer   = flagging table version number
    timeRange = pair of floats giving the beginning and end time in days,
                inclusive, of the data to be flagged
    Source    = Source name, "Any" => all sources.
    Chans     = pair of ints giving first and last spectral channel numbers
                (1-rel) to be flagged; 0s => all
    IFs       = pair of ints giving first and last IF numbers
                (1-rel) to be flagged; 0s => all
    Ants      = first and second antenna  numbers for a baseline, 0=>all
    Stokes    = String giving stokes to be flagged, 
                "FFFF"  where F is '1' to flag corresponding Stokes, '0' not.
                Stokes order 'R', 'L', 'RL' 'LR' or 'X', 'Y', 'XY', 'YX'
    subA      = Subarray
    freqID    = Frequency ID
    Reason    = reason string for flagging (max 24 char)
    """
    ################################################################
    # Checks
    if not PIsA(inUV):
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return
    #
    # Set flagging parameters 
    inInfo = inUV.List
    dim = InfoList.dim
    dim[0] = 1; dim[1] = 1; dim[2] = 1
    InfoList.PAlwaysPutInt(inInfo, "flagVer", dim, [flagVer])
    InfoList.PAlwaysPutInt(inInfo, "subA",    dim, [subA])
    InfoList.PAlwaysPutInt(inInfo, "freqID",  dim, [freqID])
    dim[0] = 2
    InfoList.PAlwaysPutInt(inInfo, "Ants",  dim, Ants)
    InfoList.PAlwaysPutInt(inInfo, "IFs",   dim, IFs)
    InfoList.PAlwaysPutInt(inInfo, "Chans", dim, Chans)
    dim[0] = 2
    InfoList.PAlwaysPutFloat(inInfo, "timeRange", dim, timeRange)
    dim[0] = len(Source)
    InfoList.PAlwaysPutString(inInfo, "Source", dim, [Source])
    dim[0] = len(Stokes)
    InfoList.PAlwaysPutString(inInfo, "Stokes", dim, [Stokes])
    dim[0] = len(Reason)
    InfoList.PAlwaysPutString(inInfo, "Reason", dim, [Reason])
    #
    Obit.UVUtilFlag (inUV.me, err.me)
    # end PFlag

def PTableCLGetDummy (inUV, outUV, ver, err, solInt=10.):
    """ Create dummy CL table table (applying will not modify data)

    Create and return a dummy CL table
    inUV    = input Python Obit UV
    outUV   = output Python Obit UV, must be previously defined
    ver     = version number of new table, 0=> create new
    err     = Python Obit Error/message stack
    solint  = time interval (sec) of table
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not outUV.UVIsA():
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"

    if err.isErr: # existing error?
        return
    
    # Interval
    inInfo = PGetList(inUV)    # 
    dim = [1,1,1,1,1]
    InfoList.PAlwaysPutFloat (inInfo, "solInt", dim, [solInt])
    #
    outTable    = Table.Table("None")
    outTable.me = Obit.TableCLGetDummy(inUV.me, outUV.me, ver, err.me)
    return outTable
    # end PTableCLGetDummy

def PTableCLfromNX(outUV, nant, err, outVer=1, calInt=1.0):
    """
    Create a CL table from an iNdeX table

     * outUV    = Obit UV object
     * nant     = Maximum antenna number
     * err      = Python Obit Error/message stack to init
     * outVer   = output CL table version
     * calInt   = calibration table interval in min
    """
    ################################################################
    # If an old table exists, delete it
    if outUV.GetHighVer("AIPS CL")>=outVer:
        zz = outUV.ZapTable("AIPS CL", outVer, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error zapping old FQ table")
    noif   = outUV.Desc.Dict["inaxes"][outUV.Desc.Dict["jlocif"]]
    npol   = outUV.Desc.Dict["inaxes"][outUV.Desc.Dict["jlocs"]]
    calI = (calInt-(1.0/60))/1440  # Imcrement in days
    cltab = outUV.NewTable(Table.READWRITE, "AIPS CL",outVer,err,numIF=noif,numPol=npol,numTerm=1)
    nxtab = outUV.NewTable(Table.READONLY, "AIPS NX", 1,err)
    if err.isErr:
        OErr.printErrMsg(err, "Error with table")
    cltab.Open(Table.READWRITE, err)
    nxtab.Open(Table.READONLY, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening table")
    # Update header
    d = cltab.Desc.List.Dict
    d["numAnt"] = nant
    cltab.Desc.List.Dict = d
    # Create row
    ia = [];     fa1 = []; fa0 = []
    for iif in range(0,noif):
        ia. append(0)
        fa1.append(1.0)
        fa0.append(0.0)
    row = {'Table name':'AIPS CL', '_status':[1], 'NumFields':34,
           'TIME':[0.0],  'TIME INTERVAL': [0.0], 'ANTENNA NO.':[1], 'SOURCE ID':[1], 'SUBARRAY':[0], 'FREQ ID':[0],
           'ATMOS':[0.0], 'DATMOS':[0.0], 'DOPPOFF':fa0,  'I.FAR.ROT':[0.0],'GEODELAY':[0.0],  
           'MBDELAY1':[0.0], 'DISP 1':[0.0], 'DDISP 1':[0.0], 'CLOCK 1':[0.0],  'DCLOCK 1':[0.0],
           'MBDELAY2':[0.0], 'DISP 2':[0.0], 'DDISP 2':[0.0], 'CLOCK 2':[0.0],  'DCLOCK 2':[0.0],
           'REAL1':fa1, 'IMAG1':fa0,  'DELAY 1':fa0, 'RATE 1':fa0, 'WEIGHT 1':fa1, 'REFANT 1':ia,
           'REAL2':fa1, 'IMAG2':fa0,  'DELAY 2':fa0, 'RATE 2':fa0,'WEIGHT 2':fa1, 'REFANT 2':ia,
            }
    # Loop over NX rows
    nrows = nxtab.Desc.Dict["nrow"]
    irow = -1
    for inxrow in range (1,nrows+1):
        nxrow = nxtab.ReadRow(inxrow,err)
        if err.isErr:
            OErr.printErrMsg(err, "Error reading NX table")
        OErr.printErr(err)
        # Divvy up scans
        ntime = max (1, int(0.5+nxrow['TIME INTERVAL'][0]/calI))
        delta = nxrow['TIME INTERVAL'][0]/ntime
        time = []
        for itime in range (0,ntime):
            time.append(nxrow['TIME'][0]-0.5*nxrow['TIME INTERVAL'][0]+itime*delta)
        # end of scan
        ntime += 1
        time.append(nxrow['TIME'][0]+0.5*nxrow['TIME INTERVAL'][0])
        row['SOURCE ID'][0] = nxrow['SOURCE ID'][0]   # source number
        # Loop over times in scan
        for itime in range (0,ntime):
            row['TIME']          = [time[itime]]
            row['TIME INTERVAL'] = [delta]
            # Loop over antennas
            for iant in range(1,nant+1):
                row['ANTENNA NO.'] = [iant]
                cltab.WriteRow(irow, row,err)
                if err.isErr:
                    OErr.printErrMsg(err, "Error writing CL table")
    cltab.Close(err)
    nxtab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing table")
    # end PTableCLfromNX

def PTableSNGetZeroFR (inUV, outUV, ver, err, solInt=10., timeInt = 1.0):
    """ Create SN tablethat will counter-rotate the data to a zero fringe rate 

    Create and return an SN table
    After this operation, all terrestial sources should be constant.
    Amplitudes reflect the effect of the difference in fringe rate and delay
    for the integration time and observed bandwidth.
    inUV    = input Python Obit UV
    outUV   = output Python Obit UV, must be previously defined
    ver     = version number of new table, 0=> create new
    err     = Python Obit Error/message stack
    solint  = time interval (min) of table
    timeInt = Integration time (sec) of data.
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError,"inUV MUST be a Python Obit UV"
    if not outUV.UVIsA():
        raise TypeError,"outUV MUST be a Python Obit UV"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"

    if err.isErr: # existing error?
        return
    
    # Interval
    inInfo = PGetList(inUV)    # 
    dim = [1,1,1,1,1]
    InfoList.PAlwaysPutFloat (inInfo, "solInt", dim, [solInt])
    # Integration time
    inInfo = PGetList(inUV)    # 
    dim = [1,1,1,1,1]
    InfoList.PAlwaysPutFloat (inInfo, "timeInt", dim, [timeInt])
    #
    outTable    = Table.Table("None")
    outTable.me = Obit.TableSNGetZeroFR(inUV.me, outUV.me, ver, err.me)
    return outTable
    # end PTableSNGetZeroFR

