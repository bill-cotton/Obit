""" Python Obit "On-the-fly" (OTF) single dish data class

This class contains single dish data and allows access.
An ObitOTF is the front end to a persistent disk resident structure.
There maybe (usually are) associated tables which either describe
the data or contain calibration and/or editing information.

OTF Members with python interfaces:
List      - used to pass instructions to processing
Desc      - Astronomical labeling of the image 
TableList - List of tables attached
RecBuf    - memory pointer into I/O Buffer

Additional Functions are available in OTFUtil, OTFSoln2Cal, OTFGetSoln,
OTFGetAtmCor,  CleanOTF

There are a number of utility routines in this module which take
control parameters in the form of python dictionaries
(e.g. AtmCal, Clean, Concat, Image, ResidCal, Soln2Cal, Split)
which each have defined dictionaries with default values and names of the
routine and "Input" appended.
Care should he taken not to change the data types of the entries in these
dictionaries.
These dictionaries can be listed in semi human readable form using the OTF.input
function.

Data selection, calibration and editing parameters on List member
  "doCalSelect" bool (1,1,1) Select/calibrate/edit data?
  "doCalib"     int  (1,1,1) >0 -> calibrate,
  "gainUse"     int  (1,1,1) SN/CL table version number, 0-> use highest
  "flagVer"     int  (1,1,1) Flag table version, 0-> use highest, <0-> none
  "BChan"       int  (1,1,1) First spectral channel selected. [def all]
  "EChan"       int  (1,1,1) Highest spectral channel selected. [def all]
  "Targets"     string (?,?,1) Target names selected. [def all]
  "timeRange"   float (2,1,1) Selected timerange in days. [def all]
  "Scans"       int  (2,1,1) Lowest and highest selected scan numbers. [def all]
  "Feeds"       int  (?,1,1) a list of selected feed numbers, [def all.]
  "keepCal"     bool (1,1,1) If true keep cal-on data, otherwise drop [def True.]
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2008
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

# Obit On the Fly (OTF) calibration and Imaging
import Obit, OErr, Image, ImageDesc, FArray, Table, InfoList, OTFDesc, types
import CleanOTF, OTFUtil, OTFRec, TableList, string
import OData

# Python shadow class to ObitOTF class
# class name in C
myClass = "ObitOTF"

class OTF(OData.OData):
    """ Python Obit "On-the-fly" (OTF) single dish data class

    This class contains single dish data and allows access.
    An ObitOTF is the front end to a persistent disk resident structure.
    There maybe (usually are) associated tables which either describe
    the data or contain calibration and/or editing information.

    OTF Members with python interfaces:
    List      - used to pass instructions to processing
    Desc      - Astronomical labeling of the image 
    TableList - List of tables attached
    RecBuf    - memory pointer into I/O Buffer
   """
    def __init__(self,name) :
        self.this = Obit.new_OTF(name)
        self.thisown = 1
        self.myClass = myClass
    def __del__(self):
        if Obit!=None:
            Obit.delete_OTF(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.OTFUnref(Obit.OTF_me_get(self.this))
            # In with the new
            Obit.OTF_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.OTF_me_get(self.this)
        # Functions to return members
        if name=="List":
            if not self.OTFIsA():
                raise TypeError,"input MUST be a Python Obit OTF"
            out    = InfoList.InfoList()
            out.me = Obit.InfoListUnref(out.me)
            out.me = Obit.OTFGetList(self.cast(myClass))
            return out
        if name=="TableList":
            if not self.OTFIsA():
                raise TypeError,"input MUST be a Python Obit OTF"
            out    = TableList.TableList("TL")
            out.me = Obit.TableListUnref(out.me)
            out.me = Obit.OTFGetTableList(self.cast(myClass))
            return out
        if name=="Desc":
            if not self.OTFIsA():
                raise TypeError,"input MUST be a Python Obit OTF"
            out    = OTFDesc.OTFDesc("None")
            out.me = Obit.OTFGetDesc(self.cast(myClass))
            return out
        if name=="RecBuf":
            if not self.OTFIsA():
                raise TypeError,"input MUST be a Python Obit OTF"
            return Obit.OTFGetRecBuf(self.cast(myClass))
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != OTF:
            return
        if self==None:
            return "None"
        return "<C OTF instance> " + Obit.OTFGetName(self.me)

    def cast(self, toClass):
        """ Casts object pointer to specified class
        
        self     = object whose cast pointer is desired
        toClass  = Class string to cast to ("ObitOTF")
        """
        # Get pointer with type of this class
        out =  self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
            
    def NewTable (self, access, tabType, tabVer, err,
                  numDet=1, numPoly=0, numParm=0):
        """ Return the specified associated table
        
        self      = Python OTF object
        access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
        tabType   = Table type, e.g. "OTFSoln"
        tabVer    = table version, if > 0 on input that table returned,
                    if 0 on input, the highest version is used.
        err       = Python Obit Error/message stack
        Optional parameters, values only used if table created
        numDet    = Number of Detectors (OTFCal, OTFSoln, OTFScanData)
        numPoly   = Number of polynomial terms (OTFCal, OTFSoln)
        numParm   = Number of model parameters (OTFModel)
        """
        inOData = self
        # Checks
        if not self.ODataIsA():
            raise TypeError,"input MUST be a Python Obit OData"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        outTable    = Table.Table("None")
        id  = inOData.cast(OData.myClass)  # Cast to OData
        outTab = None
        if tabType=="OTFArrayGeom":
            outTab = Obit.TableOTFArrayGeom(id, [tabVer], access, tabType, err.me)
        elif tabType=="OTFCal":
            outTab = Obit.TableOTFCal(id, [tabVer], access, tabType, numDet, numPoly, err.me)
        elif tabType=="OTFFlag":
            outTab = Obit.TableOTFFlag(id, [tabVer], access, tabType, err.me)
        elif tabType=="OTFIndex":
            outTab = Obit.TableOTFIndex(id, [tabVer], access, tabType, err.me)
        elif tabType=="OTFModel":
            outTab = Obit.TableOTFModel(id, [tabVer], access, tabType, numParm, err.me)
        elif tabType=="OTFScanData":
            outTab = Obit.TableOTFScanData(id, [tabVer], access, tabType, numDet, err.me)
        elif tabType=="OTFSoln":
            outTab = Obit.TableOTFSoln(id, [tabVer], access, tabType, numDet, numPoly, err.me)
        elif tabType=="SkyModel":
            outTab = Obit.TableSkyModel(id, [tabVer], access, tabType, numDet, numPoly, err.me)
        elif tabType=="OTFTarget":
            outTab = Obit.TableOTFTarget(id, [tabVer], access, tabType, err.me)
        else:  # Generic
            ret = Obit.newODataTable (inOData.cast(myClass), access, tabType, [tabVer], err.me)
            # Table and version returned in a list
            outTab = ret[0]
        # Did it work? error may not be set
        if (outTab.__class__!=str) or (string.find(outTab,"_ObitTable_p") <0):
            OErr.printErrMsg(err, "Error getting OData data Table")
            raise RuntimeError,"Failed to extract "+tabType+" table from "+inOData.GetName()
        # Error?
        if err.isErr:
            OErr.printErrMsg(err, "Error getting OData data Table")
        # Test validity of outTab
        if not Obit.TableIsA(outTab):
            OErr.printErrMsg(err, "Error getting OData data Table")
            raise RuntimeError,"Failed to extract "+tabType+" table from "+inOData.GetName()
        # Open and close to fully instantiate - should exist
        outTable.me = outTab
        Table.PFullInstantiate (outTable, access, err)
        # Make sure that it worked - the output should be a table
        if not Table.PIsA(outTable):
            raise RuntimeError,"Failed to extract "+tabType+" table from "+inOData.GetName()
        return outTable
        # end NewTable

    def Open (self, access, err):
        """ Open a OTF data persistent (disk) form

        Returns 0 on success, else failure
        self   = Python OTF object
        access = access READONLY (1), WRITEONLY (2), READWRITE(3)
        err    = Python Obit Error/message stack
        """
        inOTF = self
        # Checks
        if not self.OTFIsA():
            raise TypeError,"input MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        ret = Obit.OTFOpen(inOTF.cast(myClass), access, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Opening OTF data")
        return ret
        # end Open
        
    def Read (self, err):
        """ Read a OTF  persistent (disk) form
        
        Reads into buffer attached to OTF data, use VisBuf for access
        Returns 0 on success, else failure
        self   = Python OTF object
        err    = Python Obit Error/message stack
        """
        inOTF = self
        # Checks
        if not self.OTFIsA():
            raise TypeError,"input MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        ret = Obit.OTFRead (inOTF.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Reading OTF data")
        return ret
        # end Read
        
    def Write (self, err):
        """ Write a OTF  persistent (disk) form

        Writes buffer attached to OTF data, use VisBuf for access
        returns 0 on success, else failure
        self      = Python OTF object
        err       = Python Obit Error/message stack
        """
        inOTF = self
        # Checks
        if not self.OTFIsA():
            raise TypeError,"input MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        ret = Obit.OTFWrite (inOTF.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Writing OTF data")
        return ret
        # end Write
        
    def ReadRec (self, err):
        """ Read a OTF  persistent (disk) form
        
        Returns OTFRec structure from next record
        self   = Python OTF object
        err    = Python Obit Error/message stack
        """
        # Checks
        if not self.OTFIsA():
            raise TypeError,"input MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        ret = OTFRec.PGet(self, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error Reading OTF data")
        return ret
        # end ReadRec
        
    def WriteRec (self, outRec, err):
        """ Write a OTF  persistent (disk) form

        Writes buffer attached to OTF data, use VisBuf for access
        returns 0 on success, else failure
        self      = Python OTF object
        outRec    = OTFRec structure to write
        err       = Python Obit Error/message stack
        """
        # Checks
        if not self.OTFIsA():
            raise TypeError,"input MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        OTFRec.PSet(outRec, self, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error Writing OTF data")
        # end Write
        
    def Close (self, err):
        """ Close a OTF  persistent (disk) form
        
        returns 0 on success, else failure
        self      = Python OTF object
        err       = Python Obit Error/message stack
        """
        inOTF = self
        # Checks
        if not self.OTFIsA():
            raise TypeError,"input MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        ret = Obit.OTFClose (inOTF.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Closing OTF data")
        return ret
        # end Close
        
    def Copy (self, outOTF, err):
        """ Make a deep copy of input object.
        
        Makes structure the same as self, copies data, tables
        self   = Python OTF object to copy
        outOTF  = Output Python OTF object, must be defined
        err    = Python Obit Error/message stack
        """
        # Checks
        if not self.OTFIsA():
            raise TypeError,"self MUST be a Python Obit OTF"
        if not outOTF.OTFIsA():
            raise TypeError,"outOTF MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        Obit.OTFCopy (self.cast(myClass), outOTF.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error copying OTF data")
    # end Copy

    def Clone (self, outOTF, err):
        """ Make a copy of a object but do not copy the actual data
        
        This is useful to create an OTF similar to the input one.
        self   = Python OTF object
        outOTF  = Output Python OTF object, must be defined
        err    = Python Obit Error/message stack
        """
        # Checks
        if not self.OTFIsA():
            raise TypeError,"self MUST be a Python Obit OTF"
        if not outOTF.OTFIsA():
            raise TypeError,"outOTF MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        Obit.OTFClone (self.cast(myClass), outOTF.cast(myClass), err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error copying OTF data")
    # end Clone

    def Scratch (self, err):
        """ Create a scratch file suitable for accepting the data to be read from self
        
        A scratch OTF is more or less the same as a normal OTF except that it is
        automatically deleted on the final unreference.
        self      = Python OTF object
        err       = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not self.OTFIsA():
            raise TypeError,"self MUST be a Python Obit OTF"
        if not OErr.OErrIsA(err):
            raise TypeError,"err MUST be an OErr"
        #
        outOTF    = OTF("None")
        outOTF.me = Obit.OTFScratch (self.cast(myClass), err.me);
        if err.isErr:
            OErr.printErrMsg(err, "Error creating scratch file")
        outOTF.Info(err)          # Update info
        return outOTF
    # end Scratch

    def Header (self, err):
        """ Write image header on output
        
        self   = Python Obit OTF object
        err    = Python Obit Error/message stack
        """
        PHeader (self, err)
        # end Header
        
    def Info (self, err):
        """ Get underlying data file info
        
        self   = Python Obit OTF object
        err    = Python Obit Error/message stack
        """
        POTFInfo(self, err)
        # end Info
        
    def UpdateDesc (self, err, Desc=None):
        """ Update any disk resident structures about descriptor
        
        self      = Python OTF object
        err       = Python Obit Error/message stack
        Desc      = Descriptor, if None then use current descriptor
                    Contents can be accessed throuth the Dict member
        """
        # Checks
        inOTF = self
        if not self.OTFIsA():
            raise TypeError,"input MUST be a Python Obit OTF"
        #
        # if Desc=None make copy of current contents
        if Desc == None:
            d = inOTF.Desc.Dict
        else:
            d = Desc.Dict
        # Open for write
        inOTF.Open(READWRITE,err)           # Open
        inOTF.Desc.Dict = d                 # Update header
        Obit.OTFDirty(inOTF.cast(myClass))   # force update
        inOTF.Close(err)                    # Close to update
        # end UpdateDesc
        
    def OTFIsA (self):
        """ Tells if input really a Python Obit OTF
        
        return true, false (1,0)
        self   = Python OTF object
        """
        ################################################################
        # Allow derived types
        return Obit.OTFIsA(self.cast(myClass))
    # end PIsA
    # End of class member functions (i.e. invoked by x.func())
 
err=OErr.OErr()

# Commonly used, dangerous variables
dim=[1,1,1,1,1]
blc=[1,1,1,1,1,1,1]
trc=[0,0,0,0,0,0,0]

# Symbolic names for access codes
READONLY  = OData.READONLY  # 1
WRITEONLY = OData.WRITEONLY # 2
READWRITE = OData.READWRITE # 3
def ObitName(ObitObject):
    """Return name of an Obit object or input if not an Obit Object
    """
    ################################################################
    out = ObitObject    # in case
    print "\n type ",ObitObject.me
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

def newPOTF(name, filename, disk, exists, err, nrec=1000):
    """ Create and initialize an OTF structure

    Create, set initial access information (nrec records)
    and if exists verifies the file.
    Returns the Python OTF object
    name     = name desired for object (labeling purposes)
    filename = name of FITS file
    disk     = FITS directory number
    exists   = if true then the file is opened and closed to verify
    err      = Python Obit Error/message stack
    nrec     = Number of records read/written per call
    """
    ################################################################
    out = OTF(name)
    if err.isErr: # existing error?
        return out
    Obit.OTFSetFITS(out.me, nrec, disk, filename, err.me)
    if exists:
        Obit.OTFfullInstantiate (out.me, READWRITE, err.me)
    # show any errors 
    OErr.printErrMsg(err, "newPOTF: Error verifying file")
    out.FileType = 'FITS'
    out.FileName = filename
    out.Fname    = filename
    out.Disk     = disk
    out.Otype    = "OTF"
    return out      # seems OK
    # end newPOTF
 
def ClearCal(inOTF, err):
    """ Delete calibration tables on an OTF

    Removes all OTFSoln and OTFCal tables
    inOTF    = Extant Python OTF
    err      = Python Obit Error/message stack
    """
    ################################################################
    if err.isErr: # existing error?
        return 
    #Obit.OTFfullInstantiate (inOTF.cast(myClass), READWRITE, err.me)
    inOTF.Open(READWRITE, err) 
    inOTF.Close(err) 
    OErr.printErrMsg(err, "ClearCal: Error verifying file")
    ver = Obit.OTFGetHighVer(inOTF.cast(myClass), "OTFCal")
    while (ver>0):
        Obit.OTFZapTable (inOTF.cast(myClass), 'OTFCal', ver, err.me)
        ver = ver-1
        OErr.printErrMsg(err, "ClearCal: Error removing OTFCal")
    ver = Obit.OTFGetHighVer(inOTF.cast(myClass), "OTFSoln")
    while (ver>0):
        Obit.OTFZapTable (inOTF.cast(myClass), 'OTFSoln', ver, err.me)
        ver = ver-1
        OErr.printErrMsg(err, "ClearCal: Error removing OTFSoln")
    
    # end ClearCal

# Define AtmCal input dictionary
AtmCalInput={'structure':['AtmCal',[('InData','Input OTF'),
                                    ('solInt','Solution interval(sec)'),
                                    ('Tau0','Zenith opacity'),
                                    ('minEl','Min elev. (deg)'),
                                    ('aTemp','Atm. temperature (K)'),
                                    ('tRx','Recvr. Temp. (K)'),
                                    ('calJy','Cal. signal in Jy'),
                                    ('raOff','RA pointing offset (deg)'),
                                    ('decOff','Dec pointing offset (deg)')]],
             'InData':None, 'solInt':10000.0, 'Tau0':0.0, 'minEl':0.0,
             'aTemp':[0.0,0.0], 'tRx':[0.0,0.0], 'calJy':[1.0,1.0],
             'raOff':0.0, 'decOff':0.0}
def AtmCal (err, input=AtmCalInput):
    """ Basic atmospheric calibration.

    Applies Atmospheric calibration and optionally gross pointing offsets
    Returns the version number of the Soln Table on success.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData = input Python OTF to calibrate
    solInt = solution interval (sec)
    Tau0   = zenith opacity (nepers)
    minEl  = minimum elevation (deg)
    tTemp  = effective atmospheric temperature (per detector)
    tRx    = Receiver temperature per detector (K)
    calJy  = Noise cal value in Jy per detector
    raOff  = RA pointing offset (deg)
    decOff = Dec pointing offset (deg)
    """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    solInt  = input["solInt"]/86400.0   # convert to days
    Tau0    = input["Tau0"]
    minEl   = input["minEl"]
    aTemp   = input["aTemp"]
    tRx     = input["tRx"]
    calJy   = input["calJy"]
    RAoff   = input["raOff"]
    Decoff  = input["decOff"]
    # Checks
    if not PIsA(inData):
        raise TypeError,'AtmCal: Bad input OTF'
    #
    # Set calibration parameters
    inInfo = inData.List;
    inInfo.set("solInt",  solInt)
    inInfo.set("Tau0",    Tau0)
    inInfo.set("minEl",   minEl)
    inInfo.set("RAoff",   RAoff)
    inInfo.set("Decoff",  Decoff)
    inInfo.set("aTemp",   aTemp)
    inInfo.set("tRx",     tRx)
    inInfo.set("calJy",   calJy)
    #
    # Determine calibration (Obit object)
    solnTable = Obit.OTFGetAtmCor (inData.me, inData.me, err.me);
    #
    # show any errors 
    OErr.printErrMsg(err, "AtmCal: Error determining calibration")
    #
    # Get table version number
    tabVer = Obit.TableGetVer(solnTable)
    #
    # Cleanup Obit objects
    solnTable = Obit.TableUnref (solnTable)
    #
    return tabVer
    # end AtmCal

# Define PolyBLCal input dictionary
PolyBLCalInput={'structure':['PolyBLCal',[('InData','Input OTF'),
                                    ('solInt','Solution interval(sec)'),
                                    ('order','polynomial order'),
                                    ('gainUse','cal. table version, -1=none'),
                                    ('flagVer','flag table version, -1=none'),
                                    ('minEl','Min elev. (deg)')]],
             'InData':None, 'solInt':10.0, 'order':1, 'minEl':0.0,
             'gainUse':-1, 'flagVer':-1}
def PolyBLCal (err, input=PolyBLCalInput):
    """ Polynomial baseline fit to residual data

    Each solution interval in a scan is median averaged
    (average of 9 points around the median) and then a polynomial fitted.
    Returns the version number of the Soln Table on success.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData = input Python OTF to calibrate
    solInt = solution interval (sec)
    order  = polynomial order
    minEl  = minimum elevation (deg)
    gainUse = version number of prior table (Soln or Cal) to apply, -1 is none
    flagVer = version number of flagging table to apply, -1 is none
    """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    solInt  = input["solInt"]/86400.0   # convert to days
    order   = input["order"]
    minEl   = input["minEl"]
    gainUse = input["gainUse"]
    flagVer = input["flagVer"]
    # Checks
    if not PIsA(inData):
        raise TypeError,'PolyBLCal: Bad input OTF'
    if err.isErr: # existing error?
        return None
    #
    # Set calibration parameters
    dim[0] = 1; dim[1] = 1
    inInfo = inData.List;
    InfoList.PAlwaysPutFloat(inInfo, "solInt",  dim, [solInt]);
    InfoList.PAlwaysPutInt(inInfo,   "order",   dim, [order]);
    InfoList.PAlwaysPutFloat(inInfo, "minEl",   dim, [minEl]);
    InfoList.PAlwaysPutInt  (inInfo, "flagVer", dim, [flagVer])
    InfoList.PAlwaysPutInt  (inInfo, "gainUse", dim, [gainUse])
    if gainUse>=0:
        itemp = 1
    else:
        itemp = -1
    InfoList.PAlwaysPutInt (inInfo, "doCalib", dim, [itemp])
    doCalSelect = (gainUse >= 0) or (flagVer>0)
    doClaSelect = True
    InfoList.PAlwaysPutBoolean (inInfo, "doCalSelect", dim, [doCalSelect])
    #
    # Determine calibration (Obit object)
    solnTable = Obit.OTFGetSolnPolyBL (inData.me, inData.me, err.me);
    #
    # show any errors 
    OErr.printErrMsg(err, "PolyBLCal: Error determining calibration")
    #
    # Get table version number
    tabVer = Obit.TableGetVer(solnTable)
    #
    # Cleanup Obit objects
    solnTable = Obit.TableUnref (solnTable)
    #
    return tabVer
    # end PolyBL

# Define MBBaseCal input dictionary
MBBaseCalInput={'structure':['MBBaseCal',[('InData','Input OTF'),
                                    ('solInt','Solution interval(sec)'),
                                    ('order','polynomial order'),
                                    ('gainUse','cal. table version, -1=none'),
                                    ('flagVer','flag table version, -1=none'),
                                    ('clipSig','data outside of +/- clipsig ignored [def large]'),
                                    ('plotDet','Detector number (1-rel) to plot per scan [def =-1 = none]'),
                                    ('minEl','Min elev. (deg)')]],
             'InData':None, 'solInt':10.0, 'order':1, 'clipSig':1.0e20, 'plotDet':-1,'minEl':0.0,
             'gainUse':-1, 'flagVer':-1}

def MBBaseCal (err, input=MBBaseCalInput):
    """ Continuum baseline fitting for multibeam instrument.

    Fit one term, time variable common, atmospheric polynomial and a single offset
    per detector.
    Since the different detectors each have an individual multiplicative term, the 
    Atmospheric + offset are places in the the detector's additive term and the
    polynomial is set to zero.
    Scans in excess of 5000 samples will be broken into several.
    Returns the version number of the Soln Table on success.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData = input Python OTF to calibrate
    solInt = solution interval (sec), entries 4 times per SolInt
    order  = polynomial order
    clipSig = Data outside of +/- clipsig ignored [def large]
    plotDet = Detector number (1-rel) to plot per scan [def =-1 = none]
    minEl   = minimum elevation (deg)
    gainUse = version number of prior table (Soln or Cal) to apply, -1 is none
    flagVer = version number of flagging table to apply, -1 is none
    """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    Solint  = input["Solint"]/86400.0   # convert to days
    order   = input["order"]
    minEl   = input["minEl"]
    clipSig = input["clipSig"]
    plotDet = input["plotDet"]
    gainUse = input["gainUse"]
    flagVer = input["flagVer"]
    # Checks
    if not PIsA(inData):
        raise TypeError,'MBBaseCal: Bad input OTF'
    if err.isErr: # existing error?
        return None
    #
    # Set calibration parameters
    dim[0] = 1; dim[1] = 1
    inInfo = inData.List;
    InfoList.PAlwaysPutFloat(inInfo, "solInt",   dim, [solInt]);
    InfoList.PAlwaysPutFloat(inInfo, "clipSig",  dim, [clipSig]);
    InfoList.PAlwaysPutFloat(inInfo, "minEl",    dim, [minEl]);
    InfoList.PAlwaysPutInt(inInfo, "plotDet", dim, [plotDet]);
    InfoList.PAlwaysPutInt(inInfo, "Order",   dim, [order]);
    InfoList.PAlwaysPutInt(inInfo, "flagVer", dim, [flagVer])
    InfoList.PAlwaysPutInt(inInfo, "gainUse", dim, [gainUse])
    if gainUse>=0:
        itemp = 1
    else:
        itemp = -1
    InfoList.PAlwaysPutInt (inInfo, "doCalib", dim, [itemp])
    doCalSelect = (gainUse >= 0) or (flagVer>0)
    doClaSelect = True
    InfoList.PAlwaysPutBoolean (inInfo, "doCalSelect", dim, [doCalSelect])
    #
    # Determine calibration (Obit object)
    solnTable = Obit.OTFGetSolnMBBase (inData.me, inData.me, err.me);
    #
    # show any errors 
    OErr.printErrMsg(err, "MBBaseCal: Error determining calibration")
    #
    # Get table version number
    tabVer = Obit.TableGetVer(solnTable)
    #
    # Cleanup Obit objects
    solnTable = Obit.TableUnref (solnTable)
    #
    return tabVer
    # end MBBase

# Define ResidCal input dictionary
ResidCalInput={'structure':['ResidCal',[('InData','Input OTF'),
                                        ('Model','Model FArray'),
                                        ('ModelDesc','Model Descriptor'),
                                        ('solType','Soln type, Filter,Offset,Gain'),
                                        ('solInt','Solution interval (sec)'),
                                        ('minEl','Minimum elev (deg)'),
                                        ('minRMS','min. RMS residual (gain)'),
                                        ('minFlux','Minimum flux density in Model to use'),
                                        ('Clip','Clipping level for residuals'),
                                        ('calJy','Cal. signal in Jy'),
                                        ('gainUse','cal. table version, -1=none'),
                                        ('flagVer','flag table version, -1=none')]],
               'InData':None, 'Model':None, 'ModelDesc':None, 'solType':"Filter",
               'solInt':10000.0, 'minEl':0.0, 'minRMS':0.0, 'minFlux':-10000.0,
               'Clip':1.0e20,'calJy':[1.0,1.0], 'gainUse':-1, 'flagVer':-1};
def ResidCal (err, input=ResidCalInput):
    """ Determine residual calibration for an OTF.

    Determines a solution table for an OTF by one of a number of techniques using
    residuals from a model image.
    Returns the version number of the Soln Table on success.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData  = Python input OTF to calibrate
    Model   = Python input model FArray, "None" means do not subtract model image
    ModelDesc= Python input model ImageDesc
    minFlux = Minimum brightness in model
    solInt = solution interval (sec)
    solType = solution type:
       "Gain" solve for multiplicative term from "cals" in data.
             (solInt, minRMS, minEl, calJy)
       "Offset" Solve for additive terms from residuals to the model.
             (solInt, minEl)
       "GainOffset" Solve both gain and offset
             (solInt, minRMS, minEl, calJy)
       "Filter"  Additive terms from filters residuals to the model.
             (solInt, minEl)
       "MultiBeam" Multibeam solution
             (solInt, minEl)
    minEl  = minimum elevation (deg)
    minRMS = Minimum RMS residual to solution
    calJy  = Noise cal value in Jy per detector 
    gainUse = version number of prior table (Soln or Cal) to apply, -1 is none
    flagVer = version number of flagging table to apply, -1 is none
   """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    model   = input["Model"]
    modDesc = input["ModelDesc"]
    calType = input["solType"]
    solInt  = input["solInt"]/86400.0   # convert to days
    minEl   = input["minEl"]
    minRMS  = input["minRMS"]
    minFlux = input["minFlux"]
    Clip    = input["Clip"]
    calJy   = input["calJy"]
    gainUse = input["gainUse"]
    flagVer = input["flagVer"]
    # Checks
    if not PIsA(inData):
        raise TypeError,'ResidCal: Bad input OTF'
    if err.isErr: # existing error?
        return None
    # 
    # Apply prior calibration as requested 
    inInfo = inData.List
    inInfo.set("flagVer", flagVer)
    inInfo.set("gainUse", gainUse)
    if gainUse>=0:
        itemp = 1
    else:
        itemp = -1
    inInfo.set("doCalib", itemp)
    doCalSelect = (gainUse >= 0) or (flagVer>0)
    doClaSelect = True
    inInfo.set("doCalSelect", doCalSelect)
    #
    if model != None:
        zapIt = True     #  A scratch file - delete
        # clip image below minFlux
        print "Clip model below ", minFlux
        FArray.PClip (model, minFlux, 1.0e20, 0.0)
        # Scratch file for residual data 
        scrData = PScratch (inData, err)
        OErr.printErrMsg(err, "ResidCal: Error creating scratch file")
        #  
        # Subtract image from inData to scrData
        OTFUtil.PSubImage(inData, scrData, model, modDesc, err)
        # error check
        OErr.printErrMsg(err, "ResidCal: Error subtracting model")
    else:
        scrData = inData       # use input data
        zapIt = False          # Not a scratch file - don't delete
    #
    # Set calibration parameters 
    scrInfo = scrData.List
    scrInfo.set("solInt", solInt)
    scrInfo.set("minRMS", minRMS)
    scrInfo.set("minEl", minEl)
    scrInfo.set("Clip", Clip)
    scrInfo.set("calJy", calJy)
    scrInfo.set("calType", calType)
    #
    # Determine calibration by type 
    if calType == "Filter": # Time filtering 
        solnTable = Obit.OTFGetSolnFilter (scrData.me, inData.me, err.me)
    elif calType == "MultiBeam": # Multibeam
        solnTable = Obit.OTFGetSolnCal (scrData.me, inData.me, err.me)
    else: # offset or gain 
        solnTable = Obit.OTFGetSolnGain (scrData.me, inData.me, err.me)
    #
    # show any errors 
    OErr.printErrMsg(err, "ResidCal: Error determining calibration")
    #
    # Get table version number
    tabVer = Obit.TableGetVer(solnTable)
    #
    # Cleanup Obit Objects
    if zapIt:
        scrData   = PZap(scrData, err)   # Zap scratch file
    solnTable = Obit.TableUnref (solnTable)
    #
    return tabVer
    # end ResidCal

# Define Soln2Cal input dictionary
Soln2CalInput={'structure':['Soln2Cal',[('InData','Input OTF'),
                                        ('soln','input soln table version'),
                                        ('oldCal','input cal table version, -1=none'),
                                        ('newCal','output cal table')]],
               'InData':None, 'soln':0, 'oldCal':-1, 'newCal':0}
def Soln2Cal (err, input=Soln2CalInput):
    """ Apply a Soln (solution) table to a Cal (calibration) table.
    
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData = Python input OTF to calibrate
    soln   = Soln table version number to apply, 0-> high
    oldCal = input Cal table version number, -1 means none, 0->high
    newCal = output Cal table version number, 0->new
    """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    soln    = input["soln"]
    oldCal  = input["oldCal"]
    newCal  = input["newCal"]
    # Checks
    if not PIsA(inData):
        raise TypeError,'Soln2Cal: Bad input OTF'
    if err.isErr: # existing error?
        return 
    # Default table versions
    if soln == 0:
        soln = Obit.OTFGetHighVer(inData.me, "OTFSoln")
    if oldCal == 0:
        oldCal = Obit.OTFGetHighVer(inData.me, "OTFCal")
    if oldCal == 0:   # Must not be one
        oldCal = -1
    if newCal == 0:
        newCal = oldCal + 1;
    #
    # Specify desired tables 
    inInfo = inData.List
    inInfo.set("solnUse", soln)
    inInfo.set("calIn",   oldCal)
    inInfo.set("calOut",  newCal)
    #
    print "Soln2Cal: Soln ",soln," applied to cal",oldCal," write cal",newCal
    # Update calibration return Obit object
    calTable = Obit.OTFSoln2Cal (inData.me, inData.me, err.me);
    #
    # error check
    OErr.printErrMsg(err, "Soln2Cal: Error updating calibration")
    #
    # Cleanup
    calTable = Obit.TableUnref(calTable)
    # end Soln2Cal

# Define Split input dictionary
SplitInput={'structure':['Split',[('InData','Input OTF'),
                                  ('OutData','Extant output OTF'),
                                  ('Average','if true (1) average in frequency'),
                                  ('gainUse','cal. table version, -1=none'),
                                  ('flagVer','flag table version, -1=none')]],
            'InData':None, 'OutData':None, 'average':0, 'gainUse':-1, 'flagVer':-1};
def Split (err, input=SplitInput):
    """ Select and calibrate an OTF writing a new one.

    Applies calibration and editing/selection to inData and writes outData.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData  = input Python OTF to calibrate
    OutData = output Python OTF, must be previously defined
    Average = if true average in frequency
    gainUse = version number of prior table (Soln or Cal) to apply, -1 is none
    flagVer = version number of flagging table to apply, -1 is none
    """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    outData = input["OutData"]
    Average = input["Average"]
    gainUse = input["gainUse"]
    flagVer = input["flagVer"]
    # Checks
    if not PIsA(inData):
        raise TypeError,'Split: Bad input OTF'
    if not PIsA(outData):
        raise TypeError,'Split: Bad output OTF'
    if err.isErr: # existing error?
        return 
    #
    # Apply calibration as requested 
    dim[0] = 1; dim[1] = 1
    inInfo = inData.List
    InfoList.PAlwaysPutInt (inInfo, "flagVer", dim, [flagVer])
    InfoList.PAlwaysPutInt (inInfo, "gainUse", dim, [gainUse])
    if gainUse >= 0:
        itemp = 1
    else:
        itemp=0
    InfoList.PAlwaysPutInt (inInfo, "doCalib", dim, [itemp])
    doCalSelect = (gainUse >= 0) or (flagVer>0)
    doClaSelect = True
    InfoList.PAlwaysPutBoolean (inInfo, "doCalSelect", dim, [doCalSelect])
    #
    # Copy image from inData to outData 
    if average:
        Obit.OTFAver  (inData.me, outData.me, err.me);   # Average in frequency
    else:
        Obit.OTFCopy(inData.me, outData.me, err.me);     # Simple copy
    #
    # error check
    OErr.printErrMsg(err, "Split: Error copying data")
    #
    # end Split

# Define Concat input dictionary
ConcatInput={'InData':None, 'OutData':None};
def Concat (err, input=ConcatInput):
    """ Concatenates OTFs.

    Applies Copies InData to the end of  OutData.
    The files must be compatable (not checked)
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData  = Python input OTF to calibrate
    OutData = Python output OTF, must be previously defined
    """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    outData = input["OutData"]
    #
    # Copy image from inData to outData 
    Obit.OTFConcat(inData.me, outData.me, err.me);     # Simple append
    #
    # error check
    OErr.printErrMsg(err, "Concat: Error copying data")
    #
    # end Concat

# Define Image input dictionary
ImageInput={'structure':['Image',[('InData','Input OTF'),
                                  ('OutName','Output image file name'),
                                  ('OutWeight','Output gridding weight file name'),
                                  ('Disk','disk number for output image file'),
                                  ('ra','center RA (deg)'),
                                  ('dec','center Dec (deg)'),
                                  ('nx','number of pixels in x = RA'),
                                  ('ny','number of pixels in y = de'),
                                  ('xCells','Cell spacing in x (asec)'),
                                  ('yCells','Cell spacing in y (asec)'),
                                  ('minWt','minimum summed weight in gridded image wrt max '),
                                  ('Clip','flag data with abs. value grweater than Clip'),
                                  ('ConvType','Conv. fn type, 0=pillbox,3=Gaussian,4=exp*sinc,5=Sph wave'),
                                  ('ConvParm','Conv. fn parameters'),
                                  ('gainUse','cal. table version, -1=none'),
                                  ('doFilter','Filter out of band noise? [True]'),
                                  ('doBeam','Convolved Beam image desired? [def True]'),
                                  ('Beam','Instrumental response beam [def None]'),
                                  ('Wt','Image to save gridding weight array [def None], overrides OutWeight'),
                                  ('flagVer','flag table version, -1=none')]],
            'InData':None, 'OutName':None, 'OutWeight':None, 'Disk':1,
            'ra':0.0, 'dec':0.0, 'nx':100, 'ny':100,
            'xCells':0.001, 'yCells':0.001, 'minWt':0.01, 'Clip':1.0e19,
            'ConvParm':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'ConvType':3, 'doFilter':True,
            'gainUse':-1, 'doBeam':True, 'Beam':None, 'Wt':None, 'flagVer':-1};
def makeImage (err, input=ImageInput):
    """ Image an OTF.

    Data is convolved and resampled onto the specified grid.
    Image is created and returned on success.
    err     = Python Obit Error/message stack
    input   = input parameter dictionary

    Input dictionary entries:
    InData  = input Python OTF to image
             Additional optional parameters on InfoList member:
            "beamNx" int scalar "X" size of Beam (pixels)
            "beamNy" int scalar "Y" size of Beam(pixels)
            "doScale" bool scalar If true, convolve/scale beam [def True]
                      Only use False if providing a dirty beam which already
                      includes the effects of gridding.
            "deMode"  bool scalar Subtract image mode from image? [def False]
            "deBias"  bool scalar Subtract calibration bias from image? [def False]
                      Note, this doesn't really work the way you would like
    OutName = name of output image file
    OutWeight =  Output gridding weight file name
    Disk    = disk number for output image file
    ra      = center RA (deg)
    dec     = center Dec (deg)
    nx      = number of pixels in "x" = RA
    ny      = number of pixels in 'Y' = dec
    xCells  = Cell spacing in x (asec)
    yCells  = Cell spacing in y (asec)
    minWt   = minimum summed weight in gridded image wrt max [def 0.01]
    Clip    = data values with abs. value larger are set zero weight
    ConvType= Convolving function Type 0=pillbox,3=Gaussian,4=exp*sinc,5=Sph wave
    ConvParm= Convolving function parameters depends on ConvType
      Type 2 = Sinc, (poor function - don't use)
        Parm[0] = halfwidth in cells,
        Parm[1] = Expansion factor
      Type 3 = Gaussian,
        Parm[0] = halfwidth in cells,[def 3.0]
        Parm[1] = Gaussian with as fraction or raw beam [def 1.0]
      Type 4 = Exp*Sinc
        Parm[0] = halfwidth in cells, [def 2.0]
        Parm[1] = 1/sinc factor (cells) [def 1.55]
        Parm[2] = 1/exp factor (cells) [def 2.52]
        Parm[3] = exp power [def 2.0]
      Type 5 = Spherodial wave 
        Parm[0] = halfwidth in cells [def 3.0]
        Parm[1] = Alpha [def 5.0]
        Parm[2] = Expansion factor [not used]

    gainUse = version number of prior table (Soln or Cal) to apply, -1 is none
    flagVer = version number of flagging table to apply, -1 is none
    doFilter= Filter out of band noise?
    doBeam  = Beam convolved with convolving Fn image desired? [def True]
    Beam    = Actual instrumental Beam to use, else Gaussian [def None]
    Wt      = Image to save gridding weight array [def None], overrides OutWeight
    """
    ################################################################
    # Get input parameters
    inData  = input["InData"]
    outname = input["OutName"]
    outwt   = input["OutWeight"]
    disk    = input["Disk"]
    RA      = input["ra"]
    Dec     = input["dec"]
    nx      = input["nx"]
    ny      = input["ny"]
    xCells  = input["xCells"]
    yCells  = input["yCells"]
    minWt   = input["minWt"]
    Clip    = input["Clip"]
    ConvType= input["ConvType"]
    ConvParm= input["ConvParm"]
    gainUse = input["gainUse"]
    flagVer = input["flagVer"]
    doBeam  = input["doBeam"]
    Beam    = input["Beam"] 
    Wt      = input["Wt"]
    doFilter=  input["doFilter"]
   # Default table versions (the Obit routines will do this as well)
    if gainUse == 0:   # Get highest numbered OTFCal table
        gainUse = Obit.OTFGetHighVer(inData.me, "OTFCal")
    if gainUse == 0:   # Doesn't seem to be one, try OTFSoln
        gainUse = Obit.OTFGetHighVer(inData.me, "OTFSoln")
    if gainUse == 0:   # Must not be one
        gainUse = -1
    if Beam:
        lbeam = Beam   # Actual beam given
    else:
        lbeam = Image.Image("NoBeam")  # Not given
    #print 'Image: using cal table ',gainUse
    # Checks
    if not PIsA(inData):
        raise TypeError,'Image: Bad input OTF'
    if err.isErr: # existing error?
        return None
    #
    # Set imaging/calibration parameters
    inInfo = inData.List
    dim[0] = 1; dim[1] = 1
    InfoList.PAlwaysPutFloat(inInfo, "RA",     dim, [RA])
    InfoList.PAlwaysPutFloat(inInfo, "Dec",    dim, [Dec])
    InfoList.PAlwaysPutInt(inInfo, "nx",       dim, [nx])
    InfoList.PAlwaysPutInt(inInfo, "ny",       dim, [ny])
    InfoList.PAlwaysPutInt(inInfo, "ConvType", dim, [ConvType])
    xCells = -abs(xCells) # RA goes backwards
    InfoList.PAlwaysPutFloat(inInfo, "xCells",      dim, [xCells])
    InfoList.PAlwaysPutFloat(inInfo, "yCells",      dim, [yCells])
    InfoList.PAlwaysPutFloat(inInfo, "minWt",       dim, [minWt])
    InfoList.PAlwaysPutFloat(inInfo, "Clip",        dim, [Clip])
    docal = gainUse >= 0
    InfoList.PAlwaysPutBoolean(inInfo, "doCalSelect", dim, [True])
    InfoList.PAlwaysPutBoolean(inInfo, "doFilter",    dim, [doFilter])
    dim[0] = len(ConvParm)
    InfoList.PAlwaysPutFloat(inInfo, "ConvParm",  dim, ConvParm)
    dim[0] = 1
    if docal:
        dcal = 1
    else:
        dcal = 0
    InfoList.PAlwaysPutInt(inInfo, "doCalib",  dim, [dcal])
    InfoList.PAlwaysPutInt(inInfo, "gainUse",  dim, [gainUse])
    InfoList.PAlwaysPutInt(inInfo, "flagVer",  dim, [flagVer])
    #
    # Create output image object from OTF 
    outImage = Obit.OTFUtilCreateImage (inData.me, err.me)
    #
    # Define output 
    Obit.ImageSetFITS(outImage, 2, disk, outname, blc, trc, err.me)
    Obit.ImagefullInstantiate (outImage, 2, err.me)
    # error check 
    OErr.printErrMsg(err, "Image: Error verifying output file")
    # Define Weight image
    if Wt:
        lWt = Wt   # Actual weight image given
    else:
        # Create output Weight image?
        if outwt:
            # Define output 
            lWt = Image.Image("Gridding Weight")
            Obit.ImageSetFITS(lWt.me, 2, disk, outwt, blc, trc, err.me)
            Obit.ImageClone(outImage, lWt.me, err.me)
            Obit.ImagefullInstantiate (lWt.me, 2, err.me)
            # error check 
            OErr.printErrMsg(err, "Image: Error verifying weight file")
            input["Wt"] = lWt  # save
        else:
            lWt = Image.Image("NoWt")  # Not given
    # Form image 
    Obit.OTFUtilMakeImage (inData.me, outImage, doBeam, lbeam.me, lWt.me, err.me);
    # error check 
    OErr.printErrMsg(err, "Image: Error imaging OTF")
    #
    # Wrap output image in a Python object
    out = Image.Image(" ")
    out.me = outImage
    #
    return out
    # end makeImage

def SelfCal (err, ImageInp=ImageInput, CleanInp=None,
             ResidCalInp=ResidCalInput, Soln2CalInp=Soln2CalInput, ):
    """ Self calibrate an OTF

    Image an OTF, optionally Clean, determine residual calibration,
    apply to Soln to Cal table.  If the Clean is done, then the CLEAN result is
    used as the model in the ResidCal, otherwise the dirty image from Image is.
    err         = Python Obit Error/message stack
    ImageInp    = input parameter dictionary for Image
    CleanInp    = input parameter dictionary for Clean, "None"-> no Clean requested
                  May be modified to point to the result of the Image step
    ResidCalInp = input parameter dictionary for ResidCal
                  Will be modified to give correct derived model image
    Soln2CalInp = input parameter dictionary for Soln2Cal
    """
    ################################################################
    #
    if err.isErr: # existing error?
        return 
    # Dirty Image
    image = makeImage(err, ImageInp)
    image.Open(1,err)
    image.Read(err)
    OErr.printErrMsg(err, "OTF:SelfCal: Error reading model")
    model     = image.FArray
    modelDesc = image.Desc
    image.Close(err)
    OErr.printErrMsg(err, "OTF:SelfCal: Error imaging OTF")
    # Clean if requested
    if CleanInp != "None":
        CleanObj  = CleanInp["CleanOTF"]       # Clean object
        resid = CleanObj.Clean                 # Copy image just produced
        Image.PCopy(image, resid, err)         # to clean(resid)        
        #input(CleanInp)
        CleanOTF.PClean (err, CleanInp)         # Do Clean
        OErr.printErrMsg(err, "OTF:SelfCal: Error Cleaning")
        CCimage = CleanObj.Clean                # Get Clean for image model
        CCTab    = CCimage.NewTable(1,"AIPS CC", 0, err)
        # Use all
        dim[0] = 1; dim[1] = 1
        InfoList.PAlwaysPutInt(CCTab.List, "BComp", dim, [1]);
        InfoList.PAlwaysPutInt(CCTab.List, "EComp", dim, [0]);
        OErr.printErrMsg(err, "OTF:SelfCal: Error getting CCTable")
        # Make model image, CC convolved with beam
        model = OTFUtil.PConvBeam (CCTab, CleanObj.Beam, model, err);
        OErr.printErrMsg(err, "OTF:SelfCal: Error making model")
    #
    # Residual calibration
    ResidCalInp['Model']     = model                # Set reference image array
    ResidCalInp['ModelDesc'] = modelDesc            # Set reference image array
    #input(ResidCalInp)
    ResidCal (err, ResidCalInp)                 # Determine calibration
    OErr.printErrMsg(err, "OTF:SelfCal: Error Calibrating")
    #
    # Apply solution to calibration table
    Soln2Cal (err, Soln2CalInp)
    #
    # end SelfCal

    
def PScratch (inOTF, err):
    """ Create a scratch file suitable for accepting the data to be read from inOTF

    A scratch OTF is more or less the same as a normal OTF except that it is
    automatically deleted on the final unreference.
    inOTF     = Python OTF object
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inOTF.Scratch(err)
    # end PScratch

def PZap (inOTF, err):
    """ Delete underlying files and the basic object.

    inOTF     = Python OTF object
    err       = Python Obit Error/message stack
    """
    ################################################################
    inOTF.Zap(err)
    # end PZap

def PRename (inOTF, err, newFITSName=None):
    """ Rename underlying files

    inOTF   = Python OTF object
    err       = Python Obit Error/message stack
    For FITS files:
    newFITSName = new name for FITS file
    """
    ################################################################
    inOTF.Rename(err,newFITSName=newFITSName)
    # end PRename

def PCopy (inOTF, outOTF, err):
    """ Make a deep copy of input object.

    Makes structure the same as inOTF, copies data, tables
    inOTF     = Python OTF object to copy
    outOTF    = Output Python OTF object, must be defined
    err       = Python Obit Error/message stack
    """
    ################################################################
    inOTF.Copy (outOTF, err)
    # end PCopy

def PClone (inOTF, outOTF, err):
    """ Make a copy of a object but do not copy the actual data

    This is useful to create an OTF similar to the input one.
    inOTF   = Python OTF object
    outOTF  = Output Python OTF object, must be defined
    err     = Python Obit Error/message stack
    """
    ################################################################
    inOTF.Clone (outOTF, err)
    # end PClone

def PConcat (inOTF, outOTF, err):
    """ Copy data from inOTF to the end of outOTF

    inOTF   = Python OTF object
    outOTF  = Output Python OTF object, must be defined
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not PIsA(outOTF):
        raise TypeError,"outOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.OTFConcat (inOTF.cast(myClass), outOTF.cast(myClass), err.me)
    # end PConcat

def PNewOTFTable (inOTF, access, tabType, tabVer, err,
                  numDet=1, numPoly=0, numParm=0):
    """ Return the specified associated table
    
    inOTF     = Python OTF object
    access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
    tabType   = Table type, e.g. "OTFSoln"
    tabVer    = table version, if > 0 on input that table returned,
                if 0 on input, the highest version is used.
    err       = Python Obit Error/message stack
        Optional parameters, values only used if table created
        numDet    = Number of Detectors (OTFCal, OTFSoln, OTFScanData)
        numPoly   = Number of polynomial terms (OTFCal, OTFSoln)
        numParm   = Number of model parameters (OTFModel)
    """
    ################################################################
    return inOTF.NewTable (access, tabType, tabVer, err, \
                           numDet=numDet, numPoly=numPoly, numParm=numParm)
    # end PNewOTFTable
    

def POpen (inOTF, access, err):
    """ Open an image persistent (disk) form

    Returns 0 on success, else failure
    inOTF     = Python OTF object
    access    = access 1=READONLY, 2=WRITEONLY, 3=READWRITE
    err       = Python Obit Error/message stack
    """
    ################################################################
    return inOTF.Open (access, err)
    # end POpen

def PDirty (inOTF):
    """ Mark OTF as needing a header update to disk file

    inOTF     = Python OTF object
    """
    ################################################################
    inOTF.Dirty()
    # end PDirty

def PClose (inOTF, err):
    """ Close an image  persistent (disk) form

    inOTF     = Python OTF object
    err       = Python Obit Error/message stack
    """
    ################################################################
    inOTF.Close (err)
    # end PClose

def PZapTable (inOTF, tabType, tabVer, err):
    """ Destroy specified table

    inOTF     = Python OTF object
    tabType   = Table type, e.g. "OTFSoln"
    tabVer    = table version, integer
    err       = Python Obit Error/message stack
    """
    ################################################################
    inOTF.ZapTable (tabType, tabVer, err)
    # end PZapTable

def PCopyTables (inOTF, outOTF, exclude, include, err):
    """ Copy Tabels from one image to another

    inOTF     = Python OTF object
    outOTF    = Output Python OTF object, must be defined
    exclude   = list of table types to exclude (list of strings)
                has priority
    include   = list of table types to include (list of strings)
    err       = Python Obit Error/message stack
    """
    ################################################################
    inOTF.CopyTables (outOTF, exclude, include, err)
    # end PCopyTables

def PUpdateTables (inOTF, err):
    """ Update any disk resident structures about the current tables

    inOTF     = Python OTF object
    err       = Python Obit Error/message stack
    """
    ################################################################
    inOTF.UpdateTables (err)
    # end PUpdateTables

def PFullInstantiate (inOTF, access, err):
    """ Fully instantiate an OTF by opening and closing

    return 0 on success, else failure
    inOTF     = Python OTF object
    access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return None
    #
    return Obit.OTFfullInstantiate (inOTF.cast(myClass), access, err.me)
    # end PFullInstantiate

def PSetTarget (inOTF, Target, Flux, RA, Dec, err):
    """ Set target flux density and position

    inOTF     = Python OTF object
    Target    = Target name
    Flux      = Target Flux density
    RA        = RA in deg at mean equinox and epoch
    Dec       = Dec in deg at mean equinox and epoch
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    return Obit.OTFSetTarget (inOTF.cast(myClass), Target, Flux, RA, Dec, err.me)
    # end PSetTarget

def PGetList (inOTF):
    """ Return the member InfoList

    returns InfoList
    inOTF     = Python OTF object
    """
    ################################################################
    return inOTF.List
    # end PGetList

def PGetTableList (inOTF):
    """ Return the member tableList

    returns tableList
    inOTF   = Python OTF object
    """
    ################################################################
    return inOTF.TableList
    # end PGetTableList

def PHeader (inOTF, err):
    """ Print data descriptor

    inOTF      = Python Obit OTF object
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF data"
    #
    # Fully instantiate
    PFullInstantiate (inOTF, READONLY, err)
    # File info
    if inOTF.FileType=="AIPS":
        print "AIPS OTF Data Name: %12s Class: %6s seq: %d8 disk: %4d" % \
              (inOTF.Aname, inOTF.Aclass, inOTF.Aseq, inOTF.Disk)
    elif inOTF.FileType=="FITS":
        print "FITS OTF Data Disk: %5d File Name: %s " % \
              (inOTF.Disk, inOTF.FileName)
    # print in OTFDesc
    OTFDesc.PHeader(inOTF.Desc)
    # Tables
    TL = inOTF.TableList
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
    
def PGetDesc (inOTF):
    """ Return the member OTFDesc

    returns OTFDesc as a Python Dictionary
    inOTF   = Python OTF object
    """
    ################################################################
    return inOTD.Desc
    # end PGetDesc

def PUpdateDesc (inOTF, err, Desc=None):
    """ Update external representation of descriptor

    inOTF   = Python OTF object
    err     = Python Obit Error/message stack
    Desc    = OTF descriptor, if None then use current descriptor
    """
    ################################################################
    inOTF.UpdateDesc (err, Desc=Desc)
    # end PUpdateDesc

def POTFInfo (inOTF, err):
    """ Get file info for extant OTF data object

    Fills in information on object, useful for scratch files
    inOTF   = Python OTF object
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOTF):
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # file info
    info = Obit.OTFInfo (inOTF.cast(myClass), err.me);
    if err.isErr:
        OErr.printErrMsg(err, "Error creating scratch file")

    if info["type"]=="FITS":
        inOTF.FileType = 'FITS'
        inOTF.Disk   = info["disk"]
        inOTF.Otype  = "OTF"
        inOTF.FileName = info["filename"]
        inOTF.Fname    = info["filename"]
    # end POTFInfo

def PGetRecBuf (inOTF):
    return inOTF.RecBuf

def PGetHighVer (inOTF, tabType):
    """ Get highest version number of a specified Table

    returns highest tabType version number, 0 if none.
    inOTF     = Python OTF object
    tabType   = Table type, e.g. "OTFSoln"
    """
    ################################################################
    return  inOTF.GetHighVer (tabType)
    # end PGetHighVer

def PIsScratch (inOTF):
    """ Tells if OTF is a scratch object

    return true, false (1,0)
    inOTF   = Python OTF object
    """
    ################################################################
    return inOTF.IsScratch ()
    # end PIsScratch

def PIsA (inOTF):
    """ Tells if input really a Python Obit OTF

    return true, false (1,0)
    inOTF     = Python OTF object
    """
    ################################################################
    return Obit.OTFIsA(inOTF.cast(myClass))
    # end PIsA

def PGetName (inOTF):
    """ Tells OTF object name (label)

    returns name as character string
    inOTF     = Python OTF object
    """
    ################################################################
    return inOTF.GetName()
    # end PGetName

