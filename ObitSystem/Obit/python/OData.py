""" Python Obit base Data class

This class Is the base class for Obit data classes
An ObitData is the front end to a persistent disk resident structure.
There maybe (usually are) associated tables which either describe
the data or contain calibration and/or editing information.
Both FITS (as Tables) and AIPS cataloged data are supported.
Most access to data is through functions as the volume of the data is
inappropriate to be processed directly in python.

OData Members with python interfaces:

=====  =======================
exist  True if object previously existed prior to object creation
List   used to pass instructions to processing
=====  =======================
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007-2019
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
 
# Python shadow class to ObitOData class
from __future__ import absolute_import
import OErr, InfoList, Table, History, Image, UV
import Obit, _Obit, TableList,  string

# class name in C
myClass = "ObitData"
 
class OData(Obit.OData):
    """
    Python ObitData (OData) class
    
    OData Members with python interfaces:

    InfoList  - used to pass instructions to processing
    """
    def __init__(self, name) :
        super(OData, self).__init__()
        Obit.CreateOData(self.this, name)
        self.myClass = myClass
    def __del__(self, DeleteOData=_Obit.DeleteOData):
        if _Obit!=None:
            DeleteOData(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            if value==None:
                raise TypeError("None given for OData object")
            # Out with the old
            if self.this!=None:
                Obit.ODataUnref(Obit.OData_Get_me(self.this))
            # In with the new
            Obit.OData_Set_me(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, OData):
            return "Bogus dude",str(self.__class__)
        if name == "me" : 
            return Obit.OData_Get_me(self.this)
        # Functions to return members
        if name=="List":
            if not self.ODataIsA():
                raise TypeError("input MUST be a Python Obit OData")
            out    = InfoList.InfoList()
            out.me = Obit.ODataGetList(self.cast(myClass).me)
            return out
        if name=="TableList":
            if not self.ODataIsA():
                raise TypeError("input MUST be a Python Obit OData")
            out    = InfoList.InfoList()
            out.me = Obit.ODataGetTableList(self.cast(myClass).me)
            return out
        if name=="Desc":
            raise AttributeError("Stubbed, virtual class")
        raise AttributeError(str(name))  # Unknown
    def __repr__(self):
        if self.__class__ != OData:
            return
        return "<C OData instance> " + Obit.ODataGetName(self.me)
    def cast(self, toClass):
        """ Casts object to specified class

        * self     = object whose cast pointer is desired
        * toClass  = Class string to cast to
        """
        # by type, first UV
        try:
            if Obit.UVIsA(self.me)!=0:
                if (toClass=='ObitData'):
                    recast = OData("None")
                    recast.me = Obit.UVCastData(self.me)
                else:
                    recast = self
                return recast
        except:
            pass
        # Try Image
        try:
            if Obit.ImageIsA(self.me)!=0:
                if (toClass=='ObitData'):
                    recast = OData("None")
                    recast.me = Obit.ImageCastData(self.me)
                else:
                    recast = self
                return recast
        except:
                pass
        # Try ImageMF
        try:
            if Obit.ImageMFIsA(self.me)!=0:
                if (toClass=='ObitImage'):
                    recast = Image.Image("None")
                    recast.me = Obit.ImageMFCastImage(self.me)
                else:
                    recast = OData("None")
                    recast.me = Obit.ImageMFCastData(self.me)
                return recast
        except Exception as exception:
                print ('Fail',exception)
                pass
         # Already OData?
        try:
            if Obit.ODataIsA(self.me)!=0:
                recast = self
                return recast
        except:
                pass
        # If it gets here bail
        raise TypeError("Unknown type to cast")
    # end cast
            
    def Zap (self, err):
        """
        Delete underlying files and the basic object.

        * self      = Python OData object
        * err       = Python Obit Error/message stack
        """
        inOData = self
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        Obit.ODataZap (inOData.cast(myClass).me, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Zapping OData data")
        # end Zap
        
    def Open (self, access, err):
        """
        Open a OData data persistent (disk) form
        
        Returns 0 on success

        * self   = Python OData object
        * access = access READONLY (1), WRITEONLY (2), READWRITE(3)
        * err    = Python Obit Error/message stack
        """
        inOData = self
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        ret = Obit.ODataOpen(inOData.cast(myClass).me, access, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Opening OData data")
        return ret
        # end Open
        
    def Close (self, err):
        """
        Close a OData  persistent (disk) form
        
        Returns 0 on success

        * self      = Python OData object
        * err       = Python Obit Error/message stack
        """
        inOData = self
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        ret = Obit.ODataClose (inOData.cast(myClass).me, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Closing OData data")
        return ret
        # end Close
        
    def History (self, access, err):
        """
        Return the associated History

        * self      = Python OData object
        * access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
        * err       = Python Obit Error/message stack
        """
        out    = History.History("History", self.List, err)
        #out.me = newODataHistory (self, access, err)
        # Error?
        if err.isErr:
            OErr.printErrMsg(err, "Error getting History")
        return out
        # end History
        
    def Stalin (self, startr, endr, err):
        """
        Edit history
        
        return 0 on success, else failure

        * self   = input Python History
        * startr = first (1-rel) history record to delete
        * endr   = highest (1-rel) history record to delete, 0->to end
        * err    = Python Obit Error/message stack
        """
        hi = self.History (READWRITE, err)
        return History.PEdit (hi, startr, endr, err)
    # end Stalin

    def NewTable (self, access, tabType, tabVer, err, \
                  numOrb=0, numPCal=2, numIF=1, numPol=1, \
                  numTerm=0, numChan=1, numTones=1, numBand=1, \
                  numTabs=1, npoly=1, numCoef=5, noParms=0):
        """
        Return the specified associated table
        
        Table will be created if necessary.
        * self      = Python OData object
        * access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
        * tabType   = Table type, e.g. "AIPS AN"
        * tabVer    = table version, if > 0 on input that table returned,
          if 0 on input, the highest version is used.
        * err       = Python Obit Error/message stack

        Optional parameters, values only used if table created
        * numOrb    = Number of orbital parameters (AN)
        * numPCal   = Number of polarization parameters per IF (AN)
        * numIF     = Number of IFs (CD, FQ, SN, CL, BP, BL, SY, TY, CQ, CP, PD)
        * numPol    = Number of Stokes' (AN, SN, CD, CL, BP, BL, PC, SY, TY, GC, MC, IM, PD)
        * numTerm   = Number of terms in model polynomial (CL)
        * numChan   = Number of spectral channels (BP, CP, PD, FS)
        * numTomes  = Number of Phase cal tones (PC)
        * numTabs   = Number of ??? (GC)
        * numCoef   = Number of polynomial coefficents (NI)
        * numBand   = Number of Bands(?) (IM, GC, MC)
        * npoly     = number of polynomial terms (IM)
        * noParms   = Number of parameters in CC table model
        * maxis1-5  = Dimension of axes of IDI data matrix
        """
        inOData = self
        # Checks
        if not (self.cast(myClass)).ODataIsA():
            raise TypeError("input MUST be a Python Obit OData")
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        outTable    = Table.Table("None")
        tver = [int(tabVer)]  # Make sure long
        id  = inOData.cast(myClass).me
        if tabType=="AIPS AN":
            outTable.me = Obit.TableAN(id, tver, access, tabType, numIF, numOrb, numPCal, err.me)
        elif tabType=="AIPS SU":
            outTable.me = Obit.TableSU(id, tver, access, tabType, numIF, err.me)
        elif tabType=="AIPS NX":
            outTable.me = Obit.TableNX(id, tver, access, tabType, err.me)
        elif tabType=="AIPS FQ":
            outTable.me = Obit.TableFQ(id, tver, access, tabType, numIF, err.me)
        elif tabType=="AIPS FG":
            outTable.me = Obit.TableFG(id, tver, access, tabType, err.me)
        elif tabType=="AIPS SN":
            outTable.me = Obit.TableSN(id, tver, access, tabType, numPol, numIF, err.me)
        elif tabType=="AIPS CD":
            outTable.me = Obit.TableCD(id, tver, access, tabType, numIF, numPol, err.me)
        elif tabType=="AIPS CL":
            outTable.me = Obit.TableCL(id, tver, access, tabType,
                                  numPol, numIF, numTerm, err.me)
        elif tabType=="AIPS BP":
            outTable.me = Obit.TableBP(id, tver, access, \
                                  tabType, numPol, numIF, numChan, err.me)
        elif tabType=="AIPS PD":
            outTable.me = Obit.TablePD(id, tver, access, \
                                  tabType, numPol, numIF, numChan, err.me)
        elif tabType=="AIPS BL":
            outTable.me = Obit.TableBL(id, tver, access, tabType, numPol, numIF, err.me)
        elif tabType=="AIPS PC":
            outTable.me = Obit.TablePC(id, tver, access, tabType, numPol, numTones, err.me)
        elif tabType=="AIPS WX":
            outTable.me = Obit.TableWX(id, tver, access, tabType, err.me)
        elif tabType=="AIPS TY":
            outTable.me = Obit.TableTY(id, tver, access, tabType, numPol, numIF, err.me)
        elif tabType=="AIPS GC":
            outTable.me = Obit.TableGC(id, tver, access, tabType, \
                                numBand,  numPol, numTabs, err.me)
        elif tabType=="AIPS MC":
            outTable.me = Obit.TableMC(id, tver, access, tabType, numPol, numBand, err.me)
        elif tabType=="AIPS NI":
            outTable.me = Obit.TableNI(id, tver, access, tabType, numCoef, err.me)
        elif tabType=="AIPS PS":
            outTable.me = Obit.TablePS(id, tver, access, tabType, err.me)
        elif tabType=="AIPS CP":
            outTable.me = Obit.TableCP(id, tver, access, tabType, numIF, numChan,  err.me)
        elif tabType=="AIPS CQ":
            outTable.me = Obit.TableCQ(id, tver, access, tabType, numIF, err.me)
        elif tabType=="AIPS IM":
            outTable.me = Obit.TableIM(id, tver, access, tabType, \
                                  numPol, numBand, npoly, err.me)
        elif tabType=="AIPS CT":
            outTable.me = Obit.TableCT(id, tver, access, tabType, err.me)
        elif tabType=="AIPS OB":
            outTable.me = Obit.TableOB(id, tver, access, tabType, err.me)
        elif tabType=="AIPS OF":
            outTable.me = Obit.TableOF(id, tver, access, tabType, err.me)
        elif tabType=="AIPS CC":
            outTable.me = Obit.TableCC(id, tver, access, tabType, noParms, err.me)
        elif tabType=="AIPS VL":
            outTable.me = Obit.TableVL(id, tver, access, tabType, err.me)
        elif tabType=="AIPS VZ":
            outTable.me = Obit.TableVZ(id, tver, access, tabType, err.me)
        elif tabType=="AIPS MF":
            outTable.me = Obit.TableMF(id, tver, access, tabType, err.me)
        elif tabType=="AIPS FS":
            outTable.me = Obit.TableFS(id, tver, access, tabType, numChan, err.me)
            # IDI tables
        elif tabType=="IDI_ANTENNA":
            outTable.me = Obit.TableIDI_ANTENNA(id, tver, access, tabType,
                                           numIF, numPCal, err.me)
        elif tabType=="IDI_ARRAY_GEOMETRY":
            outTable.me = Obit.TableIDI_ARRAY_GEOMETRY(id, tver, access, tabType,
                                                  numIF, numOrb, err.me)
        elif tabType=="IDI_FREQUENCY":
            outTable.me = Obit.TableIDI_FREQUENCY(id, tver, access, tabType,
                                             numIF, err.me)
        elif tabType=="IDI_SOURCE":
            outTable.me = Obit.TableIDI_SOURCE(id, tver, access, tabType,
                                          numIF, err.me)
        elif tabType=="IDI_OData_DATA":
            outTable.me = Obit.TableIDI_OData_DATA(id, tver, access, tabType,
                                              numIF, maxis1, maxis2, maxis3, maxis4, maxis5,
                                              err.me)
        else:  # Generic
            ret = Obit.newODataTable (inOData.cast(myClass).me, access, tabType, tver, err.me)
            # Table and version returned in a list
            outTable.me = ret
        # Did it work? error may not be set
        if not outTable.IsA():
            OErr.printErrMsg(err, "Error getting OData data Table")
            raise RuntimeError("Failed to extract "+tabType+" table from "+inOData.GetName())
        # Error?
        if err.isErr:
            OErr.printErrMsg(err, "Error getting OData data Table")
        Table.PFullInstantiate (outTable, access, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error instantiating OData data Table")
        return outTable
    # end NewTable

    def Header (self, err):
        """
        Write image header on output
        
        Virtual

        * self   = Python Obit OData object
        * err    = Python Obit Error/message stack
        """
        # Stubbed
        raise RuntimeError("Header: Not Defined for virtual base class OData")
    # end Header
        
    def Info (self, err):
        """
        Get underlying data file info
        
        Virtual

        * self   = Python Obit OData object
        * err    = Python Obit Error/message stack
        """
        # Stubbed
        raise RuntimeError("Info: Not Defined for virtual base class OData")
        # end Info
        
    def ZapTable (self, tabType, tabVer, err):
        """
        Destroy specified table
        
        Returns 0 on success
        * self      = Python OData object
        * tabType   = Table type, e.g. "AIPS CC"
        * tabVer    = table version, integer
        * err       = Python Obit Error/message stack
        """
        inOData = self
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        # Open/Close to (re)Read header to be sure it's up to date
        inOData.Open(READONLY,err)
        inOData.Close(err)
        
        # delete table
        ret = Obit.ODataZapTable(inOData.cast(myClass).me, tabType, tabVer, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Zapping OData data Table")
        # Update header
        inOData.UpdateTables(err)
        return ret
        # end ZapTable

    def UpdateTables (self, err):
        """
        Update any disk resident structures about the current tables
        
        Returns 0 on success

        * self      = Python Image object
        * err       = Python Obit Error/message stack
        """
        inOData = self
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        ret = Obit.ODataUpdateTables (inOData.cast(myClass).me, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error updating Tables info")
        return ret
        # end UpdateTables
        
    def UpdateDesc (self, err, Desc=None):
        """
        Update any disk resident structures about descriptor
        
        Virtual

        * self      = Python OData object
        * err       = Python Obit Error/message stack
        * Desc      = Descriptor, if None then use current descriptor
                    Contents can be accessed throuth the Dict member
        """
        raise RuntimeError("UpdateDesc: Not Defined for virtual base class OData")
        # end UpdateDesc
        
    def Scratch (self, err):
        """
        Create a scratch file suitable for accepting the data to be read from self
        
        A scratch OData is more or less the same as a normal OData except that it is
        automatically deleted on the final unreference.

        * self      = Python OData object
        * err       = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        outOData    = OData("None")
        outOData.me = Obit.ODataScratch (self.cast(myClass).me, err.me);
        if err.isErr:
            OErr.printErrMsg(err, "Error creating scratch file")
            
        return outOData
    # end Scratch

    def Rename (self, err, newFITSName=None, \
                newAIPSName="            ", \
                newAIPSClass="      ", newAIPSSeq=0):
        """
        Rename underlying files

        * self   = Python OData object
        * err       = Python Obit Error/message stack

        For FITS files:
        * newFITSName = new name for FITS file
        
        For AIPS:
        * newAIPSName  = New AIPS Name (max 12 char) Blank => don't change.
        * newAIPSClass = New AIPS Class (max 6 char) Blank => don't change.
        * newAIPSSeq   = New AIPS Sequence number, 0 => unique value
        """
        ################################################################
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        if len(newAIPSName)>12:
            raise RuntimeError("New AIPS Name too long")
        if len(newAIPSClass)>6:
            raise RuntimeError("New AIPS Class too long")
        #
        # Set controls
        inInfo = self.List    # 
        dim = [1,1,1,1,1]
        InfoList.PAlwaysPutInt      (inInfo, "newSeq",       dim, [newAIPSSeq])
        dim[0] = 12
        InfoList.PAlwaysPutString  (inInfo, "newName",    dim, [newAIPSName])
        dim[0] = 6
        InfoList.PAlwaysPutString  (inInfo, "newClass",   dim, [newAIPSClass])
        if newFITSName:
            dim[0] = len(newFITSName)
            InfoList.PAlwaysPutString  (inInfo, "newFileName",dim, [newFITSName])
        # Rename
        Obit.ODataRename (self.cast(myClass).me, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error Renaming OData data")
    # end PRename

    def Copy (self, outOData, err):
        """
        Make a deep copy of input object.
        
        Makes structure the same as self, copies data, tables

        * self   = Python OData object to copy
        * outOData  = Output Python OData object, must be defined
        * err       = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        Obit.ODataCopy (self.cast(myClass).me, outOData.cast(myClass).me, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error copying OData data")
    # end PCopy

    def Clone (self, outOData, err):
        """
        Make a copy of a object but do not copy the actual data
        
        This is useful to create an OData similar to the input one.

        * self   = Python OData object
        * outOData  = Output Python OData object, must be defined
        * err    = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        Obit.ODataClone (self.cast(myClass).me, outOData.cast(myClass).me, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error cloning OData data")
    # end PClone

    def Dirty (self):
        """
        Mark OData as needing a header update to disk file

        * self     = Python OData object
        """
        ################################################################
        Obit.ODataDirty (self.cast(myClass).me)
    # end PDirty

    def CopyTables (self, outOData, exclude, include, err):
        """
        Copy Tables from one OData to another

        * self      = Python OData object
        * outOData     = Output Python OData object, must be defined
        * exclude   = list of table types to exclude (list of strings)
          has priority
        * include   = list of table types to include (list of strings)
        * err       = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        ret = Obit.ODataCopyTables  (self.cast(myClass).me, outOData.cast(myClass).me, \
                                     exclude, include, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error copying Tables")
        return ret
    # end PCopyTables

    def FullInstantiate (self, access, err):
        """
        Fully instantiate an OData by opening and closing
        
        return 0 on success, else failure

        * self   = Python OData object
        * access    = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
        * err       = Python Obit Error/message stack
        """
        ################################################################
        # Checks
        if not OErr.OErrIsA(err):
            raise TypeError("err MUST be an OErr")
        #
        ret = Obit.ODataFullInstantiate (self.cast(myClass).me, access, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error verifying OData data")
        return ret
    # end PfullInstantiate

    def GetHighVer (self, tabType):
        """
        Get highest version number of a specified Table
                returns highest tabType version number, 0 if none.

        * self   = Python OData object
        * tabType   = Table type, e.g. "OTFSoln"
        """
        ################################################################
        return Obit.ODataGetHighVer(self.cast(myClass).me, tabType)
    # end PGetHighVer

    def IsScratch (self):
        """
        Tells if OData is a scratch object
        
        return true, false (1,0)

        * self   = Python OData object
        """
        ################################################################
        return Obit.ODataisScratch(self.cast(myClass).me)
    # end PIsScratch

    def ODataIsA (self):
        """
        Tells if input really a Python Obit OData
        
        return True, False
        * self   = Python OData object
        """
        ################################################################
        # Allow derived types
        try:
            return Obit.ODataIsA(self.me)!=0
        except:
            return False
        else:
            return False
    # end ODataIsA

    def GetName (self):
        """
        Tells OData object name (label)
        
        returns name as character string

        * self   = Python OData object
        """
        ################################################################
        return Obit.ODataGetName(self.cast(myClass).me)
    # end GetName
    
# End of class member functions (i.e. invoked by x.func())


# Symbolic names for access codes
READONLY  = 1
WRITEONLY = 2
READWRITE = 3


