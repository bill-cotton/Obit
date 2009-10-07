""" Python Obit Table class

This class contains tabular data and allows access.
An ObitTable is the front end to a persistent disk resident structure.
Both FITS (as Tables) and AIPS cataloged data are supported.

Table Members with python interfaces:
InfoList  - used to pass instructions to processing

Table header keywords for specific table types are available in the keys
member of a Table after the table has been opened.  These will be updated
to disk when the table is closed.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004,2005,2007
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

# Python shadow class to ObitTable class
import Obit, OErr, InfoList, TableDesc

class TablePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.TableUnref(Obit.Table_me_get(self.this))
            # In with the new
            Obit.Table_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != Table:
            return
        if name == "me" : 
            return Obit.Table_me_get(self.this)
         # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="IOList":
            return PGetIOList(self)
        if name=="Desc":
            return PGetDesc(self)
        if name=="IODesc":
            return PGetIODesc(self)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != Table:
            return
        return "<C Table instance> " + Obit.TableGetName(self.me)
class Table(TablePtr):
    def __init__(self,name) :
        self.this = Obit.new_Table(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_Table(self.this)

    def Zap (self, err):
        """ Delete underlying files and the basic object.
        
        self      = Python Table object
        err       = Python Obit Error/message stack
        """
        PZap(self,err)
        # end Zap

    def Open (self, access, err):
        """ Open an table persistent (disk) form

        Specific table type keywords are written to the "keys" dict member
        self   = Python Table object
        access    = access READONLY (1), WRITEONLY (2), READWRITE(3)
        err       = Python Obit Error/message stack
        """
        POpen(self, access, err)
        # end Open

    def Close (self, err):
        """ Close an table  persistent (disk) form
        
        Specific table type keywords are written from the "keys" dict member
        self      = Python Table object
        err       = Python Obit Error/message stack
        """
        PClose (self, err)
        # end Close

    def ReadRow (self, rowno, err):
        """ Read a specified row in a table and returns as a python Dict

        self   = Python Image object
        rowno     = row number (1-rel) to read
        err    = Python Obit Error/message stack
        """
        return PReadRow (self, rowno, err)
        # end ReadFA

    def WriteRow (self, rowno, rowDict, err):
        """ Write an image  persistent (disk) form from a specified Dict

        Writes a single row
        self      = Python Image object
        rowno     = row number (1-rel) to write
        rowDict   = Python Dict of same form as returned by PReadRow
        err       = Python Obit Error/message stack
        """
        PWriteRow (self, rowno, rowDict, err)
        # end WriteRow
        # End Table class definitions

# Symbolic names for access codes
READONLY  = 1
WRITEONLY = 2
READWRITE = 3

def PZap (inTab, err):
    """ Destroy the persistent form of a Table

    inTab    = input Python Obit Table
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
   #
    Obit.TableZap (inTab.me, err.me)
    # end PZap 

def PCopy (inTab, outTab, err):
    """ Copy a Table including persistent forms

    inTab    = input Python Obit Table
    outTab   = extant output Python Obit Table
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    if not PIsA(outTab):
        raise TypeError,"outTab MUST be a Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.TableCopy(inTab.me, outTab.me, err.me)
    # end PCopy

def PClone (inTab, outTab):
    """ Copy the structure of a Table

    inTab    = input Python Table
    outTab   = extant output Python Obit Table or None
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    #
    if not outTab:
        out = Table("None")
    else:
        out = outTab
    out.me = Obit.TableClone(inTab.me, out.me)
    return out
    # end PClone


def PConcat (inTab, outTab, err):
    """ Copy row data from inTab to the end of outTab

    inTab    = input Python Obit Table
    outTab   = extant output Python Obit Table
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    if not PIsA(outTab):
        raise TypeError,"outTab MUST be a Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.TableConcat(inTab.me, outTab.me, err.me)
    # end PConcat


def PFullInstantiate (inTab, access, err):
    """ Open and close to fully instantiate

    return 0 on success, else failure
    inTab    = input Python Table
    access   = access code 1=READONLY, 2=WRITEONLY, 3=READWRITE
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    if err.isErr: # existing error?
        return None
    return Obit.TablefullInstantiate(inTab.me, access, err.me)
    # end PFullInstantiate

def POpen (inTab, access, err):
    """ Open a table persistent (disk) form

    Specific table type keywords are written to the "keys" dict member
    inTab     = Python Table object
    access    = access READONLY (1), WRITEONLY (2), READWRITE(3)
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.TableOpen(inTab.me, access, err.me)
    #
    # Get specific type keywords as dict
    tabtype = inTab.Desc.Dict["Table name"]
    if tabtype=="AIPS AN":
        inTab.keys = Obit.TableANGetHeadKeys(inTab.me)
    elif tabtype=="AIPS AT":
        inTab.keys = Obit.TableATGetHeadKeys(inTab.me)
    elif tabtype=="AIPS BL":
        inTab.keys = Obit.TableBLGetHeadKeys(inTab.me)
    elif tabtype=="AIPS BP":
        inTab.keys = Obit.TableBPGetHeadKeys(inTab.me)
    elif tabtype=="AIPS CC":
        inTab.keys = Obit.TableCCGetHeadKeys(inTab.me)
    elif tabtype=="AIPS CL":
        inTab.keys = Obit.TableCLGetHeadKeys(inTab.me)
    elif tabtype=="AIPS CQ":
        inTab.keys = Obit.TableCQGetHeadKeys(inTab.me)
    elif tabtype=="AIPS CT":
        inTab.keys = Obit.TableCTGetHeadKeys(inTab.me)
    elif tabtype=="AIPS FG":
        inTab.keys = Obit.TableFGGetHeadKeys(inTab.me)
    elif tabtype=="AIPS FQ":
        inTab.keys = Obit.TableFQGetHeadKeys(inTab.me)
    elif tabtype=="AIPS GC":
        inTab.keys = Obit.TableGCGetHeadKeys(inTab.me)
    elif tabtype=="AIPS IM":
        inTab.keys = Obit.TableIMGetHeadKeys(inTab.me)
    elif tabtype=="AIPS MC":
        inTab.keys = Obit.TableMCGetHeadKeys(inTab.me)
    elif tabtype=="AIPS MF":
        inTab.keys = Obit.TableMFGetHeadKeys(inTab.me)
    elif tabtype=="AIPS NI":
        inTab.keys = Obit.TableNIGetHeadKeys(inTab.me)
    elif tabtype=="AIPS NX":
        inTab.keys = Obit.TableNXGetHeadKeys(inTab.me)
    elif tabtype=="AIPS OB":
        inTab.keys = Obit.TableOBGetHeadKeys(inTab.me)
    elif tabtype=="AIPS OF":
        inTab.keys = Obit.TableOFGetHeadKeys(inTab.me)
    elif tabtype=="AIPS PC":
        inTab.keys = Obit.TablePCGetHeadKeys(inTab.me)
    elif tabtype=="AIPS PS":
        inTab.keys = Obit.TablePSGetHeadKeys(inTab.me)
    elif tabtype=="AIPS SN":
        inTab.keys = Obit.TableSNGetHeadKeys(inTab.me)
    elif tabtype=="AIPS SU":
        inTab.keys = Obit.TableSUGetHeadKeys(inTab.me)
    elif tabtype=="AIPS TY":
        inTab.keys = Obit.TableTYGetHeadKeys(inTab.me)
    elif tabtype=="AIPS VL":
        inTab.keys = Obit.TableVLGetHeadKeys(inTab.me)
    elif tabtype=="AIPS VZ":
        inTab.keys = Obit.TableVZGetHeadKeys(inTab.me)
    elif tabtype=="AIPS WX":
        inTab.keys = Obit.TableWXGetHeadKeys(inTab.me)
    # end POpen

def PDirty (inTable):
    """ Mark Table as needing a header update to disk file

    inTable     = Python Table object
    """
    ################################################################
    # Checks
    if not PIsA(inTable):
        raise TypeError,"inTable MUST be a Python Obit Table"
    #
    Obit.TableDirty (inTable.me)
    # end PDirty

def PClose (inTab, err):
    """ Close a table  persistent (disk) form

    Specific table type keywords are written from the "keys" dict member
    inTab     = Python Table object
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    if err.isErr: # existing error?
        return
    #
    # Set specific type keywords from dict
    tabtype = inTab.Desc.Dict["Table name"]
    if tabtype=="AIPS AN" and inTab.keys:
        Obit.TableANSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS AT" and inTab.keys:
        Obit.TableATSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS BL" and inTab.keys:
        Obit.TableBLSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS BP" and inTab.keys:
        Obit.TableBPSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS CC" and inTab.keys:
        Obit.TableCCSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS CL" and inTab.keys:
        Obit.TableCLSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS CQ" and inTab.keys:
        Obit.TableCQSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS CT" and inTab.keys:
        Obit.TableCTSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS FG" and inTab.keys:
        Obit.TableFGSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS FQ" and inTab.keys:
        Obit.TableFQSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS GC" and inTab.keys:
        Obit.TableGCSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS IM" and inTab.keys:
        Obit.TableIMSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS MC" and inTab.keys:
        Obit.TableMCSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS MF" and inTab.keys:
        Obit.TableMFSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS NI" and inTab.keys:
        Obit.TableNISetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS NX" and inTab.keys:
        Obit.TableNXSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS OB" and inTab.keys:
        Obit.TableOBSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS OF" and inTab.keys:
        Obit.TableOFSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS PC" and inTab.keys:
        Obit.TablePCSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS PS" and inTab.keys:
        Obit.TablePSSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS SN" and inTab.keys:
        Obit.TableSNSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS SU" and inTab.keys:
        Obit.TableSUSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS TY" and inTab.keys:
        Obit.TableTYSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS VL" and inTab.keys:
        Obit.TableVLSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS VZ" and inTab.keys:
        Obit.TableVZSetHeadKeys(inTab.me, inTab.keys)
    elif tabtype=="AIPS WX" and inTab.keys:
        Obit.TableWXSetHeadKeys(inTab.me, inTab.keys)
    #
    Obit.TableClose (inTab.me, err.me)
    # end PClose

def PReadRow (inTab, rowno, err):
    """ Read a specified row in a table and returns as a python Dict

    Dict has keys:
       "Table name"  to give the name of the table
       Field named   (column labels)
    data are returned as a list of the field data type.
    inTab     = Python Table object
    rowno     = row number (1-rel) to read
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return None
    #
    return Obit.TableReadRow (inTab.me, rowno, err.me)
    # end PReadRow

def PWriteRow (inTab, rowno, rowDict, err):
    """ Write an image  persistent (disk) form from a specified Dict

    Writes a single row
    inTab     = Python Table object
    rowno     = row number (1-rel) to write
    rowDict   = Python Dict of same form as returned by PReadRow
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Python Obit Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.TableWriteRow (inTab.me, rowno, rowDict, err.me)
    OErr.printErr(err)
    # end PWriteRow

def PUnref (inTab):
    """ Decrement reference count

    Decrement reference count which will destroy object if it goes to zero
    Python object stays defined.
    inTab   = Python Table object
    """
    ################################################################
     # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Python Obit Table"

    inTab.me = Obit.TableUnref(inTab.me)
    # end PUnref

def PGetList (inTab):
    """ Return the InfoList from a Table
    
    returns InfoList
    inTab    = input Python Table
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.TableGetList(inTab.me)
    return out
    # end PGetList

def PGetIOList (inTab):
    """ Return the InfoList from a Table's IO member
    
    returns InfoList from IO member (disk resident version)
    if the IO member is not defined a None is returned.
    For most reliable results, this routine should be called when
    the table is opened with Write allowed.
    inTab    = input Python Table
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    #
    out    = InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.TableGetIOList(inTab.me)
    return out
    # end PGetIOList

def PGetDesc (inTab):
    """ Return the TableDesc from a Table
    
    returns TableDesc
    inTab    = input Python Table
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    #
    out    = TableDesc.TableDesc("TableDesc")
    out.me = Obit.TableDescUnref(out.me)
    out.me = Obit.TableGetDesc(inTab.me)
    return out
    # end PGetDesc

def PGetIODesc (inTab):
    """ Return the TableDesc from a Table's IO member
    
    returns TableDesc from IO member (disk resident version)
    if the IO member is not defined a None is returned.
    For most reliable results, this routine should be called when
    the table is opened with Write allowed.
    inTab    = input Python Table
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    #
    out    = TableDesc.TableDesc("TableDesc")
    out.me = Obit.TableDescUnref(out.me)
    out.me = Obit.TableGetIODesc(inTab.me)
    return out
    # end PGetDesc

def PGetVer (inTab):
    """ Get table version number

    returns table version number
    inTab    = input Python Table
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    #
    return Obit.TableGetVer(inTab.me)
    # end PGetVer

def PIsA (inTab):
    """ Tells if object thinks it's a Python Obit Table

    return true, false (1,0)
    inTab    = input Python Table
    """
    ################################################################
    # Checks
    if inTab.__class__ != Table:
        return 0
    #
    try:
        return Obit.TableIsA(inTab.me)
    except:
        return False
    # end PIsA

def PGetName (inTab):
    """ Returns object name (label)

    return name string
    inTab    = input Python Table
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    #
    return Obit.TableGetName(inTab.me)
    # end PGetName

def PSort (inTab, colName, desc, err):
    """ Sort a table of a column

    inTab    = input Python Obit Table to sort
    colName  = Column name (e.g. "Time")
    desc     = if true sort in descending order, else ascending
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTab):
        raise TypeError,"inTab MUST be a Table"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    #
    Obit.TableUtilSort(inTab.me, colName, desc, err.me)
    # end PSort

