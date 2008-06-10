# $Id: TableDesc.py,v 1.3 2007/12/16 20:26:15 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2005
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
 
# Python shadow class to ObitTableDesc class
import Obit, InfoList, OErr

class TableDescPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.TableDesc_me_set(self.this,value)
            return
        if name=="Dict":
            return PSetDict(self,value)
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.TableDesc_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="Dict":
            return PGetDict(self)
        raise AttributeError,str(name)
    def __repr__(self):
        return "<C TableDesc instance>"
class TableDesc(TableDescPtr):
    """ Python Obit Image descriptor class

    This contains information about the structure of a table

    Image Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List (readonly)
    Dict      - (virtual) Python dictionary with contents of descriptor
                Member Dict
    """
    def __init__(self, name) :
        self.this = Obit.new_TableDesc(name)
    def __del__(self):
        Obit.delete_TableDesc(self.this)

def PGetDict (inTD):
    """ Returns the contents of an TableDesc as a Python Dictionary

    returns dictionary
    inTD = Python TableDesc to read
    """
    ################################################################
    # Checks
    if not PIsA(inTD):
        raise TypeError,"inTD MUST be a Python Obit TableDesc"
    #
    return Obit.TableDescGetDict(inTD.me)
    # end PGetDict


def PSetDict (inTD, inDict):
    """ Copies the contents of a Python Dictionary to an TableDesc

    Only write descriptive, not structural values
    It's best if this was created by PGetDict.
    inTD   = Python TableDesc to update
    inDict = Python dictionary with values
    """
    ################################################################
    # Checks
    if not PIsA(inTD):
        raise TypeError,"inTD MUST be a Python Obit TableDesc"
    #
    Obit.TableDescSetDict(inTD.me, inDict)
    # end PSetDict

def PDef (inDict):
    """  Create TableDesc from the contents of a Python Dictionary

    Returns new Table Descriptor
    inDict = Python dictionary with values, must be in the form produced
             by PGetDict
    """
    ################################################################
    #
    outTD = TableDesc(inDict["Table name"])
    outTD.me = Obit.TableDescDef(inDict)
    # Check
    if len(outTD.Dict) <= 0:
        raise RuntimeError,"Failed to create valid Table Descriptor"
    return outTD
    # end PDef

def PGetList (inDesc):
    """  Get InfoList from TableDesc

    returns InfoList
    inDesc  = Python Obit input TableDesc
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit TableDesc"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.TableDescGetList(inDesc.me)
    return out
    # end PGetList 

def PIsA (inTD):
    """ Tells if the input really is a Python Obit TableDesc

    returns true or false (1,0)
    inTD = Python TableDesc to test
    """
    ################################################################
     # Checks
    if inTD.__class__ != TableDesc:
        return 0
    #
    return Obit.TableDescIsA(inTD.me)
    # end  PIsA

def PUnref (inTD):
    """ Decrement reference count

    Decrement reference count which will destroy object if it goes to zero
    Python object stays defined.
    inTD   = Python TableDesc object
    """
    ################################################################
     # Checks
    if not PIsA(inTD):
        raise TypeError,"inTD MUST be a Python Obit TableDesc"

    inTD.me = Obit.TableDescUnref(inTD.me)
    # end PUnref


