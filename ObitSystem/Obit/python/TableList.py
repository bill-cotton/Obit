# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2005,2019
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
 
# Python shadow class to ObitTableList class
from __future__ import absolute_import
import Obit, _Obit, InfoList, OErr

class TableList(Obit.TableList):
    """ Python Obit TableList class

    This contains information about the Tables associated with an image or dataset

    Image Members with python interfaces:
    List      - (virtual) Python list of table names and numbers
    """
    def __init__(self, name) :
        super(TableList, self).__init__()
        Obit.CreateTableList(self.this, name)
    def __del__(self, DeleteTableList=_Obit.DeleteTableList):
        DeleteTableList(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.TableList_Set_me(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.TableList_Get_me(self.this)
        raise AttributeError(str(name))
    def __repr__(self):
        return "<C TableList instance>"

def PGetList (inTL, err):
    """ Returns the contents of an TableList as a Python list

    returns list
    inTL = Python TableList to read
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTL):
        raise TypeError("inTL MUST be a Python Obit TableList")
    if err.isErr:
        return None  # existing error 
    #
    return Obit.TableListGetList(inTL.me, err.me)
    # end PGetList

def PPrint (inTL, err):
    """ Print the contents of an TableList on stderr

    inTL = Python TableList to read
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTL):
        raise TypeError("inTL MUST be a Python Obit TableList")
    if err.isErr:
        return None  # existing error 
    #
    return Obit.TableListPrint(inTL.me, err.me)
    # end PPrint

def PCheck (inTL, err):
    """ Check the contents of an TableList

    Any errors lister on err
    inTL = Python TableList to check
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTL):
        raise TypeError("inTL MUST be a Python Obit TableList")
    if err.isErr:
        return   # existing error 
    #
    Obit.TableListCheck(inTL.me, err.me)
    # end PCheck

def PGetHigh (inTL, tabType):
    """ Find highest version of a table of a given type

    returns list
    inTL    = Python TableList
    tabType = Table type, e.g. "AIPS CC"
    """
    ################################################################
    # Checks
    if not PIsA(inTL):
        raise TypeError("inTL MUST be a Python Obit TableList")
    #
    return Obit.TableListGetHigh(inTL.me, tabType)
    # end PGetHigh


def PPutHi (inTL, err):
    """ Adds History to Table List

    inTL    = Python TableList
    err  = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inTL):
        raise TypeError("inTL MUST be a Python Obit TableList")
    #
    Obit.TableListPutHi(inTL.me, err.me)
    # end PPutHi


def PIsA (inTL):
    """ Tells if the input really is a Python Obit TableList

    returns true or false (1,0)
    inTL = Python TableList to test
    """
    ################################################################
     # Checks
    if not isinstance(inTL, TableList):
        return False
    #
    return Obit.TableListIsA(inTL.me)
    # end  PIsA


