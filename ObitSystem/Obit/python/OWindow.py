# $Id: OWindow.py,v 1.3 2005/10/17 12:05:51 bcotton Exp $
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
 
# Python shadow class to ObitDConCleanWindow class
import Obit, InfoList, OErr 

class OWindowPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.OWindow_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.OWindow_me_get(self.this)
        raise AttributeError,str(name)
    def __repr__(self):
        return "<C OWindow instance>"
class OWindow(OWindowPtr):
    """ Python Obit Image descriptor class

    This contains information about the Tables associated with an image or dataset

    Image Members with python interfaces:
    List      - (virtual) Python list of table names and numbers
    """
    def __init__(self) :
        self.this = Obit.new_OWindow()
    def __del__(self):
        if Obit!=None:
            Obit.delete_OWindow(self.this)

# Window type definitions (from ObitDConCleanWindow.h)
RectangleType = 0
RoundType = 1
        
def PCreate (name, mosaic, err):
    """ Create OWindow from an ImageMosaic

    name      = Name to be given to object
    mosaic    = Python ImageMosaic to attach
    err   = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not ImageMosaic.PIsA(mosaic):
        raise TypeError,"uvData MUST be a Python Obit UV"
    if err.isErr:
        return None  # existing error 
    #
    out = OWindow();
    out.me = Obit.OWindowCreate(name, mosaic.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error creating Window")
    return out;
    # end PCreate

def PCreate1 (name, naxis, err):
    """ Create single field OWindow

    name     = Name to be given to object
    naxis    = [x_dim,y_dim] size of image
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if err.isErr:
        return None  # existing error 
    #
    out = OWindow();
    out.me = Obit.OWindowCreate1(name, naxis, err.me)
    if err.isErr:
        printErrMsg(err, "Error creating Window")
    return out;
    # end PCreate1

def PGetList (inOW, field, err):
    """ Returns the contents of an OWindow field as a Python list

    returns list of lists, one per window
    window list for each window is ID, type, followed by parameters
    types: RectangleType, RoundType
    inOW  = Python OWindow to read
    field = Which field if mosaic
    err   = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOW):
        raise TypeError,"inOW MUST be a Python Obit OWindow"
    if err.isErr:
        return None  # existing error 
    #
    return Obit.OWindowGetList(inOW.me, field, err.me)
    # end PGetList


def PSetList (inOW, list, field, err):
    """ Copies list as from PGetList to an OWindow field

    Previous contents are deleted
    inOW  = Python OWindow to read
    field = Which field
    list  = list of window lists
            window list for each window is ID, type, followed by parameters
            types: RectangleType, RoundType
    err   = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOW):
        raise TypeError,"inOW MUST be a Python Obit OWindow"
    if err.isErr:
        return None  # existing error 
    #
    Obit.OWindowSetList(inOW.me, list, field, err.me)
    if err.isErr:
        printErrMsg(err, "Error setting Window InfoList")
    # end PSetList


def PAdd (inOW, field, type, window, err):
    """ Adds new window, 

    returns iD of new window
    inOW   = Python OWindow to read
    field  = Which field if mosaic
    type   = RectangleType or RoundType
    window = type dependent parameters
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOW):
        raise TypeError,"inOW MUST be a Python Obit OWindow"
    if err.isErr:
        return None  # existing error 
    #
    return Obit.OWindowAdd(inOW.me, field, type, window, err.me)
    # end PAdd

def PUpdate (inOW, field, iD, type, window, err):
    """ Updates window, creates new if iD doesn't exist

    inOW   = Python OWindow to read
    field  = Which field if mosaic
    iD     = iD of window to update
    type   = RectangleType or RoundType
    window = type dependent parameters
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOW):
        raise TypeError,"inOW MUST be a Python Obit OWindow"
    if err.isErr:
        return None  # existing error 
    #
    Obit.OWindowUpdate(inOW.me, field, iD, type, window, err.me)
    if err.isErr:
        printErrMsg(err, "Error updating Window")
    # end PUpdate

def PDel (inOW, field, iD, err):
    """ Deletes window, creates new if iD doesn't exist

    inOW   = Python OWindow to read
    field  = Which field if mosaic
    iD     = iD of window to delete, -1 => all
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inOW):
        raise TypeError,"inOW MUST be a Python Obit OWindow"
    if err.isErr:
        return None  # existing error 
    #
    Obit.OWindowDel(inOW.me, field, iD, err.me)
    if err.isErr:
        printErrMsg(err, "Error deleting Window")
    # end PDel

def PGetMaxID (inOW, field):
    """ Find highest version of a table of a given type

    returns list
    inOW    = Python OWindow
    field   = Which field
    tabType = Table type, e.g. "AIPS CC"
    """
    ################################################################
    # Checks
    if not PIsA(inOW):
        raise TypeError,"inOW MUST be a Python Obit OWindow"
    #
    return Obit.OWindowGetMaxID(inOW.me, field)
    # end PGetMaxID


def PIsA (inOW):
    """ Tells if the input really is a Python Obit OWindow

    returns true or false (1,0)
    inOW = Python OWindow to test
    """
    ################################################################
     # Checks
    if inOW.__class__ != OWindow:
        return 0
    #
    return Obit.OWindowIsA(inOW.me)
    # end  PIsA


