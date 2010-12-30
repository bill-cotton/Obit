""" Python Obit History class

This class contains processing history and allows access.
An ObitHistory is the front end to a persistent disk resident structure.
Both FITS  and AIPS cataloged data are supported.
To create a History, instantiate with the InfoList from the relevant object
to which the history is attached.
History(name, info, err)
History Members with python interfaces:
List  - used to pass instructions to processing (file info)
"""
# Python/Obit History class
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

# Python shadow class to ObitHistory class
import Obit, InfoList, OErr

class HistoryPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.HistoryUnref(Obit.History_me_get(self.this))
            # In with the new
            Obit.History_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != History:
            return
        if name == "me" : 
            return Obit.History_me_get(self.this)
        if name=="List":
            return PGetList(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != History:
            return
        return "<C History instance> " + Obit.HistoryGetName(self.me)
class History(HistoryPtr):
    """ Python Obit History class
    
    This class contains processing history and allows access.
    An ObitHistory is the front end to a persistent disk resident structure.
    Both FITS  and AIPS cataloged data are supported.
    To create a History, instantiate with the InfoList from the relevant object
    to which the history is attached.
    History(name, info, err)
    History Members with python interfaces:
    List  - used to pass instructions to processing (file info)
    """
    def __init__(self,name,info,err) :
        self.this = Obit.new_History(name, info.me, err.me)
    def __del__(self):
        if Obit!=None:
            Obit.delete_History(self.this)

    def Zap (self, err):
        """ Destroy the persistent form of a History
        
        self = input Python Obit History
        err  = Python Obit Error/message stack
        """
        PZap (self, err)
        # end Zap 

    def Open (self, access, err):
        """ Open History
        
        return 0 on success, else failure
        self    = input Python History
        access   = access code READONLY(1), WRITEONLY(2), READWRITE(3)
        err      = Python Obit Error/message stack
        """
        return POpen (self, access, err)
        # end POpen

    def Close (self, err):
        """ Close History
        
        return 0 on success, else failure
        self  = input Python History
        err   = Python Obit Error/message stack
        """
        return PClose (self, err)
    # end PClose

    def ReadRec (self, recno, err):
        """ Read History record
        
        returns string
        self  = input Python History
        recno = desired record
        err   = Python Obit Error/message stack
        """
        return PReadRec (self, recno, err)
    # end PReadRec

    def WriteRec (self, recno, hiCard, err):
        """ Write History record
        
        return 0 on success, else failure
        self   = input Python History
        recno  = desired record, -1 => end of table
        hiCard = input history record
        err    = Python Obit Error/message stack
        """
        return PWriteRec (self, recno, hiCard, err)
    # end PWriteRec

    def Stalin (self, startr, endr, err):
        """ Edit history
        
        return 0 on success, else failure
        self   = input Python History
        startr = first (1-rel) history record to delete
        endr   = highest (1-rel) history record to delete, 0->to end
        err    = Python Obit Error/message stack
        """
        return PEdit (self, startr, endr, err)
    # end Stalin

    def TimeStamp (self, label, err):
        """ Write timestamp and label to History
        
        return 0 on success, else failure
        self  = input Python History
        label = character string for label
        err   = Python Obit Error/message stack
        """
        return PTimeStamp (self, label, err)
    # end TimeStamp


# Symbolic names for access codes
READONLY  = 1
WRITEONLY = 2
READWRITE = 3

def PZap (inHis, err):
    """ Destroy the persistent form of a History

    inHis    = input Python Obit History
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.HistoryZap (inHis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error deleting History")
    # end PZap 

def PCopy (inHis, outHis, err):
    """ Copy a History including persistent forms

    inHis    = input Python Obit History
    outHis   = extant output Python Obit History
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not PIsA(outHis):
        raise TypeError,"outHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.HistoryCopy(inHis.me, outHis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying History")
    # end PCopy

def PCopyHeader (inHis, outHis, err):
    """ Copy a History from a (FITS) header

    inHis    = input Python Obit History (FITS)
    outHis   = extant output Python Obit History
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not PIsA(outHis):
        raise TypeError,"outHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.HistoryCopyHeader(inHis.me, outHis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying History from header")
    # end PCopyHeader

def PCopy2Header (inHis, outHis, err):
    """ Copy a History to a (FITS) header

    inHis    = input Python Obit History 
    outHis   = extant output Python Obit (FITS) History in header
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        print "inHis really",inHis.__class__ 
        raise TypeError,"inHis MUST be a History"
    if not PIsA(outHis):
        print "outHis really",outHis.__class__ 
        raise TypeError,"outHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.HistoryCopy2Header(inHis.me, outHis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying History to header")
    # end PCopy2Header


def PHeader2Header (inHis, outHis, err):
    """ Copy a History from a FITS header to a (FITS) header

    inHis    = input Python Obit History 
    outHis   = extant output Python Obit (FITS) History in header
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not PIsA(outHis):
        raise TypeError,"outHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.HistoryHeader2Header(inHis.me, outHis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying History header to header")
    # end PHeader2Header


def POpen (inHis, access, err):
    """ Open History

    return 0 on success, else failure
    inHis    = input Python History
    access   = access code READONLY(1), WRITEONLY(2), READWRITE(3)
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.HistoryOpen(inHis.me, access, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening History")
    return ret
    # end POpen

def PClose (inHis, err):
    """ Close History

    return 0 on success, else failure
    inHis    = input Python History
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.HistoryClose(inHis.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error closing History")
    return ret
    # end PClose

def PReadRec (inHis, recno, err):
    """ Read History record

    returns string
    inHis    = input Python History
    recno    = desired record
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.HistoryReadRec(inHis.me, recno, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error reading history record")
    return ret
    # end PReadRec

def PWriteRec (inHis, recno, hiCard, err):
    """ Write History record

    return 0 on success, else failure
    inHis    = input Python History
    recno    = desired record
    hiCard   = input history record
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.HistoryWriteRec(inHis.me, recno, hiCard, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing History record")
    return ret
    # end PWriteRec

def PEdit (inHis, startr, endr, err):
    """ Edit History

    Deletes a range of history records.
    return 0 on success, else failure
    inHis    = input Python History
    startr   = first (1-rel) history record to delete
    endr     = highest (1-rel) history record to delete, 0=>to end
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = POpen (inHis, READWRITE, err)
    ret = Obit.HistoryEdit(inHis.me, startr, endr, err.me)
    OErr.printErr(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error Editing History")
    ret = PClose (inHis, err)
    return ret
    # end PEdit

def PTimeStamp (inHis, label, err):
    """ Write timestamp and label to History

    return 0 on success, else failure
    inHis    = input Python History
    label    = character string for label
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    ret = Obit.HistoryTimeStamp(inHis.me, label, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error writing History time stamp")
    # end PTimeStamp

def PGetList (inHis):
    """ Return the InfoList from a History
    
    returns InfoList
    inHis    = input Python History
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    #
    out    = InfoList.InfoList()
    out.me = Obit.HistoryGetList(inHis.me)
    return out
    # end PGetList

def PIsA (inHis):
    """ Tells if object thinks it's a Python Obit History

    return true, false (1,0)
    inHis    = input Python History
    """
    ################################################################
    # Checks
    if inHis.__class__ != History:
        return 0
    #
    return Obit.HistoryIsA(inHis.me)
    # end PIsA

def PGetName (inHis):
    """ Returns object name (label)

    return name string
    inHis    = input Python History
    """
    ################################################################
    # Checks
    if not PIsA(inHis):
        raise TypeError,"inHis MUST be a History"
    #
    return Obit.HistoryGetName(inHis.me)
    # end PGetName
