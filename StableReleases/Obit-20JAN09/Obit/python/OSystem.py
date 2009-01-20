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

# Python shadow class to ObitSystem class
import Obit

# I don't exist until created
#ObitSys = None
#

class OSystemPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" :
            return Obit.OSystem_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C OSystem instance>"
class OSystem(OSystemPtr):
    """ Obit System Object class

    An Obit system is needed to access persistent forms of data objects.
    In particular, the ObitSystem keeps track of and scratch files creatd and deletes
    any left when the system is shutdown (OSystem object deleted)
    Currently Obit can access FITS files or AIPS (some classes) and the locations
    of data directories are specified at OSystem startup.  All external files will
    be in one of the directories specified and indicated by the "disk" variable
    inside Python Obit software.
    The system is started by creating an OSystem,
    the constructor has the following arguments:
    pgmName        = A name for the program, used in creating scratch files
                     Available from python object as member pgmName
    pgmNumber      = A version number of the program (POPS number for AIPS)
                     Available from python object as member pgmNumber
    AIPSuser       = AIPS user number if accessing AIPS data structures
    numberAIPSdisk = Number of AIPS disk directories
                     if -1 use values of $DA01, $DA02...
    AIPSdir        = List of numberAIPSdisk directories for AIPS data areas
                     For default case, give a list with "Def"
    numberFITSdisk = Number of FITS disk directories
                     if -1 use values of $DA01, $DA02...
    FITSdir        = List of numberFITSdisk directories for FITS data areas
                     For default case, give a list with "Def"
    F_TRUE         = Value of boolean TRUE (needed for Fortran interface)
    F_FALSE        = Value of boolean FALSE (needed for Fortran interface)
    err            = Obit Error stack on which to report problems.

    An example defining one FITS directory and no AIPS directories:
    # Init Obit
    err=OErr.OErr()
    ObitSys=OSystem.OSystem ("OTFSub", 1, 103, 1, ["None"], 1, ["../FITSdata/"], 1, 0, err)

    """
    def __init__(self,pgmName, pgmNumber, AIPSuser,
                 numberAIPSdisk, AIPSdir, numberFITSdisk, FITSdir,
                 F_TRUE, F_FALSE, err) :
        self.this = Obit.new_OSystem(pgmName, pgmNumber, AIPSuser,
                                     numberAIPSdisk, AIPSdir,
                                     numberFITSdisk, FITSdir,
                                     F_TRUE, F_FALSE, err.me)
        if err.isErr:
            printErrMsg(err, "Error in Obit initialization")
        # save info visible to python
        self.pgmName   = pgmName
        self.pgmNumber = pgmNumber
        # Rember who I am - can be only one
        OSystem.ObitSys = self
    def __del__(self):
        # Better not shutdown?
        if Obit!=None and self!=None:
            Obit.delete_OSystem(self.this)
    # End class definitions

        
def Shutdown (inObj=None):
    """ Shuts down Obit in all known IO related modules

    inObj   = Python Obit System object, defaults to OSystem.ObitSys
    """
    ################################################################
    #  Use default if none given
    if inObj.__class__ != OSystem:
        obj = OSystem.ObitSys
    else:
        obj = inObj
    Obit.Shutdown(obj.me)
    del (obj)
    # end Shutdown

def PIsInit ():
    """ Tells if Obit initialized

    returns True if Obit initialized else False
    """
    ################################################################
    retval = Obit.SystemIsInit();
    return retval!=0;
    # end PIsInit

def PGetAIPSuser ():
    """ Tells AIPS user number

    returns AIPS user number
    """
    ################################################################
    return Obit.SystemGetAIPSuser()
    # end PGetAIPSuser

def PSetAIPSuser (user):
    """ Sets  AIPS user number

    user = AIPS user number
    """
    ################################################################
    Obit.SystemSetAIPSuser(user)
    # end PSetAIPSuser

def PGetPgmName ():
    """ Tells Program name 

    returns name as character string
    """
    ################################################################
    return Obit.SystemGetPgmName()
    # end PGetPgmName

def PToday ():
    """ Get today's date in string suitable for descriptors

    returns date as character string
    """
    ################################################################
    return Obit.SystemToday()
    # end PToday

def PSetPgmName (pgmName):
    """ Sets Program name

    pgmName = new program name
    """
    ################################################################
    Obit.SystemSetPgmName(pgmName)
    # end PSetPgmName

def PGetPgmNumber ():
    """ Tells Program number 

    returns number
    """
    ################################################################
    return Obit.SystemGetPgmNumber()
    # end PGetPgmNumber

def PSetPgmNumber (pgmNumber):
    """ Sets Program number

    pgmNumber = new program number
    """
    ################################################################
    Obit.SystemSetPgmNumber(pgmNumber)
    # end PSetPgmNumber

def PMemPrint ():
    """ Prints contents of Obit Memory allocation on stdout
    """
    ################################################################
    #
    Obit.MemPrint()
    # end PMemPrint

def PAllowThreads (nThreads):
    """ Sets maximum number of threads in an Obit thread pool
    
    nThreads maximum number of threads in an Obit thread pool
    """
    ################################################################
    Obit.SystemAllowThreads(nThreads)
    # end PAllowThreads

def PGetNoThreads ():
    """ Tells Number of Threads enabled in Obit

    returns number
    """
    ################################################################
    return Obit.SystemGetNoThreads()
    # end PGetNoThreads

