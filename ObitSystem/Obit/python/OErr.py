# $Id: OErr.py,v 1.11 2005/12/08 02:01:20 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2004
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

# Python shadow class to ObitErr class
import Obit

class OErrPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.OErr_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.OErr_me_get(self.this)
        if name == "isErr" : 
            return PIsErr(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C OErr instance>"
    def __str__(self):
        messages = ''
        msg = Obit.OErrMsg(self.me)
        while msg:
            messages += '%s\n' % msg
            msg = Obit.OErrMsg(self.me)
            continue
        Obit.ObitErrClear (self.me);  # Clear stack
        return messages
class OErr(OErrPtr):
    """ Python ObitErr message and error stack
    
    This is an error stack class for obtaining tracebacks for error conditions.
    This is also the mechanism for passing informative messages.
    No messages, error or informative, are displayed until the contents
    of the stack are explicitly printed
    """
    def __init__(self) :
        self.this = Obit.new_OErr()
    def __del__(self):
        if Obit!=None:
            Obit.delete_OErr(self.this)

    def Clear(self):
        """ Clear Obit error stack """
        PClear(self)
        # end Clear

# Error levels
NoErr     = 0
Info      = 1
Warn      = 2
Traceback = 3
MildError = 4
Error     = 5
StrongError = 6
Fatal     = 7

def PIsErr(err):
    """ Tells if an error condition exists

    Returns True if error condition exists, else False
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    return Obit.isError(err.me) != 0
    # end PIsErr

def PClear(err):
    """ Clear Obit error stack

    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    Obit.ObitErrClear(err.me)
    #end PClear

def PSet(err):
    """ Set Obit error flag

    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    Obit.SetError(err.me)
    #end PSet

def PLog(err, eCode, message):
    """ Add message To Obit Error/message stack

    err      = Python Obit Error/message stack
    eCode    = error code defined above:
               NoErr, Info, Warn, Traceback,
               MildError, Error, StrongError, Fatal
    """
    ################################################################
    # Checks
    if not OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    Obit.LogError(err.me, eCode, message)
    #end PLog

def printErr(err):
    """ Prints Obit error stack
    
    err      = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    Obit.ObitErrLog(err.me)
    # end PrintErr
     
def printErrMsg(err, message="Error"):
    """ Prints Obit error stack and throws runtime exception on error

    err     = Python Obit Error/message stack
    message = message string for exception
    """
    ################################################################
    # Checks
    if not OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    ierr = Obit.isError(err.me)
    Obit.ObitErrLog(err.me)
    if ierr:
        print message
        raise OErr
    # end printErrMsg
     
def OErrIsA (err):
    """ Tells if object thinks it's a Python ObitErr

    return true, false (1,0)
    err    = input Python ObitErr stack
    """
    ################################################################
    # Checks
    if err.__class__ != OErr:
        return 0
    #
    return Obit.ObitErrIsA(err.me)
    # end OErrIsA

def Bomb ():
    """ Throws an exception to stop the debugger
    """
    ################################################################
    #
    Obit.Bomb()
    # end Bomb
