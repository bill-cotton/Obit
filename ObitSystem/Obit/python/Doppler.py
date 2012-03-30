""" Python Obit Doppler correction class/utility

    Image Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List 
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C)2012
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

# Doppler correction class
import Obit, OErr, InfoList, UV, UVDesc, types

# Python shadow class to ObitDoppler class
 
class DopplerPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.DopplerUnref(Obit.Doppler_me_get(self.this))
            # In with the new
            Obit.Doppler_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != Doppler:
            return
        if name == "me" : 
            return Obit.Doppler_me_get(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        raise AttributeError,str(name)
    def __repr__(self):
        if self.__class__ != Doppler:
            return
        return "<C Doppler instance> " + Obit.DopplerGetName(self.me)
#
class Doppler(DopplerPtr):
    """ Python Obit Image class

    This class handles Doppler corrections

    Doppler Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List
                (readonly)
    """
    def __init__(self, name) :
        self.this = Obit.new_Doppler(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_Doppler(self.this)


def newObit(name, err):
    """ Create and initialize an Doppler structure

    Create sky model object
    Returns the Python Doppler object
    name     = name desired for object (labeling purposes)
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = Doppler (name)
    return out      # seems OK
    # end newObit

def PCopy (inDoppler, outDoppler, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inDoppler, copies pointers
    inDoppler  = Python Doppler object to copy
    outDoppler = Output Python Doppler object, must be defined
    err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inDoppler):
        raise TypeError,"inDoppler MUST be a Python Obit Doppler"
    if not PIsA(outDoppler):
        raise TypeError,"outDoppler MUST be a Python Obit Doppler"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.DopplerCopy (inDoppler.me, outDoppler.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error copying Doppler")
    # end PCopy

def PGetList (inDoppler):
    """ Return the member InfoList

    returns InfoList
    inDoppler  = Python Doppler object
    """
    ################################################################
     # Checks
    if not PIsA(inDoppler):
        raise TypeError,"inDoppler MUST be a Python Obit Doppler"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.DopplerGetList(inDoppler.me)
    return out
    # end PGetList


def PCVel (inUV, outUV, RestFreq, err, scratch=False, \
           VLSR=0.0, refDate="0000-00-00", refChan=-9999.):
    """ Doppler correct a UV data set

    inUV     = Input UV data, any calibration/editing/selection parameters
               should be entered on the List member.
    outUV    = Output UV data, should be defined but not instantiated
               if not scratch
    RestFreq = Rest Frequency (GHz)
    err      = Obit error/message object
    scratch  = True if output a scratch file (destroyed on object deletion)
    VLSR     = desired center LSR velocity (km/s)
    refDate  = reference date for reference channel as 'yyyy-mm-dd'
               defaults to observation date.
    refChan  = reference channel (-9999=> compute)
    """
    ################################################################
    # Checks
    if not UV.PIsA(inUV):
        raise TypeError,"inUV MUST be a Python ObitUV"
    if not UV.PIsA(outUV):
        raise TypeError,"outUV MUST be a Python ObitUV"
    #
    if scratch:
        lscratch = 1
    else:
        lscratch = 0
        # set parameters
    info = inUV.List
    info.set("RestFreq", RestFreq*1.0e9, ttype="double")
    info.set("VelLSR",   VLSR*1000)
    info.set("JDref",    UVDesc.PDate2JD (refDate), ttype="double")
    info.set("refChan",  refChan)
    # Correct data
    Obit.DopplerCVel(inUV.me, lscratch, outUV.me, err.me)
    # end PCVel

def PDopplerFreqLSR(rest, vlsr, ra, dec, date, ut, \
                    x=-1.601225e06, y= -5.041980e06, z=3.554856e06):
    """ Return sky frequency of a given line

    Some rest frequencies
    HI 1.420406 GHz
    OH 1.612231, 1.665402, 1.667359. 1.72053 GHz
    CH3OH 6.668518, 36.169265 GHz
    H2O 22.235 GHz
    NH3 23.694, 23.723, 23.87 GHz
    SiO 42.820582, 43.122079, 214.0885, 215.5959 GHz
    CO 115.271 GHz
    returns frequency (GHz)
    rest    = Rest frequency (GHz)
    vlsr    = LSR velocity (km/s)
    ra      = Right ascension (deg)
    dec     = Declination (deg)
    date    = Date string as 'yyyy-mm-dd'
    ut      = Time (hr) wrt date
    x       = Earth centered X coordinate (m), default VLA W2
    y       = Earth centered Y coordinate (m)
    z       = Earth centered Z coordinate (m)
    """
    ################################################################
    return 1.0e-9*Obit.DopplerFreqLSR (rest*1.0e9, vlsr, ra, dec, date, ut, x, y, z)
    # end PDopplerFreqLSR

def PGetName (inDoppler):
    """ Tells Image object name (label)

    returns name as character string
    inDoppler  = Python Doppler object
    """
    ################################################################
     # Checks
    if not PIsA(inDoppler):
        raise TypeError,"inDoppler MUST be a Python Obit Doppler"
    #
    return Obit.DopplerGetName(inDoppler.me)
    # end PGetName

def PIsA (inDoppler):
    """ Tells if input really a Python Obit Doppler

    return true, false (1,0)
    inDoppler   = Python Doppler object
    """
    ################################################################
    # Checks
    if inDoppler.__class__ != Doppler:
        return 0
    return Obit.DopplerIsA(inDoppler.me)
    # end PIsA
