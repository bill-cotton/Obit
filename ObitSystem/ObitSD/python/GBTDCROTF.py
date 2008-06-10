# $Id: GBTDCROTF.py,v 1.5 2005/12/19 00:37:09 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2004,2005
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
 
# Python shadow class to ObitGBTDCROTF class
import Obit, OTF, OErr, types

class GBTDCROTFPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
     #def __del__(self):
     #    if self.thisown == 1 :
     #       # If Obit has been unloaded don't bother
     #       if Obit.__class__ == Obit:
     #           Obit.delete_GBTDCROTF(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.GBTDCROTF_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.GBTDCROTF_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        return "<C GBTDCROTF instance>"
class GBTDCROTF(GBTDCROTFPtr):
    """ Python Obit GBTDCROTF  GBT DCR to OTF conversion class

    This class converts from GBT DCR format data to OTF format a scan
    at a time appending the data to the member OTF.
    """
    def __init__(self,name) :
        self.this = Obit.new_GBTDCROTF(name)
        self.thisown = 1
    def __del__(self):
        if Obit!=None:
            Obit.delete_GBTDCROTF(self.this)

def newPGBTDCROTF(name, outOTF, err):
    """ Create and initialize an OGBTDCROTF structure

    Create and save output OTF member.
    Returns the Python PBTDCROTF object
    name     = name desired for object (labeling purposes)
    outOTF   = Extant OTF (possible uninitialized) to write
    err      = Python Obit Error/message stack
    """
    ################################################################
    out = GBTDCROTF("None")
    out.me = Obit.GBTDCROTFValue(name, outOTF.me, err.me)
    return out 
    # end newPGBTDCROTF
 
def PConvert (inGDO, inDisk, scanName, err):
    """  Convert a scan of DCR data and appent to OTF

    inGDO   = Python Obit input GBTDCROTF with attached OTF
    inDisk  = FITS disk number of base of input
    scanName= Name of scan (e.g. (2004_04_16_05:31:47")
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inGDO):
        raise TypeError,"inGDO MUST be a Python Obit GBTDCROTF"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    Obit.GBTDCROTFConvert(inGDO.me, inDisk, scanName, err.me)
    # end PConvert

 
def PIsA (inGDO):
    """ Tells if the input really is a Python Obit GBTDCROTF

    returns true or false (1,0)
    inGDO = Python GBTDCROTF to test
    """
    ################################################################
     # Checks
    if inGDO.__class__ != GBTDCROTF:
        return 0
    #
    return Obit.GBTDCROTFIsA(inGDO.me)
    # end  PIsA

