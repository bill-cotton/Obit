""" Python Obit Source class

    This contains information about an astronomical source
    Source Members with python interfaces:
      Name    Source Name, first SU table entry matching Name returned.
      SourID  Source ID
      Qual    Source qualifier
      RAMean  Mean RA in deg
      DecMean Mean Dec in deg
      RAApp   Apparent RA in deg
      DecApp  Apparent Dec in deg
      Equinox Equinox (years) of mean positions

"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2013
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
 
# Python shadow class to ObitSource class
import Obit, OErr, string, math

class SourcePtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.Source_me_set(self.this,value)
            return
        if name=="RAMean":
            Obit.SourceSetRAMean(self.me,value);
            return
        if name=="DecMean":
            Obit.SourceSetDecMean(self.me,value);
            return
        if name=="RAApp":
            Obit.SourceSetRAApp(self.me,value);
            return
        if name=="DecApp":
            Obit.SourceSetDecApp(self.me,value);
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.Source_me_get(self.this)
        # Functions to return members
        if name=="RAMean":
            return Obit.SourceGetRAMean(self.me);
        if name=="DecMean":
            return Obit.SourceGetDecMean(self.me);
        if name=="RAApp":
            return Obit.SourceGetRAApp(self.me);
        if name=="DecApp":
            return Obit.SourceGetDecApp(self.me);
        if name=="Name":
            return Obit.SourceGetName(self.me)
        if name=="SourID":
            return Obit.SourceGetSID(self.me)
        if name=="Qual":
            return Obit.SourceGetQual(self.me)
        if name=="Equinox":
            return Obit.SourceGetEquinox(self.me)
        raise AttributeError,str(name)
    def __repr__(self):
        return "<C Source instance>"
class Source(SourcePtr):
    """ Python Obit Source class

    """
    def __init__(self, name):
        self.this = Obit.new_Source(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_Source(self.this)


def PIsA (inSou):
    """ Tells if the input really is a Python Obit Source

    returns true or false (1,0)
    inUD = Python Source to test
    """
    ################################################################
     # Checks
    if inSou.__class__ != Source:
        return 0
    #
    return Obit.SourceIsA(inSou.me)
    # end  PIsA

def PCreateByNumber (name, inUV, SouID, OErr):
    """ Create a Source from a UV for a specified Source ID

    If inUV has a SoUrce table the Source is extracted from it, otherwise
    from the information in the descriptor
     * name   = Name for object
     * inUV   = UV data for source info
     * SouID  = Source identifier in inUV
     * err    = Python Obit Error/message stack
      returns  new Source object
    """
    ################################################################
    out    = Source(name)
    out.me = Obit.SourceCreateByNumber(name, inUV.me, SouID, OErr.me)
    return out
    # end PCreateByNumber

def PCreateByName (name, inUV, Souce, Qual, OErr):
    """ Create a Source from a UV for a specified Source name/qual

    If inUV has a SoUrce table the Source is extracted from it, otherwise
    from the information in the descriptor
     * name   = Name for object
     * inUV   = UV data for source info
     * Source = Source name in inUV
     * Qual   =  Source qualifier  in inUV
     * err    = Python Obit Error/message stack
      returns  new Source object
    """
    ################################################################
    out    = Source(name)
    out.me = Obit.SourceCreateByName(name, inUV.me, Source, Qual, OErr.me)
    return out
    # end PCreateByName




