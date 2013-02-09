""" Python Obit AntennaList class

    This contains information about an antenna list/ time information for 
    a UV data.
    Functions are provided to calculate the azimith, elevation and 
    parallactic angles of a given source for an antenna at a given time.
    Source Members with python interfaces:
        RefJD    Julian data of reference day
        ArrName  Array name

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
 
# Python shadow class to ObitAntennaList class
import Obit, OErr, string, math

class AntennaListPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.AntennaList_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.AntennaList_me_get(self.this)
        # Functions to return members
        if name=="JDRef":
            return Obit.AntennaListGetRefJD(self.me);
        if name=="ArrName":
            return Obit.AntennaListGetArrName(self.me);
        raise AttributeError,str(name)
    def __repr__(self):
        return "<C AntennaList instance>"
class AntennaList(AntennaListPtr):
    """ Python Obit AntennaList class

    """
    def __init__(self, name, inUV, subA, err):
        self.this = Obit.new_AntennaList(name, inUV.me, subA, err.me)
    def __del__(self):
        if Obit!=None:
            Obit.delete_AntennaList(self.this)
    def Elev(self, ant, time, Source):
        """ Returns elevation (rad) of a source at a given antenna and time
        
        * self    AntennaList
        * ant     Antenna number in self
        * time    Time (days) wrt reference day
        * Source  Source in question
        """
        return Obit.AntennaListGetElev (self.me, ant, time, Source.me)
    # end Elev

    def Az(self, ant, time, Source):
        """ Returns azimuth (rad) of a source at a given antenna and time
        
        * self    AntennaList
        * ant     Antenna number in self
        * time    Time (days) wrt reference day
        * Source  Source in question
        """
        return Obit.AntennaListGetAz (self.me, ant, time, Source.me)
    # end Az

    def ParAng(self, ant, time, Source):
        """ Returns parallactic angle (rad) of a source at a given antenna and time
        
        * self    AntennaList
        * ant     Antenna number in self
        * time    Time (days) wrt reference day
        * Source  Source in question
        """
        return Obit.AntennaListGetParAng (self.me, ant, time, Source.me)
    # end ParAng

def PIsA (inSou):
    """ Tells if the input really is a Python Obit AntennaList

    returns true or false (1,0)
    inUD = Python AntennaList to test
    """
    ################################################################
     # Checks
    if inSou.__class__ != AntennaList:
        return 0
    #
    return Obit.AntennaListIsA(inSou.me)
    # end  PIsA



