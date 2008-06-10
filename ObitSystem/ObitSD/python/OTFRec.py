""" Python Obit OTFRec class

    This contains data from a single OTF measurement

    OTFRec Members with python interfaces:
    time     = Time of sample in days
    timei    = Integration time of sample in days
    scan     = Scan index
    target   = Target index
    ra       = RA (deg) of sample
    dec      = Declination (deg) of sample
    rot      = Rotation angle (parallactic angle) of sample
    cal      = <> 0 => cal on
    data     = data as list of tuples (data, wt)
"""
# $Id: OTFRec.py,v 1.1 2008/01/14 12:15:40 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2008
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
 
# Python class for Obit OTF measurements
import Obit, OTF, OErr, string, math

class OTFRecPtr :
    def __init__(self):
        self.EOF    = False
        self.time   = 0.0
        self.timei  = 0.0
        self.scan   = 0
        self.target = 0
        self.ra     = 0.0
        self.dec    = 0.0
        self.rot    = 0.0
        self.cal    = 0.0
        self.data=[(0.0,1.0)]
    def __setattr__(self,name,value):
        self.__dict__[name] = value
    def __getattr__(self,name):
        raise AttributeError,str(name)
    def __repr__(self):
        return "<C OTFRec instance>"
class OTFRec(OTFRecPtr):
    """ Python Obit OTF record class

    """

def PGet (inOTF, err):
    """ Read next buffer, Generate a record from an ObitOTF buffer

    Returns OTFRec object, EOF member set to True if all data read
    inOTF     = Python OTF object, file should be opened read enabled
                and with 1 record per read.
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inOTF.OTFIsA():
        raise TypeError,"inOTF MUST be a Python ObitOTF"
    #
    out = OTFRec()
    out.__dict__ = Obit.OTFRecGet(inOTF.me, err.me)
    return out
    # end PGet


def PSet (inRec, inOTF, err):
    """ Copy an OTF record to a ObitOTF buffer

    inRec     = Python OTFRec to write
    inOTF     = Python OTF object, file should be opened write enabled
                and with 1 vis per write.
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inOTF.OTFIsA():
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    #
    Obit.OTFRecSet(inRec.__dict__, inOTF.me, err.me)
    # end PSet

def PGeomProj (inRec, inOTF):
    """ Get Ra and Declination for each detector for a OTFRec

    Return a dict with:
    "xpos" => array of detector right ascensions
    "ypos" => array of detector declinations
    inRec     = Python OTFRec 
    inOTF     = Python OTF object, file should be opened.
    """
    ################################################################
    # Checks
    if not inOTF.OTFIsA():
        raise TypeError,"inOTF MUST be a Python Obit OTF"
    #
    return Obit.OTFRecGeomProj(inRec.__dict__, inOTF.me)
    # end PSet


