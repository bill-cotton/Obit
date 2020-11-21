""" Python Obit UVVis class

    This contains data from a single UV visibility measurement

    UVVis Members with python interfaces:
    u    = u coordinate (lambda)
    v    = v coordinate (lambda)
    w    = w coordinate (lambda)
    time = Visibility time in days since 0 h on reference day
    ant1 = antenna 1 of baseline
    ant2 = antenna 2 of baseline
    suba = subarray number
    suid = Source ID number
    vis  = visibilities as list of (complex, float) (vis, wt)
"""
# $Id: UVVis.py 452 2013-06-03 14:40:12Z bill.cotton $
#-----------------------------------------------------------------------
#  Copyright (C) 2007-2020
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
 
# Python class for Obit visibility measurements
from __future__ import absolute_import
import Obit, UV, OErr, string, math

class UVVis():
    def __init__(self):
        self.EOF = False
        self.u = 0.0
        self.v = 0.0
        self.w = 0.0
        self.time = 0.0
        self.ant1 = 0
        self.ant2 = 0
        self.suid = 0
        self.vis=[(complex(0.0,00),1.0)]
    def __setattr__(self,name,value):
        self.__dict__[name] = value
    def __getattr__(self,name):
        raise AttributeError(str(name))
    def __repr__(self):
        return "<C UVVis instance>"

def PGet (inUV, err):
    """ Read next buffer, Generate a visibility from a ObitUV buffer

    Returns UVVis object, EOF member set to True if all data read
    inUV      = Python UV object, file should be opened read enabled
                and with 1 vis per read.
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError("inUV MUST be a Python Obit UV")
    #
    out = UVVis()
    vis = Obit.UVVisGet(inUV.me, err.me)
    if 'EOF' in vis:
        out.EOF = vis['EOF']
    out.u = vis['u']
    out.v = vis['v']
    out.w = vis['w']
    out.time = vis['time']
    out.ant1 = vis['ant1']
    out.ant2 = vis['ant2']
    if 'suid' in vis:
        out.suid = vis['suid']
    out.vis  = vis['vis']
    return out
    # end PGet


def PSet (inVis, inUV, err):
    """ Copy a visibility to a ObitUV buffer

    inUVis    = Python UVVis to write
    inUV      = Python UV object, file should be opened write enabled
                and with 1 vis per write.
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not inUV.UVIsA():
        raise TypeError("inUV MUST be a Python Obit UV")
    #
    Obit.UVVisSet(inVis.__dict__, inUV.me, err.me)
    # end PSet


