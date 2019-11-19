""" Python Obit interface to display server

This class is for creating and using the interface to an image display server
ONLY ONE MAY EXIST
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2019
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

# Python shadow class to ObitDisplay class
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, Image, ImageMosaic, OWindow

class ODisplay(Obit.ODisplay):
    """
    Python Obit interface to display server
    
    This class is for creating and using the interface to an image display server
    ONLY ONE MAY EXIST
    """
    def __init__(self, name, serverURL, err):
        super(ODisplay, self).__init__()
        Obit.CreateODisplay (self.this, name, serverURL, err.me)
        self.serverURL = serverURL
    def __del__(self, DeleteODisplay=_Obit.DeleteODisplay):
        if _Obit!=None:
            DeleteODisplay(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            if self.this!=None:
                Obit.ODisplayUnref(Obit.ODisplay_Get_me(self.this))
            # In with the new
            Obit.ODisplay_Set_me(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, ODisplay):
            return "Bogus dude"+str(self.__class__)
        if name == "me" : 
            return Obit.ODisplay_Get_me(self.this)
        if name == "url" : 
            return self.serverURL
        raise AttributeError(name)
    def __repr__(self):
        if not isinstance(self, ODisplay):
            return "Bogus Dude!"
        return "<C ODisplay instance> " + Obit.ODisplayGetName(self.me)

    def ping (self):
        """ See if Display server present
    
        Returns True if server present

        * self      = Display object
        """
        from six.moves.xmlrpc_client import ServerProxy
        url = self.serverURL
        if url == "ObitView":
            url = "http://localhost:8765/RPC2"
        server = ServerProxy(url)
        try:
            answer = server.ping(42)
        except:
            answer = False
            pass
        else:
            pass
        if answer:
            print("Display Server "+url+" present")
            return True
        else:
            print("Display Server "+url+" NOT present")
            return False
        # end ping

def PImage (disp, image, err, window=None) :
    """
    Display an image on the display server

    * disp     = display server
    * image    = Image to display
    * err      = Python Obit Error/message stack
    * window   = in not None, then edit this OWindow in the server

    returns True if worked
    """
    ################################################################
    # Checks
    if not PIsA(disp):
        print("Actually ",disp.__class__)
        raise TypeError("disp MUST be a Python Obit Display")
    if OWindow.PIsA(window):
        ret = Obit.ODisplayImageEdit(disp.me, image.me, window.me, err.me)
    else:
        ret = Obit.ODisplayImage(disp.me, image.me, err.me)
    return ret
    # end PImage


def PMosaic (disp, mosaic, field, err, window=None) :
    """
    Display an image mosaic on the display server

    * disp     = display server
    * mosaic   = Image Mosaic to display
    * field    = which field in mosaic (1-rel)
    * err      = Python Obit Error/message stack
    * window   = in not None, then edit this OWindow in the server

    returns True if worked
    """
    ################################################################
    # Checks
    if not PIsA(disp):
        print("Actually ",disp.__class__)
        raise TypeError("disp MUST be a Python Obit Display")
    if OWindow.PIsA(window):
        ret = Obit.ODisplayMosaic(disp.me, mosaic.me, field, window.me, err.me)
    else:
        ret = Obit.ODisplayMosaicEdit(disp.me, mosaic.me, field, err.me)
    return ret
    # end PImage

def PMarkPos (disp, pos, err) :
    """
    Display an image mosaic on the display server

    * disp     = display server
    * pos      = celestial ppositin as ("hh mm ss.s dd mm ss.s")
    * err      = Python Obit Error/message stack

    returns True if worked
    """
    ################################################################
    # Checks
    if not PIsA(disp):
        print("Actually ",disp.__class__)
        raise TypeError("disp MUST be a Python Obit Display")
    ret = Obit.ODisplayMarkPos(disp.me, pos, err.me)
    return ret
    # end PMarkPos


def PIsA (disp):
    """
    Tells if the input is a Python ObitDisplay
    
    returns True or False 
    * disp = Python Obit Display to test
    """
    ################################################################
    # Checks
    if not isinstance(disp, ODisplay):
        return False
    return Obit.ODisplayIsA(disp.me)!=0
    # end PIsA


