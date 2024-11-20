# Script to blank an image using an Image Window
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2024
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
 
# Python class for image slices
import Obit, OErr, Image, FArray, OWindow, ODisplay, ImageDesc
import os,pickle

class Blank():
    """ Python Obit Image blanking class
    
    Selection regions of an image
    1) Use set_mask to define the mask using the CLEAN window editing
       Pixels in unboxes will be blanked even if inside a box.
       Unless cleared, the window list is reused in subsequent calls.
    2) The CLEAN window list can be saved to/restored from a pickle file
       using save_window and fetch_window.
    3) Use make_mask to apply the window list to make a blanking mask
       in member Mask.
    4) Use blank to apply the mask to a specified plane in an image 
       compatible (in size) with the Mask.
       NB: this modifies the image, make a copy before running this.
    
    Members:
        Name   Name (label) for object
        Image  ObitImage to be blanked
        Window OWindow CLEAN window
        Mask   FArray to use for blanking image Image
    """
    def __init__(self, name, img, err):
        """
        Initialize blanking object for ObitImage img.
        """
        naxis = img.Desc.Dict['inaxes'][0:2]
        win = OWindow.PCreate1('win',naxis, err)
        mask = FArray.FArray("Mask",naxis=naxis)
        FArray.PFill(mask, FArray.fblank) # Init masks with blanks
        self.__dict__ = {"Name":name, "Image":img, "Window":win, "Mask":mask}; 
    
    def set_mask (self, disp, err):
        """
        Set blanking mask list using TV window edit
        disp should have already been defined in ObitTalk
        """
        print("Define blanking masks as windows, can use both boxes and unboxes")
        ODisplay.PImage(disp, self.Image, err, window=self.Window)
        OErr.printErr(err)
    # end set_mask
    def save_window (self, pfile, err):
        """
        Save window list to pickle file pfile
        """
        fd = open(pfile, "wb")
        # Can't save swig object - save list instead
        lwin = OWindow.PGetList(self.Window,1,err)
        pickle.dump(lwin, fd, protocol=2)
        fd.close()
    # end save_window
    def fetch_window (self, pfile, err):
        """
        Retrieve window list from pickle file pfile
        """
        # does it exist?
        if not os.path.lexists(pfile):
            print("Pickle jar",pfile,"does not exist")
            return
        fd = open(pfile, "rb")
        lwin = pickle.load(fd)
        fd.close()
        OWindow.PSetList(self.Window,lwin,1,err)
    # end fetch_window
    def make_mask (self,err):
        """
        Generate blanking mask
        """
        FArray.PFill(self.Mask, FArray.fblank) # Init masks with blanks
        lwin = OWindow.PGetList(self.Window,1,err)
        # Unblank Rectangle and Round windows
        for w in lwin:
            if w[1]==Rectangle:
                FArray.PRectFill(self.Mask, w[2:], 1.0)
            if w[1]==Round:
                FArray.PRoundFill(self.Mask, w[2:], 1.0)
        # Reblank UnRectangle and UnRound windows
        for w in lwin:
            if w[1]==UnRectangle:
                FArray.PRectFill(self.Mask, w[2:], FArray.fblank)
            if w[1]==UnRound:
                FArray.PRoundFill(self.Mask, w[2:], FArray.fblank)
    # end make_mask
    def blank (self,img,err,plane=1):
        """
        Blank plane plane in image img using Mask
        """
        img.GetPlane(None, [plane,1,1,1,1],err)
        FArray.PBlank(img.FArray, self.Mask, img.FArray)
        img.PutPlane(None, [plane,1,1,1,1],err)
    # end blank
# end class Blank
#Box types
Rectangle   = 0
Round       = 1
UnRectangle = 2
UnRound     = 3
