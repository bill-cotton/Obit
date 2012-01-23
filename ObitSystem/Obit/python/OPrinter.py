""" Python Obit interface to printer

This class is for creating and using the interface to a printer
This formats text files for printing.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2012
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

# Python shadow class to ObitPrinter class
import Obit, InfoList

class OPrinterPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            Obit.OPrinterUnref(Obit.OPrinter_me_get(self.this))
            # In with the new
            Obit.OPrinter_me_set(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if self.__class__ != OPrinter:
            return
        if name == "me" : 
            return Obit.OPrinter_me_get(self.this)
        raise AttributeError,name
    def __repr__(self):
        if self.__class__ != OPrinter:
            return
        return "<C OPrinter instance> " + Obit.OPrinterGetName(self.me)

class OPrinter(OPrinterPtr):
    """
    Python Obit interface to printer server
    
    This class is for creating and using the interface to a printer
    to format text.

    OPrinter Members with python interfaces:

    * List - InfoList used to pass instructions to processing

    """
    def __init__(self, name="Printer", isInteractive=True, \
                 streamname="stdout", fileName="None"):
        """
        Create Printer object
        
        isInteractive = logical, True if interactive
        streamname    = name of output stream, "stdout", or "stderr"
        fileName      = file name if output to file desired., "None" => no file
        """
        self.this = Obit.new_OPrinter(name, isInteractive, streamname, fileName)
    def __del__(self):
        if Obit!=None:
            Obit.delete_OPrinter(self.this)
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.OPrinter_me_get(self.this)
        # Functions to return members
        if name=="List":
            out    = InfoList.InfoList()
            out.me = Obit.InfoListUnref(out.me)
            out.me = Obit.OPrinterGetList(self.me)
            return out
        raise AttributeError,str(name)
    def Open (self, err, LinesPerPage=50, Title1="Page title", Title2="    "):
        """
        Open a printer givint titles and page size

        * self   = Python Image object
        * access    = access READONLY (1), WRITEONLY (2), READWRITE(3)
        * err       = Python Obit Error/message stack
        * blc       = if given and a list of integers (min 2) giving
          bottom left corner (1-rel) of subimage
        * trc       = if given and a list of integers (min 2) giving
          top right corner (1-rel) of subimage
        """
        POpen(self, Title1, Title2, err, LinesPerPage=LinesPerPage)
        # end Open
    def Write (self, line, err):
        """
        Write a text line to a printer
        
        * self      = Python printer object
        * line      = line of text
        * err       = Python Obit Error/message stack
        """
        PWrite (self, err)
        # end Write
    def Close (self, err):
        """
        Close printer

        * self      = Python Printer object
        * err       = Python Obit Error/message stack
        """
        PClose (self, err)
        # end Close

    # end class OPrinter
    
def POpen (printer, Title1, Title2, err, LinesPerPage=50) :
    """
    Open printer and set page titles.

    * printer      = printer, the following is in the List member
        pageLimit long  Maximum number of pages of output def [50]
    * Title1       = first line of page title
    * Title2       = second line of page title
    * err          = Python Obit Error/message stack
    * LinesPerPage = how many lines per page.
    """
    ################################################################
    # Checks
    if not PIsA(printer):
        print "Actually ",printer.__class__
        raise TypeError,"printer MUST be a Python Obit Printer"
    Obit.OPrinterOpen (printer.me, LinesPerPage, Title1, Title2, err.me)
    # end POpen


def PWrite (printer, line, err) :
    """
    Write output

    * printer      = printer
    * line         = line of text
    * err          = Python Obit Error/message stack

    returns quit to indicate if interactive user bailing out
    """
    ################################################################
    # Checks
    if not PIsA(printer):
        print "Actually ",printer.__class__
        raise TypeError,"printer MUST be a Python Obit Printer"
    return Obit.OPrinterWrite (printer.me, line, err.me)
    # end PWrite

def PClose (printer, err) :
    """
    Close printer and flush buffers.

    * printer      = printer
    * err          = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(printer):
        print "Actually ",printer.__class__
        raise TypeError,"printer MUST be a Python Obit Printer"
    Obit.OPrinterClose (printer.me, err.me)
    # end PClose


def PSetTitle (printer, Title1, Title2, err) :
    """
    Update page titles

    * printer      = printer
    * Title1       = first line of page title
    * Title2       = second line of page title
    * err          = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(printer):
        print "Actually ",printer.__class__
        raise TypeError,"printer MUST be a Python Obit Printer"
    Obit.OPrinterSetTitle (printer.me, Title1, Title2, err.me)
    # end PSetTitle

def PNewPage (printer, err) :
    """
    Start new page of output

    * printer      = printer
    * err          = Python Obit Error/message stack

    returns quit to indicate if interactive user bailing out
    """
    ################################################################
    # Checks
    if not PIsA(printer):
        print "Actually ",printer.__class__
        raise TypeError,"printer MUST be a Python Obit Printer"
    return Obit.OPrinterNewPage (printer.me, err.me)
    # end PNewPage

def PIsA (printer):
    """
    Tells if the input is a Python ObitPrinter
    
    returns true or false (1,0)

    * printer = Python Obit Printer to test
    """
    ################################################################
      # Checks
    if printer.__class__ != OPrinter:
        return 0
    return Obit.OPrinterIsA(printer.me)
    # end PIsA


