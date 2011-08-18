"""
Python utility package for parameter file parser
"""

# Python utility package for parameter file parser
import Obit, InfoList, OErr
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007
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
#  Correspondence about concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

def PParse(infile, list, err):
    """ 
    Parse text file

    Parse parameter format text file and copy to list
    The input file is basically free format with a keyword=value syntax.  
    Comments follow a "#" symbol.  String keywords should have no leading or
    trailing blanks before the end of the line or a comment delimiter.
    If a numeric value is followed by a comment there should be at least one
    blank before the comment delimiter.  Allowed types are:
    Str (String), Flt (float), Dbl (double), Int (int), Boo (boolean)

    Examples::

        # Comment
        $Key = SomeStr   Str (9)   # Some string
        Something
        $Key =  FArr     Flt (2)   # Array of 2 floats
        0.0 0.0
        $Key = CLEANBox  Int(4,1)  # Clean box with 4 ints
        0 0 0 0
    
    returns  0 on success, else 1

    * infile   = Name of the input text file to parse
    * list     = ObitInfoList to accept values.
    * err      = ObitErr for reporting errors.
    """
    ################################################################
    # Checks
    if not InfoList.PIsA(list):
        raise TypeError,"list MUST be a Python Obit InfoList"
    ret = Obit.Parse(infile, list.me, err.me)
    if ret or err.isErr:
        OErr.printErrMsg(err, "Error Parsing "+infile)
    return ret
    # end PParse

def PDump(infile, list, err):
    """
    Dump list to text file
    
    Dump parameter format text file and copy to list
    returns  0 on success, else 1

    * outfile  = Name of the output text file to write
    * list     = ObitInfoList to dump
    * err      = ObitErr for reporting errors.
    """
    ################################################################
    # Checks
    if not InfoList.PIsA(list):
        raise TypeError,"list MUST be a Python Obit InfoList"
    ret = Obit.Dump(infile, list.me, err.me)
    if ret or err.isErr:
        OErr.printErrMsg(err, "Error Dumping to "+infile)
    return ret
    # end PDump

