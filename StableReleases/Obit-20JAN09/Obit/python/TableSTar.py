""" AIPS STar table

Due to the funky nature of the AIPS STar table it cannot be made in the usual
Obit fashion.  This class allows doing this from python.

Symbol type codes
   1: Plus sign (default)   12: Five pointed star
   2: Cross (X)             13: Star of David
   3: Circle                14: Seven-pointed star
   4: Box                   15: Eight-pointed star
   5: Triangle              16: Nine-pointed star
   6: Diamond               17: Ten-pointed star
   7: Pentagon              18: 11-pointed star
   8: Hexagon               19: 12-pointed star
   9: Septagon              20: 13-pointed star
   10: Octagon              21: 14-pointed star
   11: Nine-gon             22: Plus with gap
                            23: Vertical line
                            24: Cross (X) with gap
"""
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
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
import Obit, Table, TableDesc, OErr, Image, ImageDesc

class TableSTar(Table.Table):
    pass
    # end class TableSTar


# Data type codes
OBIT_double = 10
OBIT_float = 9
OBIT_string = 13
OBIT_int = 2

# Non class functions
def PCreate(im, err, ver=0):
    """
    New AIPS STars table
    
    Create a ST table on input image im
    im  = Obit Image on which to attach ST Table
    err = Python Obit Error/message stack
    ver = version, 0=> new
    """
    ################################################################
    # Check
    if not im.ImageIsA():
        raise TypeError,'im MUST be a Python Obit Image'
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return None
    
    # Get image descriptor
    id = im.Desc.Dict
    # Set descriptor dict
    dd = {"FieldName":[id["ctype"][0], id["ctype"][1], "MAJOR AX", "MINOR AX", \
                       'POSANG', 'STARTYPE',  'LABEL', \
                       "_status"], \
          "FieldUnit":["DEGREES", "DEGREES", "DEGREES", "DEGREES", \
                       "DEGREES", "INDEX   ", "        ", "        "], \
          "repeat":[1,1,1,1,1,1,24,1], \
          "dim0":[1,1,1,1,1,1,24,1], \
          "dim1":[1,1,1,1,1,1,1,1], \
          "dim2":[1,1,1,1,1,1,1,1], \
          "type":[OBIT_double,OBIT_double,OBIT_float,OBIT_float,OBIT_float,OBIT_float,\
                  OBIT_string,OBIT_int], \
          "sortOrder1":0, "sortOrder2":0, "Table name":"AIPS ST", "version":1 \
          }
    # Table descriptor
    tabDesc = TableDesc.PDef(dd)
    # Table
    st = im.NewTable(Table.WRITEONLY,"AIPS ST",ver,err)
    Obit.TableSetDesc(st.me, tabDesc.me)
    # Instantiate
    Table.PFullInstantiate(st, Table.WRITEONLY, err)
    return st
    # end PCreate
    
def newRow (im):
    """ Create new row structure for writing ST Table
    
    im   = Obit Image on which to attach ST Table
    returns row:
    Position columns have labelws of first two axes of image
    (e.g. 'RA---SIN', 'DEC--SIN')
    'MAJOR AX' major axis of symbol
    'MINOR AX  Minor axis of symbol (deg)
    'POSANG'   Position angle in deg
    'STARTYPE' symbol code
    1: Plus sign (default)   12: Five pointed star
    2: Cross (X)             13: Star of David
    3: Circle                14: Seven-pointed star
    4: Box                   15: Eight-pointed star
    5: Triangle              16: Nine-pointed star
    6: Diamond               17: Ten-pointed star
    7: Pentagon              18: 11-pointed star
    8: Hexagon               19: 12-pointed star
    9: Septagon              20: 13-pointed star
    10: Octagon              21: 14-pointed star
    11: Nine-gon             22: Plus with gap
    23: Vertical line
    24: Cross (X) with gap
    'LABEL'    Label string for symbol, up to 24 char.
    """
    # Get image descriptor
    id = im.Desc.Dict
    out = {id["ctype"][0]:[0.0], id["ctype"][1]:[0.0], \
           'MINOR AX': [0.0], 'MAJOR AX': [0.0], 'POSANG': [0.0], 'STARTYPE':[3.0], \
           'LABEL': ['                        '], \
           'NumFields': 8, 'Table name': 'AIPS ST', '_status': [0]}
    return out
    # end newRow

def PWriteCirc (sttab, im, center, radius, err):
    """ Write an entry for drawing a circle
    
    sttab  = Python Table object, must be open with write enabled
    im     = Obit Image on which to attach ST Table
    center = [x,y] pixels
    radius = radius in pixels
    err    = Python Obit Error/message stack
    """
    ################################################################
    # Check
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return None
    # Get image descriptor
    id = im.Desc.Dict
    # Get row
    row = newRow(im)
    # Convert pixels to positions
    pos = ImageDesc.PGetPos(im.Desc, center, err)
    if err.isErr:
        printErrMsg(err, "Error converting pixel location to position")
    row[id["ctype"][0]] = [pos[0]]
    row[id["ctype"][1]] = [pos[1]]
    row['MAJOR AX']  = [radius * abs(id["cdelt"][0])]
    row['MINOR AX']  = row['MAJOR AX']
    row['POSANG']    = [0.0]
    row['STARTYPE']  = [3.0]
    row['LABEL']     = ["    "]
    # Write
    sttab.WriteRow(-1,row, err)
    if err.isErr:
        printErrMsg(err, "Error Writing ST table")
    # end PWriteCirc

