# Make mask file from image VL table
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2019
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
import Image, Table, ImageDesc, OErr

def MakeMask(im, vlver, maskFile, err, minSNR=5.0):
    """
    Create a CLEAN mask file for Obit Imaging tasks CLEANFile
    Radius of entries = major axis size.    

    im       = Image with VL table
    vlver    = VL Table version, generate with FndSou
    maskFile = name for output mask file
    err      = Obit error/message stack
    minSNR   = Min. SNR to use
    """
    scale=1.0
    # Open output file
    fd = open(maskFile,'w')

    vltab=im.NewTable(Table.READONLY,'AIPS VL',vlver,err)
    maxscale = scale/abs(im.Desc.Dict['cdelt'][0])
    nrow=vltab.Desc.Dict['nrow']
    vltab.Open(Table.READONLY,err)
    OErr.printErr(err)
    for irow in range(1,nrow+1):
        row = vltab.ReadRow(irow,err)
        if row['PEAK INT'][0]/row['I RMS'][0]>minSNR:
            ra =  ImageDesc.PRA2HMS(row['RA(2000)'][0])
            dec = ImageDesc.PDec2DMS(row['DEC(2000)'][0])
            size = min(50, int(row['MAJOR AX'][0]*maxscale+0.5))
            line = "%s %s %d\n"%(ra,dec,size)
            fd.write(line)
            
    vltab.Close(err)
    OErr.printErr(err)
    fd.close()
# end MakeMask
