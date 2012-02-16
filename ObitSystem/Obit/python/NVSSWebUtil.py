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

# Python package for Web Utilities 
import Obit, Image, ImageDesc, SkyGeom, ImageInterp, OErr
# Number of images per 4 deg declination strip starting as (dec, #)
NVSSnumMap = [
     ( 0., 90), ( 4., 90), ( 8., 90), (12., 90), (16., 90), (20., 90), \
     (24., 90), (28., 90), (32., 80), (36., 80), (40., 72), (44., 72), \
     (48., 72), (52., 60), (56., 60), (60., 48), (64., 48), (68., 40), \
     (72., 32), (76., 32), (80., 24), (84., 16), (88., 8)
     ]

def NVSSFindFile(RA, Dec, equinox, err, stokes='I'):
    """ Determine the NVSS image best suited for a given position

    Returns name of image
    RA       = Right ascension as string ("HH MM SS.SS")
    Dec      = Declonation as string ("sDD MM SS.SS")
    equinox  = equinox of RA,Dec, 1950 or 2000
    err      = Python Obit Error/message stack
    stokes   = Stokes desired, 'I', 'Q' or 'U'
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return
    if not type(RA)==str:
        raise TypeError,"RA must be a string (HH MM SS.S)"
    if not type(Dec)==str:
        raise TypeError,"Dec must be a string (HH MM SS.S)"
    #
    ra  = ImageDesc.PHMS2RA(RA, sep=' ')
    dec = ImageDesc.PDMS2Dec(Dec, sep=' ')
    # Precess 1950 to 2000?
    if equinox==1950:
        (ra,dec) = SkyGeom.PBtoJ(ra,dec)

    # Must be north of -41 dec
    if dec<-41.:
        OErr.PLog(err, OErr.Error, "Dec MUST be > -41")
        return "No Image"

    # Maximum distance to "belong"
    tolDec = 2.0
    idec   = -1
    adec   = abs(dec)
    for i in range(0,23):
        if abs(adec-NVSSnumMap[i][0])<=tolDec:
            dd = int(NVSSnumMap[i][0]+0.01)
            idec = i
            break;
    # Find it?
    if idec<0:
        OErr.PLog(err, OErr.Error, "Failed on declination "+Dec)
        return "No Image"
    # Width of images in RA
    deltRa = 360.0 / NVSSnumMap[idec][1]
    racell = int(0.5 + ra/deltRa)  # Cell
    if racell>=NVSSnumMap[idec][1]: # Wrap?
        racell = 0
    # Derive name
    raCen = racell*deltRa / 15.0
    iRAH = int(raCen + 0.5)
    iRAM = int(60.0 * (raCen-iRAH) + 0.1)
    if iRAM<0.0:
        iRAH -= 1
        iRAM += 59
    if iRAH == -1:
        iRAH = 23
    # Sign of declination
    if dec>=-2.0:
        dsign = 'P'
    else:
        dsign = 'M'
    # Stokes type
    if stokes in ["Q","U"]:
        itype = 'C'
    else:
        itype = 'I'
    # put it together
    file = "%s%2.2i%2.2i%s%2.2i"%(itype,iRAH,iRAM,dsign,dd)
    return file
    # end NVSSFindFile

def NVSSPtFlux(RA, Dec, equinox, err, stokes='I',dir="/home/ftp/nvss/MAPS"):
    """ Determine the NVSS Flux density

    Returns flux density, None on failure
    RA       = Right ascension as string ("HH MM SS.SS")
    Dec      = Declonation as string ("sDD MM SS.SS")
    equinox  = equinox of RA,Dec, 1950 or 2000
    err      = Python Obit Error/message stack
    stokes   = Stokes desired, 'I', 'Q' or 'U'
    dir      = directory or url of directory
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if err.isErr: # existing error?
        return None
    #
    file = NVSSFindFile(RA, Dec, equinox, err, stokes=stokes)
    if err.isErr:
        OErr.printErr(err)
        return None
    path = dir+"/"+file
    #print "Image ",path
    img = Image.newPFImage("Image", path, 0, True, err, verbose=False)
    if err.isErr:
        OErr.printErr(err)
        return None
    interp = ImageInterp.PCreate("Interpolator", img, err)
    if err.isErr:
        OErr.printErr(err)
        return None
    ra  = ImageDesc.PHMS2RA(RA, sep=' ')
    dec = ImageDesc.PDMS2Dec(Dec, sep=' ')
    # Precess 1950 to 2000?
    if equinox==1950:
        (ra,dec) = SkyGeom.PBtoJ(ra,dec)
    # Plane
    plane = 1
    if stokes=='Q':
        plane = 2
    if stokes=='U':
        plane = 3
    # Interpolate
    value = interp.Value(ra, dec, err, plane=plane)
    if err.isErr:
        OErr.printErr(err)
        return None
    return value;
    # end NVSSPtFlux
