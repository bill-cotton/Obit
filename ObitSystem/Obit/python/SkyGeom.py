""" Python Obit Sky Geometry Utility

This module contains utilities for Sky Geometry calculations
Also primary beam calculations
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007,2014
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

# Python utility package for  Sky Geometry 
import Obit, UVDesc
import math


def PShiftXY (ra, dec, rotate, shiftRA, shiftDec):
    """ Determine shift between two positions 

    Determine the shift in RA and Dec between two celestial positions.
    The shift is in (possibly) rotated coordinates.
     ra       Initial Right Ascension in deg.
     dec      Initial declination in deg.
     rotate   Rotation of field, to E from N, deg.
     shiftRA  Shifted Right Ascension in deg.
     shiftDec Shifted declination in deg.
     returns [xShift, yShift] shift in RA, dec in deg
    """
    ################################################################
    # Dummy final arguments for Swig interface.
    return Obit.SkyGeomShiftXY (ra, dec, rotate, shiftRA, shiftDec, [0.0], [0.0])
    # end PShiftXY


def PXYShift (ra, dec, xShift, yShift, rotate):
    """ Determine result of a shift to a position 

    Determine result of a shift applied to a celestial position.
    The shift is in (possibly) rotated coordinates.
     ra       Initial Right Ascension in deg.
     dec      Initial declination in deg.
     xShift   Shift from ra to shiftRA in deg.
     yShift   Shift from dec to shiftDec in deg.
     rotate   Rotation of field, to E from N, deg.
     return [shiftRA, shiftDec] Shifted RA, Dec in deg.
    """
    ################################################################
    return Obit.SkyGeomXYShift (ra, dec, xShift, yShift, rotate, [0.0], [0.0])
    # end PXYShift

    

def PNewPos (type, ra0, dec0, l, m):
    """ Returns astronomical coordinates given direction cosines, projection 

    Determines the coordinates (raout,decout) corresponding to a 
    displacement (l,m) given by the direction cosines from coordinates 
    (ra0,dec0).  the direction cosine l is assumed to be positive to 
    the east; m is positive to the north.  the routine works for 
    4 kinds of projective geometries and for celestial, ecliptic, or 
    galactic coordinate systems. 
    this subroutine always uses an accurate computation. 
    All angles in this subroutine are in radians. 
      type    projection type code e.g. "-SIN"
            Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
            projections anything else is linear 
      ra0    coordinate reference right ascension (longitude) 
      dec0   coordinate reference declination (latitude) 
      l      cosine angle of displacement to east 
      m      cosine angle of displacement to north
    Return [raout, decout, ierr]
      raout  Right ascension or longitude at (l,m) 
      decout declination or latitude at (l,m) 
      ierr   Error condition: 0 = ok, 1 = l,m crazy, 
             2 = bad type,  3 = answer undefined 
    """
    ################################################################
    return Obit.SkyGeomNewPos (type, ra0, dec0, l, m,[0.0], [0.0], [0])
    # end PNewPos

    

def PWorldPos(xpix, ypix, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, \
              type, xpos, ypos):
    """ accurate position for pixel coordinates 
    
    Determine accurate position for pixel coordinates.
      xpix    x pixel number  (RA or long without rotation)
      ypix    y pixel number  (dec or lat without rotation)
      xref    x reference coordinate value (deg)
      yref    y reference coordinate value (deg)
      xrefpix x reference pixel
      yrefpix y reference pixel
      xinc    x coordinate increment (deg)
      yinc    y coordinate increment (deg)
      rot     rotation (deg)  (from N through E)
      type    projection type code e.g. "-SIN"
    Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
    projections anything else is linear 
    Return [ierr, xpos, ypos]
      ierr  0 if successful otherwise: 1 = angle too large for projection;
      xpos  x (RA) coordinate (deg)
      ypos  y (dec) coordinate (deg)
    """
    ################################################################
    return Obit.SkyGeomWorldPos(xpix, ypix, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, \
              type,[0.0], [0.0])
    # end PWorldPos


def PCDpos(xpix, ypix, xref, yref, xrefpix, yrefpix, xinc, yinc, rot,\
           cd1, cd2, type):
    """ Position for pixel coordinates from IRAF style CD matrix 
    
    Determine accurate position for pixel coordinates from IRAF
    style CD matrix.  
    Note: xinc, yinc, and rot can be derived from cd1 and cd2 and 
    should be compatible with them.
      xpix    x pixel number  (RA or long without rotation)
      ypix    y pixel number  (dec or lat without rotation)
      xref    x reference coordinate value (deg)
      yref    y reference coordinate value (deg)
      xrefpix x reference pixel
      yrefpix y reference pixel
      xinc    x coordinate increment (deg)
      yinc    y coordinate increment (deg)
      rot     rotation (deg)  (from N through E)
      type    projection type code e.g. "-SIN"
    Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
    projections anything else is linear 
      cd1     first column of CD matrix
      cd2     second column of CD matrix
    Return [ierr, xpos, ypos]
      ierr  0 if successful otherwise: 1 = angle too large for projection;
      xpos  x (RA) coordinate (deg)
      ypos  y (dec) coordinate (deg)
    """
    ################################################################
    return Obit.SkyGeomCDpos(xpix, ypix, xref, yref, xrefpix, yrefpix, xinc, yinc, rot,\
           cd1, cd2, type, [0.0], [0.0])
    # end PCDpos



def PXYpix(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc,  \
           rot, type):
    """ Pixel coordinates for an RA and Dec
    
    Determine pixel for given coordinate
      xpos    x (RA) coordinate (deg)
      ypos    y (dec) coordinate (deg)
      xref    x reference coordinate value (deg)
      yref    y reference coordinate value (deg)
      xrefpix x reference pixel
      yrefpix y reference pixel
      xinc    x coordinate increment (deg)
      yinc    y coordinate increment (deg)
      rot     rotation (deg)  (from N through E)
      type    projection type code e.g. "-SIN"
    Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
    projections anything else is linear 
    Return [ierr, xpix, ypix]
      ierr   0 if successful otherwise: 1 = angle too large for projection
             2 = bad values
      xpix    x pixel number  (RA or long without rotation)
      ypix    y pixel number  (dec or lat without rotation)
    """
    ################################################################
    return Obit.SkyGeomXYpix(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, \
           rot, type, [0.0], [0.0])
    # end PXYpix

    

def PCDpix(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, \
           icd1, icd2, type):
    """ ixel coordinates for an RA and Dec from IRAF  style CD matrix. 
    
    Determine accurate pixel coordinates for an RA and Dec uses IRAF  
    style CD matrix. 
    Note: xinc, yinc, and rot can be derived from cd1 and cd2 and 
    should be compatible with them.
      xpos    x (RA) coordinate (deg)
      ypos    y (dec) coordinate (deg)
      xref    x reference coordinate value (deg)
      yref    y reference coordinate value (deg)
      xrefpix x reference pixel
      yrefpix y reference pixel
      xinc    x coordinate increment (deg)
      yinc    y coordinate increment (deg)
      rot     rotation (deg)  (from N through E)
      type    projection type code e.g. "-SIN"
            Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
            projections anything else is linear 
      cd1     first column of CD matrix
      cd2     second column of CD matrix
     Return [ierr, xpix, ypix]
      ierr   0 if successful otherwise: 1 = angle too large for projection
             2 = bad values
      xpix    x pixel number  (RA or long without rotation)
      ypix    y pixel number  (dec or lat without rotation)
    """
    ################################################################
    return Obit.SkyGeomCDpix(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, \
           icd1, icd2, type, [0.0], [0.0])
    # end PCDpix



def PWorldPosLM(dx, dy, xref, yref, xinc, yinc, rot, type):
    """ Position for pixel coordinates from  offsets from the reference position.
    
    Determine accurate position for pixel coordinates from offsets 
    from the reference position.
      dx      x coordinate offset  (RA or long) 
      dy      y coordinate offset  (dec or lat)
      xref    x reference coordinate value (deg)
      yref    y reference coordinate value (deg
      xinc    x coordinate increment (deg)
      yinc    y coordinate increment (deg)
      rot     rotation (deg)  (from N through E)
      type    projection type code e.g. "-SIN"
    Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
    projections anything else is linear 
    Return [ierr, xpos, ypos]
      ierr   0 if successful otherwise: 1 = angle too large for projection
      xpos    x (RA) coordinate (deg)
      ypos    y (dec) coordinate (deg)
    """
    ################################################################
    return Obit.SkyGeomWorldPosLM(dx, dy, xref, yref, xinc, yinc, rot, type, \
                                  [0.0], [0.0])
    # end PWorldPosLM



def PXYPixLM(xpos, ypos, xref, yref, xinc, yinc, rot, type, dx, dy):
    """ Coordinate offsets for an RA and Dec   
    
    Determine accurate coordinate offsets for an RA and Dec 
      xpos    x (RA) coordinate (deg)
      ypos    y (dec) coordinate (deg)
      xref    x reference coordinate value (deg)
      yref    y reference coordinate value (deg)
      xrefpix x reference pixel
      yrefpix y reference pixel
      xinc    x coordinate increment (deg)
      yinc    y coordinate increment (deg)
      rot     rotation (deg)  (from N through E)
      type    projection type code e.g. "-SIN"
    Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
    projections anything else is linear 
    Return [ierr, xpos, ypos]
      ierr   0 if successful otherwise: 1 = angle too large for projection
             2 = bad values 
      dx     x (RA) coordinate offset (deg)
      dy     y (dec) coordinate offset (deg)
    """
    ################################################################
    return Obit.SkyGeomXYPixLM(xpos, ypos, xref, yref, xinc, yinc, rot, type, \
                               [0.0], [0.0])
    # end PXYPixLM



def PBtoJ (ra, dec):
    """ Precess B1950 to J2000 coordinates
    
    Converts B1950 RA, Dec to J2000 
    Using method on page B42 of The Astronomical Almanac (1990 ed.)
      ra    in/out Right Ascension in degrees
      dec   in/out Declination in degrees
    Return [ra,dec]  J2000 position in deg.
    """
    ################################################################
    return Obit.SkyGeomBtoJ ([ra], [dec])
    # end PBtoJ



def PJtoB (ra, dec):
    """ Precess J2000 to B1950 coordinates 

    Converts J2000 RA, Dec to B1950  
    Using method on page B42 of The Astronomical Almanac (1990 ed.)
      ra    in/out Right Ascension in degrees
      dec   in/out Declination in degrees
    Return [ra,dec]  B1950 position in deg.
    """
    ################################################################
    return Obit.SkyGeomJtoB ([ra], [dec])
    # end PJtoB



def PEq2Gal (RALong, DecLat):
    """ Convert Equatorial (B1950) to Galactic coordinates  

    Converts Convert Equatorial (B1950)to Galactic coordinates
      RALong    in/out Right Ascension/longitude in degrees
      DecLat    in/out Declination.latitude in degrees
    Return [RALong, DecLat]  Galactic coordinates in deg.
    """
    ################################################################
    return Obit.SkyGeomEq2Gal ([RALong], [DecLat])
    # end PEq2Gal



def PGal2Eq (RALong, DecLat):
    """ Convert Galactic to Equatorial (B1950)  

    Converts Convert Galactic to Equatorial (B1950) coordinates
      RALong    in/out Right Ascension/longitude in degrees
      DecLat    in/out Declination.latitude in degrees
    Return [RALong, DecLat] Equatorial coordinates (B1950)  in deg.
    """
    ################################################################
    return Obit.SkyGeomGal2Eq ([RALong], [DecLat])
    # end PGal2Eq



def PEq2Ec (RALong, DecLat, epoch):
    """ Convert Equatorial to Ecliptic coordinates 

    Converts Convert Equatorial to Ecliptic coordinates
      RALong    in/out Right Ascension/longitude in degrees
      DecLat    in/out Declination.latitude in degrees
      epoch     Epoch of the coordinates to transform
    Return [RALong, DecLat] Ecliptic coordinates in deg.
    """
    ################################################################
    return Obit.SkyGeomEq2Ec ([RALong], [DecLat], epoch)
    # end PEq2Ec



def PEc2Eq (RALong, DecLat, epoch):
    """ Convert Ecliptic to Equatorial 

    Converts Convert Ecliptic to Equatorial coordinates
      RALong    in/out Right Ascension/longitude in degrees
      DecLat    in/out Declination.latitude in degrees
      epoch     Epoch of the coordinates to transform
    Return [RALong, DecLat] Equatorial coordinates in deg.
    """
    ################################################################
    return Obit.SkyGeomEc2Eq ([RALong], [DecLat], epoch)
    # end PEc2Eq



def PRADec2Zern (ra, dec, xshift, yshift):
    """ Projection to Zernike plane 

    Converts from celestial coordinates, expressed in terms of  
    a reference position and a shift from this position to coordinates  
    in a plane for fitting an Ionospheric phase screen model  
    consisting of Zernike polynomials.  The output coordinates are  
    normalized to unity at a 10 deg radius from the reference position.  
    The coordinates are projected onto a plane tangent to the sky at  
    the reference position.  
     ra      Right Ascention of reference position (deg) 
     dec     Declination of reference position (deg) 
     xshift  Shift in X (RA) to desired position (deg) 
     yshift  Shift in Y (Dec) to desired position (deg) 
    Return [xzer, yzer, ierr]
     xzer    x-coordinate on Zernike plane 
     yzer    y-coordinate on Zernike plane 
     ierr    0 ok, 1 out of range 
    """
    ################################################################
    return Obit.SkyGeomRADec2Zern (ra, dec, xshift, yshift, [0.0], [0.0], [0])
    # end PRADec2Zern

def PPBPoly (Desc, RA, Dec, cutoff=None):
    """ Compute VLA beam shape from a fitted polynomial
    
    Compute primary beam at RA, Dec with pointing given in Desc
      Desc   = Data (UV or Image) descriptor with pointing position
               and frequency
      Ra     = Right Ascension in degrees
      Dec    = Declination in degrees
      cutoff = distance from pointing beyond which to set primary
               beam to zero
      Return Primary beam  power factor [0.0, 1.0].
    """
    ################################################################
    d = Desc.Dict    # Descriptor as dictionary
    freq = d["crval"][d["jlocf"]]   # Frequency
    rotate = d["crota"][d["jlocd"]] # Rotation
    # Get pointing position, obsra, obsdec if possible
    if "obsra" in d:
        obsra  = d["obsra"]
    else:
        obsra  = 0.0
    if "obsdec" in d:
        obsdec = d["obsdec"]
    else:
        obsdec  = 0.0
    if (obsra == 0.0) and (obsdec == 0.0):
        obsra  = d["crval"][d["jlocr"]]
        obsdec = d["crval"][d["jlocd"]]
    # Angle from obs posn
    (dra, ddec) =  PShiftXY (obsra, obsdec, rotate, RA, Dec)
    angle = math.sqrt (dra*dra+ddec*ddec)
    if cutoff and (angle>cutoff):
        pbfact = 0.0
    else:
        pbfact =  Obit.PBUtilPoly(angle, freq, 0.0)
    return pbfact
    # end  PPBPoly

def PPBKAT (Desc, RA, Dec, cutoff=None):
    """ Compute KAT-7 beam shape from a fitted polynomial
    
    Compute primary beam at RA, Dec with pointing given in Desc
      Desc   = Data (UV or Image) descriptor with pointing position
               and frequency
      Ra     = Right Ascension in degrees
      Dec    = Declination in degrees
      cutoff = distance from pointing beyond which to set primary
               beam to zero
      Return Primary beam  power factor [0.0, 1.0].
    """
    ################################################################
    d = Desc.Dict    # Descriptor as dictionary
    freq = d["crval"][d["jlocf"]]   # Frequency
    rotate = d["crota"][d["jlocd"]] # Rotation
    # Get pointing position, obsra, obsdec if possible
    if "obsra" in d:
        obsra  = d["obsra"]
    else:
        obsra  = 0.0
    if "obsdec" in d:
        obsdec = d["obsdec"]
    else:
        obsdec  = 0.0
    if (obsra == 0.0) and (obsdec == 0.0):
        obsra  = d["crval"][d["jlocr"]]
        obsdec = d["crval"][d["jlocd"]]
    # Angle from obs posn
    (dra, ddec) =  PShiftXY (obsra, obsdec, rotate, RA, Dec)
    angle = math.sqrt (dra*dra+ddec*ddec)
    if cutoff and (angle>cutoff):
        pbfact = 0.0
    else:
        pbfact =  Obit.PBUtilKAT(angle, freq, 0.0)
    return pbfact
    # end  PPBKAT

def PPBJinc (Desc, RA, Dec, antSize=25.0, cutoff=None):
    """ Compute Antenna beam assuming uniform illumination of an antenna
    
    Compute primary beam at RA, Dec with pointing given in Desc 
      Desc    = Data (UV or Image) descriptor with pointing position
                and frequency
      Ra      = Right Ascension in degrees
      Dec     = Declination in degrees
      antSize = Antenna diameter (m)
      cutoff  = distance from pointing beyond which to set primary
               beam to zero
      Return Primary beam power factor [0.0, 1.0].
    """
    ################################################################
    d = Desc.Dict    # Descriptor as dictionary
    freq = d["crval"][d["jlocf"]]   # Frequency
    rotate = d["crota"][d["jlocd"]] # Rotation
    # Get pointing position, obsra, obsdec if possible
    if "obsra" in d:
        obsra  = d["obsra"]
    else:
        obsra  = 0.0
    if "obsdec" in d:
        obsdec = d["obsdec"]
    else:
        obsdec  = 0.0
    if (obsra == 0.0) and (obsdec == 0.0):
        obsra  = d["crval"][d["jlocr"]]
        obsdec = d["crval"][d["jlocd"]]
    # Angle from obs posn
    (dra, ddec) =  PShiftXY (obsra, obsdec, rotate, RA, Dec)
    angle = math.sqrt (dra*dra+ddec*ddec)
    if cutoff and (angle>cutoff):
        pbfact = 0.0
    else:
        pbfact =  Obit.PBUtilJinc(angle, freq, antSize, 0.0)
    return pbfact
    # end  PPBJinc



