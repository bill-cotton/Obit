"""
Python utility package for Zernike polynomials

This utility package supplies Zernike terms and derivatives.
"""

# Python utility package for Zernike polynomials
import Obit
# $Id: ZernikeUtil.py,v 1.1 2006/10/26 21:36:42 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2006
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

def PZernike(n, x, y):
    """ Determine Zernike coefficients for position (x,y)

    returns  array of n efficients (by order)
    n       = highest order number (1-rel) between 1 and 36 supported 
    x       = X coordinate on unit circle
    y       = Y coordinate on unit circle
    """
    ################################################################
    # Error check
    if n<1 or n>36:
        raise RuntimeError,"unsupported Zernike order:"+str(n)
    out = []   # initialize output
    for i in range(0,n):
        val = Obit.Zernike(i, x, y);
        out.append(val)
    return out
    # end PZernike

def PZernikeGradX(n, x, y):
    """ Determine Zernike gradient in x coefficients for position (x,y)

    returns  array of n efficients (by order)
    n       = highest order number (1-rel) between 1 and 36 supported 
    x       = X coordinate on unit circle
    y       = Y coordinate on unit circle
    """
    ################################################################
    # Error check
    if n<1 or n>36:
        raise RuntimeError,"unsupported Zernike order:"+str(n)
    out = []   # initialize output
    for i in range(0,n):
        val = Obit.Zernike(i, x, y);
        out.append(val)
    return out
    # end PZernikeGradX

def PZernikeGradY(n, x, y):
    """ Determine Zernike gradient in y coefficients for position (x,y)

    returns  array of n efficients (by order)
    n       = highest order number (1-rel) between 1 and 36 supported 
    x       = X coordinate on unit circle
    y       = Y coordinate on unit circle
    """
    ################################################################
    # Error check
    if n<1 or n>36:
        raise RuntimeError,"unsupported Zernike order:"+str(n)
    out = []   # initialize output
    for i in range(0,n):
        val = Obit.Zernike(i, x, y);
        out.append(val)
    return out
    # end PZernikeGradY

def PZernikePolar(n, rho, phi):
    """ Determine Zernike coefficients for Polar coordinates rho and phi

    returns  array of n efficients (by order)
    n       = highest order number (1-rel) between 1 and 36 supported 
    rho     = radial distance on unit circle
    phi     = azimuthal coordinate on unit circle (radian)
    """
    ################################################################
    # Error check
    if n<1 or n>36:
        raise RuntimeError,"unsupported Zernike order:"+str(n)
    out = []   # initialize output
    for i in range(0,n):
        val = Obit.Zernike(i, rho, phi);
        out.append(val)
    return out
    # end PZernikePolar

