/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your option) any later version.              */
/*;                                                                   */
/*;  This program is distributed in the hope that it will be useful,  */
/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*;  GNU General Public License for more details.                     */
/*;                                                                   */
/*;  You should have received a copy of the GNU General Public        */
/*;  License along with this program; if not, write to the Free       */
/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*;  MA 02139, USA.                                                   */
/*;                                                                   */
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUV structure               */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitUVGridDef.h
 * ObitUVGrid structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Weighting sums */
odouble wtSums[3];
/** Width of convolving kernel in cells */
olong convWidth;
/** Number of of tabulated points per cell in convfn */
olong convNperCell;
/** Size of uv grid (pixels) */
olong nuGrid, nvGrid;
/** Weighting box radius */
olong WtBox;
/** Weighting function index */
olong WtFunc;
/** Number of visibilities out of the inner 90% of the grid */
olong numberBad;
/** Robust weighting parameter */
ofloat Robust;
/**  Weighting power */
ofloat WtPower;
/** Scaling from wavelength to cells for u, v at reference frequency */
ofloat UScale, VScale;
/** -sigma,uu,vv,uv for taper (in cells) */
ofloat sigma1, sigma2, sigma3;
/** max, min baseline lengths (wavelengths) */
ofloat blmax, blmin;
/** Robust temperance value */
ofloat temperance;
/** Weight scaling (normalization) factor */
ofloat wtScale;
/** Noise increase factor */
ofloat noiseFactor;
/** weight grid  */
ObitFArray *wtGrid;
/** count grid, as ofloat  */
ObitFArray *cntGrid;
/** Gridding convolution table */
ObitFArray *convfn;
