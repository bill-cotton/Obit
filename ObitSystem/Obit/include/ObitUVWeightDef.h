/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2017                                          */
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
/** Number of IFs */
olong numIF;
/** Separate Robust factors per IF? */
gboolean RobustperIF;
/** Robust weighting parameter array per IF */
ofloat *Robust;
/**  Weighting power */
ofloat WtPower;
/** Scaling from wavelength to cells for u, v at reference frequency */
ofloat UScale, VScale;
/** -sigma,uu,vv,uv for taper (in cells) array per IF */
ofloat *sigma1, *sigma2, *sigma3;
/** Minimum weight for inner taper */
ofloat minInnerWt;
/** -sigma,uu,vv,uv for inverse taper (in cells) array per IF */
ofloat *isigma1, *isigma2, *isigma3;
/** max, min baseline lengths (wavelengths) */
ofloat blmax, blmin;
/** Robust temperance value array */
ofloat *temperance;
/** Weight scaling (normalization) factor array per IF */
ofloat *wtScale;
/** Noise increase factor */
ofloat noiseFactor;
/** weight grid one per IF  */
ObitFArray **wtGrid;
/** count grid, as ofloat one per IF */
ObitFArray **cntGrid;
/** Gridding convolution table */
ObitFArray *convfn;
