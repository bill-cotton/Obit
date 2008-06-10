/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*  Define the basic components of the ObitOTFGrid structure          */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitOTFGridDef.h
 * ObitOTFGrid structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** I/O status */
ObitIOStatus myStatus;
/** convolving function type */
olong convType;
/** Width of convolving kernel in cells */
olong convWidth;
/** Number of of tabulated points per cell in convfn */
olong convNperCell;
/** Center pixel (1-rel) Image */
olong icenxImage, icenyImage;
/** Size of image plane (pixels) */
olong nxImage, nyImage;
/** Scaling from position offset to cells */
ofloat XScale, YScale;
/** Minimum summed gridding weight */
ofloat minWt;
/** Convolving function parameters */
ofloat convParm[10];
/** Additional Gaussian equivalent beam size (FWHM) from convolution, in pixels */
ofloat addBM;
/** data grid */
ObitFArray *grid;
/** Weight grid */
ObitFArray *gridWt;
/** Gridding convolution table */
ObitFArray *convfn;
/** Working array of x (RA-like)  positions, one per detector in array */
ofloat *xpos;
/** Working array of y (Dec-like)  positions, one per detector in array */
ofloat *ypos;
/** Project type code */
ObitOTFProj Proj;
/** RA of projection tangent point (deg) */
ofloat raProj;
/** Dec of projection tangent point (deg) */
ofloat decProj;
/** Beam size (Gaussian sigma) in pixels */
ofloat beamSize;
/** Fitted Beam size FWHM in deg */
ofloat fitBeamSize;
/** Clipping level */
ofloat clip;
