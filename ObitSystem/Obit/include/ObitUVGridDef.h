/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
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
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** I/O status */
ObitIOStatus myStatus;
/** Width of convolving kernel in cells */
olong convWidth;
/** Number of of tabulated points per cell in convfn */
olong convNperCell;
/** Center pixel (1-rel) Beam */
olong icenxBeam, icenyBeam;
/** Center pixel (1-rel) Image */
olong icenxImage, icenyImage;
/** Size of beam plane (pixels) */
olong nxBeam, nyBeam;
/** Size of image plane (pixels) */
olong nxImage, nyImage;
/** Start channel (1-rel) in uv data to grid */
olong startChann;
/* Number of channels in uv data to grid */
olong numberChann;
/** Scaling from wavelength to cells for u, v, w at reference frequency */
ofloat UScale, VScale, WScale;
/** -2pi*Position shift parameters for x, y,z */
ofloat dxc, dyc, dzc;
/** rotation to apply in imaging */
ofloat rotate;
/** guardband as fraction of the grid size */
ofloat guardband;
/** max, min baseline lengths (wavelengths) */
ofloat blmax, blmin;
/** Prenormalization peak of Beam, used to normalize images */
ofloat BeamNorm;
/** Is this to be a Beam image? */
gboolean doBeam;
/** Need to apply 3D imaging rotations? */ 
gboolean do3Dmul;
/** 3D rotation matrix for u,v,w */
ofloat URot3D[3][3];
/** 3D rotation matrix for x,y,z */
ofloat PRot3D[3][3];
/** uvdata grid, typed as ofloat but complex as (real,imag) */
ObitCArray *grid;
/** Gridding convolution table */
ObitFArray *convfn;
/** ObitFFT for making a beam */
ObitFFT *FFTBeam;
/** ObitFFT for making image */
ObitFFT *FFTImage;
/** "X" gridding correction function for beam */
ObitFArray *xCorrBeam;
/** "Y" gridding correction function for beam */
ObitFArray *yCorrBeam;
/** "X" gridding correction function for image */
ObitFArray *xCorrImage;
/** "Y" gridding correction function for image */
ObitFArray *yCorrImage;
/** Number of threads (elements in threadArgs)  */
olong nThreads;
/** Array of FT Function structures  */
gpointer **threadArgs;
/** nThreads-1 work gridding arrays */
ObitCArray **workGrids;
