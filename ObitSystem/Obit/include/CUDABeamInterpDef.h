/* $Id: $                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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
/*  Define the basic components of the CUDABeamInterp structure     */
/*  This is intended to be included in a class structure definition   */
/**
 * \file CUDABeamInterpDef.h
 */
/* Fooey - may need local definition */
#if HAVE_GPU==1  /* Have a GPU? */
#ifndef CUDAIMAGEDESCDEF_H // Prevent multiple definitions
#define CUDAIMAGEDESCDEF_H
#define IM_MAXDIM 7       /* maximum array dimension */
#define IMLEN_VALUE 41    /* Maximum length of descriptor string value */
#define IMLEN_KEYWORD 21  /* Maximum length of descriptor keyword  */
typedef struct {
#include "CUDAImageDescDef.h"
} CUDAImageDesc;
#endif /* CUDAIMAGEDESCDEF_H */
#endif /* HAVE_GPU */

/** Structure in GPU memory */
float *d_BeamInterp;
/** Dimension of plane in myArray */
int nx, ny;
/** Half width of interpolation kernal */
int hwidth;
/** Reciprocals of Lagrangian denominators */
float denom[10];
/** CUDA stuff */
#if HAVE_GPU==1  /* Have a GPU? */
/** Array to be interpolated (Host) */
CUDAFArray *h_myArray;
/** Array to be interpolated (GPU) */
CUDAFArray *d_myArray;
/** Image descriptor */
CUDAImageDesc *d_myDesc;
#else  /* No GPU */
/** Array to be interpolated (Host) */
void* *h_myArray;
/** Array to be interpolated (GPU) */
void* *d_myArray;
/** Image descriptor */
void* *d_myDesc;
#endif /* HAVE_GPU */
