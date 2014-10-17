/* $Id:  $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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
#ifndef CUDAFINTERPOLATE_H 
#define CUDAFINTERPOLATE_H 

#include "CUDAFArray.h"

/**
 * \file CUDAFInterpolate.h
 * CUDAFInterpolate does Lagrangian interpolation of positions in an CUDAFArray
 * This is limited CUDA implementation of ObitFInterpolate
 */

/*--------------Class definitions-------------------------------------*/
/** CUDAFInterpolate Class structure. */
typedef struct {
#include "CUDAFInterpolateDef.h"   /* this class definition */
} CUDAFInterpolate;


/*---------------Public functions---------------------------*/
#if IS_CUDA==1  /* CUDA code */
/* CUDA ONLY (not C callable) */
/** Public: Constructor from value. */
extern "C"
CUDAFInterpolate* 
newCUDAFInterpolateCreate (CUDAFArray *array, int hwidth, int nstream);

/** Public: Interpolate Image */
extern "C"
void CUDAFInterpolateImage (CUDAFInterpolate *in, 
			    CUDAFArray *XPix, CUDAFArray *YPix,
			    float *outBuffer);

/** Public: Destroy */
extern "C"
void CUDAFInterpolateZap(CUDAFInterpolate *in);

/** Public: Interpolate Pixel in 2D array */
__device__ float CUDAFInterpolatePixel (CUDAFInterpolate *in, float *pixel);

/** Public: Interpolate value in 1- array */
__device__ float CUDAFInterpolate1D (CUDAFInterpolate *in, float pixel);
#else /* Not CUDA */
/** Public: Constructor from value. */
CUDAFInterpolate* 
newCUDAFInterpolateCreate (CUDAFArray *array, int hwidth, int nstream);

/** Public: Interpolate Image */
void CUDAFInterpolateImage (CUDAFInterpolate *in, 
			    CUDAFArray *XPix, CUDAFArray *YPix,
			    float *outBuffer);

/** Public: Destroy */
void CUDAFInterpolateZap(CUDAFInterpolate *in);

#endif /* IS_CUDA */
#endif /* CUDAFINTERPOLATE_H */ 
