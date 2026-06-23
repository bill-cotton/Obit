/* $Id:  $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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
#ifndef CUDABEAMINTERP_H 
#define CUDABEAMINTERP_H

#include "CUDASkyGeom.h"
#include "CUDAFArray.h"

/**
 * \file CUDABeamInterp.h
 * CUDABeamInterp does Lagrangian interpolation of positions in an CUDAFArray
 * This is limited CUDA implementation of ObitFInterpolate
 */

/*--------------Class definitions-------------------------------------*/
/** CUDABeamInterp Class structure. */
#if HAVE_GPU==1  /* Have a GPU? */
#ifndef CUDABEAMINTERPDEF_H // Prevent multiple definitions
#define CUDABEAMINTERPDEF_H
typedef struct {
#include "CUDABeamInterpDef.h"   /* this class definition */
} CUDABeamInterp;
#endif /* CUDABEAMINTERPDEF_H */ 
#endif /* HAVE_GPU */

#include "ObitCUDAUtil.h"
#include "CUDAFArray.h"
#include "CUDASkyGeom.h"


/*---------------Public functions---------------------------*/
/* Have to include code here to get it into the same file as kernals */
#if IS_CUDA==1  /* CUDA code */
/* CUDA ONLY  */
// CUDA includes
#include <helper_cuda.h>
#include <helper_functions.h> 
/*---------------Private function prototypes----------------*/
/** Private: Set convolution kernal */
extern "C"
__device__ void SetConvKernal (float Target, int naxis, int hwidth, 
	   		       const float *denom, int *Start, float *Kernal);

/*----------------------Public C functions---------------------------*/

/*----------------------Public CUDA functions---------------------------*/
/**
 * Interpolate value at requested pixel in a plane of an n(>=2)-D array.
 * Interpolation between planes is not supported.
 * CUDA callable
 * \param in    The object to interpolate
 * \param pixel Pixel location (1-rel) in planes and which plane.
 *              Should have number of dimensions equal to in.
 * \return value, magic blanked if invalid
 */
__device__ float CUDABeamInterpPixel (const CUDABeamInterp __restrict__ *in, float *pixel)
{
  CUDAFArray *myArray = (CUDAFArray *)in->d_myArray;
  float fblank = myArray->fblank;
  float value = fblank;
  float sum, sumwt, wty, wt;
  float xKernal[10], yKernal[10], *data;
  int i, j, xStart, yStart, iwid, indx, planeOff, iplane, iprod;

  /* Must be inside array */
  iplane = 1;
  iprod = 1;
  for (i=0; i<myArray->ndim; i++) {
    if (i>1) {
      iplane *= pixel[i] * iprod;  /* How many planes deep? */
      iprod *= myArray->naxis[i];
    }
    if ((pixel[i]<1.0) || (pixel[i] > myArray->naxis[i])) {
      /* out of bounds - return blank */
      return fblank;
    }
  }

  /* Set convolving x, y kernals */
  SetConvKernal (pixel[0], in->nx, in->hwidth, in->denom, &xStart, xKernal);
  SetConvKernal (pixel[1], in->ny, in->hwidth, in->denom, &yStart, yKernal);

  /* Local versions of things */
  data = myArray->d_array;
  iwid = 1 + 2 * in->hwidth;

  /* Offset to start of plane */
  planeOff = (iplane-1) * in->nx * in->ny;
    
  /* Zero sums */
  sum   = 0.0;
  sumwt = 0.0;
 
  /* Loop over data summing values times convolving weights */
  for (j=0; j<iwid; j++) {
    wty = yKernal[j];
    indx = planeOff + xStart + ((yStart + j) * in->nx);
    for (i=0; i<iwid; i++) {
      if (data[indx] != fblank) {
	wt = xKernal[i] * wty;
	sumwt += wt;
	sum   += data[indx] * wt;
      } 
      indx++;
    }
  }

  /* normalize sum if not excessive blanking */
  if (sumwt > 0.90) {
    value = sum / sumwt;
    return value;
  } else return fblank; 

} /* end CUDABeamInterpPixel */

/**
 * Interpolate value at image position in a plane of an n(>=2)-D array.
 * Interpolation between planes is not supported.
 * MUST have attached image descriptor (CUDABeamInterpDesc)
 * CUDA callable
 * \param in    The object to interpolate
 * \param rotar Rotation angle (rad)
 * \param pos   Position (deg) in image, higher dimensions give plane number
 *              Should have number of dimensions equal to in.
 * \return value, magic blanked if invalid
 */
extern "C"
__device__ float CUDABeamInterpPos (CUDABeamInterp *in, float rotar, double *pos)
{
  CUDAFArray *myArray = (CUDAFArray *)in->d_myArray;
  float fblank = myArray->fblank;
  float pixel[5]={1,1,1,1,1};

  /* Determine pixel */
  if (CUDAPositionXYpix(pos, in->d_myDesc, rotar, pixel)!=0) {
    return fblank;
  }

  pixel[3] = (int)(pos[2]+0.5);  /* plane number */

  /* Interpolate */		 
  return CUDABeamInterpPixel (in, pixel);
} /* end CUDABeamInterpPos */

/**
 * Interpolate value at requested pixel in 1-D array.
 * CUDA callable
 * \param in    The object to interpolate
 * \param pixel Pixel location (1-rel) in array
 * \return value, blanked if invalid
 */
extern "C"
__device__ float CUDABeamInterp1D (CUDABeamInterp *in, float pixel)
{
  CUDAFArray *myArray = (CUDAFArray *)in->d_myArray;
  float fblank = myArray->fblank;
  float value = fblank;
  float sum, sumwt, wt;
  float xKernal[10], *data;
  int i, xStart, iwid, indx;

  /* Must be inside 1-D array */
  if (myArray->ndim>2) return fblank; /* too many dimensions */
  if (pixel<0.5) return fblank; /* not in array */
  if (pixel>(myArray->naxis[0]+0.5)) return fblank; /* not in array */

  SetConvKernal (pixel, in->nx, in->hwidth, in->denom, &xStart, xKernal);
    
  /* Local versions of things */
  data    = myArray->d_array;
  iwid    = 1 + 2 * in->hwidth;

  /* Zero sums */
  sum   = 0.0;
  sumwt = 0.0;
 
  /* Loop over data summing values times convolving weights */
  for (i=0; i<iwid; i++) {
    indx = xStart + i;
    if (data[indx] != fblank) {
      wt = xKernal[i];
      sumwt += wt;
      sum   += data[indx] * wt;
    } 
  }

  /* normalize sum if not excessive blanking */
  if (sumwt > 0.50) {
    value = sum / sumwt;
    return value;
  } 

  /* No luck - return fblank */
  return fblank;
} /* end CUDABeamInterp1D */

/*---------------Private functions--------------------------*/
/**
 * Set Lagrangian interpolation kernal taking into account ends of the grid.
 * \param  Pixel  Which is the desired pixel (1-rel)?
 * \param  naxis  Number of pixels on axis 
 * \param  hwidth Half width of convolution kernal
 * \param  denom  Reciprocals of Lagrangian denominators
 * \param  Start  [out] first pixels in array for convolution
 * \param  Kernal [out] convolving kernal.
 */
extern "C"
__device__  void SetConvKernal (float Pixel, int naxis, int hwidth, 
			       const float *denom, int *Start, float *Kernal)
{
  float prod, xx;
  int ipos, i, j, cen, iwid;

  /* fractional pixel */
  ipos = Pixel + 0.5;
  iwid = hwidth*2 + 1;

  /* set first pixel */
  cen = ipos - hwidth;
  cen = MAX (1, MIN (cen, (naxis-iwid+1)));
  /* make 0 rel */
  cen = cen - 1;
  *Start = cen; /* returned version */

  /* set "x" at first pixel to 1.0 */
  xx = Pixel - cen;

  /* compute interpolating kernal */
  for (j= 0; j<iwid; j++) { /* loop 50 */
    prod = denom[j];
    for (i= 0; i<iwid; i++) { /* loop 30 */
      if (i != j) prod *= (xx - (i+1));
    } /* end loop  L30:  */;
    Kernal[j] = prod;
  } /* end loop  L50:  */;
} /* end SetConvKernal */

#else /* Not CUDA */
#endif /* IS_CUDA */
#endif /* CUDABEAMINTERP_H */ 
