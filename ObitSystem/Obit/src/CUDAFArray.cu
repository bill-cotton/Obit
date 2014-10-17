/* $Id: $   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
/*  Define the CUDAFArray (GPU version of ObitFArray) class            */
/**
 * \file CUDAFArray.c
 * implementation of limited GPU version of ObitFArray class
 */
#define IS_CUDA 1
#include "CUDAFArray.h"
#include "ObitCUDAUtil.h"
#include <stdlib.h>
#include <helper_cuda.h>
#include <helper_functions.h> 

/**
 * Creates an CUDAFArray of a specified geometry.
 * Called from c
 * \param ndim  Number of dimensions desired, if <=0 data array not created.
 *              maximum value = MAXFARRAYDIM.
 * \param naxis Dimensionality along each axis
 * \return the new object.
 */
extern "C"
CUDAFArray* CUDAFArrayCreate (int ndim, int *naxis)
{
  CUDAFArray* out;
  int memsize, i, size;

  /* Create basic structure - first host */
  memsize = sizeof(CUDAFArray);
  out = (CUDAFArray*) ObitCUDAUtilAllocHost(memsize);
  /* GPU structure */
  memsize = sizeof(CUDAFArray);
  out->d_FArray = ObitCUDAUtilAllocGPU(memsize);

  /* copy geometry */
  out->ndim = ndim;
  size = 1; /* total size */
  for (i=0; i<ndim; i++) {
    out->naxis[i] = max (1, min(naxis[i],524288));  /* Not too big */
    size *= out->naxis[i]; /* total size */
  }

  /* create GPU data array - add a bit extra */
  memsize = (size+10)*sizeof(float);
  out->d_array = (float*)ObitCUDAUtilAllocGPU(memsize);
  /* Host version */
  out->h_array = (float*)ObitCUDAUtilAllocHost(memsize);

  /* Blanking value */
  out->fblank = CUDAMagicF();

  /* copy basic structure to device */
  memsize = sizeof(CUDAFArray);
  ObitCUDAUtilHost2GPU ((float*)out->d_FArray, (float*)out, memsize, NULL);

  return out;
} /* end CUDAFArrayCreate */

/**
 * Destroys a CUDAFArray
 * Called from c
 * \param in   Object to delete
 */
extern "C"
void CUDAFArrayZap (CUDAFArray *in)
{
  if (in) {
    if (in->h_array)  ObitCUDAUtilFreeHost((float*)in->h_array);
    if (in->d_array)  ObitCUDAUtilFreeGPU((float*)in->d_array);
    if (in->d_FArray) ObitCUDAUtilFreeGPU((float*)in->d_FArray);
    ObitCUDAUtilFreeHost((float*)in);
  }
} /* end CUDAFArrayZap */

