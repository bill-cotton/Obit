/* $Id$    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2021                                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITCUDAGRID_H 
#define OBITCUDAGRID_H 
#include "ObitCUDAGridInfoDef.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDAGrid.h
 *
 * Define interface to (primitive) CUDA routines
 * UV data gridding using a GPU.
 * Real GPU/CUDA versions are called if HAVE_GPU==1, else stubbed versions
 * 
 */

/*---------------Public functions---------------------------*/
#if HAVE_GPU==1  /* GPU? Real versions */
/* cuda version */
#if IS_CUDA==1
/* Public: Allocate Structures */
extern "C"
void ObitCUDAGridAlloc (CUDAGridInfo *gpuInfo, long nfacet, long nchan, long nif, 
			long nVisPIO, long lenvis, long nplane);

/* Public: Initialize gridding */
extern "C"
void ObitCUDAGridInit (CUDAGridInfo *gpuInfo);

/* Public: Grid a bufferload of visibilities */
extern "C"
void ObitCUDAGridGrid (CUDAGridInfo *in);

/* Public: Fold negative u to positive */
extern "C"
void ObitCUDAGridFlip (CUDAGridInfo *in);

/* Public: Copy vis grid to locked host memory */
extern "C"
void ObitCUDAGrid2CPU (CUDAGridInfo *in, long ifacet);

/* Public: Shutdown gridding */
extern "C"
void ObitCUDAGridShutdown (CUDAGridInfo *gpuinfo);
#else /* non cuda */
/* Public: Allocate Structures */
void ObitCUDAGridAlloc (CUDAGridInfo *gpuInfo, long nfacet, long nchan, long nif, 
			long nVisPIO, long lenvis, long nplane);

/* Public: Initialize gridding */
void ObitCUDAGridInit (CUDAGridInfo *gpuInfo);

/* Public: Grid a bufferload of visibilities */
void ObitCUDAGridGrid (CUDAGridInfo *in);

/* Public: Fold negative u to positive */
void ObitCUDAGridFlip (CUDAGridInfo *in);

/* Public: Copy vis grid to locked host memory */
void ObitCUDAGrid2CPU (CUDAGridInfo *in, long ifacet);

/* Public: Shutdown gridding */
void ObitCUDAGridShutdown (CUDAGridInfo *gpuinfo);
#endif

#else  /* No GPU - stub */
/* Public: Initialize gridding */
/*extern "C"*/
void ObitCUDAGridInit (CUDAGridInfo *gpuInfo, long nfacet, long nchan, long nVisPIO)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridInit */

/* Public: Grid a bufferload of visibilities */
/*extern "C"*/
void ObitCUDAGridGrid (CUDAGridInfo *in, ObitErr *err)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridInit */

/* Public: Fold negative u to positive */
/*extern "C"*/
void ObitCUDAGridFlip (CUDAGridInfo *in, ObitErr *err)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridFlip */

/* Public: Copy vis grid to FFT buffer via locked memory */
/*extern "C"*/
void ObitCUDAGrid2CPU (CUDAGridInfo *in, olong ispec, ObitErr *err)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGrid2CPU */

/* Public: Shutdown GPU gridding */
/*extern "C"*/
void ObitCUDAGridShutdown (CUDAGridInfo *gpuinfo)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridDFTShutdown */
#endif /* HAVE_GPU */
#endif /* OBITCUDAGRID_H */ 
