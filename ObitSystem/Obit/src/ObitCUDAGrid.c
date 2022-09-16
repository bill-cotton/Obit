/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2020,2021                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "Obit.h"
#include "ObitCUDAGrid.h"


/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDAGrid.c
 * ObitCUDAGrid stub versions for non GPU use
 * This class is derived from the Obit base class.
 */
#if HAVE_GPU!=1  /* Not GPU? Stub versions */


/* Public: Initialize gridding */
void ObitCUDAGridInit (CUDAGridInfo *gpuInfo)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridInit */

/* Public: Grid a bufferload of visibilities */
void ObitCUDAGridGrid (CUDAGridInfo *in)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridGrid */

/* Public: Fold negative u to positive */
void ObitCUDAGridFlip (CUDAGridInfo *in)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridFlip */

/* Public: Copy vis grid to locked host memory */
void ObitCUDAGrid2CPU (CUDAGridInfo *in, long ifacet)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGrid2CPU */

/* Public: Shutdown gridding */
void ObitCUDAGridShutdown (CUDAGridInfo *gpuinfo)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAGridDFTShutdown */

#endif 


