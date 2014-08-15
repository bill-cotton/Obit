/* $Id$        */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITCUDAUTIL_H 
#define OBITCUDAUTIl_H 

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDAUtil.h
 *
 * C callable CUDA GPU routines
 * Routines abort on error
 * Real GPU/CUDA versions are called if HAVE_GPU==1, else stubbed versions
 */

/*---------------Public functions---------------------------*/
#if HAVE_GPU==1  /* GPU? Real versioins */
/* Public: Set device */
void ObitCUDASetGPU (int cuda_device);

/* Public: Reset device */
void ObitCUDAResetGPU ();

/* Public: synchronize GPU */
void ObitCUDADeviceSynchronize (int* event);

/* Public: Create stream */
int* ObitCUDAStreamCreate ();

/* Public: Destroy stream */
void ObitCUDAStreamDestroy (int* stream);

/* Public: Create event */
int* ObitCUDAEventCreate ();

/* Public: Destroy event */
void ObitCUDAEventDestroy (int* event);

/* Public: record event */
void ObitCUDAEventRecord (int* event, int *stream);

/* Public: synchronize event */
void ObitCUDAEventSynchronize (int* event);

/* Public: Allocate locked host memory */
float* ObitCUDAUtilAllocHost (int memsize);

/* Public: Deallocate locked host memory */
void ObitCUDAUtilFreeHost (float *host);

/* Public: Allocate locked Device memory */
float* ObitCUDAUtilAllocGPU (int memsize);

/* Public: Deallocate Device memory */
void ObitCUDAUtilFreeGPU (float *GPU);

/* Public: Copy Host to GPU memory */
void ObitCUDAUtilHost2GPU(float *GPU, float *host, int memsize, int* stream);

/* Public: Copy GPU to Host memory */
void ObitCUDAUtilGPU2Host(float *host, float *GPU, int memsize, int* stream);

#else  /* No GPU - stubb */
/* Public: Set device */
void ObitCUDASetGPU (int cuda_device)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDASetGPU */

/* Public: Reset device */
void ObitCUDAResetGPU ()
{
  g_error("GPU/CUDA not implemented");
} /* end  ObitCUDAResetGPU */

/* Public: synchronize GPU */
void ObitCUDADeviceSynchronize (int* event)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDADeviceSynchronize */

/* Public: Create stream */
int* ObitCUDAStreamCreate ()
{
  g_error("GPU/CUDA not implemented");
  return NULL;
} /* end  ObitCUDAStreamCreate */

/* Public: Destroy stream */
void ObitCUDAStreamDestroy (int* stream)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAStreamDestroy */

/* Public: Create event */
int* ObitCUDAEventCreate ()
{
  g_error("GPU/CUDA not implemented");
  return NULL;
} /* end ObitCUDAEventCreate */

/* Public: Destroy event */
void ObitCUDAEventDestroy (int* event)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAEventDestroy */

/* Public: record event */
void ObitCUDAEventRecord (int* event, int *stream)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAEventRecord */

/* Public: synchronize event */
void ObitCUDAEventSynchronize (int* event)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAEventSynchronize */

/* Public: Allocate locked host memory */
float* ObitCUDAUtilAllocHost (int memsize)
{
  g_error("GPU/CUDA not implemented");
  return NULL;
} /* end ObitCUDAUtilAllocHost */

/* Public: Deallocate locked host memory */
void ObitCUDAUtilFreeHost (float *host)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAUtilFreeHost */

/* Public: Allocate locked Device memory */
float* ObitCUDAUtilAllocGPU (int memsize)
{
  g_error("GPU/CUDA not implemented");
  return NULL;
} /* end ObitCUDAUtilAllocGPU  */

/* Public: Deallocate Device memory */
void ObitCUDAUtilFreeGPU (float *GPU)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAUtilFreeGPU */

/* Public: Copy Host to GPU memory */
void ObitCUDAUtilHost2GPU(float *GPU, float *host, int memsize, int* stream)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAUtilHost2GPU */

/* Public: Copy GPU to Host memory */
void ObitCUDAUtilGPU2Host(float *host, float *GPU, int memsize, int* stream)
{
  g_error("GPU/CUDA not implemented");
} /* end ObitCUDAUtilGPU2Host */
#endif /* HAVE_GPU */
#endif /* OBITFCUDAUtil_H */ 
