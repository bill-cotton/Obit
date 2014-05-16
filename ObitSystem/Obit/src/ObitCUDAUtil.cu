/* $Id:  $        */
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

#define IS_CUDA 1

// includes, project
#include <helper_cuda.h>
#include <helper_functions.h>  // helper for shared that are common to CUDA SDK samples

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDAUtil.h
 *
 * C callable CUDA GPU routines
 * Routines abort on error
 */

/*---------------Public functions---------------------------*/
/* Public: Set device */
/**
 * Assign a GPU
 * \param cuda_device GPU number to use
 */
extern "C"
void ObitCUDASetGPU (int cuda_device)
{
#if HAVE_GPU==1  /* CUDA code */
  checkCudaErrors(cudaSetDevice(cuda_device));
#endif /* HAVE_GPU */
} /* end ObitCUDASetGPU */

/**
 * Reset device 
 */
extern "C"
void ObitCUDAResetGPU ()
{
#if HAVE_GPU==1  /* CUDA code */
  cudaDeviceReset();
#endif /* HAVE_GPU */
} /* end ObitCUDAResetGPU */

/**
 * Wait for an event
 * \param event  that to be waited for
 */
extern "C"
void ObitCUDADeviceSynchronize (int* event)
{
#if HAVE_GPU==1  /* CUDA code */
  cudaDeviceSynchronize();
#endif /* HAVE_GPU */
} /* end ObitCUDADeviceSynchronize */

/**
 * Create a processing stream
 * \return stream as int*
 */
extern "C"
int* ObitCUDAStreamCreate ()
{
  int* out=NULL;
#if HAVE_GPU==1  /* CUDA code */
  checkCudaErrors(cudaStreamCreate((cudaStream_t*)&out));
#endif /* HAVE_GPU */
  return out;
} /* end ObitCUDAStreamCreate */

/**
 * Destroy a stream
 * \param stream stream to destroy (as int*)
 */
extern "C"
void ObitCUDAStreamDestroy (int* stream)
{
#if HAVE_GPU==1  /* CUDA code */
  cudaStreamDestroy((cudaStream_t)stream);
#endif /* HAVE_GPU */
} /* end ObitCUDAStreamDestroy */

/**
 * Create an event
 * \return event as int*
 */
extern "C"
int* ObitCUDAEventCreate ()
{
  int* out=NULL;
#if HAVE_GPU==1  /* CUDA code */
  checkCudaErrors(cudaEventCreate((cudaEvent_t*)&out));
#endif /* HAVE_GPU */
  return out;
} /* end ObitCUDAEventCreate */

/**
 * Destroy an event
 * \param event to be destroyed (as int*)
 */
extern "C"
void ObitCUDAEventDestroy (int* event)
{
#if HAVE_GPU==1  /* CUDA code */
  cudaEventDestroy((cudaEvent_t)event);
#endif /* HAVE_GPU */
} /* end ObitCUDAEventDestroy */

/**
 * Associate an event with a stream, waiting for completion
 * \param event  to wait for
 * \param stream stream
 */
extern "C"
void ObitCUDAEventRecord (int* event, int* stream)
{
#if HAVE_GPU==1  /* CUDA code */
  cudaEventRecord((cudaEvent_t)event, (cudaStream_t)stream);
#endif /* HAVE_GPU */
} /* end ObitCUDAEventRecord */

/**
 * Wait for an event defined by ObitCUDAEventRecord
 * \param event  to wait for
 */
extern "C"
void ObitCUDAEventSynchronize (int* event)
{
#if HAVE_GPU==1  /* CUDA code */
  cudaEventSynchronize((cudaEvent_t)event);
#endif /* HAVE_GPU */
} /* end ObitCUDAEventSynchronize */

/**
 * Allocate locked host memory 
 * \param memsize  size in bytes
 * \return pointer to memory
 */
extern "C"
float* ObitCUDAUtilAllocHost (int memsize)
{
  float *out=NULL;
#if HAVE_GPU==1  /* CUDA code */
  checkCudaErrors(cudaMallocHost(&out, memsize));
#endif /* HAVE_GPU */
  return out;
} /* end ObitCUDAUtilAllocHost */

/**
 * Deallocate locked host memory
 * \param host memory pointer to free
 */
extern "C"
void ObitCUDAUtilFreeHost (float *host)
{
#if HAVE_GPU==1  /* CUDA code */
  cudaFreeHost(host);
#endif /* HAVE_GPU */
} /* end ObitCUDAUtilFreeHost */

/**
 * Allocate device memory 
 * \param  size in bytes
 * \return pointer to memory
 */
extern "C"
float* ObitCUDAUtilAllocGPU (int memsize)
{
  float *out=NULL;
#if HAVE_GPU==1  /* CUDA code */
  checkCudaErrors(cudaMalloc(&out, memsize));
#endif /* HAVE_GPU */
  return out;
} /* end ObitCUDAUtilAllocGPU */

/**
 * Deallocate Device memory
 * \param GPU memory pointer to free
 */
extern "C"
void ObitCUDAUtilFreeGPU (float *GPU)
{
#if HAVE_GPU==1  /* CUDA code */
  cudaFree(GPU);
#endif /* HAVE_GPU */
} /* end ObitCUDAUtilFreeGPU */

/**
 * Copy data from host to GPU memory
 * \param GPU      GPU memory
 * \param host     locked host memory
 * \param memsize  size in bytes
 * \param stream   If non-NULL then stream pointer
 */
extern "C"
void ObitCUDAUtilHost2GPU(float *GPU, float *host, int memsize, int* stream)
{
#if HAVE_GPU==1  /* CUDA code */
  if (stream!=NULL) {
    checkCudaErrors(cudaMemcpyAsync(GPU, host, memsize, 
      cudaMemcpyHostToDevice, (cudaStream_t)stream));
  } else {
    checkCudaErrors(cudaMemcpyAsync(GPU, host, memsize, cudaMemcpyHostToDevice, 0));
  }
#endif /* HAVE_GPU */
} /* end ObitCUDAUtilHost2GPU */

/**
 * Copy data from GPU to host memory
 * \param host     locked host memory
 * \param GPU      GPU memory
 * \param memsize  size in bytes
 * \param stream   If non-NULL then stream pointer
 */
extern "C"
void ObitCUDAUtilGPU2Host(float *host, float *GPU, int memsize, int* stream)
{
#if HAVE_GPU==1  /* CUDA code */
  if (stream!=NULL) {
    checkCudaErrors(cudaMemcpyAsync(host,GPU,  memsize, 
      cudaMemcpyDeviceToHost, (cudaStream_t)stream));
  } else {
    checkCudaErrors(cudaMemcpyAsync(host, GPU, memsize, cudaMemcpyDeviceToHost, 0));
  }
#endif /* HAVE_GPU */
} /* end ObitCUDAUtilGPU2Host */

#endif /* OBITFCUDAUtil_H */ 
