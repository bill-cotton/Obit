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
#if HAVE_GPU==1  /* CUDA code */
/* Public: Set device */
/**
 * Assign a GPU
 * \param cuda_device GPU number to use
 */
extern "C"
void ObitCUDASetGPU (int cuda_device)
{
  checkCudaErrors(cudaSetDevice(cuda_device));
} /* end ObitCUDASetGPU */

/**
 * Reset device 
 */
extern "C"
void ObitCUDAResetGPU ()
{
  cudaDeviceReset();
} /* end ObitCUDAResetGPU */

/**
 * Synchronize the device
 * \param event  that to be waited for, not really used
 */
extern "C"
void ObitCUDADeviceSynchronize (int* event)
{
  cudaDeviceSynchronize();
} /* end ObitCUDADeviceSynchronize */

/**
 * Create a processing stream
 * Native CUDA types
 * \return stream
 */
cudaStream_t ObitCUDAStreamCreateCUDA ()
{
  cudaStream_t out=NULL;
  checkCudaErrors(cudaStreamCreate(&out));
  return out;
} /* end ObitCUDAStreamCreateCUDA */

/**
 * Destroy a stream
 * Native CUDA types
 * \param stream stream to destroy
 */
void ObitCUDAStreamDestroyCUDA (cudaStream_t stream)
{
  cudaStreamDestroy(stream);
} /* end ObitCUDAStreamDestroyCUDA */

/**
 * Create an event
 * Native CUDA types
 * \return event
 */
cudaEvent_t ObitCUDAEventCreateCUDA ()
{
  cudaEvent_t out=NULL;
  checkCudaErrors(cudaEventCreate(&out));
  return out;
} /* end ObitCUDAEventCreateCUDA */

/**
 * Destroy an event
 * Native CUDA types
 * \param event to be destroyed
 */
void ObitCUDAEventDestroyCUDA (cudaEvent_t event)
{
  cudaEventDestroy(event);
} /* end ObitCUDAEventDestroyCUDA */

/**
 * Associate an event with a stream, waiting for completion
 * Native CUDA types
 * \param event  to wait for
 * \param stream stream
 */
void ObitCUDAEventRecordCUDA (cudaEvent_t event, cudaStream_t stream)
{
  cudaEventRecord(event, stream);
} /* end ObitCUDAEventRecordCUDA */

/**
 * Wait for an event defined by ObitCUDAEventRecordCUDA
 * Native CUDA types
 * \param event  to wait for
 */
void ObitCUDAEventSynchronizeCUDA (cudaEvent_t event)
{
  cudaEventSynchronize(event);
} /* end ObitCUDAEventSynchronizeCUDA */

/**
 * Create a processing stream
 * \return stream as int*
 */
extern "C"
int* ObitCUDAStreamCreate ()
{
  int* out=NULL;
  checkCudaErrors(cudaStreamCreate((cudaStream_t*)&out));
  return out;
} /* end ObitCUDAStreamCreate */

/**
 * Destroy a stream
 * \param stream stream to destroy (as int*)
 */
extern "C"
void ObitCUDAStreamDestroy (int* stream)
{
  cudaStreamDestroy((cudaStream_t)stream);
} /* end ObitCUDAStreamDestroy */

/**
 * Create an event
 * \return event as int*
 */
extern "C"
int* ObitCUDAEventCreate ()
{
  int* out=NULL;
  checkCudaErrors(cudaEventCreate((cudaEvent_t*)&out));
  return out;
} /* end ObitCUDAEventCreate */

/**
 * Destroy an event
 * \param event to be destroyed (as int*)
 */
extern "C"
void ObitCUDAEventDestroy (int* event)
{
  cudaEventDestroy((cudaEvent_t)event);
} /* end ObitCUDAEventDestroy */

/**
 * Associate an event with a stream, waiting for completion
 * \param event  to wait for
 * \param stream stream
 */
extern "C"
void ObitCUDAEventRecord (int* event, int* stream)
{
  cudaEventRecord((cudaEvent_t)event, (cudaStream_t)stream);
} /* end ObitCUDAEventRecord */

/**
 * Wait for an event defined by ObitCUDAEventRecord
 * \param event  to wait for
 */
extern "C"
void ObitCUDAEventSynchronize (int* event)
{
  cudaEventSynchronize((cudaEvent_t)event);
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
  checkCudaErrors(cudaMallocHost(&out, memsize));
  return out;
} /* end ObitCUDAUtilAllocHost */

/**
 * Deallocate locked host memory
 * \param host memory pointer to free
 */
extern "C"
void ObitCUDAUtilFreeHost (float *host)
{
  cudaFreeHost(host);
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
  checkCudaErrors(cudaMalloc(&out, memsize));
  return out;
} /* end ObitCUDAUtilAllocGPU */

/**
 * Deallocate Device memory
 * \param GPU memory pointer to free
 */
extern "C"
void ObitCUDAUtilFreeGPU (float *GPU)
{
  cudaFree(GPU);
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
  if (stream!=NULL) {
    checkCudaErrors(cudaMemcpyAsync(GPU, host, memsize, 
      cudaMemcpyHostToDevice, (cudaStream_t)stream));
  } else {
    checkCudaErrors(cudaMemcpyAsync(GPU, host, memsize, cudaMemcpyHostToDevice, 0));
  }
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
  if (stream!=NULL) {
    checkCudaErrors(cudaMemcpyAsync(host,GPU,  memsize, 
      cudaMemcpyDeviceToHost, (cudaStream_t)stream));
  } else {
    checkCudaErrors(cudaMemcpyAsync(host, GPU, memsize, cudaMemcpyDeviceToHost, 0));
  }
} /* end ObitCUDAUtilGPU2Host */
#endif /* HAVE_GPU */


/**
 * Returns Obit magic blanking float value
 * This is adopted from AIPS and correcponds to the string 'INDE'
 * \return float magic value
 */
float CUDAMagicF (void)
{
  static union FBLANKequiv {
    char string[4];
    float fblank;
  } FBLANK;
  FBLANK.string[0] = 'I'; 
  FBLANK.string[1] = 'N'; 
  FBLANK.string[2] = 'D'; 
  FBLANK.string[3] = 'E'; 
  
  return FBLANK.fblank;
} /* end CUDAMagicF */
#endif /* OBITFCUDAUtil_H */ 

