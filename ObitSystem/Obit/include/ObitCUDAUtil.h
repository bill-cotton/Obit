/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014-2026                                          */
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
#define OBITCUDAUTIL_H 

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDAUtil.h
 *
 * C callable CUDA GPU routines
 * Routines abort on error
 * Real GPU/CUDA versions are called if HAVE_GPU==1, else stubbed versions
 */

/*---------------Public functions---------------------------*/
/* Magic floating value */
float CUDAMagicF (void);

#if HAVE_GPU==1  /* GPU? Real versions */
/* Public: Check device number >0=OK, <0 = invalid*/
int ObitCUDACheckGPU (int cuda_device);

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

#if IS_CUDA==1  /* CUDA code */
/* Public: Create stream */
cudaStream_t ObitCUDAStreamCreateCUDA ();

/* Public: Destroy stream */
void ObitCUDAStreamDestroyCUDA (cudaStream_t stream);

/* Public: Create event */
cudaEvent_t ObitCUDAEventCreateCUDA ();

/* Public: Destroy event */
void ObitCUDAEventDestroyCUDA (cudaEvent_t event);

/* Public: record event */
void ObitCUDAEventRecordCUDA (cudaEvent_t event, cudaStream_t stream);

/* Public: synchronize event */
void ObitCUDAEventSynchronizeCUDA (cudaEvent_t event);

/* Public: Allocate locked host memory */
extern "C"
float* ObitCUDAUtilAllocHost (int memsize);

/* Public: Deallocate locked host memory */
extern "C"
void ObitCUDAUtilFreeHost (float *host);

/* Public: Allocate locked Device memory */
extern "C"
float* ObitCUDAUtilAllocGPU (int memsize);

/* Public: Deallocate Device memory */
extern "C"
void ObitCUDAUtilFreeGPU (float *GPU);

/* Public: Copy Host to GPU memory */
extern "C"
void ObitCUDAUtilHost2GPU(float *GPU, float *host, int memsize, int* stream);

/* Public: Copy Host anything to GPU memory */
extern "C"
void ObitCUDAUtilHostAny2GPU(void *GPU, void *host, int memsize, int* stream);

/* Public: Copy GPU to Host memory */
extern "C"
void ObitCUDAUtilGPU2Host(float *host, float *GPU, int memsize, int* stream);

/* Public: Copy anything GPU to Host memory */
extern "C"
void ObitCUDAUtilGPUAny2Host(void *host, void *GPU, int memsize, int* stream);

/* Public: How much device global memory in bytes */
extern "C"
size_t ObitCUDAUtilMemory(int cuda_device);

#else /* not CUDA */
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

/* Public: pointer info type */
void ObitCUDAUtilPointerType(void *ptr, char *label);

#endif /* IS_CUDA */
#else  /* No GPU - stubb */
/* Don't need to even define */
#endif /* HAVE_GPU */

#endif /* OBITFCUDAUtil_H */ 
