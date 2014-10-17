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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#define IS_CUDA 1
#include "ObitCUDAUtil.h"
#include "CUDAFInterpolate.h"
#include "CUDAFArray.h"
#include <stdlib.h>

/**
 * \file CUDAFInterpolate.c
 * ObitFInterpolate class function definitions.
 * This class supports 1 and 2-D interpolation in ObitFArrays using 
 * Lagrange interpolation.
 * This is a limited subset of ObitFInterpolate.
 * Parallelism is invokes by multiple cally to interpolation functions 
 * rather than multiple interpolators, i.e. the functions maintain state 
 * rather than the interpolator.
 */


/*-------------------- CUDA Kernals ------------------------------------*/
#if HAVE_GPU==1  /* CUDA code */
// CUDA includes
#include <helper_cuda.h>
#include <helper_functions.h> 
/**
 * Interpolate pixels in a 2D float Array
 * block =  nx/16, ny/16, threads in block = 16 x 16 [x][y]
 * \param  g_finterp    interpolator
 * \param  g_Xpix       array of x pixels for g_out_arr
 * \param  g_Ypix       array of y pixels for g_out_arr
 * \param  yoff         offset in y
 * \param  g_out_arr    output 2D ObitFArray
 */
__global__ void floatInterpolateKernel(CUDAFInterpolate* __restrict__ g_fInterp, 
           CUDAFArray* g_XPix,  CUDAFArray* g_YPix, int Yoff,
	   float* __restrict__ g_out_arr)
{
    int i        = threadIdx.x + blockIdx.x*blockDim.x;       // i (x=fastest) index
    int j        = Yoff + threadIdx.y + blockIdx.y*blockDim.y; // j (y=slowest) index
    int nx       = g_XPix->naxis[0];
    int ny       = g_XPix->naxis[1];
    int indx     = i + j*nx;                                     // g_Xpix/g_Ypix index
    int jndx     = i + (threadIdx.y + blockIdx.y*blockDim.y)*nx; // output array index
    float pixel[2]; 

    // No more than actual number of output pixels
    if ((i>nx) || (j>ny)) return;

    // pixel numbers
    pixel[0] = g_XPix->d_array[indx];
    pixel[1] = g_YPix->d_array[indx];
    g_out_arr[jndx] = CUDAFInterpolatePixel (g_fInterp, pixel);
    //g_out_arr[jndx] = i;
 } // end floatInterpolateKernel

/*---------------Private function prototypes----------------*/
/** Private: Set convolution kernal */
__device__ void SetConvKernal (float Target, int naxis, int hwidth, 
			   float *denom, int *Start, float *Kernal);
static void InterpolateWithStreams(
            CUDAFInterpolate *interp, CUDAFArray *XPix, CUDAFArray *YPix,
            float *outBuffer, int prtLv);

/*----------------------Public C functions---------------------------*/
/**
 * Constructor from values.
 * C callable
 * \param array   The CUDAFarray to be interpolated.
 * \param hwidth  half-width in pixels of interpolation kernal
 * \param nstream Number of streams to use
 * \return the new object.
 */
extern "C"
CUDAFInterpolate* newCUDAFInterpolateCreate (CUDAFArray *array, int hwidth, int nstream)
{
  CUDAFInterpolate* out;
  int iwid, i, j, memsize;
  float prod;

  /* Create/init output structure */
  memsize = sizeof(CUDAFInterpolate);
  out = (CUDAFInterpolate*) ObitCUDAUtilAllocHost(memsize);
  /* GPU structure */
  memsize = sizeof(CUDAFInterpolate);
  out->d_FInterpolate = ObitCUDAUtilAllocGPU(memsize);

  /* Attach array */
  out->h_myArray = array;
  out->d_myArray = array->d_FArray;

  /* Kernal width */
  out->hwidth  = max (1, min (4, hwidth));

  /* Get array size info */
  out->nx     = array->naxis[0];
  out->ny     = array->naxis[1];

  /* Init Lagrangian denominators for hwidth */
  iwid = 1 + (2*out->hwidth);
  for (j= 1; j<=iwid; j++) {
    prod = 1.0;
    for (i= 1; i<=iwid; i++) {
      if (i != j) prod = prod * (j - i);
    } 
    out->denom[j-1] = 1.0 / prod;
  } 

  /* Allocate sStream stuff */
  out->nstream   = nstream;
  memsize        = out->nstream*sizeof(cudaStream_t*);
  out->stream    = (cudaStream_t*)ObitCUDAUtilAllocHost(memsize);
  memsize        = out->nstream*sizeof(cudaEvent_t*);
  out->cycleDone = (cudaEvent_t*)ObitCUDAUtilAllocHost(memsize);
  memsize        = out->nstream*sizeof(float*);
  out->d_data    = (float**)ObitCUDAUtilAllocHost(memsize);
  out->nxBuff = out->nyBuff = 0;  /* No GPU buffers yet */
  for (i=0; i<out->nstream; ++i)   {
    out->stream[i]    = ObitCUDAStreamCreateCUDA();
    out->cycleDone[i] = ObitCUDAEventCreateCUDA();
    ObitCUDAEventRecordCUDA (out->cycleDone[i], out->stream[i]);
  }

  /* copy basic structure to device */
  memsize = sizeof(CUDAFInterpolate);
  ObitCUDAUtilHost2GPU (out->d_FInterpolate, (float*)out, memsize, NULL);

 return out;
} /* end newCUDAFInterpolateCreate */

/**
 * Interpolate value at requested pixel in a plane of an n(>=2)-D array.
 * Interpolation between planes is not supported.
 * C callable
 * \param in       The object to interpolate
 * \param XPix     Array of input x pixels in outImage
 * \param YPix     Array of input y pixels in outImage
 * \param outImage Output pixel array
 */
extern "C"
void CUDAFInterpolateImage (CUDAFInterpolate *in, CUDAFArray *XPix, CUDAFArray *YPix,
			   float *outBuffer)
{
#if HAVE_GPU==1  /* CUDA code */
   // Put available fast memory in L1 cache - no apparent effect
   checkCudaErrors(cudaFuncSetCacheConfig(floatInterpolateKernel, cudaFuncCachePreferL1));

   /* Process with streams */
   InterpolateWithStreams(in, XPix, YPix, outBuffer, 0);
#endif /* HAVE_GPU */
    return;
} /* end CUDAInterpolateImage */

/**
 * Deallocates member objects.  Also Zaps FArray being interpolated.
 * C callable
 * \param  in The object to deallocate.
 */
extern "C"
void CUDAFInterpolateZap (CUDAFInterpolate *in)
{
	int i;
  if (in) {
    /* CUDA stuff */
    for (i=0; i<in->nstream; ++i)    {
      if ((in->stream)&&(in->stream[i]))       ObitCUDAStreamDestroyCUDA(in->stream[i]);
      if ((in->cycleDone)&&(in->cycleDone[i])) ObitCUDAEventDestroyCUDA(in->cycleDone[i]);
      if ((in->d_data)&&(in->d_data[i]))       ObitCUDAUtilFreeGPU(in->d_data[i]);
    }
    if (in->stream)       ObitCUDAUtilFreeHost((float*)in->stream);
    if (in->cycleDone)    ObitCUDAUtilFreeHost((float*)in->cycleDone);
    if (in->d_data)       ObitCUDAUtilFreeHost((float*)in->d_data);
    if (in->h_myArray)    CUDAFArrayZap(in->h_myArray);
    if (in->d_FInterpolate) ObitCUDAUtilFreeGPU((float*)in->d_FInterpolate);
    ObitCUDAUtilFreeHost((float*)in);
  }

} /* end CUDAFInterpolateZap */

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
__device__ float CUDAFInterpolatePixel (CUDAFInterpolate *in, float *pixel)
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

} /* end CUDAFInterpolatePixel */

/**
 * Interpolate value at requested pixel in 1-D array.
 * CUDA callable
 * \param in    The object to interpolate
 * \param pixel Pixel location (1-rel) in array
 * \return value, blanked if invalid
 */
__device__ float CUDAFInterpolate1D (CUDAFInterpolate *in, float pixel)
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
} /* end CUDAFInterpolate1D */

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
__device__  void SetConvKernal (float Pixel, int naxis, int hwidth, 
			       float *denom, int *Start, float *Kernal)
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

/*-------------------- Processing loops ------------------------------------*/
/**
 * Interpolation using a GPU.
 * Multiple streams are used to overlap I/O and computation,
 * each call divides the data into streams_used pieces.
 * \param  interp        Interpolator
 * \param  XPix          Array of input x pixels in outImage
 * \param  YPix          Array of input y pixels in outImage
 * \param  outBuffer     Host array for output
 * \param  prtLv         Print level, 5=>much directly printed
 */
static void InterpolateWithStreams(
            CUDAFInterpolate *interp, CUDAFArray *XPix, CUDAFArray *YPix,
            float *outBuffer, int prtLv)
{

    int last_stream, current_stream = 0;
    int nrowPass, npass = interp->nstream;
    int nx = XPix->naxis[0];   // Number of cells in a row
    int ny = XPix->naxis[1];   // Number of rows
    int memsize, dorow, off, Yoff;
    dim3 numBlocks, thPerBlock;

    // Divide rows (ny) by npass
    nrowPass = 1 + (ny-1)/npass;
    if (prtLv>=5) printf ("Start\n");
 
    // Do processing in a loop
    //
    // Note: All memory commands are processed in the order  they are issued,
    // independent of the stream they are enqueued in. Hence the pattern by
    // which the copy and kernel commands are enqueued in the stream
    // has an influence on the achieved overlap.

    for (int i=0; i<npass; ++i) {
        int next_stream = (current_stream + 1) % interp->nstream;
	int prev_stream = current_stream - 1;
	if (prev_stream<0) prev_stream = interp->nstream-1;
	off = next_stream*nrowPass*nx;       // Offset in data buffers 
        memsize = nrowPass*nx*sizeof(float); // Amount of data to read
  	if (prtLv>=5) printf ("\n\nLoop %d prev %d current %d next %d\n", 
            i, prev_stream,current_stream,next_stream );

	// Ensure that processing and copying of the previous cycle has finished
	if (i>0) {
  	  off = prev_stream*nrowPass*nx;  // Offset in data buffers 
	  if (prtLv>=5) printf ("download prev_stream %d off %d\n",prev_stream,off);
	  checkCudaErrors(cudaMemcpyAsync(&outBuffer[off],
					  interp->d_data[prev_stream],
					  memsize,
					  cudaMemcpyDeviceToHost,
					  interp->stream[prev_stream]));
	  if (prtLv>=5) printf ("sync prev_stream %d loop %d\n",prev_stream, i);
	  cudaEventSynchronize(interp->cycleDone[prev_stream]);
	}

        // Process current
	if (prtLv>=5) printf ("Process, nrow, %d stream %d \n",
           nrowPass, current_stream);
	// make sure to do all rows
	if (i==npass-1) dorow = ny-i*nrowPass;
	else            dorow = nrowPass;

	// package work
	numBlocks.x  = 1+(nx-1)/16; numBlocks.y  = 1 + (dorow-1)/16;
	thPerBlock.x = 16;          thPerBlock.y = 16;
	Yoff = i*nrowPass;

         floatInterpolateKernel<<<numBlocks, thPerBlock, 0, interp->stream[current_stream]>>>(
	      (CUDAFInterpolate *)interp->d_FInterpolate, 
              (CUDAFArray*)XPix->d_FArray, (CUDAFArray*)YPix->d_FArray, 
              Yoff, interp->d_data[current_stream]);

	// make sure previous frame done
	if (i>0) cudaEventSynchronize(interp->cycleDone[prev_stream]);

	last_stream    = current_stream;
        current_stream = next_stream;
    } /* end loop */

    /* Data from last pass */
    if (prtLv>=5) printf ("sync last_stream %d \n",last_stream);
    cudaEventSynchronize(interp->cycleDone[last_stream]);
    
    // last piece may be uneven size
    memsize = nx*dorow*sizeof(float);
    off = last_stream*nrowPass*nx;       // Offset in data buffers
    if (prtLv>=5) printf ("download last_stream %d off %d\n",last_stream,off);
    checkCudaErrors(cudaMemcpyAsync(&outBuffer[off],
				    interp->d_data[last_stream],
				    memsize,
				    cudaMemcpyDeviceToHost,
				    interp->stream[last_stream]));
   if (prtLv>=5) printf ("Finish\n");
    cudaDeviceSynchronize();

    return;
} /* end processWithStreams */
#endif /* HAVE_GPU */
