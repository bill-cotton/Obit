/* $Id: $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2021,2022                                          */
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

/*#include "ObitCUDAGrid.h"*/
#include "ObitCUDAUtil.h"
/* Public: pointer info type DAMN*/
void ObitCUDAUtilPointerType(void *ptr, char *label);

#include "ObitCUDAGrid.h"
#include "ObitCUDAGridInfoDef.h"
#include "ObitErr.h"
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDAGrid.cu
 * Primitive CUDA routines
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling routines.
 */

/*--------------- CUDA setup, GPU kernal  ----------------*/
/* This is a CUDA routine */
#if HAVE_GPU==1  /* CUDA code */
#define IS_CUDA 1

// includes, project
#include <helper_cuda.h>
#include <helper_functions.h>  // helper for shared that are common to CUDA SDK samples
#include "ObitCUDAGrid.h"
 
// functions
static void grid_processWithStreams(int streams_used, CUDAGridInfo* gridInfo, 
            cudaStream_t* stream, cudaEvent_t* cycleDone);
#endif /* HAVE_GPU */

// Atomic add on less capable systems
__device__ float atomicAddFloat(float* address, float val)
{
#if __CUDA_ARCH__ >=200
    return atomicAdd(address, val);
#else
    unsigned int* address_as_ul = (unsigned int*)address;
    union floatEquiv {
      unsigned int       iarg;
      float              farg;
    };
    union floatEquiv arg, assumed, old;
    old.iarg = *address_as_ul;
    do {
        assumed.iarg = old.iarg;	
	arg.farg = val+assumed.farg;
        old.iarg = atomicCAS(address_as_ul, assumed.iarg, arg.iarg);
    } while (assumed.iarg != old.iarg);
    return old.farg;
# endif
}

// use constant memory for some arrays
#define MAX_CM_CONVFN 2000  // Maximum size of convolution array
#define MAX_CM_FREQARR 8192 // Maximum size of channel frequency scaling array

__device__ __constant__ float cm_convfn[MAX_CM_CONVFN];  // hardcoded convolution function max
__device__ __constant__ float cm_freqArr[MAX_CM_FREQARR]; // hardcoded channel max      

// Grid has halfWidth extra columns in u added to avoid complication near u=0
//  The extra rows in the conjugate region need to be added to half plane before FFT.
__global__ void gridKernel(int current, int ifacet, int nvis, CUDAGridInfo* cudaInfo, float *d_grid_debug)
{
    GridInfo* gridInfo = cudaInfo->d_gridInfo;                // grid specifications
    FacetInfo* facetInfo = cudaInfo->d_facetInfo[ifacet];     // facet specifications
    long kvis         = threadIdx.x + blockIdx.x*blockDim.x;  // vis number
    long ichan        = threadIdx.y + blockIdx.y*blockDim.y;  // channel number
    long iplane       = gridInfo->d_freqPlane[ichan];         // plane number
    float *grid       = facetInfo->d_grid + iplane*facetInfo->sizeGrid;
    //float *rot        = facetInfo->rotUV;
    //float *shift      = facetInfo->shift;
    //float *convfn     = NULL; //gridInfo->d_convfn;
    float *freqArr    = NULL; //gridInfo->d_freqArr;
    float *vis_in     = gridInfo->d_data_in[current];
    long ivis         = kvis * gridInfo->lenvis;              // beginning of visibility
    int halfWidth     = gridInfo->convWidth/2;                // half width of convolution kernal
    int fullWidth     = gridInfo->convWidth;                  // full width of convolution kernal
    int convNperCell  = gridInfo->convNperCell;               // resolution of kernal
    int  nxo2         = facetInfo->nx/2;
    int  nyo2         = facetInfo->ny/2;
    //int  doBeam       = facetInfo->doBeam;
    //int  nrparm       = gridInfo->nrparm;
    //float guardu      = gridInfo->guardu*facetInfo->nx;  /* Guardband in cells */
    //float guardv      = gridInfo->guardv*facetInfo->ny;
    //float uscale      = facetInfo->uscale;
    //float vscale      = facetInfo->vscale;
    //float maxBL2      = facetInfo->maxBL*facetInfo->maxBL;
    //float minBL2      = facetInfo->minBL*facetInfo->minBL;
    //float bmTaper     = facetInfo->bmTaper;
    float sigma1      = gridInfo->d_sigma1[ichan];
    float sigma2      = gridInfo->d_sigma2[ichan];
    float sigma3      = gridInfo->d_sigma3[ichan];
    float *gp;
    float u,v,w, uu, vv, ww, bl2;
    float vr, vi, vw, vvr, vvi;
    float phase, c, s, phaseSign=1.0, freqFact;
    float *convu, *cconvu, *convv, cu, cv, tapeFact; //, ftemp;
    long iu, iv, icu, icv, jvis, lrow, icufn, icvfn;

    // in bounds?
    if ((ifacet>=gridInfo->nfacet) || (kvis>=nvis) || (ichan>=gridInfo->nchan)
         || (iplane>=gridInfo->nplane)) return;

    // pointers into constant memory
    //convfn  = cm_convfn;
    if (gridInfo->nchan<MAX_CM_FREQARR) freqArr = cm_freqArr;
    else                                freqArr = gridInfo->d_freqArr;

    jvis = ivis + gridInfo->nrparm + ichan*3;
    vw  = vis_in[jvis+2];
    if (vw<=0.0) return;  // datum flagged?

    //  Assume random parameters start with u,v,w
    u = vis_in[ivis];
    v = vis_in[ivis+1];
    w = vis_in[ivis+2];
    
    // rotate u,v,w for facet
    uu = u*facetInfo->rotUV[0] + v*facetInfo->rotUV[1] + w*facetInfo->rotUV[2];
    vv = u*facetInfo->rotUV[3] + v*facetInfo->rotUV[4] + w*facetInfo->rotUV[5];
    ww = u*facetInfo->rotUV[6] + v*facetInfo->rotUV[7] + w*facetInfo->rotUV[8];

    // Only gridding half plane, need to flip to other side?
    if (uu<0.0) {
      phaseSign = -1.0;
    } else { // no flip
      phaseSign = 1.0;
    }

    freqFact = phaseSign * freqArr[ichan];
    u = uu * freqFact;  // Scale u,v,w to channel
    v = vv * freqFact;
    w = ww * freqFact;

    // baseline limit
    bl2 = u*u+v*v+w*w;
    if ((bl2<facetInfo->minBL*facetInfo->minBL) || (bl2>facetInfo->maxBL*facetInfo->maxBL)) return;

    // Tapering
     if (sigma1!=0.0) {
 	tapeFact = __expf(u*u*(sigma2+facetInfo->bmTaper)+v*v*(sigma1+facetInfo->bmTaper)+u*v*(sigma3));
     } else if (facetInfo->bmTaper!=0.0) {
   	tapeFact = __expf(facetInfo->bmTaper*bl2);
     } else tapeFact=1.0;
     vw *= tapeFact;
 
    // visibility or (1,0) for beam
    if (facetInfo->doBeam) {
      vvr = vw; vvi = 0.0;  // for beams use (1 (weight),0) 
    } else {  // image
      vvr = vw*vis_in[jvis];
      vvi = phaseSign*vw*vis_in[jvis+1];
    }

    // position shift
    phase =  (u*facetInfo->shift[0] + v*facetInfo->shift[1] + w*facetInfo->shift[2]);
    // convert to cells
    u = u * facetInfo->uscale;
    v = v * facetInfo->vscale;

   // Data valid? within guardband
   if ((u<gridInfo->guardu*facetInfo->nx) && (fabs(v)<gridInfo->guardv*facetInfo->ny)) {

    // rotate vis to facet position
    __sincosf (phase, &s, &c);
    vr = c*vvr - s*vvi;
    vi = s*vvr + c*vvi;
    iu = lroundf(u);   // to grid cell number, u=zero in [halfWidth]
    iv = lroundf(v);

// DEBUG
//atomicAddFloat(&d_grid_debug[ifacet], sigma1); //DEBUG
//if ((ichan==1)) {
//   d_grid_debug[0]=sigma1; d_grid_debug[1]=sigma2; d_grid_debug[2]=sigma3; d_grid_debug[3]=bmTaper;
//}
//    atomicAddFloat(&d_grid_debug[100+iu+halfWidth], 1.0);
//    atomicAddFloat(&d_grid_debug[500+iv+nyo2], 1.0);
//    for (icv=0; icv<200; icv++) {d_grid_debug[icv*2]=gridInfo->d_convfn[icv]; d_grid_debug[1+icv*2]=cm_convfn[icv];}
//} // end  DEBUG

    // start in convolution function
    lrow = 2*(1 + nxo2 + halfWidth);   // length of grid row in floats
    icufn = lroundf((convNperCell*(iu + 0.5 - u)));
    convu = cm_convfn + icufn;
    icvfn = lroundf((convNperCell*(iv + 0.5 - v)));
    convv = cm_convfn + icvfn;
    cv = (*convv);
    iv += nyo2-halfWidth;
    for (icv=0; icv<fullWidth; icv++) {  // loop in v
      gp = grid + iu*2 + (iv+icv)*lrow;
      cconvu = convu;
      for (icu=0; icu<fullWidth; icu++) { // loop in u
        cu = (*cconvu)*cv;
	// atomic update of grid
 	atomicAddFloat(gp++, vr*cu);
	atomicAddFloat(gp++, vi*cu);
	cconvu += convNperCell;
      } // end u convolution loop
      convv += convNperCell;
      cv = (*convv);
      }  // end v convolution loop
   } // end data valid
} // end gridKernel

// Flip/add conjugate rows in grid
// block.x = vrow, threads in block = facet
__global__ void gridFlipKernel(CUDAGridInfo *cudaInfo, int ifacet, float *d_grid_debug)
{
   GridInfo* gridInfo   = cudaInfo->d_gridInfo;              // grid specifications
   FacetInfo* facetInfo = cudaInfo->d_facetInfo[ifacet];     // facet specifications
   long irow       = threadIdx.x + blockIdx.x*blockDim.x;    // row number
   long iplane     = threadIdx.y + blockIdx.y*blockDim.y;    // plane
   int  nx         = facetInfo->nx;
   int  ny         = facetInfo->ny;
   long halfWidth  = gridInfo->convWidth/2;    // half width of convolution kernal
   float *grid     = facetInfo->d_grid + iplane*facetInfo->sizeGrid;
   float *gxi, *gxo, *gci, *gco, xxo[2], cjo[2], xxi[2], cji[2];
   long iu, vc, lrow;

   // Check out of range 
   if ((irow>ny/2) || (iplane>=gridInfo->nplane)) return;

   vc = ny - irow;         // conjugate row number
   lrow = 2 + nx + 2*halfWidth;  // length of grid row 
   gci = grid + vc*lrow   + halfWidth*2;
   gxo = grid + irow*lrow + halfWidth*2;
   gxi = gxo; gco = gci;
   // loop over u columns
   for (iu=0; iu<=halfWidth; iu++) {
     // Read initial values from both sides of both rows 
     cji[0] = gci[0]; cji[1] = gci[1];
     cjo[0] = gco[0]; cjo[1] = gco[1];
     xxi[0] = gxi[0]; xxi[1] = gxi[1];
     xxo[0] = gxo[0]; xxo[1] = gxo[1];
     // update both row and conjugate row 
     gxo[0] = xxo[0] + cji[0];
     gxo[1] = xxo[1] - cji[1];
     gco[0] = cjo[0] + xxi[0];
     gco[1] = cjo[1] - xxi[1];
     gxo += 2; gco += 2; gxi -= 2; gci -= 2;
   }
} // end gridFlipKernel

#if HAVE_GPU==1  /* CUDA code */
/**
 * Copy cuda base structures to GPU
 * \param cudaInfo GPU gridding information
 */
extern "C"
void UpdateGPUBase (CUDAGridInfo *cudaInfo) {
  CUDAGridInfo *h_cudaInfo=NULL;
  GridInfo* h_gridInfo    = cudaInfo->h_gridInfo;      // host grid specifications
  FacetInfo** h_facetInfo = cudaInfo->h_facetInfo;     // host array of facet specifications
  GridInfo* d_gridInfo    = cudaInfo->d_gridInfo;      // device grid specifications
  FacetInfo** d_facetInfo = cudaInfo->d_facetInfo;     // device array of facet specifications
  int i, nfacet=h_gridInfo->nfacet, nchanIF=h_gridInfo->nchan*h_gridInfo->nif;
  size_t memsize;

  // Load constant memory
  memsize = ((h_gridInfo->convWidth+1)*h_gridInfo->convNperCell)*sizeof(float);
  cudaMemcpyToSymbol(cm_convfn,  h_gridInfo->d_convfn,  memsize);
  memsize = MIN(MAX_CM_FREQARR,h_gridInfo->nchan)*sizeof(float);
  cudaMemcpyToSymbol(cm_freqArr, h_gridInfo->d_freqArr, memsize);

  /* copy base structure to device */
  memsize = sizeof(CUDAGridInfo);
  checkCudaErrors(cudaMallocHost(&h_cudaInfo, memsize));
  memcpy (h_cudaInfo, cudaInfo, memsize);
  checkCudaErrors(cudaMemcpy(cudaInfo->d_base, h_cudaInfo, memsize, cudaMemcpyHostToDevice));  
  // DEBUG - check
 // memsize = 100;
 // checkCudaErrors(cudaMemcpy(cudaInfo->h_grid_debug, cudaInfo->d_base, memsize, cudaMemcpyDeviceToHost));  

  /* copy gridInfo to device */
  memsize = sizeof(GridInfo);
  checkCudaErrors(cudaMemcpy(d_gridInfo, h_gridInfo, memsize, cudaMemcpyHostToDevice));  

  /* arrays */
  memsize = nchanIF*sizeof(long);
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_freqPlane, h_gridInfo->h_freqPlane, memsize, cudaMemcpyHostToDevice));  
  memsize = nchanIF*sizeof(float);
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_freqArr, h_gridInfo->h_freqArr, memsize, cudaMemcpyHostToDevice));  
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_sigma1, h_gridInfo->h_sigma1, memsize, cudaMemcpyHostToDevice));  
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_sigma2, h_gridInfo->h_sigma2, memsize, cudaMemcpyHostToDevice));  
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_sigma3, h_gridInfo->h_sigma3, memsize, cudaMemcpyHostToDevice));  
  memsize = 1000*sizeof(float);
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_convfn, h_gridInfo->h_convfn, memsize, cudaMemcpyHostToDevice));  
  /* Saved data buffer pointers */
  memsize = (GRID_STREAM_COUNT*sizeof(float*));
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_data_in, h_gridInfo->h_d_data_in, memsize, cudaMemcpyHostToDevice));  

  /* copy to device */
  memsize = nfacet*sizeof(FacetInfo*);
  checkCudaErrors(cudaMemcpy(d_facetInfo, cudaInfo->h_d_facetInfo, memsize, cudaMemcpyHostToDevice));  
  for (i=0; i<nfacet; i++) {
    memsize = sizeof(FacetInfo);
    checkCudaErrors(cudaMemcpy(cudaInfo->h_d_facetInfo[i], h_facetInfo[i], memsize, cudaMemcpyHostToDevice));  
  }
} /* end UpdateGPUBase */

/**
 * Grids UV data using a GPU.
 * Multiple streams are used to overlap I/O and computation,
 * each call divides the data into streams_used pieces.
 * \param  streams_used  Number of streams to use
 * \param  cudaInfo      Information structure
 */
void grid_processWithStreams(int streams_used, CUDAGridInfo *cudaInfo, 
     cudaStream_t* stream, cudaEvent_t* cycleDone)
{
    GridInfo* gridInfo    = cudaInfo->h_gridInfo;       // grid specifications
    //FacetInfo** facetInfo = cudaInfo->h_facetInfo;     // array of facet specifications
    int nvis   = gridInfo->nvis;
    int nfacet = gridInfo->nfacet;
    int nchan  = gridInfo->nchan;
    int last_stream, current_stream = 0;
    int npass = streams_used;
    int nvisPass = (nvis+npass-1)/npass;  // nearly round up
    int lenvis = gridInfo->lenvis;
    int off, dovis,  nleft, i, ms, ifacet, nChBlock, nVisBlock;
    size_t lmemsize, memsize = (lenvis*nvisPass)*sizeof(float);
    dim3 numBlocks, thPerBlock;
    float *h_grid_debug=cudaInfo->h_grid_debug, *d_grid_debug=cudaInfo->d_grid_debug;
    // float *debug=NULL; // DEBUG
    gboolean initdebug=h_grid_debug==NULL;
    ms = 100*sizeof(float);
    // Initialize debug arrays 
    if (h_grid_debug==NULL) checkCudaErrors(cudaMallocHost(&h_grid_debug, ms));
    if (d_grid_debug==NULL) checkCudaErrors(cudaMalloc(&d_grid_debug, ms));
    //checkCudaErrors(cudaMallocHost(&debug, ms)); // DEBUG
    cudaInfo->h_grid_debug = h_grid_debug;
    cudaInfo->d_grid_debug = d_grid_debug;
    if (initdebug) {
      //memset (h_grid_debug, 0, ms);
      //checkCudaErrors(cudaMemcpy(d_grid_debug, h_grid_debug, ms, cudaMemcpyHostToDevice));
      // Update GPU structures
      UpdateGPUBase (cudaInfo);
    }

    // Do processing in a loop
    //
    // Note: All memory commands are processed in the order  they are issued,
    // independent of the stream they are enqueued in. Hence the pattern by
    // which the copy and kernel commands are enqueued in the stream
    // has an influence on the achieved overlap.

    int next_stream = (current_stream + 1) % streams_used;
    int prev_stream = current_stream - 1;
    if (prev_stream<0) prev_stream = streams_used-1;

   // Upload first frame
   memsize = (lenvis*nvisPass)*sizeof(float);
   checkCudaErrors(cudaMemcpyAsync(gridInfo->h_d_data_in[0],
			  	   gridInfo->h_data_in,
				   memsize,
				   cudaMemcpyHostToDevice,
				   gridInfo->stream[0]));
    nleft = nvis - nvisPass;  // How many left to copy?*/
    cudaEventSynchronize(cycleDone[0]);
    for (i=0; i<npass; ++i) {
        next_stream = (current_stream + 1) % streams_used;
	prev_stream = current_stream - 1;
	if (prev_stream<0) prev_stream = streams_used-1;
	off = next_stream*lenvis*nvisPass;  /* Offset in data buffers */
	if (nleft<nvisPass) lmemsize = (lenvis*nleft)*sizeof(float);
	else                lmemsize = memsize;
        // start copying next set of data
	if (nleft>0)
	  checkCudaErrors(cudaMemcpyAsync(gridInfo->h_d_data_in[next_stream],
	  		      &gridInfo->h_data_in[off],
                              lmemsize,
                              cudaMemcpyHostToDevice,
                              gridInfo->stream[next_stream]));

        // Process current
	// make sure to do all visibilities
	if (i==npass-1) dovis = nvis-i*nvisPass;
	else            dovis = nvisPass;

	checkCudaErrors(cudaEventSynchronize(cycleDone[current_stream]));
	if (dovis>0) { // More to do?
	  int maxCh = 16, maxVis=64;
          // Process current, vis (maxVis/block), chan(maxCh/block), loop over facet
 	  nVisBlock = (olong)(0.9999+dovis/(float)maxVis); 
	  nChBlock  = (olong)(0.9999+nchan/(float)maxCh); 
	  numBlocks.x = nVisBlock; numBlocks.y = nChBlock; thPerBlock.x = maxVis; thPerBlock.y = maxCh;
	  if (dovis<=maxVis) {numBlocks.x = 1; thPerBlock.x = dovis;}
	  if (nchan<=maxCh)  {numBlocks.y = 1; thPerBlock.y = nchan;}
	  numBlocks.x = MAX(1,numBlocks.x); numBlocks.y = MAX(1,numBlocks.y); 
	  thPerBlock.x = MAX(1,thPerBlock.x); thPerBlock.y = MAX(1,thPerBlock.y); 
//fprintf (stderr,"Grid X=%d %d, Y=%d %d,nVisBlock=%d, nChBlock=%d, dovis=%d, nchan=%d, nleft=%d\n", 
//	  numBlocks.x,thPerBlock.x, numBlocks.y,thPerBlock.y, nVisBlock,nChBlock,dovis,nchan,nleft); // DEBUG
          for (ifacet=0; ifacet<nfacet; ifacet++) {
	     gridKernel<<<numBlocks, thPerBlock, 0, gridInfo->stream[current_stream]>>>
                (current_stream, ifacet, dovis, (CUDAGridInfo *)cudaInfo->d_base, d_grid_debug);
          }
        } // end more data

	// make sure previous frame done
	if (i>0) checkCudaErrors(cudaEventSynchronize(cycleDone[prev_stream]));
	last_stream = current_stream;
        current_stream = next_stream;
	nleft -= nvisPass;
    } /* end loop */

    /* Last pass */
    checkCudaErrors(cudaEventSynchronize(cycleDone[last_stream]));
    checkCudaErrors(cudaDeviceSynchronize());
   
// Get debugging info 
//memsize = 1000*sizeof(float);
//checkCudaErrors(cudaMemcpy(h_grid_debug, d_grid_debug, memsize, cudaMemcpyDeviceToHost));
//fprintf(stderr,"Facet Counts %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g  \n", 
//  h_grid_debug[0],h_grid_debug[1],h_grid_debug[2],h_grid_debug[3], h_grid_debug[4],h_grid_debug[5],h_grid_debug[6],h_grid_debug[7]);

    return;
} // end grid_processWithStreams

/* Public: Allocate basic structures */
/**
 * Allocate most gridding structures, creates locked memory structures,
 * \param gridInfo gridding information
 * \param nfacet   Number of facets
 * \param nchan    Number of channels
 * \param nif      Number of IFs
 * \param nVisPIO  Number of vis records per transaction
 * \param lenvis   Length in floats of a vis record
 * \param nplane   Number of planes in image
 */
extern "C"
void ObitCUDAGridAlloc (CUDAGridInfo *cudaInfo, long nfacet, long nchan, long nif, 
     long nVisPIO, long lenvis, long nplane) {
  int i, ms;
  float *alloc;

  /* Debugging array pointers */
  cudaInfo->h_grid_debug = NULL;
  cudaInfo->d_grid_debug = NULL;

  // Number of Vis per IO
  cudaInfo->nVisPIO = nVisPIO;

  // Allocate arrays - General Griding 
  /* base address in GPU */

  ms = sizeof (GridInfo);
  checkCudaErrors(cudaMallocHost(&alloc, ms));
  cudaInfo->h_gridInfo = (GridInfo*)alloc;
  GridInfo* gridInfo    = cudaInfo->h_gridInfo;      // grid specifications
  gridInfo->nfacet = nfacet;
  gridInfo->nchan  = nchan*nif;
  gridInfo->nvis   = nVisPIO;
  gridInfo->lenvis = lenvis;
  gridInfo->nplane = nplane;
  ms = sizeof (GridInfo);
  checkCudaErrors(cudaMalloc(&alloc, ms));
  cudaInfo->d_gridInfo = (GridInfo*)alloc;

  ms = nVisPIO*lenvis*sizeof(float);
  checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_data_in, ms));
 
  ms = nchan*nif*sizeof(long);
  checkCudaErrors(cudaMallocHost(&alloc, ms));
  cudaInfo->h_gridInfo->h_freqPlane = (int*)alloc;
  checkCudaErrors(cudaMalloc(&alloc, ms));
  cudaInfo->h_gridInfo->d_freqPlane = (int*)alloc;
  memset (cudaInfo->h_gridInfo->h_freqPlane, 0, ms);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_gridInfo->d_freqPlane, 
      cudaInfo->h_gridInfo->h_freqPlane, ms, cudaMemcpyHostToDevice));

  ms = nchan*nif*sizeof(float);
  checkCudaErrors(cudaMallocHost(&alloc, ms));
  cudaInfo->h_gridInfo->h_sigma1 = (float*)alloc;
  checkCudaErrors(cudaMallocHost(&alloc, ms));
  cudaInfo->h_gridInfo->h_sigma2 = (float*)alloc;
  checkCudaErrors(cudaMallocHost(&alloc, ms));
  cudaInfo->h_gridInfo->h_sigma3 = (float*)alloc;
  checkCudaErrors(cudaMalloc(&alloc, ms));
  cudaInfo->h_gridInfo->d_sigma1 = (float*)alloc;
  checkCudaErrors(cudaMalloc(&alloc, ms));
  cudaInfo->h_gridInfo->d_sigma2 = (float*)alloc;
  checkCudaErrors(cudaMalloc(&alloc, ms));
  cudaInfo->h_gridInfo->d_sigma3 = (float*)alloc;
  memset (cudaInfo->h_gridInfo->h_sigma1, 0, ms);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_gridInfo->d_sigma1, 
      cudaInfo->h_gridInfo->h_sigma1, ms, cudaMemcpyHostToDevice));
  memset (cudaInfo->h_gridInfo->h_sigma2, 0, ms);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_gridInfo->d_sigma2, 
      cudaInfo->h_gridInfo->h_sigma2, ms, cudaMemcpyHostToDevice));
  memset (cudaInfo->h_gridInfo->h_sigma3, 0, ms);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_gridInfo->d_sigma3, 
      cudaInfo->h_gridInfo->h_sigma3, ms, cudaMemcpyHostToDevice));

  ms = nchan*nif*sizeof(float);
  checkCudaErrors(cudaMalloc    (&cudaInfo->h_gridInfo->d_freqArr, ms));
  checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_freqArr, ms));
  memset (cudaInfo->h_gridInfo->h_freqArr, 0, ms);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_gridInfo->d_freqArr, 
      cudaInfo->h_gridInfo->h_freqArr, ms, cudaMemcpyHostToDevice));

  ms = 2000*sizeof(float);
  checkCudaErrors(cudaMalloc    (&cudaInfo->h_gridInfo->d_convfn, ms));
  checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_convfn, ms));
  memset (cudaInfo->h_gridInfo->h_convfn, 0, ms);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_gridInfo->d_convfn, 
      cudaInfo->h_gridInfo->h_convfn, ms, cudaMemcpyHostToDevice));

  /* Streams */
  cudaInfo->h_gridInfo->nStream   = GRID_STREAM_COUNT;
  cudaInfo->h_gridInfo->stream    = (cudaStream_t*)g_malloc0(GRID_STREAM_COUNT*sizeof(cudaStream_t*));
  cudaInfo->h_gridInfo->cycleDone = (cudaEvent_t*) g_malloc0(GRID_STREAM_COUNT*sizeof(cudaEvent_t*)); 
  ms = (GRID_STREAM_COUNT*sizeof(float*));
  checkCudaErrors(cudaMalloc    (&alloc, ms));
  cudaInfo->h_gridInfo->d_data_in = (float**)alloc;
  // One in host locked memory to store pointers and copy later
  checkCudaErrors(cudaMallocHost (&alloc, ms));
  cudaInfo->h_gridInfo->h_d_data_in = (float**)alloc;
  memset (cudaInfo->h_gridInfo->h_d_data_in, 0, ms);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_gridInfo->d_data_in, 
      cudaInfo->h_gridInfo->h_d_data_in, ms, cudaMemcpyHostToDevice));

  for (i=0; i<GRID_STREAM_COUNT; i++) {
    checkCudaErrors(cudaStreamCreate(&cudaInfo->h_gridInfo->stream[i]));
    checkCudaErrors(cudaEventCreate(&cudaInfo->h_gridInfo->cycleDone[i]));
    ms = nVisPIO*lenvis*sizeof(float)/GRID_STREAM_COUNT;
    // device data buffers per stream, save pointer for now
    checkCudaErrors(cudaMalloc(&cudaInfo->h_gridInfo->h_d_data_in[i], ms));
  }

  /* Facet stuff */
  ms = nfacet*sizeof(FacetInfo*);
  checkCudaErrors(cudaMallocHost(&alloc, ms));
  cudaInfo->h_facetInfo = (FacetInfo**)alloc;
  FacetInfo** h_facetInfo = cudaInfo->h_facetInfo;     // array of facet specifications
  checkCudaErrors(cudaMallocHost(&alloc, ms));
  cudaInfo->h_d_facetInfo = (FacetInfo**)alloc;
  checkCudaErrors(cudaMalloc(&alloc, ms));
  cudaInfo->d_facetInfo = (FacetInfo**)alloc;
  FacetInfo** d_facetInfo = cudaInfo->h_d_facetInfo;
  for (i=0; i<nfacet; i++) {
    ms = sizeof(FacetInfo);
    checkCudaErrors(cudaMalloc    (&alloc, ms));
    d_facetInfo[i] = (FacetInfo*)alloc;
    checkCudaErrors(cudaMallocHost(&alloc, ms));
    h_facetInfo[i] = (FacetInfo*)alloc;
    memset (h_facetInfo[i], 0, ms);
    h_facetInfo[i]->rotUV[0]=1.0; h_facetInfo[i]->rotUV[4]=1.0; h_facetInfo[i]->rotUV[8]=1.0; 
    checkCudaErrors(cudaMemcpy(d_facetInfo[i], h_facetInfo[i], ms, cudaMemcpyHostToDevice));
  }
 
  // packet of info to pass to GPU
  CUDAGridInfo *d_base;
  ms = sizeof (CUDAGridInfo);
  checkCudaErrors(cudaMalloc(&d_base, ms));
  cudaInfo->d_base = d_base;

} /* end CUDAGridAlloc */

/**
 * Copy structures to GPU, allocate grid buffers
 * \param gridInfo gridding information
 */
extern "C"
void ObitCUDAGridInit (CUDAGridInfo *cudaInfo) {
  CUDAGridInfo *h_cudaInfo=NULL;
  GridInfo* h_gridInfo    = cudaInfo->h_gridInfo;      // host grid specifications
  FacetInfo** h_facetInfo = cudaInfo->h_facetInfo;     // host array of facet specifications
  GridInfo* d_gridInfo    = cudaInfo->d_gridInfo;      // device grid specifications
  FacetInfo** d_facetInfo = cudaInfo->d_facetInfo;     // device array of facet specifications
  int i, nfacet=h_gridInfo->nfacet, nplane=h_gridInfo->nplane,
    nchanIF=h_gridInfo->nchan*h_gridInfo->nif;
  size_t memsize;

  /* copy base structure to device */
  memsize = sizeof(CUDAGridInfo);
  checkCudaErrors(cudaMallocHost(&h_cudaInfo, memsize));
  memcpy (h_cudaInfo, cudaInfo, memsize);
  checkCudaErrors(cudaMemcpy(cudaInfo->d_base, h_cudaInfo, memsize, cudaMemcpyHostToDevice));  
  // DEBUG
  //memsize = 100;
  //checkCudaErrors(cudaMemcpy(cudaInfo->h_grid_debug, cudaInfo->d_base, memsize, cudaMemcpyDeviceToHost));  

  /* copy gridInfo to device */
  memsize = sizeof(GridInfo);
  checkCudaErrors(cudaMemcpy(d_gridInfo, h_gridInfo, memsize, cudaMemcpyHostToDevice));  

  /* arrays */
  memsize = nchanIF*sizeof(long);
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_freqPlane, h_gridInfo->h_freqPlane, memsize, cudaMemcpyHostToDevice));  
  memsize = nchanIF*sizeof(float);
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_sigma1,    h_gridInfo->h_sigma1,    memsize, cudaMemcpyHostToDevice));  
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_sigma2,    h_gridInfo->h_sigma2,    memsize, cudaMemcpyHostToDevice));  
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_sigma3,    h_gridInfo->h_sigma3,    memsize, cudaMemcpyHostToDevice));  
  memsize = nchanIF*sizeof(float);
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_freqArr, h_gridInfo->h_freqArr, memsize, cudaMemcpyHostToDevice));  
  memsize = 2000*sizeof(float);
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_convfn, h_gridInfo->h_convfn, memsize, cudaMemcpyHostToDevice));  
  /* Saved data buffer pointers */
  memsize = (GRID_STREAM_COUNT*sizeof(float*));
  checkCudaErrors(cudaMemcpy(h_gridInfo->d_data_in, h_gridInfo->h_d_data_in, memsize, cudaMemcpyHostToDevice));  

  /* Facet arrays 
      Allocate Grid buffers per facet */
  for (i=0; i<nfacet; i++) {
    memsize = h_facetInfo[i]->sizeGrid*nplane*sizeof(float);
    checkCudaErrors(cudaMallocHost(&h_facetInfo[i]->h_grid, memsize));
    checkCudaErrors(cudaMalloc    (&h_facetInfo[i]->d_grid, memsize));
    // initialize to zero
    memset (h_facetInfo[i]->h_grid, 0, memsize);
    checkCudaErrors(cudaMemcpy(h_facetInfo[i]->d_grid, h_facetInfo[i]->h_grid, memsize, cudaMemcpyHostToDevice));  
}
  /* copy to device */
  memsize = nfacet*sizeof(FacetInfo*);
  checkCudaErrors(cudaMemcpy(d_facetInfo, cudaInfo->h_d_facetInfo, memsize, cudaMemcpyHostToDevice));  
  for (i=0; i<nfacet; i++) {
    memsize = sizeof(FacetInfo);
    checkCudaErrors(cudaMemcpy(cudaInfo->h_d_facetInfo[i], h_facetInfo[i], memsize, cudaMemcpyHostToDevice));  
  }


} /* end ObitCUDAGridInit */

/* Public: Grid a bufferload of visibilities */
extern "C"
void ObitCUDAGridGrid (CUDAGridInfo *cudaInfo) {
  grid_processWithStreams(GRID_STREAM_COUNT, cudaInfo, 
                          cudaInfo->h_gridInfo->stream, cudaInfo->h_gridInfo->cycleDone);
} /* end ObitCUDAGridGrid   */

/* Public: Fold negative u to positive */
extern "C"
void ObitCUDAGridFlip (CUDAGridInfo *cudaInfo) {
   GridInfo*   h_gridInfo  = cudaInfo->h_gridInfo;
   FacetInfo** h_facetInfo = cudaInfo->h_facetInfo;
   int i, nrow, nrowBlock=32, nplane, nplaneBlock=16;
   dim3 numBlocks, thPerBlock;
   // Add conjugate rows

   // loop over facets, split by blocks of rows and planes
   nplane = h_gridInfo->nplane;
   nrow = h_facetInfo[0]->ny;

   for (i=0; i<h_gridInfo->nfacet; i++) {
     nrow = 1+h_facetInfo[i]->ny/2;
     numBlocks.x = (olong)(0.9999+nrow/((float)nrowBlock)); 
     numBlocks.y = (olong)(0.9999+nplane/((float)nplaneBlock)); 
     thPerBlock.x = nrowBlock;thPerBlock.y = nplaneBlock;
     if (nrow<=nrowBlock)     {numBlocks.x = 1; thPerBlock.x = nrow;}
     if (nplane<=nplaneBlock) {numBlocks.y = 1; thPerBlock.y = nplane;}
     numBlocks.x = MAX(1,numBlocks.x); numBlocks.y = MAX(1,numBlocks.y);
//fprintf (stderr,"Flip X=%d %d, Y=%d %d\n", numBlocks.x,thPerBlock.x, numBlocks.y,thPerBlock.y); // DEBUG
     gridFlipKernel<<<numBlocks,thPerBlock>>>
          ((CUDAGridInfo *)cudaInfo->d_base, i, cudaInfo->d_grid_debug);
        checkCudaErrors(cudaDeviceSynchronize());
    } // end facet loop

     cudaDeviceSynchronize();
} /* end ObitCUDAGridFlip */

/* Public: Copy  vis grid to locked host memory */
extern "C"
void ObitCUDAGrid2CPU (CUDAGridInfo *cudaInfo, long ifacet) {
   GridInfo* gridInfo    = cudaInfo->h_gridInfo;
   FacetInfo** facetInfo = cudaInfo->h_facetInfo;
   size_t memsize = facetInfo[ifacet]->sizeGrid*gridInfo->nplane*sizeof(float);
   checkCudaErrors(cudaMemcpy(facetInfo[ifacet]->h_grid, facetInfo[ifacet]->d_grid, 
                   memsize, cudaMemcpyDeviceToHost));
} /* end ObitCUDAGrid2CPU */

/* Public: Shutdown gridding */
extern "C"
void ObitCUDAGridShutdown (CUDAGridInfo *cudaInfo) {
  GridInfo* h_gridInfo    = cudaInfo->h_gridInfo;
  //GridInfo* d_gridInfo    = cudaInfo->d_gridInfo;
  FacetInfo** h_facetInfo = cudaInfo->h_facetInfo;
  int i;
  // Free allocated memory
  if (cudaInfo->h_grid_debug) cudaFreeHost(cudaInfo->h_grid_debug);
  cudaInfo->h_grid_debug = NULL;
  if (cudaInfo->d_grid_debug) cudaFree(cudaInfo->d_grid_debug);
  cudaInfo->d_grid_debug = NULL;

  // streams
  for (i=0; i<GRID_STREAM_COUNT; i++) {
    //if ((d_gridInfo->h_d_data_in) && (d_gridInfo->h_d_data_in[i]))
    //  {cudaFree(d_gridInfo->h_d_data_in[i]); d_gridInfo->h_d_data_in[i] = NULL;}
    if ((h_gridInfo->h_d_data_in) && (h_gridInfo->h_d_data_in[i]))
      {cudaFree(h_gridInfo->h_d_data_in[i]); h_gridInfo->h_d_data_in[i] = NULL;}
    cudaStreamDestroy((cudaStream_t)h_gridInfo->stream[i]);
    cudaEventDestroy(h_gridInfo->cycleDone[i]);
  }

  // data buffers
  //cudaFree    (d_gridInfo->h_d_data_in); d_gridInfo->d_data_in   = NULL;
  cudaFreeHost(h_gridInfo->h_data_in);   h_gridInfo->h_data_in   = NULL;

  // Others
  if (h_gridInfo->h_freqArr)   {cudaFreeHost(h_gridInfo->h_freqArr);    h_gridInfo->h_freqArr = NULL;}
  if (h_gridInfo->h_freqPlane) {cudaFreeHost(h_gridInfo->h_freqPlane);  h_gridInfo->h_freqPlane = NULL;}
  if (h_gridInfo->h_sigma1)    {cudaFreeHost(h_gridInfo->h_sigma1);     h_gridInfo->h_sigma1 = NULL;}
  if (h_gridInfo->h_sigma2)    {cudaFreeHost(h_gridInfo->h_sigma2);     h_gridInfo->h_sigma2 = NULL;}
  if (h_gridInfo->h_sigma3)    {cudaFreeHost(h_gridInfo->h_sigma3);     h_gridInfo->h_sigma3 = NULL;}
  if (h_gridInfo->h_convfn)    {cudaFreeHost(h_gridInfo->h_convfn);     h_gridInfo->h_convfn = NULL;}
  if (h_gridInfo->d_freqArr)   {cudaFree    (h_gridInfo->d_freqArr);    h_gridInfo->d_freqArr = NULL;}
  if (h_gridInfo->d_freqPlane) {cudaFree    (h_gridInfo->d_freqPlane);  h_gridInfo->d_freqPlane = NULL;}
  if (h_gridInfo->d_sigma1)    {cudaFree    (h_gridInfo->h_sigma1);     h_gridInfo->d_sigma1 = NULL;}
  if (h_gridInfo->d_sigma2)    {cudaFree    (h_gridInfo->h_sigma2);     h_gridInfo->d_sigma2 = NULL;}
  if (h_gridInfo->d_sigma3)    {cudaFree    (h_gridInfo->h_sigma3);     h_gridInfo->d_sigma3 = NULL;}
  if (h_gridInfo->d_convfn)    {cudaFree    (h_gridInfo->d_convfn);     h_gridInfo->d_convfn = NULL;}
  // Facet stuff 
  for (i=0; i<h_gridInfo->nfacet; i++) {
     if (h_facetInfo[i]->d_grid) {cudaFree    (h_facetInfo[i]->d_grid); h_facetInfo[i]->d_grid = NULL;}
     if (h_facetInfo[i]->h_grid) {cudaFreeHost(h_facetInfo[i]->h_grid); h_facetInfo[i]->h_grid = NULL;}
     if (h_facetInfo[i]) {cudaFreeHost(h_facetInfo[i]); h_facetInfo[i]=NULL;}
   }
   cudaFree    (cudaInfo->d_gridInfo);  cudaInfo->d_gridInfo  = NULL;
   cudaFreeHost(cudaInfo->h_gridInfo);  cudaInfo->h_gridInfo  = NULL;
   cudaFreeHost(cudaInfo->h_facetInfo); cudaInfo->h_facetInfo = NULL;

   // Others
   if (cudaInfo->d_base)      {cudaFree    (cudaInfo->d_base);      cudaInfo->d_base = NULL;}
   if (cudaInfo->d_facetInfo) {cudaFree    (cudaInfo->d_facetInfo); cudaInfo->d_facetInfo = NULL;}

   // Reset GPU
   cudaDeviceReset();

  } /* end ObitCUDAGridShutdown */
#endif // CUDA code

