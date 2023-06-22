/* $Id: $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2021,2023                                          */
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

#include "ObitCUDAGrid.h"
#include "ObitCUDAUtil.h"
/* Public: pointer info type DAMN*/
void ObitCUDAUtilPointerType(void *ptr, char *label);

#include "ObitCUDAGrid.h"
#include "ObitCUDAGridInfoDef.h"
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
//#include "ObitCUDAGrid.h"
 
// functions
static void grid_processWithStreams(int streams_used, CUDAGridInfo* gridInfo, int gpu,
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
__global__ void gridKernel(int current, int off, int nvis, CUDAGridInfo* cudaInfo, float *d_grid_debug)
{
    long kvis         = threadIdx.x + blockIdx.x*blockDim.x;  // vis number
    long ichan        = threadIdx.y + blockIdx.y*blockDim.y;  // channel number
    long ifacet       = threadIdx.z + blockIdx.z*blockDim.z;  // facet number
    GridInfo* gridInfo = cudaInfo->d_gridInfo;                // grid specifications
    FacetInfo* facetInfo = cudaInfo->d_facetInfo[ifacet];     // facet specifications
    long iplane       = gridInfo->d_freqPlane[ichan];         // plane number
    float *grid       = facetInfo->d_grid + iplane*facetInfo->sizeGrid;
    float *freqArr    = NULL; //gridInfo->d_freqArr;
    float *vis_in     = &gridInfo->d_data_in[off];
    long ivis         = kvis * gridInfo->lenvis;              // beginning of visibility
    int halfWidth     = gridInfo->convWidth/2;                // half width of convolution kernal
    int fullWidth     = gridInfo->convWidth;                  // full width of convolution kernal
    int convNperCell  = gridInfo->convNperCell;               // resolution of kernal
    int  nxo2         = facetInfo->nx/2;
    int  nyo2         = facetInfo->ny/2;
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
    if (kvis>=nvis) return;
    if (ifacet>=gridInfo->nfacet) return;
    if ((ichan>=gridInfo->nchan) || (iplane>=gridInfo->nplane)) return;

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
//atomicAddFloat(&d_grid_debug[ifacet], sigma1); 
//end DEBUG

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
__global__ void gridFlipKernel(int gpu, CUDAGridInfo *cudaInfo, int ifacet, float *d_grid_debug)
{
   long irow          = threadIdx.x + blockIdx.x*blockDim.x;  // row number
   long iplane        = threadIdx.y + blockIdx.y*blockDim.y;  // plane
   GridInfo* gridInfo = cudaInfo->d_gridInfo;                 // grid specifications
   long halfWidth     = gridInfo->convWidth/2;    // half width of convolution kernal

   FacetInfo* facetInfo = cudaInfo->d_facetInfo[ifacet];      // facet specifications
   int  nx            = facetInfo->nx;
   int  ny            = facetInfo->ny;
   float *grid        = facetInfo->d_grid + iplane*facetInfo->sizeGrid;
   float *gxi, *gxo, *gci, *gco, xxo[2], cjo[2], xxi[2], cji[2];
   long iu, vc, lrow;
   // Check out of range 
   if ((irow>ny/2) || (iplane>=gridInfo->nplane)) return;

   vc = ny - irow;               // conjugate row number
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
/* Public: Allocate basic structures */
/**
 * Allocate most gridding structures, creates locked host memory structures,
 * \param gridInfo gridding information
 * \param nfacet   Number of facets
 * \param nchan    Number of channels
 * \param nif      Number of IFs
 * \param nVisPIO  Number of vis records per transaction
 * \param lenvis   Length in floats of a vis record
 * \param nplane   Number of planes in image
 * \param iGPU     Which GPU index in use
 * \param doBuff   If TRUE allocate buffer
 */
extern "C"
void ObitCUDAGridAlloc (CUDAGridInfo *cudaInfo, long nfacet, long nchan, long nif, 
     long nVisPIO, long lenvis, long nplane, long iGPU, int doBuff) {
  int i, ms;
  long nGPU = (long)cudaInfo->nGPU;
  FacetInfo **t_facetArr=NULL;     // Temporary host version of device specific facetInfo Structures
  CUDAGridInfo *t_cudaInfo=NULL;   // Host version of device specific cudaInfo Structure

  /* Debugging array pointers */
  cudaInfo->h_grid_debug = NULL;
  cudaInfo->d_grid_debug = NULL;

  // copies of cudaInfo to pass to GPUs
  //  Device
  ms = nGPU * sizeof (CUDAGridInfo**);
  if (!cudaInfo->d_base) {  // device version
    checkCudaErrors(cudaMallocHost(&cudaInfo->d_base, nGPU*sizeof(CUDAGridInfo**))); 
    for (i=0; i<nGPU; i++) cudaInfo->d_base[i] = NULL;
  }
  if (!cudaInfo->d_base[iGPU]) {  // copy for each gpu - is there already one?
    checkCudaErrors(cudaMalloc(&cudaInfo->d_base[iGPU], sizeof (CUDAGridInfo)));
  }

  //  Host
  if (!cudaInfo->h_base) {  // host version
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_base, nGPU*sizeof(CUDAGridInfo**)));
    for (i=0; i<nGPU; i++) cudaInfo->h_base[i] = NULL;
  }
  ms = sizeof (CUDAGridInfo);
  if (!cudaInfo->h_base[iGPU]) { // extant version?
    checkCudaErrors(cudaMallocHost(&t_cudaInfo, sizeof(CUDAGridInfo))); // for this device
    cudaInfo->h_base[iGPU] = t_cudaInfo;
    // copy current host version of cudaInfo
    memcpy(t_cudaInfo, cudaInfo, sizeof (CUDAGridInfo));
  } else {  // Use old one
    t_cudaInfo = (CUDAGridInfo*)cudaInfo->h_base[iGPU];
  }

  cudaInfo->nVisPIO = nVisPIO;  // Number of Vis per IO
  cudaInfo->nfacet  = nfacet;   // Number of facets

  // *******************gridInfo*********************************************
  // Host array to keep device pointers
  if (!cudaInfo->h_d_gridInfo) {
    ms = nGPU* sizeof(GridInfo*);
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_d_gridInfo, ms));
    memset (cudaInfo->h_d_gridInfo, 0, ms); // zero
  }  // end create h_d_gridInfo

  // Allocate arrays - General Griding 
  if (!cudaInfo->h_gridInfo) {
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo, sizeof (GridInfo)));
    memset (cudaInfo->h_gridInfo, 0, sizeof(GridInfo)); // zero
    cudaInfo->h_gridInfo->nStream = GRID_STREAM_COUNT;
    // stream, events arrays per GPU 
    cudaInfo->h_gridInfo->stream    = (cudaStream_t**)malloc(nGPU*sizeof(cudaStream_t*));
    cudaInfo->h_gridInfo->cycleDone = (cudaEvent_t**)malloc(nGPU*sizeof(cudaEvent_t*));
  }
  GridInfo* gridInfo    = cudaInfo->h_gridInfo;      // grid specifications
  gridInfo->nfacet = nfacet;
  gridInfo->nchan  = nchan*nif;
  gridInfo->nif    = nif;
  gridInfo->nvis   = nVisPIO;
  gridInfo->lenvis = lenvis;
  gridInfo->nplane = nplane;
  gridInfo->nGPU   = nGPU;
  // device version for this GPU
  if (!t_cudaInfo->d_gridInfo) {
    checkCudaErrors(cudaMalloc(&t_cudaInfo->d_gridInfo,sizeof (GridInfo)));
    // Copy host version to where it's going in device
    checkCudaErrors(cudaMemcpy(t_cudaInfo->d_gridInfo, gridInfo, sizeof(GridInfo), cudaMemcpyHostToDevice));
   }

  GridInfo* t_gridInfo;              // temp host copy of device grid specifications
  checkCudaErrors(cudaMallocHost(&t_gridInfo, sizeof(GridInfo)));
  // Init with host version
  memcpy(t_gridInfo, gridInfo, sizeof(GridInfo));

  // Host vis data buffer
  if (doBuff && !cudaInfo->h_gridInfo->h_data_in) {
    ms = nVisPIO*lenvis*sizeof(float);
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_data_in, ms));
  }
  // save pointer in device specific version
  t_gridInfo->h_data_in = cudaInfo->h_gridInfo->h_data_in;
 
  // Host storage for device buffer per GPU
  if (!cudaInfo->h_gridInfo->h_d_data_in) 
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_d_data_in, nGPU*sizeof(float*)));
 

  // This device vis data_in buffer
  if (doBuff) {
    checkCudaErrors(cudaMalloc(&t_gridInfo->d_data_in, nVisPIO*lenvis*sizeof(float)));
    cudaInfo->h_gridInfo->h_d_data_in[iGPU] = t_gridInfo->d_data_in; // Save device pointer in host copy
    t_cudaInfo->d_data_in                   = t_gridInfo->d_data_in; // Save device pointer in host copy
  }
 
  // Frequency planes
  // host
  if (!cudaInfo->h_gridInfo->h_freqPlane) {
    ms = nchan*nif*sizeof(long);
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_freqPlane, ms));
    memset (cudaInfo->h_gridInfo->h_freqPlane, 0, ms);  // zero
  }
  // this device
  ms = nchan*nif*sizeof(long);
  checkCudaErrors(cudaMalloc(&t_gridInfo->d_freqPlane, ms));
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_freqPlane, 
        cudaInfo->h_gridInfo->h_freqPlane, ms, cudaMemcpyHostToDevice));

  //  grid info arrays
  // Host
  ms = nchan*nif*sizeof(float);
  if (!cudaInfo->h_gridInfo->h_sigma1) {
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_sigma1, ms));
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_sigma2, ms));
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_sigma3, ms));
    //  zero
    memset (cudaInfo->h_gridInfo->h_sigma1, 0, ms);
    memset (cudaInfo->h_gridInfo->h_sigma2, 0, ms);
    memset (cudaInfo->h_gridInfo->h_sigma3, 0, ms);
  }
  // This device
  checkCudaErrors(cudaMalloc(&t_gridInfo->d_sigma1, ms));
  checkCudaErrors(cudaMalloc(&t_gridInfo->d_sigma2, ms));
  checkCudaErrors(cudaMalloc(&t_gridInfo->d_sigma3, ms));
  //  zero
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_sigma1, 
      cudaInfo->h_gridInfo->h_sigma1, ms, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_sigma2, 
    cudaInfo->h_gridInfo->h_sigma2, ms, cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_sigma3, 
    cudaInfo->h_gridInfo->h_sigma3, ms, cudaMemcpyHostToDevice));

  // frequency arrays
  // Host
  if (!cudaInfo->h_gridInfo->h_freqArr) {
    ms = nchan*nif*sizeof(float);
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_freqArr, ms));
    // zero
    memset (cudaInfo->h_gridInfo->h_freqArr, 0, ms);
  }
  // This device
  ms = nchan*nif*sizeof(float);
  checkCudaErrors(cudaMalloc (&t_gridInfo->d_freqArr, ms));
  // zero
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_freqArr, 
        cudaInfo->h_gridInfo->h_freqArr, ms, cudaMemcpyHostToDevice));

  // convolution functions
  // Host
  if (!cudaInfo->h_gridInfo->d_convfn) {
  ms = MAX_CM_CONVFN*sizeof(float);  // should always be big enough
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_gridInfo->h_convfn, ms));
    memset (cudaInfo->h_gridInfo->h_convfn, 0, ms); // zero
  }
  // This device
  ms = MAX_CM_CONVFN*sizeof(float);  // should always be big enough
  checkCudaErrors(cudaMalloc (&t_gridInfo->d_convfn, ms));
  // zero
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_convfn, 
      cudaInfo->h_gridInfo->h_convfn, ms, cudaMemcpyHostToDevice));

  // Streams - only on host
  cudaInfo->h_gridInfo->nStream   = GRID_STREAM_COUNT;
  cudaInfo->h_gridInfo->stream[iGPU]    = (cudaStream_t*)malloc(GRID_STREAM_COUNT*sizeof(cudaStream_t));
  cudaInfo->h_gridInfo->cycleDone[iGPU] = (cudaEvent_t*) malloc(GRID_STREAM_COUNT*sizeof(cudaEvent_t)); 

  // create streams, events this device
  for (i=0; i<GRID_STREAM_COUNT; i++) {
    checkCudaErrors(cudaStreamCreate(&cudaInfo->h_gridInfo->stream[iGPU][i]));
    checkCudaErrors(cudaEventCreate(&cudaInfo->h_gridInfo->cycleDone[iGPU][i]));
  }

  // Finished gridInfo; copy t_GridInfo to device
  checkCudaErrors(cudaMemcpy(t_cudaInfo->d_gridInfo, t_gridInfo, sizeof(GridInfo), cudaMemcpyHostToDevice));
  cudaInfo->h_d_gridInfo[iGPU] = t_cudaInfo->d_gridInfo;  // save pointer in host memory
  
  if (t_gridInfo) cudaFreeHost(t_gridInfo);  // Free temp array

  // *******************facetInfo*********************************************
  // Host array to keep device pointers
  if (!cudaInfo->h_d_facetInfo) {
    ms = nGPU* sizeof(FacetInfo**);
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_d_facetInfo, ms));
    memset (cudaInfo->h_d_facetInfo, 0, ms); // zero
  }  // end create h_d_facetInfo
  // This device
  checkCudaErrors(cudaMallocHost(&cudaInfo->h_d_facetInfo[iGPU], nfacet*sizeof(FacetInfo*)));
  

  // Facet stuff - host version 
  if (!cudaInfo->h_facetInfo) {
    ms = nfacet*sizeof(FacetInfo*);
    checkCudaErrors(cudaMallocHost(&cudaInfo->h_facetInfo, ms));
    for (i=0; i<nfacet; i++) {
      checkCudaErrors(cudaMallocHost(&cudaInfo->h_facetInfo[i], sizeof(FacetInfo)));
      memset (cudaInfo->h_facetInfo[i], 0, sizeof(FacetInfo)); // zero
      cudaInfo->h_facetInfo[i]->GPU_num = iGPU;
      cudaInfo->h_facetInfo[i]->rotUV[0]=1.0; 
      cudaInfo->h_facetInfo[i]->rotUV[4]=1.0; 
      cudaInfo->h_facetInfo[i]->rotUV[8]=1.0; 
    }
  }

  // This device
  ms = nfacet*sizeof(FacetInfo*);
  checkCudaErrors(cudaMalloc(&t_cudaInfo->d_facetInfo, ms));
  // Allocate entries in array using t_facetArr,  temp. host array for d_facetInfo array
  checkCudaErrors(cudaMallocHost(&t_facetArr, ms));
  // Facet entries in device
  for (i=0; i<nfacet; i++) {
    ms = sizeof(FacetInfo);  
    checkCudaErrors(cudaMalloc(&t_facetArr[i], ms));
    // copy host version (zeroed)
    checkCudaErrors(cudaMemcpy(t_facetArr[i], cudaInfo->h_facetInfo[i], ms, cudaMemcpyHostToDevice));
  }
  // copy array of pointers to device
  ms = nfacet*sizeof(FacetInfo*);
  checkCudaErrors(cudaMemcpy(t_cudaInfo->d_facetInfo, t_facetArr, ms, cudaMemcpyHostToDevice));
  //??cudaInfo->h_d_facetInfo[iGPU] = t_cudaInfo->d_facetInfo;  // save pointer in host

  // Copy working host version of device cudaInfo to device
  ms  = sizeof (CUDAGridInfo);
  checkCudaErrors(cudaMemcpy(cudaInfo->d_base[iGPU], t_cudaInfo, ms, cudaMemcpyHostToDevice));
  memcpy(cudaInfo->h_base[iGPU], t_cudaInfo, ms);  // also to host copy
  if (t_facetArr) cudaFreeHost(t_facetArr);  // Free temp array
} // end CUDAGridAlloc

/**
 * Copy structures to GPU, allocate grid buffers, fill in info
 * \param gridInfo gridding information
 */
extern "C"
void ObitCUDAGridInit (CUDAGridInfo *cudaInfo, long iGPU) {
  CUDAGridInfo *t_cudaInfo=NULL;
  GridInfo* h_gridInfo = cudaInfo->h_gridInfo;       // host grid specifications
  FacetInfo** h_facetInfo;                           // host array of facet specifications
  FacetInfo **t_facetArr;                            // Temp host version of device array
  FacetInfo *t_facetInfo;                            // Temp host version of a device facetInfo
  int i, nfacet=h_gridInfo->nfacet, nplane=h_gridInfo->nplane,
    nchanIF=h_gridInfo->nchan*h_gridInfo->nif;
  size_t ms;

  // Get base structure (cudaInfo) for device 
  t_cudaInfo = (CUDAGridInfo*)cudaInfo->h_base[iGPU];

  GridInfo* t_gridInfo;              // temp host copy of device grid specifications
  checkCudaErrors(cudaMallocHost(&t_gridInfo, sizeof(GridInfo)));
  // get host copy of device gridInfo
  checkCudaErrors(cudaMemcpy(t_gridInfo, t_cudaInfo->d_gridInfo, sizeof(GridInfo), cudaMemcpyDeviceToHost));

  // copy info from host version
  t_gridInfo->nfacetPerGPU = h_gridInfo-> nfacetPerGPU;
  t_gridInfo->nplane       = h_gridInfo->nplane;
  t_gridInfo->nchan        = h_gridInfo->nchan;
  t_gridInfo->nif          = h_gridInfo->nif;
  t_gridInfo->nstok        = h_gridInfo->nstok;
  t_gridInfo->nrparm       = h_gridInfo->nrparm;
  t_gridInfo->nvis         = h_gridInfo->nvis;
  t_gridInfo->convWidth    = h_gridInfo->convWidth;
  t_gridInfo->convNperCell = h_gridInfo->convNperCell;
  t_gridInfo->rotate       = h_gridInfo->rotate;
  t_gridInfo->guardu       = h_gridInfo->guardu;
  t_gridInfo->guardv       = h_gridInfo->guardv;
  t_gridInfo->oldnVisPIO   = h_gridInfo->oldnVisPIO;
  t_gridInfo->nStream      = h_gridInfo->nStream;
  t_gridInfo->prtLv        = h_gridInfo->prtLv;
 
  /* arrays */
  ms = nchanIF*sizeof(long);
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_freqPlane, h_gridInfo->h_freqPlane, ms, cudaMemcpyHostToDevice));  
  ms = nchanIF*sizeof(float);
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_sigma1,    h_gridInfo->h_sigma1,    ms, cudaMemcpyHostToDevice));  
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_sigma2,    h_gridInfo->h_sigma2,    ms, cudaMemcpyHostToDevice));  
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_sigma3,    h_gridInfo->h_sigma3,    ms, cudaMemcpyHostToDevice));  
  ms = nchanIF*sizeof(float);
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_freqArr, h_gridInfo->h_freqArr, ms, cudaMemcpyHostToDevice));  
  ms = MAX_CM_CONVFN*sizeof(float);
  checkCudaErrors(cudaMemcpy(t_gridInfo->d_convfn, h_gridInfo->h_convfn, ms, cudaMemcpyHostToDevice));  
  // Put it backt to device
  checkCudaErrors(cudaMemcpy(t_cudaInfo->d_gridInfo, t_gridInfo, sizeof(GridInfo), cudaMemcpyHostToDevice));
  if (t_gridInfo) cudaFreeHost(t_gridInfo);  // Free temp array

  /* Facet arrays - Allocate Grid buffers per facet */
  h_facetInfo = cudaInfo->h_facetInfo;     // host array of facet specifications
  for (i=0; i<nfacet; i++) {
    if (iGPU==cudaInfo->FacetGPU[i]) {  // This one in this GPU?
      ms = h_facetInfo[i]->sizeGrid*nplane*sizeof(float);
      checkCudaErrors(cudaMallocHost(&h_facetInfo[i]->h_grid, ms));
      checkCudaErrors(cudaMalloc    (&h_facetInfo[i]->d_grid, ms));
      // initialize to zero
      memset (h_facetInfo[i]->h_grid, 0, ms);
      checkCudaErrors(cudaMemcpy(h_facetInfo[i]->d_grid, h_facetInfo[i]->h_grid, ms, cudaMemcpyHostToDevice));
      }
  }
  // copy all to device - fetch pointer array from device
  ms = nfacet*sizeof(FacetInfo*);
  checkCudaErrors(cudaMallocHost(&t_facetArr, ms));
  checkCudaErrors(cudaMallocHost(&t_facetInfo, sizeof(FacetInfo)));
  checkCudaErrors(cudaMemcpy(t_facetArr, t_cudaInfo->d_facetInfo, ms, cudaMemcpyDeviceToHost));
  // update with pointer
  ms = sizeof(FacetInfo);
  for (i=0; i<nfacet; i++) {
    if (h_facetInfo[i]->d_grid) {  // need something
      checkCudaErrors(cudaMemcpy(t_facetInfo, t_facetArr[i], ms, cudaMemcpyDeviceToHost));
      t_facetInfo->d_grid = h_facetInfo[i]->d_grid;
      checkCudaErrors(cudaMemcpy(t_facetArr[i], t_facetInfo, ms, cudaMemcpyHostToDevice));
      cudaInfo->h_d_facetInfo[iGPU][i] = t_facetArr[i];  // save pointer in host
    }
  }

  cudaFreeHost(t_facetArr); cudaFreeHost(t_facetInfo);  // Free temp arrays
  //* no cudaInfo->d_facetInfo[iGPU] = t_facetArr;
  //NO if (t_cudaInfo) cudaFreeHost(t_cudaInfo);  // Free temp device cudaInfo array

} /* end ObitCUDAGridInit */

/**
 * Copy cuda base structures to GPU
 * \param gpu      GPU index
 * \param cudaInfo GPU gridding information
 */
extern "C"
void UpdateGPUBase (long gpu, CUDAGridInfo *cudaInfo) {
  CUDAGridInfo *t_cudaInfo=NULL;
  GridInfo* h_gridInfo    = cudaInfo->h_gridInfo;      // host grid specifications
  FacetInfo** h_facetInfo;                             // host array of facet specifications
  FacetInfo** d_facetInfo;                             // device array of facet specifications
  FacetInfo **t_facetArr;                              // Temp host version of device array
  int nfacet=h_gridInfo->nfacet;
  size_t ms;
  long i, iGPU, lastGPU=-1;

  // Load constant memory
  ms = ((h_gridInfo->convWidth+1)*h_gridInfo->convNperCell)*sizeof(float);
  checkCudaErrors(cudaMemcpyToSymbol(cm_convfn,  h_gridInfo->h_convfn,  ms));
  ms = MIN(MAX_CM_FREQARR,h_gridInfo->nchan)*sizeof(float);
  checkCudaErrors(cudaMemcpyToSymbol(cm_freqArr, h_gridInfo->h_freqArr, ms));

  // Get base structure (cudaInfo) for device 
  t_cudaInfo = (CUDAGridInfo*)cudaInfo->h_base[gpu];
  // copy to GPU
  checkCudaErrors(cudaMemcpy(cudaInfo->d_base[gpu], t_cudaInfo, sizeof(CUDAGridInfo), cudaMemcpyHostToDevice));

  GridInfo* t_gridInfo;              // temp host copy of device grid specifications
  checkCudaErrors(cudaMallocHost(&t_gridInfo, sizeof(GridInfo)));
  // get host copy of device gridInfo
  checkCudaErrors(cudaMemcpy(t_gridInfo, t_cudaInfo->d_gridInfo, sizeof(GridInfo), cudaMemcpyDeviceToHost));

  // copy to device - fetch faceInfo pointer array from device
  checkCudaErrors(cudaMallocHost(&t_facetArr, nfacet*sizeof(FacetInfo**)));
  h_facetInfo = cudaInfo->h_facetInfo; // host array of facet specifications
  for (i=0; i<nfacet; i++) {
    iGPU = cudaInfo->FacetGPU[i]; // which gpu index am I
    if (iGPU!=gpu) continue;  // This GPU?
    if (iGPU!=lastGPU) {
      // copy facetInfo* array from h_d_facetInfo[iGPU] to device
      ms = nfacet*sizeof(FacetInfo*);
      d_facetInfo = cudaInfo->h_d_facetInfo[iGPU]; // device array of facet specifications
      checkCudaErrors(cudaMemcpy(t_facetArr, d_facetInfo, ms, cudaMemcpyDeviceToHost));
    }
    lastGPU = iGPU;
    ms = sizeof(FacetInfo);
    checkCudaErrors(cudaMemcpy(t_facetArr[i], h_facetInfo[i], ms, cudaMemcpyHostToDevice));  
  } // end facet loop

  // save pointer array where host can see it
  ms = nfacet*sizeof(FacetInfo**);
  checkCudaErrors(cudaMemcpy(cudaInfo->h_d_facetInfo[gpu], t_facetArr, ms, cudaMemcpyHostToDevice));
  cudaFreeHost(t_facetArr);  // Free temp array
} // end UpdateGPUBase 

/**
 * Grids UV data using a GPU.
 * Multiple streams are used to overlap I/O and computation,
 * each call divides the data into streams_used pieces.
 * \param  streams_used  Number of streams to use
 * \param  cudaInfo      Information structure
 * \param  gpu           GPU index
 * \param  stream        Streams for processing
 * \param  cycleDone     Events for processing
 */
void grid_processWithStreams(int streams_used, CUDAGridInfo *cudaInfo, 
     int gpu, cudaStream_t* stream, cudaEvent_t* cycleDone)
{
    GridInfo* gridInfo    = cudaInfo->h_gridInfo;       // grid specifications
    //FacetInfo** facetInfo = cudaInfo->h_facetInfo;     // array of facet specifications
    int nvis   = gridInfo->nvis;
    int nfacet = gridInfo->nfacet;
    int last_stream, current_stream = 0;
    int npass = streams_used;
    int nvisPass = (nvis+npass-1)/npass;  // nearly round up
    int lenvis = gridInfo->lenvis;
    int off, dovis,  nleft, i, ms, nVisBlock, nFacetBlock;
    size_t lmemsize, memsize = (lenvis*nvisPass)*sizeof(float);
    dim3 numBlocks, thPerBlock;
    float *h_grid_debug=cudaInfo->h_grid_debug, *d_grid_debug=cudaInfo->d_grid_debug;
    // float *debug=NULL; // DEBUG
    int initdebug=h_grid_debug==NULL;
    ms = 1000*sizeof(float);
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
   checkCudaErrors(cudaMemcpyAsync(gridInfo->h_d_data_in[gpu],
			  	   gridInfo->h_data_in,
				   memsize,
				   cudaMemcpyHostToDevice,
				   stream[0]));
    nleft = nvis - nvisPass;  // How many left to copy?*/
    cudaEventSynchronize(cycleDone[0]);
    for (i=0; i<npass; ++i) {
        next_stream = (current_stream + 1) % streams_used;
	prev_stream = current_stream - 1;
	if (prev_stream<0) prev_stream = streams_used-1;
	off = next_stream*lenvis*nvisPass;  // Offset in data buffers
	if (nleft<nvisPass) lmemsize = (lenvis*nleft)*sizeof(float);
	else                lmemsize = memsize;
        // start copying next set of data
	if (nleft>0)
	  checkCudaErrors(cudaMemcpyAsync(&gridInfo->h_d_data_in[gpu][off],
	  		      &gridInfo->h_data_in[off],
                              lmemsize,
                              cudaMemcpyHostToDevice,
                              stream[next_stream]));

        // Process current
	// make sure to do all visibilities
	if (i==npass-1) dovis = nvis-i*nvisPass;
	else            dovis = nvisPass;

	//checkCudaErrors(cudaEventSynchronize(cycleDone[current_stream]));
	if (dovis>0) { // More to do?
	  // divide work depending on number of facets.
	  int nChBlock, nchan,  maxCh, maxVis, maxFacet;
	  nchan = gridInfo->nchan;
	  if (nfacet<8) {maxCh = 32; maxVis=32; maxFacet=1;}
	  else          {maxCh = 16; maxVis=8;  maxFacet=8;}
          // Process current, vis (maxVis/block), chan(maxCh/block), loop over facet
 	  nVisBlock = (long)(0.9999+dovis/(float)maxVis); 
	  nChBlock  = (long)(0.9999+nchan/(float)maxCh); 
	  nFacetBlock =  (long)(0.9999+nfacet/(float)maxFacet);
	  numBlocks.x = nVisBlock; numBlocks.y = nChBlock; thPerBlock.x = maxVis; thPerBlock.y = maxCh;
	  numBlocks.z = nFacetBlock; thPerBlock.z = maxFacet; 
	  if (dovis<=maxVis)  {numBlocks.x = 1; thPerBlock.x = dovis;}
	  if (nchan<=maxCh)   {numBlocks.y = 1; thPerBlock.y = nchan;}
	  if (nfacet<=maxFacet)  {numBlocks.z = 1; thPerBlock.z = nfacet;}
	  numBlocks.x =  MAX(1,numBlocks.x);  numBlocks.y = MAX(1,numBlocks.y); 
	  thPerBlock.x = MAX(1,thPerBlock.x); thPerBlock.y = MAX(1,thPerBlock.y); 
	  numBlocks.z  = MAX(1,numBlocks.z);  thPerBlock.z = MAX(1,thPerBlock.z); 
//fprintf (stderr,"Grid X=%d %d, Y=%d %d,nVisBlock=%d, nChBlock=%d, dovis=%d, nchan=%d, nleft=%d\n", 
//	  numBlocks.x,thPerBlock.x, numBlocks.y,thPerBlock.y, nVisBlock,nChBlock,dovis,nchan,nleft); // DEBUG
	   gridKernel<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>
              (gpu, off, dovis, (CUDAGridInfo *)cudaInfo->d_base[gpu], d_grid_debug);
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
//memsize = sizeof(float); // just to give a continuation of the routine

    return;
} // end grid_processWithStreams

/* Public: Grid a bufferload of visibilities */
extern "C"
void ObitCUDAGridGrid (CUDAGridInfo *cudaInfo, long iGPU, int init) {
  if (init) UpdateGPUBase (iGPU, cudaInfo);
  grid_processWithStreams(GRID_STREAM_COUNT, cudaInfo, iGPU,
                          cudaInfo->h_gridInfo->stream[iGPU], cudaInfo->h_gridInfo->cycleDone[iGPU]);
} /* end ObitCUDAGridGrid   */

/* Public: Fold negative u to positive */
extern "C"
void ObitCUDAGridFlip (CUDAGridInfo *cudaInfo, long iGPU) {
   GridInfo*   h_gridInfo  = cudaInfo->h_gridInfo;
   FacetInfo** h_facetInfo = cudaInfo->h_facetInfo;
   int i, nrow, nrowBlock=32, nplane, nplaneBlock=16;
   dim3 numBlocks, thPerBlock;
   // Add conjugate rows

   // loop over facets, split by blocks of rows and planes
   nplane = h_gridInfo->nplane;
   nrow = h_facetInfo[0]->ny;

   for (i=0; i<h_gridInfo->nfacet; i++) {
     // Is this facet on this GPU?
     if (iGPU==cudaInfo->FacetGPU[i]) {  // This one in this GPU?
       nrow = 1+h_facetInfo[i]->ny/2;
       numBlocks.x = (long)(0.9999+nrow/((float)nrowBlock)); 
       numBlocks.y = (long)(0.9999+nplane/((float)nplaneBlock)); 
       thPerBlock.x = nrowBlock;thPerBlock.y = nplaneBlock;
       if (nrow<=nrowBlock)     {numBlocks.x = 1; thPerBlock.x = nrow;}
       if (nplane<=nplaneBlock) {numBlocks.y = 1; thPerBlock.y = nplane;}
       numBlocks.x = MAX(1,numBlocks.x); numBlocks.y = MAX(1,numBlocks.y);
       gridFlipKernel<<<numBlocks,thPerBlock>>>
          ((int)iGPU, (CUDAGridInfo *)cudaInfo->d_base[iGPU], i, cudaInfo->d_grid_debug);
       } // end if this GPU
    } // end facet loop

     cudaDeviceSynchronize();
} // end ObitCUDAGridFlip 

/* Public: Copy  vis grid to locked host memory */
extern "C"
void ObitCUDAGrid2CPU (CUDAGridInfo *cudaInfo, long ifacet) {
   GridInfo* gridInfo    = cudaInfo->h_gridInfo;
   //long iGPU =  cudaInfo->FacetGPU[ifacet]; // which gpu index am I
   FacetInfo** facetInfo = cudaInfo->h_facetInfo;
   size_t ms = facetInfo[ifacet]->sizeGrid*gridInfo->nplane*sizeof(float);
   checkCudaErrors(cudaMemcpy(facetInfo[ifacet]->h_grid, facetInfo[ifacet]->d_grid, 
                   ms, cudaMemcpyDeviceToHost));
} /* end ObitCUDAGrid2CPU */

/* Public: Shutdown gridding */
extern "C"
void ObitCUDAGridShutdown (CUDAGridInfo *cudaInfo, long iGPU) {
  GridInfo* h_gridInfo    = cudaInfo->h_gridInfo;
  FacetInfo** h_facetInfo = cudaInfo->h_facetInfo;
  CUDAGridInfo *t_cudaInfo=NULL;
  int nGPU = cudaInfo->h_gridInfo->nGPU;
  int i, last = (iGPU==(nGPU-1));               // Last GPU? get host stuff

  // Free allocated debug memory
  if (cudaInfo->h_grid_debug) {cudaFreeHost(cudaInfo->h_grid_debug); cudaInfo->h_grid_debug = NULL;}
  if (cudaInfo->d_grid_debug) {cudaFree(cudaInfo->d_grid_debug);     cudaInfo->d_grid_debug = NULL;}

  // Get base structure (cudaInfo) for device 
  t_cudaInfo = (CUDAGridInfo*)cudaInfo->h_base[iGPU];

  GridInfo* t_gridInfo;              // temp host copy of device grid specifications
  checkCudaErrors(cudaMallocHost(&t_gridInfo, sizeof(GridInfo)));
  // get host copy of device gridInfo
  checkCudaErrors(cudaMemcpy(t_gridInfo, t_cudaInfo->d_gridInfo, sizeof(GridInfo), cudaMemcpyDeviceToHost));

  // free gridInfo arrays
  if (t_gridInfo->d_freqArr)    {cudaFree (t_gridInfo->d_freqArr);   t_gridInfo->d_freqArr = NULL;}
  if (t_gridInfo->d_freqPlane)  {cudaFree (t_gridInfo->d_freqPlane); t_gridInfo->d_freqPlane = NULL;}
  if (t_gridInfo->d_sigma1)     {cudaFree (t_gridInfo->d_sigma1);    t_gridInfo->d_sigma1 = NULL;}
  if (t_gridInfo->d_sigma2)     {cudaFree (t_gridInfo->d_sigma2);    t_gridInfo->d_sigma2 = NULL;}
  if (t_gridInfo->d_sigma3)     {cudaFree (t_gridInfo->d_sigma3);    t_gridInfo->d_sigma3 = NULL;}
  if (t_gridInfo->d_convfn)     {cudaFree (t_gridInfo->d_convfn);    t_gridInfo->d_convfn = NULL;}
  if (t_cudaInfo->d_gridInfo)   {cudaFree(t_cudaInfo->d_gridInfo);   t_cudaInfo->d_gridInfo = NULL;}  // device GridInfo

  // device data buffer
  if (t_gridInfo->d_data_in)  {cudaFree    (t_gridInfo->d_data_in); t_gridInfo->d_data_in = NULL;}
  if (t_gridInfo) cudaFreeHost(t_gridInfo);  // Free temp array

  // streams - only on host but gpu dependent
  for (i=0; i<GRID_STREAM_COUNT; i++) {
    cudaStreamDestroy((cudaStream_t)h_gridInfo->stream[iGPU][i]);
    cudaEventDestroy(h_gridInfo->cycleDone[iGPU][i]);
  }

  // this GPU base
  if (cudaInfo->d_base[iGPU]) {cudaFree (cudaInfo->d_base[iGPU]); cudaInfo->d_base[iGPU]=NULL;}

  // this GPU d_facetInfo
  // Facet stuff 
  for (i=0; i<h_gridInfo->nfacet; i++) {
     if (cudaInfo->FacetGPU[i]==iGPU) {  // Is this facet in this GPU? If so free d_grid array 
       if (h_facetInfo[i]->d_grid) {cudaFree (h_facetInfo[i]->d_grid); h_facetInfo[i]->d_grid = NULL;}
       }
   }
  //Not visible on host if (cudaInfo->d_facetInfo[iGPU])   {cudaFree (cudaInfo->d_facetInfo[iGPU]);       cudaInfo->d_facetInfo[iGPU]=NULL;}
  if (cudaInfo->h_d_facetInfo[iGPU]) {cudaFreeHost (cudaInfo->h_d_facetInfo[iGPU]); cudaInfo->h_d_facetInfo[iGPU]=NULL;}
  
  // Get device arrays 
  if (t_cudaInfo->d_facetInfo) { // d_facet array
      cudaFree (t_cudaInfo->d_facetInfo); t_cudaInfo->d_facetInfo = NULL;
  }

  if (t_cudaInfo->d_gridInfo) {  // d_gridInfo array
    for (i=0; i<cudaInfo->h_gridInfo->nGPU; i++) cudaFree (cudaInfo->h_d_gridInfo[i]); 
    cudaFree (t_cudaInfo->d_gridInfo); cudaInfo->d_gridInfo  = NULL;
  }
  if (last) {  
    // Free host arrays on last call
    // host data buffer
    if (h_gridInfo->h_data_in)    {cudaFreeHost(h_gridInfo->h_data_in);    h_gridInfo->h_data_in = NULL;}
    // Other arrays
    if (h_gridInfo->h_freqArr)    {cudaFreeHost(h_gridInfo->h_freqArr);    h_gridInfo->h_freqArr = NULL;}
    if (h_gridInfo->h_freqPlane)  {cudaFreeHost(h_gridInfo->h_freqPlane);  h_gridInfo->h_freqPlane = NULL;}
    if (h_gridInfo->h_sigma1)     {cudaFreeHost(h_gridInfo->h_sigma1);     h_gridInfo->h_sigma1 = NULL;}
    if (h_gridInfo->h_sigma2)     {cudaFreeHost(h_gridInfo->h_sigma2);     h_gridInfo->h_sigma2 = NULL;}
    if (h_gridInfo->h_sigma3)     {cudaFreeHost(h_gridInfo->h_sigma3);     h_gridInfo->h_sigma3 = NULL;}
    if (h_gridInfo->h_convfn)     {cudaFreeHost(h_gridInfo->h_convfn);     h_gridInfo->h_convfn = NULL;}
    if (h_gridInfo->GPU_device_no){cudaFreeHost(h_gridInfo->GPU_device_no);h_gridInfo->GPU_device_no = NULL;}
    if (h_gridInfo->stream)       {free(h_gridInfo->stream);             h_gridInfo->stream    = NULL;}
    if (h_gridInfo->cycleDone)    {free(h_gridInfo->cycleDone);          h_gridInfo->cycleDone = NULL;}

    // Facet stuff 
    for (i=0; i<h_gridInfo->nfacet; i++) {
      if (cudaInfo->FacetGPU[i]==iGPU) {  /* Is this facet in this GPU? */
        if (h_facetInfo[i]->h_grid) {cudaFreeHost(h_facetInfo[i]->h_grid); h_facetInfo[i]->h_grid = NULL;}
        if (h_facetInfo[i]) {cudaFreeHost(h_facetInfo[i]); h_facetInfo[i]=NULL;}
      }
     } // end facet loop

     // host facetInfo array and gridInfo
     if (cudaInfo->h_gridInfo)    {cudaFreeHost(cudaInfo->h_gridInfo);    cudaInfo->h_gridInfo    = NULL;}
     if (cudaInfo->h_facetInfo)   {cudaFreeHost(cudaInfo->h_facetInfo);   cudaInfo->h_facetInfo   = NULL;}
     if (cudaInfo->h_d_facetInfo) {cudaFreeHost(cudaInfo->h_d_facetInfo); cudaInfo->h_d_facetInfo = NULL;}
     if (cudaInfo->h_d_gridInfo)  {cudaFreeHost(cudaInfo->h_d_gridInfo);  cudaInfo->h_d_gridInfo = NULL;}

     // Others
     if (cudaInfo->d_base) {   // d_base array in host
      cudaFreeHost (cudaInfo->d_base); cudaInfo->d_base = NULL;
     }

     // Reset GPU
     cudaDeviceSynchronize();
     cudaDeviceReset();

     }  // end if last  
  } /* end ObitCUDAGridShutdown */
#endif // CUDA code

