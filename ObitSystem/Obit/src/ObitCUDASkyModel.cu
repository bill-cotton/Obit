/* Still need 
5) Grid
*/
/* $Id: $        */
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

/*#include "ObitCUDASkyModel.h"*/
#include "ObitCUDASkyModelInfoDef.h"
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDASkyModel.cu
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

/* DEBUG */
float *h_debug, *d_debug;
__global__ void debugKernal( float* __restrict__ d_debug, GPUVisInfo* visInfo)
{
    int iprod      = threadIdx.x+256*blockIdx.y;     // product number
    // No more than actual number of products
    if (iprod>=visInfo->nprod) return;
    d_debug[0] = blockIdx.x;
    d_debug[1] = blockIdx.y;
    d_debug[2] = threadIdx.x;
    d_debug[3] = threadIdx.y;
} /* end debugKernel */

/**
 * Point DFT GPU kernal.
 * block = vis, threads in block = data product = channel/stokes/IF
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 */
__global__ void dftPointKernel(float* __restrict__ g_data, 
	   GPUModelInfo* modelInfo,  GPUVisInfo* visInfo)
{
    int lenvis     = visInfo->lenvis;
    int idx        = blockIdx.x * lenvis;        // beginning of a visibility
    int iprod      = threadIdx.x+256*blockIdx.y; // product number
    int nrparm     = visInfo->nrparm;
    float *FreqArr = visInfo->freqScale;
    int modelSize  = modelInfo->size;
    int nModel;
    float *Model;
    int ichan, istok, iIF, ivis, iMod;
    float arg, amp, s, c, sumR, sumI; 
    float u, v, w, tr, ti;
    float freqFact;

    // No more than actual number of products
    if (iprod>=visInfo->nprod) return;

    // get channel,stokes, IF from data product
    ichan = (iprod / visInfo->incf)  % visInfo->nchan;
    istok = (iprod / visInfo->incs)  % visInfo->nstok;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;

    // real part of vis
    ivis = idx + nrparm + iprod*3;
 
    // This one desired?
    if ((ichan<visInfo->chanb) || (ichan>visInfo->chane) ||
        (istok<visInfo->stokb) || (istok>visInfo->stoke) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) || (g_data[ivis+2]<=0.0)) {
	return;
    }

    // Model per poln 
    nModel = modelInfo->nmodel[istok];
    Model  = modelInfo->model[istok];

    // frequency scaling factor	  
    freqFact = FreqArr[ichan*visInfo->kincf+iIF*visInfo->kincif];

    // get scaled u,v,w factors
    u = g_data[idx+visInfo->ilocu];
    v = g_data[idx+visInfo->ilocv];
    w = g_data[idx+visInfo->ilocw];
    u *= freqFact;
    v *= freqFact;
    w *= freqFact;
    sumR = sumI = 0.0;
 
    // model = flux, x,y,z factors
    iMod = 0;
    for (int i=0; i<nModel; i++) {
	amp = Model[iMod];
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
        sumR += amp * c;
        sumI += amp * s;
        iMod += modelSize;
    } // end loop over model comps

    // Multiply by Factor
    sumR *= modelInfo->stokFact[istok];
    sumI *= modelInfo->stokFact[istok];
    // Swap R/-I for LR?
    if ((modelInfo->doSwap) && (istok==3)) {
      arg  =  sumR;
      sumR = -sumI;
      sumI =  arg;
    }
    // Correct Data by operation type 
    if (modelInfo->opType==GPUOpTypeRepl) {      // Replace
      g_data[ivis]   = sumR;
      g_data[ivis+1] = sumI;
    } else if (modelInfo->opType==GPUOpTypeDiv) { // Divide
      amp = sumR*sumR + sumI*sumI;
      sumR /= amp;
      sumI /= amp;
      amp = sqrt(amp);
      tr = g_data[ivis];
      ti = g_data[ivis+1];
      g_data[ivis]   = sumR*tr + sumI*ti;
      g_data[ivis+1] = sumR*ti - sumI*tr;
      g_data[ivis+2] *= amp;
    } else {                                     // Subtract
      g_data[ivis]   -= sumR;
      g_data[ivis+1] -= sumI;
    } // end opType
 } // end dftPointKernel

/**
 * Point DFT with spectrum GPU kernal.
 * block = vis, threads in block = data product = channel/stokes/IF
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 */
__global__ void dftPointSpecKernel(float* __restrict__ g_data, 
	   GPUModelInfo* modelInfo,  GPUVisInfo* visInfo)
{
    int lenvis       = visInfo->lenvis;
    int idx          = blockIdx.x * lenvis;        // beginning of a visibility
    int iprod        = threadIdx.x+256*blockIdx.y; // product number
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->freqScale;
    float freqRat    = visInfo->d_freqRat[0];
    int modelSize    = modelInfo->size;
    int iterm, nterm = modelInfo->nterm;
    int nModel;
    float *Model;
    int ichan, istok, iIF, ivis, iMod;
    float arg, amp, s, c, sumR, sumI; 
    float u, v, w, tr, ti;
    float freqFact, lll, lnspecFreqFact;

    // No more than actual number of products
    if (iprod>=visInfo->nprod) return;

    // get channel,stokes, IF from data product
    ichan = (iprod / visInfo->incf)  % visInfo->nchan;
    istok = (iprod / visInfo->incs)  % visInfo->nstok;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;

    // real part of vis
    ivis = idx + nrparm + iprod*3;
 
   // This one desired?
    if ((ichan<visInfo->chanb) || (ichan>visInfo->chane) ||
        (istok<visInfo->stokb) || (istok>visInfo->stoke) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) || (g_data[ivis+2]<=0.0)) {
	return;
    }

    // Model per poln 
    nModel = modelInfo->nmodel[istok];
    Model  = modelInfo->model[istok];

    // frequency scaling factor	  
    freqFact = FreqArr[ichan*visInfo->kincf+iIF*visInfo->kincif];
    // log of ratio to reference freq
    lnspecFreqFact = __logf(freqFact*freqRat);

    // get scaled u,v,w factors
    u = g_data[idx+visInfo->ilocu];
    v = g_data[idx+visInfo->ilocv];
    w = g_data[idx+visInfo->ilocw];
    u *= freqFact;
    v *= freqFact;
    w *= freqFact;
    sumR = sumI = 0.0;

    // model = flux, x,y,z factors, spectral terms
    iMod = 0;
    for (int i=0; i<nModel; i++) {
        // amplitude
	amp = Model[iMod];
        // phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
        // Frequency dependent spectral term */
        lll = lnspecFreqFact;
        arg = 0.0;
        for (iterm=0; iterm<nterm; iterm++) {
	  arg += Model[iMod+4+iterm] * lll;
	  lll *= lnspecFreqFact;
	}
	amp *= __expf(-arg);
        sumR += amp * c;
        sumI += amp * s;
        iMod += modelSize;
    } // end loop over model comps

    // Multiply by Factor
    sumR *= modelInfo->stokFact[istok];
    sumI *= modelInfo->stokFact[istok];
    // Swap R/-I for LR?
    if ((modelInfo->doSwap) && (istok==3)) {
      arg  =  sumR;
      sumR = -sumI;
      sumI =  arg;
    }
    // Correct Data by operation type 
    if (modelInfo->opType==GPUOpTypeRepl) {      // Replace
      g_data[ivis]   = sumR;
      g_data[ivis+1] = sumI;
    } else if (modelInfo->opType==GPUOpTypeDiv) { // Divide
      amp = sumR*sumR + sumI*sumI;
      sumR /= amp;
      sumI /= amp;
      amp = sqrt(amp);
      tr = g_data[ivis];
      ti = g_data[ivis+1];
      g_data[ivis]   = sumR*tr + sumI*ti;
      g_data[ivis+1] = sumR*ti - sumI*tr;
      g_data[ivis+2] *= amp;
    } else {                                     // Subtract
      g_data[ivis]   -= sumR;
      g_data[ivis+1] -= sumI;
    } // end opType
 } // end dftPointSpecKernel

/**
 * Point DFT with tabulated spectrum GPU kernal.
 * block = vis, threads in block = data product = channel/stokes/IF
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  debug      debug array
 */
__global__ void dftPointTSpecKernel(float* __restrict__ g_data,
	   GPUModelInfo* modelInfo,  GPUVisInfo* visInfo)
//, float* __restrict__ debug)
{
    int lenvis       = visInfo->lenvis;
    int idx          = blockIdx.x * lenvis;        // beginning of a visibility
    int iprod        = threadIdx.x+256*blockIdx.y; // product number
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->freqScale;
    float *freqRat   = visInfo->d_freqRat;
    int modelSize    = modelInfo->size;
    int *specIndex   = modelInfo->specIndex;
    int nModel;
    float *Model;
    int ichan, istok, iIF, ivis, itab, iMod, ifq;
    float arg, amp, s, c, sumR, sumI; 
    float u, v, w, tr, ti;
    float freqFact, lnspecFreqFact;

    // No more than actual number of products
    if (iprod>=visInfo->nprod) return;

    // get channel,stokes, IF from data product
    ichan = (iprod / visInfo->incf)  % visInfo->nchan;
    istok = (iprod / visInfo->incs)  % visInfo->nstok;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;

    // real part of vis
    ivis = idx + nrparm + iprod*3;
 
    // This one desired?
    if ((ichan<visInfo->chanb) || (ichan>visInfo->chane) ||
        (istok<visInfo->stokb) || (istok>visInfo->stoke) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) || (g_data[ivis+2]<=0.0)) {
	return;
    }

    // Model per poln 
    nModel = modelInfo->nmodel[istok];
    Model  = modelInfo->model[istok];

    // frequency scaling factor	  
    ifq = ichan*visInfo->kincf+iIF*visInfo->kincif;
    itab = 5 + specIndex[ifq]; // which coarse channel?
    freqFact = FreqArr[ifq];   // channel frequency scaling factor
    // log of ratio to reference freq
    lnspecFreqFact = -__logf(freqFact*freqRat[itab-5]);

    // get scaled u,v,w factors
    u = g_data[idx+visInfo->ilocu];
    v = g_data[idx+visInfo->ilocv];
    w = g_data[idx+visInfo->ilocw];
    u *= freqFact;
    v *= freqFact;
    w *= freqFact;
    sumR = sumI = 0.0;

    // model = [flux], x,y,z factors, spectral index, coarse channel fluxes
    iMod = 0;
    for (int i=0; i<nModel; i++) {
        // phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
        // include spectral index
	arg = lnspecFreqFact*Model[iMod+4];
	amp =  Model[iMod+itab] * __expf(arg);
        sumR += amp * c;
        sumI += amp * s;
        iMod += modelSize;
    } // end loop over model comps

    // Multiply by Factor
    sumR *= modelInfo->stokFact[istok];
    sumI *= modelInfo->stokFact[istok];
    // Swap R/-I for LR?
    if ((modelInfo->doSwap) && (istok==3)) {
      arg  =  sumR;
      sumR = -sumI;
      sumI =  arg;
    }
    // Correct Data by operation type 
    if (modelInfo->opType==GPUOpTypeRepl) {      // Replace
      g_data[ivis]   = sumR;
      g_data[ivis+1] = sumI;
    } else if (modelInfo->opType==GPUOpTypeDiv) { // Divide
      amp = sumR*sumR + sumI*sumI;
      sumR /= amp;
      sumI /= amp;
      amp = sqrt(amp);
      tr = g_data[ivis];
      ti = g_data[ivis+1];
      g_data[ivis]   = sumR*tr + sumI*ti;
      g_data[ivis+1] = sumR*ti - sumI*tr;
      g_data[ivis+2] *= amp;
    } else {                                     // Subtract
      g_data[ivis]   -= sumR;
      g_data[ivis+1] -= sumI;
    } // end opType
 } // end dftPointTSpecKernel

/**
 * Gauss DFT GPU kernal.
 * block = vis, threads in block = data product = channel/stokes/IF
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 */
__global__ void dftGaussKernel(float* __restrict__ g_data,
	   GPUModelInfo* modelInfo,  GPUVisInfo* visInfo)
{
    int lenvis     = visInfo->lenvis;
    int idx        = blockIdx.x * lenvis;        // beginning of a visibility
    int iprod      = threadIdx.x+256*blockIdx.y; // product number
    int nrparm     = visInfo->nrparm;
    float *FreqArr = visInfo->freqScale;
    int modelSize  = modelInfo->size;
    int nModel;
    float *Model;
    int ichan, istok, iIF, ivis, iMod;
    float arg, amp, s, c, sumR, sumI; 
    float u, v, w, tr, ti;
    float freqFact;

    // No more than actual number of products
    if (iprod>=visInfo->nprod) return;

    // get channel,stokes, IF from data product
    ichan = (iprod / visInfo->incf)  % visInfo->nchan;
    istok = (iprod / visInfo->incs)  % visInfo->nstok;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;

    // real part of vis
    ivis = idx + nrparm + iprod*3;
 
    // This one desired?
    if ((ichan<visInfo->chanb) || (ichan>visInfo->chane) ||
        (istok<visInfo->stokb) || (istok>visInfo->stoke) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) || (g_data[ivis+2]<=0.0)) {
	return;
    }

    // Model per poln 
    nModel = modelInfo->nmodel[istok];
    Model  = modelInfo->model[istok];

    // frequency scaling factor	  
    freqFact = FreqArr[ichan*visInfo->kincf+iIF*visInfo->kincif];

    // get scaled u,v,w factors
    u = g_data[idx+visInfo->ilocu];
    v = g_data[idx+visInfo->ilocv];
    w = g_data[idx+visInfo->ilocw];
    u *= freqFact;
    v *= freqFact;
    w *= freqFact;
    sumR = sumI = 0.0;
 
    // model = flux, x,y,z factors, uu, vv, uv factors
    iMod = 0;
    for (int i=0; i<nModel; i++) {
	amp = Model[iMod];
  	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
        arg = u*u*Model[iMod+4] + v*v*Model[iMod+5] + u*v*Model[iMod+6];
	amp *= __expf(arg);
        sumR += amp * c;
        sumI += amp * s;
        iMod += modelSize;
    } // end loop over model comps

    // Multiply by Factor
    sumR *= modelInfo->stokFact[istok];
    sumI *= modelInfo->stokFact[istok];
    // Swap R/-I for LR?
    if ((modelInfo->doSwap) && (istok==3)) {
      arg  =  sumR;
      sumR = -sumI;
      sumI =  arg;
    }
    // Correct Data by operation type 
    if (modelInfo->opType==GPUOpTypeRepl) {      // Replace
      g_data[ivis]   = sumR;
      g_data[ivis+1] = sumI;
    } else if (modelInfo->opType==GPUOpTypeDiv) { // Divide
      amp = sumR*sumR + sumI*sumI;
      sumR /= amp;
      sumI /= amp;
      amp = sqrt(amp);
      tr = g_data[ivis];
      ti = g_data[ivis+1];
      g_data[ivis]   = sumR*tr + sumI*ti;
      g_data[ivis+1] = sumR*ti - sumI*tr;
      g_data[ivis+2] *= amp;
    } else {                                     // Subtract
      g_data[ivis]   -= sumR;
      g_data[ivis+1] -= sumI;
    } // end opType
 } // end dftGaussKernel

/**
 * Gauss DFT with spectrum GPU kernal.
 * block = vis, threads in block = data product = channel/stokes/IF
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 */
__global__ void dftGaussSpecKernel(float* __restrict__ g_data,
	   GPUModelInfo* modelInfo,  GPUVisInfo* visInfo)
{
    int lenvis       = visInfo->lenvis;
    int idx          = blockIdx.x * lenvis;        // beginning of a visibility
    int iprod        = threadIdx.x+256*blockIdx.y; // product number
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->freqScale;
    float freqRat    = visInfo->d_freqRat[0];
    int modelSize    = modelInfo->size;
    int iterm, nterm = modelInfo->nterm;
    int nModel;
    float *Model;
    int ichan, istok, iIF, ivis, iMod;
    float arg, amp, s, c, sumR, sumI; 
    float u, v, w, tr, ti;
    float freqFact, lll, lnspecFreqFact;

    // No more than actual number of products
    if (iprod>=visInfo->nprod) return;

    // get channel,stokes, IF from data product
    ichan = (iprod / visInfo->incf)  % visInfo->nchan;
    istok = (iprod / visInfo->incs)  % visInfo->nstok;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;

    // real part of vis
    ivis = idx + nrparm + iprod*3;
 
    // This one desired?
    if ((ichan<visInfo->chanb) || (ichan>visInfo->chane) ||
        (istok<visInfo->stokb) || (istok>visInfo->stoke) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) || (g_data[ivis+2]<=0.0)) {
	return;
    }

    // Model per poln 
    nModel = modelInfo->nmodel[istok];
    Model  = modelInfo->model[istok];

    // frequency scaling factor	  
    freqFact = FreqArr[ichan*visInfo->kincf+iIF*visInfo->kincif];
    // log of ratio to reference freq
    lnspecFreqFact = __logf(freqFact*freqRat);

   // get scaled u,v,w factors
    u = g_data[idx+visInfo->ilocu];
    v = g_data[idx+visInfo->ilocv];
    w = g_data[idx+visInfo->ilocw];
    u *= freqFact;
    v *= freqFact;
    w *= freqFact;
    sumR = sumI = 0.0;

    // model = flux, x,y,z factors, uu, vv, uv factors, spectral terms
    iMod = 0;
    for (int i=0; i<nModel; i++) {
        // amplitude
	amp = Model[iMod];
  	// Phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
        arg  = u*u*Model[iMod+4] + v*v*Model[iMod+5] + u*v*Model[iMod+6];
	amp *= __expf(arg);
        // Frequency dependent spectral term */
        lll = lnspecFreqFact;
        arg = 0.0;
        for (iterm=0; iterm<nterm; iterm++) {
	  arg += Model[iMod+7+iterm] * lll;
	  lll *= lnspecFreqFact;
	}
	amp *= __expf(-arg);
        // Gaussian factor
        sumR += amp * c;
        sumI += amp * s;
        iMod += modelSize;
    } // end loop over model comps

    // Multiply by Factor
    sumR *= modelInfo->stokFact[istok];
    sumI *= modelInfo->stokFact[istok];
    // Swap R/-I for LR?
    if ((modelInfo->doSwap) && (istok==3)) {
      arg  =  sumR;
      sumR = -sumI;
      sumI =  arg;
    }
    // Correct Data by operation type 
    if (modelInfo->opType==GPUOpTypeRepl) {      // Replace
      g_data[ivis]   = sumR;
      g_data[ivis+1] = sumI;
    } else if (modelInfo->opType==GPUOpTypeDiv) { // Divide
      amp = sumR*sumR + sumI*sumI;
      sumR /= amp;
      sumI /= amp;
      amp = sqrt(amp);
      tr = g_data[ivis];
      ti = g_data[ivis+1];
      g_data[ivis]   = sumR*tr + sumI*ti;
      g_data[ivis+1] = sumR*ti - sumI*tr;
      g_data[ivis+2] *= amp;
    } else {                                     // Subtract
      g_data[ivis]   -= sumR;
      g_data[ivis+1] -= sumI;
    } // end opType
 } // end dftGaussSpecKernel

/**
 * Gauss DFT with tabulated spectrum GPU kernal.
 * block = vis, threads in block = data product = channel/stokes/IF
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  debug      debug array
 */
__global__ void dftGaussTSpecKernel(float* __restrict__ g_data,
	   GPUModelInfo* modelInfo,  GPUVisInfo* visInfo)
//, float* __restrict__ debug)
{
    int lenvis       = visInfo->lenvis;
    int idx          = blockIdx.x * lenvis;        // beginning of a visibility
    int iprod        = threadIdx.x+256*blockIdx.y; // product number
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->freqScale;
    float *freqRat   = visInfo->d_freqRat;
    int modelSize    = modelInfo->size;
    int *specIndex   = modelInfo->specIndex;
    int nModel;
    float *Model;
    int ichan, istok, iIF, ivis, itab, iMod, ifq;
    float arg, amp, s, c, sumR, sumI; 
    float u, v, w, tr, ti;
    float freqFact, lnspecFreqFact;

    // No more than actual number of products
    if (iprod>=visInfo->nprod) return;

    // get channel,stokes, IF from data product
    ichan = (iprod / visInfo->incf)  % visInfo->nchan;
    istok = (iprod / visInfo->incs)  % visInfo->nstok;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;

    // real part of vis
    ivis = idx + nrparm + iprod*3;
 
    // This one desired?
    if ((ichan<visInfo->chanb) || (ichan>visInfo->chane) ||
        (istok<visInfo->stokb) || (istok>visInfo->stoke) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) || (g_data[ivis+2]<=0.0)) {
	return;
    }

    // Model per poln 
    nModel = modelInfo->nmodel[istok];
    Model  = modelInfo->model[istok];

    // frequency scaling factor	  
    ifq = ichan*visInfo->kincf+iIF*visInfo->kincif;
    itab = 8 + specIndex[ifq]; // which coarse channel?
    freqFact = FreqArr[ifq];   // channel frequency scaling factor
    // log of ratio to reference freq
    lnspecFreqFact = -__logf(freqFact*freqRat[itab-8]);

    // get scaled u,v,w factors
    u = g_data[idx+visInfo->ilocu];
    v = g_data[idx+visInfo->ilocv];
    w = g_data[idx+visInfo->ilocw];
    u *= freqFact;
    v *= freqFact;
    w *= freqFact;
    sumR = sumI = 0.0;

    // model = [flux], x,y,z factors, uu, vv, uv factors, spectral index, coarse channel fluxes
    iMod = 0;
    for (int i=0; i<nModel; i++) {
        // Phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
        // Gaussian factor
        arg = u*u*Model[iMod+4] + v*v*Model[iMod+5] + u*v*Model[iMod+6];
	amp  =__expf(arg);
        // include spectral index
	arg = lnspecFreqFact*Model[iMod+7];
	amp *= Model[iMod+itab] * __expf(arg);
        sumR += amp * c;
        sumI += amp * s;
        iMod += modelSize;
    } // end loop over model comps

    // Multiply by Factor
    sumR *= modelInfo->stokFact[istok];
    sumI *= modelInfo->stokFact[istok];
    // Swap R/-I for LR?
    if ((modelInfo->doSwap) && (istok==3)) {
      arg  =  sumR;
      sumR = -sumI;
      sumI =  arg;
    }
    // Correct Data by operation type 
    if (modelInfo->opType==GPUOpTypeRepl) {      // Replace
      g_data[ivis]   = sumR;
      g_data[ivis+1] = sumI;
    } else if (modelInfo->opType==GPUOpTypeDiv) { // Divide
      amp = sumR*sumR + sumI*sumI;
      sumR /= amp;
      sumI /= amp;
      amp = sqrt(amp);
      tr = g_data[ivis];
      ti = g_data[ivis+1];
      g_data[ivis]   = sumR*tr + sumI*ti;
      g_data[ivis+1] = sumR*ti - sumI*tr;
      g_data[ivis+2] *= amp;
    } else {                                     // Subtract
      g_data[ivis]   -= sumR;
      g_data[ivis+1] -= sumI;
    } // end opType
 } // end dftGaussTSpecKernel
#endif /* HAVE_GPU */

#if HAVE_GPU==1  /* CUDA code */
static void DFTprocessWithStreams(int streams_used, int nvis, 
            GPUModelInfo* h_modelInfo, GPUModelInfo* d_modelInfo, 
            GPUVisInfo* h_visInfo, GPUVisInfo* d_visInfo, 
            cudaStream_t* stream, cudaEvent_t* cycleDone,
            float *h_data_source, float *h_data_sink, 
            float *d_data[], int prtLv);
#endif /* HAVE_GPU */

/*----------------------Public functions---------------------------*/
/**
 * Initialize an ObitCUDASkyModel 
 * Currently only DFT point supported
 * \param gpuInfo    processing info
 */
extern "C"
void ObitCUDASkyModelDFTInit (GPUModelInfo *gpuInfo)
{
} /* end  ObitCUDASkyModelDFTInit */


/**
 * Setup model for ObitCUDASkyModel 
 * Currently only DFT point supported
 * \param gpuInfo    processing info
 */
extern "C"
void ObitCUDASkyModelDFTSetMod (GPUModelInfo *gpuInfo)
{
#if HAVE_GPU==1  /* CUDA code */
#endif /* HAVE_GPU */
    return;
} /* end ObitCUDASkyModelDFTSetMod */

/**
 * Calculate an ObitCUDASkyModel 
 * Currently only DFT point supported
 * \param gpuInfo    processing info
 */
extern "C"
void ObitCUDASkyModelDFTCalc (GPUInfo *gpuInfo)
{
#if HAVE_GPU==1  /* CUDA code */
   // Put available fast memory in L1 cache - no apparent effect
   checkCudaErrors(cudaFuncSetCacheConfig(dftPointKernel, cudaFuncCachePreferL1));

   /* Process with streams */
   DFTprocessWithStreams(gpuInfo->nstream, gpuInfo->nvis,gpuInfo->h_modelInfo, gpuInfo->d_modelInfo, 
                        gpuInfo->h_visInfo, gpuInfo->d_visInfo, 
                        (cudaStream_t *)gpuInfo->stream, (cudaEvent_t *)gpuInfo->cycleDone, 
                        gpuInfo->h_data, gpuInfo->h_data, gpuInfo->d_data, 0);
#endif /* HAVE_GPU */
    return;
} /* end ObitCUDASkyModelDFTCalc */


/**
 * Shutdown  ObitCUDASkyModel 
 * Currently only DFT point supported
 * \param gpuInfo    processing info
 */
extern "C"
void ObitCUDASkyModelDFTShutdown (GPUInfo *gpuInfo)
{

#if HAVE_GPU==1  /* CUDA code */
#endif /* HAVE_GPU */

} /*end ObitCUDASkyModelDFTShutdown */
#if HAVE_GPU==1  /* CUDA code */
/**
 * Calculates DFT model using a GPU.
 * Multiple streams are used to overlap I/O and computation,
 * each call divides the data into streams_used pieces.
 * \param  streams_used  Number of streams to use
 * \param  nvis          Number of visibilities
 * \param  nmodel        Number of sky model components
 * \param  h_modelInfo   Host resident modelInfo
 * \param  d_modelInfo   GPU resident modelInfo
 * \param  h_visInfo     Host resident visInfo
 * \param  d_visInfo     GPU resident visInfo
 * \param  stream        GPU stream array [streams_used]
 * \param  cycleDone     GPU event array [streams_used]
 * \param  h_data_source Host resident input buffer, should be locked
 * \param  h_data_sink   Host resident output buffer, should be locked
 *                       may be h_data_source
 * \param  d_data        GPU resident data buffer (nvis/streams_used)
 * \param  prtLv         Print level, 5=>much directly printed
 */
static void DFTprocessWithStreams(int streams_used, int nvis, 
            GPUModelInfo* h_modelInfo, GPUModelInfo* d_modelInfo, 
            GPUVisInfo* h_visInfo, GPUVisInfo* d_visInfo, 
            cudaStream_t* stream, cudaEvent_t* cycleDone,
            float *h_data_source, float *h_data_sink, 
            float *d_data[], int prtLv)
{

    int  last_stream, current_stream = 0;
    int npass = streams_used;
    int nvisPass = (nvis+npass-1)/npass;  // nearly round up
    int lenvis = h_visInfo->lenvis;
    int off, nprod, dovis;
    int memsize = (lenvis*nvisPass)*sizeof(float);
    dim3 numBlocks, thPerBlock;

    // DEBUG  create debug arrays
    //int ms = 10000*sizeof(float);
    //float *d_debug, *h_debug;
    //checkCudaErrors(cudaMalloc(&d_debug, ms));
    //checkCudaErrors(cudaMallocHost(&h_debug, ms));
    //h_debug[0] = 111.111;  h_debug[2] = 123456.;
    //checkCudaErrors(cudaMemcpyAsync(d_debug,h_debug,  ms, cudaMemcpyHostToDevice,0));

    // Number of data products
    nprod = h_visInfo->nchan * h_visInfo->nIF * h_visInfo->nstok;

    if (prtLv>=5) printf ("Start\n");
 
    // Do processing in a loop
    //
    // Note: All memory commands are processed in the order  they are issued,
    // independent of the stream they are enqueued in. Hence the pattern by
    // which the copy and kernel commands are enqueued in the stream
    // has an influence on the achieved overlap.


    // Upload first frame
    if (prtLv>=5) printf ("upload current_stream %d off %d\n",0,0);
    checkCudaErrors(cudaMemcpyAsync(d_data[0],
			  	   &h_data_source[0],
				   memsize,
				   cudaMemcpyHostToDevice,
				   stream[0]));
    //?cudaEventSynchronize(cycleDone[0]);

    for (int i=0; i<npass; ++i) {
        int next_stream = (current_stream + 1) % streams_used;
	int prev_stream = current_stream - 1;
	if (prev_stream<0) prev_stream = streams_used-1;
	off = next_stream*lenvis*nvisPass;  /* Offset in data buffers */
  	if (prtLv>=5) printf ("\n\nLoop %d prev %d current %d next %d\n",i, prev_stream,current_stream,next_stream );

	// Upload next frame
	if (prtLv>=5) printf ("upload next_stream %d off %d\n",next_stream,off);
	checkCudaErrors(cudaMemcpyAsync(d_data[next_stream],
			    &h_data_source[off],
                            memsize,
                            cudaMemcpyHostToDevice,
                            stream[next_stream]));

	// Ensure that processing and copying of the previous cycle has finished
	if (i>0) {
	  off = prev_stream*lenvis*nvisPass;  /* Offset in data buffers */
	  if (prtLv>=5) printf ("download prev_stream %d off %d\n",prev_stream,off);
	  checkCudaErrors(cudaMemcpyAsync(&h_data_sink[off],
					  d_data[prev_stream],
					  memsize,
					  cudaMemcpyDeviceToHost,
					  stream[prev_stream]));
	  if (prtLv>=5) printf ("sync prev_stream %d loop %d\n",prev_stream, i);
	  cudaEventSynchronize(cycleDone[prev_stream]);
	}

        // Process current
	if (prtLv>=5) printf ("Process, nvis, %d nch %d stream %d \n",nvisPass, h_visInfo->nchan, current_stream);
	// make sure to do all visibilities
	if (i==npass-1) dovis = nvis-i*nvisPass;
	else            dovis = nvisPass;

	// package work
	numBlocks.x  = dovis; numBlocks.y = (nprod+255)/256;
	thPerBlock.x = 256;   thPerBlock.y = 1;

	//debugKernal<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>(d_debug, d_visInfo);

        // by model type
        switch (h_modelInfo->type) {
        case GPUModTypePoint:
  	  dftPointKernel<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo);
          break;
        case GPUModTypePointSpec:
  	  dftPointSpecKernel<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo);
          break;
        case GPUModTypePointTSpec:
  	  dftPointTSpecKernel<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              //d_modelInfo, d_visInfo);
              d_modelInfo, d_visInfo);
//, d_debug);
          break;
        case GPUModTypeGauss:
  	  dftGaussKernel<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo);
          break;
        case GPUModTypeGaussSpec:
  	  dftGaussSpecKernel<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo);
          break;
        case GPUModTypeGaussTSpec:
  	  dftGaussTSpecKernel<<<numBlocks, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo);
          break;
          default:
              printf("Model type not supported in GPU",d_modelInfo->type);
          };  // End switch by model type 

	// make sure previous frame done
	if (i>0) cudaEventSynchronize(cycleDone[prev_stream]);

	last_stream = current_stream;
        current_stream = next_stream;
    } /* end loop */

    /* Data from last pass */
    if (prtLv>=5) printf ("sync last_stream %d \n",last_stream);
    cudaEventSynchronize(cycleDone[last_stream]);
    
    // last piece may be uneven size
    dovis   = nvis-nvisPass*(npass-1);
    memsize = (lenvis*dovis)*sizeof(float);
    off = last_stream*lenvis*nvisPass;  /* Offset in data buffers */
    if (prtLv>=5) printf ("download last_stream %d off %d\n",last_stream,off);
    checkCudaErrors(cudaMemcpyAsync(&h_data_sink[off],
				    d_data[last_stream],
				    memsize,
				    cudaMemcpyDeviceToHost,
				    stream[last_stream]));
    if (prtLv>=5) printf ("Finish\n");
    cudaDeviceSynchronize();

    // DEBUG fetch debug arrays from GPU memory
    //checkCudaErrors(cudaMemcpyAsync(h_debug, d_debug, ms, cudaMemcpyDeviceToHost,0));
    //cudaFreeHost(h_debug);
    //cudaFree(d_debug);		 

    return;
} /* end processWithStreams */
#endif /* HAVE_GPU */

