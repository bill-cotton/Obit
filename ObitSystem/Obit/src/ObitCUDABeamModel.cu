/* $Id: $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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

/*--------------- CUDA setup, GPU kernal  ----------------*/
/* This is a CUDA routine */
#if HAVE_GPU==1  /* Only if have a GPU */
#define IS_CUDA 1
#include "CUDASkyGeom.h"
#include "ObitCUDAUtil.h"
#include "ObitCUDABeamModel.h"
#include "ObitCUDABeamModelInfoDef.h"
#include "CUDALinCorSum.h"
#include "CUDAFArray.h"
#include "CUDABeamInterp.h"
#define NCHBLK 256  /* Number of channels per block, 256 OK
	               512 too few resources */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDABeamModel.cu
 * Primitive CUDA routines for skyModel with Beam corrections
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling routines.
 */


// includes, project
#include <helper_cuda.h>
#include <helper_functions.h>  // helper for shared that are common to CUDA SDK samples


/****   Some utilities   ****/
/**
 * Calculate inverse square root of beam gain
 * Uses cos2Beam adapted fromObitBeamShape.c 
 * \param fact   antSize*fudge*DG2RAD*VELIGHT*, fudge=1.052
 * \param freq   Frequency (Hz)
 * \param ang    Offset from pointing in deg.
 * \return 1/sqrt(gain), 100 past 10% point
 */
__device__ __forceinline__ float InvBeam(float fact, float freq, float ang) {
  float arg = fact * freq * ang;
  if (arg>1.02) return 100.0;  
  float div = (1.-4.*(arg*arg));
  float igain = div/(__cosf(3.141592653589793*arg));
  return igain;
} /* end  InvBeam */

/**
 * Allows access to either data with "BASELINE" random parameters or
 * "ANTENNA1" and "ANTENNA2" (with 0.01*(subarray-1)).
 * \param uvDesc  Data descriptor
 * \param buffer  UV data buffer for single visibility
 * \param ant1    First antenna number (0-rel)
 * \param ant2    Second antenna number (0-rel)
 */
__device__ __forceinline__ void GetAnts(int ilocb, int iloca1, int iloca2, int ilocsa,
	   float *buffer, int *ant1, int *ant2)  {
  int cbase;
  if (ilocb>=0) {  /* Baseline */
    cbase = buffer[ilocb];
    *ant1 = (cbase / 256.0) + 0.001;
    *ant2 = (cbase - *ant1 * 256) + 0.001;
  } else { /* Antennas */
    *ant1 = (long)(buffer[iloca1]+0.5);
    *ant2 = (long)(buffer[iloca2]+0.5);
  }
  // zero relative
  *ant1 -= 1;  *ant2 -= 1;
}  /*  end GetAnts */

/* Thoughts
  - use flux density + spectral index per subband, tractable memory usage, faster
    and probably good enough
    Tabulated spectrum Model:
      +0 amp (?)               = single subband/channel flux
      +1,2,3                   = u,v,w
      +4,5,6                   = Gaussian terms
      +5 (Pt) + 2*subband      = flux
      +8 (Gauss) + 2*subband   = flux
      +5 (Pt) + 2*subband+1    = SI
      +8 (Gauss) + 2*subband+1 = SI

  Order (slow to fast) of Jones matrices:ant, component, channel
*/

#include "ObitCUDAUtil.h"
/**
 * Shutdown  ObitCUDABeamModel 
 * Currently only DFT supported
 * \param gpuInfo    processing info
 */
//extern "C"
void ObitCUDABeamModelDFTShutdown (GPUBeamInfo *gpuInfo)
{

   // Reset GPU
   cudaDeviceSynchronize();
   cudaDeviceReset();
} /*end ObitCUDABeamModelDFTShutdown */

/****   Data processing kernels   ****/
/**
 * Point, single flux density DFT with correction GPU kernal.
 * block = vis, threads in block = data product = channel
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  antInto    Antenna beam information
 * \param  nvis       number of visibilities
 */
 extern "C"
__global__ void beamPointKernel(float* __restrict__ g_data,
	   GPUBeamModelInfo* modelInfo,  GPUBeamVisInfo* visInfo,
	   GPUBeamAntInfo* antInfo, int nvis)
{
    int lenvis       = visInfo->lenvis;
    int iprod        = threadIdx.x+NCHBLK*blockIdx.x; // product (channel) number
    int idx          = blockIdx.y * lenvis;           // beginning of a visibility 
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->d_freqScale;
    int modelSize    = modelInfo->size;
    int nModel=modelInfo->nmodel;
    float *Model=modelInfo->d_model;
    long ichan, jchan, kchan, iIF, ivis, iMod, ifq, jindx1, jindx2;
    int mChan  = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    // number of channels in beams
    int ia1, ia2, it1, it2, ilocb=antInfo->ilocb,  iloca1=antInfo->iloca1, iloca2=antInfo->iloca2;
    int ilocsa=antInfo->ilocsa, allBad;
    int notRepl=modelInfo->opType!=GPUBeamOpTypeRepl, gotSome = 0;  // any valid data?
    cudaMatx *Jones = antInfo->d_Jones;
    float arg, amp, s, c; 
    float u, v, w, freqFact;
    cudaMatx SumVis, visMatx;

    // No more than actual number of channels, or visibilities
    if (iprod>=visInfo->nchan)  return;
    if (blockIdx.y>=nvis)       return;

    // get channel, IF from data product
    kchan = (iprod / visInfo->incf)  % visInfo->nchan;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;
    ichan = iprod;

    // real part of first vis (XX or RR)
    ivis = idx + nrparm + iprod*visInfo->incf;
 
     // If all data flagged and not replacing, return
    allBad = (g_data[ivis+2]<=0.0) && (g_data[ivis+5]<=0.0);
    if (visInfo->incf>=12) { // Also XY,YX
       allBad = allBad && (g_data[ivis+8]<=0.0) && (g_data[ivis+11]<=0.0);
    }
    if (allBad && (modelInfo->opType!=GPUBeamOpTypeRepl)) return;
 
    // This one desired?
    if ((kchan<visInfo->chanb) || (kchan>visInfo->chane) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) ||
	((g_data[ivis+2]<=0.0)&&notRepl)) {
	return;
    }

    // frequency scaling factor	  
    ifq = ichan*visInfo->kincf+iIF*visInfo->kincif;
    freqFact = FreqArr[ifq];     // channel frequency scaling factor

    // get scaled u,v,w factors
    u = freqFact * g_data[idx+visInfo->ilocu];
    v = freqFact * g_data[idx+visInfo->ilocv];
    w = freqFact * g_data[idx+visInfo->ilocw];

    // Antenna numbers and types
    GetAnts(ilocb, iloca1, iloca2, ilocsa,  &g_data[idx], &ia1, &ia2); // antenna numbers
    it1 = antInfo->d_AntType[ia1];  // Type of ant1
    it2 = antInfo->d_AntType[ia2];  // Type of ant2
    cudaMatxZero2C(&SumVis);  // Zero model accumulator

    // model = [flux], x,y,z factors, spectral index, coarse channel fluxes
    jchan = antInfo->d_visJChann[ichan];  /* Jones channel for this vis channel */
    iMod = -modelSize;
    for (int i=0; i<nModel; i++) {
        iMod += modelSize;
        if (Model[iMod]==0.0) continue;  // Anything to add?
       // model phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
	// model amplitude
        amp = Model[iMod];	
	// Jones matrices
	// Visibility contribution (IPol for now) - indexed by component number and type
	jindx1 = jchan + i*mChan + it1*mChan*nModel;
	jindx2 = jchan + i*mChan + it2*mChan*nModel;
	if ((jindx1<0) && (jindx2<0)) continue;  // Should not happen - but does???
	// accumulate by Stokes type
	switch (modelInfo->stokType) {
	case 1: // Stokes I
	  cudaJonesCor1_Lin_I (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	case 2: // Stokes Q
	  cudaJonesCor1_Lin_Q (amp, s, c, &Jones[jindx1], &Jones[jindx2], antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 3: // Stokes U
	  cudaJonesCor1_Lin_U (amp, s, c, &Jones[jindx1], &Jones[jindx2],  antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 4: // Stokes V
	  cudaJonesCor1_Lin_V (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	default:
	  printf ("Unsupported Stokes type %d\n",modelInfo->stokType);
	  // should never get here
	}; /* end switch */
	gotSome = 1;
    } // end loop over model comps

    if (!gotSome) return;  // Any valid data?

    // Multiply by Factor
    CUDA_COMPLEX_SMULT(&SumVis.cpx[0], modelInfo->stokFact[0], SumVis.cpx[0]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[1], modelInfo->stokFact[1], SumVis.cpx[1]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[2], modelInfo->stokFact[2], SumVis.cpx[2]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[3], modelInfo->stokFact[3], SumVis.cpx[3]);

    // extract visibility data reordering, if only XX,YY XY,YX not really used
    if (visInfo->incf<12) { // Only XX,YY
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1], 0., 0., 0., 0.,
                              g_data[ivis+3], g_data[ivis+4]);
    } else {
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1],
                              g_data[ivis+6], g_data[ivis+7],
                              g_data[ivis+9], g_data[ivis+10],
                              g_data[ivis+3], g_data[ivis+4]);
    }

    // Correct Data by operation type
    switch (modelInfo->opType) {
    case GPUBeamOpTypeSub: // Subtract
      cudaMatxSub(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
     break;
    case GPUBeamOpTypeDiv: // Divide
      cudaMatxDiv(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
      break;
    case GPUBeamOpTypeRepl: // Replace, even if flagged
      // Replace in Vis data reordering
      if (visInfo->incf<12) { // Only XX,YY
         g_data[ivis+0] = SumVis.cpx[0].real; g_data[ivis+1] = SumVis.cpx[0].imag;
         g_data[ivis+3] = SumVis.cpx[3].real; g_data[ivis+4] = SumVis.cpx[3].imag;
      } else {
        cudaMatxGet2C(&SumVis, &g_data[ivis+0], &g_data[ivis+1],
                               &g_data[ivis+6], &g_data[ivis+7],
                               &g_data[ivis+9], &g_data[ivis+10],
                               &g_data[ivis+3], &g_data[ivis+4]);
      }
      break;
    default:
      printf ("Unsupported op type %d\n",modelInfo->opType);
      // should never get here
    }; /* end switch over operation */

 } // end beamPointKernel

/**
 * Point DFT with parameterized spectrum and beam correction GPU kernal.
 * block = vis, threads in block = data product = channel
 * x dimension = NCHBLK product (channel), y is visibility
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  antInto    Antenna beam information
 * \param  nvis       number of visibilities
 * \param  prtLv      if >=6 print diagnostics
 */
 extern "C"
__global__ void beamPointSpecKernel(float* __restrict__ g_data,
	   GPUBeamModelInfo* modelInfo,  GPUBeamVisInfo* visInfo, GPUBeamAntInfo* antInfo, int nvis, int prtLv)
{
    int lenvis       = visInfo->lenvis;
    int iprod        = threadIdx.x+NCHBLK*blockIdx.x; // product (channel) number
    int idx          = blockIdx.y * lenvis;           // beginning of a visibility 
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->d_freqScale;
    float *freqRat   = visInfo->d_freqRat;
    int modelSize    = modelInfo->size;
    int *specIndex   = modelInfo->d_specIndex;
    int iterm, nterm = modelInfo->nterm;
    int nModel=modelInfo->nmodel;
    float *Model=modelInfo->d_model;
    long ichan, jchan,  kchan, iIF, ivis, iMod, ifq,jindx1, jindx2;
    int ilocb=antInfo->ilocb,  iloca1=antInfo->iloca1, iloca2=antInfo->iloca2;
    int ilocsa=antInfo->ilocsa, ia1, ia2, it1, it2, allBad;
    int mChan  = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    // number of channels in beams
    int notRepl=modelInfo->opType!=GPUBeamOpTypeRepl, gotSome = 0;  // any valid data?
    cudaMatx *Jones = antInfo->d_Jones;
    float arg, amp, s, c; 
    float u, v, w, freqFact, lll, lnspecFreqFact;
    cudaMatx SumVis, visMatx;

    // No more than actual number of channels, or visibilities
    if (iprod>=visInfo->nchan) return;
    if (blockIdx.y>=nvis)      return;

    // get channel, IF from data product
    kchan = (iprod / visInfo->incf)  % visInfo->nchan;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;
    ichan = iprod;

    // real part of vis
    ivis = idx + nrparm + iprod*visInfo->incf;
 
    // If all data flagged and not replacing, return
    allBad = (g_data[ivis+2]<=0.0) && (g_data[ivis+5]<=0.0);
    if (visInfo->incf>=12) { // Also XY,YX
       allBad = allBad && (g_data[ivis+8]<=0.0) && (g_data[ivis+11]<=0.0);
    }
    if (allBad && (modelInfo->opType!=GPUBeamOpTypeRepl)) return;
 
    // This one desired?
    if ((kchan<visInfo->chanb) || (kchan>visInfo->chane) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) ||
	((g_data[ivis+2]<=0.0)&&notRepl)) {
	return;
    }

    // frequency scaling factor	  
    ifq = ichan*visInfo->kincf+iIF*visInfo->kincif;
    //??itab = 7 + 2*specIndex[ifq]; // which coarse channel (subband)?
    freqFact = FreqArr[ifq];     // channel frequency scaling factor
    // log of ratio to reference freq
    lnspecFreqFact = -__logf(freqFact*freqRat[specIndex[ifq]]);

    // get scaled u,v,w factors
    u = freqFact * g_data[idx+visInfo->ilocu];
    v = freqFact * g_data[idx+visInfo->ilocv];
    w = freqFact * g_data[idx+visInfo->ilocw];

    // Antenna numbers and types
    GetAnts(ilocb, iloca1, iloca2, ilocsa,  &g_data[idx], &ia1, &ia2); // antenna numbers (0-rel)
    it1 = antInfo->d_AntType[ia1];  // Type of ant1
    it2 = antInfo->d_AntType[ia2];  // Type of ant2
    cudaMatxZero2C(&SumVis);  // Zero model accumulator

    // model = [flux], x,y,z factors, spectral index, coarse channel fluxes
    jchan = antInfo->d_visJChann[ichan];  /* Jones channel for this vis channel */
    iMod = -modelSize;
    for (int i=0; i<nModel; i++) {
        iMod += modelSize;
        if (Model[iMod]==0.0) continue;  // Anything to add?
        // model phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);

        // amp, Frequency dependent spectral term 
	amp = Model[iMod];
        lll = lnspecFreqFact;
        arg = 0.0;
        for (iterm=0; iterm<nterm; iterm++) {
	  arg += Model[iMod+4+iterm] * lll;
	  lll *= lnspecFreqFact;
	}
	amp *= __expf(-arg);

	// Jones matrices
	// Visibility contribution (IPol for now) - indexed by component number and type
	jindx1 = jchan + i*mChan + it1*mChan*nModel;
	jindx2 = jchan + i*mChan + it2*mChan*nModel;
	if ((jindx1<0) && (jindx2<0)) continue;  // Should not happen - but does???
	// accumulate by Stokes type
	switch (modelInfo->stokType) {
	case 1: // Stokes I
	  cudaJonesCor1_Lin_I (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	case 2: // Stokes Q
	  cudaJonesCor1_Lin_Q (amp, s, c, &Jones[jindx1], &Jones[jindx2], antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 3: // Stokes U
	  cudaJonesCor1_Lin_U (amp, s, c, &Jones[jindx1], &Jones[jindx2],  antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 4: // Stokes V
	  cudaJonesCor1_Lin_V (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	default:
	  printf ("Unsupported Stokes type %d\n",modelInfo->stokType);
	  // should never get here
	}; /* end switch */
	gotSome = 1;
     } // end loop over model comps

    if (!gotSome) return;  // Any valid data?

    // Multiply by Factor
    CUDA_COMPLEX_SMULT(&SumVis.cpx[0], modelInfo->stokFact[0], SumVis.cpx[0]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[1], modelInfo->stokFact[1], SumVis.cpx[1]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[2], modelInfo->stokFact[2], SumVis.cpx[2]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[3], modelInfo->stokFact[3], SumVis.cpx[3]);

    // extract visibility data reordering, if only XX,YY XY,YX not really used
    if (visInfo->incf<12) { // Only XX,YY
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1], 0., 0., 0., 0.,
                              g_data[ivis+3], g_data[ivis+4]);
    } else {
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1],
                              g_data[ivis+6], g_data[ivis+7],
                              g_data[ivis+9], g_data[ivis+10],
                              g_data[ivis+3], g_data[ivis+4]);
    }

    // Correct Data by operation type
    switch (modelInfo->opType) {
    case GPUBeamOpTypeSub: // Subtract
      cudaMatxSub(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
     break;
    case GPUBeamOpTypeDiv: // Divide
      cudaMatxDiv(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
      break;
    case GPUBeamOpTypeRepl: // Replace, even if flagged
      // Replace in Vis data reordering
      if (visInfo->incf<12) { // Only XX,YY
         g_data[ivis+0] = SumVis.cpx[0].real; g_data[ivis+1] = SumVis.cpx[0].imag;
         g_data[ivis+3] = SumVis.cpx[3].real; g_data[ivis+4] = SumVis.cpx[3].imag;
      } else {
        cudaMatxGet2C(&SumVis, &g_data[ivis+0], &g_data[ivis+1],
                               &g_data[ivis+6], &g_data[ivis+7],
                               &g_data[ivis+9], &g_data[ivis+10],
                               &g_data[ivis+3], &g_data[ivis+4]);
      }
      break;
    default:
      printf ("Unsupported op type %d\n",modelInfo->opType);
      // should never get here
    }; /* end switch over operation */

    // diagnostics
    if ((prtLv>=6) && (ia1==0) && (ia2==1) && (ichan==0)) {
      float r1, i1, r2, i2, r3, i3, r4, i4;
      cudaMatxGet2C(&SumVis, &r1, &i1, &r2, &i2, &r3, &i3, &r4, &i4);
      printf("bl %d-%d ch %ld, mod vis %f %f, %f %f, %f %f, %f %f\n",
        ia1+1, ia2+1, ichan, r1, i1, r2, i2, r3, i3, r4, i4);
    }  // end diagnostics

 } // end beamPointSpecKernel

/**
 * Point DFT with tabulated spectrum and beam correction GPU kernal.
 * block = vis, threads in block = data product = channel
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  antInto    Antenna beam information
 * \param  nvis       number of visibilities
 * \param  prtLv      if >=6 print diagnostics
 */
 extern "C"
__global__ void beamPointTSpecKernel(float* __restrict__ g_data,
	   GPUBeamModelInfo* modelInfo,  GPUBeamVisInfo* visInfo, GPUBeamAntInfo* antInfo,
	   int nvis, int prtLv)
{
    int lenvis       = visInfo->lenvis;
    int iprod        = threadIdx.x+NCHBLK*blockIdx.x; // product (channel) number
    int idx          = blockIdx.y * lenvis;           // beginning of a visibility 
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->d_freqScale;
    float *freqRat   = visInfo->d_freqRat;
    int modelSize    = modelInfo->size;
    int *specIndex   = modelInfo->d_specIndex;
    int nModel       = modelInfo->nmodel;
    float *Model     = modelInfo->d_model;
    int mChan        = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    // number of channels in beams
    long ichan, jchan, kchan, iIF, ivis, itab, iMod, ifq, nchpIF,jindx1, jindx2;
    int ilocb=antInfo->ilocb,  iloca1=antInfo->iloca1, iloca2=antInfo->iloca2;
    int ilocsa=antInfo->ilocsa, ia1, ia2, it1, it2, allBad;
    int notRepl=modelInfo->opType!=GPUBeamOpTypeRepl, gotSome = 0;  // any valid data?
    cudaMatx *Jones = antInfo->d_Jones;
    float arg, amp, s, c; 
    float u, v, w, freqFact, lnspecFreqFact;
    cudaMatx SumVis,visMatx;

    // No more than actual number of channels, or visibilities
    if (iprod>=visInfo->nchan) return;
    if (blockIdx.y>=nvis)      return;

    // real part of vis
    ivis = idx + nrparm + iprod*visInfo->incf;

    // If all data flagged and not replacing, return
    allBad = (g_data[ivis+2]<=0.0) && (g_data[ivis+5]<=0.0);
    if (visInfo->incf>=12) { // Also XY,YX
       allBad = allBad && (g_data[ivis+8]<=0.0) && (g_data[ivis+11]<=0.0);
    }
    if (allBad && (modelInfo->opType!=GPUBeamOpTypeRepl)) return;
 
    // get channel, IF from data product
    nchpIF = visInfo->nchan/visInfo->nIF;
    kchan = iprod  % nchpIF;  // ch in IF
    iIF   = iprod / nchpIF;
    ichan = iprod;  // total channel number

    // This one desired?
    if (((kchan+1)<visInfo->chanb) || (kchan>visInfo->chane) ||
 	((iIF+1)<visInfo->IFb)     || (iIF>visInfo->IFe) ||
	((g_data[ivis+2]<=0.0)&&notRepl)) {
 	return;
     }

    // frequency scaling factor	  
    ifq = ichan;
    itab = 7 + 2*specIndex[ifq]; // which coarse channel (subband)?
    freqFact = FreqArr[ifq];     // channel frequency scaling factor
    // log of ratio channel frequency to subband frequency
    lnspecFreqFact = __logf(freqFact*freqRat[specIndex[ifq]]); 

    // get scaled u,v,w factors
    u = freqFact * g_data[idx+visInfo->ilocu];
    v = freqFact * g_data[idx+visInfo->ilocv];
    w = freqFact * g_data[idx+visInfo->ilocw];

    // Antenna numbers  and types
    GetAnts(ilocb, iloca1, iloca2, ilocsa,  &g_data[idx], &ia1, &ia2); // antenna numbers
    it1 = antInfo->d_AntType[ia1];  // Type of ant1
    it2 = antInfo->d_AntType[ia2];  // Type of ant2
    cudaMatxZero2C(&SumVis);  // Zero model accu,ulator

    // model = [flux], x,y,z factors, spectral index, coarse channel fluxes
    jchan = antInfo->d_visJChann[ichan];  /* Jones channel for this vis channel */
    iMod = -modelSize;
    for (int i=0; i<nModel; i++) {
        iMod += modelSize;
        if (Model[iMod+itab]==0.0) continue;  // Anything to add?
        // model phase
 	arg = u*Model[iMod+4] + v*Model[iMod+5] + w*Model[iMod+6];//Fixed
	__sincosf(arg, &s, &c);
        // Amplitude, include spectral index
	arg = lnspecFreqFact*Model[iMod+itab+1];
	amp =  Model[iMod+itab] * __expf(arg);
	//amp *= -6.28;  //DEBUG
	// Get Jones matrix indices
	jindx1 = jchan + i*mChan + it1*mChan*nModel;
	jindx2 = jchan + i*mChan + it2*mChan*nModel;
	if ((jindx1<0) && (jindx2<0)) continue;  // Should not happen - but does???

	// accumulate by Stokes type
	switch (modelInfo->stokType) {
	case 1: // Stokes I
	  cudaJonesCor1_Lin_I (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	case 2: // Stokes Q
	  cudaJonesCor1_Lin_Q (amp, s, c, &Jones[jindx1], &Jones[jindx2], antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 3: // Stokes U
	  cudaJonesCor1_Lin_U (amp, s, c, &Jones[jindx1], &Jones[jindx2],  antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 4: // Stokes V
	  cudaJonesCor1_Lin_V (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	default:
	  printf ("Unsupported Stokes type %d\n",modelInfo->stokType);
	  // should never get here
	}; /* end switch */
	gotSome = 1;
    } // end loop over model comps

    if (!gotSome) return;  // Any valid data?

    // Multiply by Factor
    CUDA_COMPLEX_SMULT(&SumVis.cpx[0], modelInfo->stokFact[0], SumVis.cpx[0]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[1], modelInfo->stokFact[1], SumVis.cpx[1]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[2], modelInfo->stokFact[2], SumVis.cpx[2]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[3], modelInfo->stokFact[3], SumVis.cpx[3]);

    // extract visibility data reordering, if only XX,YY XY,YX not really used
    if (visInfo->incf<12) { // Only XX,YY
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1], 0., 0., 0., 0.,
                              g_data[ivis+3], g_data[ivis+4]);
    } else {
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1],
                              g_data[ivis+6], g_data[ivis+7],
                              g_data[ivis+9], g_data[ivis+10],
                              g_data[ivis+3], g_data[ivis+4]);
    }

    // Correct Data by operation type
    switch (modelInfo->opType) {
    case GPUBeamOpTypeSub: // Subtract
      cudaMatxSub(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
     break;
    case GPUBeamOpTypeDiv: // Divide
      cudaMatxDiv(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
      break;
    case GPUBeamOpTypeRepl: // Replace, even if flagged
      // Replace in Vis data reordering
      if (visInfo->incf<12) { // Only XX,YY
         g_data[ivis+0] = SumVis.cpx[0].real; g_data[ivis+1] = SumVis.cpx[0].imag;
         g_data[ivis+3] = SumVis.cpx[3].real; g_data[ivis+4] = SumVis.cpx[3].imag;
      } else {
        cudaMatxGet2C(&SumVis, &g_data[ivis+0], &g_data[ivis+1],
                               &g_data[ivis+6], &g_data[ivis+7],
                               &g_data[ivis+9], &g_data[ivis+10],
                               &g_data[ivis+3], &g_data[ivis+4]);
      }
      break;
    default:
      printf ("Unsupported op type %d\n",modelInfo->opType);
      // should never get here
    }; /* end switch over operation */


} // end beamPointTSpecKernel

/************************************** Gaussian model ******************************************************/
/**
 * Gauss single channel DFT with correction GPU kernal.
 * block = vis, threads in block = data product = channel
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  antInto    Antenna beam information
 * \param  nvis       number of visibilities
 */
 extern "C"
__global__ void beamGaussKernel(float* __restrict__ g_data,
	   GPUBeamModelInfo* modelInfo,  GPUBeamVisInfo* visInfo, GPUBeamAntInfo* antInfo, int nvis)
{
    int lenvis       = visInfo->lenvis;
    int iprod        = threadIdx.x+NCHBLK*blockIdx.x; // product (channel) number
    int idx          = blockIdx.y * lenvis;           // beginning of a visibility 
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->d_freqScale;
    int modelSize    = modelInfo->size;
    int nModel       = modelInfo->nmodel;
    float *Model     = modelInfo->d_model;
    int mChan       = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    // number of channels in beams
    long ichan, jchan,  kchan, iIF, ivis, iMod, ifq, jindx1, jindx2;
    int ilocb=antInfo->ilocb,  iloca1=antInfo->iloca1, iloca2=antInfo->iloca2;
    int ilocsa=antInfo->ilocsa, ia1, ia2, it1, it2, allBad;
    int notRepl=modelInfo->opType!=GPUBeamOpTypeRepl, gotSome = 0;  // any valid data?
    cudaMatx *Jones = antInfo->d_Jones;
    float arg, amp, s, c; 
    float u, v, w, freqFact;
    cudaMatx SumVis, visMatx;

    // No more than actual number of channels
    if (iprod>=visInfo->nchan) return;
    if (blockIdx.y>=nvis)      return;

    // get channel, IF from data product
    kchan = (iprod / visInfo->incf)  % visInfo->nchan;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;
    ichan = iprod;

    // real part of vis
    ivis = idx + nrparm + iprod*visInfo->incf;
 
    // If all data flagged and not replacing, return
    allBad = (g_data[ivis+2]<=0.0) && (g_data[ivis+5]<=0.0);
    if (visInfo->incf>=12) { // Also XY,YX
       allBad = allBad && (g_data[ivis+8]<=0.0) && (g_data[ivis+11]<=0.0);
    }
    if (allBad && (modelInfo->opType!=GPUBeamOpTypeRepl)) return;
 
    // This one desired?
    if ((kchan<visInfo->chanb) || (kchan>visInfo->chane) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) ||
	((g_data[ivis+2]<=0.0)&&notRepl)) {
	return;
    }

    // frequency scaling factor	  
    ifq = ichan*visInfo->kincf+iIF*visInfo->kincif;
    freqFact = FreqArr[ifq];     // channel frequency scaling factor

    // get scaled u,v,w factors
    u = freqFact * g_data[idx+visInfo->ilocu];
    v = freqFact * g_data[idx+visInfo->ilocv];
    w = freqFact * g_data[idx+visInfo->ilocw];

    // Antenna numbers and types
    GetAnts(ilocb, iloca1, iloca2, ilocsa,  &g_data[idx], &ia1, &ia2); // antenna numbers
    it1 = antInfo->d_AntType[ia1];  // Type of ant1
    it2 = antInfo->d_AntType[ia2];  // Type of ant2
    cudaMatxZero2C(&SumVis);  // Zero model accumulator

    // model = [flux], x,y,z factors, spectral index, coarse channel fluxes
    jchan = antInfo->d_visJChann[ichan];  /* Jones channel for this vis channel */
    iMod = -modelSize;
    for (int i=0; i<nModel; i++) {
        iMod += modelSize;
        if (Model[iMod]==0.0) continue;  // Anything to add?
       // model phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);

        // model amplitude
        amp = Model[iMod];	
        arg = u*u*Model[iMod+4] + v*v*Model[iMod+5] + u*v*Model[iMod+6];  // Gaussian term
	amp *= __expf(arg);

        // Jones matrices
	// Get Jones matrix indices
	jindx1 = jchan + i*mChan + it1*mChan*nModel;
	jindx2 = jchan + i*mChan + it2*mChan*nModel;
	if ((jindx1<0) && (jindx2<0)) continue;  // Should not happen - but does???
	// accumulate by Stokes type
	switch (modelInfo->stokType) {
	case 1: // Stokes I
	  cudaJonesCor1_Lin_I (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	case 2: // Stokes Q
	  cudaJonesCor1_Lin_Q (amp, s, c, &Jones[jindx1], &Jones[jindx2], antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 3: // Stokes U
	  cudaJonesCor1_Lin_U (amp, s, c, &Jones[jindx1], &Jones[jindx2],  antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 4: // Stokes V
	  cudaJonesCor1_Lin_V (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	default:
	  printf ("Unsupported Stokes type %d\n",modelInfo->stokType);
	  // should never get here
	}; /* end switch */
        gotSome = 1;
    } // end loop over model comps

    if (!gotSome) return;  // Any valid data?

    // Multiply by Factor
    CUDA_COMPLEX_SMULT(&SumVis.cpx[0], modelInfo->stokFact[0], SumVis.cpx[0]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[1], modelInfo->stokFact[1], SumVis.cpx[1]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[2], modelInfo->stokFact[2], SumVis.cpx[2]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[3], modelInfo->stokFact[3], SumVis.cpx[3]);

    // extract visibility data reordering, if only XX,YY XY,YX not really used
    if (visInfo->incf<12) { // Only XX,YY
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1], 0., 0., 0., 0.,
                              g_data[ivis+3], g_data[ivis+4]);
    } else {
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1],
                              g_data[ivis+6], g_data[ivis+7],
                              g_data[ivis+9], g_data[ivis+10],
                              g_data[ivis+3], g_data[ivis+4]);
    }

    // Correct Data by operation type
    switch (modelInfo->opType) {
    case GPUBeamOpTypeSub: // Subtract
      cudaMatxSub(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
     break;
    case GPUBeamOpTypeDiv: // Divide
      cudaMatxDiv(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
      break;
    case GPUBeamOpTypeRepl: // Replace, even if flagged
      // Replace in Vis data reordering
      if (visInfo->incf<12) { // Only XX,YY
         g_data[ivis+0] = SumVis.cpx[0].real; g_data[ivis+1] = SumVis.cpx[0].imag;
         g_data[ivis+3] = SumVis.cpx[3].real; g_data[ivis+4] = SumVis.cpx[3].imag;
      } else {
        cudaMatxGet2C(&SumVis, &g_data[ivis+0], &g_data[ivis+1],
                               &g_data[ivis+6], &g_data[ivis+7],
                               &g_data[ivis+9], &g_data[ivis+10],
                               &g_data[ivis+3], &g_data[ivis+4]);
      }
      break;
    default:
      printf ("Unsupported op type %d\n",modelInfo->opType);
      // should never get here
    }; /* end switch over operation */
 } // end beamGaussKernel

/**
 * Gauss DFT with parameterized spectrum and beam correction GPU kernal.
 * block = vis, threads in block = data product = channel
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  antInto    Antenna beam information
 * \param  nvis       number of visibilities
 */
 extern "C"
__global__ void beamGaussSpecKernel(float* __restrict__ g_data,
	   GPUBeamModelInfo* modelInfo,  GPUBeamVisInfo* visInfo, GPUBeamAntInfo* antInfo, int nvis)
{
    int lenvis       = visInfo->lenvis;
    int iprod        = threadIdx.x+NCHBLK*blockIdx.x; // product (channel) number
    int idx          = blockIdx.y * lenvis;           // beginning of a visibility 
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->d_freqScale;
    float *freqRat   = visInfo->d_freqRat;
    int modelSize    = modelInfo->size;
    int *specIndex   = modelInfo->d_specIndex;
    int iterm, nterm = modelInfo->nterm;
    int nModel       = modelInfo->nmodel;
    float *Model     = modelInfo->d_model;
    int mChan        = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    // number of channels in beams
    long ichan, jchan,  kchan, iIF, ivis, iMod, ifq, jindx1, jindx2;
    int ilocb=antInfo->ilocb,  iloca1=antInfo->iloca1, iloca2=antInfo->iloca2;
    int ilocsa=antInfo->ilocsa, ia1, ia2, it1, it2, allBad;
    int notRepl=modelInfo->opType!=GPUBeamOpTypeRepl, gotSome = 0;  // any valid data?
    cudaMatx *Jones = antInfo->d_Jones;
    float arg, amp, s, c; 
    float u, v, w, freqFact, lll, lnspecFreqFact;
    cudaMatx SumVis, visMatx;

    // No more than actual number of channels, or visibilities
    if (iprod>=visInfo->nchan) return;
    if (blockIdx.y>=nvis)      return;

    // get channel, IF from data product
    kchan = (iprod / visInfo->incf)  % visInfo->nchan;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;
    ichan = iprod;

    // real part of vis
    ivis = idx + nrparm + iprod*visInfo->incf;
 
    // If all data flagged and not replacing, return
    allBad = (g_data[ivis+2]<=0.0) && (g_data[ivis+5]<=0.0);
    if (visInfo->incf>=12) { // Also XY,YX
       allBad = allBad && (g_data[ivis+8]<=0.0) && (g_data[ivis+11]<=0.0);
    }
    if (allBad && (modelInfo->opType!=GPUBeamOpTypeRepl)) return;
 
    // This one desired?
    if ((kchan<visInfo->chanb) || (kchan>visInfo->chane) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) ||
	((g_data[ivis+2]<=0.0)&&notRepl)) {
	return;
    }

    // frequency scaling factor	  
    ifq = ichan*visInfo->kincf+iIF*visInfo->kincif;
    //??itab = 8 + 2*specIndex[ifq]; // which coarse channel (subband)?
    freqFact = FreqArr[ifq];     // channel frequency scaling factor
    // log of ratio to reference freq
    lnspecFreqFact = -__logf(freqFact*freqRat[specIndex[ifq]]);

    // get scaled u,v,w factors
    u = freqFact * g_data[idx+visInfo->ilocu];
    v = freqFact * g_data[idx+visInfo->ilocv];
    w = freqFact * g_data[idx+visInfo->ilocw];

    // Antenna numbers and types
    GetAnts(ilocb, iloca1, iloca2, ilocsa,  &g_data[idx], &ia1, &ia2); // antenna numbers
    it1 = antInfo->d_AntType[ia1];  // Type of ant1
    it2 = antInfo->d_AntType[ia2];  // Type of ant2
    cudaMatxZero2C(&SumVis);  // Zero model accumulator

    // model = [flux], x,y,z factors, spectral index, coarse channel fluxes
    jchan = antInfo->d_visJChann[ichan];  /* Jones channel for this vis channel */
    iMod = -modelSize;
    for (int i=0; i<nModel; i++) {
        iMod += modelSize;
        if (Model[iMod]==0.0) continue;  // Anything to add?
       // model phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);

        // amp, Frequency dependent spectral term 
	amp = Model[iMod];
        lll = lnspecFreqFact;
        arg = 0.0;
        for (iterm=0; iterm<nterm; iterm++) {
	  arg += Model[iMod+4+iterm] * lll;
	  lll *= lnspecFreqFact;
	}
	amp *= __expf(-arg);
        arg  = u*u*Model[iMod+4] + v*v*Model[iMod+5] + u*v*Model[iMod+6];  // Gaussian term
	amp *= __expf(arg);

	// Jones matrices
	// Visibility contribution (IPol for now) - indexed by component number and type
	jindx1 = jchan + i*mChan + it1*mChan*nModel;
	jindx2 = jchan + i*mChan + it2*mChan*nModel;
	if ((jindx1<0) && (jindx2<0)) continue;  // Should not happen - but does???
	// accumulate by Stokes type
	switch (modelInfo->stokType) {
	case 1: // Stokes I
	  cudaJonesCor1_Lin_I (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	case 2: // Stokes Q
	  cudaJonesCor1_Lin_Q (amp, s, c, &Jones[jindx1], &Jones[jindx2], antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 3: // Stokes U
	  cudaJonesCor1_Lin_U (amp, s, c, &Jones[jindx1], &Jones[jindx2],  antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 4: // Stokes V
	  cudaJonesCor1_Lin_V (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	default:
	  printf ("Unsupported Stokes type %d\n",modelInfo->stokType);
	  // should never get here
	}; /* end switch */
	gotSome = 1;
    } // end loop over model comps

    if (!gotSome) return;  // Any valid data?

    // Multiply by Factor
    CUDA_COMPLEX_SMULT(&SumVis.cpx[0], modelInfo->stokFact[0], SumVis.cpx[0]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[1], modelInfo->stokFact[1], SumVis.cpx[1]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[2], modelInfo->stokFact[2], SumVis.cpx[2]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[3], modelInfo->stokFact[3], SumVis.cpx[3]);

    // extract visibility data reordering, if only XX,YY XY,YX not really used
    if (visInfo->incf<12) { // Only XX,YY
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1], 0., 0., 0., 0.,
                              g_data[ivis+3], g_data[ivis+4]);
    } else {
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1],
                              g_data[ivis+6], g_data[ivis+7],
                              g_data[ivis+9], g_data[ivis+10],
                              g_data[ivis+3], g_data[ivis+4]);
    }

    // Correct Data by operation type
    switch (modelInfo->opType) {
    case GPUBeamOpTypeSub: // Subtract
      cudaMatxSub(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
     break;
    case GPUBeamOpTypeDiv: // Divide
      cudaMatxDiv(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
      break;
    case GPUBeamOpTypeRepl: // Replace, even if flagged
      // Replace in Vis data reordering
      if (visInfo->incf<12) { // Only XX,YY
         g_data[ivis+0] = SumVis.cpx[0].real; g_data[ivis+1] = SumVis.cpx[0].imag;
         g_data[ivis+3] = SumVis.cpx[3].real; g_data[ivis+4] = SumVis.cpx[3].imag;
      } else {
        cudaMatxGet2C(&SumVis, &g_data[ivis+0], &g_data[ivis+1],
                               &g_data[ivis+6], &g_data[ivis+7],
                               &g_data[ivis+9], &g_data[ivis+10],
                               &g_data[ivis+3], &g_data[ivis+4]);
      }
      break;
    default:
      printf ("Unsupported op type %d\n",modelInfo->opType);
      // should never get here
    }; /* end switch over operation */
 } // end beamGaussSpecKernel

/**
 * Gauss DFT with tabulated spectrum and beam correction GPU kernal.
 * block = vis, threads in block = data product = channel
 * does full model for one data product
 * \param  g_data     vis data
 * \param  modelInfo  model information
 * \param  visInfo    visibility information
 * \param  antInto    Antenna beam information
 * \param  nvis       number of visibilities
 */
 extern "C"
__global__ void beamGaussTSpecKernel(float* __restrict__ g_data,
	   GPUBeamModelInfo* modelInfo,  GPUBeamVisInfo* visInfo, GPUBeamAntInfo* antInfo, int nvis)
{
    int lenvis       = visInfo->lenvis;
    int iprod        = threadIdx.x+NCHBLK*blockIdx.x; // product (channel) number
    int idx          = blockIdx.y * lenvis;           // beginning of a visibility 
    int nrparm       = visInfo->nrparm;
    float *FreqArr   = visInfo->d_freqScale;
    float *freqRat   = visInfo->d_freqRat;
    int modelSize    = modelInfo->size;
    int *specIndex   = modelInfo->d_specIndex;
    int nModel       = modelInfo->nmodel;
    float *Model     = modelInfo->d_model;
    int mChan        = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    // number of channels in beams
    long ichan, jchan,  kchan, iIF, ivis, itab, iMod, ifq, jindx1, jindx2;
    int ilocb=antInfo->ilocb,  iloca1=antInfo->iloca1, iloca2=antInfo->iloca2;
    int ilocsa=antInfo->ilocsa, ia1, ia2, it1, it2, allBad;
    int notRepl=modelInfo->opType!=GPUBeamOpTypeRepl, gotSome = 0;  // any valid data?
    cudaMatx *Jones = antInfo->d_Jones;
    float arg, amp, s, c; 
    float u, v, w, freqFact, lnspecFreqFact;
    cudaMatx SumVis, visMatx;

    // No more than actual number of channels, or visibilities
    if (iprod>=visInfo->nchan) return;
    if (blockIdx.y>=nvis)      return;

    // get channel, IF from data product
    kchan = (iprod / visInfo->incf)  % visInfo->nchan;
    iIF   = (iprod / visInfo->incif) % visInfo->nIF;
    ichan = iprod;

    // real part of vis
    ivis = idx + nrparm + iprod*visInfo->incf;
 
    // If all data flagged and not replacing, return
    allBad = (g_data[ivis+2]<=0.0) && (g_data[ivis+5]<=0.0);
    if (visInfo->incf>=12) { // Also XY,YX
       allBad = allBad && (g_data[ivis+8]<=0.0) && (g_data[ivis+11]<=0.0);
    }
    if (allBad && (modelInfo->opType!=GPUBeamOpTypeRepl)) return;
 
    // This one desired?
    if ((kchan<visInfo->chanb) || (kchan>visInfo->chane) ||
	(iIF<visInfo->IFb)     || (iIF>visInfo->IFe) ||
	((g_data[ivis+2]<=0.0)&&notRepl)) {
	return;
    }

    // frequency scaling factor	  
    ifq = ichan*visInfo->kincf+iIF*visInfo->kincif;
    itab = 10 + 2*specIndex[ifq]; // which coarse channel (subband)?
    freqFact = FreqArr[ifq];     // channel frequency scaling factor
    // log of ratio to reference freq
    lnspecFreqFact = __logf(freqFact*freqRat[specIndex[ifq]]);

    // get scaled u,v,w factors
    u = freqFact * g_data[idx+visInfo->ilocu];
    v = freqFact * g_data[idx+visInfo->ilocv];
    w = freqFact * g_data[idx+visInfo->ilocw];

    // Antenna numbers and types
    GetAnts(ilocb, iloca1, iloca2, ilocsa,  &g_data[idx], &ia1, &ia2); // antenna numbers
    it1 = antInfo->d_AntType[ia1];  // Type of ant1
    it2 = antInfo->d_AntType[ia2];  // Type of ant2
    cudaMatxZero2C(&SumVis);  // Zero model accumulator

    // model = [flux], x,y,z factors, spectral index, coarse channel fluxes
    jchan = antInfo->d_visJChann[ichan];  /* Jones channel for this vis channel */
    iMod = -modelSize;
    for (int i=0; i<nModel; i++) {
        iMod += modelSize;
        if (Model[iMod+itab]==0.0) continue;  // Anything to add?
        // model phase
 	arg = u*Model[iMod+1] + v*Model[iMod+2] + w*Model[iMod+3];
	__sincosf(arg, &s, &c);
	
        // Amplitude, include spectral index
	arg = lnspecFreqFact*Model[iMod+itab+1];
	amp =  Model[iMod+itab] * __expf(arg);
        arg = u*u*Model[iMod+4] + v*v*Model[iMod+5] + u*v*Model[iMod+6];  // Gaussian term
	amp  =__expf(arg);
	
	// Jones matrices
	// Visibility contribution (IPol for now) - indexed by component number and type
	jindx1 = jchan + i*mChan + it1*mChan*nModel;
	jindx2 = jchan + i*mChan + it2*mChan*nModel;
	if ((jindx1<0) && (jindx2<0)) continue;  // Should not happen - but does???
	// accumulate by Stokes type
	switch (modelInfo->stokType) {
	case 1: // Stokes I
	  cudaJonesCor1_Lin_I (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	case 2: // Stokes Q
	  cudaJonesCor1_Lin_Q (amp, s, c, &Jones[jindx1], &Jones[jindx2], antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 3: // Stokes U
	  cudaJonesCor1_Lin_U (amp, s, c, &Jones[jindx1], &Jones[jindx2],  antInfo->c2p, antInfo->s2p, &SumVis);
	  break;
	case 4: // Stokes V
	  cudaJonesCor1_Lin_V (amp, s, c, &Jones[jindx1], &Jones[jindx2], 0.0, 0.0, &SumVis);
	  break;
	default:
	  printf ("Unsupported Stokes type %d\n",modelInfo->stokType);
	  // should never get here
	}; /* end switch */
        gotSome = 1;
    } // end loop over model comps

    if (!gotSome) return;  // Any valid data?

    // Multiply by Factor
    CUDA_COMPLEX_SMULT(&SumVis.cpx[0], modelInfo->stokFact[0], SumVis.cpx[0]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[1], modelInfo->stokFact[1], SumVis.cpx[1]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[2], modelInfo->stokFact[2], SumVis.cpx[2]);
    CUDA_COMPLEX_SMULT(&SumVis.cpx[3], modelInfo->stokFact[3], SumVis.cpx[3]);

    // extract visibility data reordering, if only XX,YY XY,YX not really used
    if (visInfo->incf<12) { // Only XX,YY
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1], 0., 0., 0., 0.,
                              g_data[ivis+3], g_data[ivis+4]);
    } else {
      cudaMatxSet2C(&visMatx, g_data[ivis+0], g_data[ivis+1],
                              g_data[ivis+6], g_data[ivis+7],
                              g_data[ivis+9], g_data[ivis+10],
                              g_data[ivis+3], g_data[ivis+4]);
    }

    // Correct Data by operation type
    switch (modelInfo->opType) {
    case GPUBeamOpTypeSub: // Subtract
      cudaMatxSub(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
     break;
    case GPUBeamOpTypeDiv: // Divide
      cudaMatxDiv(&visMatx, &SumVis, &visMatx);
      // Replace in Vis data reordering, only if not flagged
      if (g_data[ivis+2]>0.0) {g_data[ivis+0] = visMatx.cpx[0].real; g_data[ivis+1] = visMatx.cpx[0].imag;}
      if (g_data[ivis+5]>0.0) {g_data[ivis+3] = visMatx.cpx[3].real; g_data[ivis+4] = visMatx.cpx[3].imag;}
      if (visInfo->incf>=12) { // Also XY,YX
        if (g_data[ivis+8]>0.0)  {g_data[ivis+6] = visMatx.cpx[1].real; g_data[ivis+7]  = visMatx.cpx[1].imag;}
        if (g_data[ivis+11]>0.0) {g_data[ivis+9] = visMatx.cpx[2].real; g_data[ivis+10] = visMatx.cpx[2].imag;}
      }
      break;
    case GPUBeamOpTypeRepl: // Replace, even if flagged
      // Replace in Vis data reordering
      if (visInfo->incf<12) { // Only XX,YY
         g_data[ivis+0] = SumVis.cpx[0].real; g_data[ivis+1] = SumVis.cpx[0].imag;
         g_data[ivis+3] = SumVis.cpx[3].real; g_data[ivis+4] = SumVis.cpx[3].imag;
      } else {
        cudaMatxGet2C(&SumVis, &g_data[ivis+0], &g_data[ivis+1],
                               &g_data[ivis+6], &g_data[ivis+7],
                               &g_data[ivis+9], &g_data[ivis+10],
                               &g_data[ivis+3], &g_data[ivis+4]);
      }
      break;
    default:
      printf ("Unsupported op type %d\n",modelInfo->opType);
      // should never get here
    }; /* end switch over operation */
 } // end beamGaussTSpecKernel

/**
 * Compose Jones matrices by interpolating a pixel in beam images
 * per antenna, component, channel slowest to fastest
 * where "channel" is the channelization of beam
 * This routine loops over channel for an antenna and model component
 * Output matrices are divided by sqrt(symmetric_beam_gain) for reference ant. size
* Uses unit matrices where "perfect" beam gain < 0.10 or high beam Jones matrices
 *
 * \param  antInfo    Antenna information
 * \param  modelInfo  CLEAN model information
 * \param  parAng     Parallactic Angle (radians)
 * \param  prtLv      if >=6 print diagnostic values for ant 2, model 2.
 */
 extern "C" __global__ 
void GetJonesKernel (GPUBeamAntInfo* __restrict__ antInfo, 
                     const GPUBeamModelInfo* __restrict__ modelInfo, 
		     float parAng, int prtLv)
{
    long imod   = threadIdx.x + blockIdx.x*blockDim.x;  // model number
    long itype  = threadIdx.y + blockIdx.y*blockDim.y;  // type number
    long nchan  = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    //number of channels in beams
    int  nType  = antInfo->numAntType;     // Number of antenna types
    long ichan;
    long indx;
    float refSizeFact = antInfo->refSizeFact;  // Antenna size factor for PB gain
    float P1r, P1i, P2r, P2i, X1r, X1i, X2r, X2i;
    float *Model, pixel[3];
    float maxIGain=10.0;  // maximum inverse "perfect" gain
    int lenmodel, nModel = modelInfo->nmodel;
    cudaMatx *Jones;
    double pos[2];
    float ang, iPBCor;

    // test if values in range
    if ((imod>=nModel) || (itype>=nType)) return;

    Model    = modelInfo->d_model;
    lenmodel = modelInfo->size;

    // Calculate pixel
    pos[0] = -(double)Model[imod*lenmodel+2];
    pos[1] = +(double)Model[imod*lenmodel+3];
    indx = itype*8;  // a set of  8 for each type 
    if (CUDAPositionXYpix(pos, (CUDAImageDesc*)antInfo->d_BeamInterp[indx]->d_myDesc, parAng, pixel)!=0) {
      //  Pixel not in beam, use Unit matrix
      for (ichan=0; ichan<nchan; ichan++) {
        indx = ichan + imod*nchan + itype*(long)nModel*nchan;
        Jones = &antInfo->d_Jones[indx];
        cudaMatxSet2C (Jones, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	} // end channel loop
      return;
    }

    // Interpolate, looping over channels in beams
    ang = __fsqrt_rn((float)(pos[0]*pos[0]+pos[1]*pos[1])); /* Offset angle of comp. */
    for (ichan=0; ichan<nchan; ichan++) {
      // Get reference beam gain = sqrt(1/PB gain)
      iPBCor = InvBeam (refSizeFact, antInfo->d_BeamFreq[ichan], ang);
      pixel[2] = ichan+1;
      indx = itype*8;  // a set of  8 for each type 
      P1r = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+0], pixel);
      P2r = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+2], pixel);
      // Sanity check, use unit matrix if bad, parallel reals >1.3
      if ((iPBCor<maxIGain) && (fabs(P1r)<=1.3) && (fabs(P1r)<=1.3)) {
        P1i = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+1], pixel);
        P2i = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+3], pixel);
        X1r = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+4], pixel);
        X1i = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+5], pixel);
        X2r = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+6], pixel);
        X2i = iPBCor*CUDABeamInterpPixel (antInfo->d_BeamInterp[indx+7], pixel);
      } else {  // Use unit matrix
        P1r = P2r = 1.0;  P1i = P2i = 0.0; 
        X1r = X2r = X1i = X2i = 0.0; 
      }
      indx = ichan + imod*nchan + itype*(long)nchan*nModel;
      Jones = &antInfo->d_Jones[indx];
// debug printing
    if (prtLv>=6 && itype==0 && imod==0) {
        printf("in GetJonesKernel chan %ld type %ld imod %ld pos %lg %lg pixel %f %f %f gain %f freq %f ang %f parAng %f Jones %f %f %f %f %f %f %f %f\n",
	 ichan, itype, imod, pos[0], pos[1], pixel[0], pixel[1],pixel[2],1./iPBCor, antInfo->d_BeamFreq[ichan], ang, parAng,P1r, P1i, X1r, X1i, X2r, X2i, P2r, P2i);
        printf("     Raw Jones  %f %f %f %f %f %f %f %f\n",
	P1r/iPBCor, P1i/iPBCor, X1r/iPBCor, X1i/iPBCor, X2r/iPBCor, X2i/iPBCor, P2r/iPBCor, P2i/iPBCor);
    }
      cudaMatxSet2C (Jones, P1r, P1i, X1r, X1i, X2r, X2i, P2r, P2i);
    } // end channel loop
} /* End GetJonesKernel  */

/**
 * Set Unit Jones matrices
 * per antenna, component, channel slowest to fastest
 * where "channel* is the channelization of beam
 * This routine loop over channel for an antenna and model component
 *
 * \param  antInfo    Antenna information
 * \param  modelInfo  CLEAN model information
 */
 extern "C" __global__ 
void UnitJonesKernel (GPUBeamAntInfo* __restrict__ antInfo,
                     const GPUBeamModelInfo* __restrict__ modelInfo)
{
    long imod   = threadIdx.x + blockIdx.x*blockDim.x;  // model number
    long itype  = threadIdx.y + blockIdx.y*blockDim.y;  // type number
    long nchan  = antInfo->d_BeamInterp[0]->d_myDesc->inaxes[2];    //number of channels in beams
    long ichan;
    int nType  = antInfo->numAntType;     // Number of antenna types
    int nModel = modelInfo->nmodel;       // Number of models
    int indx;
    cudaMatx *Jones;

    // test if values in range
    if ((imod>=nModel) || (itype>=nType)) return;

    indx = itype*8;  // a set of  8 for each type 
    //  Unit matrix
    for (ichan=0; ichan<nchan; ichan++) {
      indx = ichan + imod*nchan + itype*nModel*nchan;
      Jones = &antInfo->d_Jones[indx];
      cudaMatxSet2C (Jones, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
      } // end channel loop
      return;
} /* End UnitJonesKernel  */
/**
 * Update parallactic angle info on antInfo
 * a single execution should do it
 *
 * \param  antInfo    Antenna information
 * \param  parAng     Parallactic angle (rad)
 */
 extern "C" __global__ 
void ParAngKernel (GPUBeamAntInfo* __restrict__ antInfo, float parAng)
{
    // Save parallactic angle sin, cos
    __sincosf(2*parAng, &antInfo->s2p, &antInfo->c2p);
    //DEBUG - NO__sincosf(-2*parAng, &antInfo->s2p, &antInfo->c2p);

      return;
} /* End ParAngKernel  */

static void BeamProcessWithStreams(int streams_used, 
            int doJones, int doUnit, float parAng, long lvisRange[2],
            GPUBeamModelInfo* h_modelInfo, GPUBeamModelInfo* d_modelInfo, 
            GPUBeamVisInfo* h_visInfo, GPUBeamVisInfo* d_visInfo, 
            GPUBeamAntInfo* h_antInfo, GPUBeamAntInfo* d_antInfo, 
            cudaStream_t* stream, cudaEvent_t* cycleDone,
            float *h_data_source, float *h_data_sink, 
            float *d_data[], int prtLv);
#endif /* HAVE_GPU */

/*----------------------Public functions---------------------------*/
/**
 * Initialize an ObitCUDABeamModel 
 * Currently only DFT point supported
 * \param gpuInfo    processing info
 */
extern "C"
void ObitCUDABeamModelDFTInit (GPUBeamInfo *gpuInfo)
{
  /* Not much to do here */
} /* end  ObitCUDABeamModelDFTInit */


/**
 * Setup model for ObitCUDABeamModel 
 * Currently only DFT point supported
 * \param gpuInfo    processing info
 */
extern "C"
void ObitCUDABeamModelDFTSetMod (GPUBeamInfo *gpuInfo)
{
  /* Not much to do here */
  return;
} /* end ObitCUDABeamModelDFTSetMod */

/**
 * Calculate an ObitCUDABeamModel with beam corrections
 * Only DFT supported
 * \param  gpuInfo   processing info
 * \param  doJones   1/0 to signal if Jones matrices need updating
 * \param  doUnit    1/0 use unit matrices as Jones
 * \param  parAng    Parallactic angle (rad)
 */
extern "C"
void ObitCUDABeamModelDFTCalc (GPUBeamInfo *gpuInfo, int doJones, int doUnit, float parAng, 
     long visRange[2])
{
   // Put available fast memory in L1 cache - no apparent effect
   //What??? checkCudaErrors(cudaFuncSetCacheConfig(dftPointKernel, cudaFuncCachePreferL1));

   /* Process with streams */
   int prtLv = 0;  // >=5 -> Debugging messages
  BeamProcessWithStreams(gpuInfo->nstream,
     doJones, doUnit, (float)parAng, visRange,
     gpuInfo->h_modelInfo, gpuInfo->d_modelInfo, 
     gpuInfo->h_visInfo, gpuInfo->d_visInfo, 
     gpuInfo->h_antInfo, gpuInfo->d_antInfo, 
     (cudaStream_t *)gpuInfo->stream, (cudaEvent_t *)gpuInfo->cycleDone, 
     gpuInfo->h_data, gpuInfo->h_data, gpuInfo->d_data, prtLv);
    return;
} /* end ObitCUDABeamModelDFTCalc */

#if HAVE_GPU==1  /* have a GPU? */
/**
 * Calculates DFT model using a GPU.
 * Multiple streams are used to overlap I/O and computation,
 * each call divides the data into streams_used pieces.
 * \param  streams_used  Number of streams to use
 * \param  doJones       1/0 to signal if Jones matrices need updating
 * \param  doUnit        1/0 use unit matrices as Jones
 * \param  parAng        Parallactic angle (rad)
 * \param  visRange      Range of visibility numbers to process in buffer (0-rel)
 * \param  nmodel        Number of sky model components
 * \param  h_modelInfo   Host resident modelInfo
 * \param  d_modelInfo   GPU resident modelInfo
 * \param  h_visInfo     Host resident visInfo
 * \param  d_visInfo     GPU resident visInfo
 * \param  h_antInfo     Host resident antInfo (Beam model)
 * \param  d_antInfo     GPU resident antInfo
 * \param  stream        GPU stream array [streams_used]
 * \param  cycleDone     GPU event array [streams_used]
 * \param  h_data_source Host resident input buffer, should be locked
 * \param  h_data_sink   Host resident output buffer, should be locked
 *                       may be h_data_source
 * \param  d_data        GPU resident data buffer (nvis/streams_used)
 * \param  prtLv         Print level, 5=>much directly printed
 */
static void BeamProcessWithStreams(int streams_used, 
            int doJones, int doUnit, float parAng, long visRange[2],
            GPUBeamModelInfo* h_modelInfo, GPUBeamModelInfo* d_modelInfo, 
            GPUBeamVisInfo* h_visInfo, GPUBeamVisInfo* d_visInfo, 
            GPUBeamAntInfo* h_antInfo, GPUBeamAntInfo* d_antInfo, 
            cudaStream_t* stream, cudaEvent_t* cycleDone,
            float *h_data_source, float *h_data_sink, 
            float *d_data[], int prtLv)
{
    int last_stream, current_stream = 0;
    int npass = streams_used;
    int nvisPass, lenvis = h_visInfo->lenvis;
    int nprod, dovis, nleft, nmodel, ntype= h_antInfo->numAntType;
    size_t lmemsize, memsize, off;
    int nload, nvis = (int)(visRange[1]-visRange[0]+1);
    size_t voff, curvindx, vindx = visRange[0];    /* Start in input vis buffer */
    size_t vmem, ocurvindx, ovindx = visRange[0];  /* Start in output vis buffer */
    dim3 numBlock, thPerBlock;
    int maxComp, maxType;
    nvisPass = (nvis+npass-1)/npass;  // nearly round up
    memsize = (lenvis*nvisPass)*sizeof(float);

    /* Update parallactic angle on antInfo */
    int nb=1, th=1;
    ParAngKernel <<<nb,th>>>
        (d_antInfo, parAng);
    // Wait for kernel to finish and check for execution errors
    checkCudaErrors(cudaDeviceSynchronize());

    if (prtLv>=5) printf ("\n\nEnter processing doJones %d doUnit %d parAng %f visRange %d %d lenvis %d nvisPass %d\n",
      doJones,doUnit,parAng,visRange[0],visRange[1],lenvis,nvisPass);

// Do the Jones/Unit matrices need updating?
    if (doJones||doUnit) {
      // looping over , type
      nmodel = h_modelInfo->nmodel;
      maxComp=MIN(NCHBLK, nmodel); maxType = 2;
      // x = chan, y = type
      numBlock.x = MAX(1,nmodel); numBlock.y = MAX(1,ntype); numBlock.z = 1;
      thPerBlock.x = MAX(1, maxComp); thPerBlock.y = MAX(1, maxType); thPerBlock.z = 1;
      if (doJones) {
        if (prtLv>=5) printf ("run GetJonesKernel numBlock %d %d %d thPerBlock %d %d %d\n",numBlock.x,numBlock.y,numBlock.z,thPerBlock.x,thPerBlock.y,thPerBlock.z);
        GetJonesKernel <<<numBlock, thPerBlock, 0, stream[0]>>>
          (d_antInfo, d_modelInfo,  parAng, prtLv);
      } else if (doUnit) {
        if (prtLv>=5) printf ("run UnitJonesKernel numBlock %d %d %d thPerBlock %d %d %d\n",numBlock.x,numBlock.y,numBlock.z,thPerBlock.x,thPerBlock.y,thPerBlock.z);
        UnitJonesKernel <<<numBlock, thPerBlock, 0, stream[0]>>>
          (d_antInfo, d_modelInfo);
      }
      // Check for launch errors (like invalid configuration)
      cudaError_t err = cudaGetLastError();
      if (err != cudaSuccess) printf("GetJones/UnitKernel Error: %s\n", cudaGetErrorString(err));

      // Wait for kernel to finish and check for execution errors
      checkCudaErrors(cudaDeviceSynchronize());
    } // End Jones update

    // Number of data products (channels) all Stokes
    nprod = h_visInfo->nchan;

    if (prtLv>=5) printf ("Start\n");
    if (prtLv>=5) printf ("nvis %d nvisPass %d lenvis %d memsize %d npass %d nchan %d nIF %d\n",nvis, nvisPass, lenvis, memsize, npass,h_visInfo->nchan,h_visInfo->nIF );
     
    // Do processing in a loop
    //
    // Note: All memory commands are processed in the order  they are issued,
    // independent of the stream they are enqueued in. Hence the pattern by
    // which the copy and kernel commands are enqueued in the stream
    // has an influence on the achieved overlap.


    // Upload first frame
    // curvindx = current start input vis index, vindx = next start input vis index
    curvindx = vindx; vindx = MIN(vindx+nvisPass, visRange[1]); 
    nload = (vindx-curvindx+1);
    lmemsize = (lenvis*nload)*sizeof(float);
    off = curvindx*lenvis;
    if (prtLv>=5) {
       voff = off/lenvis; vmem = (lmemsize/sizeof(float))/lenvis;
       printf ("upload current_stream %d off %d (%d) curvindx %d  vindx %d memsize %d (%d) nload %d nvisPass %d\n",
          0,off,voff,curvindx,vindx,lmemsize,vmem,nload, nvisPass);
    }
    checkCudaErrors(cudaMemcpyAsync(d_data[0],
			  	   &h_data_source[off],
				   lmemsize,
				   cudaMemcpyHostToDevice,
				   stream[0]));
    checkCudaErrors(cudaEventSynchronize(cycleDone[0]));
    nleft = nvis - nvisPass;  // How many left to copy?
    for (int i=0; i<npass; ++i) {
        int next_stream = (current_stream + 1) % streams_used;
	int prev_stream = current_stream - 1;
	if (prev_stream<0) prev_stream = streams_used-1;
	curvindx = vindx; vindx = MIN(vindx+nvisPass, visRange[1]); 
        nload = (vindx-curvindx+1);
        lmemsize = (lenvis*nload)*sizeof(float);
	off = curvindx*lenvis;  /* Offset in data buffers */
  	if (prtLv>=5) printf ("\n\nLoop %d prev %d current %d next %d curvindx %d vindx %d\n",i+1, prev_stream,current_stream,next_stream,curvindx,vindx);

	// Upload next frame - last may be less - or none
	if (nleft>0) {
  	  if (prtLv>=5) {
  	    voff = off/lenvis; vmem = (lmemsize/sizeof(float))/lenvis;
	    printf ("upload next_stream %d off %d (%d) memsize %d (%d) nleft %d nload %d curvindx %d vindx %d \n",
	      next_stream,off,voff,lmemsize,vmem,nleft,nload,curvindx,vindx);
	  }
	    checkCudaErrors(cudaMemcpyAsync(d_data[next_stream],
	  		      &h_data_source[off],
                              lmemsize,
                              cudaMemcpyHostToDevice,
                              stream[next_stream]));
         }

	// Ensure that processing and copying of the previous cycle has finished
	if (i>0) {
	  // ocurvindx = current start output vis index, ovindx = next start output vis index
	  ocurvindx = ovindx;  ovindx = MIN(ovindx+nvisPass, visRange[1]);
	  off = ocurvindx*lenvis;  /* Offset in data buffers */
	  if (prtLv>=5) {
	    voff = off/lenvis; vmem = (memsize/sizeof(float))/lenvis;
	    printf ("download prev_stream %d off %d (%d) memsize %d (%d) ocurvindx %d ovindx %d\n",
	      prev_stream,off,voff,memsize,vmem,ocurvindx,ovindx);
	  }
	  checkCudaErrors(cudaMemcpyAsync(&h_data_sink[off],
					  d_data[prev_stream],
					  memsize,
					  cudaMemcpyDeviceToHost,
					  stream[prev_stream]));
	  if (prtLv>=5) printf ("sync prev_stream %d loop %d\n",prev_stream, i);
	  cudaEventSynchronize(cycleDone[prev_stream]);
	}

        // Process current
	if (prtLv>=5) printf ("Process, nvis %d nch %d stream %d \n",nvisPass, h_visInfo->nchan, current_stream);
	// make sure to do all visibilities
	if (i==npass-1) dovis = nvis-i*nvisPass;
	else            dovis = nvisPass;

	// package work:
	// NCHBLK products per block (channel)  threadIdx.x+NCHBLK*blockIdx.x = prod number
	// 1 vis per block, blockIdx.y = vis number,
	numBlock.x = (nprod+(NCHBLK-1))/NCHBLK; numBlock.y  = dovis; numBlock.z = 1;
	thPerBlock.x = NCHBLK;                  thPerBlock.y = 1;    thPerBlock.z = 1;

        if (prtLv>=5) printf ("run model calc numBlock %d %d %d thPerBlock %d %d %d\n",numBlock.x,numBlock.y,numBlock.z,thPerBlock.x,thPerBlock.y,thPerBlock.z);
        // by model type
        switch (h_modelInfo->type) {
        case GPUBeamModTypePoint:
  	  beamPointKernel<<<numBlock, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo, d_antInfo, dovis);
          break;
        case GPUBeamModTypePointSpec:
  	  beamPointSpecKernel<<<numBlock, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo, d_antInfo, dovis, prtLv);
          break;
        case GPUBeamModTypePointTSpec:
  	  beamPointTSpecKernel<<<numBlock, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              //d_modelInfo, d_visInfo);
              d_modelInfo, d_visInfo, d_antInfo, dovis, prtLv);
          break;
        case GPUBeamModTypeGauss:
  	  beamGaussKernel<<<numBlock, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo, d_antInfo, dovis);
          break;
        case GPUBeamModTypeGaussSpec:
  	  beamGaussSpecKernel<<<numBlock, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo, d_antInfo, dovis);
          break;
        case GPUBeamModTypeGaussTSpec:
  	  beamGaussTSpecKernel<<<numBlock, thPerBlock, 0, stream[current_stream]>>>(
              d_data[current_stream],
              d_modelInfo, d_visInfo, d_antInfo, dovis);
          break;
          default:
              printf("Model type %d not supported in GPU",d_modelInfo->type);
          };  // End switch by model type  

       cudaError_t err = cudaGetLastError();  // Did it run?
       if (err != cudaSuccess) printf("Run Error: %s\n", cudaGetErrorString(err));

	// make sure previous frame done
	if (i>0) cudaEventSynchronize(cycleDone[prev_stream]);

	last_stream = current_stream;
        current_stream = next_stream;
	nleft -= nvisPass;
	if (prtLv>=5) printf ("End loop %d/%d curvindx %d vindx %d nleft %d\n",i+1,npass,curvindx,vindx,nleft);
    } /* end loop */

    /* Data from last pass */
    if (prtLv>=5) printf ("sync last_stream %d \n",last_stream);
    cudaEventSynchronize(cycleDone[last_stream]);
    
    // last piece may be uneven size
    dovis   = nvis-nvisPass*(npass-1);
    memsize = (lenvis*dovis)*sizeof(float);
    ocurvindx = ovindx;  ovindx = MIN(ovindx+nvisPass, visRange[1]);
    off = ocurvindx*lenvis;  /* Offset in data buffers */
    if (prtLv>=5) {
      voff = off/lenvis; vmem = (memsize/sizeof(float))/lenvis;
      printf ("download last_stream %d off %d (%d) memsize %d (%d) ocurvindx %d ovindx %d dovis %d\n",
        last_stream,off,voff,memsize,vmem,ocurvindx,ovindx,dovis);
    }
    checkCudaErrors(cudaMemcpyAsync(&h_data_sink[off],
				    d_data[last_stream],
				    memsize,
				    cudaMemcpyDeviceToHost,
				    stream[last_stream]));
    if (prtLv>=5) printf ("Finish\n");
    cudaDeviceSynchronize();

    return;
} /* end processWithStreams */
#endif /* HAVE_GPU */
