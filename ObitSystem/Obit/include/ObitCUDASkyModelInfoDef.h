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
#ifndef OBITCUDASKYMODELINFODEF_H 
#define OBITCUDASKYMODELINFODEF_H 
#define STREAM_COUNT 4
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDASkyModelInfoDef.cuh
 *
 * Define interface structures with primitive CUDA code
 *
 */

/*----------------------------------- enumerations -------------------------*/
/**
 * \enum gpuModType
 * enum for sky model types supported in GPU
 */
enum gpuModType {
  GPUModTypePoint=0,        /* Simple point model */
  GPUModTypePointSpec,      /* Point model with parameterized spectrum */
  GPUModTypePointTSpec,     /* Point model with tabulated spectrum */
  GPUModTypeGauss,          /* Simple Gaussian model */
  GPUModTypeGaussSpec,      /* Gaussian model with parameterized spectrum */
  GPUModTypeGaussTSpec,     /* Gaussian model with tabulated spectrum */
}; /* end enum  gpuModType */
typedef enum gpuModType GPUModType;

/**
 * \enum gpuOpType
 * enum for sky model operation types supported in GPU
 */
enum gpuOpType {
  GPUOpTypeSub=0,          /* Subtract model */
  GPUOpTypeDiv,            /* Divide model */
  GPUOpTypeRepl,           /* Replace data with model */
}; /* end enum  gpuOpType */
typedef enum gpuOpType GPUOpType;

/*---------------------------------- Structures ----------------------------*/
typedef struct {
  /* Number of visibilities in buffer */
  int nvis;
  /* Number of data products */
  int nprod;
  /* Number random parameters */
  int nrparm;
  /* Length of visibilities in buffer */
  int lenvis;
  /* Number channels */
  int nchan;
  /* offset in vis of u */
  int ilocu;
  /* offset in vis of v */
  int ilocv;
  /* offset in vis of w */
  int ilocw;
  /* offset in vis of time */
  int iloct;
  /* offset in vis of baseline */
  int ilocb;
  /* offset in vis of source */
  int ilocsu;
  /* First channel (0-rel) */
  int chanb;
  /* highest channel (0-rel) */
  int chane;
  /* vis (not float) increment between channels */
  int incf;
  /* Number IFs */
  int nIF;
  /* First IF (0-rel) */
  int IFb;
  /* highest IF (0-rel) */
  int IFe;
  /* vis increment between IFs */
  int incif;
  /* Number Stokes */
  int nstok;
  /* First Stokes (0-rel) */
  int stokb;
  /* highest Stokes (0-rel) */
  int stoke;
  /* vis increment between Stokes */
  int incs;
  /* stokes type 0=>IQUV, 1=>RR,LL,RL,LR, 2=XX,YY,XY,YX */
  int stype;
  /* increment in freqScale for IF */
  int kincif;
  /* increment in freqScale for frequency */
  int kincf;
  /* Ratio of UV ref. freq. to image ref. freq. (per coarse channel for TSpec) */
  float *freqRat;
  /* GPU pointer for freqRat */
  float *d_freqRat;
  /* GPU Frequency factor array dim nchan*nif */
  float *freqScale;
 } GPUVisInfo;

typedef struct {
  /* Number of spectral terms */
  int nterm;
  /* Size of components */
  int size;
  /* Swap R/I? */
  int doSwap;
  /* Type of components */
  GPUModType type;
  /* Type of operation */
  GPUOpType opType;
  /* array of coarse frequency index per frequency channel */
  int *specIndex;
  /* GPU address of specIndex */
  int *d_specIndex;
  /* Number of components per Stokes */
  int nmodel[4];
  /* GPU Model Data array per Stokes */
  float *model[4];
  /* Factors per Stokes */
  float stokFact[4];
 } GPUModelInfo;

typedef struct {
  /* Number of visibilities in buffer */
  int nvis;
  /* Number channels */
  int nchan;
  /* Number random parameters */
  int nrparm;
  /* Length of visibilities in buffer */
  int lenvis;
  /* size of model components */
  int modelSize;
  /* number of streams */
  int nstream;
  /* Print level */
  int prtLv;
  /* Old nVisPIO */
  int oldnVisPIO;
  /* Which GPU */
  int cuda_device;
  /* data buffer size in floats */
  int bufferSize;
  /* number of model components  per poln */
  int nmodel[4];
  /* device model component array per poln */
  float *d_model[4];
  /* Frequency scaling array  */
  float *d_freq;
  /* host data buffer */
  float *h_data;
  /* device data buffer, per stream  */
  float **d_data;
  /* description of visibilities */
  GPUVisInfo *d_visInfo, *h_visInfo;
  /* description of model */
  GPUModelInfo *d_modelInfo, *h_modelInfo;
  /* CUDA dependent values */
 #if IS_CUDA==1
   /* stream structures */
   cudaStream_t *stream;
   /* Stream events  */
   cudaEvent_t *cycleDone;
 #else  /* Not CUDA */
   /* stream structures (cudaStream_t) */
   int **stream;
   /* Stream events (cudaEvent_t) */
   int **cycleDone;
 #endif /* IS_CUDA */
 } GPUInfo;

#endif /* OBITCUDASKYMODELINFODEF_H */
