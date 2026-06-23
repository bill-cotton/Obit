/* $Id$        */
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
#ifndef OBITCUDABEAMMODELINFODEF_H 
#define OBITCUDABEAMMODELINFODEF_H 
#define BEAM_STREAM_COUNT 2
#include "CUDAMath.h"
//#include "ObitGPUFInterpolate.h"
#include "CUDAFInterpolate.h"
//NO CAN DO #include "ObitAntennaList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDABeamModelInfoDef.h
 *
 * CUDA implemention of skyModel with Beam corrections
 * Define interface structures with primitive CUDA code
 *
 */

/* Fooey - need local definition */
#if HAVE_GPU==1  /* Have a GPU? */
#ifndef CUDABEAMINTERPDEF_H // Prevent multiple definitions
#define CUDABEAMINTERPDEF_H
typedef struct {
#include "CUDABeamInterpDef.h"   /* CUDA stuff */
} CUDABeamInterp;
#endif /* CUDABEAMINTERPDEF_H */ 
#endif /* HAVE_GPU */

/*----------------------------------- enumerations -------------------------*/
/**
 * \enum gpuModType
 * enum for sky model types supported in GPU
 * Should be the equivalents of what's in ObitCUDASkyModelInfoDef.h
 */
enum gpuBeamModType {
  GPUBeamModTypePoint=0,        /* Simple point model */
  GPUBeamModTypePointSpec,      /* Point model with parameterized spectrum */
  GPUBeamModTypePointTSpec,     /* Point model with tabulated spectrum */
  GPUBeamModTypeGauss,          /* Simple Gaussian model */
  GPUBeamModTypeGaussSpec,      /* Gaussian model with parameterized spectrum */
  GPUBeamModTypeGaussTSpec,     /* Gaussian model with tabulated spectrum */
}; /* end enum  gpuBeamModType */
typedef enum gpuBeamModType GPUBeamModType;

/**
 * \enum gpuOpType
 * enum for sky model operation types supported in GPU
 */
enum gpuBeamOpType {
  GPUBeamOpTypeSub=0,          /* Subtract model */
  GPUBeamOpTypeDiv,            /* Divide model */
  GPUBeamOpTypeRepl,           /* Replace data with model */
}; /* end enum  gpuBeamOpType */
typedef enum gpuBeamOpType GPUBeamOpType;

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
  /* number of subbands */
  int nSpec;
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
  /* increment in d_freqScale for IF */
  int kincif;
  /* increment in d_freqScale for frequency */
  int kincf;
  /* Ratio of UV ref. freq. to image subband frequencies to  (per coarse channel for TSpec) */
  float *h_freqRat;
  /* GPU pointer for freqRat */
  float *d_freqRat;
  /* GPU Frequency factor array dim nchan*nif */
  float *d_freqScale;
 } GPUBeamVisInfo;

typedef struct {
  /* Number of spectral terms */
  int nterm;
  /* number of subbands */
  int nSpec;
  /* Size of components */
  int size;
  /* Swap R/I? */
  int doSwap;
  /* Stokes Type of components , 1,2,3,4=>I,Q,U,V */
  GPUBeamModType type;
  /* Type of operation */
  GPUBeamOpType opType;
  /* Stokes parameter I,Q,U,V = 1,2,3,4 */
  int stokType;
  /* array of coarse frequency (subband) index per frequency channel */
  int *h_specIndex;
  /* GPU address of specIndex */
  int *d_specIndex;
  /* Number of components */
  int nmodel;
  /* GPU Model Data array  
    Tabulated spectrum Model:
      +0 amp (?)               = single subband/channel flux
      +1,2,3                   = x,y,z offsets
      +4,5,6                   = Gaussian terms
      +5 (Pt) + 2*subband      = flux
      +8 (Gauss) + 2*subband   = flux
      +5 (Pt) + 2*subband+1    = SI
      +8 (Gauss) + 2*subband+1 = SI
      +4..nterm                = nterm frequency expansion terms
  */
  float *d_model;
  /* Factors per Stokes product xx,yy,xy,yx ot RR,LL,RL,LR */
  float stokFact[4];
 } GPUBeamModelInfo;

// Antenna/beam information
typedef struct {
  /* Number of antennas */
  int numAnt;
  /* Number of antenna types */
  int numAntType;
  /* Number of channels for Jones arrays = planes in beam images */
  int numJChan;
  /* Baseline random parameter if given (nonnegative) */
  int ilocb;
  /* Else antenna & subaray random parameters */
  int iloca1, iloca2, ilocsa;
  /* sine and cosine of twice the parallacticangles used in Jones */
  float s2p, c2p;
  /* Reference antenna diameter (m), factor */
  float refantSize, refSizeFact;
  /* Beam channel frequencies */
  float *h_BeamFreq, *d_BeamFreq;
  /* Antenna type index per antenna */
  int *h_AntType, *d_AntType;
  /* Table of Jones channels per vis channel */
  int *h_visJChann, *d_visJChann;
  /* Array of Jones matrices ant, component, channel
     numAnt*nmodel*numJChan - device only */
 #if HAVE_GPU==1  /* Have a GPU? */
 cudaMatx *d_Jones;
  /* Beam interpolators by type, antenna, Stokes type, r/i - slowest to fastest */
  CUDABeamInterp **d_BeamInterp;
#endif /* HAVE_GPU */
} GPUBeamAntInfo;

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
  /* number of model components */
  int nmodel;
  /* Current time range (day) */
  float timeRange[2];
  /* Allowed range of parallactic angle (deg) */
  float delta_ParAng;
  /* Current Parallactic angle */
  float curParAng;
  /* Current range of visibilities in buffer (0-rel) */
  long visRange[2];
  /* Use Unit matrices rather than from beam images? */
  int doUnit;
  /* Have Unit matrices been created */
  int haveUnit;
  /* device model component array */
  float *d_model;
  /* Frequency scaling array  */
  float *d_freq;
  /* host data buffer */
  float *h_data;
  /* device data buffer, per stream  */
  float **d_data;
  /* description of visibilities */
  GPUBeamVisInfo *d_visInfo, *h_visInfo;
  /* description of model */
  GPUBeamModelInfo *d_modelInfo, *h_modelInfo;
  /* description of antenna beam model */
  GPUBeamAntInfo *d_antInfo, *h_antInfo;
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
 } GPUBeamInfo;

#endif /* OBITCUDABEAMMODELINFODEF_H */
