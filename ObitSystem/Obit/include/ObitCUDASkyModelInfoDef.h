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
  /* GPU Frequency factor array dim nchan*nif */
  float *freqScale;
 } GPUVisInfo;

typedef struct {
  /* Number of components */
  int nmodel;
  /* Size of components */
  int size;
  /* Type of components */
  int type;
  /* GPU Model Data array */
  float *model;
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
  /* number of model components */
  int nmodel;
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
  /* model component array */
  float *d_model;
  /* Frequency scaling array  */
  float *d_freq;
  /* host data buffer */
  float *h_data;
  /* device input data buffer, per stream  */
  float **d_data_in;
  /* device output data buffer, per stream */
  float **d_data_out;
  /* description of visibilities */
  GPUVisInfo *d_visInfo, *h_visInfo;
  /* description of model */
  GPUModelInfo *d_modelInfo;
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
