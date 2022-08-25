/* $Id$        */
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
#ifndef OBITCUDAGRIDINFODEF_H 
#define OBITCUDAGRIDINFODEF_H 
#define GRID_STREAM_COUNT 2
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDAGridInfoDef.h
 *
 * Define interface structures with primitive CUDA code
 *
 */

/* Gridding info for GPU */
/*---------------------------------- Structures ----------------------------*/
/** Common info  */
typedef struct  {
/** Number of facets */
int nfacet;
/** Number of planes per facets */
int nplane;
/** Number of channels */
int nchan;
/** Number of IFs */
int nif;
/** Number of Stokes correlations */
int nstok;
/** Number of random parameters in visibilities */
int nrparm;
/** number of visibilities */
int nvis;
/** length in floats of visibilities */
int lenvis;
/** Convolution kernal width */
int convWidth;
/** Number of convolution entries per cell */
int convNperCell;
/** rotation to apply in imaging */
float rotate;
/** guardband in u, v in lambda */
float guardu, guardv;
/* Original nVisPIO on uv data */
int oldnVisPIO;
/** How many streams to use */
int nStream;
/** Print level */
int prtLv;
/* 0-rel Frequency plane per channel host*/
int *h_freqPlane;
/* 0-rel Frequency plane per channel device */
int *d_freqPlane;
/** Taper sigmas per channel/IF host */
float *h_sigma1, *h_sigma2, *h_sigma3;
/** Taper sigmas per channel/IF device */
float *d_sigma1, *d_sigma2, *d_sigma3;
/** Frequency scaling array (nchan*nif) host */
float *h_freqArr;
/** Frequency scaling array (nchan*nif) device */
float *d_freqArr;
/** gridding kernal order, given fractional cell fastest - host */
float *h_convfn;
/** gridding kernal order, given fractional cell fastest - device */
float *d_convfn;
/** Vis data locked buffer */
float *h_data_in;
/** Vis data device buffer array of GRID_STREAM_COUNT buffers. */
float **d_data_in;
/** d_data_in Array in locked host memory. */
float **h_d_data_in;
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
} GridInfo; 

/** Info per facet */
typedef struct  {
/* Making Beam (replace data w/ (1.0)? */
int  doBeam;
/** Size of facets (cells) per Facet */
int nx, ny;
/** Scaling for u,v,w to cells */
float uscale, vscale, wscale;
/** shift parameters per facet (dxc, dyc, dzc) */
float shift[3];
/** UVW, Rotation matrix per facet [3][3] */
float rotUV[9];
/** Max/min baseline */
float maxBL, minBL;
/** Beam Taper */
float bmTaper;
/** Size of each grid (nx x ny/2+1 x 2 floats)per facet  */
int sizeGrid;
/** Host facet grids each nplane x nx x ny/2+1 x 2 floats */
float *h_grid;
/** device facet grids each nplane x nx x ny/2+1 x 2 floats */
float *d_grid;
} FacetInfo; 

typedef struct  {
/** which cuda enabled GPU? */
long cuda_device;
/** How much GPU global memory */
size_t gpu_memory;
/* How many visibilities per I/O? */
long nVisPIO;
/* GPU base address */
void* d_base;
/** Information in common to all facets */
GridInfo* h_gridInfo;
GridInfo* d_gridInfo;
/** Per facet information */
FacetInfo** h_facetInfo;
FacetInfo** d_facetInfo;
/** d_facetInfo Array in locked host memory. */
FacetInfo** h_d_facetInfo;
/** DEBUG stuff */
float *h_grid_debug, *d_grid_debug;
} CUDAGridInfo; 

#endif /* OBITCUDAGRIDINFODEF_H */
