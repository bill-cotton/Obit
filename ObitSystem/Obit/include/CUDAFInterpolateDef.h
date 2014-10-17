/* $Id: $                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*  Define the basic components of the CUDAFInterpolate structure     */
/*  This is intended to be included in a class structure definition   */
/**
 * \file CUDAFInterpolateDef.h
 * CUDAInterpolate structure members for derived classes.
 * This is a limited GPU implementation of ObitFInterpolate
 */
/** Structure in GPU memory */
float *d_FInterpolate;
/** Array to be interpolated (Host) */
CUDAFArray *h_myArray;
/** Array to be interpolated (GPU) */
float *d_myArray;
/** Dimension of plane in myArray */
int nx, ny;
/** Half width of interpolation kernal */
int hwidth;
/** Reciprocals of Lagrangian denominators */
float denom[10];
/** CUDA stuff */
/** Number of streams */
int nstream;
/** Size of buffers */
int nxBuff, nyBuff;
/** GPU work buffers */
float **d_data;
#if IS_CUDA==1  /* CUDA code */
cudaStream_t *stream;
/* Stream events  */
cudaEvent_t *cycleDone;
#else  /* Not CUDA */
/* stream structures (cudaStream_t) */
int **stream;
/* Stream events (cudaEvent_t) */
int **cycleDone;
#endif /* IS_CUDA */
