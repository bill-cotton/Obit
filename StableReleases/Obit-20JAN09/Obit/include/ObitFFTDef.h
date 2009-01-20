/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitFFT structure              */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitFFTDef.h
 * ObitFFT structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Requested FFT direction */
ObitFFTdir dir;
/** Requested FFT type */
ObitFFTtype type;
/** Rank of matrix (<=7)*/
olong rank;
/** dimensionality of array */
olong dim[7];
#if HAVE_FFTW3==1  /* FFTW3 version */
/** FFTW3 plan for full Complex transforms */
fftwf_plan CPlan;
/** FFTW3 plan for Half Complex/Real transforms */
fftwf_plan RPlan;

#elif HAVE_FFTW==1  /* FFTW 2 version */
/** FFTW plan for full Complex transforms*/
fftwnd_plan CPlan;
/** FFTW plan for Half Complex/Real transforms */
rfftwnd_plan RPlan;

#elif HAVE_GSL==1  /* Else try GSL version */
/* Suppport up to rank 3 */
/** Stride on each axis */
size_t stride[3];
/** Wavetable for full float Complex transforms*/
gsl_fft_complex_wavetable_float *CWavetab[3];
/** Wavetable forfloat  half plane complex to real transforms*/
gsl_fft_halfcomplex_wavetable_float *HCWavetab[3];
/** Wavetable for float real to half plane complex transforms*/
gsl_fft_real_wavetable_float *RWavetab[3];
/** Workspace for float real to half plane complex transforms*/
gsl_fft_real_workspace_float *Rwork[3];
/** Workspace for float real to half plane complex transforms*/
gsl_fft_complex_workspace_float *Cwork[3];

/** Wavetable for full double Complex transforms*/
gsl_fft_complex_wavetable *dblCWavetab[3];
/** Wavetable for double half plane complex to real transforms*/
gsl_fft_halfcomplex_wavetable *dblHCWavetab[3];
/** Wavetable for double real to half plane complex transforms*/
gsl_fft_real_wavetable *dblRWavetab[3];
/** Workspace for double real to half plane complex transforms*/
gsl_fft_real_workspace *dblRwork[3];
/** Workspace for double real to half plane complex transforms*/
gsl_fft_complex_workspace *dblCwork[3];
#endif /* HAVE_GSL */
