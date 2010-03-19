/* $Id: ObitUVGSolveDef.h 2 2008-06-10 15:32:27Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
/*  Define the basic components of the ObitUVGSolve structure          */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitUVGSolveDef.h
 * ObitUVGSolve structure members for this and any derived classes.
 */
#include "ObitUVGSolveDef.h"  /* Parent class definitions */
/** Maximum antenna number (1-rel) */
olong maxAnt;
/** Number of IFs */
olong numIF;
/** Number of poln */
olong numPoln;
/** Print level for diagnostics */
olong prtLv;
/** Reference antenna number (1-rel) */
olong refAnt;
/** Data for current scan  */
ScanData *scanData;
/** Complex gain (r,i) per antenna, IF,poln IF varies fastest */
ofloat *antGain;
/** Group delay (sec) per antenna, IF,poln IF varies fastest */
ofloat *antDelay;
/** Dispersion term (?) per antenna, IF,poln IF varies fastest */
ofloat *antDisp;
/** Weight per antenna, IF,poln IF varies fastest */
ofloat *antWeight;
/** Work complex arrays */
ObitCArray *cWork1, *cWork2, *cWork3;
/** Work float arrays */
ObitFArray *fWork1, *fWork2, *fWork3;
/** FFT for coarse fringe search  */
ObitFFT *myFFT;
/** FFT over sampling */
olong FFTOverSamp;
/** Complex work array for FFT fit */
ObitCArray *FFTFitArray;
/** Interpolator for FFT fit */
ObitCInterpolate* myCInterp;
/* GSL Stuff */
#if HAVE_GSL==1  /* GSL stuff */
/** Coarse solution fitter */
gsl_multifit_fdfsolver* myCoarseSolver;
/** Coarse solution fitter function */
gsl_multifit_function_fdf *myCoarseFunc;
#endif /* GSL stuff */
/** Number of coefficients in coarse fit */
olong ncoefCoarse;
/** Number of data points in coarse fit */
olong ndataCoarse;
/** Coarse soln initial values */
double *coarse_coef_init;
/** Coarse soln data values */
ofloat *coarseData;
/** Coarse soln data weights */
ofloat *coarseWt;
/** Coarse Frequency offset from ref. (GHz) */
ofloat *coarseFreq;

