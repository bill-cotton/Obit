/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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
/*  Define the basic components of the ObitSpectrumFit structure      */
/**
 * \file ObitSpectrumFitDef.h
 * ObitSpectrumFit structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Number of terms in spectral polynomial */
olong nterm;
/** Number of frequency data points to be fitted */
olong nfreq;
/** Size of planes in pixels */
olong nx, ny;
/** Minimum pixel SNR to fit */
ofloat minSNR;
/** Output Image descriptor */
ObitImageDesc *outDesc;
/** Array of pixel arrays for input data (nfreq) */
ObitFArray **inFArrays;
/** Array of BeamShape objects for input planes (nfreq) */
ObitBeamShape **BeamShapes;
/** Array of RMSes per inFArrays */
ofloat *RMS;
/** Array of calibration fractional error per inFArrays */
ofloat *calFract;
/** Array of pixel arrays for output data (2*nterm) */
ObitFArray **outFArrays;

/* GSL fitting stuff */
#ifdef HAVE_GSL
/** Matrix of powers of log(nu) per datum */
gsl_matrix *X;
/** Vector of log(S)=data to be fitted */
gsl_vector *y;
/** Data weights */
gsl_vector *w;
/** Vector of fitted coeffieients */
gsl_vector *coef;
/** Covariance matrix from fitting */
gsl_matrix *covar;
/** Fitting work space */
gsl_multifit_linear_workspace *work;
#endif /* HAVE_GSL */ 
