/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008-2009                                          */
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
/** Do Error analysis: */
gboolean doError;
/** make Primary beam correction? */
gboolean doPBCorr;
/** Do broken power law ? */
gboolean doBrokePow;
/** Size of planes in pixels */
olong nx, ny;
/** Maximum Chi square to accept a partial fit */
ofloat maxChi2;
/** Output Image descriptor */
ObitImageDesc *outDesc;
/** Array of pixel arrays for input data (nfreq) */
ObitFArray **inFArrays;
/** Array of BeamShape objects for input planes (nfreq) */
ObitBeamShape **BeamShapes;
/** Array of RMSes per inFArrays (nfreq) */
ofloat *RMS;
/** Array of calibration fractional error per inFArrays */
ofloat *calFract;
/** reference frequency */
odouble refFreq;
/** Array of frequencies (nfreq) */
odouble *freqs;
/** Array of pixel arrays for output data (2*nterm) */
ObitFArray **outFArrays;
