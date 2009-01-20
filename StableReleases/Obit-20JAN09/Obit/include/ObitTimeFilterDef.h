/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003,2008                                          */
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
/*  Define the basic components of the ObitTimeFilter structure              */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitTimeFilterDef.h
 * ObitTimeFilter structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Array of Time series */
ObitFArray **timeSeries;
/** Pointers to data in time arrays */
ofloat **timeData;
/** Times in series */
ofloat *times;
/** Array of Frequency series */
ObitCArray **freqSeries;
/** Pointers to data in frequency arrays */
ofloat **freqData;
/** Frequencies in series */
ofloat *freqs;
/** Forward (time->freq) FFT */
ObitFFT *FFTFor;
/** Reverse (freq->time) FFT */
ObitFFT *FFTRev;
/** Interpolation object to deblank time series */
ObitFInterpolate *interp;
/** Number of times in series */
olong nTime;
/** Number of frequencies in series */
olong nFreq;
/** Number of series */
olong nSeries;
/** Time step in series */
ofloat dTime;
/** Frequency step in series */
ofloat dFreq;
