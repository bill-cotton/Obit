/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
/*  Define the basic components of the ObitDoppler structure          */
/*  This is intended to be included in a class structure definition   */
/* and to be used as the template for generating new classes derived  */
/* from Obit.                                                         */
/**
 * \file ObitDopplerDef.h
 * ObitDoppler structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** UV data set being operated on */
ObitUV *uvdata;
/** Temporary spectrum  */
ObitCArray *Spectrum;
/** work spectrum  */
ObitCArray *Work;
/** Forward FFT */
ObitFFT *FFTFor;
/** Reverse FFT */
ObitFFT *FFTRev;
/** Antenna List */
ObitAntennaList *antList;
/** source list */
ObitSourceList *sourceList;
/** Current source structure */
ObitSource *source;
/** Number of channels/IF */
olong nchan;
/** Number of IFs */
olong nif;
/** Number of polarizations */
olong npoln;
/** Year of reference day */
olong year;
/** Day of year of reference day */
olong doy;
/**JD of reference velocity. */
odouble JDref;
/** Rest frequency (Hz) of line */
odouble RestFreq;
/** Desired LSR velocity (km/s) */
ofloat VelLSR ;
/** Desired channel for VelLSR */
ofloat refChan;
