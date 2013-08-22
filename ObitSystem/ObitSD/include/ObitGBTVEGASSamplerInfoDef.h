/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
/*  Define the basic components of the ObitGBTDCRSampler\Info structure.    */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitGBTVEGASSampler\InfoListDef.h
 * ObitGBTVEGASSampler\InfoDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
#define MAXNUMVEGASSAMPLER 10
/** This will always be the same as the BANK keyword in the Primary HDU. */
gchar  bank_a[MAXNUMVEGASSAMPLER];
/** The PORT corresponding to the ``A'' input.  This will usually be the same PORT */
gchar  port_a[MAXNUMVEGASSAMPLER];
/** This will always be the same as the BANK keyword in the Primary HDU. */
gchar  bank_b[MAXNUMVEGASSAMPLER];
/** The PORT corresponding to the ``B'' input.  This will usually be the same PORT */
gchar  port_b[MAXNUMVEGASSAMPLER];
/** 0-rel index indicationg which sub-band this SAMPLER corresponds to. */
gshort  subband[MAXNUMVEGASSAMPLER];
/** The frequency at the center of the channel given by the CRPIX1 keyword. */
odouble  crval1[MAXNUMVEGASSAMPLER];
/** The increment between adjacent channels.  Note that the total bandwidth for a given */
odouble  cdelt1[MAXNUMVEGASSAMPLER];
/** The effective channel resolution. */
odouble  freqres[MAXNUMVEGASSAMPLER];
/** One of 'REAL' or 'IMAG' These describe the type of data found at each element alon */
gchar  datatype[MAXNUMVEGASSAMPLER][4];
/** Number of samplers */
olong nVEGASSampler;
