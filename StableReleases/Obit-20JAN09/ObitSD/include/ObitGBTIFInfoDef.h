/* $Id$ */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitGBTIFInfo structure.          */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitGBTIFInfoListDef.h
 * ObitGBTIFInfoDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
#define MAXNUMIF 50
/** Backend name */
gchar *Backend;
/** Polarization (R, L, X, Y) values */
gchar poln[MAXNUMIF];
/** Number of IFs in this structure */
olong nIF;
/** Number of frequency channels */
olong nchan;
/** Bank identifier */
gchar bank[MAXNUMIF][2];
/** Port identifier */
olong port[MAXNUMIF];
/** Feed number */
olong feed[MAXNUMIF];
/** index of first FEED of a sig/ref pair */
olong srfeed1[MAXNUMIF];
/** index of second FEED of a sig/ref pair */
olong srfeed2[MAXNUMIF];
/** Channel frequency increment (Hz) */
ofloat delta[MAXNUMIF];
/** Frequency reference pixel */
ofloat refPixel[MAXNUMIF];
/** Reference frequency (Hz) */
odouble refFrequency[MAXNUMIF];
