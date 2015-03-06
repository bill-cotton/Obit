/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012-2015                                          */
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
/*  Define the basic components of the ObitSourceEphemerus structure  */
/*  This is intended to be included in a class structure definition   */
/* and to be used as the template for generating new classes derived  */
/* from Obit.                                                         */
/**
 * \file ObitSourceEphemerusDef.h
 * ObitSourceEphemerus structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Number of entries given */
olong nentry;
/** Last source number */
olong lastSrc;
/** Last entry number */
olong lastEntry;
/** Time (day) of last calculation */
ofloat lastTime;
/** Time (day) until which last calculation is valid */
ofloat validTime;
/** last uv rotation angle (rad) */
ofloat lastUVrot;
/** Last RA calculated (rad) */
odouble lastRA;
/** Last Dec calculated (rad) */
odouble lastDec;
/** Last Distance calculated */
odouble lastDist;
/** How often to update position (day) */
odouble updtime;
/* Reference day as JD */
odouble refJD;
/** List of source IDs per entry */
olong *SID;
/** Reference time wrt reference day (days) per entry */
odouble *refTime;
/** beginning of validity (days) per entry */
odouble *startTime;
/** end of validity (days) per entry */
odouble *endTime;
/** Apparent RA at reference time (deg), per entry*/
odouble *RARef;
/** number of RA derivatives per entry */
olong *numRADeriv;
/** Array of arrays of RA derivatives per entry, deg/day... */
odouble **RADeriv;
/** Apparent Declination at reference time (deg), per entry */
odouble *DecRef;
/** number of Dec derivatives per entry */
olong *numDecDeriv;
/** Array of arrays of Dec derivatives per entry, deg/day... */
odouble **DecDeriv;
/** Distance at reference time (m), per entry */
odouble *distRef;
/** number of distance derivatives per entry */
olong *numDistDeriv;
/** Array of arrays of distance derivatives per entry, m/day... */
odouble **DistDeriv;
/** Work source for precessing */
ObitSource *source;
/** UV data descriptor */
ObitUVDesc *uvDesc;
