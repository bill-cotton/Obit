/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013,2014                                          */
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
/*  Define the basic components of the ObitUVWCalc structure          */
/*  This is intended to be included in a class structure definition   */
/* and to be used as the template for generating new classes derived  */
/* from ObitUVWCalc.                                                  */
/**
 * \file ObitUVWCalcDef.h
 * ObitUVWCalc structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** UV data */
ObitUV* myData;
/** Source list */
ObitSourceList *SouList;
/** (Modified) Antenna list */
ObitAntennaList *AntList;
/** Source structure for single source data */
ObitSource *mySource;
/** Current Source ppointer */
ObitSource *curSource;
/** Antenna number index into AntList */
olong *antIndex;
/** Maximum antenna number */
olong maxAnt;
/** Current source ID */
olong curSID;
/** Number of antennas in arrays */
olong nant;
/** Flip sign of u,v,w? (GMRT, LOFAR) */
gboolean doFlip;
/* Time of current antenna celestial coordinates */
ofloat curTime;
/** Lorentz contraction factor */
ofloat LorentzFact;
/** inverse of wavelength for curent source */
ofloat ilambda;
/** Cosine, Sine of current source rotation of u-v for differential precession */
ofloat cuvrot, suvrot;
/** Earth's rotation rate, turns per day */
odouble Rate;
/** Greenwich Sidereal Time at UT=0, hours*/
odouble GSTUTC0;
/** Current JD */
odouble JD;
/** Cos and Sin of Declination of current source */
odouble cDec, sDec;
/** Cos and Sin of RA of current source */
odouble cRA, sRA;
/** Position for diurnal abberation */
odouble obsPos[3];
/** Antenna coordinates in celestial frame */
odouble *xm, *ym, *zm;
