/* $Id$  */
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
/*  Define the basic components of the ObitSource structure.          */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitSourceListDef.h
 * ObitSourceDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Source ID number */
olong SourID;
/** Source qualifier */
olong Qual;
/** Number of IFs */
olong numIF;
/** Source Name */
gchar SourceName[20];
/** Source cal code */
gchar CalCode[8];
/** Standard Equinox  */
ofloat equinox;
/** Source Apparent RA (deg) */
odouble RAApp;
/** Source Apparent Declination (deg) */
odouble DecApp;
/** Source  RA at standard equinox (deg) */
odouble RAMean;
/** Source Declination at standard equinox (deg) */
odouble DecMean;
/** Proper motion (deg/day) in RA */
odouble  PMRa;
/** Proper motion (deg/day) in declination */
odouble  PMDec;
/** Total Stokes I flux density per IF */
gfloat*  IFlux;
/** Total Stokes Q flux density per IF */
gfloat*  QFlux;
/** Total Stokes U flux density per IF */
gfloat*  UFlux;
/** Total Stokes V flux densityper IF */
gfloat*  VFlux;
/** Frequency offset (Hz) from IF nominal per IF */
gdouble*  FreqOff;
/** Bandwidth */
odouble  Bandwidth;
/** LSR velocity per IF */
gdouble*  LSRVel;
/** Line rest frequency per IF */
gdouble*  RestFreq;
