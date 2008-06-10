/* $Id: ObitOTFArrayGeomDef.h,v 1.1.1.1 2004/07/19 17:04:43 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
/*  Define the basic components of the ObitCArray structure           */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitOTFArrayGeom.h
 * ObitOTFArrayGeom structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Number of detectors */
olong numberDetect;
/** Offsets in azimuth per detector in degrees */
ofloat *azOffset;
/** Offsets in elevation per detector in degrees */
ofloat *elOffset;
/** Reference date as "YYYYMMDD"*/
gchar RefDate[12]; 
/**Time system, 'IAT' or 'UTC' */
gchar TimeSys[4];
/** Earth centered telescope coordinates */
odouble TeleX, TeleY, TeleZ;
/** Earth rotaiton rate, deg/IAT day */
odouble DegDay;
/** GST at time=0 (degrees) on the reference date */
odouble GSTiat0;
/** Polar position (x,y) on ref. date */
ofloat PolarX,PolarY;
/** UT1-UTC  (time sec.) */
ofloat ut1Utc;
/** data time-UTC  (time sec.) */
ofloat dataUtc;
/** IAT - UTC (sec). */
ofloat iatUtc;
/** telescope latitude in radians */
ofloat lat;
/** telescope longitude in radians */
ofloat lon;
/** LST at iat0 in radians */
ofloat LSTiat0;
/** Earth rotation rate in rad/day */
ofloat RadDay;
/** Data - IAT in days */
ofloat dataIat;
