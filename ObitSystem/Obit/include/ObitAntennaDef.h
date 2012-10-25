/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2012                                          */
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
/*  Define the basic components of the ObitAntenna structure.          */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitAntennaListDef.h
 * ObitAntennaDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Antenna name */
gchar AntName[8];
/** Antenna location (m) */
odouble AntXYZ[3];
/** Antenna E. Longitude (rad) */
odouble AntLong;
/** Antenna Latitude (rad) */
odouble AntLat;
/** Antenna distance from earth center (m) */
odouble AntRad;
/** Feed A feed position angle (deg) */
ofloat FeedAPA;
/** Feed B feed position angle (deg) */
ofloat FeedBPA;
/** Feed A Polarization calibration parameters */
ofloat *FeedAPCal;
/** Feed B Polarization calibration parameters */
ofloat *FeedBPCal;
/** Antenna ID number (1-rel), <=0 -> no information*/
olong AntID;
/** Antenna Mount Type 0=altaz, 1=equatorial, 2=orbiting */
olong AntMount;
/** Number of polarization calibration parameters per feed */
olong numPCal;
/** Feed A Feed Poln type 'R' 'L', 'X', 'Y' */
gchar FeedAType;
/** Feed B Feed Poln type 'R' 'L', 'X', 'Y' */
gchar FeedBType;
