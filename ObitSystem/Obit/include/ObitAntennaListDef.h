/* $Id: ObitAntennaListDef.h,v 1.2 2006/03/06 21:19:11 bcotton Exp $ */
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
/*  Define the basic components of the ObitAntennaList structure.     */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitAntennaListDef.h
 * ObitAntennaListDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Number of entries */
olong number;
/** Array of source structures */
ObitAntenna **ANlist;
/** Polarization calibration type  OBIT_UVPoln_Approx, OBIT_UVPoln_VLBI
    OBIT_UVPoln_ELORI OBIT_UVPoln_XYLin */
ObitUVPolCalType polType;
/** Array center (m) */
odouble ArrayXYZ[3];
/** GST at IAT=0 (rad) */
odouble GSTIAT0;
/** Earth rotation rate (rad/day) */
odouble RotRate;
/** Julian date of reference day **/
odouble JD;
/** Polar offset (m) */
ofloat PolarXY[2];
/** UT1 - UTC (rad) */
ofloat ut1Utc;
/** Data time - UTC (rad) */
ofloat dataUtc;
/** Data time - IAT  (rad) */
ofloat dataIat;
/** R-L phase difference  1 per IF */
ofloat *RLPhaseDiff;
/** Number of polarization terms per feed */
olong numPoln;
/** Number of IFs (sets of poln calibration) */
olong numIF;
/** Number of polarization calibration constants */
olong  numPCal;
/** Polarization Reference antenna */
olong polRefAnt;
/** The FQid for which the polarization calibration is valid */
olong FreqID;
/** Time system, e.g. "UTC" */
gchar TimSys[12];
/** Array name */
gchar ArrName[12];
/** Is this the VLA (all parallactic angles the same ) */
gboolean isVLA;
