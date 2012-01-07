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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitPolCalList structure.      */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitPolCalListDef.h
 * ObitAntennaListDef structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Polarization calibration type  OBIT_UVPoln_Approx, OBIT_UVPoln_VLBI
    OBIT_UVPoln_ELORI OBIT_UVPoln_XYLin */
ObitUVPolCalType polType;
/** Number of antenna entries */
olong numAnt;
/** Number of polarization terms per feed */
olong numPoln;
/** Number of IFs (sets of poln calibration) */
olong numIF;
/** Number of channels per IF */
olong numChan;
/** Number of polarization calibration constants */
olong  numPCal;
/** Polarization Reference antenna */
olong polRefAnt;
/** R-L phase difference  1 per channel/IF */
ofloat *RLPhaseDiff;
/** Array of antenna based lists, 
    each antenna list has numPCal entries per channel/IF in order:
    P1p1, P1p2, P2p1, P2p2
    where P is polarization, e.g. "R", "L",
    and p is parameter, (r,i) for D terms, else (ellipticity, orientation)
*/
ofloat **ANlist;
