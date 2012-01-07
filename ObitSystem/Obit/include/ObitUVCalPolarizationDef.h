/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2012                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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
#ifndef OBITUVCALPOLARIZATIONDEF_H 
#define OBITUVCALPOLARIZATIONDEF_H 
#include "ObitAntennaList.h"
#include "ObitSourceList.h"
#include "ObitPolCalList.h"


/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalPolarizationDef.h
 * ObitUVCal utilities for applying polarization calibration to uv data 
 */

/*--------------Structure definitions-------------------------------------*/
/** Polarization calibration structure */
typedef struct {
#include "ObitDef.h"  /* Parent class definitions */
/** Polarization calibration type  OBIT_UVPoln_Approx, OBIT_UVPoln_VLBI
    OBIT_UVPoln_ELORI OBIT_UVPoln_XYLin */
  ObitUVPolCalType polType;
  /** Number of antennas in calibration table (actually max antenna no). */
  olong numAnt;
  /** Number of subarrays in the data */
  olong numSubA;
  /** Number of IFs in data. */
  olong numIF;
  /** Selected Subarray number. <=0-> all */
  olong SubA;
  /** Selected Frequency ID  number. <=0-> all */
  olong FreqID;
  /** Start IF number (1-rel) selected */
  olong bIF;
  /** highest IF (1-rel) selected */
  olong eIF;
  /** Number of channels in data */
  olong numChan;
  /** Start channel number (1-rel) */
  olong bChan;
  /** highest channel (1-rel) */
  olong eChan;
  /** Number of polarizations in the calibration table (1 or 2) */
  olong numPoln;
  /** Are corrections per channel? i.e. use AIPS PD table */
  gboolean perChan;
  /** Are all stokes correlations required */
  gboolean allStokes;
  /** Current source ID, -1 => single source */
  olong curSourID;
  /** Current Subarray */
  olong curSubA;
  /** Current time */
  ofloat curTime;
  /** Array  current parallactic angles, one per antenna 
   as cosine and sine */
  ofloat *curPA, *curCosPA, *curSinPA;
  /** Subarray inverse Jones matrices (2x2 complex per IF or 
      channel/IF) per antenna */
  ofloat **Jones;
  /** Polarization Mueller matrix for one baseline, all IFs or 
      channel/IFs */
  ofloat *PolCal;
  /** Current source RA (rad) */
  odouble curRA;
  /** Cosine current source Declination */
  odouble curCosDec;
  /** Sine current source Declination */
  odouble curSinDec;
  /** Polarication calibration from AIPS PD table */
  ObitPolCalList *PCal;
} ObitUVCalPolarizationS;
#endif /* OBITUVCALPOLARIZATIONDEF_H */ 


