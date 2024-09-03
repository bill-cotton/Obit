/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2024                                               */
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
#ifndef OBITUVCALJONESDEF_H 
#define OBITUVCALJONESDEF_H 
#include "ObitAntennaList.h"
#include "ObitSourceList.h"
#include "ObitPolCalList.h"


/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalJonesDef.h
 * ObitUVCal utilities for applying polarization calibration to uv data 
 */

/*--------------Structure definitions-------------------------------------*/
/** Polarization calibration structure */
typedef struct {
#include "ObitDef.h"  /* Parent class definitions */
  /** Calibration JI table (as Obit*) Single source */
  Obit *JITable;
  /** JI Table Row (as Obit*) */
  Obit *JITableRow;
  /** Calibration JT table (as Obit*) Multi Source */
  Obit *JTTable;
  /** JT Table Row (as Obit*) */
  Obit *JTTableRow;
  /** Number of rows in calibration table */
  olong numRow;
  /** Last Row read */
  olong LastRowRead;
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
  /** Using JI table? */
  gboolean useJI;
  /** Using JT table? */
  gboolean useJT;
  /** Are array feeds circularly polarized (else linear) */
  gboolean circFeed;
  /** Is this MeerKAT data? */
  gboolean isMeerKAT;
  /** Current source ID, -1 => single source */
  olong curSourID;
  /** Current Subarray */
  olong curSubA;
  /** Prior source ID in cal table */
  olong PriorSourID;
  /** Following source ID in cal table */
  olong FollowSourID;
  /** Prior time in Jones cal table */
  ofloat PriorCalTime;
  /** Prior time per antenna */
  ofloat *PriorAntTime;
  /** Following Time in Jones cal table */
  ofloat FollowCalTime;
  /** Following time per antenna */
  ofloat *FollowAntTime;
  /** time of current Jones calibration  */
  ofloat CalTime;
  /** Current time */
  ofloat curTime;
  /** Prior calibration arrays from JI/JT table, in order
   *    By antenna {numAnt}
   *       By IF {bIF->eIF}
   *          By channel {bChan->eChan}
   *            8 elements of 2x2 complex Jones matrix */
  ofloat *CalPrior;
  /** Following calibration arrays from JI/JT table - like CalPrior */
  ofloat *CalFollow;
  /** Array  current parallactic angles, one per antenna 
   as cosine and sine */
  ofloat *curPA, *curCosPA, *curSinPA;
  /** Subarray Jones inverse matrices (2x2 complex) in order
   *    By antenna {numAnt}
   *       By IF {bIF->eIF}
   *          By channel {bChan->eChan} 
     Calibrate by: Jones1^-1 x vis x Jones2^-1^H */
  ObitMatx **Jones;
  /** Work matrix */
  ObitMatx *TempMatx, *TempMatx2;
  /** Rotation matrix to convert circular basis to linear */
  ofloat *C2L_Matrix;
  /** Current source RA (rad) */
  odouble curRA;
  /** Cosine current source Declination */
  odouble curCosDec;
  /** Sine current source Declination */
  odouble curSinDec;
} ObitUVCalJonesS;
#endif /* OBITUVCALJONESDEF_H */ 


