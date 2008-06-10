/* $Id: ObitOTFCalFlagDef.h,v 1.1.1.1 2004/07/19 17:04:44 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
#ifndef OBITOTFCALFLAGDEF_H 
#define OBITOTFCALFLAGDEF_H 

#include "Obit.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalFlagDef.h
 * ObitOTFCal utilities for applying flagging to OTF data 
 */

/*--------------Structure definitions-------------------------------------*/
/** calibration structure */
typedef struct {
  /** Calibration ObitTableOTFFlag table (as Obit*) */
  Obit *FlagTable;
  /** Flag Table Row (as Obit*) */
  Obit *FlagTableRow;
  /** Number of rows in flag table */
  olong numRow;
  /** Last Row read */
  olong LastRowRead;
  /** Number of Feeds in flag table. */
  olong numFeed;
  /** Number of Stokes' in data */
  olong numStok;
  /** First Stokes value in data */
  olong stoke0;
  /** Number of channels in data */
  olong numChan;
  /** Time (days) for which the flagging tables is current */
  ofloat flagTime;
  /** Maximum number of flag table entries  */
  olong maxFlag;
  /** Current number of flag table entries */
  olong numFlag;
  /** Target ID per flag entry  */
  olong *flagTarget;
  /** Feed number per flag entry  */
  olong *flagFeed;
  /**  First channel per flag entry */
  olong *flagBChan;
  /**  Highest channel per flag entry  */
  olong *flagEChan;
  /**  Flags for the polarizations per flag entry
   Note: there are 4*numFlag values */
  gboolean *flagPol;
  /**  End time  per flag entry  */
  ofloat *flagEndTime;
} ObitOTFCalFlagS;
#endif /* OBITOTFCALFLAGDEF_H */ 



