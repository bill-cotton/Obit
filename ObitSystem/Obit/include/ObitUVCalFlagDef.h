/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003,2010                                          */
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
#ifndef OBITUVCALFLAGDEF_H 
#define OBITUVCALFLAGDEF_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitUVCal.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalFlagDef.h
 * ObitUVCal utilities for applying flagging to uv data 
 */

/*--------------Structure definitions-------------------------------------*/
/** Amp/phase/delay/rate calibration structure */
typedef struct {
  /** Calibration FG table (as Obit*) */
  Obit *FGTable;
  /** FG Table Row (as Obit*) */
  Obit *FGTableRow;
  /** Number of rows in flag table */
  olong numRow;
  /** Last Row read */
  olong LastRowRead;
  /** Number of antennas in flag table (actually max antenna no). */
  olong numAnt;
  /** Number of subarrays in the data */
  olong numSubA;
  /** Number of IFs in data. */
  olong numIF;
  /** Number of Stokes' in data */
  olong numStok;
  /** First Stokes value in data */
  olong stoke0;
  /** Selected Subarray number. <=0-> all */
  olong SubA;
  /** Selected Frequency ID  number. <=0-> all */
  olong FreqID;
  /** Number of channels in data */
  olong numChan;
  /** Time (days) for which the flagging tables is current {TMFLST} */
  ofloat flagTime;
  /** Maximum number of flag table entries {MAXFLG} */
  olong maxFlag;
  /** Current number of flag table entries {NUMFLG} */
  olong numFlag;
  /** Source ID per flag entry {FLGSOU} */
  olong *flagSour;
  /** Antenna number per flag entry  {FLGANT} */
  olong *flagAnt;
  /** Baseline number (A1*256+A2) per flag entry  {FLGBAS} */
  olong *flagBase;
  /**  Subarray number per flag entry  {FLGSUB} */
  olong *flagSubA;
  /**  Freqid numbers per flag entry  {FLGFQD} */
  olong *flagFQID;
  /**  First IF per flag entry  {FLGBIF} */
  olong *flagBIF;
  /**  Highest IF per flag entry  {FLGEIF} */
  olong *flagEIF;
  /**  First channel per flag entry  {FLGBCH} */
  olong *flagBChan;
  /**  Highest channel per flag entry  {FLGECH} */
  olong *flagEChan;
  /**  Flags for the polarizations per flag entry  {FLGPOL} 
   Note: there are 4*numFlag values */
  gboolean *flagPol;
  /**  End time of validity per flag entry  {FLGTND} */
  ofloat *flagEndTime;
  /**  Flag all due to excessive number of flags  */
  gboolean flagAll;
  /**  Maximum number of simultaneous flags */
  olong maxSimFlag;
  /**  Number of thread argument structures in thArgArr */
  olong nThArg;
  /** Array of thread argument structures - 
      gpointer to avoid circular definition */
  gpointer thArgArr;
} ObitUVCalFlagS;
#endif /* OBITUVCALFLAGDEF_H */ 


