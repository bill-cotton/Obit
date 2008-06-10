/* $Id: ObitUVCalBaselineDef.h,v 1.2 2005/10/06 20:22:56 bcotton Exp $  */
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
#ifndef OBITUVCALBASELINEDEF_H 
#define OBITUVCALBASELINEDEF_H 

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
 * \file ObitUVCalBaselineDef.h
 * ObitUVCal utilities for applying baseline dependent calibration to
 * uv data 
 */

/*--------------Structure definitions-------------------------------------*/
/** Amp/phase/delay/rate calibration structure */
typedef struct {
  /** Calibration BL table (as Obit*) */
  Obit *BLTable;
  /** BL Table Row (as Obit*)*/
  Obit *BLTableRow;
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
  /** Number of baselines */
  olong numBase;
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
  olong numPol;
  /** Offset from the beginning an IF entry in the calibration array 
   * for a given polarization.  The first dimension is the polarization 
   * pixel number and the second is the antenna number of a baseline
   * (e.g. first or second = 1 or 2). 
   */
  olong PolOff[2][4];
  /** current source ID in cal table */
  olong CurSourID;
  /** Prior source ID in cal table */
  olong PriorSourID;
  /** Following source ID in cal table */
  olong FollowSourID;
  /** Length of calibration array entry */
  olong lenCalArrayEntry;
  /** Calibrate Weights? */
  gboolean doCalWt;
  /** Integration time of data */
  ofloat DeltaTime;
  /** Prior time in cal table {BLTIM(1)} */
  ofloat PriorCalTime;
  /** Following Time in cal table {BLLTIM(2)} */
  ofloat FollowCalTime;
  /** time of calibration in CalApply {LCALTM} */
  ofloat CalTime;
  /** Calibration array to apply to data {BLFAC} Values in order:
   *    Indexing scheme: an entry defined by ant1<ant2 starts in element:
   *    lentry * (((ant1-1)*numant-((ant1+1)*ant1)/2 + ant2) - 1) + 1
   *      where lentry = 2 * NUMPOL * (EIF-BIF+1)
   *         An entry contains the values in order:
   *           By IF (NUMIF)
   *             By Polarization (NUMPOL)
   *               Real part, imaginary part.
   *               Applied only to cross corelation data.
   */
  ofloat *CalApply;
  /** Prior calibration array from cal (BL table) {BLLTAB(*,1)} */
  ofloat *CalPrior;
  /** Following calibration array from cal (BL table) {BLLTAB(*,1)} */
  ofloat *CalFollow;
} ObitUVCalBaselineS;
#endif /* OBITUVCALBASELINEDEF_H */ 


