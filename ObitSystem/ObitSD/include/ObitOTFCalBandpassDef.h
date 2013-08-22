/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
#ifndef OBITOTFCALBANDPASSDEF_H 
#define OBITOTFCALBANDPASSDEF_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"
#include "ObitFFT.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalBandpassDef.h
 * ObitOTFCal utilities for applying bandpass calibration to uv data 
 */

/*--------------Structure definitions-------------------------------------*/
/** Amp/phase/delay/rate calibration structure */
typedef struct {
  /** Calibration OTFBP table (as Obit*) */
  Obit *BPTable;
  /** OTFBP Table Row (as Obit*)*/
  Obit *BPTableRow;
  /** Number of rows in calibration table */
  olong numRow;
  /** Last Row read */
  olong LastRowRead;
  /** Number of feeds in data */
  olong numFeed;
  /** Number of channels in data */
  olong numChan;
  /** Start channel number (1-rel) */
  olong bChan;
  /** highest channel (1-rel) */
  olong eChan;
  /** Number of polarizations in the calibration table (1 or 2) */
  olong numPol;
  /** Offset from the beginning an feed entry in the calibration array 
   * for a given polarization.  The first dimension is the polarization 
   * pixel number and the second is the antenna number of a baseline
   * (e.g. first or second = 1 or 2). 
   */
  olong PolOff[2][4];
  /** Apply Bandpass calibration?  */
  gboolean doBand;
  /** Calibrate Weights? */
  gboolean doBPWt;
  /** Multiplicative calibration array to apply to data Values in order:
   *    channel, poln, feed
   */
  ofloat *BPApply;
} ObitOTFCalBandpassS;
#endif /* OBITOTFCALBANDPASSDEF_H */ 


