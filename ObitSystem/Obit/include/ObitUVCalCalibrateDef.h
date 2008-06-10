/* $Id: ObitUVCalCalibrateDef.h,v 1.2 2005/10/06 20:22:56 bcotton Exp $ */
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
#ifndef OBITUVCALCALIBRATEDEF_H 
#define OBITUVCALCALIBRATEDEF_H 

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
 * \file ObitUVCalCalibrateDef.h
 * ObitUVCal utilities for applying amp/phase/delay/rate calibration to
 * uv data 
 */

/*--------------Structure definitions-------------------------------------*/
/** Amp/phase/delay/rate calibration structure */
typedef struct {
  /** Calibration CL table (as Obit*) */
  Obit *CLTable;
  /** CL Table Row (as Obit*)*/
  Obit *CLTableRow;
  /** Calibration SN table (as Obit*)*/
  Obit *SNTable;
  /** SN Table Row (as Obit*)*/
  Obit *SNTableRow;
  /** Number of rows in calibration table */
  olong numRow;
  /** Last Row read */
  olong LastRowRead;
  /** Use SN (else CL) table? */
  gboolean doSNTable;
  /** Calibrate Weights? */
  gboolean doCalWt;
  /** Make rate smearing correction? {DORSM} */
  gboolean doRateSmear;
  /** Make delay/rate correction? */
  gboolean doDelayRate;
  /** length of entry in CalApply/CalPrior/CalFollow (LCUCAL) */
  olong lenCalArrayEntry;
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
  olong numPol;
  /** number of entries in Lambda for each IF {NLAMBDA}*/
  olong numLambda;
  /** Offset from the beginning an IF entry in the calibration array 
   * for a given polarization.  The first dimension is the polarization 
   * pixel number and the second is the antenna number of a baseline
   * (e.g. first or second = 1 or 2). 
   */
  olong PolOff[2][4];
  /** Count of the good, the bad and ugly */
  olong countRec[3][2];
  /** current source ID in cal table */
  olong CurSourID;
  /** Prior source ID in cal table */
  olong PriorSourID;
  /** Following source ID in cal table */
  olong FollowSourID;
  /** Prior time in cal table {CALTIM(1)} */
  ofloat PriorCalTime;
  /** Prior time per antenna */
  ofloat *PriorAntTime;
  /** Following Time in cal table {CALTIM(2)} */
  ofloat FollowCalTime;
  /** Following time per antenna */
  ofloat *FollowAntTime;
  /** time of calibration in CalApply {LCALTM} */
  ofloat CalTime;
  /** Integration time of data (days) */
  ofloat DeltaTime;
  /** Calibration array to apply to data {CURCAL} Values in order:
   *    By antenna {NUMANT}
   *       By IF (EIF-BIF+1)
   *          By Polarization {NUMPOL}
   *              Real part, imaginary part, cos(delta), sin(delta), rate
   *                  Where delta is the phase change between
   *                  channels and rate is the fringe rate in
   *                  radians/day
   */
  ofloat *CalApply;
  /** Prior calibration array from cal (SN or CL table) {CALTAB(*,1)} */
  ofloat *CalPrior;
  /** Following calibration array from cal (SN or CL table) {CALTAB(*,1)} */
  ofloat *CalFollow;
  /** Current Ionospheric Faraday rotation per ant {IFR} */
  ofloat *IFR;
  /** Prior Ionospheric Faraday rotation per ant {IFRTAB(*,1)} */
  ofloat *PriorIFR;
  /** Following Ionospheric Faraday rotation per IF {IFRTAB(*,2)} */
  ofloat *FollowIFR;
  /** Dispersive delay at time of current datum, one per antenna {DDELAY} */
  ofloat *DDelay;
  /** Prior dispersive delay per ant {DDTAB(*,1)} */
  ofloat *PriorDDelay;
  /** Following dispersive delay  per ant {DDTAB(*,1)} */
  ofloat *FollowDDelay;
  /** IF scaling factor to convert s/s to rad/day */
  ofloat *RateFact;
  /** IF scaling factor to convert s to rad/channel */
  ofloat *DelayFact;
  /** Array of wavelengths for each channel and IF {LAMBDA} */
  ofloat *Lambda;
  /** VLBA Filter Taper function (eg. 'HANNING') from CQ table{LTAPER, LTPVBA} 
   * \li 1 = "HANNING"
   * \li 2 = "UNIFORM"
   */
  olong *LTaper;
  /** VLBA Spectral averaging factor per IF and Corr-ID from CQ table {NSPECA,NXDSM}*/
  olong *NSpecA;
  /** VLBA Time interval per bit in FFT (s) per IF and Corr-ID from CQ table {DELBIT, DBTVBA} */
  odouble *DelBit;
  /** VLBA FFT size per IF and Corr-ID from CQ table {NFFTSZ,NFTVBA} */
  olong *NFFTSize;
  /** VLBA Time filter type (0=>boxcar; 3=>32-4; 4=>64-8 per IF and Corr-ID from CQ table 
      (ITYPTF,ITFVBA) */
  olong *typeTimeFilt;
  /** VLBA Time filter averaging time (sec) per IF and Corr-ID from CQ table {TAVGTF} */
  ofloat *TimeFiltTime;
  /** is VLBA delay decor. correction enabled per IF and Corr-ID from CQ table {WDODSM, DODSM} */
  gboolean *doDelayDecorr;
  /** VLBA  Subarray correlator type (1=VLBA; else non-VLBA) {NETVLB, ICQVBA} */
  olong *corrType;
  /** VLBA: Print warning about missing CQ per IF and Corr-ID from CQ table {WARN,WRNVBA} */
  gboolean *warnNoCQ;
} ObitUVCalCalibrateS;
#endif /* OBITUVCALCALIBRATEDEF_H */ 
