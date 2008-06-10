/* $Id: ObitUVCalBandpassDef.h,v 1.2 2005/10/06 20:22:56 bcotton Exp $  */
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
#ifndef OBITUVCALBANDPASSDEF_H 
#define OBITUVCALBANDPASSDEF_H 

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
#include "ObitFFT.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalBandpassDef.h
 * ObitUVCal utilities for applying bandpass calibration to uv data 
 */

/*--------------Structure definitions-------------------------------------*/
/** Amp/phase/delay/rate calibration structure */
typedef struct {
  /** Calibration BP table (as Obit*) */
  Obit *BPTable;
  /** BP Table Row (as Obit*)*/
  Obit *BPTableRow;
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
  /** Type of bandpass solution 1->cross power, 2->total power, 3->mix, default 1 */
  olong BPType;
  /** Length of calibration array entry */
  olong lenBPArrayEntry;
  /** Bandpass calibration type 
   * \li if = 1 then all the bandpass data for each antenna  will be averaged to form a 
   *      composite bandpass spectrum, this will then be used to correct the data.
   * \li if = 2 the bandpass spectra nearest in time (in a weighted  sense) to the uv 
   *      data point will be used to correct the data.
   * \li if = 3 the bandpass data will be interpolated in time using the solution weights to 
   *      form a composite bandpass spectrum, this interpolated spectrum will then be used 
   *      to correct the data.
   * \li if = 4 the bandpass spectra nearest in time (neglecting weights) to the uv data point 
   *      will be used to correct the data.
   * \li if = 5 the bandpass data will be interpolated in time ignoring weights to form a 
   *      composite bandpass spectrum, this interpolated spectrum will then be used to 
   * \li if = 1 then all the bandpass data for each antenna  will be averaged to form a 
   *      composite bandpass spectrum, this will then be used to correct the data.
   * \li if = 2 the bandpass spectra nearest in time (in a weighted  sense) to the uv 
   *      data point will be used to correct the data.
   * \li if = 3 the bandpass data will be interpolated in time using the solution weights to 
   *      form a composite bandpass spectrum, this interpolated spectrum will then be used 
   *      to correct the data.
   * \li if = 4 the bandpass spectra nearest in time (neglecting weights) to the uv data point 
   *      will be used to correct the data.
   * \li if = 5 the bandpass data will be interpolated in time ignoring weights to form a 
   *      composite bandpass spectrum, this interpolated spectrum will then be used to 
   *      correct the data.
   */
  olong doBand;
  /** if this VLBA data? per subarray */
  gboolean *isVLBA;
  /** Calibrate Weights? */
  gboolean doBPWt;
  /** Integration time of data */
  ofloat DeltaTime;
  /** Prior time in cal table {BPTIM(1)} */
  ofloat PriorBPTime;
  /** Prior time per antenna */
  ofloat *PriorAntTime;
  /** Following Time in cal table {BPLTIM(2)} */
  ofloat FollowBPTime;
  /** time of calibration in CalApply {LCALTM} */
  ofloat BPTime;
  /** Following time per antenna */
  ofloat *FollowAntTime;
  /** Calibration array to apply to data {BLFAC} Values in order:
   *    Indexing scheme: an entry defined by ant1<ant2 starts in element:
   *    lentry * (((ant1-1)*numant-((ant1+1)*ant1)/2 + ant2) - 1) + 1
   *      where lentry = 2 * NUMPOL * (EIF-BIF+1)
   *         An entry contains the values in order:
   *           By IF (NUMIF)
   *             By Polarization (NUMPOL)
   *                By Channel
   *                  Real part, 
   *                  imaginary part.
   */
  ofloat *BPApply;
  /** Prior calibration array from cal (BP table) {BLLTAB(*,1)} */
  ofloat *BPPrior;
  /** Following calibration array from cal (BP table) {BLLTAB(*,1)} */
  ofloat *BPFollow;
  /** Solution Weights, per, antenna, per selected IF,  per poln, prior solution */
  ofloat *PriorSolnWt;
  /** Solution Weights, per, antenna, per selected IF,  per poln, prior solution */
  ofloat *FollowSolnWt;
  /** Spectra work arrays 2*nchan complex arrays */
  ObitCArray *spectrum1, *spectrum2;
  /** ObitFFTs for autocorrelations */
  ObitFFT *ACFFTFor, *ACFFTRev;
  /** ObitFFTs for crosscorrelations */
  ObitFFT *CCFFTFor, *CCFFTRev;
  /** Work arrays for expanding polynomial BP table entries */
  ofloat *BPWork1, *BPWork2, *BPWork3, *BPWork4;
} ObitUVCalBandpassS;
#endif /* OBITUVCALBANDPASSDEF_H */ 


