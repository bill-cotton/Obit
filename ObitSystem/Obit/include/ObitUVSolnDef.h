/* $Id: ObitUVSolnDef.h,v 1.2 2006/12/27 17:32:07 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUVSoln structure         */
/*  This is intended to be included in a class structure definition   */
/* and to be used as the template for generating new classes derived  */
/* from Obit.                                                         */
/**
 * \file ObitUVSolnDef.h
 * ObitUVSoln structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** UV data object associated with input Solution table */
ObitUV *myUV;
/** Calibration SN table */
ObitTableSN *SNTable;
/** SN Table Row */
ObitTableSNRow *SNTableRow;
/** Calibrator selector */
ObitUVSel *CalSel;
/** is SN table smoothed? (to be deleted) */
gboolean isSNSmoo;
/** Number of rows in calibration table */
olong numRow;
/** Last Row read */
olong LastRowRead;
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
/** Number of polarizations in the calibration table (1 or 2) */
olong numPol;
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
 *              Real part, imaginary part, delay, rate, weight, refant
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
/** Current Multiband delay per ant. and poln.  */
ofloat *MBDelay;
/** Prior Multiband delay per ant. and poln.  */
ofloat *PriorMBDelay;
/** Following Multiband delay per ant. and poln. */
ofloat *FollowMBDelay;
/** IF scaling factor to convert s/s to rad/day */
ofloat *RateFact;
/** Current Reference antenna per ant and poln */
olong *RefAnt;
/** Interpolation mode */
ObitUVSolnInterMode interMode;
/** Interpolation parameters */
ofloat interParm[10];
/** Max. time over which to interpolate */
ofloat maxInter;
/** number of terms in interpolation polynomial */
olong interNPoly;
/** Missing Antennas (no solutions) */
olong *MissAnt;

