/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2011                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUVCal structure            */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitUVCalDef.h
 * ObitUVCal structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** I/O status */
ObitIOStatus myStatus;
/** Obit UV data Descriptor */
ObitUVDesc* myDesc;
/** Obit UV data Selector */
ObitUVSel* mySel;
/** Drop Subarray info? */
gboolean dropSubA;
/** Apply flagging? */
gboolean doFlag;
/** Apply Spectral smoothing*/
gboolean doSmo;
/** Apply Baseline based calibration? */
gboolean doBL;
/** Apply amp/phase/delay/rate/baseline calibration */
gboolean doCal;
/** Apply Bandpass calibration? if > 0 specifies bandpass type */
olong doBP;
/** Apply Polarization calibration? */
gboolean doPol;
/** Stokes translation indices */
olong jadr[4][2];
/** Number of antennas in calibration table (actually max antenna no). */
olong numAnt;
/** Number of subarrays in the data */
olong numSubA;
/** Number of Spectral channels */
olong numChan;
/** Start channel number (1-rel) */
olong bChan;
/** highest channel (1-rel) */
olong eChan;
/** Number of IFs */
olong numIF;
/** Start IF number (1-rel) selected */
olong bIF;
/** highest IF (1-rel) selected */
olong eIF;
/** Number of Stokes in data*/
olong numStok;
/** Stokes translation factors */
ofloat selFact[4][2];
/** Spectral smoothing
 *         SMOOTH(1) = type of smoothing to apply:
 *            0 => no smoothing
 *            1 => Hanning
 *            2 => Gaussian
 *            3 => Boxcar
 *            4 => Sinc (i.e. sin(x)/x)
 *          SMOOTH(2) = the "diameter" of the function, i.e.
 *            width between first nulls of Hanning triangle
 *            and sinc function, FWHM of Gaussian, width of
 *            Boxcar. Defaults (if < 0.1) are 4, 2, 2 and 3
 *            channels for SMOOTH(1) = 1 - 4.
 *          SMOOTH(3) = the diameter over which the convolving
 *            function has value - in channels.
 *            Defaults: 1, 3, 1, 4 times SMOOTH(2) used when
 */
ofloat smooth[3];
/** Spectral index to apply to data (wrt alphaRefF) */
ofloat alpha;
/** Reference frequency (Hz) for spectral index */
odouble alphaRefF;
/** Start channel number for smoothing (1-rel) {BCHANS}*/
olong bChanSmo;
/** highest channel for smoothing (1-rel) {ECHANS} */
olong eChanSmo;
/** Width of spectral Smoothing kernal */
olong SmoothWidth;
/** Smoothing convolution function {SMTAB} */
ofloat *SmoothConvFn;
/** Spectral smoothing work array */
ofloat *SmoothWork;
/** Spectral Index work array */
ofloat *SpecIndxWork;
/** SourceList with source info */
ObitSourceList *sourceList;
/** Array of AntennaLists one per subarray */
ObitAntennaList **antennaLists;
/** data flagging structure */
ObitUVCalFlagS *flag;
/** Baseline dependent calibration structure */
ObitUVCalBaselineS *baselineCal;
/** Amp/phase/delay/rate calibration structure */
ObitUVCalCalibrateS *ampPhaseCal;
/** Bandpass calibration structure */
ObitUVCalBandpassS *bandpassCal;
/** Polarization calibration structure */
ObitUVCalPolarizationS *polnCal;
/** Number of subarrays (number of entries in ANtables) */
olong numANTable;
/** Array of AN tables for polarization calibration, one per subarray (as Obit*)*/
Obit **ANTables;
/** Baseline dependent calibration table (as Obit*) */
Obit *BLTable;
/** Bandpass calibration table (as Obit*) */
Obit *BPTable;
/** VLBA corelator parameter table */
Obit *CQTable;
/** CL table for multisource calibration */
Obit *CLTable;
/* Flag table */
Obit *FGTable;
/* Solution table for single source calibration */
Obit *SNTable;
/* Source  table, NULL -> none */
Obit *SUTable;

