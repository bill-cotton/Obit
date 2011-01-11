/* $Id$ */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUVSel structure            */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitUVSelDef.h
 * ObitUVSel structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** File Type requested (FITS/AIPS) */
ObitIOType FileType;
/** Max. number of visibilities per read/write */
olong nVisPIO;
/** Number of visibilities for next read */
olong numVisRead;
/** Size of uncompressed record */
olong lrecUC;
/** Number of random parameters in the uncompressed version */
olong nrparmUC;
/** Number of visibilities in output */
olong numberVis;
/** Number of Stokes parameters in output */
olong numberPoln;
/** Increment in visibility float array between Stokes */
olong jincs;
/** Start channel (1-rel) */
olong startChann;
/** Number of channels */
olong numberChann;
/** Increment in visibility float array between Frequencies */
olong jincf;
/** Start IF (1-rel) */
olong startIF;
/** Number of IFs */
olong numberIF;
/** Increment in visibility float array between IFs */
olong jincif;
/** Selected Subarray number. <=0 -> all */
olong SubA;
/** Selected Frequency ID  number. <=0 -> all */
olong FreqID;
/** Calibrate/edit/select Data? */
gboolean doCalSelect;
/** Input data compressed? */
gboolean Compress;
/** Translate Stokes? */
gboolean transPol;
/** Need both correlations to form output Stokes parameter? */
gboolean bothCorr;
/** Select Antennas? (deselect if FALSE) */
gboolean selectAnts;
/** Number of entries in ants 0=> all selected. */
olong numberAntList;
/** List of selected (or deselected) antennas, NULL => all selected */
olong *ants;
/** Select Sources? (deselect if FALSE) */
gboolean selectSources;
/** Number of entries in sources, 0=> all selected. */
olong numberSourcesList;
/** List of selected (or deselected) source ids, NULL => all selected */
olong *sources;
/** Selected Stokes parameter(s) */
gchar Stokes[5];
/** Start and end times in days */
ofloat timeRange[2];
/** UV range (lower, upper) in wavelengths at reference frequency */
ofloat UVRange[2];
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
/** do Polarization Calibration? */
gboolean doPolCal;
/** do Baseline dependent calibration? */
gboolean doBLCal;
/** version number of BL table for Baseline dependent calibration */
olong BLversion;
/** do Bandpass calibration? */
gboolean doBPCal;
/** Bandpass calibration type, <=0 -> none,
 *      (1) if = 1 then all the bandpass data for each antenna
 *          will be averaged to form a composite bandpass
 *          spectrum, this will then be used to correct the data.
 *      (2) if = 2 the bandpass spectra nearest in time (in a weighted
 *          sense) to the uv data point will be used to correct the data.
 *      (3) if =3 the bandpass data will be interpolated in time using
 *          the solution weights to form a composite bandpass spectrum,
 *          this interpolated spectrum will then be used to correct the
 *          data.
 *      (4) if = 4 the bandpass spectra nearest in time (neglecting
 *          weights) to the uv data point will be used to correct the
 *          data.
 *      (5) if = 5 the bandpass data will be interpolated in time ignoring
 *          weights to form a composite bandpass spectrum, this
 *          interpolated spectrum will then be used to correct the data.
 */
olong doBand;
/** version number of BP table for Bandpass calibration */
olong BPversion;
/** Desired correlation type, 0= Cross only, 1=both, 2= Auto only */
olong corrType;
/** do amp/phase/delay/rate calibration? */
gboolean doCal;
/** do amp/phase/delay/rate calibration of weights? */
gboolean doCalWt;
/** Drop Subarray info? */
gboolean dropSubA;
/** version number of SN/CL table for  calibration */
olong calVersion;
/** apply flagging? */
gboolean doFlag;
/** version number of FG table for flagging */
olong FGversion;
/** Using an index table? */
gboolean doIndex;
/** NX table (as Obit*) */
Obit *NXTable;
/** NX Table Row (as Obit*) */
Obit *NXTableRow;
/** Number of rows in flag table */
olong numRow;
/** Last Row read */
olong LastRowRead;
/** First visibility in scan */
olong scanFirstVis;
/** Last visibility in scan */
olong scanLastVis;
/** Desired length in days of subscans */
ofloat SubScanTime;
/** Interval in current scan leading to even number of subscans */
ofloat SubScanSuggest;
/** Input data averaging time for fringe rate decorrelation calc. */
ofloat InputAvgTime;
/** Keep all data whether flagged or not */
gboolean passAll;
