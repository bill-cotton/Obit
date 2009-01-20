/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
#ifndef OBITUVSEL_H 
#define OBITUVSEL_H 
#include "Obit.h"
#include "ObitUVDesc.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSel.h
 * ObitUVSel Obit uv data selector class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This contains the descriptions of data selection and calibration.
 *
 * \section ObitUVSelUsage Usage
 * Instances can be obtained using the #newObitUVSel constructor
 * the #ObitUVSelCopy copy constructor or a pointer duplicated using 
 * the #ObitUVSelRef function.
 * When an instance is no longer needed, use the #ObitUVSelUnref macro
 * to release it.
 *
 * \section ObitUVSelCalibration Data selection and Calibration
 * The ObitUVSel member of a #ObitUV is used to pass information the
 * the data selection and calibration routines.  This information is
 * stored on the ObitInfoList of the ObitUV data before it is opened with 
 * access OBIT_IO_ReadCal.  Subsequent calls to ObitUVReadSelect will apply
 * the data selection and calibration requested.  The calibration/selection 
 * paramters are described in the following list.
 * \li  "doCalSelect" OBIT_bool (1,1,1) Select/calibrate/edit data?
 * \li  "Stokes" OBIT_string (4,1,1) Selected output Stokes parameters:
 *               "    "=> no translation,"I   ","V   ","Q   ", "U   ", 
 *               "IQU ", "IQUV",  "IV  ", "RR  ", "LL  ", "RL  ", "LR  ", 
 *               "HALF" = RR,LL, "FULL"=RR,LL,RL,LR. [default "    "]
 *               In the above 'F' can substitute for "formal" 'I' (both RR+LL).
 * \li  "BChan" OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li  "EChan" OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li  "BIF"   OBIT_int (1,1,1) First "IF" selected. [def all]
 * \li  "EIF"   OBIT_int (1,1,1) Highest "IF" selected. [def all]
 * \li  "doPol"   OBIT_int (1,1,1) >0 -> calibrate polarization.
 * \li  "doCalib" OBIT_int (1,1,1) >0 -> calibrate, 2=> also calibrate Weights
 * \li  "gainUse" OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \li  "flagVer" OBIT_int (1,1,1) Flag table version, 0-> use highest, <0-> none
 * \li  "BLVer"   OBIT_int (1,1,1) BL table version, 0> use highest, <0-> none
 * \li  "BPVer"   OBIT_int (1,1,1) Band pass (BP) table version, 0-> use highest
 * \li  "Subarray" OBIT_int (1,1,1) Selected subarray, <=0->all [default all]
 * \li  "dropSubA" OBIT_bool (1,1,1) Drop subarray info?
 * \li  "FreqID"   OBIT_int (1,1,1) Selected Frequency ID, <=0->all [default all]
 * \li  "timeRange" OBIT_float (2,1,1) Selected timerange in days.
 * \li  "timeRange" OBIT_float (2,1,1) Selected timerange in days.
 * \li  "UVRange" OBIT_float (2,1,1) Selected UV range in kilowavelengths.
 * \li  "InputAvgTime" OBIT_float (1,1,1) Input data averaging time (sec).
 *                used for fringe rate decorrelation correction.
 * \li  "Sources" OBIT_string (?,?,1) Source names selected unless any starts with
 *                a '-' in which case all are deselected (with '-' stripped).
 * \li  "souCode" OBIT_string (4,1,1) Source Cal code desired, '    ' => any code selected
 *                                   '*   ' => any non blank code (calibrators only)
 *                                   '-CAL' => blank codes only (no calibrators)
 * \li  "Qual"    Obit_int (1,1,1)  Source qualifier, -1 [default] = any
 * \li  "Antennas" OBIT_int (?,1,1) a list of selected antenna numbers, if any is negative
 *                 then the absolute values are used and the specified antennas are deselected.
 * \li  "corrType" OBIT_int (1,1,1) Correlation type, 0=cross corr only, 1=both, 2=auto only.
 * \li  "passAll" OBIT_bool (1,1,1) If True, pass along all data when selecting/calibration
 *                                  even if it's all flagged, 
 *                                  data deselected by time, source, antenna etc. is not passed.
 * \li  "doBand"  OBIT_int (1,1,1) Band pass application type <0-> none
 *      (1) if = 1 then all the bandpass data for each antenna
 *          will be averaged to form a composite bandpass
 *          spectrum, this will then be used to correct the data.
 *      (2) if = 2 the bandpass spectra nearest in time (in a weighted
 *          sense) to the uv data point will be used to correct the data.
 *      (3) if = 3 the bandpass data will be interpolated in time using
 *          the solution weights to form a composite bandpass spectrum,
 *          this interpolated spectrum will then be used to correct the
 *          data.
 *      (4) if = 4 the bandpass spectra nearest in time (neglecting
 *          weights) to the uv data point will be used to correct the
 *          data.
 *      (5) if = 5 the bandpass data will be interpolated in time ignoring
 *          weights to form a composite bandpass spectrum, this
 *          interpolated spectrum will then be used to correct the data.
 * \li  "Smooth"  OBIT_float (3,1,1) specifies the type of spectral smoothing
 *         Smooth(1) = type of smoothing to apply:
 *            0 => no smoothing
 *            1 => Hanning
 *            2 => Gaussian
 *            3 => Boxcar
 *            4 => Sinc (i.e. sin(x)/x)
 *          Smooth(2) = the "diameter" of the function, i.e.
 *            width between first nulls of Hanning triangle
 *            and sinc function, FWHM of Gaussian, width of
 *            Boxcar. Defaults (if < 0.1) are 4, 2, 2 and 3
 *            channels for Smooth(1) = 1 - 4.
 *          Smooth(3) = the diameter over which the convolving
 *            function has value - in channels.
 *            Defaults: 1, 3, 1, 4 times Smooth(2) used when
 * \li "SubScanTime" Obit_float scalar [Optional] if given, this is the 
 *          desired time (days) of a sub scan.  This is used by the 
 *          selector to suggest a value close to this which will
 *          evenly divide the current scan.  See #ObitUVSelSubScan
 *          0 => Use scan average.
 *          This is only useful for ReadSelect operations on indexed ObitUVs.
 */

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVSel
 * returns a ObitUVSel* (NULL).
 * \li in = object to unreference.
 */
#define ObitUVSelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVSel.
 * returns a ObitUVSel*.
 * in = object to reference
 */
#define ObitUVSelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVSelIsA(in) ObitIsA (in, ObitUVSelGetClass())
/*--------------Class definitions-------------------------------------*/
/**
 * ObitUVSel Class structure.
 *
 * This class contains descriptions of interferometric visibility data.
 */  
typedef struct {
#include "ObitUVSelDef.h" /* actual definition */
} ObitUVSel;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVSelClassInit (void);

/** Public: Constructor. */
ObitUVSel* newObitUVSel (gchar *name);

/** Public: Return class pointer. */
gconstpointer ObitUVSelGetClass (void);

/** Public: Copy UVSel */
ObitUVSel* ObitUVSelCopy (ObitUVSel* in, ObitUVSel* out,
			  ObitErr *err);

/** Public: How big a buffer is needed for a data transfer? */
olong ObitUVSelBufferSize (ObitUVDesc* desc, 
			   ObitUVSel* sel);

/** Public: Enforces defaults in inaxes, blc, trc */
void ObitUVSelDefault (ObitUVDesc* in, 
		       ObitUVSel* sel);

/** Public: Applies selection to a Descriptor for writing */
void ObitUVSelGetDesc (ObitUVDesc* in, ObitUVSel* sel,
		       ObitUVDesc* out, ObitErr *err);

/** Public: Applies selection to a Descriptor for reading */
void ObitUVSelSetDesc (ObitUVDesc* in, ObitUVSel* sel,
		       ObitUVDesc* out, ObitErr *err);

/* Public: Initialize indexing the uv data*/
void ObitUVSelNextInit (ObitUVSel *in, ObitUVDesc* desc, ObitErr *err);

/* Public: Next visibility to read */
gboolean ObitUVSelNext (ObitUVSel *in, ObitUVDesc* desc, ObitErr *err);

/* Public: Shutdown uv data indexing */
void ObitUVSelShutdown (ObitUVSel *in, ObitErr *err);

/* Public: Setup for source selection */
void ObitUVSelSetSour (ObitUVSel* sel, gpointer inData, olong Qual, 
		       gchar *souCode, gchar *Sources, olong lsou, olong nsou, 
		       ObitErr *err);

/* Public: Setup for antenna selection */
void ObitUVSelSetAnt (ObitUVSel* sel, olong *Antennas, olong nant);

/* Public: Is a given source ID selected? */
gboolean ObitUVSelWantSour (ObitUVSel* sel, olong SourID);

/* Public: Is a given antenna ID selected? */
gboolean ObitUVSelWantAnt (ObitUVSel* sel, olong ant);

/* Public: Sugget subscan length */
ofloat ObitUVSelSubScan (ObitUVSel* sel);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVSelClassDef.h" /* Actual definition */
} ObitUVSelClassInfo; 

#endif /* OBITUVSEL_H */ 

