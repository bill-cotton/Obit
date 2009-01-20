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
#ifndef OBITOTFSEL_H 
#define OBITOTFSEL_H 
#include <glib.h>
#include "Obit.h"
#include "ObitOTFDesc.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFSel.h
 * ObitOTFSel Obit "On the Fly" data selector class definition.
 * This class is derived from the Obit class.
 * This contains the descriptions of data selection and calibration.
 *
 * \section ObitOTFSelUsage Usage
 * Instances can be obtained using the #newObitOTFSel constructor
 * the #ObitOTFSelCopy copy constructor or a pointer duplicated using 
 * the #ObitOTFSelRef function.
 * When an instance is no longer needed, use the #ObitOTFSelUnref macro
 * to release it.
 *
 * \section ObitOTFSelCalibration Data selection and Calibration
 * The ObitOTFSel member of a #ObitOTF is used to pass information the
 * the data selection and calibration routines.  This information is
 * stored on the ObitInfoList of the ObitOTF data before it is opened with 
 * access OBIT_IO_ReadCal.  Subsequent calls to ObitOTFReadSelect will apply
 * the data selection and calibration requested.  The calibration/selection 
 * paramters are described in the following list.
 * \li  "doCalSelect" OBIT_bool (1,1,1) Select/calibrate/edit data?
 * \li  "doCalib" OBIT_int (1,1,1) >0 -> calibrate,
 * \li  "gainUse" OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \li  "flagVer" OBIT_int (1,1,1) Flag table version, 0-> use highest, <0-> none
 * \li  "BChan" OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li  "EChan" OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li  "Targets" OBIT_string (?,?,1) Target names selected. [def all]
 * \li  "timeRange" OBIT_float (2,1,1) Selected timerange in days. [def all]
 * \li  "Scans" OBIT_int (2,1,1) Lowest and highest selected scan numbers. [def all]
 * \li  "Feeds" OBIT_int (?,1,1) a list of selected feed numbers, [def all.]
 * \li  "keepCal" OBIT_bool (1,1,1) If true keep cal-on data, otherwise drop [def True.]
 * \li  "replCal" OBIT_bool (1,1,1) If true replace data with cal [def False.]
*/

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOTFSel
 * returns a ObitOTFSel* (NULL).
 * \li in = object to unreference.
 */
#define ObitOTFSelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOTFSel.
 * returns a ObitOTFSel*.
 * in = object to reference
 */
#define ObitOTFSelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOTFSelIsA(in) ObitIsA (in, ObitOTFSelGetClass())
/*--------------Class definitions-------------------------------------*/
/**
 * ObitOTFSel Class structure.
 *
 * This class contains descriptions of interferometric visibility data.
 */  
typedef struct {
#include "ObitOTFSelDef.h" /* actual definition */
} ObitOTFSel;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitOTFSelClassInit (void);

/** Public: Constructor. */
ObitOTFSel* newObitOTFSel (gchar *name);

/** Public: Return class pointer. */
gconstpointer ObitOTFSelGetClass (void);

/** Public: Copy OTFSel */
ObitOTFSel* ObitOTFSelCopy (ObitOTFSel* in, ObitOTFSel* out, ObitErr *err);

/** Public: How big a buffer is needed for a data transfer? */
olong ObitOTFSelBufferSize (ObitOTFDesc* desc, ObitOTFSel* sel);

/** Public: Enforces defaults in inaxes ... */
void ObitOTFSelDefault (ObitOTFDesc* in, ObitOTFSel* sel);

/** Public: Applies selection to a Descriptor for writing */
void ObitOTFSelGetDesc (ObitOTFDesc* in, ObitOTFSel* sel,
		       ObitOTFDesc* out, ObitErr *err);

/** Public: Applies selection to a Descriptor for reading */
void ObitOTFSelSetDesc (ObitOTFDesc* in, ObitOTFSel* sel,
		       ObitOTFDesc* out, ObitErr *err);

/* Public: Initialize indexing the uv data*/
void ObitOTFSelNextInit (ObitOTFSel *in, ObitOTFDesc* desc, ObitErr *err);

/* Public: Next record to read */
gboolean ObitOTFSelNext (ObitOTFSel *in, ObitOTFDesc* desc, ObitErr *err);

/* Public: Shutdown uv data indexing */
void ObitOTFSelShutdown (ObitOTFSel *in, ObitErr *err);

/* Public: Is a given target ID selected? */
gboolean ObitOTFSelWantTarget (ObitOTFSel* sel, olong TargetID);

/* Public: Is a given Feed selected? */
gboolean ObitOTFSelWantFeed (ObitOTFSel* sel, olong feed);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOTFSelClassDef.h" /* Actual definition */
} ObitOTFSelClassInfo; 

#endif /* OBITOTFSEL_H */ 

