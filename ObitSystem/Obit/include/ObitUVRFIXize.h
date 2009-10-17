/* $Id:  $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
/*;Correspondence aboutthis software should be addressed as follows:  */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVRFIXIZE_H 
#define OBITUVRFIXIZE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitTableSN.h"
#include "ObitUVSoln.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVRFIXize.h
 * Obit class for RFI excision
 *
 * This class looks for signals at zero fringe rate (i.e. terresterial)
 * detects when present, estimates and subtracts from a uv dataset.
 *     This technique will estimate and remove a quasi-continuous interfering
 * signal from a uvdata set.
 * Estimates the contributions to UV data from RFI and removes them
 * RFI is estimated by counterrotating a set of residual UV data by the
 * interferometer field center fringe rate and the averaging.
 * This operation is performed on each channel/IF and polarization datastream.
 * The outUV data passed to PCreate will be filled with a corrected
 * copy of myUV.
 * Control parameters on object->info:
 * \li "solInt"    OBIT_float (1,1,1) Counter rotated SN table interval [def 1 min]
 * \li "doInvert"  OBIT_bool (1,1,1) If TRUE invert solution [def FALSE];
 * \li "timeInt"   OBIT_float (1,1,1) Data integration time in sec [def 10 sec].
 * \li "timeAvg"   OBIT_float  (1,1,1) Time interval over which to average 
 *                 residual data to estimate RFI (min) [def = 1 min.]
 *                 NB: this should be at least 2 integrations.
 *                 This is the time interval on which RFI will be searched
 * \li "minRFI"    OBIT_float (1,1,1) Minimum RFI amplitude (Jy) [def 50]
 *
 * \section ObitUVRFIXizeaccess Creators and Destructors
 * An ObitUVRFIXize will usually be created using ObitUVRFIXizeCreate which allows 
 * specifying a name for the object as well as other information.

 * This operation involves three steps:
 * \li Counter-rotate/average residual \\
 *  The residual data is counter-rotated by the field center fringe rate and then
 * averaged to timeAvg min and written to scratch file RFIUV.
 * \li Filter to estimate RFI \\
 * The counter-rotated and averaged residuals are filtered to produce a file which 
 * contains the estimated of the RFI contribution to the signal.
 * \li Correct input data by estimated RFI \\
 * The closest RFI visibility estimate to each visibility in myUV is subtracted to
 * write outUV.
 *
 * A copy of a pointer to an ObitUVRFIXize should always be made using the
 * #ObitUVRFIXizeRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVRFIXize or changing a pointer, the function
 * #ObitUVRFIXizeUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/*--------------Class definitions-------------------------------------*/
/** ObitUVRFIXize Class structure. */
typedef struct {
#include "ObitUVRFIXizeDef.h"   /* this class definition */
} ObitUVRFIXize;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVRFIXize
 * returns a ObitUVRFIXize*.
 * in = object to unreference
 */
#define ObitUVRFIXizeUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVRFIXize.
 * returns a ObitUVRFIXize*.
 * in = object to reference
 */
#define ObitUVRFIXizeRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVRFIXizeIsA(in) ObitIsA (in, ObitUVRFIXizeGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVRFIXizeClassInit (void);

/** Public: Default Constructor. */
ObitUVRFIXize* newObitUVRFIXize (gchar* name);

/** Public: Create/initialize ObitUVRFIXize structures */
ObitUVRFIXize* ObitUVRFIXizeCreate (gchar* name, ObitUV *inUV, ObitUV *residUV, 
				    ObitUV *outUV);
/** Typedef for definition of class pointer structure */
typedef ObitUVRFIXize* (*ObitUVRFIXizeCreateFP) (gchar* name, ObitUV *inUV, 
						 ObitUV *residUV, ObitUV *outUV);

/** Public: ClassInfo pointer */
gconstpointer ObitUVRFIXizeGetClass (void);

/** Public: Counterrotate and average residual data */
void ObitUVRFIXizeCounterRot (ObitUVRFIXize *in, ObitErr *err);
typedef void (*ObitUVRFIXizeCounterRotFP) (ObitUVRFIXize *in, ObitErr *err);

/** Public: Filter Counterrotated/averaged residual data */
void ObitUVRFIXizeFilter (ObitUVRFIXize *in, ObitErr *err);
typedef void (*ObitUVRFIXizeFilterFP) (ObitUVRFIXize *in, ObitErr *err);

/** Public: Remove estimated RFI from data */
void ObitUVRFIXizeCorrect (ObitUVRFIXize *in, ObitErr *err);
typedef void (*ObitUVRFIXizeCorrectFP) (ObitUVRFIXize *in, ObitErr *err);

/** Public: Initialize Fetching filtered data */
void ObitUVRFIXizeFetchStartUp (ObitUVRFIXize *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitUVRFIXizeFetchStartUpFP) (ObitUVRFIXize *in, ObitErr *err);

/** Public: interpolate calibration at a given time */
gboolean ObitUVRFIXizeFetch (ObitUVRFIXize *in, ofloat time, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitUVRFIXizeFetchFP) (ObitUVRFIXize *in, ofloat time, 
					  ObitErr *err);

/** Public: Shutdown interpolation */
void ObitUVRFIXizeFetchShutDown (ObitUVRFIXize *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitUVRFIXizeFetchShutDownFP) (ObitUVRFIXize *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVRFIXizeClassDef.h"
} ObitUVRFIXizeClassInfo; 

#endif /* OBITFUVRFIXIZE_H */ 
