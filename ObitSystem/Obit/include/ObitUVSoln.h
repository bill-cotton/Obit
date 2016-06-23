/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2016                                          */
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
#ifndef OBITUVSOLN_H 
#define OBITUVSOLN_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitTableSN.h"
#include "ObitTableCL.h"
#include "ObitTableBP.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSoln.h
 * Obit utilities for manipulating Solution (SN) tables
 *
 * This class manipulates, mostly interpolates, solutions tables to 
 * specified times and using a vareity of interpolation techniques.
 *
 * The following options can be entered onto the info list
 * prior to StartUp:
 * \li "interMode", OBIT_string (4,1,1)  Interpolation mode:
                    default or blank = "2PT "
 * \li "2PT " = linear vector interpolation with no SN smoothing.
 * \li "SELF" = Use only SN solution from same source which is closest in time.
 * \li "POLY" = Fit a polynomial to the SN rates and delays.
 *             Use the integral of the rate polynomial for the phases.
 * \li "SIMP" = Simple linear phase connection between SN phase
 *              entries, assumes phase difference less than 180 degrees.
 * \li "AMBG" = Linear phase connection using rates to resolve phase ambiguities.
 * \li "CUBE" = As AMBG but fit third order polynomial to phases and rates.
 * \li "MWF " = Median window filter of SN table before 2PT interpolation
 * \li "BOX " = Boxcar smoothing of SN table before 2PT interpolation,
 *              boxcar width set by adverb INTPARM.
 * \li "GAUS" = Gaussian smoothing of SN table before 2PT interpolation,
 * \li "interParm", OBIT_float (3,1,1) interpolation parameters
 * \li mode="BOX ", smoothing time in hours for amplitude, phase, delay/rate
 * \li mode="MWF ", window size in hours for amplitude, phase, delay/rate
 * 
 * \li "interNPoly", OBIT_int (1,1,1) number of terms in polynomial for
 *                   mode="POLY", default = 2
 *
 *
 * \li "maxInter", OBIT_float (1,1,1) Max. time (min) over which to interpolate.
 *                 default = 1440.0;
 *
 * \li "solnVer", OBIT_int (1,1,1) Solution (SN) table version to use.
 *                0=> highest, default 0
 *
 * Any data selection parameters on the input UV data info object will 
 * be applied.
 *
 * \section ObitUVSolnaccess Creators and Destructors
 * An ObitUVSoln will usually be created using ObitUVSolnCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitUVSoln should always be made using the
 * #ObitUVSolnRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVSoln or changing a pointer, the function
 * #ObitUVSolnUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitUVSolnInterMode
 * enum for UVSoln interpolation modes
 */
enum obitUVSolnInterMode {
  /** Undefined interpolation type */
  OBIT_UVSolnInterUnknown=0, 
  /** linear vector interpolation with no SN smoothing. */
  OBIT_UVSolnInter2PT, 
  /** Use only SN solution from same source */
  OBIT_UVSolnInterSELF, 
  /** Fit a polynomial to the SN rates and delays */
  OBIT_UVSolnInterPOLY, 
  /** Simple linear phase connection  */
  OBIT_UVSolnInterSIMP, 
  /** Linear phase connection using rates to resolve phase ambiguities. */
  OBIT_UVSolnInterAMBG, 
  /** As AMBG but fit third order polynomial to phases and rates. */
  OBIT_UVSolnInterCUBE, 
  /** Median window filter of SN table before 2PT interpolation */
  OBIT_UVSolnInterMWF, 
  /** Boxcar smoothing of SN table before 2PT interpolation */
  OBIT_UVSolnInterBOX,
  /** Gaussian smoothing of SN table before 2PT interpolation */
  OBIT_UVSolnInterGAUS
}; /* end enum obitIOStatus */
/** typedef for enum for UVSoln interpolation modes */
typedef enum obitUVSolnInterMode ObitUVSolnInterMode; 

/*--------------Class definitions-------------------------------------*/
/** ObitUVSoln Class structure. */
typedef struct {
#include "ObitUVSolnDef.h"   /* this class definition */
} ObitUVSoln;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVSoln
 * returns a ObitUVSoln*.
 * in = object to unreference
 */
#define ObitUVSolnUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVSoln.
 * returns a ObitUVSoln*.
 * in = object to reference
 */
#define ObitUVSolnRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVSolnIsA(in) ObitIsA (in, ObitUVSolnGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVSolnClassInit (void);

/** Public: Default Constructor. */
ObitUVSoln* newObitUVSoln (gchar* name);

/** Public: Create/initialize ObitUVSoln structures */
ObitUVSoln* ObitUVSolnCreate (gchar* name, ObitUV *inUV);
/** Typedef for definition of class pointer structure */
typedef ObitUVSoln* (*ObitUVSolnCreateFP) (gchar* name, ObitUV *inUV);

/** Public: ClassInfo pointer */
gconstpointer ObitUVSolnGetClass (void);

/** Public: Copy (deep) constructor. */
ObitUVSoln* ObitUVSolnCopy  (ObitUVSoln *in, ObitUVSoln *out, ObitErr *err);

/** Public: Copy structure. */
void ObitUVSolnClone (ObitUVSoln *in, ObitUVSoln *out, ObitErr *err);

/** Public: Initialize interpolation */
void ObitUVSolnStartUp (ObitUVSoln *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitUVSolnStartUpFP) (ObitUVSoln *in, ObitErr *err);

/** Public: interpolate calibration at a given time */
gboolean ObitUVSolnGetSN (ObitUVSoln *in, ObitTableSNRow *SNrow, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitUVSolnGetSNFP) (ObitUVSoln *in, 
					 ObitTableSNRow *SNrow, ObitErr *err);

/** Public: Shutdown interpolation */
void ObitUVSolnShutDown (ObitUVSoln *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitUVSolnShutDownFP) (ObitUVSoln *in, ObitErr *err);

/** Update calibration arrays. */
void ObitUVSolnUpdate (ObitUVSoln *in, ofloat time, olong SourID,
		       ObitErr *err);

/** Refererence phases to a common reference antenna */
void ObitUVSolnRefAnt (ObitTableSN *SNTab, olong isuba, olong* refant, 
		       ObitErr* err);
/** Average fringe rates over IF/poln in an SN table. */
void ObitUVSolnAvgRate (ObitTableSN *SNTab, ObitErr* err); 

/** Smooth an SN table and possible interpolate blanked soln. */
void ObitUVSolnSNSmo (ObitTableSN *SNTab, olong isuba, ObitErr* err);

/** Deselect records in an SN table */
void ObitUVSolnDeselSN (ObitTableSN *SNTab, olong isuba, olong fqid, 
			olong nantf, olong *ants, olong nsou, olong *sources,
			ofloat timerange[2], ObitErr* err);

/** Deselect records in a CL table */
void ObitUVSolnDeselCL (ObitTableCL *CLTab, olong isuba, olong fqid, 
			olong nantf, olong *ants, olong nsou, olong *sources,
			ofloat timerange[2], ObitErr* err);

/** Deselect records in a BP table */
void ObitUVSolnDeselBP (ObitTableBP *BPTab, olong isuba, olong fqid, 
			olong nantf, olong *ants, olong nsou, olong *sources,
			ofloat timerange[2], ObitErr* err);

/** Smooth solutions for a given IF and subarray  */
void 
ObitUVSolnSNSmooth (ObitTableSN *SNTab, gchar* smoFunc, gchar* smoType, ofloat alpha, ofloat *smoParm, 
		    olong iif, olong sub, ofloat* gncnt, ofloat* gnsum, 
		    olong nxt, ofloat* work1, ofloat* work2, gboolean doBlank, ObitErr* err);

/** Boxcar smoothing with weighting of an irregularly spaced array */
void 
ObitUVSolnSmooBox (ofloat smoTime, ofloat* x, ofloat* y, ofloat *w, olong n, 
		   ofloat* ys, ofloat* ws, gboolean doBlank);

/** Gaussian smoothing with weighting of an irregularly spaced array */
void 
ObitUVSolnSmooGauss (ofloat smoTime, ofloat* x, ofloat* y, ofloat *w, olong n, 
		     ofloat* ys, ofloat* ws, ofloat* wtsum, gboolean doBlank);

/** Median Window smoothing with weighting of an irregularly spaced array */
void 
ObitUVSolnSmooMWF (ofloat smoTime, ofloat alpha, ofloat* x, ofloat* y, ofloat *w, 
		   olong n, ofloat* ys, ofloat* ws, ofloat* yor, ofloat* wor, gboolean doBlank);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVSolnClassDef.h"
} ObitUVSolnClassInfo; 

#endif /* OBITFUVSOLN_H */ 
