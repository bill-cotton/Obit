/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVGSOLVE_H 
#define OBITUVGSOLVE_H 

#include "ObitUV.h"
#include "ObitTableSN.h"
#include "ObitSkyModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGSolve.h
 * ObitUVGSolve Class for uv data self calibration
 *
 * This class is derived from the #Obit class.
 *
 * Routines determine calibration for an ObitUV and write an SN table.
 * Control parameters are on the info member.
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "solInt"  OBIT_float (1,1,1) Solution interval (min). (default scan)
 * \li "refAnt"  OBIT_int   (1,1,1) Ref ant to use. (default 1)
 * \li "avgPol"  OBIT_bool  (1,1,1) True if RR and LL to be averaged (false)
 * \li "avgIF"   OBIT_bool  (1,1,1) True if all IFs to be averaged (false)
 * \li "minSNR"  OBIT_float (1,1,1) Minimum acceptable SNR (5)
 * \li "doMGM"   OBIT_bool  (1,1,1) True then find the mean gain modulus (true)
 * \li "elevMGM" OBIT_float (1,1,1) Min. elevation to include in mean gain modulus
 * \li "solType" OBIT_string (4,1,1 Solution type '  ', 'L1',  (' ')
 * \li "solMode" OBIT_string (4,1,1 Solution mode: 'A&P', 'P', 'P!A', 'GCON' ('P')
 * \li "minNo"   OBIT_int   (1,1,1) Min. no. antennas. (default 4)
 * \li "antWt"   OBIT_float (*,1,1) Antenna weights. (default 1.0)
 * \li "UVR_Full"OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (kilolamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "prtLv"   OBIT_int   (1,1,1) Print level (default no print)
 * 
 * \section ObitUVGSolveaccess Creators and Destructors
 * An ObitUVGSolve will usually be created using ObitUVGSolveCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitUVGSolve should always be made using the
 * #ObitUVGSolveRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVGSolve or changing a pointer, the function
 * #ObitUVGSolveUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUVGSolve Class structure. */
typedef struct {
#include "ObitUVGSolveDef.h"   /* this class definition */
} ObitUVGSolve;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVGSolve
 * returns a ObitUVGSolve*.
 * in = object to unreference
 */
#define ObitUVGSolveUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVGSolve.
 * returns a ObitUVGSolve*.
 * in = object to reference
 */
#define ObitUVGSolveRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVGSolveIsA(in) ObitIsA (in, ObitUVGSolveGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVGSolveClassInit (void);

/** Public: Default Constructor. */
ObitUVGSolve* newObitUVGSolve (gchar* name);

/** Public: Create/initialize ObitUVGSolve structures */
ObitUVGSolve* ObitUVGSolveCreate (gchar* name);
/** Typedef for definition of class pointer structure */
typedef ObitUVGSolve* (*ObitUVGSolveCreateFP) (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitUVGSolveGetClass (void);

/** Public: Copy (deep) constructor. */
ObitUVGSolve* ObitUVGSolveCopy  (ObitUVGSolve *in, ObitUVGSolve *out, 
				   ObitErr *err);

/** Public: Copy structure. */
void ObitUVGSolveClone (ObitUVGSolve *in, ObitUVGSolve *out, ObitErr *err);

/** Public: Determine calibration from UV data divided by model. */
ObitTableSN* ObitUVGSolveCal (ObitUVGSolve *in, ObitUV *inUV, ObitUV *outUV, 
			      ObitUVSel *sel, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVGSolveClassDef.h"
} ObitUVGSolveClassInfo; 

#endif /* OBITFUVGSOLVE_H */ 
