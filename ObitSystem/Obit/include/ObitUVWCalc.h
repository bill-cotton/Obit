/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVWCALC_H 
#define OBITUVWCALC_H 

#include "Obit.h"
#include "ObitUV.h"
#include "ObitSourceList.h"
#include "ObitAntennaList.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVWCalc.h
 *
 * ObitUVWCalc calculates [u,v,w] vectors for an ObitUV
 * Uses short baseline approximations and NO relativistic corrections
 * 
 * \section ObitUVWCalcaccess Creators and Destructors
 * An ObitUVWCalc will usually be created using ObitUVWCalcCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitUVWCalc should always be made using the
 * #ObitUVWCalcRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVWCalc or changing a pointer, the function
 * #ObitUVWCalcUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUVWCalc Class structure. */
typedef struct {
#include "ObitUVWCalcDef.h"   /* this class definition */
} ObitUVWCalc;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVWCalc
 * returns a ObitUVWCalc*.
 * in = object to unreference
 */
#define ObitUVWCalcUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVWCalc.
 * returns a ObitUVWCalc*.
 * in = object to reference
 */
#define ObitUVWCalcRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVWCalcIsA(in) ObitIsA (in, ObitUVWCalcGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVWCalcClassInit (void);

/** Public: Default Constructor. */
ObitUVWCalc* newObitUVWCalc (gchar* name);

/** Public: Create/initialize ObitUVWCalc structures */
ObitUVWCalc* ObitUVWCalcCreate (gchar* name, ObitUV *inUV, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitUVWCalc* (*ObitUVWCalcCreateFP) (gchar* name, ObitUV *inUV, 
					     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitUVWCalcGetClass (void);

/** Public: Copy (deep) constructor. */
ObitUVWCalc* ObitUVWCalcCopy  (ObitUVWCalc *in, ObitUVWCalc *out, ObitErr *err);

/** Public: Copy structure. */
void ObitUVWCalcClone (ObitUVWCalc *in, ObitUVWCalc *out, ObitErr *err);

/** Public: Calculate u,v,w  */
void ObitUVWCalcUVW (ObitUVWCalc *in, ofloat time, olong SId,
		     olong subA, olong ant1, olong ant2, ofloat *uvw, 
		     ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitUVWCalcUVWFP) (ObitUVWCalc *in, ofloat time, olong SId,
				  olong subA, olong ant1, olong ant2, 
				  ofloat *uvw, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVWCalcClassDef.h"
} ObitUVWCalcClassInfo; 

#endif /* OBITFUVWCALC_H */ 
