/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
#ifndef OBITTSYS_H 
#define OBITTSYS_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitTablePC.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPCal.h
 *
 * ObitPCal VLBA Pulsed Cal interpolator class
 *
 * The ObitPCal class provides an interface to a Pulsed Cal (PC)
 * table allowing interpolation for a given time and antenna.
 * 
 * \section ObitPCalaccess Creators and Destructors
 * An ObitPCal will usually be created using ObitPCalCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitPCal should always be made using the
 * #ObitPCalRef function which updates the reference count in the object.
 * Then whenever freeing an ObitPCal or changing a pointer, the function
 * #ObitPCalUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitPCal Class structure. */
typedef struct {
#include "ObitPCalDef.h"   /* this class definition */
} ObitPCal;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitPCal
 * returns a ObitPCal*.
 * in = object to unreference
 */
#define ObitPCalUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitPCal.
 * returns a ObitPCal*.
 * in = object to reference
 */
#define ObitPCalRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitPCalIsA(in) ObitIsA (in, ObitPCalGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitPCalClassInit (void);

/** Public: Default Constructor. */
ObitPCal* newObitPCal (gchar* name);

/** Public: Create/initialize ObitPCal structures */
ObitPCal* ObitPCalCreate (gchar* name, ObitTablePC *PCTable, 
			  ObitUV *UVData, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitPCal* (*ObitPCalCreateFP) (gchar* name, ObitTablePC *PCTable, 
				       ObitUV *UVData, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitPCalGetClass (void);

/** Public: Copy (deep) constructor. */
ObitPCal* ObitPCalCopy  (ObitPCal *in, ObitPCal *out, ObitErr *err);

/** Public: Copy structure. */
void ObitPCalClone (ObitPCal *in, ObitPCal *out, ObitErr *err);

/** Public: Interpolate PCal values. */
void ObitPCalReport (ObitPCal *in, ofloat time, olong ant, olong suba,
		     odouble *CableCal,
		     odouble *Freq1, ofloat *PCal1, 
		     odouble *Freq2, ofloat *PCal2, 
		     ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitPCalClassDef.h"
} ObitPCalClassInfo; 

#endif /* OBITTSYS_H */ 
