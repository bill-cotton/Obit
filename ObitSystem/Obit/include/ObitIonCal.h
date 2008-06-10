/* $Id: ObitIonCal.h,v 1.3 2007/08/31 17:24:48 bcotton Exp $        */
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
#ifndef OBITIONCAL_H 
#define OBITIONCAL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitImageMosaic.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIonCal.h
 * ObitIonCal Ionospheric model calibration class
 *
 * This class is derived from the #Obit class.
 *
 * \section ObitIonCalaccess Creators and Destructors
 * An ObitIonCal will usually be created using ObitIonCalCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitIonCal should always be made using the
 * #ObitIonCalRef function which updates the reference count in the object.
 * Then whenever freeing an ObitIonCal or changing a pointer, the function
 * #ObitIonCalUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitIonCal Class structure. */
typedef struct {
#include "ObitIonCalDef.h"   /* this class definition */
} ObitIonCal;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIonCal
 * returns a ObitIonCal*.
 * in = object to unreference
 */
#define ObitIonCalUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIonCal.
 * returns a ObitIonCal*.
 * in = object to reference
 */
#define ObitIonCalRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIonCalIsA(in) ObitIsA (in, ObitIonCalGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIonCalClassInit (void);

/** Public: Default Constructor. */
ObitIonCal* newObitIonCal (gchar* name);

/** Public: Create/initialize ObitIonCal structures */
ObitIonCal* ObitIonCalCreate (gchar* name);
/** Typedef for definition of class pointer structure */
typedef ObitIonCal* (*ObitIonCalCreateFP) (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitIonCalGetClass (void);

/** Public: Copy (deep) constructor. */
ObitIonCal* ObitIonCalCopy  (ObitIonCal *in, ObitIonCal *out, ObitErr *err);

/** Public: Copy structure. */
void ObitIonCalClone (ObitIonCal *in, ObitIonCal *out, ObitErr *err);

/** Public: Attach uv data. */
void ObitIonCalSetData (ObitIonCal *in, ObitUV* inUV);

/** Public: Lookup calibrators for image. */
void ObitIonCalFindImage (ObitIonCal *in, ObitImage* image, ObitErr* err);

/** Public: Fit multiple calibrators in same image. */
void ObitIonCalPosMul (ObitIonCal *in, ObitImage* image, ObitErr* err);
 
/** Public: Fit single epoch Zernike model */
ofloat ObitIonCalFit1 (ObitIonCal *in, olong epoch, ofloat *coef, 
		     ObitErr* err);
 
/** Public: Determine Ionospheric calibration for a UV data */
void ObitIonCaldoCal (ObitIonCal*in, ObitErr* err);

/** Public: Fit position offsets to sources expected at centers of a mosaic */
void ObitIonCalPosMosaic (ObitIonCal *in, ObitImageMosaic* mosaic, 
			  olong epoch, ObitErr* err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIonCalClassDef.h"
} ObitIonCalClassInfo; 

#endif /* OBITFIONCAL_H */ 
