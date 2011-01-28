/* $Id$        */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITIMAGEINTERP_H
#define OBITIMAGEINTERP_H

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitFArray.h"
#include "ObitFInterpolate.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageInterp.h
 *
 * ObitImageInterp Class to interpolate  pixel values in images
 * 
 * \section ObitImageInterpaccess Creators and Destructors
 * An ObitImageInterp will usually be created using ObitImageInterpCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitImageInterp should always be made using the
 * #ObitImageInterpRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImageInterp or changing a pointer, the function
 * #ObitImageInterpUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitImageInterp Class structure. */
typedef struct {
#include "ObitImageInterpDef.h"   /* this class definition */
} ObitImageInterp;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageInterp
 * returns a ObitImageInterp*.
 * in = object to unreference
 */
#define ObitImageInterpUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitImageInterp.
 * returns a ObitImageInterp*.
 * in = object to reference
 */
#define ObitImageInterpRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageInterpIsA(in) ObitIsA (in, ObitImageInterpGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageInterpClassInit (void);

/** Public: Default Constructor. */
ObitImageInterp* newObitImageInterp (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitImageInterpGetClass (void);

/** Public: Copy (deep) constructor. */
ObitImageInterp* ObitImageInterpCopy  (ObitImageInterp *in, 
				       ObitImageInterp *out, 
				       ObitErr *err);

/** Public: Copy structure. */
void ObitImageInterpClone (ObitImageInterp *in, 
			   ObitImageInterp *out, 
			   ObitErr *err);

/** Public: Create/initialize ObitImageInterp structures */
ObitImageInterp* ObitImageInterpCreate (gchar* name, 
					ObitImage *image, olong hwidth, 
					ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ObitImageInterp* (*ObitImageInterpCreateFP) (gchar* name, 
					     ObitImage *image, 
					     olong hwidth, 
					     ObitErr *err);

/** Public: Get pixel value at a position */
ofloat ObitImageInterpValue (ObitImageInterp* in, 
			     odouble RA, odouble Dec, 
			     ofloat Angle, olong plane,
			     ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ofloat (*ObitImageInterpValueFP) (ObitImageInterp* in,
				  odouble RA, odouble Dec, 
				  ofloat Angle, olong plane,
				  ObitErr *err);

/** Public: Get beam value providing interpolator */
ofloat ObitImageInterpValueInt (ObitImageInterp* in, ObitFInterpolate* interp,
				odouble RA, odouble Dec, 
				ofloat PAngle, olong plane,
				ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ofloat (*ObitImageInterpValueIntFP) (ObitImageInterp* in, ObitFInterpolate* interp,
				     odouble RA, odouble Dec, 
				     ofloat PAngle, olong plane,
				     ObitErr *err);

/** Public: Lookup plane for frequency */
olong ObitImageInterpFindPlane (ObitImageInterp* in, odouble freq);
/** Typedef for definition of class pointer structure */
typedef 
olong (*ObitImageInterpFindPlaneFP) (ObitImageInterp* in, odouble freq);

/** Public: Get clone of interpolator */
ObitFInterpolate*  ObitImageInterpCloneInterp (ObitImageInterp* in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ObitFInterpolate* (*ObitImageInterpCloneInterpFP) (ObitImageInterp* in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageInterpClassDef.h"
} ObitImageInterpClassInfo; 

#endif /* OBITIMAGEINTERP_H */ 
