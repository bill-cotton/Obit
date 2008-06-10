/* $Id: ObitBeamShape.h,v 1.1 2008/05/06 13:20:14 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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
#ifndef OBITBEAMSHAPE_H 
#define OBITBEAMSHAPE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitBeamShape.h
 *
 * ObitBeamShape Class providing estimates of beam shapes
 *
 * \section ObitBeamShapeaccess Creators and Destructors
 * An ObitBeamShape will usually be created using ObitBeamShapeCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitBeamShape should always be made using the
 * #ObitBeamShapeRef function which updates the reference count in the object.
 * Then whenever freeing an ObitBeamShape or changing a pointer, the function
 * #ObitBeamShapeUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitBeamShape Class structure. */
typedef struct {
#include "ObitBeamShapeDef.h"   /* this class definition */
} ObitBeamShape;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitBeamShape
 * returns a ObitBeamShape*.
 * in = object to unreference
 */
#define ObitBeamShapeUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitBeamShape.
 * returns a ObitBeamShape*.
 * in = object to reference
 */
#define ObitBeamShapeRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitBeamShapeIsA(in) ObitIsA (in, ObitBeamShapeGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitBeamShapeClassInit (void);

/** Public: Default Constructor. */
ObitBeamShape* newObitBeamShape (gchar* name);

/** Public: Create/initialize ObitBeamShape structures */
ObitBeamShape* ObitBeamShapeCreate (gchar* name, ObitImage *image, 
				    ofloat pbmin, ofloat antSize, 
				    gboolean doGain);
/** Typedef for definition of class pointer structure */
typedef ObitBeamShape* (*ObitBeamShapeCreateFP) (gchar* name, ObitImage *image, 
						 ofloat pbmin, ofloat antSize, 
						 gboolean doGain);

/** Public: ClassInfo pointer */
gconstpointer ObitBeamShapeGetClass (void);

/** Public: Copy (deep) constructor. */
ObitBeamShape* ObitBeamShapeCopy  (ObitBeamShape *in, ObitBeamShape *out, 
				   ObitErr *err);

/** Public: Copy structure. */
void ObitBeamShapeClone (ObitBeamShape *in, ObitBeamShape *out, ObitErr *err);

/** Public: Calculate gain in a given direction */
ofloat ObitBeamShapeGain (ObitBeamShape *in, odouble ra, odouble dec, 
			  ofloat parAng);
/** Typedef for definition of class pointer structure */
typedef ofloat (*ObitBeamShapeGainFP) (ObitBeamShape *in, odouble ra, 
				       odouble dec, ofloat parAng);

/** Public: Calculate gain in a symmetric Beam */
ofloat ObitBeamShapeGainSym (ObitBeamShape *in, odouble Angle);
/** Typedef for definition of class pointer structure */
typedef ofloat (*ObitBeamShapeGainSymFP) (ObitBeamShape *in, odouble Angle);

/** Public: Calculate angular offset from beam center */
odouble ObitBeamShapeAngle (ObitBeamShape *in, odouble ra, 
			    odouble dec, ofloat parAng);
/** Typedef for definition of class pointer structure */
typedef odouble (*ObitBeamShapeAngleFP) (ObitBeamShape *in, odouble ra, 
					 odouble dec, ofloat parAng);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitBeamShapeClassDef.h"
} ObitBeamShapeClassInfo; 

#endif /* OBITFBEAMSHAPE_H */ 
