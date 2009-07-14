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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVFULLBEAM_H 
#define OBITUVFULLBEAM_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitUV.h"
#include "ObitFArray.h"
#include "ObitFInterpolate.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFullBeam.h
 *
 * ObitFullBeam Class to generate full beam corrections from beam images
 * An ObitFullBeam takes images or (hyper)cubes  of a primary beam in a given
 * Stokes parameter and assists in image plane corrections
 * 
 * \section ObitFullBeamaccess Creators and Destructors
 * An ObitFullBeam will usually be created using ObitFullBeamCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitFullBeam should always be made using the
 * #ObitFullBeamRef function which updates the reference count in the object.
 * Then whenever freeing an ObitFullBeam or changing a pointer, the function
 * #ObitFullBeamUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitFullBeam Class structure. */
typedef struct {
#include "ObitFullBeamDef.h"   /* this class definition */
} ObitFullBeam;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFullBeam
 * returns a ObitFullBeam*.
 * in = object to unreference
 */
#define ObitFullBeamUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFullBeam.
 * returns a ObitFullBeam*.
 * in = object to reference
 */
#define ObitFullBeamRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFullBeamIsA(in) ObitIsA (in, ObitFullBeamGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFullBeamClassInit (void);

/** Public: Default Constructor. */
ObitFullBeam* newObitFullBeam (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitFullBeamGetClass (void);

/** Public: Copy (deep) constructor. */
ObitFullBeam* ObitFullBeamCopy  (ObitFullBeam *in, 
				     ObitFullBeam *out, 
				     ObitErr *err);

/** Public: Copy structure. */
void ObitFullBeamClone (ObitFullBeam *in, 
			  ObitFullBeam *out, 
			  ObitErr *err);

/** Public: Create/initialize ObitFullBeam structures */
ObitFullBeam* ObitFullBeamCreate (gchar* name, ObitInfoList *myInput,
				  ObitImage *image, 
				  ObitUV *uvData, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ObitFullBeam* (*ObitFullBeamCreateFP) (gchar* name, 
				       ObitInfoList *myInput,
				       ObitImage *image, 
				       ObitUV *uvData, 
				       ObitErr *err);

/** Public: Get beam value */
ofloat ObitFullBeamValue (ObitFullBeam* in, 
			  odouble RA, odouble Dec, 
			  ofloat PAngle, olong plane,
			  ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ofloat (*ObitFullBeamValueFP) (ObitFullBeam* in,
			     odouble RA, odouble Dec, 
			     ofloat PAngle, olong plane,
			     ObitErr *err);

/** Public: Get beam value providing interpolator */
ofloat ObitFullBeamValueInt (ObitFullBeam* in, ObitFInterpolate* interp,
			  odouble RA, odouble Dec, 
			  ofloat PAngle, olong plane,
			  ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ofloat (*ObitFullBeamValueIntFP) (ObitFullBeam* in, ObitFInterpolate* interp,
			     odouble RA, odouble Dec, 
			     ofloat PAngle, olong plane,
			     ObitErr *err);

/** Public: Lookup plane for frequency */
olong ObitFullBeamFindPlane (ObitFullBeam* in, odouble freq);
/** Typedef for definition of class pointer structure */
typedef 
olong (*ObitFullBeamFindPlaneFP) (ObitFullBeam* in, odouble freq);

/** Public: Get clone of interpolator */
ObitFInterpolate*  ObitFullBeamCloneInterp (ObitFullBeam* in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef 
ObitFInterpolate* (*ObitFullBeamCloneInterpFP) (ObitFullBeam* in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFullBeamClassDef.h"
} ObitFullBeamClassInfo; 

#endif /* OBITFUVFULLBEAM_H */ 
