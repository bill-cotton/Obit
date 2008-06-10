/* $Id: ObitFitRegion.h,v 1.2 2007/08/31 17:24:48 bcotton Exp $        */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITFITREGION_H 
#define OBITFITREGION_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFitModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFitRegion.h
 * ObitFitRegion Information about model fitting a regions of the sky
 *
 * This class is derived from the #Obit class.
 * 
 * \section ObitFitRegionaccess Creators and Destructors
 * An ObitFitRegion will usually be created using ObitFitRegionCreate which 
 * allows specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitFitRegion should always be made using the
 * #ObitFitRegionRef function which updates the reference count in the object.
 * Then whenever freeing an ObitFitRegion or changing a pointer, the function
 * #ObitFitRegionUnref will decrement the reference count and destroy the 
 * object when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitFitRegion Class structure. */
typedef struct {
#include "ObitFitRegionDef.h"   /* this class definition */
} ObitFitRegion;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFitRegion
 * returns a ObitFitRegion*.
 * in = object to unreference
 */
#define ObitFitRegionUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFitRegion.
 * returns a ObitFitRegion*.
 * in = object to reference
 */
#define ObitFitRegionRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFitRegionIsA(in) ObitIsA (in, ObitFitRegionGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFitRegionClassInit (void);

/** Public: Default Constructor. */
ObitFitRegion* newObitFitRegion (gchar* name);

/** Public: Create/initialize ObitFitRegion structures */
ObitFitRegion* 
ObitFitRegionCreate (gchar* name, olong corner[2], olong dim[2],
		     ofloat peak, ofloat peakResid, ofloat RMSResid,
		     ofloat fluxResid, olong nmodel, ObitFitModel **models);
/** Typedef for definition of class pointer structure */
typedef ObitFitRegion* 
(*ObitFitRegionCreateFP) (gchar* name, olong corner[2], olong dim[2],
			  ofloat peak, ofloat peakResid, ofloat RMSResid,
			  ofloat fluxResid, olong nmodel, ObitFitModel **models);

/** Public: Generate the region name from an index. */
gchar* ObitFitRegionName(gint indx);

/** Public: ClassInfo pointer */
gconstpointer ObitFitRegionGetClass (void);

/** Public: Copy (deep) constructor. */
ObitFitRegion* ObitFitRegionCopy  (ObitFitRegion *in, ObitFitRegion *out, 
				   ObitErr *err);

/** Public: Copy structure. */
void ObitFitRegionClone (ObitFitRegion *in, ObitFitRegion *out, ObitErr *err);

/** Public: Resize. */
void ObitFitRegionResize  (ObitFitRegion *in, olong nmodel);

/** Public: Subtract from image. */
void ObitFitRegionSubtract (ObitFitRegion* reg, ObitImage *image, ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFitRegionClassDef.h"
} ObitFitRegionClassInfo; 

#endif /* OBITFFITREGION_H */ 
