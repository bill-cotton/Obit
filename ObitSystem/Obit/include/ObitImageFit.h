/* $Id$        */
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
#ifndef OBITIMAGEFIT_H 
#define OBITIMAGEFIT_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFitRegion.h"
#include "ObitImageFitData.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageFit.h
 * ObitImageFit fits models to images
 *
 * This class is derived from the #Obit class.
 *
 * Fits image models in an ObitFitRegion.
 * 
 * \section ObitImageFitaccess Creators and Destructors
 * An ObitImageFit will usually be created using ObitImageFitCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitImageFit should always be made using the
 * #ObitImageFitRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImageFit or changing a pointer, the function
 * #ObitImageFitUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitImageFit Class structure. */
typedef struct {
#include "ObitImageFitDef.h"   /* this class definition */
} ObitImageFit;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageFit
 * returns a ObitImageFit*.(NULL)
 * in = object to unreference
 */
#define ObitImageFitUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitImageFit.
 * returns a ObitImageFit*.
 * in = object to reference
 */
#define ObitImageFitRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageFitIsA(in) ObitIsA (in, ObitImageFitGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageFitClassInit (void);

/** Public: Default Constructor. */
ObitImageFit* newObitImageFit (gchar* name);

/** Public: Create/initialize ObitImageFit structures */
ObitImageFit* 
ObitImageFitCreate (gchar* name);
/** Typedef for definition of class pointer structure */
typedef ObitImageFit* 
(*ObitImageFitCreateFP) (gchar* name);

/** Public: Fit an ObitFitRegion */
olong ObitImageFitFit 
(ObitImageFit* in,  ObitImage *image, ObitFitRegion* reg, 
 ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef olong (*ObitImageFitFitFP) 
     (ObitImageFit* in,  ObitImage *image, ObitFitRegion* reg, 
      ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitImageFitGetClass (void);

/** Public: Copy (deep) constructor. */
ObitImageFit* ObitImageFitCopy  (ObitImageFit *in, ObitImageFit *out, 
				 ObitErr *err);

/** Public: Copy structure. */
void ObitImageFitClone (ObitImageFit *in, ObitImageFit *out, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageFitClassDef.h"
} ObitImageFitClassInfo; 

#endif /* OBITIMAGEFIT_H */ 
