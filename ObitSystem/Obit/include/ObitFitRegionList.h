/* $Id: ObitFitRegionList.h,v 1.2 2007/08/31 17:24:48 bcotton Exp $  */
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
#ifndef OBITFITREGIONLIST_H 
#define OBITFITREGIONLIST_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitFitRegion.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFitRegionList.h
 * ObitFitRegionList List of ObitFitRegion for an image
 *
 * This class is derived from the #Obit class.
 * 
 * \section ObitFitRegionListaccess Creators and Destructors
 * An ObitFitRegionList will usually be created using ObitFitRegionListCreate 
 * which allows specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitFitRegionList should always be made using the
 * #ObitFitRegionListRef function which updates the reference count in the object.
 * Then whenever freeing an ObitFitRegionList or changing a pointer, the function
 * #ObitFitRegionListUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitFitRegionList Class structure. */
typedef struct {
#include "ObitFitRegionListDef.h"   /* this class definition */
} ObitFitRegionList;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFitRegionList
 * returns a ObitFitRegionList*.
 * in = object to unreference
 */
#define ObitFitRegionListUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFitRegionList.
 * returns a ObitFitRegionList*.
 * in = object to reference
 */
#define ObitFitRegionListRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFitRegionListIsA(in) ObitIsA (in, ObitFitRegionListGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFitRegionListClassInit (void);

/** Public: Default Constructor. */
ObitFitRegionList* newObitFitRegionList (gchar* name);

/** Public: Create/initialize ObitFitRegionList structures */
ObitFitRegionList* ObitFitRegionListCreate (gchar* name, ObitImage *image);
/** Typedef for definition of class pointer structure */
typedef ObitFitRegionList* 
(*ObitFitRegionListCreateFP) (gchar* name, ObitImage *image);

/** Public: ClassInfo pointer */
gconstpointer ObitFitRegionListGetClass (void);

/** Public: Append an FitRegion to the list. */
void ObitFitRegionListAppend(ObitFitRegionList *in, ObitFitRegion *reg);

/** Public: Remove an ObitFitRegion from the list. */
void ObitFitRegionListRemove (ObitFitRegionList *in, ObitFitRegion *reg);

/** Public: Find item in a list */
ObitFitRegion*  ObitFitRegionListFind(ObitFitRegionList *in, gchar *name);

/** Public: Subtract regions from image. */
void ObitFitRegionListSubtract (ObitFitRegionList *in, ObitImage *outImage,
				ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFitRegionListClassDef.h"
} ObitFitRegionListClassInfo; 

#endif /* OBITFFITREGIONLIST_H */ 
