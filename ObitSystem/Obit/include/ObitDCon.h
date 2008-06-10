/* $Id: ObitDCon.h,v 1.3 2007/08/31 17:24:48 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
#ifndef OBITDCON_H 
#define OBITDCON_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImageMosaic.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/*--------------Class definitions-------------------------------------*/

/**
 * \file ObitDCon.h
 * ObitDCon virtual deconvolution base class.
 *
 * This class is derived from the #Obit class.
 *
 * Actual deconvolution classes are derived from this class
 * 
 * \section ObitDConaccess Creators and Destructors
 * An ObitDCon will usually be created using ObitDConCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDCon should always be made using the
 * #ObitDConRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDCon or changing a pointer, the function
 * #ObitDConUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */
/** 
 * ObitDCon virtual deconvolution base class.
 *
 * This class is derived from the #Obit class.
 *
 * Actual deconvolution classes are derived from this class
 * 
 */
typedef struct {
#include "ObitDConDef.h"   /* this class definition */
} ObitDCon;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDCon
 * returns a ObitDCon*.
 * in = object to unreference
 */
#define ObitDConUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDCon.
 * returns a ObitDCon*.
 * in = object to reference
 */
#define ObitDConRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConIsA(in) ObitIsA (in, ObitDConGetClass())

/*---------------Public functions---------------------------*/
/**  Public: Class initializer. */
void ObitDConClassInit (void);

/** Public: Default Constructor. */
ObitDCon* newObitDCon (gchar* name);

/** Public: Create/initialize ObitDCon structures */
ObitDCon* ObitDConCreate (gchar* name, ObitImageMosaic *mosaic, 
			  ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDCon* ObitDConCopy  (ObitDCon *in, ObitDCon *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDConClone (ObitDCon *in, ObitDCon *out, ObitErr *err);

/** Public: Get Parameters. */
void ObitDConGetParms (ObitDCon *in, ObitErr *err);
typedef void (*ObitDConGetParmsFP) (ObitDCon *in, ObitErr *err);

/** Public: Do deconvolution. */
void ObitDConDeconvolve (ObitDCon *in, ObitErr *err);
typedef void (*ObitDConDeconvolveFP) (ObitDCon *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConClassDef.h"
} ObitDConClassInfo; 

#endif /* OBITDCON_H */ 
