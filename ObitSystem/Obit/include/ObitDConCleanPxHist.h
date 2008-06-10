/* $Id: ObitDConCleanPxHist.h,v 1.3 2007/08/31 17:24:48 bcotton Exp $ */
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
#ifndef OBITDCONCLEANPXHIST_H 
#define OBITDCONCLEANPXHIST_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanPxHist.h
 * ObitDConCleanPxHist CLEAN image pixel histogram class.
 *
 * This class is derived from the #Obit class.
 * 
 * \section ObitDConCleanPxHistaccess Creators and Destructors
 * An ObitDConCleanPxHist will usually be created using ObitDConCleanPxHistCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanPxHist should always be made using the
 * #ObitDConCleanPxHistRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanPxHist or changing a pointer, the function
 * #ObitDConCleanPxHistUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanPxHist Class structure. */
typedef struct {
#include "ObitDConCleanPxHistDef.h"   /* this class definition */
} ObitDConCleanPxHist;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanPxHist
 * returns a ObitDConCleanPxHist*.
 * in = object to unreference
 */
#define ObitDConCleanPxHistUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanPxHist.
 * returns a ObitDConCleanPxHist*.
 * in = object to reference
 */
#define ObitDConCleanPxHistRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanPxHistIsA(in) ObitIsA (in, ObitDConCleanPxHistGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanPxHistClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanPxHist* newObitDConCleanPxHist (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanPxHistGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanPxHist* 
ObitDConCleanPxHistCopy  (ObitDConCleanPxHist *in, ObitDConCleanPxHist *out, 
			  ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanPxHistClone (ObitDConCleanPxHist *in, 
			       ObitDConCleanPxHist *out, 
			       ObitErr *err);

/** Public: Update with new image and window. */
void ObitDConCleanPxHistUpdate (ObitDConCleanPxHist *in, olong field, 
				olong *plane, ObitImageMosaic *mosaic,
				ObitDConCleanWindow *window, 
				ObitErr *err);

/** Public: Tell how many pixels are larger than a given abs. value. */
olong ObitDConCleanPxHistNumber (ObitDConCleanPxHist *in, ofloat value,
				 ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanPxHistClassDef.h"
} ObitDConCleanPxHistClassInfo; 

#endif /*  OBITDCONCLEANBMHIST_H  */ 
