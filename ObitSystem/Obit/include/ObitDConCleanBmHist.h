/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                               */
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
#ifndef OBITDCONCLEANBMHIST_H 
#define OBITDCONCLEANBMHIST_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanBmHist.h
 * ObitDConCleanBmHist CLEAN beam histogram class.
 *
 * This class is derived from the #Obit class.
 *
 * This class has information about the maximum sidelobe exterior to 
 * a given radius.
 * 
 * \section ObitDConCleanBmHistaccess Creators and Destructors
 * An ObitDConCleanBmHist will usually be created using ObitDConCleanBmHistCreate 
 * which allows specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanBmHist should always be made using the
 * #ObitDConCleanBmHistRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanBmHist or changing a pointer, the function
 * #ObitDConCleanBmHistUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanBmHist Class structure. */
typedef struct {
#include "ObitDConCleanBmHistDef.h"   /* this class definition */
} ObitDConCleanBmHist;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanBmHist
 * returns a ObitDConCleanBmHist*.
 * in = object to unreference
 */
#define ObitDConCleanBmHistUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanBmHist.
 * returns a ObitDConCleanBmHist*.
 * in = object to reference
 */
#define ObitDConCleanBmHistRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanBmHistIsA(in) ObitIsA (in, ObitDConCleanBmHistGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanBmHistClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanBmHist* newObitDConCleanBmHist (gchar* name);

/** Public: Create/initialize ObitDConCleanBmHist structures */
ObitDConCleanBmHist* 
ObitDConCleanBmHistCreate (gchar* name, ObitImage *Beam,  ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanBmHistGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanBmHist* 
ObitDConCleanBmHistCopy  (ObitDConCleanBmHist *in, ObitDConCleanBmHist *out, 
			  ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanBmHistClone (ObitDConCleanBmHist *in, 
			       ObitDConCleanBmHist *out, 
			       ObitErr *err);

/** Public: Update with new Beam image. */
void ObitDConCleanBmHistUpdate (ObitDConCleanBmHist *in, ObitImage *Beam,
				olong *plane, ObitErr *err);

/** Public: Give maximum abs. exterior sidelobe */
ofloat ObitDConCleanBmHistPeak (ObitDConCleanBmHist *in, olong radius,
				ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanBmHistClassDef.h"
} ObitDConCleanBmHistClassInfo; 

#endif /*  OBITDCONCLEANBMHIST_H  */ 
