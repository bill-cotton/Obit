/* $Id$    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
#ifndef OBITIMAGESEL_H 
#define OBITIMAGESEL_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitInfoList.h"
#include "ObitImageDesc.h"
#include "ObitFArray.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageSel.h
 * ObitImageSel Obit image selection class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This contains information about portions of an image selected.
 *
 * \section ObitImageSelUsage Usage
 * Instances can be obtained using the #newObitImageSel constructor
 * the #ObitImageSelCopy copy constructor or a pointer duplicated using 
 * the #ObitImageSelRef function.
 * When an instance is no longer needed, use the #ObitImageSelUnref macro
 * to release it.
 */

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageSel
 * returns a ObitImageSel* (NULL).
 * \li in = object to unreference.
 */
#define ObitImageSelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitImageSel.
 * returns a ObitImageSel*.
 * in = object to reference
 */
#define ObitImageSelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageSelIsA(in) ObitIsA (in, ObitImageSelGetClass())

/** Maximum number of dimensions */
#define IM_MAXDIM 7       /* maximum array dimension */

/*--------------Class definitions-------------------------------------*/
/**
 * ObitImageSel Class structure.
 *
 * This class contains descriptions of interferometric visibility data.
 */  
typedef struct {
#include "ObitImageSelDef.h" /* actual definition */
} ObitImageSel;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageSelClassInit (void);

/** Public: Constructor. */
ObitImageSel* newObitImageSel (gchar *name);

/** Public: Return class pointer. */
gconstpointer ObitImageSelGetClass (void);

/** Public: Copy ImageSel */
ObitImageSel* ObitImageSelCopy (ObitImageSel* in, ObitImageSel* out,
				ObitErr *err);

/** Public: Create/resize buffer for image? */
ObitFArray* 
ObitImageSelBuffer (ObitFArray *buffer, ObitImageDesc* desc, 
		    ObitImageSel* sel);

/** Public: Enforces defaults in inaxes, blc, trc */
void ObitImageSelDefault (ObitImageDesc* in, 
			  ObitImageSel* sel);

/** Public: Applies selection to a Descriptor */
void ObitImageSelSetDesc (ObitImageDesc* in, ObitImageSel* sel,
			  ObitImageDesc* out, ObitErr *err);
/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageSelClassDef.h" /* Actual definition */
} ObitImageSelClassInfo; 

#endif /* OBITIMAGESEL_H */ 

