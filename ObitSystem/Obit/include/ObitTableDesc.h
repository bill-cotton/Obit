/* $Id$                            */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITTABLEDESC_H 
#define OBITTABLEDESC_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitInfoList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableDesc.h
 * ObitTableDesc Obit table descriptor class definition.
 * This class is derived from the Obit class.
 * This contains information about the observations and the size and 
 * structure of the data.
 *
 * \section ObitTableDescUsage Usage
 * Instances can be obtained using the #newObitTableDesc constructor
 * the #ObitTableDescCopy copy constructor or a pointer duplicated using 
 * the #ObitTableDescRef function.
 * When an instance is no longer needed, use the #ObitTableDescUnref macro
 * to release it.
 */

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTableDesc
 * returns a ObitTableDesc* (NULL).
 * \li in = object to unreference.
 */
#define ObitTableDescUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTableDesc.
 * returns a ObitTableDesc*.
 * in = object to reference
 */
#define ObitTableDescRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTableDescIsA(in) ObitIsA (in, ObitTableDescGetClass())

/** Maximum number of dimensions in regular data array */
#define Table_MAXDIM 7       /* maximum array dimension */
/** Maximum number of "random" parameters */
#define Table_MAX_RANP 14       /* maximum array dimension */
/** Maximum length of descriptor string value */
#define TableLEN_VALUE 41
/** Maximum length of descriptor keyword  */
#define TableLEN_KEYWORD 21

/*--------------Class definitions-------------------------------------*/
/**
 * ObitTableDesc Class structure.
 *
 * This class contains descriptions of interferometric visibility data.
 */  
typedef struct {
#include "ObitTableDescDef.h"  /* Actual definitions */
} ObitTableDesc;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitTableDescClassInit (void);

/** Public: Constructor. */
ObitTableDesc* newObitTableDesc (gchar *name);

/** Public: Copy TableDesc */
ObitTableDesc* ObitTableDescCopy (ObitTableDesc* in, ObitTableDesc* out,
			    ObitErr *err);

/** Public: Return class pointer. */
gconstpointer ObitTableDescGetClass (void);

/** Public: Copy descriptive (nonstructural) information. */
void ObitTableDescCopyDesc (ObitTableDesc* in, ObitTableDesc* out,
			 ObitErr *err);

/** Public: Index for easier access */
void ObitTableDescIndex (ObitTableDesc* in);

/** Public: Reallocate arrays */
void ObitTableDescRealloc (ObitTableDesc* in, olong nfield);

/** Public: Check compatibility */
gboolean ObitTableDescCompatible (ObitTableDesc* in1, 
				  ObitTableDesc* in2);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTableDescClassDef.h" /* Actual definition */
} ObitTableDescClassInfo; 


#endif /* OBITTABLEDESC_H */ 

