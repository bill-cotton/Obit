/* $Id: ObitSourceList.h,v 1.5 2007/08/31 17:24:48 bcotton Exp $      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
#ifndef OBITSOURCELIST_H 
#define OBITSOURCELIST_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitSource.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSourceList.h
 * ObitSourceList class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is a list of sources.
 *
 * \section ObitSourceListUsage Usage
 * Instances can be obtained using the #newObitSourceList constructor,
 * the #ObitSourceListCopy constructor or a pointer duplicated using 
 * the #ObitSourceListRef macro.
 * When an instance is no longer needed, use the #ObitSourceListUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitSourceList Class. */
typedef struct {
#include "ObitSourceListDef.h"   /* actual definition */
} ObitSourceList;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSourceList
 * returns a ObitSourceList*.
 * in = object to unreference
 */
#define ObitSourceListUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSourceList.
 * returns a ObitSourceList*.
 * in = object to reference
 */
#define ObitSourceListRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSourceListIsA(in) ObitIsA (in, ObitSourceListGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSourceListClassInit (void);

/** Public: Constructor. */
ObitSourceList* newObitSourceList (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitSourceListGetClass (void);

/** Public: Copy  constructor. */
ObitSourceList* 
ObitSourceListCopy  (ObitSourceList *in, ObitSourceList *out, ObitErr *err);

/** Public: Create from value */
ObitSourceList* ObitSourceListCreate (gchar* name, olong nsou);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitSourceListClassDef.h" /* Actual definition */
} ObitSourceListClassInfo; 


#endif /* OBITSOURCELIST_H */ 
