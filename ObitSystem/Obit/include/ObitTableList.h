/* $Id: ObitTableList.h,v 1.7 2007/08/31 17:24:48 bcotton Exp $   */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITTABLELIST_H 
#define OBITTABLELIST_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitTable.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableList.h
 * ObitTableList class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is a list of associated tables.
 *
 * \section ObitTableListUsage Usage
 * Instances can be obtained using the #newObitTableList constructor,
 * the #ObitTableListCopy constructor or a pointer duplicated using 
 * the #ObitTableListRef macro.
 * When an instance is no longer needed, use the #ObitTableListUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitTableList Class. */
typedef struct {
#include "ObitTableListDef.h"   /* actual definition */
} ObitTableList;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTableList
 * returns a ObitTableList*.
 * in = object to unreference
 */
#define ObitTableListUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTableList.
 * returns a ObitTableList*.
 * in = object to reference
 */
#define ObitTableListRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTableListIsA(in) ObitIsA (in, ObitTableListGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitTableListClassInit (void);

/** Public: Constructor. */
ObitTableList* newObitTableList (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitTableListGetClass (void);

/** Public: Copy  constructor. */
ObitTableList* 
ObitTableListCopy  (ObitTableList *in, ObitTableList *out, ObitErr *err);

/** Public: Store item to TableList. */
void 
ObitTableListPut(ObitTableList *in, 
		gchar* name, olong *version, ObitTable *table, ObitErr *err);
typedef void (*ObitTableListPutFP) (ObitTableList *in, 
		gchar* name, olong *version, ObitTable *table, ObitErr *err);

/** Public: Retrieve item from TableList by name. */
gboolean 
ObitTableListGet(ObitTableList *in, 
		gchar* name, olong *version, ObitTable **table, ObitErr *err);
typedef gboolean (*ObitTableListGetFP) (ObitTableList *in, 
				       gchar* name, olong version, ObitTable *table, 
				       ObitErr *err);

/** Public: Retrieve item from TableList by number. */
gboolean 
ObitTableListGetNumber(ObitTableList *in, olong number,
		      gchar **name, olong *version, ObitTable **table, ObitErr *err);
typedef gboolean 
(*ObitTableListGetNumberFP) (ObitTableList *in, olong number,
			    gchar **name, olong *version, ObitTable **table, ObitErr *err);

/** Public: Return highered numbered table of a given type. */
olong ObitTableListGetHigh(ObitTableList *in, gchar *name);
typedef olong (*ObitTableListGetHighFP) (ObitTableList *in, gchar* name);

/** Public: Remove item from list. */
void 
ObitTableListRemove (ObitTableList *in, gchar *name,  olong version);
typedef void 
(*ObitTableListRemoveFP) (ObitTableList *in, gchar *name,  olong version);

/** Public: Print list to stderr. */
void ObitTableListPrint  (ObitTableList *in, ObitErr *err);

/** Public: Check validity */
void ObitTableListCheck  (ObitTableList *in, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitTableListClassDef.h" /* Actual definition */
} ObitTableListClassInfo; 


#endif /* OBITTABLELIST_H */ 
