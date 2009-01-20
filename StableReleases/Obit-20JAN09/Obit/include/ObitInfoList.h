/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2008                                          */
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
#ifndef OBITINFOLIST_H 
#define OBITINFOLIST_H 
#include <stdio.h>
#include <glib.h>
#include "ObitErr.h"
#include "ObitInfoElem.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitInfoList.h
 * ObitInfoList Linked list of labeled items class definition.
 * This facility allows storing arrays of values of the same (native) 
 * data type and retrieving them by name or order number in the list.
 * Implementation uses the glib GSList class.
 *
 * This class provides a linked list of labeled arrays of data of a given
 * type (#ObitInfoElem).
 * The type codes are defined as enum #obitInfoType.
 * This facility is used for passing information and control parameters.
 *
 * Strings are always special.
 * When a string array is passed, is should be dimensioned similar to 
 * the dim array, i.e. if dim=[10,1,1,1,1] the data should be passed as a 
 * gchar*, if dim=[10,4,1,1,1] the a gchar**.
 * String arrays ar always returned is one block with no NULL, and in 
 * column major order.
 * Strings arrays of up to 3-D (gchar***) can be handled.
 *
 * \section ObitInfoListUsage Usage
 * Instances can be obtained using the #newObitInfoList constructor,
 * the #ObitInfoListCopy constructor or a pointer duplicated using 
 * the #ObitInfoListRef function.
 * When an instance is no longer needed, use the #ObitInfoListUnref 
 * function to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitInfoList Class structure. */
typedef struct {
  /** class name for verification */
  gchar className[16];
  /** Reference count of pointers to this object */
  gint32  ReferenceCount;
  /** Number of entries */
  gint32 number;
  /** glib singly linked list */
  GSList* list;
  /** temporary storage for dim array */
  gint32 dim[MAXINFOELEMDIM];
  /** temporary olong storage */
  olong work[10];
  /** temporary ofloat storage */
  ofloat fwork[10];
} ObitInfoList;

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: constructor. */
ObitInfoList* newObitInfoList (void);

/** Public: destructor. */
ObitInfoList* freeObitInfoList (ObitInfoList *in);

/** Public: Copy constructor. */
ObitInfoList* ObitInfoListCopy (ObitInfoList* in);

/** Public: Reference object pointer. */
ObitInfoList* ObitInfoListRef (ObitInfoList* in);

/** Public: Unreference object pointer. */
ObitInfoList* ObitInfoListUnref (ObitInfoList* in);

/** Public: Copy entries from one list to another. */
ObitInfoList* ObitInfoListCopyData (ObitInfoList* in, ObitInfoList* out);

/** Public: Copy entries from one list to another controlled by a list. */
void ObitInfoListCopyList(ObitInfoList* in, ObitInfoList* out, gchar **list);

/** Public: Copy entries from one list to another controlled by a list with rename. */
void ObitInfoListCopyListRename(ObitInfoList* in, ObitInfoList* out, 
				gchar **inList, gchar **outList);

/** Public: Copy entries from one list to another adding a prefix. */
void ObitInfoListCopyAddPrefix(ObitInfoList* in, ObitInfoList* out, 
			       gchar *prefix);

/** Public: Copy entries with a given prefix from one list to another. */
void ObitInfoListCopyWithPrefix(ObitInfoList* in, ObitInfoList* out, 
				gchar *prefix, gboolean strip);

/** Public: Store item to InfoList. */
void ObitInfoListPut(ObitInfoList *in, 
		      gchar* name, ObitInfoType type, gint32 *dim, 
		      gconstpointer data, ObitErr *err);

/** Public: Store item to InfoList, redefine item if necessary */
void ObitInfoListAlwaysPut(ObitInfoList *in, 
		      gchar* name, ObitInfoType type, gint32 *dim, 
		      gconstpointer data);

/** Public: Get info about item from InfoList by name. */
gboolean ObitInfoListInfo(ObitInfoList *in, 
		      gchar *name, ObitInfoType *type, gint32 *dim, 
		      ObitErr *err);

/** Public: Retrieve item from InfoList by name. */
gboolean ObitInfoListGet(ObitInfoList *in, 
		      gchar *name, ObitInfoType *type, gint32 *dim, 
		      gpointer data, ObitErr *err);

/** Public: Return pointers to an item in InfoList by name. */
gboolean ObitInfoListGetP(ObitInfoList *in, 
			  gchar *name, ObitInfoType *type, gint32 *dim, 
			  gpointer *data);

/** Public: Test retrieve of an item from InfoList. */
gboolean ObitInfoListGetTest(ObitInfoList *in, 
		      gchar *name, ObitInfoType *type, gint32 *dim, 
		      gpointer data);

/** Public: Retrieve item from InfoList by number. */
gboolean 
ObitInfoListGetNumber(ObitInfoList *in, olong number,
		      gchar **name, ObitInfoType *type, gint32 *dim, 
		      gpointer data, ObitErr *err);

/** Public: Return pointers to an item in InfoList by number. */
gboolean 
ObitInfoListGetNumberP(ObitInfoList *in, olong number,
		       gchar **name, ObitInfoType *type, gint32 *dim, 
		       gpointer *data);

/** Public: Remove item from list. */
void ObitInfoListRemove (ObitInfoList *in, gchar *name);

/** Public: Change dimension or type of an item. */
void ObitInfoListResize(ObitInfoList *in, 
			gchar *name, ObitInfoType type, gint32 *dim);

/** Public: Print contents to file (e.g. stdout) */
void ObitInfoListPrint (ObitInfoList *in, FILE *file);

/** Public: Returns true if input is a  ObitInfoList* */
gboolean ObitInfoListIsA (ObitInfoList* in);

#endif /* OBITINFOLIST_H */ 
