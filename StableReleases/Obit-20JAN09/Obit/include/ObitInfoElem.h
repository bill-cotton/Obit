/* $Id$  */
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
#ifndef OBITINFOELEM_H 
#define OBITINFOELEM_H 
#include <stdio.h>
#include <glib.h>
#include "ObitTypes.h"

/**
 * \file ObitInfoElem.h
 * Elements of an ObitInfoList class definition.
 * Stores a label, size and type info and a data array.
 * The limit on the number of dimensions is MAXINFOELEMDIM.
 * The type coded are defined as enum ObitInfoType.
 */

#undef CLASS
/** ObitInfoElem Class */
#define CLASS "ObitInfoElem"

/** Maximum number of dimensions */
#define MAXINFOELEMDIM 5

typedef struct { 
  /**  element name */
  gchar        *iname;
  /** Data type as enum */
  ObitInfoType itype;
  /** Dimensionality array */
  gint32       idim[MAXINFOELEMDIM];
  /** Size of data in bytes  */
  gint32       size;
  /** Pointer to data array */
  gpointer     data;
}  ObitInfoElem; 
  
/** Constructor  */ 
ObitInfoElem* newObitInfoElem (gchar *label, ObitInfoType type, 
				gint32 *dim, gconstpointer data); 
/** Copy Constructor  */ 
ObitInfoElem* ObitInfoElemCopy (ObitInfoElem *in);

/** Destructor  */ 
void freeObitInfoElem(ObitInfoElem *me); 
  
/**  Compare element name with test string. true=match, else no match  */ 
gboolean ObitInfoElemTest (ObitInfoElem *me, char *testname); 
  
/**  Compare size and type. TRUE = same  */ 
gboolean ObitInfoElemComp (ObitInfoElem *me, ObitInfoType type, gint32 *dim); 
  
/**  Update contents of an info element; returns TRUE if successful.  */ 
gboolean ObitInfoElemUpdate (ObitInfoElem *me, gint32 type, 
			 gint32 *dim, gconstpointer data, gboolean warn); 

/** store data */ 
void ObitInfoElemSave (ObitInfoElem *me, gconstpointer data); 
  
/**  Change type, resize data */ 
void ObitInfoElemResize  (ObitInfoElem *me, ObitInfoType type, gint32 *dim); 

/** determine the size in bytes of an element */ 
olong ObitInfoElemSize  (ObitInfoType type, gint32 *dim); 

/** Print the contents of an object */ 
void  ObitInfoElemPrint(ObitInfoElem *me, FILE *file);

#endif /* OBITINFOELEM_H */ 
