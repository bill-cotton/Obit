/* $Id: ObitMem.h,v 1.5 2007/08/31 17:24:48 bcotton Exp $     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
#ifndef OBITMEM_H 
#define OBITMEM_H 
#include "ObitThread.h"
#include "ObitTypes.h"
#include <stdio.h>
#include <stdlib.h>
#include "memwatch.h"  /* For debugging memory problems */

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitMem.h
 * Obit memory manager utilities
 * This system allocates, deallocates and determines the usability of memory.
 * It is designed to minimize the chances of some of the more common c
 * memory problems, especially freeing bad pointers or memory multiple times.
 *
 */

/*------------------- Macroes ----------------------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitMemClassInit (void);

#ifndef MEMWATCH
/** Public: allocate memory. */
gpointer ObitMemAlloc (gulong size);

/** Public: allocate memory and zero fill. */
gpointer ObitMemAlloc0 (gulong size);

/** Public: allocate memory giving name. */
gpointer ObitMemAllocName (gulong size, gchar *name);

/** Public: allocate memory and zero fill, giving name. */
gpointer ObitMemAlloc0Name (gulong size, gchar *name);

/** Public: reallocate memory */
gpointer ObitMemRealloc (gpointer mem, gulong size);

/** Public: deallocate */
gpointer ObitMemFree (gpointer mem);

/** Public: Check if in allocated block */
gboolean ObitMemValid (gpointer mem);

#endif  /* MEMWATCH */

/** Public: Print contents to file (e.g. stdout) */
void ObitMemPrint (FILE *file);

/** Public: Summary of contents */
void ObitMemSummary (olong *number, olong *total);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */

#endif /* OBITMEM_H */ 

