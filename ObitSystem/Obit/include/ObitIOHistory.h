/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2010                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITIOHISTORY_H 
#define OBITIOHISTORY_H 
#include "Obit.h"
#include "ObitIO.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOHistory.h
 * ObitIOHistory base class definition.
 * This class is derived from the ObitIO class.
 *
 * This is a virtual base class and should never be directly instantiated,
 * However, its functions should be called and the correct version will 
 * be run.
 * This class is the base for specific History access functions.
 * Derived classes provide an I/O interface to various underlying disk
 * structures.
 * The structure is also defined in ObitIOHistoryDef.h to allow recursive 
 * definition in derived classes. 
 *
 * \section ObitIOHistoryUsage Usage
 * No instances should be created of this class but the class member 
 * functions, given a derived type, will invoke the correct function.
 */

/*---------------Class Structure---------------------------*/
/** ObitIOHistory Class. */
typedef struct {
#include "ObitIOHistoryDef.h"   /* actual definition */
} ObitIOHistory;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOHistory
 * returns a ObitIOHistory*.
 * in = object to unreference
 */
#define ObitIOHistoryUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOHistory.
 * returns a ObitIOHistory*.
 * in = object to reference
 */
#define ObitIOHistoryRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOHistoryIsA(in) ObitIsA (in, ObitIOHistoryGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOHistoryClassInit (void);

/** Public: Constructor. */
ObitIOHistory* newObitIOHistory (gchar *name, ObitInfoList *info,
		   ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOHistoryGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOHistorySame (ObitIOHistory *in, ObitInfoList *in1, 
			    ObitInfoList *in2, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOHistoryZap  (ObitIOHistory *in, ObitErr *err);

/** Public: Copy  constructor. */
ObitIOHistory* ObitIOHistoryCopy  (ObitIOHistory *in, ObitIOHistory *out, ObitErr *err);

/** Public:  Open */ 
ObitIOCode ObitIOHistoryOpen (ObitIOHistory *in, ObitIOAccess access, ObitInfoList *info, 
				     ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOHistoryClose (ObitIOHistory *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOHistorySet (ObitIOHistory *in, ObitInfoList *info, ObitErr *err);

/** Public:  Read Record */
ObitIOCode ObitIOHistoryReadRec (ObitIOHistory *in, olong recno, gchar *hiCard, ObitErr *err);
typedef ObitIOCode (*ObitIOHistoryReadRecFP) (ObitIOHistory *in, olong recno, gchar *hiCard, 
					      ObitErr *err);

/** Public:  Write Record */
ObitIOCode ObitIOHistoryWriteRec (ObitIOHistory *in, olong recno, gchar *hiCard, ObitErr *err);
typedef ObitIOCode (*ObitIOHistoryWriteRecFP) (ObitIOHistory *in, olong recno, gchar *hiCard, 
				     ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOHistoryFlush (ObitIOHistory *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOHistoryReadDescriptor (ObitIOHistory *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOHistoryWriteDescriptor (ObitIOHistory *in, ObitErr *err);

/** Public:  Get number of records */
olong ObitIOHistoryNumRec (ObitIOHistory *in);
typedef olong (*ObitIOHistoryNumRecFP) (ObitIOHistory *in);

/** Public:  Set number of records */
void ObitIOHistorySetNumRec (ObitIOHistory *in, olong current);
typedef void (*ObitIOHistorySetNumRecFP) (ObitIOHistory *in, 
					  olong current);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any base class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOHistoryClassDef.h" /* Actual definition */
} ObitIOHistoryClassInfo; 


#endif /* OBITIOHISTORY_H */ 
