/* $Id$   */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITIOHISTORYAIPS_H 
#define OBITIOHISTORYAIPS_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIOHistory.h"
#include "ObitFile.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOHistoryAIPS.h
 * ObitIOHistoryAIPS class definition.
 *
 * This class is derived from #ObitIOHistory
 *
 * This class allows history access to AIPS files
 */

/*---------------Class Structure---------------------------*/
/** ObitIOHistoryAIPS Class. */
typedef struct {
#include "ObitIOHistoryAIPSDef.h"   /* actual definition */
} ObitIOHistoryAIPS;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOHistoryAIPS
 * returns a ObitIOHistoryAIPS*.
 * in = object to unreference
 */
#define ObitIOHistoryAIPSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOHistoryAIPS.
 * returns a ObitIOHistoryAIPS*.
 * in = object to reference
 */
#define ObitIOHistoryAIPSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOHistoryAIPSIsA(in) ObitIsA (in, ObitIOHistoryAIPSGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOHistoryAIPSClassInit (void);

/** Public: Constructor. */
ObitIOHistoryAIPS* newObitIOHistoryAIPS (gchar *name, ObitInfoList *info,
		   ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOHistoryAIPSGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOHistoryAIPSSame  (ObitIO *in, ObitInfoList *in1, 
				 ObitInfoList *in2, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOHistoryAIPSZap  (ObitIOHistoryAIPS *in, ObitErr *err);

/** Public: Copy  constructor. */
ObitIOHistoryAIPS* ObitIOHistoryAIPSCopy  (ObitIOHistoryAIPS *in, ObitIOHistoryAIPS *out, 
					   ObitErr *err);

/** Public:  Open */ 
ObitIOCode ObitIOHistoryAIPSOpen (ObitIOHistoryAIPS *in, ObitIOAccess access, ObitInfoList *info, 
				     ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOHistoryAIPSClose (ObitIOHistoryAIPS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOHistoryAIPSSet (ObitIOHistoryAIPS *in, ObitInfoList *info, ObitErr *err);

/** Public:  Read Record */
ObitIOCode ObitIOHistoryAIPSReadRrec (ObitIOHistoryAIPS *in, olong recno, gchar *hiCard, 
				      ObitErr *err);
typedef ObitIOCode (*ObitIOHistoryAIPSReadRecFP) (ObitIOHistoryAIPS *in, olong recno, 
						  gchar *hiCard, ObitErr *err);

/** Public:  Write Record */
ObitIOCode ObitIOHistoryAIPSWriteRec (ObitIOHistoryAIPS *in, olong recno, gchar *hiCard, 
				      ObitErr *err);
typedef ObitIOCode (*ObitIOHistoryAIPSWriteRecFP) (ObitIOHistoryAIPS *in, olong recno, 
						   gchar *hiCard, ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOHistoryAIPSFlush (ObitIOHistoryAIPS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOHistoryAIPSReadDescriptor (ObitIOHistoryAIPS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOHistoryAIPSWriteDescriptor (ObitIOHistoryAIPS *in, ObitErr *err);

/** Public:  number of records */
olong ObitIOHistoryAIPSNumRec (ObitIOHistory *in);


/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any base class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOHistoryAIPSClassDef.h" /* Actual definition */
} ObitIOHistoryAIPSClassInfo; 


#endif /* OBITIOHISTORYAIPS_H */ 
