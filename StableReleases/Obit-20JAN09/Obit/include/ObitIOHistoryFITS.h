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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITIOHISTORYFITS_H 
#define OBITIOHISTORYFITS_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIOHistory.h"
#include "ObitTableHistory.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOHistoryFITS.h
 * ObitIOHistoryFITS class definition.
 *
 * This class is derived from the #ObitIOHistory class.
 *
 * This class allows history access to FITS files
 */

/*---------------Class Structure---------------------------*/
/** ObitIOHistoryFITS Class. */
typedef struct {
#include "ObitIOHistoryFITSDef.h"   /* actual definition */
} ObitIOHistoryFITS;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOHistoryFITS
 * returns a ObitIOHistoryFITS*.
 * in = object to unreference
 */
#define ObitIOHistoryFITSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOHistoryFITS.
 * returns a ObitIOHistoryFITS*.
 * in = object to reference
 */
#define ObitIOHistoryFITSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOHistoryFITSIsA(in) ObitIsA (in, ObitIOHistoryFITSGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOHistoryFITSClassInit (void);

/** Public: Constructor. */
ObitIOHistoryFITS* newObitIOHistoryFITS (gchar *name, ObitInfoList *info,
		   ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOHistoryFITSGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOHistoryFITSSame  (ObitIO *in, ObitInfoList *in1, 
				  ObitInfoList *in2, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOHistoryFITSZap  (ObitIOHistoryFITS *in, ObitErr *err);

/** Public: Copy  constructor. */
ObitIOHistoryFITS* ObitIOHistoryFITSCopy  (ObitIOHistoryFITS *in, ObitIOHistoryFITS *out, 
ObitErr *err);

/** Public:  Open */ 
ObitIOCode ObitIOHistoryFITSOpen (ObitIOHistoryFITS *in, ObitIOAccess access, ObitInfoList *info, 
				     ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOHistoryFITSClose (ObitIOHistoryFITS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOHistoryFITSSet (ObitIOHistoryFITS *in, ObitInfoList *info, ObitErr *err);

/** Public:  Read Record */
ObitIOCode ObitIOHistoryFITSReadRrec (ObitIOHistoryFITS *in, olong recno, gchar *hiCard, ObitErr *err);
typedef ObitIOCode (*ObitIOHistoryFITSReadRecFP) (ObitIOHistoryFITS *in, olong recno, gchar *hiCard, 
					      ObitErr *err);

/** Public:  Write Record */
ObitIOCode ObitIOHistoryFITSWriteRec (ObitIOHistoryFITS *in, olong recno, gchar *hiCard, ObitErr *err);
typedef ObitIOCode (*ObitIOHistoryFITSWriteRecFP) (ObitIOHistoryFITS *in, olong recno, gchar *hiCard, 
				     ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOHistoryFITSFlush (ObitIOHistoryFITS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOHistoryFITSReadDescriptor (ObitIOHistoryFITS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOHistoryFITSWriteDescriptor (ObitIOHistoryFITS *in, ObitErr *err);

/** Public:  number of records */
olong ObitIOHistoryFITSNumRec (ObitIOHistory *in);


/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any base class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOHistoryFITSClassDef.h" /* Actual definition */
} ObitIOHistoryFITSClassInfo; 


#endif /* OBITIOHISTORYFITS_H */ 
