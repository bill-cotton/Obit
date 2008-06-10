/* $Id$    */
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
#ifndef OBITFILEFITS_H 
#define OBITFILEFITS_H 
#include "fitsio.h"
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitFile.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFileFITS.h
 * ObitFileFITS class definition.
 *
 * This class is derived from the #ObitFile class.
 *
 * This class provides an I/O interface to FITS files.
 * This implementation uses cfitsio.
 * The structure is also defined in ObitFileFITSDef.h to allow recursive 
 * definition in derived classes. 
 *
 * \section ObitFileFITSUsage Usage
 * Instances of this class are for access to disk files and is used 
 * for access to AIPS data files.
 * Instances can be made using the #newObitFileFITS constructor,
 * or the #ObitFileFITSCopy copy constructor and pointers copied 
 * (with reference pointer update) using #ObitFileFITSRef.
 * The destructor (when reference count goes to zero) is
 * #ObitIOUnref.
 */

/*----------------- typedefs ---------------------------*/
/*---------------Class Structure---------------------------*/
/** ObitFileFITS Class. */
typedef struct {
#include "ObitFileFITSDef.h"   /* actual definition */
} ObitFileFITS;

/*----------------- Macroes ---------------------------*/

/** 
 * Macro to unreference (and possibly destroy) an ObitFileFITS
 * returns a ObitFileFITS*.
 * in = object to unreference
 */
#define ObitFileFITSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFileFITS.
 * returns a ObitFileFITS*.
 * in = object to reference
 */
#define ObitFileFITSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFileFITSIsA(in) ObitIsA (in, ObitFileFITSGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFileFITSClassInit (void);

/** Public: Constructor. */
ObitFileFITS* newObitFileFITS (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitFileFITSGetClass (void);

/** Public:  destroy */
ObitFileFITS* ObitFileFITSZap (ObitFileFITS *in, ObitErr *err);

/** Public:  Destroy current HDU */
ObitIOCode ObitFileFITSZapHDU (ObitFileFITS *in, ObitErr *err);

/** Public: Copy  constructor. */
ObitFileFITS* ObitFileFITSCopy  (ObitFileFITS *in, ObitFileFITS *out, ObitErr *err);

/** Public:  Open */ 
ObitIOCode 
ObitFileFITSOpen (ObitFileFITS *in, gchar *fileName, olong disk, 
		  ObitIOAccess access,  ObitErr *err);

/** Public:  Close */
ObitIOCode ObitFileFITSClose (ObitFileFITS *in, ObitErr *err);

/** Public: Does a given file exist? */
gboolean ObitFileFITSExist (gchar *fileName, ObitErr *err);

/** Public:  Position to given extention number */
ObitIOCode ObitFileFITSNum (ObitFileFITS *in, olong hdunum, olong *hdutype, 
			    ObitErr *err);

/** Public:  Position to given extention name */
ObitIOCode ObitFileFITSName (ObitFileFITS *in, olong hdutype, gchar *extname, 
			    olong extver, ObitErr *err);

/** Public: Read String keyword */
ObitIOCode 
ObitFileFITSReadKeyStr (ObitFileFITS *in, gchar *Name, 
			gchar *Value, gchar *Comment, ObitErr *err);

/** Public: Read float keyword */
ObitIOCode 
ObitFileFITSReadKeyFlt (ObitFileFITS *in, gchar *Name, 
			ofloat *Value, gchar *Comment, ObitErr *err);

/** Public: Read double keyword */
ObitIOCode 
ObitFileFITSReadKeyDbl (ObitFileFITS *in, gchar *Name, 
			odouble *Value, gchar *Comment, ObitErr *err);

/** Public: Read long keyword */
ObitIOCode 
ObitFileFITSReadKeyLng (ObitFileFITS *in, gchar *Name, 
			olong *Value, gchar *Comment, ObitErr *err);

/** Public: Read next HISTORY header card */
ObitIOCode 
ObitFileFITSReadHistory (ObitFileFITS *in, gchar *hiCard, ObitErr *err);

/** Public: Write String keyword */
ObitIOCode 
ObitFileFITSWriteKeyStr (ObitFileFITS *in, gchar *Name, gboolean update,
			gchar *Value, gchar *Comment, ObitErr *err);

/** Public: Write float keyword */
ObitIOCode 
ObitFileFITSWriteKeyFlt (ObitFileFITS *in, gchar *Name, gboolean update,
			ofloat Value, gchar *Comment, ObitErr *err);

/** Public: Write double keyword */
ObitIOCode 
ObitFileFITSWriteKeyDbl (ObitFileFITS *in, gchar *Name, gboolean update,
			odouble Value, gchar *Comment, ObitErr *err);

/** Public: Write long keyword */
ObitIOCode 
ObitFileFITSWriteKeyLng (ObitFileFITS *in, gchar *Name, gboolean update,
			olong Value, gchar *Comment, ObitErr *err);

/** Public: Write String HISTORY keyword */
ObitIOCode 
ObitFileFITSWriteHisKeyStr (fitsfile *inFptr, gchar *Name, gchar *Value, 
			    gchar *Comment, ObitErr *err);

/** Public: Write float HISTORY keyword */
ObitIOCode 
ObitFileFITSWriteHisKeyFlt (fitsfile *inFptr, gchar *Name, ofloat Value, 
			    gchar *Comment, ObitErr *err);

/** Public: Write double HISTORY keyword */
ObitIOCode 
ObitFileFITSWriteHisKeyDbl (fitsfile *inFptr, gchar *Name, odouble Value, 
			    gchar *Comment, ObitErr *err);

/** Public: Write long HISTORY keyword */
ObitIOCode 
ObitFileFITSWriteHisKeyLng (fitsfile *inFptr, gchar *Name, olong Value, 
			    gchar *Comment, ObitErr *err);

/** Public: Write Date keyword */
ObitIOCode 
ObitFileFITSWriteDate (ObitFileFITS *in, ObitErr *err);

/** Public: Write next HISTORY header card */
ObitIOCode 
ObitFileFITSWriteHistory (ObitFileFITS *in, gchar *hiCard, ObitErr *err);

/** Public: Expand the current HDR by a specified number of keywords*/
ObitIOCode 
ObitFileFITSAddKeys (ObitFileFITS *in, olong morekeys, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFileFITSClassDef.h" /* Actual definition */
} ObitFileFITSClassInfo; 


#endif /* OBITFILEFITS_H */ 
