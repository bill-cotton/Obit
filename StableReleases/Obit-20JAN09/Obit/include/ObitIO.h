/* $Id$              */
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
#ifndef OBITIO_H 
#define OBITIO_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitImageDesc.h"
#include "ObitUVDesc.h"
#include "ObitUVCal.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIO.h
 * ObitIO base class definition.

 * This class is derived from the #Obit class.
 *
 * This is a virtual base class and should never be directly instantiated,
 * However, its functions mshould be called and the correct version will 
 * be run.
 * Derived classes provide an I/O interface to various underlying disk
 * structures.
 * The structure is also defined in ObitIODef.h to allow recursive 
 * definition in derived classes. 
 *
 * \section ObitIOUsage Usage
 * No instances should be created of this class but the class member 
 * functions, given a derived type, will invoke the correct function.
 */

/*---------------Class Structure---------------------------*/
/** ObitIO Class. */
typedef struct {
#include "ObitIODef.h"   /* actual definition */
} ObitIO;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIO
 * returns a ObitIO*.
 * in = object to unreference
 */
#define ObitIOUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIO.
 * returns a ObitIO*.
 * in = object to reference
 */
#define ObitIORef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOIsA(in) ObitIsA (in, ObitIOGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOClassInit (void);

/** Public: Constructor. */
ObitIO* newObitIO (gchar *name, ObitInfoList *info,
		   ObitErr *err);
/** define type for ClassInfo structure */
typedef ObitIO* (*newObitIOFP) (gchar* name, ObitInfoList *info,
				ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOSame (ObitIO *in, ObitInfoList *in1, ObitInfoList *in2, 
		     ObitErr *err);
typedef gboolean (*ObitIOSameFP) (ObitIO *in, ObitInfoList *in1, 
				  ObitInfoList *in2, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitIORename  (ObitIO *in, ObitInfoList *info, ObitErr *err);
typedef void (*ObitIORenameFP) (ObitIO *in,  ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOZap  (ObitIO *in, ObitErr *err);
typedef void (*ObitIOZapFP) (ObitIO *in, ObitErr *err);

/** Public: Copy  constructor. */
ObitIO* ObitIOCopy  (ObitIO *in, ObitIO *out, ObitErr *err);

/** Public:  Open */ 
ObitIOCode ObitIOOpen (ObitIO *in, ObitIOAccess access, ObitInfoList *info, 
	     ObitErr *err);
typedef ObitIOCode (*ObitIOOpenFP) (ObitIO *in, ObitIOAccess access, 
				  ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOClose (ObitIO *in, ObitErr *err);
typedef ObitIOCode (*ObitIOCloseFP) (ObitIO *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOSet (ObitIO *in, ObitInfoList *info, ObitErr *err);
typedef ObitIOCode (*ObitIOSetFP) (ObitIO *in, ObitInfoList *info, 
				 ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIORead (ObitIO *in, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitIOReadFP) (ObitIO *in, ofloat *data, ObitErr *err);

/** Public:  Read Row */
ObitIOCode ObitIOReadRow (ObitIO *in, olong rowno, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitIOReadRowFP) (ObitIO *in, olong rowno, ofloat *data, 
				    ObitErr *err);

/** Public:  Read with selection */
ObitIOCode ObitIOReadSelect (ObitIO *in, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitIOReadSelectFP) (ObitIO *in, ofloat *data, 
					  ObitErr *err);

/** Public:  Read Row with selection */
ObitIOCode ObitIOReadRowSelect (ObitIO *in, olong rowno, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitIOReadRowSelectFP) (ObitIO *in, olong rowno, ofloat *data, 
					  ObitErr *err);

/** Public:  Write */
ObitIOCode ObitIOWrite (ObitIO *in, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitIOWriteFP) (ObitIO *in, ofloat *data, ObitErr *err);

/** Public:  Write Row */
ObitIOCode ObitIOWriteRow (ObitIO *in, olong rowno, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitIOWriteRowFP) (ObitIO *in, olong rowno, ofloat *data, 
				     ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOFlush (ObitIO *in, ObitErr *err);
typedef ObitIOCode (*ObitIOFlushFP) (ObitIO *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOReadDescriptor (ObitIO *in, ObitErr *err);
typedef ObitIOCode (*ObitIOReadDescriptorFP) (ObitIO *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOWriteDescriptor (ObitIO *in, ObitErr *err);
typedef ObitIOCode (*ObitIOWriteDescriptorFP) (ObitIO *in, ObitErr *err);

/** Public:  Create buffer */
void ObitIOCreateBuffer (ofloat **data, olong *size, ObitIO *in, 
			 ObitInfoList *info, ObitErr *err);
typedef void (*ObitIOCreateBufferFP) (ofloat **data, olong *size, ObitIO *in, 
			 ObitInfoList *info, ObitErr *err);

/** Public: Destroy buffer */
void ObitIOFreeBuffer (ofloat *buffer);
typedef void (*ObitIOFreeBufferFP) (ofloat *buffer);

/** Public: Create an associated Table 
 * Typed as base class to avoid problems. */
Obit* newObitIOTable (ObitIO *in, ObitIOAccess access, 
		      gchar *tabType, olong *tabver, ObitErr *err);
typedef Obit* (*newObitIOTableFP) (ObitIO *in, ObitIOAccess access, 
				   gchar *tabType, olong *tabver, 
				   ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitIOUpdateTables (ObitIO *in, ObitInfoList *info, ObitErr *err);
typedef ObitIOCode (*ObitIOUpdateTablesFP) (ObitIO *in, ObitInfoList *info, 
					    ObitErr *err);

/** Public: Extract information about underlying file */
void ObitIOGetFileInfo (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
			ObitInfoList *outList, ObitErr *err);
typedef void 
(*ObitIOGetFileInfoFP) (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
			ObitInfoList *outList, ObitErr *err);

/*-------------------Class Info--------------------------*/

/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any base class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOClassDef.h" /* Actual definition */
} ObitIOClassInfo; 


#endif /* OBITIO_H */ 
