/* $Id: ObitIOImageFITS.h,v 1.9 2008/04/27 20:39:29 bcotton Exp $   */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITIOIMAGEFITS_H 
#define OBITIOIMAGEFITS_H 
#include "fitsio.h"
#include "ObitIO.h"
#include "ObitImageDesc.h"


/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOImageFITS.h
 * ObitIOImageFITS class definition.
 * This class provides an interface to the cfitsio package for FITS images.
 *
 * This class is derived from the #ObitIO class.
 *
 * \section ObitIOImageFITSUsage Usage
 * Instances of this class are for access to FITS image files using
 * the cfitsio package.
 * Instances can be made using the #newObitIOImageFITS constructor,
 * or the #ObitIOImageFITSCopy copy constructor and pointers copied 
 * (with reference pointer update) using #ObitIORef.
 * The destructor (when reference count goes to zero) is
 * #ObitIOUnref.
 * This class should seldom need be accessed directly outside of the 
 * ObitIO class.
 * Parameters needed (passed via ObitInfoList) are:
 * \li "BLC" OBIT_int (?,1,1) the bottom-left corner desired as expressed 
 * in 1-rel pixel indices.  If absent, the value (1,1,1...) will be assumed.
 * dimension of this array is [IM_MAXDIM].
 * \li "TRC" OBIT_int (?,1,1) the top-right corner desired as expressed 
 * in 1-rel pixel indices.  If absent, all pixels are included.
 * dimension of this array is [IM_MAXDIM].
 * \li "IOBy" OBIT_int (1,1,1) an ObitIOSize enum defined in ObitIO.h
 *  giving values OBIT_IO_byRow or  OBIT_IO_byPlane to specify 
 * if the data transfers  are to be by row or plane at a time.
 * \li "FileName" OBIT_string (?,1,1) Name of disk file.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitIOImageFITS Class structure. */
typedef struct {
  #include "ObitIOImageFITSDef.h" /* class definition */
} ObitIOImageFITS;

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOImageFITS
 * returns a ObitIOImageImageFITS* (NULL).
 * in = object to unreference.
 */
#define ObitIOImageFITSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an IOImageFITS.
 * returns a IOImageFITS*.
 * in = object to reference
 */
#define IOImageFITSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOImageFITSIsA(in) ObitIsA (in, ObitIOImageFITSGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOImageFITSClassInit (void);

/** Public: Constructor. */
ObitIOImageFITS* newObitIOImageFITS (gchar* name, ObitInfoList *info,
				     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOImageFITSGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOImageFITSSame (ObitIO *in, ObitInfoList *in1, 
			      ObitInfoList *in2, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitIOImageFITSRename  (ObitIO *in, ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOImageFITSZap  (ObitIOImageFITS *in, ObitErr *err);

/** Public: Copy  constructor. */
ObitIOImageFITS* ObitIOImageFITSCopy  (ObitIOImageFITS *in, 
				       ObitIOImageFITS *out, ObitErr *err);
/** Public: Unconditional destructor. */
ObitIOImageFITS* freeObitIOImageFITS (ObitIOImageFITS *in);

/** Public:  Open */ 
ObitIOCode ObitIOImageFITSOpen (ObitIOImageFITS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOImageFITSClose (ObitIOImageFITS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOImageFITSSet (ObitIOImageFITS *in, ObitInfoList *info, 
			       ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIOImageFITSRead (ObitIOImageFITS *in, ofloat *data, 
				ObitErr *err);

/** Public:  Write */
ObitIOCode ObitIOImageFITSWrite (ObitIOImageFITS *in, ofloat *data, 
				 ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOImageFITSFlush (ObitIOImageFITS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOImageFITSReadDescriptor (ObitIOImageFITS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOImageFITSWriteDescriptor (ObitIOImageFITS *in, ObitErr *err);

/** Public:  Create buffer */
void
ObitIOImageFITSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOImageFITS *in, ObitInfoList *info, 
			     ObitErr *err);

/** Public: Create an associated Table 
 * Typed as base class to avoid problems. */
Obit* 
newObitIOImageFITSTable (ObitIOImageFITS *in, ObitIOAccess access, 
			 gchar *tabType, olong *tabVer, ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitIOImageFITSUpdateTables (ObitIOImageFITS *in, ObitInfoList *info,
					ObitErr *err);

/** Public: Update header BSCALE,BZERO */
void ObitIOImageFITSUpdateScale (ObitIOImageFITS *in, ofloat quant,
				       ObitErr *err);
/*---------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOImageFITSClassDef.h" /* Actual definition */
} ObitIOImageFITSClassInfo; 

#endif /* OBITIOIMAGEFITS_H */ 
