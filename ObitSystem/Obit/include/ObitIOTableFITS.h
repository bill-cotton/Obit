/* $Id$    */
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
#ifndef OBITIOTABLEFITS_H 
#define OBITIOTABLEFITS_H 
#include "fitsio.h"
#include "Obit.h"
#include "ObitIO.h"
#include "ObitTableDesc.h"
#include "ObitTableSel.h"


/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOTableFITS.h
 * ObitIOTableFITS class definition.
 * This class provides an interface to the cfitsio package for FITS images.
 *
 * This class is derived from the #ObitIO class.
 * Related functions are in the 
 * \link ObitIOTableFITSUtil.h ObitIOTableFITSUtil 
 * \endlink module.
 *
 * \section ObitIOTableFITSUsage Usage
 * Instances of this class are for access to FITS image files using
 * the cfitsio package.
 * Instances can be made using the $newObitIOTableFITS constructor,
 * or the #ObitIOTableFITSCopy copy constructor and pointers copied 
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
/** ObitIOTableFITS Class structure. */
typedef struct {
  #include "ObitIOTableFITSDef.h" /* class definition */
} ObitIOTableFITS;

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOTableFITS
 * returns a ObitIOTableImageFITS* (NULL).
 * in = object to unreference.
 */
#define ObitIOTableFITSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOTableFITS.
 * returns a ObitIOTableFITS*.
 * in = object to reference
 */
#define ObitIOTableFITSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOTableFITSIsA(in) ObitIsA (in, ObitIOTableFITSGetClass())
/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOTableFITSClassInit (void);

/** Public: Constructor. */
ObitIOTableFITS* newObitIOTableFITS (gchar* name, ObitInfoList *info,
				     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOTableFITSGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOTableFITSSame (ObitIO *in, ObitInfoList *in1, 
			      ObitInfoList *in2, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOTableFITSZap  (ObitIOTableFITS *in, ObitErr *err);

/** Public: Copy  constructor. */
ObitIOTableFITS* ObitIOTableFITSCopy  (ObitIOTableFITS *in, 
				       ObitIOTableFITS *out, ObitErr *err);

/** Public:  Open */ 
ObitIOCode ObitIOTableFITSOpen (ObitIOTableFITS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOTableFITSClose (ObitIOTableFITS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOTableFITSSet (ObitIOTableFITS *in, ObitInfoList *info, 
			       ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIOTableFITSRead (ObitIOTableFITS *in, ofloat *data, 
				ObitErr *err);

/** Public:  Read specifying start row*/
ObitIOCode ObitIOTableFITSReadRow (ObitIOTableFITS *in, olong rowno, ofloat *data, 
				ObitErr *err);
/** Public:  Write */
ObitIOCode ObitIOTableFITSWrite (ObitIOTableFITS *in, ofloat *data, 
				 ObitErr *err);

/** Public:  Write specifying start row */
ObitIOCode ObitIOTableFITSWriteRow (ObitIOTableFITS *in, olong rowno, ofloat *data, 
				 ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOTableFITSFlush (ObitIOTableFITS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOTableFITSReadDescriptor (ObitIOTableFITS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOTableFITSWriteDescriptor (ObitIOTableFITS *in, ObitErr *err);

/** Public:  Create buffer */
void
ObitIOTableFITSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOTableFITS *in, ObitInfoList *info, 
			     ObitErr *err);

/*---------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOTableFITSClassDef.h" /* Actual definition */
} ObitIOTableFITSClassInfo; 

#endif /* OBITIOTABLEFITS_H */ 
