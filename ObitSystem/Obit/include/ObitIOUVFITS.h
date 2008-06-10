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
#ifndef OBITIOUVFITS_H 
#define OBITIOUVFITS_H 
#include "fitsio.h"
#include "Obit.h"
#include "ObitIO.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"


/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOUVFITS.h
 * ObitIOUVFITS class definition.
 * This class provides an interface to the cfitsio package for FITS uv data.
 *

 * This class is derived from the #ObitIO class.
 *
 * \section ObitIOUVFITSUsage Usage
 * Instances of this class are for access to FITS uv data files using
 * the cfitsio package.
 * Instances can be made using the #newObitIOUVFITS constructor,
 * or the #ObitIOUVFITSCopy copy constructor and pointers copied 
 * (with reference pointer update) using #ObitIORef.
 * The destructor (when reference count goes to zero) is
 * #ObitIOUnref.
 * This class should seldom need be accessed directly outside of the 
 * ObitIO class.
 * Parameters needed (passed via ObitInfoList) are:
 * \li "IOBy" OBIT_int (1,1,1) an ObitIOSize enum defined in ObitIO.h
 *  giving values OBIT_IO_byRow or  OBIT_IO_byPlane to specify 
 * if the data transfers  are to be by row or plane at a time.
 * \li "FileName" OBIT_string (?,1,1) Name of disk file.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitIOUVFITS Class structure. */
typedef struct {
  #include "ObitIOUVFITSDef.h" /* class definition */
} ObitIOUVFITS;

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOUVFITS
 * returns a ObitIOUVFITS* (NULL).
 * in = object to unreference.
 */
#define ObitIOUVFITSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOUVFITS.
 * returns a ObitIOUVFITS*.
 * in = object to reference
 */
#define ObitIOUVFITSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOUVFITSIsA(in) ObitIsA (in, ObitIOUVFITSGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOUVFITSClassInit (void);

/** Public: Constructor. */
ObitIOUVFITS* newObitIOUVFITS (gchar* name, ObitInfoList *info,
			       ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOUVFITSGetClass (void);

/** Public: Copy  constructor. */
ObitIOUVFITS* ObitIOUVFITSCopy  (ObitIOUVFITS *in, 
				       ObitIOUVFITS *out, ObitErr *err);

/** Public: Are underlying structures the same. */
gboolean ObitIOUVFITSSame (ObitIO *in, ObitInfoList *in1, 
			   ObitInfoList *in2, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitIOUVFITSRename  (ObitIO *in, ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOUVFITSZap  (ObitIOUVFITS *in, ObitErr *err);


/** Public:  Open */ 
ObitIOCode ObitIOUVFITSOpen (ObitIOUVFITS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOUVFITSClose (ObitIOUVFITS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOUVFITSSet (ObitIOUVFITS *in, ObitInfoList *info, 
			       ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIOUVFITSRead (ObitIOUVFITS *in, ofloat *data, 
				ObitErr *err);

/** Public:  Write */
ObitIOCode ObitIOUVFITSWrite (ObitIOUVFITS *in, ofloat *data, 
				 ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOUVFITSFlush (ObitIOUVFITS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOUVFITSReadDescriptor (ObitIOUVFITS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOUVFITSWriteDescriptor (ObitIOUVFITS *in, ObitErr *err);

/** Public:  Create buffer */
void
ObitIOUVFITSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOUVFITS *in, ObitInfoList *info, 
			     ObitErr *err);

/** Public: Create an associated Table 
 * Typed as base class to avoid problems. */
Obit* 
newObitIOUVFITSTable (ObitIOUVFITS *in, ObitIOAccess access, 
		      gchar *tabType, olong *tabVer, ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitIOUVFITSUpdateTables (ObitIOUVFITS *in,  ObitInfoList *info,
				     ObitErr *err);

/*---------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOUVFITSClassDef.h" /* Actual definition */
} ObitIOUVFITSClassInfo; 

#endif /* OBITIOUVFITS_H */ 
