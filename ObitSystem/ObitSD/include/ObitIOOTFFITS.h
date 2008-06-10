/* $Id: ObitIOOTFFITS.h,v 1.6 2007/09/11 12:50:20 bcotton Exp $   */
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
#ifndef OBITIOOTFFITS_H 
#define OBITIOOTFFITS_H 
#include "fitsio.h"
#include "ObitIO.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"


/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOOTFFITS.h
 *
 * ObitIOOTFFITS class definition.
 * This class provides an interface to the cfitsio package for FITS images.
 * This class is derived from the ObitIO class.
 *
 * \section ObitIOOTFFITSUsage Usage
 * Instances of this class are for access to FITS format GBT/OTF data.
 * Instances can be made using the $newObitIOOTFFITS constructor,
 * or the #ObitIOOTFFITSCopy copy constructor and pointers copied 
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
/** ObitIOOTFFITS Class structure. */
typedef struct {
  #include "ObitIOOTFFITSDef.h" /* class definition */
} ObitIOOTFFITS;

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOOTFFITS
 * returns a ObitIOOTFFITS* (NULL).
 * in = object to unreference.
 */
#define ObitIOOTFFITSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOOTFFITS.
 * returns a ObitIOOTFFITS*.
 * in = object to reference
 */
#define ObitIOOTFFITSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOOTFFITSIsA(in) ObitIsA (in, ObitIOOTFFITSGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOOTFFITSClassInit (void);

/** Public: Constructor. */
ObitIOOTFFITS* newObitIOOTFFITS (gchar* name, ObitInfoList *info,
			       ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOOTFFITSGetClass (void);

/** Public: Copy  constructor. */
ObitIOOTFFITS* ObitIOOTFFITSCopy  (ObitIOOTFFITS *in, 
				       ObitIOOTFFITS *out, ObitErr *err);

/** Public: Are underlying structures the same. */
gboolean ObitIOOTFFITSSame (ObitIO *in, ObitInfoList *in1, 
			   ObitInfoList *in2, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitIOOTFFITSRename  (ObitIO *in, ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOOTFFITSZap  (ObitIOOTFFITS *in, ObitErr *err);


/** Public:  Open */ 
ObitIOCode ObitIOOTFFITSOpen (ObitIOOTFFITS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOOTFFITSClose (ObitIOOTFFITS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOOTFFITSSet (ObitIOOTFFITS *in, ObitInfoList *info, 
			       ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIOOTFFITSRead (ObitIOOTFFITS *in, ofloat *data, ObitErr *err);

/** Public:  Write */
ObitIOCode ObitIOOTFFITSWrite (ObitIOOTFFITS *in, ofloat *data, ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOOTFFITSFlush (ObitIOOTFFITS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOOTFFITSReadDescriptor (ObitIOOTFFITS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOOTFFITSWriteDescriptor (ObitIOOTFFITS *in, ObitErr *err);

/** Public:  Create buffer */
void
ObitIOOTFFITSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOOTFFITS *in, ObitInfoList *info, 
			     ObitErr *err);

/** Public: Create an associated Table 
 * Typed as base class to avoid problems. */
Obit* 
newObitIOOTFFITSTable (ObitIOOTFFITS *in, ObitIOAccess access, 
		      gchar *tabType, olong *tabVer, ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitIOOTFFITSUpdateTables (ObitIOOTFFITS *in,  ObitInfoList *info,
				     ObitErr *err);

/*---------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOOTFFITSClassDef.h" /* Actual definition */
} ObitIOOTFFITSClassInfo; 

#endif /* OBITIOOTFFITS_H */ 
