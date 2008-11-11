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
#ifndef OBITIOIMAGEAIPS_H 
#define OBITIOIMAGEAIPS_H 
#include "fitsio.h"
#include "ObitIO.h"
#include "ObitFile.h"
#include "ObitImageDesc.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOImageAIPS.h
 * ObitIOImageAIPS class definition.
 *
 * This class is derived from the #ObitIO class.
 *
 * \section ObitIOImageAIPSUsage Usage
 * Instances of this class are for access to AIPS image files.
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * Instances can be made using the #newObitIOImageAIPS constructor,
 * or the #ObitIOImageAIPSCopy copy constructor and pointers copied 
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
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitIOImageAIPS Class structure. */
typedef struct {
#include "ObitIOImageAIPSDef.h" /* This class definition */
} ObitIOImageAIPS;

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOImageAIPS
 * returns a ObitIOImageImageAIPS* (NULL).
 * \li in = object to unreference.
 */
#define ObitIOImageAIPSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an IOImageAIPS.
 * returns a IOImageAIPS*.
 * in = object to reference
 */
#define IOImageAIPSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOImageAIPSIsA(in) ObitIsA (in, ObitIOImageAIPSGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOImageAIPSClassInit (void);

/** Public: Constructor. */
ObitIOImageAIPS* newObitIOImageAIPS (gchar* name, ObitInfoList *info,
				     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOImageAIPSGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOImageAIPSSame (ObitIO *in, ObitInfoList *in1, 
			      ObitInfoList *in2, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitIOImageAIPSRename (ObitIO *in, ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOImageAIPSZap  (ObitIOImageAIPS *in, ObitErr *err);

/** Public: Copy constructor. */
ObitIOImageAIPS* ObitIOImageAIPSCopy  (ObitIOImageAIPS *in, 
				       ObitIOImageAIPS *out, ObitErr *err);

/** Public:  Open */ 
ObitIOCode ObitIOImageAIPSOpen (ObitIOImageAIPS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOImageAIPSClose (ObitIOImageAIPS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOImageAIPSSet (ObitIOImageAIPS *in, ObitInfoList *info, 
			       ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIOImageAIPSRead (ObitIOImageAIPS *in, ofloat *data, 
				ObitErr *err);

/** Public:  Write */
ObitIOCode ObitIOImageAIPSWrite (ObitIOImageAIPS *in, ofloat *data, 
				 ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOImageAIPSFlush (ObitIOImageAIPS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOImageAIPSReadDescriptor (ObitIOImageAIPS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOImageAIPSWriteDescriptor (ObitIOImageAIPS *in, ObitErr *err);

/** Public:  Create buffer */
void ObitIOImageAIPSCreateBuffer (ofloat **data, olong *size, 
				  ObitIOImageAIPS *in, ObitInfoList *info, 
				  ObitErr *err);

/** Public: Create an associated Table 
 * Typed as base class to avoid problems. */
Obit* newObitIOImageAIPSTable (ObitIOImageAIPS *in, ObitIOAccess access, 
		      gchar *tabType, olong *tabver, ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitIOImageAIPSUpdateTables (ObitIOImageAIPS *in, 
					ObitInfoList *info, 
					ObitErr *err);

/** Public: Extract information about underlying file */
void ObitIOImageAIPSGetFileInfo (ObitIO *in, ObitInfoList *myInfo, 
				 gchar *prefix, ObitInfoList *outList, 
				 ObitErr *err);
/*---------------Class Info--------------------------*/

/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOImageAIPSClassDef.h" /* Actual definition */
} ObitIOImageAIPSClassInfo; 

#endif /* OBITIOIMAGEAIPS_H */ 
