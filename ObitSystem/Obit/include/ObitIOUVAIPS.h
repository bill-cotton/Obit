/* $Id: ObitIOUVAIPS.h,v 1.7 2007/08/31 17:24:48 bcotton Exp $     */
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
#ifndef OBITIOUVAIPS_H 
#define OBITIOUVAIPS_H 
#include "fitsio.h"
#include "Obit.h"
#include "ObitIO.h"
#include "ObitFile.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOUVAIPS.h
 * ObitIOUVAIPS class definition.
 *
 * This class is derived from the #ObitIO class.
 *
 * \section ObitIOUVAIPSUsage Usage
 * Instances of this class are for access to AIPS uv data files.
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * Instances can be made using the #newObitIOUVAIPS constructor,
 * or the #ObitIOUVAIPSCopy copy constructor and pointers copied 
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
/** ObitIOUVAIPS Class structure. */
typedef struct {
#include "ObitIOUVAIPSDef.h" /* This class definition */
} ObitIOUVAIPS;

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOUVAIPS
 * returns a ObitIOUVUVAIPS* (NULL).
 * \li in = object to unreference.
 */
#define ObitIOUVAIPSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOUVAIPS.
 * returns a ObitIOUVAIPS*.
 * in = object to reference
 */
#define ObitIOUVAIPSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOUVAIPSIsA(in) ObitIsA (in, ObitIOUVAIPSGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOUVAIPSClassInit (void);

/** Public: Constructor. */
ObitIOUVAIPS* newObitIOUVAIPS (gchar* name, ObitInfoList *info,
			       ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOUVAIPSGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOUVAIPSSame (ObitIO *in, ObitInfoList *in1, 
			   ObitInfoList *in2, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitIOUVAIPSRename (ObitIO *in, ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOUVAIPSZap  (ObitIOUVAIPS *in, ObitErr *err);

/** Public: Copy constructor. */
ObitIOUVAIPS* ObitIOUVAIPSCopy  (ObitIOUVAIPS *in, 
				       ObitIOUVAIPS *out, ObitErr *err);
/** Public:  Open */ 
ObitIOCode ObitIOUVAIPSOpen (ObitIOUVAIPS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOUVAIPSClose (ObitIOUVAIPS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOUVAIPSSet (ObitIOUVAIPS *in, ObitInfoList *info, 
			       ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIOUVAIPSRead (ObitIOUVAIPS *in, ofloat *data, 
				ObitErr *err);

/** Public:  Write */
ObitIOCode ObitIOUVAIPSWrite (ObitIOUVAIPS *in, ofloat *data, 
				 ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOUVAIPSFlush (ObitIOUVAIPS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOUVAIPSReadDescriptor (ObitIOUVAIPS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOUVAIPSWriteDescriptor (ObitIOUVAIPS *in, ObitErr *err);

/** Public:  Create buffer */
void ObitIOUVAIPSCreateBuffer (ofloat **data, olong *size, 
				  ObitIOUVAIPS *in, ObitInfoList *info, 
				  ObitErr *err);

/** Public: Create an associated Table 
 * Typed as base class to avoid problems. */
Obit* newObitIOUVAIPSTable (ObitIOUVAIPS *in, ObitIOAccess access, 
		      gchar *tabType, olong *tabver, ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitIOUVAIPSUpdateTables (ObitIOUVAIPS *in, ObitInfoList *info, 
				     ObitErr *err);

/*---------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOUVAIPSClassDef.h" /* Actual definition */
} ObitIOUVAIPSClassInfo; 

#endif /* OBITIOUVAIPS_H */ 
