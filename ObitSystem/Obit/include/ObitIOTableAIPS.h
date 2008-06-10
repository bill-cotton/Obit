/* $Id$  */
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
#ifndef OBITIOTABLEAIPS_H 
#define OBITIOTABLEAIPS_H 
#include "fitsio.h"
#include "Obit.h"
#include "ObitIO.h"
#include "ObitFile.h"
#include "ObitTableDesc.h"
#include "ObitTableSel.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOTableAIPS.h
 * ObitIOTableAIPS class definition.
 *
 * This class is derived from the #ObitIO class.
 * Related functions are in the 
 * \link ObitIOTableAIPSUtil.h ObitIOTableAIPSUtil 
 * \endlink module.
 *
 * \section ObitIOTableAIPSUsage Usage
 * Instances of this class are for access to AIPS uv data files.
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * Instances can be made using the #newObitIOTableAIPS constructor,
 * or the #ObitIOTableAIPSCopy copy constructor and pointers copied 
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
/** ObitIOTableAIPS Class structure. */
typedef struct {
#include "ObitIOTableAIPSDef.h" /* This class definition */
} ObitIOTableAIPS;

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitIOTableAIPS
 * returns a ObitIOTableTableAIPS* (NULL).
 * \li in = object to unreference.
 */
#define ObitIOTableAIPSUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitIOTableAIPS.
 * returns a ObitIOTableAIPS*.
 * in = object to reference
 */
#define ObitIOTableAIPSRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitIOTableAIPSIsA(in) ObitIsA (in, ObitIOTableAIPSGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitIOTableAIPSClassInit (void);

/** Public: Constructor. */
ObitIOTableAIPS* newObitIOTableAIPS (gchar* name, ObitInfoList *info,
				     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitIOTableAIPSGetClass (void);

/** Public: Are underlying structures the same. */
gboolean ObitIOTableAIPSSame (ObitIO *in, ObitInfoList *in1, 
			      ObitInfoList *in2, ObitErr *err);

/** Public: Delete underlying structures. */
void ObitIOTableAIPSZap  (ObitIOTableAIPS *in, ObitErr *err);

/** Public: Copy constructor. */
ObitIOTableAIPS* ObitIOTableAIPSCopy  (ObitIOTableAIPS *in, 
				       ObitIOTableAIPS *out, ObitErr *err);
/** Public:  Open */ 
ObitIOCode ObitIOTableAIPSOpen (ObitIOTableAIPS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err);

/** Public:  Close */
ObitIOCode ObitIOTableAIPSClose (ObitIOTableAIPS *in, ObitErr *err);

/** Public:  Init I/O */
ObitIOCode ObitIOTableAIPSSet (ObitIOTableAIPS *in, ObitInfoList *info, 
			       ObitErr *err);

/** Public:  Read */
ObitIOCode ObitIOTableAIPSRead (ObitIOTableAIPS *in, ofloat *data, 
				ObitErr *err);

/** Public:  Read specifying start row*/
ObitIOCode ObitIOTableAIPSReadRow (ObitIOTableAIPS *in, olong rowno, ofloat *data, 
				ObitErr *err);

/** Public:  Write */
ObitIOCode ObitIOTableAIPSWrite (ObitIOTableAIPS *in, ofloat *data, 
				 ObitErr *err);

/** Public:  Write specifying start row */
ObitIOCode ObitIOTableAIPSWriteRow (ObitIOTableAIPS *in, olong rowno, ofloat *data, 
				 ObitErr *err);

/** Public:  Flush */
ObitIOCode ObitIOTableAIPSFlush (ObitIOTableAIPS *in, ObitErr *err);

/** Public:  Read Descriptor */
ObitIOCode ObitIOTableAIPSReadDescriptor (ObitIOTableAIPS *in, ObitErr *err);

/** Public:  Write Descriptor */
ObitIOCode ObitIOTableAIPSWriteDescriptor (ObitIOTableAIPS *in, ObitErr *err);

/** Public:  Create buffer */
void ObitIOTableAIPSCreateBuffer (ofloat **data, olong *size, 
				  ObitIOTableAIPS *in, ObitInfoList *info, 
				  ObitErr *err);

/*---------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitIOTableAIPSClassDef.h" /* Actual definition */
} ObitIOTableAIPSClassInfo; 

#endif /* OBITIOTABLEAIPS_H */ 
