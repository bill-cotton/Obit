/* $Id$ */
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
#ifndef OBITOTF_H 
#define OBITOTF_H 

#include "ObitData.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"
#include "ObitOTFArrayGeom.h"
#include "ObitTableList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTF.h
 * ObitOTF GBT On the fly data class definition.
 *
 * This class contains single dish data and allows access.
 * An ObitOTF is the front end to a persistent disk resident structure.
 * There maybe (usually are) associated tables which either describe
 * the data or contain calibration and/or editing information.
 * These associated tables are listed in an #ObitTableList member and
 * the #newObitOTFTable function allows access to these tables.
 * ObitOTF is a derived class from class ObitData.
 * Only FITS data are supported.
 *
 * \section ObitOTFSpecification Specifying desired data transfer parameters
 * The desired data transfers are specified in the member ObitInfoList.
 * There are separate sets of parameters used to specify the FITS or AIPS 
 * data files.
 * Data is read and written as arrays of floats.
 * In the following an ObitInfoList entry is defined by 
 * the name in double quotes, the data type code as an #ObitInfoType enum 
 * and the dimensions of the array (? => depends on application).
 * Only FITS format data are supported.
 *
 * The following apply to all types of files:
 * \li "nVisPIO", OBIT_int, Max. Number of visibilities per 
 *     "Read" or "Write" operation.  Default = 1.
 *
 * \subsection OTF FITS files
 * This implementation uses cfitsio which allows using, in addition to 
 * regular FITS idata, gzip compressed files, pipes, shared memory 
 * and a number of other input forms.
 * The convenience Macro #ObitOTFSetFITS simplifies specifying the 
 * desired data.
 * Binary tables of the type created by AIPS program FITAB are used 
 * for storing visibility data in FITS.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) Name of disk file.
 *
 * \section ObitOTFaccess Creators and Destructors
 * An ObitOTF can be created using newObitOTF which allows specifying 
 * a name for the object.  This name is used to label messages.
 * The copy constructors #ObitOTFClone and ObitOTFCopy make shallow
 * and deep copies of an extant ObitOTF.  If the output ObitOTF has
 * previously been specified, including its disk resident information,
 * then ObitOTFCopy will copy the disk resident as well as the memory 
 * resident information.  Also, any associated tables will be copied.
 *
 * A copy of a pointer to an ObitOTF should always be made using the
 * #ObitOTFRef function which updates the reference count in the object.
 * Then whenever freeing an ObitOTF or changing a pointer, the function
 * #ObitOTFUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitOTFUsage I/O
 * Visibility data is available after an input object is "Opened"
 * and "Read".
 * "Read Select" also allows specifying the data to be read as well as
 * optional calibration and editing to be applied as the data is read.
 * I/O optionally uses a buffer attached to the ObitOTF or some external
 * location.
 * To Write an ObitOTF, create it, open it, and write.
 * The object should be closed to ensure all data is flushed to disk.
 * Deletion of an ObitOTF after its final unreferencing will automatically
 * close it.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitOTF Class structure. */
typedef struct {
#include "ObitOTFDef.h"   /* this class definition */
} ObitOTF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOTF
 * returns a ObitOTF*.
 * in = object to unreference
 */
#define ObitOTFUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOTF.
 * returns a ObitOTF*.
 * in = object to reference
 */
#define ObitOTFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOTFIsA(in) ObitIsA (in, ObitOTFGetClass())

/** 
 * Convenience Macro to define OTF I/O to a FITS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitOTF to specify i/O for.
 *\li nsamp = Max. Number of visibilities per read.
 *\li disk = FITS disk number
 *\li file = Specified FITS file name.
 *\li err = ObitErr to receive error messages.
 */
#define ObitOTFSetFITS(in,nsamp,disk,file,err)  G_STMT_START{       \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_FITS;                            \
       in->info->work[1] = nsamp; in->info->work[2]= disk;          \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       ObitInfoListPut (in->info, "nRecPIO", OBIT_long,              \
		  in->info->dim, (gpointer)&in->info->work[1], err);\
       ObitInfoListPut (in->info, "IOBy", OBIT_long, in->info->dim,  \
		 (gpointer)&in->info->work[1], err);                \
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&in->info->work[2], err); \
       in->info->dim[0] = strlen(file);                             \
       ObitInfoListPut (in->info, "FileName", OBIT_string,          \
                 in->info->dim, (gpointer)file, err);               \
     }G_STMT_END  

/*---------------Public functions---------------------------*/
/**  Public: Class initializer. */
void ObitOTFClassInit (void);

/** Public: Constructor. */
ObitOTF* newObitOTF (gchar* name);

/** Public: Copy Constructor for scratch file. */
ObitOTF* newObitOTFScratch (ObitOTF *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitOTF* (*newObitOTFScratchFP) (ObitOTF *in, ObitErr *err);

/** Public: Fully instantiate. */
void ObitOTFFullInstantiate (ObitOTF *in, gboolean exist, ObitErr *err);
typedef void (*ObitOTFFullInstantiateFP) (ObitOTF *in, gboolean exist, 
					  ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitOTFGetClass (void);

/** Public: Rename underlying structures. */
void ObitOTFRename  (ObitOTF *in, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitIOOTFFITSRename  (ObitIO *in, ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
ObitOTF* ObitOTFZap  (ObitOTF *in, ObitErr *err);

/** Public: Copy (deep) constructor. */
ObitOTF* ObitOTFCopy  (ObitOTF *in, ObitOTF *out, 
		     ObitErr *err);

/** Public: Copy (deep) constructor averaging over frequency. */
ObitOTF* ObitOTFAver (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Copy structure only - no data copy. */
void ObitOTFClone (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Do two OTFs have the same underlying structures?. */
gboolean ObitOTFSame (ObitOTF *in1, ObitOTF *in2, ObitErr *err );
typedef gboolean (*ObitOTFSameFP) (ObitOTF *in1, ObitOTF *in2, 
				   ObitErr *err);

/** Public: Concatenate two OTF, (in at end of out) */
ObitIOCode ObitOTFConcat  (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Create ObitIO structures and open file */
ObitIOCode ObitOTFOpen (ObitOTF *in, ObitIOAccess access, 
			  ObitErr *err);

/** Public: Close file and become inactive */
ObitIOCode ObitOTFClose (ObitOTF *in, ObitErr *err);

/** Public: Reset IO to start of file */
ObitIOCode ObitOTFIOSet (ObitOTF *in, ObitErr *err);
typedef ObitIOCode (*ObitOTFIOSetFP) (ObitOTF *in, ObitErr *err);

/** Public: Read specified data */
ObitIOCode ObitOTFRead (ObitOTF *in, ofloat *data, ObitErr *err);

/** Public: Read select, edit, calibrate specified data */
ObitIOCode ObitOTFReadSelect (ObitOTF *in, ofloat *data, ObitErr *err);

/** Public: Write specified data */
ObitIOCode ObitOTFWrite (ObitOTF *in, ofloat *data, ObitErr *err);

/** Public: Return an associated Table */
ObitTable* newObitOTFTable (ObitOTF *in, ObitIOAccess access, 
			   gchar *tabType, olong *tabver, ObitErr *err);
typedef ObitTable* (*newObitOTFTableFP) (ObitOTF *in, ObitIOAccess access, 
					gchar *tabType, olong *tabver, 
					ObitErr *err);

/** Public: Destroy an associated Table */
ObitIOCode ObitOTFZapTable (ObitOTF *in, gchar *tabType, olong tabVer, 
			    ObitErr *err);
typedef ObitIOCode (*ObitOTFZapTableFP) (ObitOTF *in, gchar *tabType, 
					olong tabVer, ObitErr *err);

/** Public: Copy associated Tables */
ObitIOCode ObitOTFCopyTables (ObitOTF *in, ObitOTF *out, gchar **exclude,
			     gchar **include, ObitErr *err);
typedef ObitIOCode (*ObitOTFCopyTablesFP) (ObitOTF *in, ObitOTF *out, 
					    gchar **exclude, gchar **include, 
					    ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitOTFUpdateTables (ObitOTF *in, ObitErr *err);
typedef ObitIOCode (*ObitOTFUpdateTablesFP) (ObitOTF *in, ObitErr *err);

/** Public: Write header keyword */
void ObitOTFWriteKeyword (ObitOTF *in, 
			  gchar* name, ObitInfoType type, gint32 *dim, 
			  gconstpointer data, ObitErr *err);
/** Public: Read header keyword */
void ObitOTFReadKeyword (ObitOTF *in, 
			 gchar* name, ObitInfoType *type, gint32 *dim, 
			 gpointer data, ObitErr *err);
/** Public: How many records in current scan? */
olong ObitOTFNumRecScan (ObitOTF *inOTF);
/** Public: What is the highest scan number? */
olong ObitOTFHighScan (ObitOTF *inOTF, ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOTFClassDef.h"
} ObitOTFClassInfo; 

#endif /* OBITOTF_H */ 
