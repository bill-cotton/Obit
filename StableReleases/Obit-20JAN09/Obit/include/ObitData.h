/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#ifndef OBITDATA_H 
#define OBITDATA_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitTableList.h"
#include "ObitHistory.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitData.h
 * ObitData virtual base class
 *
 * This class is derived from the #Obit class.
 *
 * This class is the virtual base class for Obit data.
 * The derived classes are data access objects which allow access to
 * potentially, multiple data structures and present a uniform
 * internal representation.
 * There maybe (usually are) associated tables which either describe
 * the data or contain calibration and/or editing information.
 * These associated tables are listed in an #ObitTableList member and
 * the #newObitDataTable function allows access to these tables.
 * ObitData is a derived class from class Obit.
 *
 * Generally ObitData objects will be of derived class but a generic
 * ObitData will allow access to tables but not the main file data or 
 * descriptor.
 *
 * \section ObitDataSpecification Specifying desired data parameters
 * The desired data  are specified in the member ObitInfoList.
 * There are separate sets of parameters used to specify the FITS or AIPS 
 * data files.
 * In the following an ObitInfoList entry is defined by 
 * the name in double quotes, the data type code as an #ObitInfoType enum 
 * and the dimensions of the array (? => depends on application).
 * To specify whether the underlying data files are FITS or AIPS
 * \li "FileType" OBIT_int (1,1,1) OBIT_IO_FITS or OBIT_IO_AIPS 
 * which are values of an #ObitIOType enum defined in ObitTypes.h.
 *
 * \subsection FITS files
 * This implementation uses cfitsio which allows using, in addition to 
 * regular FITS idata, gzip compressed files, pipes, shared memory 
 * and a number of other input forms.
 * The convenience Macro #ObitUVSetFITS simplifies specifying the 
 * desired data.
 * Binary tables of the type created by AIPS program FITAB are used 
 * for storing visibility data in FITS.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) Name of disk file.
 *
 * \subsection ObitUVAIPS AIPS files
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * For accessing AIPS files, the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 *
 * \section ObitDataaccess Creators and Destructors
 * There should only be instances of derived classes (
 * which do not include "Data" in the class name).
 * An object derived from ObitData can be created using newObit?
 * which allows specifying  a name for the object.  
 * This name is used to label messages.
 * The copy constructors #Obit?Clone and Obit?Copy make shallow
 * and deep copies of an extant Obit?.  If the output ObitU? has
 * previously been specified, including its disk resident information,
 * then Obit?Copy will copy the disk resident as well as the memory 
 * resident information.  Also, any associated tables will be copied.
 *
 * A copy of a pointer to an ObitData should always be made using the
 * #Obit?Ref function which updates the reference count in the object.
 * Then whenever freeing an ObitUV or changing a pointer, the function
 * #Obit?Unref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitData Class structure. */
typedef struct {
#include "ObitDataDef.h"   /* this class definition */
} ObitData;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitData
 * returns a ObitData*.
 * in = object to unreference
 */
#define ObitDataUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitData.
 * returns a ObitData*.
 * in = object to reference
 */
#define ObitDataRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDataIsA(in) ObitIsA (in, ObitDataGetClass())

/** 
 * Convenience Macro to define ObitData I/O to a FITS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitData to specify i/O for.
 *\li disk = FITS disk number
 *\li file = Specified FITS file name.
 *\li err = ObitErr to receive error messages.
 */
#define ObitDataSetFITS(in,disk,file,err)  G_STMT_START{       \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_FITS;                            \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&in->info->work[2], err); \
       in->info->dim[0] = strlen(file);                             \
       ObitInfoListPut (in->info, "FileName", OBIT_string,          \
                 in->info->dim, (gpointer)file, err);               \
     }G_STMT_END  

/** 
 * Convenience Macro to define ObitData I/O to an AIPS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitData to specify i/O for.
 *\li disk = AIPS disk number
 *\li cno  = catalog slot number
 *\li user = User id number
 *\li err = ObitErr to receive error messages.
 */
#define ObitDataSetAIPS(in,disk,cno,user,err)  G_STMT_START{   \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_AIPS;                            \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&disk, err);              \
       ObitInfoListPut (in->info, "DISK", OBIT_long,                 \
                 in->info->dim, (gpointer)&disk, err);              \
       ObitInfoListPut (in->info, "CNO", OBIT_long,                  \
                 in->info->dim, (gpointer)&cno, err);               \
       ObitInfoListPut (in->info, "User", OBIT_long,                 \
                 in->info->dim, (gpointer)&user, err);              \
     }G_STMT_END   


/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDataClassInit (void);

/** Public: Constructor. */
ObitData* newObitData (gchar* name);

/** Public: Create Data object from description in an ObitInfoList */
ObitData* ObitDataFromFileInfo (gchar *prefix, ObitInfoList *inList, 
				ObitErr *err);
typedef ObitData*
(*ObitDataFromFileInfoFP) (gchar *prefix, ObitInfoList *inList, 
			   ObitErr *err);

/** Public: Copy Constructor for scratch file. */
ObitData* newObitDataScratch (ObitData *in, ObitErr *err);
typedef ObitData* (*newObitDataScratchFP) (ObitData *in, ObitErr *err);

/** Public: Fully instantiate. */
void ObitDataFullInstantiate (ObitData *in, gboolean exist, ObitErr *err);
typedef void (*ObitDataFullInstantiateFP) (ObitData *in, gboolean exist, 
					   ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDataGetClass (void);

/** Public: Rename underlying structures. */
void ObitDataRename  (ObitData *in, ObitErr *err);
typedef void (*ObitDataRenameFP) (ObitData *in, ObitErr *err);

/** Public: Delete underlying structures. */
ObitData* ObitDataZap  (ObitData *in, ObitErr *err);
typedef ObitData* (*ObitDataZapFP) (ObitData *in, ObitErr *err);

/** Public: Copy (deep) constructor. */
ObitData* ObitDataCopy  (ObitData *in, ObitData *out, 
		     ObitErr *err);
typedef ObitData* (*ObitDataCopyFP) (ObitData *in, ObitData *out, 
				     ObitErr *err);

/** Public: Copy structure. */
void ObitDataClone (ObitData *in, ObitData *out, ObitErr *err);
typedef void (*ObitDataCloneFP) (ObitData *in, ObitData *out, 
				 ObitErr *err);

/** Public: Do two ObitDatas have the same underlying structures?. */
gboolean ObitDataSame (ObitData *in1, ObitData *in2, ObitErr *err );
typedef gboolean (*ObitDataSameFP) (ObitData *in1, ObitData *in2, 
				  ObitErr *err);

/** Public: Assign/Initialize IO member  */
void ObitDataSetupIO (ObitData *in, ObitErr *err );
typedef void (*ObitDataSetupIOFP) (ObitData *in, ObitErr *err);

/** Public: Create ObitIO structures and open file */
ObitIOCode ObitDataOpen (ObitData *in, ObitIOAccess access, 
			 ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitIOCode (*ObitDataOpenFP) (ObitData *in, ObitIOAccess access, 
				      ObitErr *err);

/** Public: Close file and become inactive */
ObitIOCode ObitDataClose (ObitData *in, ObitErr *err);
typedef ObitIOCode (*ObitDataCloseFP) (ObitData *in, ObitErr *err);

/** Public: Reset IO to start of file */
ObitIOCode ObitDataIOSet (ObitData *in, ObitErr *err);
typedef ObitIOCode (*ObitDataIOSetFP) (ObitData *in, ObitErr *err);

/** Public: Return an associated Table */
ObitTable* newObitDataTable (ObitData *in, ObitIOAccess access, 
			     gchar *tabType, olong *tabver, ObitErr *err);
typedef ObitTable* (*newObitDataTableFP) (ObitData *in, ObitIOAccess access, 
					  gchar *tabType, olong *tabver, 
					  ObitErr *err);

/** Public: Return an associated History */
ObitHistory* newObitDataHistory (ObitData *in, ObitIOAccess access, ObitErr *err);
typedef ObitHistory* (*newObitDataHistoryFP) (ObitData *in, ObitIOAccess access, 
					      ObitErr *err);

/** Public: Destroy an associated Table */
ObitIOCode ObitDataZapTable (ObitData *in, gchar *tabType, olong tabVer, 
			   ObitErr *err);
typedef ObitIOCode (*ObitDataZapTableFP) (ObitData *in, gchar *tabType, 
					  olong tabVer, ObitErr *err);

/** Public: Copy associated Tables */
ObitIOCode ObitDataCopyTables (ObitData *in, ObitData *out, 
			       gchar **exclude, gchar **include, 
			       ObitErr *err);
typedef ObitIOCode 
(*ObitDataCopyTablesFP) (ObitData *in, ObitData *out, 
			 gchar **exclude, gchar **include, 
			 ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitDataUpdateTables (ObitData *in, ObitErr *err);
typedef ObitIOCode (*ObitDataUpdateTablesFP) (ObitData *in, ObitErr *err);

/** Public: Copy a given table from one ObitData to another */
void ObitDataCopyTable (ObitData *in, ObitData *out, gchar *tabType, 
			olong *inver, olong *outver, ObitErr *err);
typedef void (*ObitDataCopyTableFP) (ObitData *in, ObitData *out, 
				     gchar *tabType, 
				     olong *inver, olong *outver, ObitErr *err);

/** Public: Write header keyword */
void ObitDataWriteKeyword (ObitData *in, 
			   gchar* name, ObitInfoType type, gint32 *dim, 
			   gconstpointer data, ObitErr *err);
typedef void 
(*ObitDataWriteKeywordFP) (ObitData *in, 
			   gchar* name, ObitInfoType type, gint32 *dim, 
			   gconstpointer data, ObitErr *err);

/** Public: Read header keyword */
void ObitDataReadKeyword (ObitData *in, 
			  gchar* name, ObitInfoType *type, gint32 *dim, 
			  gpointer data, ObitErr *err);
typedef void 
(*ObitDataReadKeywordFP) (ObitData *in, 
			  gchar* name, ObitInfoType *type, gint32 *dim, 
			  gpointer data, ObitErr *err);

/** Public: Extract information about underlying file */
void ObitDataGetFileInfo (ObitData *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err);
typedef void 
(*ObitDataGetFileInfoFP) (ObitData *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err);

/*----------- ClassInfo Structure -------------------------------------*/

/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDataClassDef.h"
} ObitDataClassInfo; 

#endif /* OBITDATA_H */ 
