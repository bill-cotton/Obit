/* $Id: ObitUV.h,v 1.15 2007/08/31 17:24:48 bcotton Exp $       */
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
#ifndef OBITUV_H 
#define OBITUV_H 

#include "ObitData.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitTableList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUV.h
 * ObitUV uv data class definition.
 *
 * This class is derived from the #ObitData class.
 * Related functions are in the 
 * \link ObitUVUtil.h ObitUVUtil 
 * \endlink ,
 * \link ObitUVEdit.h ObitUVEdit 
 * \endlink and
 * \link ObitUVPeelUtil.h ObitUVPeelUtil 
 * \endlink modules.
 *
 * This class contains interoferometric data and allows access.
 * An ObitUV is the front end to a persistent disk resident structure.
 * There maybe (usually are) associated tables which either describe
 * the data or contain calibration and/or editing information.
 * These associated tables are listed in an #ObitTableList member and
 * the #newObitUVTable function allows access to these tables.
 * Both FITS (as Tables) and AIPS cataloged data are supported.
 * The knowledge of underlying classes should be limited to private 
 * function #ObitUVSetupIO in ObitUV.c
 *
 * \section ObitUVSpecification Specifying desired data transfer parameters
 * The desired data transfers are specified in the member ObitInfoList.
 * There are separate sets of parameters used to specify the FITS or AIPS 
 * data files.
 * Data is read and written as arrays of floats, data compressed on the 
 * disk is compressed/uncompressed on the fly.
 * In the following an ObitInfoList entry is defined by 
 * the name in double quotes, the data type code as an #ObitInfoType enum 
 * and the dimensions of the array (? => depends on application).
 * To specify whether the underlying data files are FITS or AIPS
 * \li "FileType" OBIT_int (1,1,1) OBIT_IO_FITS or OBIT_IO_AIPS 
 * which are values of an #ObitIOType enum defined in ObitTypes.h.
 *
 * The following apply to both types of files:
 * \li "nVisPIO", OBIT_int, Max. Number of visibilities per 
 *     "Read" or "Write" operation.  Default = 1.
 *
 * \subsection UVFITS FITS files
 * This implementation uses cfitsio which allows using, in addition to 
 * regular FITS data, gzip compressed files, pipes, shared memory 
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
 * The convenience macro #ObitUVSetAIPS simplifies specifying the 
 * desired data.
 * For accessing AIPS files, the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 *
 * \section ObitUVaccess Creators and Destructors
 * An ObitUV can be created using newObitUV which allows specifying 
 * a name for the object.  This name is used to label messages.
 * The copy constructors #ObitUVClone and ObitUVCopy make shallow
 * and deep copies of an extant ObitUV.  If the output ObitUV has
 * previously been specified, including its disk resident information,
 * then ObitUVCopy will copy the disk resident as well as the memory 
 * resident information.  Also, any associated tables will be copied.
 *
 * A copy of a pointer to an ObitUV should always be made using the
 * #ObitUVRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUV or changing a pointer, the function
 * #ObitUVUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitUVUsage I/O
 * Visibility data is available after an input object is "Opened"
 * and "Read".
 * "Read Select" also allows specifying the data to be read as well as
 * optional calibration and editing to be applied as the data is read.
 * I/O optionally uses a buffer attached to the ObitUV or some external
 * location.
 * Data consists of a set of "random parameters" (u,v,w time, baseline, etc)
 * and a rectangular data array of complex visibilities with a weight.
 * The order, presence and size of components of the data are described
 * in an #ObitUVDesc object which also tells which visibility numbers are in 
 * the buffer.
 * To Write an ObitUV, create it, open it, and write.
 * The object should be closed to ensure all data is flushed to disk.
 * Deletion of an ObitUV after its final unreferencing will automatically
 * close it.
 *
 * \section Data Selection, Editing and Calibration
 * All IO supports (where appropriate) data selection, editing an calibration.
 * These are controlled by information on the ObitUV data object's info member,
 * details are given in the #ObitUVSel class documentation.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUV Class structure. */
typedef struct {
#include "ObitUVDef.h"   /* this class definition */
} ObitUV;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUV
 * returns a ObitUV*.
 * in = object to unreference
 */
#define ObitUVUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUV.
 * returns a ObitUV*.
 * in = object to reference
 */
#define ObitUVRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVIsA(in) ObitIsA (in, ObitUVGetClass())

/** 
 * Convenience Macro to define UV I/O to a FITS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitUV to specify i/O for.
 *\li nvis = Max. Number of visibilities per read.
 *\li disk = FITS disk number
 *\li file = Specified FITS file name.
 *\li err = ObitErr to receive error messages.
 */
#define ObitUVSetFITS(in,nvis,disk,file,err)  G_STMT_START{         \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_FITS;                            \
       in->info->work[1] = nvis; in->info->work[2]= disk;           \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       ObitInfoListPut (in->info, "nVisPIO", OBIT_long,              \
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

/** 
 * Convenience Macro to define UV I/O to an AIPS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitUV to specify i/O for.
 *\li nvis = Max. Number of visibilities per read.
 *\li disk = AIPS disk number
 *\li cno  = catalog slot number
 *\li user = User id number
 *\li err = ObitErr to receive error messages.
 */
#define ObitUVSetAIPS(in,nvis,disk,cno,user,err)  G_STMT_START{     \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_AIPS;                            \
       in->info->work[1]= nvis;                                     \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       ObitInfoListPut (in->info, "nVisPIO", OBIT_long, in->info->dim,\
		 (gpointer)&in->info->work[1], err);                \
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


/** 
 * Divide one complex number by another.
 *\li in1  = Numerator complex (real,imaginary)
 *\li in2  = Denominator complex
 *\li out  = Output complex value, can be in1
 *           (0,0) on zero divide
 *\li work = Array of 3 elements like in...
 */
#define ObitUVCpxDivide(in1,in2,out,work)  G_STMT_START{      \
       work[2] = in2[0]*in2[0] + in2[1]*in2[1];               \
       if (work[2]==0.0) work[2] = 1;                         \
       work[0] = in1[0]/work[2]; work[1] = in1[1]/work[2];    \
       out[0] = work[0]*in2[0] + work[1]*in2[1];              \
       out[1] = work[1]*in2[0] - work[0]*in2[1];              \
     }G_STMT_END  

/** 
 * Divide one complex number with weight by another.
 * Sets values on ObitInfoList on input object.
 *\li in1  = Numerator complex (real,imaginary,weight)
 *\li in2  = Denominator complex
 *\li out  = Output complex value, can be in1, 
 *           (0,0,0) on zero divide
 *\li work = Array of 3 elements like in...
 */
#define ObitUVWtCpxDivide(in1,in2,out,work)  G_STMT_START{    \
       if ((in1[2]<=0.0) || (in2[2]<=0.0)) { /* bad */        \
         out[0] = out[1] = out[2] = 0.0;                      \
       } else { /* do division */                             \
         work[2] = in2[0]*in2[0] + in2[1]*in2[1];             \
         if (work[2]==0.0) {out[0] = out[1] = out[2] = 0.0;   \
         } else {  /* OK */                                   \
           work[0] = in1[0]/work[2]; work[1] = in1[1]/work[2];\
           out[0] = work[0]*in2[0] + work[1]*in2[1];          \
           out[1] = work[1]*in2[0] - work[0]*in2[1];          \
           out[2] *= sqrt(work[2]);                           \
         }                                                    \
       }                                                      \
     }G_STMT_END  

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVClassInit (void);

/** Public: Constructor. */
ObitUV* newObitUV (gchar* name);

/** Public: Copy Constructor for scratch file. */
ObitUV* newObitUVScratch (ObitUV *in, ObitErr *err);
typedef ObitUV* (*newObitUVScratchFP) (ObitUV *in, ObitErr *err);

/** Public: Fully instantiate. */
void ObitUVFullInstantiate (ObitUV *in, gboolean exist, ObitErr *err);
typedef void (*ObitUVFullInstantiateFP) (ObitUV *in, gboolean exist, 
					 ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitUVGetClass (void);

/** Public: Delete underlying structures. */
ObitUV* ObitUVZap  (ObitUV *in, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitUVRename  (ObitUV *in, ObitErr *err);

/** Public: Copy (deep) constructor. */
ObitUV* ObitUVCopy  (ObitUV *in, ObitUV *out, ObitErr *err);

/** Public: Copy structure. */
void ObitUVClone (ObitUV *in, ObitUV *out, ObitErr *err );

/** Public: Do two UVs have the same underlying structures?. */
gboolean ObitUVSame (ObitUV *in1, ObitUV *in2, ObitErr *err );
typedef gboolean (*ObitUVSameFP) (ObitUV *in1, ObitUV *in2, 
				  ObitErr *err);

/** Public: Create ObitIO structures and open file */
ObitIOCode ObitUVOpen (ObitUV *in, ObitIOAccess access, 
			  ObitErr *err);

/** Public: Close file and become inactive */
ObitIOCode ObitUVClose (ObitUV *in, ObitErr *err);

/** Public: Reset IO to start of file */
ObitIOCode ObitUVIOSet (ObitUV *in, ObitErr *err);

/** Public: Read specified data */
ObitIOCode ObitUVRead (ObitUV *in, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitUVReadFP) (ObitUV *in, ofloat *data, 
				    ObitErr *err);

/** Public: Read select, edit, calibrate specified data */
ObitIOCode ObitUVReadSelect (ObitUV *in, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitUVReadSelectFP) (ObitUV *in, ofloat *data, 
					  ObitErr *err);

/** Public: Write specified data */
ObitIOCode ObitUVWrite (ObitUV *in, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitUVWriteFP) (ObitUV *in, ofloat *data, ObitErr *err);

/** Public: Rewrite specified data */
ObitIOCode ObitUVRewrite (ObitUV *in, ofloat *data, ObitErr *err);
typedef ObitIOCode (*ObitUVRewritefp) (ObitUV *in, ofloat *data, ObitErr *err);

/** Public: Return an associated Table */
ObitTable* newObitUVTable (ObitUV *in, ObitIOAccess access, 
			   gchar *tabType, olong *tabver, ObitErr *err);
typedef ObitTable* (*newObitUVTableFP) (ObitUV *in, ObitIOAccess access, 
					gchar *tabType, olong *tabver, 
					ObitErr *err);

/** Public: Destroy an associated Table */
ObitIOCode ObitUVZapTable (ObitUV *in, gchar *tabType, olong tabVer, 
			   ObitErr *err);
typedef ObitIOCode (*ObitUVZapTableFP) (ObitUV *in, gchar *tabType, 
					olong tabVer, ObitErr *err);

/** Public: Copy associated Tables */
ObitIOCode ObitUVCopyTables (ObitUV *in, ObitUV *out, gchar **exclude,
			     gchar **include, ObitErr *err);
typedef ObitIOCode (*ObitUVCopyTablesFP) (ObitUV *in, ObitUV *out, 
					   gchar **exclude, gchar **include, 
					   ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitUVUpdateTables (ObitUV *in, ObitErr *err);
typedef ObitIOCode (*ObitUVUpdateTablesFP) (ObitUV *in, ObitErr *err);

/** Public: Get Frequency arrays */
void ObitUVGetFreq (ObitUV* in, ObitErr *err);

/** Public: Obtains Subarray info for an ObitUV */
ObitIOCode ObitUVGetSubA (ObitUV *in, ObitErr *err);

/** Public: Get source position */
void ObitUVGetRADec (ObitUV *uvdata, odouble *ra, odouble *dec, 
		     ObitErr *err);

/** Public: Get single source info */
void ObitUVGetSouInfo (ObitUV *uvdata, ObitErr *err);

/** Public: Write header keyword */
void ObitUVWriteKeyword (ObitUV *in, 
			 gchar* name, ObitInfoType type, gint32 *dim, 
			 gconstpointer data, ObitErr *err);
/** Public: Read header keyword */
void ObitUVReadKeyword (ObitUV *in, 
			gchar* name, ObitInfoType *type, gint32 *dim, 
			gpointer data, ObitErr *err);

/** Public: Channel selection in FG table */
olong ObitUVChanSel (ObitUV *in, gint32 *dim, olong *IChanSel, ObitErr *err);
typedef olong (*ObitUVChanSelFP) (ObitUV *in, gint32 *dim, olong *IChanSel, 
				 ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVClassDef.h"
} ObitUVClassInfo; 

#endif /* OBITUV_H */ 
