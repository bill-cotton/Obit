/* $Id$           */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2016                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITIMAGE_H 
#define OBITIMAGE_H 

#include "ObitData.h"
#include "ObitImageDesc.h"
#include "ObitImageSel.h"
#include "ObitTableList.h"
#include "ObitFArray.h"
#include "ObitUVGrid.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImage.h
 * ObitImage class definition.
 *
 * This class is derived from the #ObitData class.
 * Related functions are in the 
 * \link ObitImageUtil.h ObitImageUtil 
 * \endlink ,
 * \link ObitConvUtil.h ObitConvUtil 
 * \endlink and
 * \link ObitFeatherUtil.h ObitFeatherUtil 
 * \endlink modules.
 *
 * This class contains an astronomical image and allows access.
 * An ObitImage is the front end to a persistent disk resident structure.
 * Magic value blanking is supported, blanked pixels have the value
 * returned by ObitMagicF().
 * There may be associated tables (e.g. "AIPS CC" tables).
 * These associated tables are listed in an #ObitTableList member and
 * the #newObitUVTable function allows access to these tables.
 * Both FITS and AIPS cataloged images are supported.
 *
 * \section ObitImageSpecification Specifying desired data transfer parameters
 * The desired data transfers are specified in the member ObitInfoList.
 * There are separate sets of parameters used to specify the FITS or AIPS 
 * data files.
 * In the following an ObitInfoList entry is defined by 
 * the name in double quotes, the data type code as an #ObitInfoType enum 
 * and the dimensions of the array (? => depends on application).
 * To specify whether the underlying data files are FITS or AIPS:
 * \li "FileType" OBIT_int (1,1,1) OBIT_IO_FITS or OBIT_IO_AIPS 
 * which are values of an #ObitIOType enum defined in ObitIO.h.
 *
 * The following apply to both types of files:
 * \li "BLC" OBIT_int (?,1,1) the bottom-left corner desired as expressed 
 * in 1-rel pixel indices.  If absent, the value (1,1,1...) will be assumed.
 * dimension of this array is [IM_MAXDIM].
 * \li "TRC" OBIT_int (?,1,1) the top-right corner desired as expressed 
 * in 1-rel pixel indices.  If absent, all pixels are included.
 * dimension of this array is [IM_MAXDIM].
 * \li "IOBy" OBIT_int (1,1,1) an ObitIOSize enum defined in ObitIO.h
 *  giving values OBIT_IO_byRow or  OBIT_IO_byPlane to specify 
 * if the data transfers  are to be by row or plane at a time.
 * Default is OBIT_IO_byRow.
 *
 * \subsection ImageFITS FITS files
 * This implementation uses cfitsio which allows using, in addition to 
 * regular FITS images, gzip compressed files, pipes, shared memory 
 * and a number of other input forms.
 * The convenience Macro #ObitImageSetFITS simplifies specifying the 
 * desired data.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) Name of disk file.
 *
 * The #ObitImageDesc member may contain:
 * \li "Quant" OBIT_float (1,1,1) Quantization level
 *           If given and > 0.0 and an integer output (Bitpix 16, 32) 
 *           is specified then the output will be quantized at this level.  
 *
 * \subsection ObitImageAIPS AIPS files
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * The convenience macro #ObitImageSetAIPS simplifies specifying the 
 * desired data.
 * For accessing AIPS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 *
 * \section ObitImageaccess Creators and Destructors
 * An ObitImage can be created using newObitImage which allows specifying 
 * a name for the object.  This name is used to label messages.
 * The copy constructors #ObitImageClone and ObitImageCopy make shallow
 * and deep copies of an extant #ObitImage.  If the output ObitImage has
 * previously been specified, including its disk resident information,
 * then #ObitImageCopy will copy the disk resident as well as the memory 
 * resident information.
 *
 * A copy of a pointer to an ObitImage should always be made using the
 * ObitImageRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImage or changing a pointer, the function
 * ObitImageUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 *
 * \section ObitImageUsage I/O
 * Pixel data in an image is available after an input object is "Opened"
 * and "Read".
 * I/O optionally uses a buffer attached to the ObitImage or some external
 * location.
 * To Write an ObitImage, create it, open it, and write.
 * The object should be closed to ensure all data is flushed to disk.
 * Deletion of an ObitImage after its final unreferencing will automatically
 * close it.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitImage Class structure. */
typedef struct {
#include "ObitImageDef.h"   /* this class definition */
} ObitImage;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImage
 * returns a ObitImage*.
 * in = object to unreference
 */
#define ObitImageUnref(in) ObitUnref ( in)

/** 
 * Macro to reference (update reference count) an ObitImage.
 * returns a ObitImage*.
 * in = object to reference
 */
#define ObitImageRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageIsA(in) ObitIsA (in, ObitImageGetClass())

/** 
 * Convenience Macro to define Image I/O to a FITS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitImage to specify i/O for.
 *\li size = size of I/O (OBIT_IO_byPlane or OBIT_IO_byRow).
 *\li disk = fits disk number
 *\li file = Specified FITS file name.
 *\li blc  = gint[IM_MAXDIM] giving bottom left corner (1-rel)
 *\li trc  = gint[IM_MAXDIM] giving top right corner (1-rel)
 *        0s => whole image
 *\li err = ObitErr to receive error messages.
 */
#define ObitImageSetFITS(in,size,disk,file,blc,trc,err) G_STMT_START{ \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_FITS;                            \
       in->info->work[1]= size; in->info->work[2]= disk;            \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       ObitInfoListPut (in->info, "IOBy", OBIT_long, in->info->dim,  \
		 (gpointer)&in->info->work[1], err);                \
       in->info->dim[0] = IM_MAXDIM;                                \
       ObitInfoListPut (in->info, "BLC", OBIT_long, in->info->dim,   \
		 (gpointer)blc, err);                               \
       ObitInfoListPut (in->info, "TRC", OBIT_long, in->info->dim,   \
		 (gpointer)trc, err);                               \
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&in->info->work[2], err); \
       in->info->dim[0] = strlen(file);                             \
       ObitInfoListPut (in->info, "FileName", OBIT_string,          \
                 in->info->dim, (gpointer)file, err);               \
     }G_STMT_END   

/** 
 * Convenience Macro to define Image I/O to an AIPS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitImage to specify i/O for.
 *\li size = size of I/O (OBIT_IO_byPlane or OBIT_IO_byRow).
 *\li disk = AIPS disk number
 *\li cno  = catalog slot number
 *\li user = User id number
 *\li blc  = gint[IM_MAXDIM] giving bottom left corner (1-rel)
 *\li trc  = gint[IM_MAXDIM] giving top right corner (1-rel)
 *        0s => whole image
 *\li err = ObitErr to receive error messages.
 */
#define ObitImageSetAIPS(in,size,disk,cno,user,blc,trc,err)  G_STMT_START{  \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_AIPS;                            \
       in->info->work[1]= size;                                     \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       ObitInfoListPut (in->info, "IOBy", OBIT_long, in->info->dim,  \
		 (gpointer)&in->info->work[1], err);                \
       in->info->dim[0] = IM_MAXDIM;                                \
       ObitInfoListPut (in->info, "BLC", OBIT_long, in->info->dim,   \
		 (gpointer)blc, err);                               \
       ObitInfoListPut (in->info, "TRC", OBIT_long, in->info->dim,   \
		 (gpointer)trc, err);                               \
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&disk, err);              \
       ObitInfoListPut (in->info, "CNO", OBIT_long,                  \
                 in->info->dim, (gpointer)&cno, err);               \
       ObitInfoListPut (in->info, "User", OBIT_long,                 \
                 in->info->dim, (gpointer)&user, err);              \
     }G_STMT_END   


/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageClassInit (void);

/** Public: Constructor. */
ObitImage* newObitImage (gchar* name);

/** Public: Create Image object from description in an ObitInfoList */
ObitImage* ObitImageFromFileInfo (gchar *prefix, ObitInfoList *inList, 
                                  ObitErr *err);
typedef ObitImage*
(*ObitImageFromFileInfoFP) (gchar *prefix, ObitInfoList *inList, 
                            ObitErr *err);

/** Public: Copy Constructor for scratch file. */
ObitImage* newObitImageScratch (ObitImage *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitImage* (*newObitImageScratchFP) (ObitImage *in, ObitErr *err);

/** Public: Fully instantiate. */
void ObitImageFullInstantiate (ObitImage *in, gboolean exist, ObitErr *err);
typedef void (*ObitImageFullInstantiateFP) (ObitImage *in, gboolean exist, 
					    ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitImageGetClass (void);

/** Public: Rename underlying structures. */
void ObitImageRename  (ObitImage *in, ObitErr *err);

/** Public: Delete underlying structures. */
ObitImage* ObitImageZap  (ObitImage *in, ObitErr *err);

/** Public: Rename underlying structures. */
void ObitImageRename  (ObitImage *in, ObitErr *err);

/** Public: Copy (deep) constructor. */
ObitImage* ObitImageCopy  (ObitImage *in, ObitImage *out, 
			   ObitErr *err);

/** Public: Copy structure. */
void ObitImageClone (ObitImage *in, ObitImage *out, ObitErr *err);

/** Public: Copy structure of in1 on grid of in2. */
void ObitImageClone2 (ObitImage *in1, ObitImage *in2, ObitImage *out, 
		      ObitErr *err);

/** Public: Do two Imagess have the same underlying structures?. */
gboolean ObitImageSame (ObitImage *in1, ObitImage *in2, ObitErr *err );
typedef gboolean (*ObitImageSameFP) (ObitImage *in1, ObitImage *in2, 
				  ObitErr *err);

/** Public: Copy structure of in to memory resident Image out */
void ObitImageCloneMem (ObitImage *in, ObitImage *out, ObitErr *err);

/** Public: Create ObitIO structures and open file */
ObitIOCode ObitImageOpen (ObitImage *in, ObitIOAccess access, 
			  ObitErr *err);

/** Public: Close file and become inactive */
ObitIOCode ObitImageClose (ObitImage *in, ObitErr *err);

/** Public: Read specified data */
ObitIOCode ObitImageRead (ObitImage *in, ofloat *data, ObitErr *err);

/** Public: Write specified data */
ObitIOCode ObitImageWrite (ObitImage *in, ofloat *data, ObitErr *err);

/** Public: Read specified image plane */
ObitIOCode ObitImageGetPlane (ObitImage *in, ofloat *data, olong plane[5], ObitErr *err);
typedef ObitIOCode (*ObitImageGetPlaneFP) (ObitImage *in, ofloat *data, 
					   olong plane[5], ObitErr *err);

/** Public: Write specified image plane */
ObitIOCode ObitImagePutPlane (ObitImage *in, ofloat *data, olong plane[5], ObitErr *err);
typedef ObitIOCode (*ObitImagePutPlaneFP) (ObitImage *in, ofloat *data, 
					   olong plane[5], ObitErr *err);

/** Public: Return an associated Table */
ObitTable* newObitImageTable (ObitImage *in, ObitIOAccess access, 
			      gchar *tabType, olong *tabver, ObitErr *err);
typedef ObitTable* (*newObitImageTableFP) (ObitImage *in, ObitIOAccess access, 
					   gchar *tabType, olong *tabver, 
					   ObitErr *err);

/** Public: Destroy an associated Table */
ObitIOCode ObitImageZapTable (ObitImage *in, gchar *tabType, olong tabVer, 
			   ObitErr *err);
typedef ObitIOCode (*ObitImageZapTableFP) (ObitImage *in, gchar *tabType, 
					olong tabVer, ObitErr *err);

/** Public: Copy associated Tables */
ObitIOCode ObitImageCopyTables (ObitImage *in, ObitImage *out, gchar **exclude,
			     gchar **include, ObitErr *err);
typedef ObitIOCode (*ObitImageCopyTablesFP) (ObitImage *in, ObitImage *out, 
					      gchar **exclude, gchar **include, 
					      ObitErr *err);

/** Public: Update disk resident tables information */
ObitIOCode ObitImageUpdateTables (ObitImage *in, ObitErr *err);
typedef ObitIOCode (*ObitImageUpdateTablesFP) (ObitImage *in, ObitErr *err);

/** Public: Set name etc for a beam associated with an image */
void ObitImageSetBeamName (ObitImage *image, ObitErr *err);
typedef void (*ObitImageSetBeamNameFP) (ObitImage *in, ObitErr *err);

/** Public: Write header keyword */
void ObitImageWriteKeyword (ObitImage *in, 
			    gchar* name, ObitInfoType type, gint32 *dim, 
			    gconstpointer data, ObitErr *err);
/** Public: Read header keyword */
void ObitImageReadKeyword (ObitImage *in, 
			   gchar* name, ObitInfoType *type, gint32 *dim, 
			   gpointer data, ObitErr *err);
/** Public: Set selection */
void ObitImageSetSelect (ObitImage *in, ObitIOSize IOBy, 
			 olong blc[IM_MAXDIM], olong trc[IM_MAXDIM],
			 ObitErr *err);
/** Public: Return image Beam. */
ObitImage* ObitImageGetBeam (ObitImage *in, olong beamNo, olong plane[5], 
			     ObitErr *err);
typedef ObitImage* (*ObitImageGetBeamFP) (ObitImage *in, olong beamNo,
					  olong plane[5], ObitErr *err);
/** Public: Return image Beam order. */
olong ObitImageGetBeamOrder (ObitImage *in);
typedef olong (*ObitImageGetBeamOrderFP) (ObitImage *in);

/** Public: Return Frequency of current plane */
odouble ObitImageGetPlaneFreq (ObitImage *in);
typedef double (*ObitImageGetPlaneFreqFP) (ObitImage *in);

/*----------- ClassInfo Structure -----------------------------------*/

/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageClassDef.h"
} ObitImageClassInfo; 

#endif /* OBITIMAGE_H */ 
