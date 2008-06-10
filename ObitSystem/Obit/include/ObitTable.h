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
#ifndef OBITTABLE_H 
#define OBITTABLE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitTableDesc.h"
#include "ObitTableSel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTable.h
 * ObitTable class definition.
 *
 * This class contains tabular data and allows access.
 * An ObitTable is the front end to a persistent disk resident structure.
 * Both FITS (as Tables) and AIPS cataloged data are supported.
 *
 * This class is derived from the #Obit class.
 * Related functions are in the 
 * \link ObitImageUtil.h ObitImageUtil 
 * \endlink module and table type specific classes and utility modules.
 *
 * \section TableDataStorage Table data storage
 * In memory tables are stored in a fashion similar to how they are 
 * stored on disk - in large blocks in memory rather than structures.
 * Due to the word alignment requirements of some machines, they are 
 * stored by order of the decreasing element size: 
 * double, float long, int, short, char rather than the logical order.
 * The details of the storage in the buffer are kept in the 
 * #ObitTableDesc.
 *
 * In addition to the normal tabular data, a table will have a "_status"
 * column to indicate the status of each row.
 * The status value is read from and written to (some modification) AIPS 
 * tables but are not written to externally generated FITS tables which
 * don't have these colummns.  It will be written to Obit generated tables
 * which will have these columns.
 * Status values:
 * \li status = 0 => normal
 * \li status = 1 => row has been modified (or created) and needs to be
 *                   written.
 * \li status = -1 => row has been marked invalid.
 *
 * \section ObitTableSpecification Specifying desired data transfer parameters
 * The desired data transfers are specified in the member ObitInfoList.
 * There are separate sets of parameters used to specify the FITS or AIPS 
 * data files.
 * In the following an ObitInfoList entry is defined by 
 * the name in double quotes, the data type code as an #ObitInfoType enum 
 * and the dimensions of the array (? => depends on application).
 * To specify whether the underlying data files are FITS or AIPS
 * \li "FileType" OBIT_int (1,1,1) OBIT_IO_FITS or OBIT_IO_AIPS 
 * which are values of an #ObitIOType enum defined in ObitIO.h.
 *
 * The following apply to both types of files:
 * \li "nRowPIO", OBIT_int, Max. Number of visibilities per 
 *     "Read" or "Write" operation.  Default = 1.
 *
 * \subsection TableFITS FITS files
 * This implementation uses cfitsio which allows using, in addition to 
 * regular FITS images, gzip compressed files, pipes, shared memory 
 * and a number of other input forms.
 * The convenience Macro #ObitTableSetFITS simplifies specifying the 
 * desired data.
 * Binary tables are used for storing visibility data in FITS.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk"     OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) FITS file name.
 * \li "TabName"  OBIT_string (?,1,1) Table name (e.g. "AIPS CC").
 * \li "Ver"      OBIT_int    (1,1,1) Table version number
 *
 * \subsection ObitTableAIPS AIPS files
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * The convenience Macro ObitTableSetAIPS simplifies specifying the 
 * desired data.
 * For accessing AIPS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 * \li "TableType" OBIT_string (2,1,1) AIPS Table type
 * \li "Ver"  OBIT_int    (1,1,1) AIPS table version number.
 *
 * \section ObitTableaccess Creators and Destructors
 * An ObitTable can be created using newObitTable which allos specifying 
 * a name for the object.  This name is used to label messages.
 * The copy constructors ObitTableClone and ObitTableCopy make shallow
 * and deep copies of an extant ObitTable.  If the output ObitTable has
 * previously been specified, including its disk resident information,
 * then ObitTableCopy will copy the disk resident as well as the memory 
 * resident information.
 *
 * A copy of a pointer to an ObitTable should always be made using the
 * ObitTableRef function which updates the reference count in the object.
 * Then whenever freeing an ObitTable or changing a pointer, the function
 * ObitTableUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 *
 * \section ObitTableUsage I/O
 * Visibility data is available after an input object is "Opened"
 * and "Read".
 * I/O optionally uses a buffer attached to the ObitTable or some external
 * location.
 * To Write an ObitTable, create it, open it, and write.
 * The object should be closed to ensure all data is flushed to disk.
 * Deletion of an ObitTable after its final unreferencing will automatically
 * close it.
 *
 * \section ObitTableRow ObitTableRow
 * An ObitTableRow is used to pass the data for a single table row.
 * It is most useful in derived classes where it includes entries by 
 * symbolic names. ObitTableRow class definitions and functions are 
 * included in the files defining the associated ObitTable.
 * ObitTableRows are derived from basal class Obit.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitTable Class structure. */
typedef struct {
#include "ObitTableDef.h"   /* this class definition */
} ObitTable;

/** ObitTableRow Class structure. */
typedef struct {
#include "ObitTableRowDef.h"   /* this class row definition */
} ObitTableRow;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTable
 * returns an ObitTable*.
 * in = object to unreference
 */
#define ObitTableUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTable.
 * returns an ObitTable*.
 * in = object to reference
 */
#define ObitTableRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTableIsA(in) ObitIsA (in, ObitTableGetClass())

/** 
 * Macro to unreference (and possibly destroy) an ObitTableRow
 * returns an ObitTableRow*.
 * in = object to unreference
 */
#define ObitTableRowUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTableRow.
 * returns an ObitTableRow*.
 * in = object to reference
 */
#define ObitTableRowRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTableRowIsA(in) ObitIsA (in, ObitTableRowGetClass())

/** 
 * Convenience Macro to define Table I/O to a FITS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitTable to specify i/O for.
 *\li disk = FITS disk number
 *\li file = Specified FITS file name.
 *\li tab  = table type (e.g. 'AIPS CC')
 *\li ver  = table version number
 *\li nrow = Maximum number of Rows per read/write.
 *\li err = ObitErr to receive error messages.
 */
#define ObitTableSetFITS(in,disk,file,tab,ver,nrow,err) G_STMT_START{ \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_FITS;                            \
       in->info->work[1] = ver; in->info->work[2]= nrow;            \
       in->info->work[2] = disk;                                    \
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&in->info->work[2], err); \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       ObitInfoListPut (in->info, "Ver", OBIT_long,                  \
                 in->info->dim, (gpointer)&in->info->work[1], err); \
       ObitInfoListPut (in->info, "nRowPIO", OBIT_long,              \
                 in->info->dim, (gpointer)&in->info->work[2], err); \
       in->info->dim[0] = strlen(tab);                              \
       ObitInfoListPut (in->info, "TabName", OBIT_string,           \
                 in->info->dim, (gpointer)tab, err);                \
       in->info->dim[0] = strlen(file);                             \
       ObitInfoListPut (in->info, "FileName", OBIT_string,          \
                 in->info->dim, (gpointer)file, err);               \
     }G_STMT_END  

/** 
 * Convenience Macro to define Table I/O to an AIPS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = ObitTable to specify i/O for.
 *\li disk = AIPS disk number
 *\li cno  = catalog slot number
 *\li tab  = table type (e.g. 'CC')
 *\li ver  = table version number
 *\li user = User id number
 *\li nrow = Maximum number of Rows per read/write.
 *\li err  = ObitErr to receive error messages.
 */
#define ObitTableSetAIPS(in,disk,cno,tab,ver,user,nrow,err)  G_STMT_START{  \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_AIPS;                            \
       in->info->work[1]= ver;  in->info->work[2]= nrow;            \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&disk, err);              \
       ObitInfoListPut (in->info, "Ver", OBIT_long,                  \
                 in->info->dim, (gpointer)&in->info->work[1], err); \
       ObitInfoListPut (in->info, "CNO", OBIT_long,                  \
                 in->info->dim, (gpointer)&cno, err);               \
       ObitInfoListPut (in->info, "User", OBIT_long,                 \
                 in->info->dim, (gpointer)&user, err);              \
       ObitInfoListPut (in->info, "nRowPIO", OBIT_long,              \
                 in->info->dim, (gpointer)&in->info->work[2], err); \
       in->info->dim[0] = 2;                                        \
       ObitInfoListPut (in->info, "TableType", OBIT_string,         \
                 in->info->dim, (gpointer)&tab, err);               \
     }G_STMT_END   


/*---------------Public functions---------------------------*/
/*----------------Table Row Functions ----------------------*/
/** Public: Row Class initializer. */
void ObitTableRowClassInit (void);

/** Public: Constructor. */
ObitTableRow* newObitTableRow (ObitTable *table);
/** Typedef for definition of class pointer structure */
typedef ObitTableRow* (*newObitTableRowFP) (ObitTable *in);

/** Public: ClassInfo pointer */
gconstpointer ObitTableRowGetClass (void);

/*------------------Table Functions ------------------------*/
/** Public: Class initializer. */
void ObitTableClassInit (void);

/** Public: Default constructor. */
ObitTable* newObitTable (gchar* name);

/** Public: Fully instantiate. */
void ObitTableFullInstantiate (ObitTable *in, gboolean exist, ObitErr *err);
typedef void (*ObitTableFullInstantiateFP) (ObitTable *in, gboolean exist, 
					    ObitErr *err);
/** Public: Remove any previous entries - also forces instantiate. */
void ObitTableClearRows (ObitTable *in, ObitErr *err);
typedef void (*ObitTableClearRowsFP) (ObitTable *in, ObitErr *err);
/** Public: ClassInfo pointer */
gconstpointer ObitTableGetClass (void);

/** Public: Delete underlying structures. */
ObitTable* ObitTableZap  (ObitTable *in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitTable* (*ObitTableZapFP) (ObitTable *in, ObitErr *err);

/** Public: Copy (deep) constructor. */
ObitTable* ObitTableCopy  (ObitTable *in, ObitTable *out, 
			   ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitTable* ObitTableClone (ObitTable *in, ObitTable *out);

/** Public: Concatenate two tables. */
void ObitTableConcat (ObitTable *in, ObitTable *out, ObitErr *err);

/** Public: Convert an ObitTable to a derived type*/
ObitTable* ObitTableConvert  (ObitTable *in);
typedef ObitTable* (*ObitTableConvertFP) (ObitTable *in);

/** Public: Create ObitIO structures and open file */
ObitIOCode ObitTableOpen (ObitTable *in, ObitIOAccess access, 
			  ObitErr *err);
typedef ObitIOCode (*ObitTableOpenFP) (ObitTable *in, ObitIOAccess access, 
			  ObitErr *err);

/** Public: Close file and become inactive */
ObitIOCode ObitTableClose (ObitTable *in, ObitErr *err);
typedef ObitIOCode (*ObitTableCloseFP) (ObitTable *in, ObitErr *err);

/** Public: Read specified data */
ObitIOCode ObitTableRead (ObitTable *in, olong row, ofloat *data, 
			  ObitErr *err);
typedef ObitIOCode (*ObitTableReadFP) (ObitTable *in, olong row, 
				       ofloat *data, ObitErr *err);

/** Public: Read/select data */
ObitIOCode ObitTableReadSelect (ObitTable *in, olong row, ofloat *data, 
				ObitErr *err);
typedef ObitIOCode (*ObitTableReadSelectFP) (ObitTable *in, olong row, 
					     ofloat *data, ObitErr *err);

/** Public: Attach Row to buffer */
void ObitTableSetRow  (ObitTable *in, ObitTableRow *row, ObitErr *err);
typedef void (*ObitTableSetRowFP) (ObitTable *in, ObitTableRow *row, 
				   ObitErr *err);

/** Public: Write specified data */
ObitIOCode ObitTableWrite (ObitTable *in, olong row, ofloat *data, 
			   ObitErr *err);
typedef ObitIOCode (*ObitTableWriteFP) (ObitTable *in, olong row, 
					ofloat *data, ObitErr *err);

/** Public: Read specified Row */
ObitIOCode ObitTableReadRow (ObitTable *in, olong rowno, ObitTableRow *row, 
			     ObitErr *err);
typedef ObitIOCode (*ObitTableReadRowFP) (ObitTable *in, olong rowno, 
					  ObitTableRow *row, ObitErr *err);

/** Public: Write specified Row */
ObitIOCode ObitTableWriteRow (ObitTable *in, olong rowno, ObitTableRow *row, 
			      ObitErr *err);
typedef ObitIOCode (*ObitTableWriteRowFP) (ObitTable *in, olong rowno, 
					   ObitTableRow *row, ObitErr *err);

/** Public: Return table type */
gchar* ObitTableGetType (ObitTable *in, ObitErr *err);
typedef gchar* (*ObitTableGetTypeFP) (ObitTable *in, ObitErr *err);

/** Public: Return table version */
olong ObitTableGetVersion (ObitTable *in, ObitErr *err);
typedef olong (*ObitTableGetVersionFP) (ObitTable *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure For Table.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTableClassDef.h"
} ObitTableClassInfo; 

/**
 * ClassInfo Structure For TableRow.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTableRowClassDef.h"
} ObitTableRowClassInfo; 

#endif /* OBITTABLE_H */ 
