/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2024                                              */
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
#ifndef OBITTABLEJI_H 
#define OBITTABLEJI_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTable.h"
#include "ObitData.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableJI.h
 * ObitTableJI class definition.
 *
 * This class is derived from the #ObitTable class.
 *
 * This class contains tabular data and allows access.
 * "AIPS JI" Jones matrix calibration information derived from a calibration solution 
 * and can be used either for "self calibration" or correction to a multisource 
 * Total Jones (JT) calibration table.  Jones matrices are 2x2 complex matrices
 * stored in a linearized form as floats.
 * 
 * Calibration of the visibility on a baseline 1-2 can use the following:\hfill\break
 * J1 = Jones matrix for first antenna,\hfill\break
 * J2 = conjugate transpose of the Jones matrix for second antenna,\hfill\break
 * vis = 
 * \beginverbatim
 *
 * This class contains tabular data and allows access.
 * "AIPS JI" table
 * An ObitTableJI is the front end to a persistent disk resident structure.
 * Both FITS (as Tables) and AIPS cataloged data are supported.
 *
 * \section TableDataStorage Table data storage
 * In memory tables are stored in a fashion similar to how they are 
 * stored on disk - in large blocks in memory rather than structures.
 * Due to the word alignment requirements of some machines, they are 
 * stored by order of the decreasing element size: 
 * double, float long, int, short, char rather than the logical order.
 * The details of the storage in the buffer are kept in the 
 * #ObitTableJIDesc.
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
 * \section ObitTableJISpecification Specifying desired data transfer parameters
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
 * The convenience Macro #ObitTableJISetFITS simplifies specifying the 
 * desired data.
 * Binary tables are used for storing visibility data in FITS.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "FileName" OBIT_string (?,1,1) FITS file name.
 * \li "TabName"  OBIT_string (?,1,1) Table name (e.g. "AIPS CC").
 * \li "Ver"      OBIT_int    (1,1,1) Table version number
 *
 * subsection ObitTableJIAIPS AIPS files
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * For accessing AIPS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 * \li "TableType" OBIT_string (2,1,1) AIPS Table type
 * \li "Ver"  OBIT_int    (1,1,1) AIPS table version number.
 *
 * \section ObitTableJIaccess Creators and Destructors
 * An ObitTableJI can be created using newObitTableJIValue which attaches the 
 * table to an ObitData for the object.  
 * If the output ObitTableJI has previously been specified, including file information,
 * then ObitTableJICopy will copy the disk resident as well as the memory 
 * resident information.
 *
 * A copy of a pointer to an ObitTableJI should always be made using the
 * ObitTableJIRef function which updates the reference count in the object.
 * Then whenever freeing an ObitTableJI or changing a pointer, the function
 * ObitTableJIUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 *
 * \section ObitTableJIUsage I/O
 * Visibility data is available after an input object is "Opened"
 * and "Read".
 * I/O optionally uses a buffer attached to the ObitTableJI or some external
 * location.
 * To Write an ObitTableJI, create it, open it, and write.
 * The object should be closed to ensure all data is flushed to disk.
 * Deletion of an ObitTableJI after its final unreferencing will automatically
 * close it.
 */

/*--------------Class definitions-------------------------------------*/

/** Number of characters for Table keyword */
// Defensive move, 24 likely to cause trouble, fixed by hand
 #define MAXKEYCHARTABLEJI 9

/** ObitTableJI Class structure. */
typedef struct {
#include "ObitTableJIDef.h"   /* this class definition */
} ObitTableJI;

/** ObitTableJIRow Class structure. */
typedef struct {
#include "ObitTableJIRowDef.h"   /* this class row definition */
} ObitTableJIRow;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTableJI
 * returns an ObitTableJI*.
 * in = object to unreference
 */
#define ObitTableJIUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTableJI.
 * returns an ObitTableJI*.
 * in = object to reference
 */
#define ObitTableJIRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTableJIIsA(in) ObitIsA (in, ObitTableJIGetClass())

/** 
 * Macro to unreference (and possibly destroy) an ObitTableJIRow
 * returns an ObitTableJIRow*.
 * in = object to unreference
 */
#define ObitTableJIRowUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTableJIRow.
 * returns an ObitTableJIRow*.
 * in = object to reference
 */
#define ObitTableJIRowRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTableJIRowIsA(in) ObitIsA (in, ObitTableJIRowGetClass())

/*---------------Public functions---------------------------*/
/*----------------Table Row Functions ----------------------*/
/** Public: Row Class initializer. */
void ObitTableJIRowClassInit (void);

/** Public: Constructor. */
ObitTableJIRow* newObitTableJIRow (ObitTableJI *table);

/** Public: ClassInfo pointer */
gconstpointer ObitTableJIRowGetClass (void);

/*------------------Table Functions ------------------------*/
/** Public: Class initializer. */
void ObitTableJIClassInit (void);

/** Public: Constructor. */
ObitTableJI* newObitTableJI (gchar* name);

/** Public: Constructor from values. */
ObitTableJI* 
newObitTableJIValue (gchar* name, ObitData *file, olong *ver,
  		     ObitIOAccess access,
                     oint numIF, oint numChan,
		     ObitErr *err);

/** Public: Class initializer. */
void ObitTableJIClassInit (void);

/** Public: ClassInfo pointer */
gconstpointer ObitTableJIGetClass (void);

/** Public: Copy (deep) constructor. */
ObitTableJI* ObitTableJICopy  (ObitTableJI *in, ObitTableJI *out, 
			   ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitTableJI* ObitTableJIClone (ObitTableJI *in, ObitTableJI *out);

/** Public: Convert an ObitTable to an ObitTableJI */
ObitTableJI* ObitTableJIConvert  (ObitTable *in);

/** Public: Create ObitIO structures and open file */
ObitIOCode ObitTableJIOpen (ObitTableJI *in, ObitIOAccess access, 
			  ObitErr *err);

/** Public: Read a table row */
ObitIOCode 
ObitTableJIReadRow  (ObitTableJI *in, olong iJIRow, ObitTableJIRow *row,
		     ObitErr *err);

/** Public: Init a table row for write */
void 
ObitTableJISetRow  (ObitTableJI *in, ObitTableJIRow *row,
		     ObitErr *err);

/** Public: Write a table row */
ObitIOCode 
ObitTableJIWriteRow  (ObitTableJI *in, olong iJIRow, ObitTableJIRow *row,
		     ObitErr *err);

/** Public: Close file and become inactive */
ObitIOCode ObitTableJIClose (ObitTableJI *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTableJIClassDef.h"
} ObitTableJIClassInfo; 

/**
 * ClassInfo Structure For TableJIRow.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTableJIRowClassDef.h"
} ObitTableJIRowClassInfo; 
#endif /* OBITTABLEJI_H */ 