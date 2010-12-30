/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2010                                          */
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
#ifndef OBITTABLEHISTORY_H 
#define OBITTABLEHISTORY_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIOHistory.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitHistory.h
 * ObitHistory class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class contains a processing history in tabular form.
 * An ObitHistory is the front end to a persistent disk resident structure.
 * Both FITS (as Tables) and AIPS cataloged data are supported.
 * This class is derived from the Obit class.
 * The AIPS conventions for history records are give the program name followed
 * by parameters in keyword=value form and to preceed any non parsable text
 * by the FITS comment header comment delimiter '/'.
 *
 * \section TableDataStorage HistoryTable data storage
 * History tables are accessed as tables although the AIPS implementation 
 * is a pre-table version.  History records are blocked into 70 character 
 * fixed strings althought AIPS internally uses 72.
 *
 * \section ObitHistorySpecification Access to History
 * History records are stored in a system dependent fashion.
 * AIPS history records are stored in an AIPS HI* file with 72 characters 
 * per record (Obit only uses 70).
 * FOR FITS files, history records are normally kept in A History binary table
 * but can be read or written to tyhe more traditional HISTORY keywords using
 * #ObitHistoryCopyHeader, or #ObitHistoryCopy2Header for entire collections
 * or #ObitFileFITS functions #ObitFileFITSReadHistory and #ObitFileFITSWriteHistory
 * for individual records.
 * Access to the history component of an object (e.g. Image, UV data) can
 * be obtained using the ObitInfoList containing the information defining the 
 * underlying file.
 * This uses routine #newObitHistoryValue.
 * Then history lines can be read or written one at a time using #ObitHistoryOpen, 
 * #ObitHistoryReadRec, #ObitHistoryWriteRec, #ObitHistoryTimeStamp, #ObitHistoryClose.
 * The contents of entire history files may be copied  using #ObitHistoryCopy, or
 * #ObitHistoryCopyHeader to copy HISTORY records from a FITS header or
 * #ObitHistoryCopy2Header to copy a history file to a FITS header.
 *
 * \subsection TableFITS FITS files
 * The Obit FITS implementation uses a HISTORY table unlike the standard FITS 
 * practice of keeping histories in the main HDU header
 * regular FITS images, gzip compressed files, pipes, shared memory 
 * and a number of other input forms.
 * The convenience Macro #ObitHistorySetFITS simplifies specifying the 
 * desired data.
 * \li "Disk"     OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) FITS file name.
 *
 * \subsection ObitHistoryAIPS AIPS files
 * The ObitAIPS class must be initialized before accessing AIPS files; 
 * this uses #ObitAIPSClassInit.
 * The convenience Macro ObitHistorySetAIPS simplifies specifying the 
 * desired data.
 * For accessing AIPS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 *
 * \section ObitHistoryaccess Creators and Destructors
 *
 * A copy of a pointer to an ObitHistory should always be made using the
 * ObitHistoryRef function which updates the reference count in the object.
 * Then whenever freeing an ObitHistory or changing a pointer, the function
 * ObitHistoryUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 *
 * \section ObitHistoryUsage I/O
 * Visibility data is available after an input object is "Opened"
 * and "Read".
 * I/O optionally uses a buffer attached to the ObitHistory or some external
 * location.
 * To Write an ObitHistory, create it, open it, and write.
 * The object should be closed to ensure all data is flushed to disk.
 * Deletion of an ObitHistory after its final unreferencing will automatically
 * close it.
 *
 * \section ObitHistoryRow ObitHistoryRow
 * An ObitHistoryRow is used to pass the data for a single table row.
 * It is most useful in derived classes where it includes entries by 
 * symbolic names. ObitHistoryRow class definitions and functions are 
 * included in the files defining the associated ObitHistory.
 * ObitHistoryRows are derived from basal class Obit.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitHistory Class structure. */
typedef struct {
#include "ObitHistoryDef.h"   /* this class definition */
} ObitHistory;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitHistory
 * returns an ObitHistory*.
 * in = object to unreference
 */
#define ObitHistoryUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitHistory.
 * returns an ObitHistory*.
 * in = object to reference
 */
#define ObitHistoryRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitHistoryIsA(in) ObitIsA (in, ObitHistoryGetClass())

/** 
 * Convenience Macro to define History I/O to a FITS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = Obit Object with InfoList member to specify i/O for.
 *\li disk = FITS disk number
 *\li file = Specified FITS file name.
 *\li err = ObitErr to receive error messages.
 */
#define ObitHistorySetFITS(in,disk,file,err) G_STMT_START{ \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_FITS;                            \
       in->info->work[2] = disk;                                    \
       in->info->dim[0] = 1;                                        \
       ObitInfoListPut (in->info, "Disk", OBIT_long,                 \
                 in->info->dim, (gpointer)&in->info->work[2], err); \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       in->info->dim[0] = strlen(file);                             \
       ObitInfoListPut (in->info, "FileName", OBIT_string,          \
                 in->info->dim, (gpointer)file, err);               \
     }G_STMT_END  

/** 
 * Convenience Macro to define History I/O to an AIPS file.
 * Sets values on ObitInfoList on input object.
 *\li in   = Obit Object with InfoList member to specify i/O for.
 *\li disk = AIPS disk number
 *\li cno  = catalog slot number
 *\li user = User id number
 *\li err  = ObitErr to receive error messages.
 */
#define ObitHistorySetAIPS(in,disk,cno,user,err)  G_STMT_START{  \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0] = OBIT_IO_AIPS;                            \
       ObitInfoListPut (in->info, "FileType", OBIT_long,             \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
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
void ObitHistoryClassInit (void);

/** Public: Default constructor. */
ObitHistory* newObitHistory (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitHistoryGetClass (void);

/** Public: Constructor from object infoList */
ObitHistory* 
newObitHistoryValue (gchar* name, ObitInfoList *info, ObitErr *err);

/** Public: Delete underlying structures. */
ObitHistory* ObitHistoryZap (ObitHistory *in, ObitErr *err);
typedef ObitHistory*(*ObitHistoryZapFP) (ObitHistory *in, ObitErr *err);

/** Public: Deep copy. */
ObitHistory* ObitHistoryCopy  (ObitHistory *in, ObitHistory *out, 
			       ObitErr *err);

/** Public: Copy history from header (FITS). */
ObitIOCode ObitHistoryCopyHeader  (ObitHistory *in, ObitHistory *out, 
				   ObitErr *err);

/** Public: Copy history to header (FITS). */
ObitIOCode ObitHistoryCopy2Header  (ObitHistory *in, ObitHistory *out, 
				    ObitErr *err);

/** Public: Copy history from header (FITS) to header (FITS). */
ObitIOCode ObitHistoryHeader2Header  (ObitHistory *in, ObitHistory *out, 
				      ObitErr *err);

/** Public: Create ObitIO structures and open file */
ObitIOCode ObitHistoryOpen (ObitHistory *in, ObitIOAccess access, 
			    ObitErr *err);
typedef ObitIOCode(*ObitHistoryOpenFP) (ObitHistory *in, ObitIOAccess access,
					ObitErr *err);

/** Public: Close file and become inactive */
ObitIOCode ObitHistoryClose (ObitHistory *in, ObitErr *err);
typedef ObitIOCode(*ObitHistoryCloseFP) (ObitHistory *in,  ObitErr *err);

/** Public: Read specified Record */
ObitIOCode 
ObitHistoryReadRec (ObitHistory *in, olong recno, gchar hiCard[73],
		    ObitErr *err);

/** Public: Write specified Record */
ObitIOCode 
ObitHistoryWriteRec (ObitHistory *in, olong rowno, gchar hiCard[73], 
		     ObitErr *err);

/** Public: Add time stamp and label */
ObitIOCode 
ObitHistoryTimeStamp (ObitHistory *in, gchar *label,
			      ObitErr *err);

/** Public: Copy a list of values from an InfoList to a History */
ObitIOCode 
ObitHistoryCopyInfoList (ObitHistory *out, gchar *pgmName, gchar *list[], 
			 ObitInfoList *info, ObitErr *err);

/** Public: Tell number of history records */
olong ObitHistoryNumRec (ObitHistory *in);

/** Public: Remove a range of records */
ObitIOCode ObitHistoryEdit (ObitHistory *in, olong startr, olong endr,
			    ObitErr *err);
typedef ObitIOCode(*ObitHistoryEditFP) (ObitHistory *in, 
					olong startr, olong endr,
					ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure For Table.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitHistoryClassDef.h"
} ObitHistoryClassInfo; 

#endif /* OBITTABLEHISTORY_H */ 
