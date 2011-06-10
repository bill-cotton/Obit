/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2011                                          */
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
#ifndef OBITFILE_H 
#define OBITFILE_H 
#include <sys/types.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFile.h
 * ObitFile class definition for disk file access.
 *
 * This class is derived from the #Obit class.
 *
 * This class provides an I/O interface to disk files.
 *
 * \section ObitFileUsage Usage
 * Instances of this class are for access to disk files and is used 
 * for access to AIPS data files.
 * Instances can be made using the #newObitFile constructor,
 * or the #ObitFileCopy copy constructor and pointers copied 
 * (with reference pointer update) using #ObitFileRef.
 * The destructor (when reference count goes to zero) is
 * #ObitIOUnref.
 */

/*----------------- typedefs ---------------------------*/
/** Type for file position */
typedef off_t ObitFilePos;
/*typedef gint64 ObitFilePos;*/

/*---------------Class Structure---------------------------*/
/** ObitFile Class. */
typedef struct {
#include "ObitFileDef.h"   /* actual definition */
} ObitFile;

/*----------------- Macroes ---------------------------*/

/** 
 * Macro to unreference (and possibly destroy) an ObitFile
 * returns a ObitFile*.
 * in = object to unreference
 */
#define ObitFileUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFile.
 * returns a ObitFile*.
 * in = object to reference
 */
#define ObitFileRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFileIsA(in) ObitIsA (in, ObitFileGetClass())

/* Private functions are only defined in the .c file */

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFileClassInit (void);

/** Public: Constructor. */
ObitFile* newObitFile (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitFileGetClass (void);

/** Public:  destroy */
ObitFile* ObitFileZap (ObitFile *in, ObitErr *err);
typedef ObitFile* (*ObitFileZapFP) (ObitFile *in, ObitErr *err);

/** Public:  destroy file */
void ObitFileZapFile (gchar *fileName, ObitErr *err);
typedef void (*ObitFileZapFileFP) (gchar *fileName, ObitErr *err);

/** Public:  rename */
void ObitFileRename (gchar *oldName, gchar *newName, ObitErr *err);
typedef void (*ObitFileRenameFP) (gchar *oldName, gchar *newName,  ObitErr *err);

/** Public: Copy  constructor. */
ObitFile* ObitFileCopy  (ObitFile *in, ObitFile *out, ObitErr *err);

/** Public:  Open */ 
ObitIOCode 
ObitFileOpen (ObitFile *in, gchar *fileName, ObitIOAccess access,  
	      ObitIOType type, olong blockSize, ObitErr *err);
typedef ObitIOCode (*ObitFileOpenFP) (ObitFile *in, gchar *fileName, 
				      ObitIOAccess access,  
				      ObitIOType type, olong blockSize, 
				      ObitErr *err);

/** Public:  Close */
ObitIOCode ObitFileClose (ObitFile *in, ObitErr *err);
typedef ObitIOCode (*ObitFileCloseFP) (ObitFile *in, ObitErr *err);

/** Public:  Position  file */
ObitIOCode ObitFileSetPos (ObitFile *in, ObitFilePos filePos, ObitErr *err);
typedef ObitIOCode (*ObitFileSetPosFP) (ObitFile *in, ObitFilePos filePos, 
					ObitErr *err);

/** Public:  Position at end of file */
ObitIOCode ObitFileEnd (ObitFile *in, ObitErr *err);
typedef ObitIOCode (*ObitFileEndFP) (ObitFile *in, ObitErr *err);

/** Public:  Read */
ObitIOCode 
ObitFileRead (ObitFile *in, ObitFilePos filePos, olong size, gchar *buffer, 
	      ObitErr *err);
typedef ObitIOCode (*ObitFileReadFP) (ObitFile *in, ObitFilePos filePos, 
				      olong size,  gchar *buffer, 
				      ObitErr *err);

/** Public:  Read next line of text file */
ObitIOCode 
ObitFileReadLine (ObitFile *in, gchar *line, olong lineMax, ObitErr *err);
typedef ObitIOCode (*ObitFileReadLineFP) (ObitFile *in, gchar *line, olong lineMax, ObitErr *err);

/** Public:  Read next segment of XML file */
ObitIOCode 
ObitFileReadXML (ObitFile *in, gchar *line, olong lineMax, ObitErr *err);
typedef ObitIOCode (*ObitFileReadXMLFP) (ObitFile *in, gchar *line, olong lineMax, ObitErr *err);

/** Public:  Write */
ObitIOCode 
ObitFileWrite (ObitFile *in, ObitFilePos filePos, olong size, gchar *buffer, 
	       ObitErr *err);
typedef ObitIOCode (*ObitFileWriteFP) (ObitFile *in, ObitFilePos filePos, 
				       olong size, gchar *buffer, 
				       ObitErr *err);

/** Public:  Write  next line to text file */
ObitIOCode 
ObitFileWriteLine (ObitFile *in, gchar *line, ObitErr *err);
typedef ObitIOCode (*ObitFileWriteLineFP) (ObitFile *in, gchar *line, ObitErr *err);

/** Public: Pad remainder of a block. */
ObitIOCode 
ObitFilePad (ObitFile *in, olong padTo, gchar *buffer, olong bsize, 
	     ObitErr *err);
typedef ObitIOCode (*ObitFilePadFP) (ObitFile *in, olong padTo, 
				     gchar *buffer, 
				     olong bsize, ObitErr *err);

/** Public: Pad end of file to an integral of a given size. */
ObitIOCode 
ObitFilePadFile (ObitFile *in, olong blksize, ObitErr *err);
typedef ObitIOCode (*ObitFilePadFileFP) (ObitFile *in, olong blksize, 
					 ObitErr *err);

/** Public:  Flush Buffer */
ObitIOCode ObitFileFlush (ObitFile *in, ObitErr *err);
typedef ObitIOCode (*ObitFileFlushFP) (ObitFile *in, ObitErr *err);

/** Public: Does a given file exist? */
gboolean ObitFileExist (gchar *fileName, ObitErr *err);
typedef gboolean (*ObitFileExistFP) (gchar *fileName, ObitErr *err);

/** Public: What is the current size of a file */
ObitFilePos ObitFileSize (gchar *fileName, ObitErr *err);
typedef ObitFilePos (*ObitFileSizeFP) (ObitFile *in, ObitErr *err);

/** Public: What is the name (without path) if a file */
gchar* ObitFileName (gchar *fileName);
typedef gchar* (*ObitFileNameFP) (ObitFile *in);


/** Public: Utility, Error message for errno */
gboolean ObitFileErrMsg (ObitErr *err);
typedef gboolean (*ObitFileErrMsgFP) (ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFileClassDef.h" /* Actual definition */
} ObitFileClassInfo; 


#endif /* OBITFILE_H */ 
