/* $Id$        */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "ObitFile.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFile.c
 * ObitFile class function definitions.
 *
 * This Class allows access to the underlying file system.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitFile";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitFileClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFileClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFileInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFileClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFileClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitFile* newObitFile (gchar* name)
{
  ObitFile* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFileClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFile));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFileInit((gpointer)out);

  return out;
} /* end newObitFile */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitFileGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFileClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Delete the file.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return NULL as value of pointer to in, on failure returns in.
 */
ObitFile* ObitFileZap (ObitFile *in, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  olong status;
  gchar *routine = "ObitFileZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return in;
  if (in->fileName==NULL) return in;
  errno = 0;  /* reset any system error */

 /* Close if still open */
  if ((in->status==OBIT_Modified) || (in->status==OBIT_Active)) {
    retCode = ObitFileClose (in, err);
    if (err->error) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, in);
  }

  /* Make sure it exists */
  if (!ObitFileExist(in->fileName, err)) return in;

  /* Zap file */
  status = remove (in->fileName);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR deleting file %s", in->name);
    ObitFileErrMsg(err); /* Existing system error? */
    return in;
  }
  
  /* delete the rest of the structure */
  in = ObitFileUnref(in); 
  return in;
} /* end ObitFileZap */

/**
 * Rename a file.
 * \param oldName Current name for file
 * \param newName New name for file
 * \param err     ObitErr for reporting errors.
 */
void ObitFileRename (gchar *oldName, gchar *newName, ObitErr *err)
{
  olong status;
  gchar *routine = "ObitFileRename";

  /* error checks */
  if (err->error) return;
  g_assert (oldName!=NULL);
  g_assert (newName!=NULL);
  errno = 0;  /* reset any system error */

  /* Make sure it exists */
  if (!ObitFileExist(oldName, err)){
    Obit_log_error(err, OBIT_Error, 
		   "%s: File %s does not exist",
		   routine, oldName);
    return;
  }

  /* rename file */
  status = rename (oldName, newName);
    if (status) {             /* it went wrong */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR renaming file %s to %s", oldName, newName);
      ObitFileErrMsg(err);     /* system error message*/
      return;
    }
    errno = 0; /* In case */
 
} /* end ObitFileRename */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * The contents of any files are not modified.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitFile* ObitFileCopy  (ObitFile *in, ObitFile *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));
  errno = 0;  /* reset any system error */


  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitFile(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->access    = OBIT_IO_None; /* Not currently active */
  out->type      = in->type;
  out->exist     = in->exist;
  out->blockSize = in->blockSize;
  out->filePos   = in->filePos;
  out->status    = in->status;
  out->myFile    = in->myFile;
  if (out->fileName) g_free(out->fileName);
  out->fileName  = g_strdup(in->fileName);

  return out;
} /* end ObitFileCopy */

/**
 * Determine if a given file name exists.
 * \param fileName  Name of file to test.
 * \param err       ObitErr for reporting errors.
 * \return TRUE if exists, else FALSE.
 */
gboolean ObitFileExist (gchar *fileName, ObitErr *err)
{
  struct stat stbuf;
  olong status;  
  gboolean exist;

  /* error checks */
  if (err->error) return FALSE;
  g_assert (fileName!=NULL);
  errno = 0;  /* reset any system error */

  /* does it exist ? */
  status = stat(fileName, &stbuf);
  exist = (status==0);

  /* Clear status as nonexistance of a file is considered an error */
  if (!exist) errno = 0;

  return exist;
} /* end ObitFileExist */

/**
 * Determine file size
 * \param fileName  Name of file to test.
 * \param err       ObitErr for reporting errors.
 * \return File size in bytes
 */
ObitFilePos ObitFileSize (gchar *fileName, ObitErr *err)
{
  ObitFilePos out;
  struct stat stbuf;
  olong status;  

  /* error checks */
  if (err->error) return FALSE;
  g_assert (fileName!=NULL);
  errno = 0;  /* reset any system error */

  /* Get file info */
  status = stat(fileName, &stbuf);
  out = stbuf.st_size;

  return out;
} /* end ObitFileSize */

/**
 * Determine file name without path
 * \param fileName  Name of file
 * \return File size in bytes, g_free when done
 */
gchar* ObitFileName (gchar *fileName)
{
  gchar *out = NULL;
  olong i, last;  

  /* error checks */
  g_assert (fileName!=NULL);

  /* Find last occurance of directory seperator */
  last = -1;
  for (i=0; i<strlen(fileName); i++) 
    if (fileName[i]==G_DIR_SEPARATOR) last = i;

  /* Extract name */
  out = g_strdup(&fileName[last+1]);

  return out;
} /* end ObitFileName */

/**
 * Initialize structures and open file.
 * The file will be positioned at the beginning.
 * \param in        Pointer to object to be opened.
 * \param fileName  Name of file to open.
 * \param access    access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param type      type of file (OBIT_IO_Binary, OBIT_IO_Text).
 * \param blockSize Size of any blocking (used for AIPS files) 0->none.
 * \param err       ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileOpen (ObitFile *in, gchar *fileName, ObitIOAccess access,  
	      ObitIOType type, olong blockSize, ObitErr *err)
{
  ObitIOCode status, retCode = OBIT_IO_SpecErr;
  gchar ft[5] = {"a"};
  struct stat stbuf;
  olong next, temp;
  
  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (fileName!=NULL);
  errno = 0;  /* reset any system error */

  /* Save call arguments */
  in->access    = access;
  in->type      = type;
  in->blockSize = blockSize;
  if (in->fileName) g_free(in->fileName);
  in->fileName  = g_strdup(fileName);

  /* If currently open, flush (if needed) and  close first */
  if ((in->status==OBIT_Modified) || (in->status==OBIT_Active)) {
    retCode = ObitFileClose (in, err);
    if (err->error) /* add traceback on error */
      Obit_traceback_val (err, "ObitFileOpen", in->name, retCode);
  }

  in->status = OBIT_ErrorExist; /* in case something goes wrong */

  /* open file by access type */
  /*------------------------ Read Only --- --------------------------*/
  if ((access == OBIT_IO_ReadOnly) || (access == OBIT_IO_ReadCal) ) {
    in->exist = TRUE;  /* it better */
    ft[0] = 'r'; ft[1]=0; ft[2]=0;
    if (type==OBIT_IO_Binary) ft[1]='b';
    in->myFile = fopen (in->fileName, ft);
    if (in->myFile==NULL) {
      Obit_log_error(err, OBIT_Error,
		     "ERROR opening file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return OBIT_IO_OpenErr;
    }
    
    /*------------------------ Read/Write ---------------------------------*/
  } else if (access == OBIT_IO_ReadWrite) {
    /* does it exist ? */
    temp = stat(in->fileName, &stbuf);
    in->exist = (temp==0);
    ft[0] = 'r'; ft[1]='+'; ft[2]=0; ft[3]=0;
    if (!in->exist) ft[0] = 'w'; 
    if (type==OBIT_IO_Binary) ft[2]='b';
    in->myFile = fopen (in->fileName, ft);
    if (in->myFile==NULL) {
      Obit_log_error(err, OBIT_Error,
		     "ERROR opening file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return OBIT_IO_OpenErr;
    }

    /*------------------------ Write Only ---------------------------------*/
  } else if (access == OBIT_IO_WriteOnly) {
    /* does it exist ? */
    temp = stat(in->fileName, &stbuf);
    in->exist = (temp==0);
    /* Clear status as nonexistance of a file is considered an error */
    if (!in->exist) errno = 0;
    ft[0] = 'w'; ft[1]='+'; ft[2]=0; ft[3]=0;
    next = 2;
    if (type==OBIT_IO_Binary) ft[next++]='b';
    if (in->exist) ft[0]='r'; 
    in->myFile = fopen (in->fileName, ft);
    if (in->myFile==NULL) {
      Obit_log_error(err, OBIT_Error,
		     "ERROR opening file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return OBIT_IO_OpenErr;
    }
    /* unknown */      
  } else {
    /* should never get here */
    g_assert_not_reached(); 
  }

  /* position at the beginning of the file */
  in->filePos = 0;
  status = fseeko (in->myFile, in->filePos, SEEK_SET);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR Positioning file %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
    return OBIT_IO_OpenErr;
  }

  in->status = OBIT_Active;  /* seems to be OK */
  errno = 0; /* clear error code */
  return OBIT_IO_OK;
} /* end ObitFileOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitFileClose (ObitFile *in, ObitErr *err)
{
  ObitIOCode  status;
  gboolean   flush;

  /* error checks */
  /* attempt close even if error exists */
  g_assert (ObitIsA(in, &myClassInfo));

  /* don't bother if it's not open (shutdown even if in error) */
  if ((in->status!=OBIT_Modified) && (in->status!=OBIT_Active) 
      && (in->status!=OBIT_ErrorExist)) 
    return OBIT_IO_OK;

  /* if writing flush buffer */
  flush = (in->status==OBIT_Modified);
  if (flush) {
    status = fflush (in->myFile);
    if (status) {             /* it went wrong */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR Flushing file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return OBIT_IO_ReadErr;
    }
  }
  /* close file */
  status = fclose(in->myFile);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR Flushing file %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
    return OBIT_IO_ReadErr;
  }

  /* mark file structure */
  in->myFile = NULL;
  in->status = OBIT_Inactive;

  return OBIT_IO_OK;
} /* end ObitFileClose */

/** Position file
 * \param in      Pointer to object to be positioned.
 * \param filePos File position in bytes of beginning of file
 * \param err     ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitFileSetPos (ObitFile *in, ObitFilePos filePos, ObitErr *err)
{
  ObitIOCode status, retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitFileSetPos";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* position the file if needbe */
  if (in->filePos != filePos) {
    in->filePos = filePos;
    status = fseeko (in->myFile, in->filePos, SEEK_SET);
    if (status) {             /* it went wrong */
      Obit_log_error(err, OBIT_Error, 
		     "%s: ERROR Positioning file %s", routine, in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return OBIT_IO_OpenErr;
    }
  }

  return OBIT_IO_OK;
} /* end ObitFileSetPos */

/**
 * Position file at EOF
 * Note on return, in->filePos is invalid
 * \param in Pointer to object to be positioned.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitFileEnd (ObitFile *in, ObitErr *err)
{
  ObitIOCode status, retCode = OBIT_IO_SpecErr;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* position at the end of the file */
  in->filePos = 0;
  status = fseeko (in->myFile, in->filePos, SEEK_END);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR Positioning file %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
    return OBIT_IO_OpenErr;
  }

  /* reset file position - the SEEK_END confuses ftello */
  in->filePos = -1L;

  return OBIT_IO_OK;
} /* end ObitFileEnd */

/**
 * Read data from disk.
 * \param in      Pointer to object to be read.
 * \param filePos File position in bytes of beginning of read
 *                <0 => start at current position.
 * \param size    number of bytes to read.
 * \param buffer  pointer to buffer to accept results.
 * \param err     ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitFileRead (ObitFile *in, ObitFilePos filePos, olong size, gchar *buffer, 
	      ObitErr *err)
{
  ObitIOCode status, retCode = OBIT_IO_SpecErr;
  size_t nRead;
  gchar *routine = "ObitFileRead";
  /* olong tellpos; DEBUG  AIPS SUCKS BIG TIME */

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (buffer != NULL);
  /* A previous error condition? */
  if (in->status == OBIT_ErrorExist) return retCode;
  errno = 0;  /* reset any system error */

  /* If last operation was a write, flush */
  if (in->status == OBIT_Modified) {
    retCode = ObitFileFlush (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    in->status = OBIT_Active;
  }

  /* Position file if needbe */
  if ((filePos>=0) && (in->filePos != filePos)) {
    in->filePos = filePos;
    /* DEBUG 
       tellpos = ftell(in->myFile);
       fprintf (stdout, "%s read fseek %d was %d positioned %d \n",
       in->fileName, (olong)filePos, (olong)in->filePos, tellpos);
       END DEBUG */
    status = fseeko (in->myFile, in->filePos, SEEK_SET);
    if (status) {             /* it went wrong */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR Positioning file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return OBIT_IO_ReadErr;
    }
  } /* end of position file */

  /* Do read */
  /* DEBUG 
     tellpos = ftell(in->myFile);
     fprintf (stdout, "%s read %d at %d positioned  %d\n",
     in->fileName, size, (olong)in->filePos, tellpos);
     END DEBUG */
  retCode = OBIT_IO_OK;
  nRead = fread (buffer, size, 1, in->myFile);
  if (nRead!=1) {		/* it went wrong */
     /* set type of error */
    retCode = OBIT_IO_ReadErr;
    if (feof(in->myFile)) retCode = OBIT_IO_EOF;
    if (retCode!=OBIT_IO_EOF) { /* EOF not an error */
      in->status = OBIT_ErrorExist;
      Obit_log_error(err, OBIT_Error, 
		     "ERROR reading file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return retCode;
    }
  }

  /* update position */
  if (in->filePos>=0) in->filePos += size;
  
  return retCode;
} /* end ObitFileRead */

/**
 * Read line of text from file.
 * \param in      Pointer to object to be read.
 * \param line    pointer to memory to accept the text
 *                may be newline terminated.
 * \param lineMax size in characters of line
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileReadLine (ObitFile *in, gchar *line, olong lineMax, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *OK;
  gchar *routine = "ObitFileReadLine";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (line != NULL);
  /* A previous error condition? */
  if (in->status == OBIT_ErrorExist) return retCode;
  errno = 0;  /* reset any system error */

  /* This must be opened as a text file */
  if (in->type!=OBIT_IO_Text) {
      Obit_log_error(err, OBIT_Error, 
		     "%s ERROR NOT text file %s", routine, in->name);
      return retCode;
  }

  /* This must be allowed to read */
  if ((in->access!=OBIT_IO_ReadOnly) && (in->access!=OBIT_IO_ReadWrite)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s ERROR NO read access for %s", routine, in->name);
      return retCode;
  }

  /* If last operation was a write, flush */
  if (in->status == OBIT_Modified) {
    retCode = ObitFileFlush (in, err);
    if (err->error) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->status = OBIT_Active;
  }

  /* Do read */
  retCode = OBIT_IO_OK;
  OK = fgets (line, lineMax, in->myFile);
  if (!OK) {			/* it went wrong */
    /* set type of error */
    retCode = OBIT_IO_ReadErr;
    if (feof(in->myFile)) retCode = OBIT_IO_EOF;
    if (retCode!=OBIT_IO_EOF) { /* EOF not an error */
      in->status = OBIT_ErrorExist;
      Obit_log_error(err, OBIT_Error, "ERROR reading file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return retCode;
    }
  }

  /* update position */
  in->filePos = ftello(in->myFile);
  
  return retCode;
} /* end ObitFileReadLine */

/**
 * Write information to disk.
 * \param in     Pointer to object to be written.
 * \param filePos File position in bytes of beginning of write
 *                <0 => start at current position.
 * \param size    number of bytes to write.
 * \param buffer Pointer to buffer to receive output data.
 * \param err    ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitFileWrite (ObitFile *in, ObitFilePos filePos, olong size, gchar *buffer, 
	       ObitErr *err)
{
  ObitIOCode status, retCode = OBIT_IO_SpecErr;
  size_t nWrit;
  /*glong tellpos;  DEBUG  */
  gchar *routine = "ObitFileWrite";
  
  if (size<=0) return OBIT_IO_OK; /* Nothing to do? */

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (buffer != NULL);
  /* A previous error condition? */
  if (in->status== OBIT_ErrorExist) return retCode;
  errno = 0;  /* reset any system error */

  /* Trap request to write to end of file positioned at the end */
  if ((filePos<0) && (in->filePos<0)) {
    filePos = ObitFileSize(in->fileName, err);
    if (err->error) Obit_traceback_val (err, routine, in->name,  retCode);
  }
  
  /* Position file if needbe */
  if ((filePos>=0) && (in->filePos != filePos)) {
    /* DEBUG 
       tellpos = ftell(in->myFile);
       fprintf (stdout, "%s write fseek %d was %d positioned %d \n",
       in->fileName, (olong)filePos, (olong)in->filePos, tellpos);
       END DEBUG */
    in->filePos = filePos;
    status = fseeko (in->myFile, in->filePos, SEEK_SET);
    if (status) {             /* it went wrong */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR Positioning file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      return OBIT_IO_WriteErr;
    }
  } /* end of position file */
  
  /* Do write */
  nWrit = fwrite (buffer, size, 1, in->myFile);
  if (nWrit!=1) {		/* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR writing file %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/

    /* set type of error */
    retCode = OBIT_IO_WriteErr;
    if (feof(in->myFile)) retCode = OBIT_IO_EOF;
    in->status = OBIT_ErrorExist;
    return retCode;
  }

  /* DEBUG 
     status = fflush (in->myFile);
     tellpos = ftell(in->myFile);
     fprintf (stdout, "%s write %d at %d positioned  %d\n",
     in->fileName, size, (olong)in->filePos, tellpos);
     END DEBUG */

  /* update position / status */
  if (in->filePos>=0) in->filePos += size;
  in->status = OBIT_Modified;

  return OBIT_IO_OK;
} /* end ObitFileWrite */

/**
 * Write line of text to file.
 * \param in      Pointer to object to be written
 * \param line    pointer to memory to text to be written, 
 *                should be newline terminated.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileWriteLine (ObitFile *in, gchar *line, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong OK;
  gchar *routine = "ObitFileWriteLine";
  
  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (line != NULL);
  /* A previous error condition? */
  if (in->status== OBIT_ErrorExist) return retCode;
  errno = 0;  /* reset any system error */
 
 /* This must be opened as a text file */
  if (in->type!=OBIT_IO_Text) {
      Obit_log_error(err, OBIT_Error, 
		     "%s ERROR NOT text file %s", routine, in->name);
      return retCode;
  }

  /* This must be allowed to write */
  if ((in->access!=OBIT_IO_WriteOnly) && (in->access!=OBIT_IO_ReadWrite)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s ERROR NO write access for %s", routine, in->name);
      return retCode;
  }

  /* Do write */
  OK = fputs (line, in->myFile);
  if (OK<0) {			/* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR writing file %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
 
    /* set type of error */
    retCode = OBIT_IO_WriteErr;
    if (feof(in->myFile)) retCode = OBIT_IO_EOF;
    in->status = OBIT_ErrorExist;
    return retCode;
  }

  /* update position / status */
  in->filePos = ftello(in->myFile);
  in->status = OBIT_Modified;

  return OBIT_IO_OK;
} /* end ObitFileWriteLine */

/**
 * Fill the file from the current position to padTo using the 
 * contents of buffer.
 * This is mostly useful for AIPS files for which ancient
 * concepts of file blocking are enforced.
 * \param in     Pointer to object to be written.
 * \param padTo  Final (0-rel) file position to be written.
 * \param buffer buffer whose contents are to be used to pad the file.
 * \param bsize  Size in bytes of buffer.
 * \param err    ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode
ObitFilePad (ObitFile *in, olong padTo, gchar *buffer, olong bsize, 
	     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  size_t size, nWrit;
  /* olong tellpos; DEBUG  */

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (buffer != NULL);
  /* A previous error condition? */
  if (in->status== OBIT_ErrorExist) return retCode;
  errno = 0;  /* reset any system error */

  /* DEBUG
     tellpos = ftell(in->myFile);
     fprintf (stdout, "%s Pad to %d was %d positioned %d \n",
     in->fileName, (olong)padTo, (olong)in->filePos, tellpos);
     END DEBUG */

  /* Is anything needed? */
  if (padTo < in->filePos) return OBIT_IO_OK;

  /* pad it */
  size = bsize;
  while (in->filePos < padTo) { /* pad 'em */
    size = MIN (size, (padTo - in->filePos));

    /* write */
    nWrit = fwrite (buffer, size, 1, in->myFile);
    /* status = fflush (in->myFile); DEBUG */
    if (nWrit!=1) {		/* it went wrong */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR writing file %s", in->name);
      ObitFileErrMsg(err);     /* system error message*/
      
      /* set type of error */
      retCode = OBIT_IO_WriteErr;
      if (feof(in->myFile)) retCode = OBIT_IO_EOF;
      in->status = OBIT_ErrorExist;
      return retCode;
    }

    if (in->filePos>=0) in->filePos += size;  /* update file position */
  } /* end loop padding */

  return OBIT_IO_OK;
} /* end ObitFilePad */

/**
 * Assure that the total length of a file is an integral multiple
 * of blksize bytes.  Any unwritten bytes are zero filled.
 * This is mostly useful for AIPS files for which ancient
 * concepts of file blocking are enforced.
 * \param in      Pointer to object to be padded.
 * \param blksize Size in bytes of blocking size.
 * \param err     ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode
ObitFilePadFile (ObitFile *in, olong blksize, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *tmpBuff=NULL;
  size_t size, nBlocks, nWrit;
  ObitFilePos wantPos;
  struct stat stbuf;
  olong status, need;
  /* olong tellpos; DEBUG  */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (blksize > 0);
  /* A previous error condition? */
  if (in->status == OBIT_ErrorExist) return retCode;
  /* Better have write enabled */
  if (in->access==OBIT_IO_ReadOnly) {
    Obit_log_error(err, OBIT_Error, 
		   "Attempt to pad ReadOnly file for %s", in->name);
    return retCode;
  }

  /* DEBUG
     tellpos = ftell(in->myFile);
     fprintf (stdout, "%s PadFile was %d positioned %d \n",
     in->fileName, (olong)in->filePos, tellpos);
     END DEBUG */

  /* get file info */
  status = stat(in->fileName, &stbuf);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR stating file for %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
  }

  /* is padding needed? Total number of blocks needed */
  nBlocks = 1 + ((stbuf.st_size-1)/blksize);
  size = nBlocks * blksize;
  if (size <= stbuf.st_size) return OBIT_IO_OK;

  /* how much padding needed? */
  need = size - stbuf.st_size - 1;
  tmpBuff = g_malloc0(need);

  /* DEBUG
     tellpos = ftell(in->myFile);
     fprintf (stdout, "%s PadFile was %d positioned %d size  %d\n",
     in->fileName, (olong)in->filePos, tellpos, (olong)stbuf.st_size);
     END DEBUG */

  /* pad it */
  size = need;
  wantPos = stbuf.st_size+1;
  /* position */
  status = fseeko (in->myFile, wantPos, SEEK_SET);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR Positioning file %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
    if (tmpBuff) g_free(tmpBuff); tmpBuff = NULL;
    return OBIT_IO_OpenErr;
  }

  /* DEBUG 
     tellpos = ftell(in->myFile);
     fprintf (stdout, "%s PadFile middle %d positioned %d size %d  %d\n",
     in->fileName, (olong)wantPos, tellpos, (olong)size, (olong)stbuf.st_size);
     END DEBUG */
  /* write */
  nWrit = fwrite (tmpBuff, size, 1, in->myFile);
  if (nWrit!=1) {      /* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR padding file for %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
    
    /* set type of error */
    retCode = OBIT_IO_WriteErr;
    if (feof(in->myFile)) retCode = OBIT_IO_EOF;
    in->status = OBIT_ErrorExist;
    if (tmpBuff) g_free(tmpBuff); tmpBuff = NULL;
    return retCode;
  }
  
  /* DEBUG
     status = fflush (in->myFile);
     status = stat(in->fileName, &stbuf);
     tellpos = ftell(in->myFile);
     fprintf (stdout, "%s PadFile after positioned %d size  %d\n",
     in->fileName, tellpos, (olong)stbuf.st_size);
     END DEBUG */

  if (in->filePos>=0) in->filePos += size;  /* update file position */
  if (tmpBuff) g_free(tmpBuff); tmpBuff = NULL;

  return OBIT_IO_OK;
} /* end ObitFilePadFile */

/**
 * Force the transfer of the I/O buffer to disk.
 * Note: this is done automatically on Close.
 * \param in  Pointer to object to flush.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitFileFlush (ObitFile *in,  ObitErr *err)
{
  ObitIOCode status, retCode = OBIT_IO_SpecErr;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* A previous error condition? */
  if (in->status== OBIT_ErrorExist) return retCode;

  /* Do Flush */
  status = fflush (in->myFile);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR Flushing file %s", in->name);
    ObitFileErrMsg(err);     /* system error message*/
    return OBIT_IO_ReadErr;
  }

  in->status = OBIT_Active;  
  return OBIT_IO_OK;
} /* end ObitFileFlush */

/**
 * Adds error message to err if (errno.h) errno non zero
 * Note: set errno (in <errno.h> to zero before a system call 
 * that might set it.
 * errno is reset to 0.
 * \param err ObitErr for reporting errors.
 * \return True if error
 */
gboolean ObitFileErrMsg (ObitErr *err)
{
  gchar *errMsg;
  olong temp;

  /* error checks */
  temp = errno;
  if (temp==0) return FALSE;  /* Error? */

  errMsg = strerror(errno);

  Obit_log_error(err, OBIT_Error, errMsg); /* Set message */

  errno = 0;  /* reset errno */

  return TRUE;
} /* end ObitFileErrMsg */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFileClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFileClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFileClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFileClassInfoDefFn (gpointer inClass)
{
  ObitFileClassInfo *theClass = (ObitFileClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFileClassInit;
  theClass->newObit       = (newObitFP)newObitFile;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFileClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFileGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitFileCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFileClear;
  theClass->ObitInit      = (ObitInitFP)ObitFileInit;
  theClass->ObitFileZap   = (ObitFileZapFP)ObitFileZap;
  theClass->ObitFileRename= (ObitFileRenameFP)ObitFileRename;
  theClass->ObitCopy      = (ObitCopyFP)ObitFileCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFileClear;
  theClass->ObitInit      = (ObitInitFP)ObitFileInit;
  theClass->ObitFileExist = (ObitFileExistFP)ObitFileExist;
  theClass->ObitFileSize  = (ObitFileSizeFP)ObitFileSize;
  theClass->ObitFileName  = (ObitFileNameFP)ObitFileName;
  theClass->ObitFileOpen  = (ObitFileOpenFP)ObitFileOpen;
  theClass->ObitFileClose = (ObitFileCloseFP)ObitFileClose;
  theClass->ObitFileSetPos= (ObitFileSetPosFP)ObitFileSetPos;
  theClass->ObitFileEnd   = (ObitFileEndFP)ObitFileEnd;
  theClass->ObitFileRead  = (ObitFileReadFP)ObitFileRead;
  theClass->ObitFileReadLine  = 
    (ObitFileReadLineFP)ObitFileReadLine;
  theClass->ObitFileWrite = (ObitFileWriteFP)ObitFileWrite;
  theClass->ObitFileWriteLine = 
    (ObitFileWriteLineFP)ObitFileWriteLine;
  theClass->ObitFilePad   = (ObitFilePadFP)ObitFilePad;
  theClass->ObitFilePadFile = 
    (ObitFilePadFileFP)ObitFilePadFile;
  theClass->ObitFileFlush = (ObitFileFlushFP)ObitFileFlush;

} /* end ObitFileClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitFileInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFile *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->access    = OBIT_IO_None;
  in->type      = OBIT_IO_Binary;
  in->status    = OBIT_Inactive;
  in->exist     = FALSE;
  in->blockSize = 0;
  in->filePos   = 0;
  in->fileName  = NULL;
} /* end ObitFileInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitFileClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFile *in = inn;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->status==OBIT_Active) || (in->status==OBIT_Modified)) {
    err = newObitErr();
    ObitFileClose (in, err); 
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* free this class members */
  if (in->thread) in->thread  = ObitThreadUnref(in->thread);
  if (in->fileName) g_free(in->fileName); in->fileName = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitFileClear */

