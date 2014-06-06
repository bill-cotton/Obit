/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2014                                          */
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

#include "ObitIO.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIO.c
 * ObitIO class function definitions.
 *
 * This is a virtual base class and should never be directly instantiated.
 * Derived classes provide an I/O interface to various underlying disk
 * structures.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIO";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitIOClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitIO* newObitIO (gchar* name, ObitInfoList *info,
		       ObitErr *err)
{
  ObitIO* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitIO));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitIOInit((gpointer)out);

  return out;
} /* end newObitIO */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Check if underlying files are the same.
 * This test is done using values entered into the #ObitInfoList
 * in case the object has not yet been opened.
 * \param in1 ObitInfoList for first object to be tested
 * \param in2 ObitInfoList for second object to be tested
 * \param err ObitErr for reporting errors.
 * \return TRUE if to objects have the same underlying structures
 * else FALSE
 */
gboolean ObitIOSame (ObitIO *in, ObitInfoList *in1, ObitInfoList *in2, ObitErr *err)
{
  const ObitIOClassInfo *myClass;
  ObitIOType FileType1,FileType2 ;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOSame";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return same;

  /* Check if the underlying types are the same */
  if (!ObitInfoListGet(in1, "FileType", &type, dim, 
		       (gpointer)&FileType1, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }
  if (!ObitInfoListGet(in2, "FileType", &type, dim, 
		       (gpointer)&FileType2, err)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }
  if (err->error) return FALSE;

  /* Compare */
  if (FileType1 != FileType2) return FALSE;

  /* call IO function */
  myClass = in->ClassInfo;

  /* Don't call yourself */
  Obit_retval_if_fail((myClass!=(const ObitIOClassInfo*)&myClassInfo), 
		      err, same,
		      "%s: recursive call for %s", routine, in->name);
  same = myClass->ObitIOSame (in, in1, in2, err);

  return same;
} /* end ObitIOSame */

/**
 * Delete underlying files.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOZap (ObitIO *in, ObitErr *err)
{
  const ObitIOClassInfo *myClass;

  /* Don't bother if NULL */
  if (!in) return;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;

  /* call actual function */
  myClass->ObitIOZap (in, err);

  return;
} /* end ObitIOZap */

/**
 * Rename underlying files.
 * New name information depends on the underlying file type and is
 * given on the info member.
 * \param in Pointer to object to be zapped.
 * \param info Associated ObitInfoList
 * For FITS files:
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 *
 * For AIPS:
 * \li "newName" OBIT_string (12,1,1) New AIPS Name 
 *      absent or Blank = don't change
 * \li "newClass" OBIT_string (6,1,1) New AIPS Class
 *      absent or Blank = don't changeO
 * \li "newSeq" OBIT_int (1,1,1) New AIPS Sequence
 *      0 => unique value
 * \param err ObitErr for reporting errors.
 */
void ObitIORename (ObitIO *in, ObitInfoList *info, ObitErr *err)
{
  const ObitIOClassInfo *myClass;

  /* Don't bother if NULL */
  if (!in) return;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;

  /* call actual function */
  myClass->ObitIORename (in, info, err);

  return;
} /* end ObitIORename */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIO* ObitIOCopy  (ObitIO *in, ObitIO *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitIO(outName, NULL, err);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->access = OBIT_IO_None; /* not currently reading file */

  return out;
} /* end ObitIOCopy */

/**
 * Initialize structures and open file.
 * The file and selection info member should have been stored in the ObitInfoList
 * prior to calling.  See derived classes for details.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOpen (ObitIO *in, ObitIOAccess access, ObitInfoList *info, 
	     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOOpen != NULL);

  /* call actual function */
  retCode = myClass->ObitIOOpen (in, access, info, err);
  in->access = access; /* just in case */

  return retCode;
} /* end ObitIOOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOClose (ObitIO *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  if (in==NULL) return OBIT_IO_OK;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOClose != NULL);

  /* call actual function */
  retCode = myClass->ObitIOClose (in, err);

  return retCode;
} /* end ObitIOClose */

/**
 * initialize I/O
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOSet (ObitIO *in, ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOSet != NULL);

  /* call actual function */
  retCode = myClass->ObitIOSet (in, info, err);

  return retCode;
} /* end ObitIOInit */

/**
 * Read data from disk.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIORead (ObitIO *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIORead != NULL);

  /* call actual function */
  retCode = myClass->ObitIORead (in, data, err);

  return retCode;
} /* end ObitIORead */

/**
 * Read data from disk to multiple buffers .
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOReadMulti (olong nBuff, ObitIO **in, ofloat **data, 
			    ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  if (err->error) return retCode;

  /* this is a virtual function, see if actual one defined */
  myClass = in[0]->ClassInfo;
  g_assert (myClass->ObitIOReadMulti != NULL);

  /* call actual function */
  retCode = myClass->ObitIOReadMulti (nBuff, in, data, err);

  return retCode;
} /* end ObitIOReadMulti */

/**
 * Reread data from disk to multiple buffers.
 * Retreives data read in a previous call to ObitIOReadMulti
 * NOTE: this depends on retreiving the data from the first element in in
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOReReadMulti (olong nBuff, ObitIO **in, ofloat **data, 
			    ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  if (err->error) return retCode;

  /* this is a virtual function, see if actual one defined */
  myClass = in[0]->ClassInfo;
  g_assert (myClass->ObitIOReadMulti != NULL);

  /* call actual function */
  retCode = myClass->ObitIOReReadMulti (nBuff, in, data, err);

  return retCode;
} /* end ObitIOReReadMulti */

/**
 * Read data from disk specifying starting row.
 * \param in Pointer to object to be read.
 * \param rowno Starting row number (1-rel) -1=> next.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOReadRow (ObitIO *in, olong rowno, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOReadRow != NULL);

  /* call actual function */
  retCode = myClass->ObitIOReadRow (in, rowno, data, err);

  return retCode;
} /* end ObitIOReadRow */

/**
 * Read data from disk applying selection.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOReadSelect (ObitIO *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOReadSelect != NULL);

  /* call actual function */
  retCode = myClass->ObitIOReadSelect (in, data, err);

  return retCode;
} /* end ObitIOReadSelect */

/**
 * Read data from disk applying selection to multiple buffers .
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOReadMultiSelect (olong nBuff, ObitIO **in, ofloat **data, 
				  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  if (err->error) return retCode;

  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in[0]->ClassInfo;
  g_assert (myClass->ObitIOReadMultiSelect != NULL);

  /* call actual function */
  retCode = myClass->ObitIOReadMultiSelect (nBuff, in, data, err);

  return retCode;
} /* end ObitIOReadMultiSelect */

/**
 * Reread data from disk applying selection to multiple buffers .
 * Retreives data read in a previous call to ObitIOReadMultiSelect
 * possibly applying new calibration.
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOReReadMultiSelect (olong nBuff, ObitIO **in, ofloat **data, 
				    ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  if (err->error) return retCode;
  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in[0]->ClassInfo;
  g_assert (myClass->ObitIOReReadMultiSelect != NULL);

  /* call actual function */
  retCode = myClass->ObitIOReReadMultiSelect (nBuff, in, data, err);

  return retCode;
} /* end ObitIOReReadMultiSelect */

/**
 * Read data from disk specifying start row and applying selection.
 * \param in Pointer to object to be read.
 * \param rowno Starting row number (1-rel) -1=> next.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOReadRowSelect (ObitIO *in, olong rowno, ofloat *data, 
			     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOReadRowSelect != NULL);

  /* call actual function */
  retCode = myClass->ObitIOReadRowSelect (in, rowno, data, err);

  return retCode;
} /* end ObitIOReadRowSelect */

/**
 * Write information to disk.
 * Writes row in->myDesc->row + 1; plane in->myDesc->plane + 1
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOWrite (ObitIO *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOWrite != NULL);

  /* call actual function */
  retCode = myClass->ObitIOWrite (in, data, err);

  return retCode;
} /* end ObitIOWrite */

/**
 * Write information to disk specifying start row.
 * Writes row in->myDesc->row + 1; plane in->myDesc->plane + 1
 * \param in Pointer to object to be written.
 * \param rowno Starting row number (1-rel) -1=> next.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOWriteRow (ObitIO *in, olong rowno, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOWriteRow != NULL);

  /* call actual function */
  retCode = myClass->ObitIOWriteRow (in, rowno, data, err);

  return retCode;
} /* end ObitIOWriteRow */

/**
 * Read image Descriptor data from disk.
 * \param in Pointer to object  with ObitImageDescto be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOReadDescriptor (ObitIO *in,  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;
  gchar *routine = "ObitIOReadDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  Obit_retval_if_fail((in->myStatus!=OBIT_Inactive), err, retCode,
		      "%s: IO inactive for %s", routine, in->name);
  

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOReadDescriptor !=NULL);

  /* call actual function */
  retCode = myClass->ObitIOReadDescriptor (in, err);

  return retCode;
} /* end ObitIOReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in Pointer to object with ObitImageDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOWriteDescriptor (ObitIO *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;
  gchar *routine = "ObitIOWriteDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  Obit_retval_if_fail((in->myStatus!=OBIT_Inactive), err, retCode,
		      "%s: IO inactive for %s", routine, in->name);
  
  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOWriteDescriptor != NULL);

  /* call actual function */
  retCode = myClass->ObitIOWriteDescriptor (in, err);

  return retCode;
} /* end ObitIOWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOFlush (ObitIO *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOFlush != NULL);

  /* call actual function */
  retCode = myClass->ObitIOFlush (in, err);

  return retCode;
} /* end ObitIOFlush */

/**
 * Create buffer approptiate for I/O request
 * \param data (output) pointer to data array
 * \param size (output) size of data array in floats.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 */
void ObitIOCreateBuffer (ofloat **data, ollong *size, ObitIO *in, 
			    ObitInfoList *info, ObitErr *err)
{
  const ObitIOClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  g_assert (data != NULL);
  g_assert (size != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOCreateBuffer != NULL);

  /* call actual function */
  myClass->ObitIOCreateBuffer (data, size, in, info, err);

  return;
} /* end ObitIOCreateBuffer */

/**
 * Destroy buffer
 * \param buffer Pointer to buffer to destroy.
 */
void ObitIOFreeBuffer (ofloat *buffer)
{
  /* error checks */
  if (buffer==NULL) return;

  ObitMemFree (buffer);

} /* end ObitIOFreeBuffer */

/**
 * Return a ObitTable Object to a specified table associated with
 * the input ObitIO.  
 * If such an object exists, a reference to it is returned,
 * else a new object is created and entered in the ObitTableList.
 * Returned object is typed an Obit to prevent circular definitions
 * between the ObitTable and the ObitIO classes.
 * \param in       Pointer to object with associated tables.
 *                 This MUST have been opened before this call.
 * \param access   access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *                 or OBIT_IO_WriteOnly).
 *                 This is used to determine defaulted version number
 *                 and a different value may be used for the actual 
 *                 Open.
 * \param tabType  The table type (e.g. "AIPS CC").
 * \param tabVer   Desired version number, may be zero in which case
 *                 the highest extant version is returned for read
 *                 and the highest+1 for write.
 * \param err      ObitErr for reporting errors.
 * \return pointer to created ObitTable, NULL on failure.
 */
Obit* 
newObitIOTable (ObitIO *in, ObitIOAccess access, 
		   gchar *tabType, olong *tabVer, ObitErr *err)
{
  Obit *out;
  const ObitIOClassInfo *myClass;
  gchar *routine = "newObitIOTable";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  g_assert(tabType!=NULL);
  g_assert(tabVer!=NULL);

  /* the Tablelist object must be present */
  if (in->tableList==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "my tableList member is NULL, open %s first", 
		     in->name);
      return NULL;
  }

  /* details depend on underlying file type,
     pass it down to ObitIO to call relevant routine */
  myClass = in->ClassInfo;
  g_assert (myClass->newObitIOTable != NULL);

  /* call actual function */
  out = myClass->newObitIOTable (in, access, tabType, tabVer, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  return out;
} /* end newObitIOTable */

/**
 * Update any disk resident structures about the current tables.
 * \param in   Pointer to object to be updated.
 * \param info ObitInfoList of parent object.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOUpdateTables (ObitIO *in, ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOClassInfo *myClass;
  gboolean openClose;
  gchar *routine="ObitIOUpdateTables";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;

  /* If I/O not defined - ignore */
  if (!ObitIsA(in, &myClassInfo)) return OBIT_IO_OK;

  /* Need to open and close? */
  openClose = !((in->myStatus==OBIT_Active) || (in->myStatus==OBIT_Modified));

  /* Open if needed */
  if (openClose) {
    retCode = ObitIOOpen (in, OBIT_IO_ReadWrite, info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);

    /* read descriptor */
    retCode = ObitIOReadDescriptor (in, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOUpdateTables != NULL);

  /* call actual function */
  retCode = myClass->ObitIOUpdateTables (in, info, err);

  /* Close if needed */
  if (openClose) {
    retCode = ObitIOClose (in, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  }
 
  return retCode;
} /* end ObitImageUpdateTables */

/**
 * Get underlying file information in entries to an ObitInfoList
 * Virtual function - should always be a call to an explicit class instance
 * Following entries for AIPS files ("xxx" = prefix):
 * \param in      Object of interest.
 * \param myInfo  InfoList on basic object with selection
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 * \param outList InfoList to write entries into
 * Following entries for AIPS files ("xxx" = prefix):
 * \li xxxName  OBIT_string  AIPS file name
 * \li xxxClass OBIT_string  AIPS file class
 * \li xxxDisk  OBIT_oint    AIPS file disk number
 * \li xxxSeq   OBIT_oint    AIPS file Sequence number
 * \li AIPSuser OBIT_oint    AIPS User number
 * \li xxxCNO   OBIT_oint    AIPS Catalog slot number
 * \li xxxDir   OBIT_string  Directory name for xxxDisk
 *
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFileName OBIT_string  FITS file name
 * \li xxxDisk     OBIT_oint    FITS file disk number
 * \li xxxDir      OBIT_string  Directory name for xxxDisk
 *
 * For all File types:
 * \li xxxFileType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxDataType OBIT_string "AIPS", "FITS"
 *    
 * For xxxDataType = "Table"
 * \li xxxTab   OBIT_string  (Tables only) Table type (e.g. "AIPS CC")
 * \li xxxVer   OBIT_oint    (Tables Only) Table version number
 *    
 * For xxxDataType = "MA"
 * \li xxxBLC   OBIT_oint[7] (Images only) 1-rel bottom-left corner pixel
 * \li xxxTRC   OBIT_oint[7] (Images Only) 1-rel top-right corner pixel
 *    
 * For xxxDataType = "OTF"
 * \li xxxnRecPIO OBIT_int (1,1,1) Number of vis. records per IO call
 *
 * For xxxDataType = "UV"
 * \li xxxnVisPIO  OBIT_int (1,1,1) Number of vis. records per IO call
 * \li xxxdoCalSelect OBIT_bool (1,1,1) Select/calibrate/edit data?
 * \li xxxStokes OBIT_string (4,1,1) Selected output Stokes parameters:
 *              "    "=> no translation,"I   ","V   ","Q   ", "U   ", 
 *              "IQU ", "IQUV",  "IV  ", "RR  ", "LL  ", "RL  ", "LR  ", 
 *              "HALF" = RR,LL, "FULL"=RR,LL,RL,LR. [default "    "]
 *              In the above 'F' can substitute for "formal" 'I' (both RR+LL).
 * \li xxxBChan OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li xxxEChan OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li xxxBIF   OBIT_int (1,1,1) First "IF" selected. [def all]
 * \li xxxEIF   OBIT_int (1,1,1) Highest "IF" selected. [def all]
 * \li xxxdoPol OBIT_int (1,1,1) >0 -> calibrate polarization.
 * \li xxxdoCalib OBIT_int (1,1,1) >0 -> calibrate, 2=> also calibrate Weights
 * \li xxxgainUse OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \li xxxflagVer OBIT_int (1,1,1) Flag table version, 0-> use highest, <0-> none
 * \li xxxBLVer   OBIT_int (1,1,1) BL table version, 0> use highest, <0-> none
 * \li xxxBPVer   OBIT_int (1,1,1) Band pass (BP) table version, 0-> use highest
 * \li xxxSubarray OBIT_int (1,1,1) Selected subarray, <=0->all [default all]
 * \li xxxdropSubA OBIT_bool (1,1,1) Drop subarray info?
 * \li xxxFreqID   OBIT_int (1,1,1) Selected Frequency ID, <=0->all [default all]
 * \li xxxtimeRange OBIT_float (2,1,1) Selected timerange in days.
 * \li xxxUVRange  OBIT_float (2,1,1) Selected UV range in kilowavelengths.
 * \li xxxInputAvgTime OBIT_float (1,1,1) Input data averaging time (sec).
 *               used for fringe rate decorrelation correction.
 * \li xxxSources OBIT_string (?,?,1) Source names selected unless any starts with
 *               a '-' in which case all are deselected (with '-' stripped).
 * \li xxxsouCode OBIT_string (4,1,1) Source Cal code desired, '    ' => any code selected
 *                                  '*   ' => any non blank code (calibrators only)
 *                                  '-CAL' => blank codes only (no calibrators)
 * \li xxxQual    Obit_int (1,1,1)  Source qualifier, -1 [default] = any
 * \li xxxAntennas OBIT_int (?,1,1) a list of selected antenna numbers, if any is negative
 *                 then the absolute values are used and the specified antennas are deselected.
 * \li xxxcorrType OBIT_int (1,1,1) Correlation type, 0=cross corr only, 1=both, 2=auto only.
 * \li xxxpassAl l OBIT_bool (1,1,1) If True, pass along all data when selecting/calibration
 *                                 even if it's all flagged, 
 *                                 data deselected by time, source, antenna etc. is not passed.
 * \li xxxdoBand  OBIT_int (1,1,1) Band pass application type <0-> none
 *     (1) if = 1 then all the bandpass data for each antenna
 *         will be averaged to form a composite bandpass
 *         spectrum, this will then be used to correct the data.
 *     (2) if = 2 the bandpass spectra nearest in time (in a weighted
 *         sense) to the uv data point will be used to correct the data.
 *     (3) if = 3 the bandpass data will be interpolated in time using
 *         the solution weights to form a composite bandpass spectrum,
 *         this interpolated spectrum will then be used to correct the
 *         data.
 *     (4) if = 4 the bandpass spectra nearest in time (neglecting
 *         weights) to the uv data point will be used to correct the
 *         data.
 *     (5) if = 5 the bandpass data will be interpolated in time ignoring
 *         weights to form a composite bandpass spectrum, this
 *         interpolated spectrum will then be used to correct the data.
 * \li xxxSmooth  OBIT_float (3,1,1) specifies the type of spectral smoothing
 *        Smooth(1) = type of smoothing to apply:
 *           0 => no smoothing
 *           1 => Hanning
 *           2 => Gaussian
 *           3 => Boxcar
 *           4 => Sinc (i.e. sin(x)/x)
 *         Smooth(2) = the "diameter" of the function, i.e.
 *           width between first nulls of Hanning triangle
 *           and sinc function, FWHM of Gaussian, width of
 *           Boxcar. Defaults (if < 0.1) are 4, 2, 2 and 3
 *           channels for Smooth(1) = 1 - 4.
 *         Smooth(3) = the diameter over which the convolving
 *           function has value - in channels.
 *           Defaults: 1, 3, 1, 4 times Smooth(2) used when
 * \li xxxSubScanTime Obit_float scalar [Optional] if given, this is the 
 *          desired time (days) of a sub scan.  This is used by the 
 *          selector to suggest a value close to this which will
 *          evenly divide the current scan.  See #ObitUVSelSubScan
 *          0 => Use scan average.
 *          This is only useful for ReadSelect operations on indexed ObitUVs.
 * \param err     ObitErr for reporting errors.
 */
void ObitIOGetFileInfo (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
			ObitInfoList *outList, ObitErr *err)
{
  gchar *routine = "ObitIOGetFileInfo";

  Obit_log_error(err, OBIT_Error, 
		 "%s: Virtual Function called for %s", routine, in->name);
} /* end ObitIOGetFileInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOClassInfoDefFn (gpointer inClass)
{
  ObitIOClassInfo *theClass = (ObitIOClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOGetClass;
  theClass->newObit       = NULL;
  theClass->newObitIO     = (newObitIOFP)newObitIO;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOCopy;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOSame;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOZap;
  theClass->ObitIORename  = (ObitIORenameFP)ObitIORename;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitIOClear;
  theClass->ObitInit      = (ObitInitFP)ObitIOInit;
  theClass->ObitIOOpen    = (ObitIOOpenFP)ObitIOOpen;
  theClass->ObitIOClose   = (ObitIOCloseFP)ObitIOClose;
  theClass->ObitIOSet     = (ObitIOSetFP)ObitIOSet;
  theClass->ObitIORead    = (ObitIOReadFP)ObitIORead;
  theClass->ObitIOReadMulti = (ObitIOReadMultiFP)ObitIOReadMulti;
  theClass->ObitIOReReadMulti = (ObitIOReReadMultiFP)ObitIOReReadMulti;
  theClass->ObitIOReadRow = (ObitIOReadRowFP)ObitIOReadRow;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOReadSelect;
  theClass->ObitIOReadMultiSelect = 
    (ObitIOReadMultiSelectFP)ObitIOReadMultiSelect;
  theClass->ObitIOReReadMultiSelect = 
    (ObitIOReReadMultiSelectFP)ObitIOReReadMultiSelect;
  theClass->ObitIOReadRowSelect = 
    (ObitIOReadRowSelectFP)ObitIOReadRowSelect;
  theClass->ObitIOWriteRow= (ObitIOWriteRowFP)ObitIOWriteRow;
  theClass->ObitIOWrite   = (ObitIOWriteFP)ObitIOWrite;
  theClass->ObitIOFlush   = (ObitIOFlushFP)ObitIOFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOCreateBuffer;
  theClass->ObitIOFreeBuffer   = 
    (ObitIOFreeBufferFP)ObitIOFreeBuffer;
  theClass->newObitIOTable   = 
    (newObitIOTableFP)newObitIOTable;
  theClass->ObitIOUpdateTables   = 
    (ObitIOUpdateTablesFP)ObitIOUpdateTables;
} /* end ObitIOClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIO *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->access   = OBIT_IO_None;
  in->myStatus = OBIT_Inactive;
  in->myDesc   = NULL;
  in->mySel    = NULL;
  in->myCal    = NULL;
  in->tableList= NULL;
} /* end ObitIOInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIO *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->myCal)  in->myCal  = ObitUVCalUnref((ObitUVCal*)in->myCal);
  if (in->myDesc) in->myDesc = ObitUVDescUnref(in->myDesc);
  if (in->mySel)  in->mySel  = ObitUVSelUnref(in->mySel);
  if (in->tableList) in->tableList = ObitUnref(in->tableList);
  
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOClear */


