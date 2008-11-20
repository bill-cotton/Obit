/* $Id$   */
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

#include "ObitIOOTFFITS.h"
#include "ObitOTFCal.h"
#include "ObitTableList.h"
#include "ObitFile.h"
#include "ObitFITS.h"
#include "ObitMem.h"
#include "ObitOTFArrayGeom.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOOTFFITS.c
 * ObitIOOTFFITS class function definitions.
 * This class is derived from the ObitIO class.
 *
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOOTFFITS";

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOOTFFITSClassInfo myClassInfo = {FALSE};

/** Function to obtain parent ClassInfo - ObitIO */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/*----------------- Union definitions ----------------------*/
/** Used for byte swapping shorts */
 union sequiv { 
   gshort full; 
   gchar parts[2];
 }; 

/** Used for byte swapping floats */
 union fequiv { 
   ofloat full; 
   gchar parts[4]; 
 }; 

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOOTFFITSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOOTFFITSClear (gpointer in);

/** Private: Determine next visibility to read */
static gboolean ObitIOOTFFITSNext (ObitIOOTFFITS *in, ObitErr *err);

/** Private: Copy Floats with byte swap to FITS order */
static void ObitIOOTFFITSfH2F (olong n, ofloat *in, ofloat *out);

/** Private: Copy Floats with byte swap to host order */
static void ObitIOOTFFITSfF2H (olong n, ofloat *in, ofloat *out);

/** Private: Read Sort Order. */
static void  ObitIOOTFFITSSortRead(ObitIOOTFFITS *in, olong *lstatus);

/** Private: Write  Sort Order.*/
static void  ObitIOOTFFITSSortWrite (ObitIOOTFFITS *in, olong *lstatus);

/** Private: Copy other header keywords. */
void  ObitIOOTFKeysOtherRead(ObitIOOTFFITS *in, olong *lstatus, 
			     ObitErr *err);

/** Private: Fix bug in cfitsio keyword parsing */
void ObitIOOTFFITSFixBug(gchar *out, gchar *in, olong maxn);

/** Private: Set Class function pointers. */
static void ObitIOOTFFITSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class on the first call.
 * \param name An optional name for the object.
 * \param info if non-NULL it is used to initialize the new object.
 * \param err  ObitErr for error messages.
 * \return the new object.
 */
ObitIOOTFFITS* newObitIOOTFFITS (gchar *name, ObitInfoList *info,
			       ObitErr *err)
{
  ObitIOOTFFITS* out;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  gchar tempStr[201];
  gchar *routine = "newObitIOOTFFITS";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitIOOTFFITSClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitIOOTFFITS));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOOTFFITSInit((gpointer)out);

  /* Get any info from info input */
  if (info!=NULL) {
    if(!ObitInfoListGet(info, "Disk", &type, (gint32*)dim, 
			&out->disk, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "FileName", &type, (gint32*)dim, 
			tempStr, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    /* form file name for file */
    /* fetch file name from temporary buffer, null terminate. */
    tempStr[dim[0]] = 0;
    if (out->FileName) g_free(out->FileName); /* release old */
    out->FileName = ObitFITSFilename (out->disk, tempStr, err);
    if (err->error) Obit_traceback_val (err, routine, name, out);
  }

  return out;
} /* end newObitIOOTFFITS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOOTFFITSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOOTFFITSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIOOTFFITSGetClass */

/**
 * Check if underlying files are the same.
 * This test is done using values entered into the #ObitInfoList
 * in case the object has not yet been opened.
 * \param in  ObitIO for test
 * \param in1 ObitInfoList for first object to be tested
 * \param in2 ObitInfoList for second object to be tested
 * \param err ObitErr for reporting errors.
 * \return TRUE if to objects have the same underlying structures
 * else FALSE
 */
gboolean ObitIOOTFFITSSame (ObitIO *in, ObitInfoList *in1, 
			   ObitInfoList *in2, ObitErr *err)
{
  olong disk1, disk2;
  gchar *filename1, *filename2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOOTFFITSSame";

  /* error checks */
  if (err->error) return same;

  /* get file from info */
  if(!ObitInfoListGet(in1, "Disk", &type, dim, &disk1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in1, "FileName", &type, (gint32*)&dim, 
		       (gpointer)&filename1)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }

  if(!ObitInfoListGet(in2, "Disk", &type, dim, &disk2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in2, "FileName", &type, (gint32*)&dim, 
		       (gpointer)&filename2)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }

  /* Compare */
  same = (disk1==disk2) && 
    !strncmp (filename1,filename2, 200);
 
  return same;
} /* end ObitIOOTFFITSSame */

/**
 * Rename underlying files.
 * New name information is given on the info member:
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 * \param in   Pointer to object to be renamed
 * \param info Associated ObitInfoList
 * \param err   ObitErr for reporting errors.
 */
void ObitIOOTFFITSRename (ObitIO *in, ObitInfoList *info, 
			  ObitErr *err)
{
  gchar *routine = "ObitIOOTFFITSRename";

  /* Don't bother if NULL */
  if (!in) return;

  /* error checks */
  if (err->error) return;

  /* Rename */
  ObitFITSRename (in, info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitIOOTFFITSRename */
/**
 * Delete underlying files.
 * Delete the whole FITS file.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOOTFFITSZap (ObitIOOTFFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar tempStr[201];
  gchar *routine = "ObitIOOTFFITSZap";

   /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Close if still open */
  if ((in->myStatus==OBIT_Modified) || (in->myStatus==OBIT_Active)) {
    retCode = ObitIOOTFFITSClose (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Destroy using cfitsio */
  if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
  else strncpy (tempStr, in->FileName, 200);
  fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status);
  fits_delete_file(in->myFptr, &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR %d deleting FITS file %s", 
		   status, in->FileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
  }

  /* FITS tables are in the same file - delete table list */
  in->tableList = ObitTableListUnref(in->tableList);

 return;
} /* end ObitIOOTFFITSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOOTFFITS* ObitIOOTFFITSCopy  (ObitIOOTFFITS *in, 
				       ObitIOOTFFITS *out, ObitErr *err)
{
  const ObitIOClassInfo *ParentClass;
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
    out = newObitIOOTFFITS(outName, NULL, err);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy this class */
  if (out->FileName!=NULL) g_free(out->FileName);
  out->FileName = g_strdup(in->FileName);

  return out;
} /* end ObitIOOTFFITSCopy */

/**
 * Initialize structures and open file.
 * The file etc. info should have been stored in the ObitInfoList.
 * The descriptor is read if ReadOnly or ReadOnly and
 * written to disk if opened WriteOnly.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk"     OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) FITS file name.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSOpen (ObitIOOTFFITS *in, ObitIOAccess access, 
			     ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  gchar tempStr[201];
  ObitOTFDesc* desc;
  gchar *routine = "ObitIOOTFFITSOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  g_assert (in->myDesc != NULL);
  g_assert (in->mySel != NULL);

  desc = in->myDesc;
  if(!ObitInfoListGet(info, "Disk", &type, (gint32*)dim, 
		      &in->disk, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
 
  if(!ObitInfoListGet(info, "FileName", &type, (gint32*)dim, 
		      tempStr, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* form file name for file */
  /* fetch file name from temporary buffer, null terminate. */
  tempStr[dim[0]] = 0;
  if (in->FileName) g_free(in->FileName); /* release old */
  in->FileName = ObitFITSFilename (in->disk, tempStr, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  /* open file by access type */
  /*------------------------ Read Only ---------------------------------*/
  if (access == OBIT_IO_ReadOnly) {
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);

    /* cfitsio refuses to open a file readwrite after it has been opened
       readonly so first try opening readwrite even if requested ReadOnly.
       If that fails try readonly. */
    /* Test open readwrite */
    fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status);
    if ((status==FILE_NOT_OPENED) || (status==READONLY_FILE))  { 
      /* Failed - try readonly */
      status = 0;
      fits_clear_errmsg();   /* Clear cfitsio error stack */
      if (fits_open_file(&(in->myFptr), (char*)tempStr, READONLY, &status) ) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR %d opening input FITS file %s", status, 
		       in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
    }

  /*------------------------ Read/Write ---------------------------------*/
  } else if ((access == OBIT_IO_ReadWrite) || (access == OBIT_IO_ReadCal)) {
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);
    /* Initialize output file */
    if ( fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status) ) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening output FITS file %s", 
		     in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      retCode = OBIT_IO_OpenErr;
      return retCode;
    }

  /*------------------------ Write Only ---------------------------------*/
  } else if (access == OBIT_IO_WriteOnly) {
    /* Initialize output file */
    /* Output file may already exist test open */
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);

    /* Open read/write to see if it's there */
    fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status);

    if (status==FILE_NOT_OPENED) { /* Not there - initialize output file */
      status = 0;
      fits_clear_errmsg();   /* Clear error stack */
      if (fits_create_file(&(in->myFptr), (char*)in->FileName, &status) ) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS file %s", 
		       in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
    } else if (status!=0) { /* error */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening output FITS file %s", 
		     in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      retCode = OBIT_IO_OpenErr;
    return retCode;
    }
  } else {
    /* should never get here */
    g_assert_not_reached(); 
  }

  /* Input table positioned in ReadDescriptor function */
  /* Output table created in WriteDescriptor function  */

  /* save information */
  in->access = access;
  in->myStatus = OBIT_Active;

  /* initialize location in data */
  desc->firstRec   = 0;
  desc->numRecBuff = 0;
  
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOOTFFITSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSClose (ObitIOOTFFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* don't bother if it's not open */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) 
    return OBIT_IO_OK;

  fits_close_file (in->myFptr, &status);
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR closing FITS file");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    retCode = OBIT_IO_CloseErr;
    return retCode;
  }

 in->myStatus = OBIT_Inactive;
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOOTFFITSClose */

/**
 * initialize I/O, position to beginning of uvdata.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSSet (ObitIOOTFFITS *in, ObitInfoList *info, 
			       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  int status = 0;

  /* Position to "OTFScanData" table version 1 */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, "OTFScanData", 1, &status);
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR positioning FITS file");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    if ((in->access==OBIT_IO_ReadOnly) ||
	(in->access==OBIT_IO_ReadWrite)) retCode = OBIT_IO_ReadErr;
    else  retCode = OBIT_IO_WriteErr;
    return retCode;
  }

  return retCode;
} /* end ObitIOOTFFITSSet */

/**
 * Read data from disk.
 * The number read will be mySel->nRecPIO (until the end of the selected
 * range of records in which case it will be smaller).
 * The first row number after a read is myDesc->firstRec
 * and the number of rows attempted is mySel->numRecRead; 
 * actual value saved as myDesc->numRecBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSRead (ObitIOOTFFITS *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitOTFDesc* desc;
  ObitOTFSel* sel;
  long size;
  int status = 0;
  olong len, i, ip, op, numRecBuff;
  gboolean done;
  ofloat *IOBuff = data;
  gchar *routine = "ObitIOOTFFITSRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  desc = in->myDesc; /* OTF descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  desc->numRecBuff =  0; /* No data yet */

  /* what next ? */
  done = ObitIOOTFFITSNext (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* check if done - all records read */
  if (done) {
    ObitIOOTFFITSClose (in, err); /* Close */
    return OBIT_IO_EOF;
  }

  len = desc->lrec; /* How big is a visibility */

  /* read block of sel->numRecRead visibilities at a time  */
  /* transfer size in bytes */
  size = sel->numRecRead * len * sizeof(ofloat); 

  /* Read */
  fits_read_tblbytes (in->myFptr, desc->firstRec, 1, size, (guchar*)IOBuff, 
		      &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR reading FITS OTF data for %s", in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* Byte swap */
  ip = op = 0; /* array pointers */
  numRecBuff = 0; /* How many valid visibilities */
  for (i=0; i<sel->numRecRead; i++) {
    /* Copy/byteswap descriptive parameters */
    ObitIOOTFFITSfF2H (desc->numDesc, &IOBuff[ip], &data[op]);
    /* Copy/byteswap data array */
    ObitIOOTFFITSfF2H (desc->colRepeat[desc->ilocdata], &IOBuff[ip+desc->ilocdata],
		      &data[op+desc->ilocdata]);
    ip += desc->lrec;  /* index in i/O array */
    op += desc->lrec;  /* index in output array */
    numRecBuff++;      /* count number */
 }
  
  desc->numRecBuff =  numRecBuff; /* How many good */
  return  OBIT_IO_OK;
} /* end ObitIOOTFFITSRead */

/**
 * Read data from disk applying selection and any calibration.
 * The number read will be mySel->nRecPIO (until the end of the selected
 * range of records in which case it will be smaller).
 * The first record number after a read is myDesc->firstRec
 * and the number of records is myDesc->numRecBuff (which
 * may be zero); number attempted is mySel->numRecRead.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param in   Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err  ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSReadSelect (ObitIOOTFFITS *in, ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitOTFDesc* desc;
  ObitOTFSel* sel;
  olong size;
  int status = 0;
  olong len, i, ip, op, numRecBuff;
  gboolean done, OK;
  ofloat *IOBuff = data;
  gchar *routine = "ObitIOOTFFITSReadSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Make sure access was set correctly */
  if (in->access!=OBIT_IO_ReadCal) {
    Obit_log_error(err, OBIT_Error, "%s: access not ReadCal for %s", 
		   routine, in->name);
    return retCode;
  }
  g_assert (ObitIsA(in->myCal, ObitOTFCalGetClass()));
  g_assert (data != NULL);

  desc = in->myDesc; /* OTF descriptor pointer */
  sel  = in->mySel;  /* selector pointer */
  desc->numRecBuff =  0; /* No data yet */

  /* what next ? */
  done = ObitIOOTFFITSNext (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* check if done - all records read */
  if (done) {
    ObitIOOTFFITSClose (in, err); /* Close */
    return OBIT_IO_EOF;
  }
 
  len = desc->lrec; /* How big is a visibility */

  /* read block of sel->numRecRead visibilities at a time  */
  /* transfer size in bytes */
  size = sel->numRecRead * len * sizeof(ofloat); 

  /* Read */
  fits_read_tblbytes (in->myFptr, desc->firstRec, 1, size, (guchar*)IOBuff, 
		      &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR reading FITS OTF data for %s", in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* calibrate/edit/select ... to output */
  ip = op = 0;          /* array pointers */
  numRecBuff = 0; /* How many valid records? */
  for (i=0; i<sel->numRecRead; i++) {

    /* Do byteswap if necessary - this assumes that the data remain the same size */
      /* descriptive  parameters */
      ObitIOOTFFITSfF2H (desc->numDesc, &IOBuff[ip],  &IOBuff[ip]);
      /* data array */
      ObitIOOTFFITSfF2H (desc->colRepeat[desc->ilocdata], &IOBuff[ip+desc->ilocdata], 
			 &IOBuff[ip+desc->ilocdata]);
  
      /* Calibrate and transform */
      OK = ObitOTFCalApply ((ObitOTFCal*)in->myCal, &IOBuff[ip], &data[op], err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      
      ip += desc->lrec;   /* index in i/O array */
      if (OK) { /* at least some of the data unflagged */
	op += desc->lrec;  /* index in output array */
	numRecBuff++;       /* count number */
      }
  } /* end calibration loop */

  desc->numRecBuff =  numRecBuff; /* How many good */
  return  OBIT_IO_OK;
} /* end ObitIOOTFFITSReadSelect */

/**
 * Write information to disk.
 * The data in the buffer will be written starting at record
 * myDesc->firstRec and the number written will be myDesc->numRecBuff
 * which should not exceed mySel->nRecPIO if the internal buffer is used.
 * myDesc->firstRec will be maintained and need not be changed for
 * sequential writing.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSWrite (ObitIOOTFFITS *in, ofloat *data, 
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitOTFDesc* desc;
  ObitOTFSel* sel;
  int status = 0;
  olong size, len, i, ip, op;
  ofloat *IOBuff = data;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  desc = in->myDesc; /* OTF descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  len = desc->lrec; /* How big is a visibility */
  size = desc->numRecBuff * len * sizeof(ofloat); /* transfer size in bytes */

  /* write block of sel->nRecPIO visibilities at a time  */

  /* if output compressed data */
  ip = op = 0; /* array pointers */
  for (i=0; i<desc->numRecBuff; i++) {
    /* byteswap whole record */
    ObitIOOTFFITSfH2F (desc->lrec, &data[op], &IOBuff[ip]);
    ip += desc->lrec;   /* index in i/O array */
    op += desc->lrec;   /* index in output array */
  } /* end byteswap loop */

  /* Write */
  fits_write_tblbytes (in->myFptr, desc->firstRec, 1, size, (guchar*)IOBuff, 
		       &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR writing FITS OTF data for %s", in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    retCode = OBIT_IO_WriteErr;
    return retCode;
  }

  /* keep track of number of records */
  desc->nrecord = MAX (desc->nrecord, desc->firstRec+desc->numRecBuff-1);

  /* where will the next write start */
  desc->firstRec += desc->numRecBuff;

  in->myStatus = OBIT_Modified; /* file has been modified */

  return  OBIT_IO_OK;
} /* end ObitIOOTFFITSWrite */

/**
 * Read image Descriptor data from disk.
 * Also reads the Array Geometry table.
 * \param in Pointer to object with ObitOTFDesc to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSReadDescriptor (ObitIOOTFFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar commnt[FLEN_COMMENT], keyword[FLEN_KEYWORD+1];
  gchar cdata[FLEN_CARD], tunit[FLEN_CARD], TString[12];
  gchar typechar;
  int i, datacol, ncol, nhdu, hdutype, naxis, status = 0;
  double scale, zero;
  long extver, repeat, inaxes[10];
  long ltemp;
  float ftemp;
  double dtemp;
  olong temp;
  ObitOTFDesc* desc;
  ObitOTFSel* sel;
  ObitTableList* tableList;
  gchar *routine = "ObitIOOTFFITSReadDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);

  retCode = OBIT_IO_OK; /* until proven otherwise */

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */
  tableList = (ObitTableList*)in->tableList;

  /* Index tables in file and update TableList if not already done*/
  if (tableList->number <= 0) {
    fits_get_num_hdus (in->myFptr, &nhdu, &status); /* how many? */
    for (i=1; i<=nhdu; i++) {
      fits_movabs_hdu (in->myFptr, i, &hdutype, &status);
      if (hdutype==BINARY_TBL) { /* If it's a table enter it in the list */
	/* table name */
	fits_read_key_str (in->myFptr, "EXTNAME", cdata, commnt, &status);
	/* version number */
	extver = 0;
	fits_read_key_lng (in->myFptr, "EXTVER", &extver, commnt, &status);
	if (status==0) { /* Add to TableList unless it's the main data */
	  if (strcmp (cdata, "OTFScanData")) {
	    temp = (olong)extver;
	    ObitTableListPut (tableList, cdata, &temp, NULL, err);
	    if (err->error)
	      Obit_traceback_val (err, routine, tableList->name, OBIT_IO_OpenErr);
	  }
	}
      }
    } /* end loop indexing file */
  } /* end update Table List */

  /* Position to "OTFScanData" table version 1 */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, "OTFScanData", 1, &status);
  if (status!=0) { /* error? */
     retCode = OBIT_IO_SpecErr;
    Obit_log_error(err, OBIT_Error, "%s ERROR locating OTFScanData table ", routine);
    return retCode;
  }

  /* How many rows in table? */
  fits_get_num_rows (in->myFptr, &ltemp, &status);
  desc->nrecord = ltemp;

  /* How many columns in table? */
  fits_get_num_cols (in->myFptr, &ncol, &status);
  if (status!=0) { /* error? */
     retCode = OBIT_IO_SpecErr;
    Obit_log_error(err, OBIT_Error, "%s ERROR reading OTFScanData table ", routine);
    return retCode;
  }

  /* Get  parameter names, "DATA" should be the last column and contain the detector data. */
  desc->ncol = ncol;  /* number of columns */
  desc->numDesc = ncol-1;
  datacol = (int)ncol;
  for (i=1; i<=desc->numDesc; i++) {
    /* Read column info */
    fits_get_bcolparms (in->myFptr, i, 
			(char*)cdata, (char*)tunit, (char*)&typechar, &repeat, &scale, &zero,
			NULL, NULL, &status);
    /* Get column label, units */
    if (status==0) {
      /* Copy strings fixing bug in cfitsio */
      ObitIOOTFFITSFixBug(desc->colType[i-1], cdata, OTFLEN_KEYWORD);
      ObitIOOTFFITSFixBug(desc->colUnit[i-1], tunit, OTFLEN_VALUE);
    }
    /* Assure that data is a scalar float */
    desc->colRepeat[i-1] = repeat; /* dimensionality */
    if ((typechar!='E') && (repeat!=1)) {
      retCode = OBIT_IO_SpecErr;
      Obit_log_error(err, OBIT_Error, "%s OTF table column %d not correct type/dim ", routine, i);
    }
    /* If this is the TIME parameter get offset */
    if (!strncmp (desc->colType[i-1], "TIME  ",  6)) desc->JDObs = zero;
  } /* end loop reading descriptive parameters */

  /* get visibility dimensionality array */
  fits_read_tdim (in->myFptr, datacol, OTF_MAXDIM, &naxis,
		  inaxes, &status);
  desc->naxis = naxis; /* type mismatch for call argument */
  for (i=0; i<naxis; i++) desc->inaxes[i] = (olong)inaxes[i];

  /* Be sure last column is "DATA  " */
  fits_get_bcolparms (in->myFptr, datacol,
		      cdata, NULL, &typechar, &repeat, &scale, &zero,
		      NULL, NULL, &status);
  ObitIOOTFFITSFixBug(desc->colType[datacol-1], cdata, OTFLEN_KEYWORD);
  desc->colRepeat[datacol-1] = repeat; /* dimensionality */
  if (strncmp (desc->colType[datacol-1], "DATA  ", 6)) {
    retCode = OBIT_IO_SpecErr;
    Obit_log_error(err, OBIT_Error, "%s Last OTF table column NOT DATA ", routine);
  }
  
  /*  Data must be TFLOAT. */
  if (typechar!='E') {
    retCode = OBIT_IO_SpecErr;
    Obit_log_error(err, OBIT_Error, "Illegal DATA type %c", typechar);
  }

  /* Data array definition */

  /* Axis labels */
  for (i=0; i<OTF_MAXDIM; i++) strncpy (desc->ctype[i], "        ", 9);
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCTYP%d", i+1, datacol);
    fits_read_key_str (in->myFptr, (char*)keyword, (char*)cdata, (char*)commnt, &status);
    if (status==0) strncpy (desc->ctype[i], cdata, 9);
    if (status==KEY_NO_EXIST) status = 0;
 }

  /* Axis increments */
  for (i=0; i<OTF_MAXDIM; i++) desc->cdelt[i] = 0.0; /* defaults */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCDLT%d", i+1, datacol);
   ftemp = (float)desc->cdelt[i];
   fits_read_key_flt (in->myFptr, (char*)keyword, &ftemp, (char*)commnt, 
		       &status);
    if (status==KEY_NO_EXIST) status = 0;
    desc->cdelt[i] = (ofloat)ftemp;
  }

  /* Axis reference pixel */
  for (i=0; i<OTF_MAXDIM; i++) desc->crpix[i] = 1.0;/* defaults */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCRPX%d", i+1, datacol);
    ftemp = 0.0;
    fits_read_key_flt (in->myFptr, (char*)keyword, &ftemp, (char*)commnt, 
		       &status);
    if (status==KEY_NO_EXIST) status = 0;
    desc->crpix[i] = (ofloat)ftemp;
  }
  
  /* Axis rotation */
  for (i=0; i<OTF_MAXDIM; i++) desc->crota[i] = 0.0;/* defaults */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCROT%d", i+1, datacol);
    ftemp = 0;
    fits_read_key_flt (in->myFptr, (char*)keyword, &ftemp, (char*)commnt, 
		       &status);
    if (status==KEY_NO_EXIST) status = 0;
    desc->crota[i] = (ofloat)ftemp;
  }

  /* Axis coordinate value at reference pixel */
  for (i=0; i<OTF_MAXDIM; i++) desc->crval[i] = 0.0;
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCRVL%d", i+1, datacol);
    dtemp = 0.0;
    fits_read_key_dbl (in->myFptr, (char*)keyword, &dtemp, (char*)commnt, 
		       &status);
    if (status==KEY_NO_EXIST) status = 0;
    desc->crval[i] = (odouble)dtemp;
 }
  /* descriptive information */
  /* Read keyword values, use default where possible */
  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "EPOCH", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->epoch = ftemp;

  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "EQUINOX", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->equinox = ftemp;
  if ((desc->equinox==0.0) && (desc->epoch>0.0)) desc->equinox = desc->epoch;
  if ((desc->epoch==0.0) && (desc->equinox>0.0)) desc->epoch = desc->equinox;
  
  strncpy (desc->teles, "        ", 8);
  fits_read_key_str (in->myFptr, "TELESCOP", cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->teles, cdata, OTFLEN_VALUE);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->object, "        ", 8);
  fits_read_key_str (in->myFptr, "OBJECT", cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->object, cdata, OTFLEN_VALUE);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->origin, "        ", 8);
  fits_read_key_str (in->myFptr, "ORIGIN", cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->origin, cdata, OTFLEN_VALUE);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->bunit, "        ", 8);
  fits_read_key_str (in->myFptr, "BUNIT", cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->bunit, cdata, OTFLEN_VALUE);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->obsdat, "        ", 9);
  fits_read_key_str (in->myFptr, "DATE-OBS", cdata, (char*)commnt, &status);
  if (status==0)  strncpy (desc->obsdat, cdata, OTFLEN_VALUE); 
  if (status==KEY_NO_EXIST) status = 0;

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "OBSRA", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->obsra = (odouble)dtemp;

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "OBSDEC", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->obsdec = (odouble)dtemp;

  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "ALTRPIX", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->altCrpix = (ofloat)ftemp;

  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "BEAMSIZE", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->beamSize = (ofloat)ftemp;

  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "DIAMETER", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->diameter = (ofloat)ftemp;

  desc->OTFType = OBIT_GBTOTF_Unknown;
  fits_read_key_str (in->myFptr, "OTFTYPE",  TString, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  else {
    desc->OTFType = ObitOTFDescString2Type (TString);
  }

  /* sort order */
  ObitIOOTFFITSSortRead (in, &status);

  /* Look for anything else and add it to the InfoList on desc */
  ObitIOOTFKeysOtherRead(in, &status, err);
  if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);

  /* was there an error? */
  if (err->error) return retCode;
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file header");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* enforce defaults */
  ObitOTFSelDefault(desc, sel);

  return retCode;
} /* end ObitIOOTFFITSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * Also writes the Array Geometry table.
 * \param in  Pointer to object with ObitOTFDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSWriteDescriptor (ObitIOOTFFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar keyword[FLEN_KEYWORD], commnt[FLEN_COMMENT+1];
  gchar *ttype[50], *tform[50], *tunit[50];
  gchar *formRP="1E      ", *DATA = "DATA  ";
  gchar keyName[FLEN_KEYWORD+1], *keyNameP;
  gchar cdata[FLEN_CARD], TString[12];
  int i, k, datacol, tfield, status = 0;
  long naxes[OTF_MAXDIM], extver, nrows;
  gboolean doFill;
  ObitOTFDesc* desc;
  ObitOTFSel* sel;
  ObitInfoType keyType;
  gint32 nkey, dim[MAXINFOELEMDIM];
  union blobEquiv {
    gchar    s[201];
    odouble   d;
    ofloat    f;
    gboolean b;
    oint     o;
    olong    i;
  } blob;
 gchar *routine = "ObitIOOTFFITSWriteDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);

  retCode = OBIT_IO_OK; /* until proven otherwise */

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* enforce defaults */
  ObitOTFSelDefault(desc, sel);

  /* Position to "OTFScanData" table version 1 if it exists */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, "OTFScanData", 1, &status);
  /* Create if it doesn't already exist */
  if (status == BAD_HDU_NUM) {
    status = 0;
    fits_clear_errmsg();   /* Clear error stack */
    /* fill descriptive arrays */
    tfield = desc->ncol;  /* number of columns */
    for (i=0; i<desc->numDesc; i++) {
      tform[i] = formRP;  /* Scalar float */
      ttype[i] = desc->colType[i];  /* column label */
      tunit[i] = desc->colUnit[i];  /* column Units */
    } /* end settingcolumns */

    /* Data column */
    datacol = tfield;
    tunit[tfield-1] = desc->bunit;
    ttype[tfield-1] = DATA;

    /* length and type*/
    g_snprintf (cdata, FLEN_CARD-1, "%dE", desc->colRepeat[tfield-1]);
    tform[tfield-1] = cdata;
    
    /* create table */
    fits_create_tbl (in->myFptr, BINARY_TBL, (long)desc->nrecord, 
		     tfield, ttype, tform, tunit,
		     "OTFScanData", &status);
    if (status!=0) Obit_cfitsio_error(err); 

    /* Add version number 1 */
    extver = 1;
    strncpy (commnt, "Table version number", FLEN_COMMENT);
    fits_update_key_lng (in->myFptr, "EXTVER", extver, (char*)commnt, &status);

    /* Data array definition */
    datacol = (int)desc->ncol; /* must be highest */
  
    /* set data dimensionality array as TDIM */
    for (i=0; i<desc->naxis; i++) naxes[i] = desc->inaxes[i];
    fits_write_tdim (in->myFptr, datacol, (olong)desc->naxis, 
      naxes, &status);
    if (status!=0) Obit_cfitsio_error(err); 

   /* Offset time to JD */
    if (desc->iloct>=0) {
      strncpy (commnt, "Offset of Date from JD", FLEN_COMMENT);
      g_snprintf (keyword, FLEN_KEYWORD-1, "TZERO%d", desc->iloct+1);
      fits_update_key_dbl (in->myFptr, (char*)keyword, (double)desc->JDObs, 12, 
			   (char*)commnt, &status);
    }

  } /* end initialize new table */
  
  /* Truncate table if needed  */
  nrows = (long)desc->nrecord;
  fits_get_num_rows (in->myFptr, &nrows, &status);
  if (nrows > desc->nrecord) {  /* Truncate */
    nrows -= desc->nrecord;
    fits_delete_rows (in->myFptr, desc->nrecord+1, nrows, &status);
  }

  /* Axis labels */
  datacol = (int)desc->ncol;
  strncpy (commnt, "        ", 9);
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCTYP%d", i+1, datacol);
    strncpy (cdata, desc->ctype[i], 9);
    fits_update_key_str (in->myFptr, (char*)keyword, (char*)cdata, (char*)commnt, &status);
  }
  
  /* Axis increments */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCDLT%d", i+1, datacol);
    fits_update_key_flt (in->myFptr, (char*)keyword, (float)desc->cdelt[i], 6, (char*)commnt, 
			 &status);
  }
  
  /* Axis reference pixel */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCRPX%d", i+1, datacol);
    fits_update_key_flt (in->myFptr, (char*)keyword, (float)desc->crpix[i], 6, (char*)commnt, 
			 &status);
  }
  
  /* Axis rotation */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCROT%d", i+1, datacol);
    fits_update_key_flt (in->myFptr, (char*)keyword, (float)desc->crota[i], 6, (char*)commnt, 
			 &status);
  }
  
  /* Axis coordinate value at reference pixel */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCRVL%d", i+1, datacol);
    fits_update_key_dbl (in->myFptr, (char*)keyword, (double)desc->crval[i], 12, (char*)commnt, 
			 &status);
  }
  
  /* descriptive information */
  /* Write keyword values */
  strncpy (commnt, "Originator of file", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "ORIGIN", (char*)desc->origin,  commnt, 
		       &status);
  strncpy (commnt, "Name of object", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "OBJECT", (char*)desc->object,  commnt, 
		       &status);
  strncpy (commnt, "Telescope used", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "TELESCOP", (char*)desc->teles,  commnt, 
		       &status);
  strncpy (commnt, "Date (yyyy-mm-dd) of observation", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "DATE-OBS", (char*)desc->obsdat, (char*)commnt, 
		       &status);
  strncpy (commnt, "Celestial coordiate equinox", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "EPOCH", (float)desc->epoch, 6,  commnt, 
		       &status);

  strncpy (commnt, "Data units", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "BUNIT", (char*)desc->bunit, (char*)commnt, &status);

  strncpy (commnt, "Observed Right Ascension", FLEN_COMMENT);
  fits_update_key_dbl (in->myFptr, "OBSRA", (double)desc->obsra, 12, (char*)commnt, 
		       &status);
  strncpy (commnt, "Observed declination ", FLEN_COMMENT);
  fits_update_key_dbl (in->myFptr, "OBSDEC", (double)desc->obsdec, 12, (char*)commnt, 
		       &status);
  strncpy (commnt, "Beam size (FWHM in deg) ", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "BEAMSIZE", (float)desc->beamSize, 6, (char*)commnt, 
		       &status);
  strncpy (commnt, "Antenna diameter (m) ", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "DIAMETER", (float)desc->diameter, 6, (char*)commnt, 
		       &status);
  strncpy (commnt, "OTF data type ", FLEN_COMMENT);
  ObitOTFDescType2String (desc->OTFType, TString);
  fits_update_key_str (in->myFptr, "OTFTYPE", TString, (char*)commnt, &status);

  /* sort order */
  ObitIOOTFFITSSortWrite (in, &status);

  /* Write other keywords from descriptor */
  nkey = desc->info->number; /* How many keywords? */
  retCode = OBIT_IO_WriteErr;
  strncpy (commnt, "             ", FLEN_COMMENT);
  for (i=0; i<nkey; i++) {
    /* Copy from ObitInfoList */
    ObitInfoListGetNumber(desc->info, i, &keyNameP, &keyType, dim, 
			  blob.s, err);
    if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);
    /* Copy, possibly truncating name */
    strncpy (keyName, keyNameP, FLEN_KEYWORD); keyName[FLEN_KEYWORD] = 0;
    /* blankfill after NULL */
    doFill = FALSE;
    for (k=0; k<8; k++) {
      if (keyName[k]==0) doFill = TRUE;
	    if (doFill) keyName[k] = ' ';
	  }
	  /* Replace any non printing characters with blanks */
	  for (k=0; k<8; k++) if (!g_ascii_isprint(keyName[k])) keyName[k]=' ';
    /* write by type */
    if (keyType==OBIT_double) {
      fits_update_key_dbl (in->myFptr, (char*)keyName, (double)blob.d, 12, (char*)commnt, 
			   &status);
    } else if (keyType==OBIT_float) { 
      fits_update_key_flt (in->myFptr, (char*)keyName, (float)blob.f, 6, (char*)commnt, 
			   &status);
    } else if (keyType==OBIT_string) { 
      blob.s[dim[0]] = 0; /* may not be null terminated */
      fits_update_key_str (in->myFptr, (char*)keyName, (char*)blob.s, (char*)commnt, 
			   &status);
    } else if (keyType==OBIT_oint) { 
      fits_update_key_lng (in->myFptr, (char*)keyName, (long)blob.o, (char*)commnt, 
			   &status);
    } else if (keyType==OBIT_long) { 
      fits_update_key_lng (in->myFptr, (char*)keyName, (long)blob.i, (char*)commnt, 
			   &status);
    } else if (keyType==OBIT_bool) { 
      fits_update_key_log (in->myFptr, (char*)keyName, (int)blob.b, (char*)commnt, 
			   &status);
    }
  } /* end loop writing additional keywords */

  /* was there an error? */
  if (err->error) return retCode;
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR updating FITS file header");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOOTFFITSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOOTFFITSFlush (ObitIOOTFFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* cfitsio does the buffer flushing on close */

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOOTFFITSFlush */

/**
 * Create buffer approptiate for I/O request.
 * Should be called after ObitIO is opened.
 * \param data (output) pointer to data array
 * \param size (output) size of data array in floats.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 */
void 
ObitIOOTFFITSCreateBuffer (ofloat **data, olong *size, 
			  ObitIOOTFFITS *in, ObitInfoList *info, ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  Obit_return_if_fail(((in->myStatus==OBIT_Modified) ||
		       (in->myStatus==OBIT_Active)), 
		      err, "Cannot define buffer, I/O not currently active");

  /* get size */
  *size = ObitOTFSelBufferSize(in->myDesc, in->mySel);

  /* allocate */
  /* (re)allocate */
  if (*data) *data = ObitMemRealloc (*data, (*size)*sizeof(ofloat));
  else *data = ObitMemAlloc0Name((*size)*sizeof(ofloat), "OTFBuffer");

} /* end ObitIOOTFFITSCreateBuffer */

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
 *                 and the highest+1 for OBIT_IO_WriteOnly.
 * \param err      ObitErr for reporting errors.
 * \return pointer to created ObitTable, NULL on failure.
 */
Obit* 
newObitIOOTFFITSTable (ObitIOOTFFITS *in, ObitIOAccess access, 
			 gchar *tabType, olong *tabVer, ObitErr *err)
{
  ObitTable *out;
  olong version;
  gboolean gotIt;
  gchar *outName, tabName[51];
  gchar *routine = "newObitIOOTFFITSTable";

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

  /* Do we already have this one? */
  version = *tabVer;
  gotIt = ObitTableListGet ((ObitTableList*)in->tableList, tabType, &version, 
			    &out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

  /* Check if we're forcing a new table */
  if ((access==OBIT_IO_WriteOnly) && (*tabVer <= 0)) {
    version++;
    out = ObitTableUnref(out);
  }

  /* Set output table version */
  *tabVer = version;

  if (gotIt && (out!=NULL)) return (Obit*)out; /* that was easy */

   /* If it doesn't exist and request is read only - return NULL */
  if ((!gotIt) && (access==OBIT_IO_ReadOnly)) return NULL;

 /* Create one - make descriptive name */
  g_snprintf (tabName, 50, "%s table %d for ",tabType, *tabVer);
  outName =  g_strconcat (tabName, in->name, NULL);
  out = newObitTable (outName);
  g_free(outName);

  /* Setup info needed for access */
  ObitTableSetFITS(out, -1, in->FileName, tabType, version, 
		  25, err);
 if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
 
 /* register it in the TableList */
 ObitTableListPut ((ObitTableList*)in->tableList, tabType, &version, 
		   out, err);
 if (err->error)   Obit_traceback_val (err, routine, in->name, NULL);

 return (Obit*)out;
} /* end newObitIOOTFFITSTable */

/**
 * Update any disk resident structures about the current tables.
 * Nothing is needed for FITS files.
 * \param in   Pointer to object to be updated.
 * \param info ObitInfoList of parent object (not used here).
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOOTFFITSUpdateTables (ObitIOOTFFITS *in, ObitInfoList *info,
				     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  return retCode;
} /* end ObitIOOTFFITSUpdateTables */

/**
 * Get underlying file information in entries to an ObitInfoList
 * Following entries for AIPS files ("xxx" = prefix):
 * \param in      Object of interest.
 * \param myInfo  InfoList on basic object with selection
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 * \param outList InfoList to write entries into
 *
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFileName OBIT_string  FITS file name
 * \li xxxDisk     OBIT_oint    FITS file disk number
 * \li xxxDir      OBIT_string  Directory name for xxxDisk
 *
 * For all File types types:
 * \li xxxFileType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxDataType OBIT_string "AIPS", "FITS"
 *
 * For xxxDataType = "OTF"
 * \li xxxnRecPIO OBIT_int (1,1,1) Number of vis. records per IO call
 * \param err     ObitErr for reporting errors.
 */
void ObitIOOTFFITSGetFileInfo (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
			       ObitInfoList *outList, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *FileType="OTF", *DataType="FITS", *dirname;
  gchar tempStr[201], *here="./";
  olong disk, i;
  gpointer listPnt;
  gchar *parm[] = {"nRecPIO", "doCalSelect", "doCalib", "gainUse", "flagVer",
		   "BChan", "EChan", "Targets", "timeRange", "Scans",
		   "Feeds", "keepCal", "replCal",
		   NULL};
  gchar *routine = "ObitIOOTFFITSGetFileInfo";

  if (err->error) return;

  /* Set basic information */
  if (prefix) keyword =  g_strconcat (prefix, "FileType", NULL);
  else keyword =  g_strdup ("FileType");
  dim[0] = strlen(FileType);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, FileType);
  g_free(keyword);
  
   /* FITS */
  if (prefix) keyword =  g_strconcat (prefix, "DataType", NULL);
  else keyword =  g_strdup ("DataType");
  dim[0] = strlen(DataType);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, DataType);
  g_free(keyword);
  
  /* Filename */
  if (!ObitInfoListGet(myInfo, "FileName", &type, dim, tempStr, err)) 
      Obit_traceback_msg (err, routine, in->name);
  if (prefix) keyword =  g_strconcat (prefix, "FileName", NULL);
  else keyword =  g_strdup ("FileName");
  ObitInfoListAlwaysPut (outList, keyword, type, dim, tempStr);
  g_free(keyword);

  /* Disk number */
  if (!ObitInfoListGet(myInfo, "Disk", &type, dim, &disk, err)) 
      Obit_traceback_msg (err, routine, in->name);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (prefix) keyword =  g_strconcat (prefix, "Disk", NULL);
  else keyword =  g_strdup ("Disk");
  ObitInfoListAlwaysPut (outList, keyword, type, dim, &disk);
  g_free(keyword);
 
  /* Disk directory if disk>0 */
  if (disk>0) dirname = ObitFITSDirname(disk, err); 
  else dirname = here;
  dim[0] = strlen(dirname);
  if (prefix) keyword =  g_strconcat (prefix, "Dir", NULL);
  else keyword =  g_strdup ("Dir");
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, dirname);
  g_free(keyword);


  /* Selection/calibration */
  i = 0;
  while (parm[i]) {
    if (prefix) keyword = g_strconcat (prefix, parm[i], NULL);
    else        keyword = g_strdup(parm[i]);
    if (ObitInfoListGetP(myInfo, parm[i], &type, dim, (gpointer*)&listPnt)) {
      ObitInfoListAlwaysPut(outList, keyword, type, dim, listPnt);
    }
    i++;
    g_free(keyword);
  }
} /* end ObitIOOTFFITSGetFileInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOOTFFITSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOOTFFITSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOOTFFITSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOOTFFITSClassInfoDefFn (gpointer inClass)
{
  ObitIOOTFFITSClassInfo *theClass = (ObitIOOTFFITSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOOTFFITSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOOTFFITSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOOTFFITSGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitIOOTFFITSClear;
  theClass->ObitInit      = (ObitInitFP)ObitIOOTFFITSInit;
  theClass->newObit    = NULL;
  theClass->newObitIO  = (newObitIOFP)newObitIOOTFFITS;
  theClass->ObitIOSame = (ObitIOSameFP)ObitIOOTFFITSSame;
  theClass->ObitIORename = (ObitIORenameFP)ObitIOOTFFITSRename;
  theClass->ObitIOZap  = (ObitIOZapFP)ObitIOOTFFITSZap;
  theClass->ObitCopy   = (ObitCopyFP)ObitIOOTFFITSCopy;
  theClass->ObitClone  = NULL;
  theClass->ObitIOOpen = (ObitIOOpenFP)ObitIOOTFFITSOpen;
  theClass->ObitIOClose= (ObitIOCloseFP)ObitIOOTFFITSClose;
  theClass->ObitIOSet  = (ObitIOSetFP)ObitIOOTFFITSSet;
  theClass->ObitIORead = (ObitIOReadFP)ObitIOOTFFITSRead;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOOTFFITSReadSelect;
  theClass->ObitIOWrite   = 
    (ObitIOWriteFP)ObitIOOTFFITSWrite;
  theClass->ObitIOFlush   = 
    (ObitIOFlushFP)ObitIOOTFFITSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOOTFFITSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOOTFFITSWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOOTFFITSCreateBuffer;
  theClass->newObitIOTable = 
    (newObitIOTableFP)newObitIOOTFFITSTable; 
  theClass->ObitIOUpdateTables   = 
    (ObitIOUpdateTablesFP)ObitIOOTFFITSUpdateTables;
  theClass->ObitIOGetFileInfo   =
    (ObitIOGetFileInfoFP)ObitIOOTFFITSGetFileInfo;

} /* end ObitIOOTFFITSClassDefFn */

/*--------------- Private functions --------------------------*/

/**
 * Creates empty member objects.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOOTFFITSInit  (gpointer inn)
{
  const ObitClassInfo *ParentClass;
  ObitIOOTFFITS *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && (ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->FileName  = NULL;

} /* end ObitIOOTFFITSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOOTFFITSClear (gpointer inn)
{
  ObitIOOTFFITS *in = inn;
  const ObitClassInfo *ParentClass;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->myStatus==OBIT_Active) ||(in->myStatus==OBIT_Modified)) {
    err = newObitErr();
    ObitIOOTFFITSClose (in, err); 
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  if (in->FileName) g_free(in->FileName);  in->FileName = NULL;

 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && (ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOOTFFITSClear */

/**
 * Uses selector member to decide which records to read next.
 * Leaves values in myDesc->firstRec and mySel->numRecRead
 * \param  in   Pointer to the object.
 * \param  err  ObitErr for reporting errors.
 */
static gboolean ObitIOOTFFITSNext (ObitIOOTFFITS *in, ObitErr *err)
{
  ObitOTFDesc* desc;
  ObitOTFSel* sel;
  gboolean done = FALSE;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitIsA(in, &myClassInfo));
 
  sel  = in->mySel;  /* selector pointer */
  desc = in->myDesc; /* OTF descriptor pointer */
 
  /* Let Selector decide */
  done = ObitOTFSelNext (sel, desc, err);

  return done;  
} /* end ObitIOOTFFITSNext */

/**
 * Swaps byte order in floats if host order differs from FITS.
 * Input data assumed in host order and will be swapped to FITS order 
 * (bigendian) if different.  
 * Converts to NaN if blanked.
 * Will work in place.
 * \param  n    Number of floats
 * \param  in   Array of input floats in host byte order.
 * \param  out  Array of output floats in FITS order
 */
static void ObitIOOTFFITSfH2F (olong n, ofloat *in, ofloat *out)
{
  olong i;
  union fequiv inu, outu;
  ofloat fblank = ObitMagicF(), *tmyNaN, myNaN;
  long iNan = ~0;

  /* jump through hoops to get a NaN (all bits on)
     The "standard" routine don't seem to work in gcc */
  tmyNaN = (ofloat*)&iNan;
  myNaN = *tmyNaN;

#if G_BYTE_ORDER==G_BIG_ENDIAN  /* no byte swap needed */
  /* if the input and output point to the same place - just return */
  if (in==out) return;
  for (i=0; i<n; i++) {
    if (in[i]!=fblank) out[i] = in[i]; /* simple copy */
    else out[i] = myNaN;
  }

#elif G_BYTE_ORDER==G_LITTLE_ENDIAN   /* byte swap */
  for (i=0; i<n; i++) {
    if (in[i]!=fblank) inu.full = in[i];
    else inu.full = myNaN;
    /*inu.full = in[i];*/
    outu.parts[0] = inu.parts[3]; 
    outu.parts[1] = inu.parts[2]; 
    outu.parts[2] = inu.parts[1]; 
    outu.parts[3] = inu.parts[0]; 
    out[i] = outu.full;
  }

#else /* unknown */
  g_error("ObitIOOTFFITSfH2F: Unsupported host byte order");
#endif
} /* end ObitIOOTFFITSfH2F */

/**
 * Swaps byte order in floats if host order differs from FITS.
 * Input data assumed in FITS byte order (bigendian) and will be
 * swapped to host order if different.
 * Converts NaNs to fblanks
 * Will work in place.
 * \param  n    Number of floats
 * \param  in   Array of input floats in FITS order
 * \param  out  Array of output floats in host byte order.
 */
static void ObitIOOTFFITSfF2H (olong n, ofloat *in, ofloat *out)
{
  olong i;
  union fequiv inu, outu;
  ofloat fblank = ObitMagicF();

#if G_BYTE_ORDER==G_BIG_ENDIAN  /* no byte swap needed */
  /* if the input and output point to the same place - just return */
  if (in==out) return;
  for (i=0; i<n; i++) {
    if (isnan(in[i]))  out[i] = fblank;
    else out[i] = in[i]; /* simple copy */
  }

#elif G_BYTE_ORDER==G_LITTLE_ENDIAN   /* byte swap */
  for (i=0; i<n; i++) {
    inu.full = in[i];
    outu.parts[0] = inu.parts[3]; 
    outu.parts[1] = inu.parts[2]; 
    outu.parts[2] = inu.parts[1]; 
    outu.parts[3] = inu.parts[0]; 
    out[i] = outu.full;
    if (isnan(out[i])) out[i] = fblank;
  }

#else /* unknown */
  g_error("ObitIOOTFFITSfF2H: Unsupported host byte order");
#endif
} /* end ObitIOOTFFITSfF2H */

/**
 * Look for keyword SORTORD for the Sort Order.
 * \param in      Pointer to ObitIOOTFFITS.
 * \param status (Output) cfitsio status.
 * \return return code, 0=> OK
 */
static void  ObitIOOTFFITSSortRead(ObitIOOTFFITS *in, olong *status)
{
  gchar commnt[FLEN_COMMENT], cdata[FLEN_COMMENT];
  ObitOTFDesc *desc;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  desc = in->myDesc; /* set descriptor */

  /* init */
  desc->isort[0] = ' '; desc->isort[1] = ' '; desc->isort[2] = 0; 

  /* Attempt rational keywords */
  fits_read_key_str (in->myFptr, "SORTORD", (char*)cdata, (char*)commnt, status);
  if (*status==0) {
    desc->isort[0]=cdata[0]; 
    desc->isort[1]=cdata[1]; 
    return;
  }
  if (*status==KEY_NO_EXIST) *status = 0;

} /* end ObitIOOTFFITSSortRead */

/**
 * Write SORTORD  keyword for  Descriptor value isort
 * \param in Pointer to ObitIOOTFFITS.
 * \param status (Output) cfitsio status.
 * \return return code, 0=> OK
 */
static void  ObitIOOTFFITSSortWrite (ObitIOOTFFITS *in, olong *status)
{
  gchar commnt[FLEN_COMMENT+1];
  ObitOTFDesc *desc;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  desc = in->myDesc; /* set descriptor */
  /* Don't bother if not given */
  if (((desc->isort[0]==' ') && (desc->isort[1]==' ')) ||
      (desc->isort[0]==0)) return;
 
  /* Rational keyword */
  strncpy (commnt, "Sort Order code ", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "SORTORD", (char*)desc->isort, (char*)commnt, 
		       status);

} /* end ObitIOOTFFITSSortWrite */

/**
 * Look for additional descriptive keywords, any that are not 
 * on the exclusion list are copied to the descriptor InfoList.
 * \param in      Pointer to ObitIOOTFFITS.
 * \param status (Output) cfitsio status.
 * \param err    ObitErr stack.
 * \return return code, 0=> OK
 */
void  ObitIOOTFKeysOtherRead(ObitIOOTFFITS *in, olong *lstatus, 
			       ObitErr *err)
{
  gchar keywrd[FLEN_KEYWORD], value[FLEN_VALUE], commnt[FLEN_COMMENT+1];
  gchar *first, *last, *anF, *aT, dtype, svalue[FLEN_VALUE];
  int i, j, k, keys, morekeys, status=(int)*lstatus;
  olong ivalue;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  double dvalue;
  ObitOTFDesc *desc;
  gchar *exclude[] = 
  {"SIMPLE", "BITPIX", "EXTEND", "HISTORY", "COMMENT", "BLANK", "        ",
   "XTENSION", "PCOUNT", "GCOUNT", "EXTNAME", "EXTVER", "SORTORD",
   "TFORM", "TTYPE", "TUNIT", "TSCAL", "TZERO", "TDIM", "TDISP",
   "1CTYP", "2CTYP", "3CTYP", "4CTYP", "5CTYP", "6CTYP", "7CTYP", 
   "1CRVL", "2CRVL", "3CRVL", "4CRVL", "5CRVL", "6CRVL", "7CRVL", 
   "1CDLT", "2CDLT", "3CDLT", "4CDLT", "5CDLT", "6CDLT", "7CDLT", 
   "1CRPX", "2CRPX", "3CRPX", "4CRPX", "5CRPX", "6CRPX", "7CRPX", 
   "1CROT", "2CROT", "3CROT", "4CROT", "5CROT", "6CROT", "7CROT", 
   "BSCALE", "BZERO", "NAXIS", "TFIELDS",
   "CTYPE", "CDELT", "CRPIX", "CROTA", "CRVAL", "OBSRA", "OBSDEC", 
   "OBJECT", "TELESCOP", "DATE", "EPOCH", "DATAMAX", "DATAMIN", "BUNIT", 
   "ALTRVAL", "ALTRPIX", "VELREF", "RESTFREQ", "XSHIFT", "YSHIFT", 
   "CLEAN", NULL};
  olong number, *len;
  gboolean bvalue, bad=FALSE;
  gchar *routine = "ObitIOOTFKeysOtherRead";

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete old InfoList and restart */
  ((ObitOTFDesc*)in->myDesc)->info = 
    ObitInfoListUnref (((ObitOTFDesc*)in->myDesc)->info);
  ((ObitOTFDesc*)in->myDesc)->info = (gpointer)newObitInfoList ();
  desc = in->myDesc; /* set descriptor */

  /* get number and length of exclusion strings */
  number = 0;
  i = 0;
  while (exclude[i]!=NULL) {
    number++;
    i++;
  }
  len = g_malloc0(number*sizeof(olong));
  for (i=0; i<number; i++) len[i] = strlen(exclude[i]);

  /* how many keywords to look at? */
  fits_get_hdrspace (in->myFptr, &keys, &morekeys, &status);
  for (k=1; k<=keys; k++) {
    fits_read_keyn (in->myFptr, k, keywrd, value, (char*)commnt, &status);
    if (status==0) {
      /* Is this on the list? */
      for (j=0; j<number; j++) {
	bad = !strncmp (keywrd, exclude[j], len[j]);
	bad = bad || (strlen(keywrd)<=0); /* blank keyword */
	if (bad) break;
      } /* end loop over exclusions */
      /* want this one? */

      if (!bad) {
	/* ask cfitsio what it is */
	fits_get_keytype (value, &dtype, &status);
	switch (dtype) { 
	case 'C':  /* Character string */
	  first = index (value,'\'')+1; /* a string? */
	  last = rindex(value,'\'')-1;
	  g_memmove(svalue, first, (last-first+1));
	  svalue[last-first+1] = 0; /* null terminate */
	  /* add to InfoList */
	  dim[0] = strlen(svalue);
	  ObitInfoListPut(desc->info, keywrd, OBIT_string, dim, 
			  (gconstpointer)svalue, err);
	  
	  break;
	case 'L':  /* logical 'T', 'F' */
	  anF   = index (value,'F'); /* Logical */
	  aT    = index (value,'T'); /* Logical */
	  bvalue = FALSE;
	  if (aT!=NULL) bvalue = TRUE;
	  /* add to InfoList */
	  dim[0] = 1;
	  ObitInfoListPut(desc->info, keywrd, OBIT_bool, dim, 
			  (gconstpointer)&bvalue, err);
	  break;
	case 'I':  /* Integer */
	  ivalue = strtol(value, NULL, 10);
	  /* add to InfoList */
	  dim[0] = 1;
	  ObitInfoListPut(desc->info, keywrd, OBIT_oint, dim, 
			  (gconstpointer)&ivalue, err);
	  break;
	case 'F':  /* Float - use double */
	  dvalue = strtod(value, &last);
	  /* add to InfoList */
	  dim[0] = 1;
	  ObitInfoListPut(desc->info, keywrd, OBIT_double, dim, 
			  (gconstpointer)&dvalue, err);
	  break;
	case 'X':  /* Complex - can't handle */
	default:
	  g_assert_not_reached(); /* unknown, barf */
	}; /* end switch on type */
	
	/* error check */
	if (err->error) Obit_traceback_msg (err, routine, in->name);
      }
    }
  } /* end loop over keywords */

} /* end ObitIOOTFKeysOtherRead */

/**
 * Work around for cfitsio bug, trailing blanks in keywords are dropped.
 * This routine blank fills out and null terminates at out[maxn-1].
 * \param out    Output string 
 * \param in     Input string from cfitsio keyword read
 * \param maxn   length of out
 */
void ObitIOOTFFITSFixBug (gchar *out, gchar *in, olong maxn)
{
  olong i, len;

  len = strlen(in);
  for (i=0; i<len; i++)    out[i] = in[i];
  for (i=len; i<maxn; i++) out[i] = ' ';
  out[maxn-1] = 0;
  
} /* end ObitIOOTFFITSFixBug */
