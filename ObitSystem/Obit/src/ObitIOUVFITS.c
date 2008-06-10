/* $Id: ObitIOUVFITS.c,v 1.38 2008/02/18 16:29:15 bcotton Exp $      */
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

#include <errno.h>
#include "ObitIOUVFITS.h"
#include "ObitTableList.h"
#include "ObitFile.h"
#include "ObitFITS.h"
#include "ObitFileFITS.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOUVFITS.c
 * ObitIOUVFITS class function definitions.
 * This class is derived from the ObitIO class.
 *
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOUVFITS";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOUVFITSClassInfo myClassInfo = {FALSE};

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
void  ObitIOUVFITSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOUVFITSClear (gpointer in);

/** Private: Determine next visibility to read */
static gboolean ObitIOUVFITSNext (ObitIOUVFITS *in, ObitErr *err);

/** Private: Compress visibilities. */
static void 
ObitIOUVFITSCompress (gint ncorr, const ofloat *visin, ofloat *wtscl, 
		      ofloat *visout);

/** Private: Uncompress visibilities. */
static void 
ObitIOUVFITSUncompress (gint ncorr, const ofloat *visin, 
			const ofloat *wtscl, ofloat *visout);

/** Private: Copy Floats with byte swap to FITS order */
static void ObitIOUVFITSfH2F (gint n, ofloat *in, ofloat *out);

/** Private: Copy Floats with byte swap to host order */
static void ObitIOUVFITSfF2H (gint n, ofloat *in, ofloat *out);

/** Private: Copy Shorts with byte swap to FITS order */
static void ObitIOUVFITSsH2F (gint n, gshort *in, gshort *out);

/** Private: Copy Shorts with byte swap to host order */
static void ObitIOUVFITSsF2H (gint n, gshort *in, gshort *out);

/** Private: Read AIPS (and other) Sort Order. */
static void  ObitIOUVFITSSortRead(ObitIOUVFITS *in, olong *lstatus);

/** Private: Write  AIPS (and other) Sort Order.*/
static void  ObitIOUVFITSSortWrite (ObitIOUVFITS *in, ObitErr *err);

/** Private: Copy other header keywords. */
void  ObitIOUVKeysOtherRead(ObitIOUVFITS *in, olong *lstatus, 
			       ObitErr *err);

/** Private: Fix bug in cfitsio keyword parsing */
void ObitIOUVFITSFixBug(gchar *out, gchar *in, olong maxn);

/** Private: Write AIPSish main file header */
static ObitIOCode WriteAIPSUVHeader (ObitIOUVFITS *in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitIOUVFITSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class on the first call.
 * \param name An optional name for the object.
 * \param info if non-NULL it is used to initialize the new object.
 * \param err  ObitErr for error messages.
 * \return the new object.
 */
ObitIOUVFITS* newObitIOUVFITS (gchar *name, ObitInfoList *info,
			       ObitErr *err)
{
  ObitIOUVFITS* out;
  gint32 dim[UV_MAXDIM];
  ObitInfoType type;
  gchar tempStr[201];
  gchar *routine = "newObitIOUVFITS";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitIOUVFITSClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitIOUVFITS));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOUVFITSInit((gpointer)out);

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
} /* end newObitIOUVFITS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOUVFITSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOUVFITSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIOUVFITSGetClass */

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
gboolean ObitIOUVFITSSame (ObitIO *in, ObitInfoList *in1, 
			   ObitInfoList *in2, ObitErr *err)
{
  olong disk1, disk2;
  gchar *filename1, *filename2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOUVFITSSame";

  /* error checks */
  if (err->error) return same;
  errno = 0;  /* reset any system error */

  /* get file from info */
  if(!ObitInfoListGet(in1, "Disk", &type, dim, &disk1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in1, "FileName", &type, dim, (gpointer)&filename1)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }

  if(!ObitInfoListGet(in2, "Disk", &type, dim, &disk2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in2, "FileName", &type, dim, (gpointer)&filename2)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }

  /* Compare */
  same = (disk1==disk2) && 
    !strncmp (filename1,filename2, 200);
 
  return same;
} /* end ObitIOUVFITSSame */

/**
 * Rename underlying files.
 * New name information is given on the info member:
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 * \param in   Pointer to object to be renamed
 * \param info Associated ObitInfoList
 * \param err   ObitErr for reporting errors.
 */
void ObitIOUVFITSRename (ObitIO *in, ObitInfoList *info, 
			    ObitErr *err)
{
  gchar *routine = "ObitIOUVFITSRename";

  /* Don't bother if NULL */
  if (!in) return;

  /* error checks */
  if (err->error) return;
  errno = 0;  /* reset any system error */

  /* Rename */
  ObitFITSRename (in, info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitIOUVFITSRename */

/**
 * Delete underlying files.
 * Delete the whole FITS file.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOUVFITSZap (ObitIOUVFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar tempStr[201];
  gchar *routine = "ObitIOUVFITSZap";

   /* error check */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* Close if still open */
  if ((in->myStatus==OBIT_Modified) || (in->myStatus==OBIT_Active)) {
    retCode = ObitIOUVFITSClose (in, err);
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
    ObitFileErrMsg(err);     /* system error message*/
  }

  /* FITS tables are in the same file - delete table list */
  in->tableList = ObitTableListUnref(in->tableList);

  return;
} /* end ObitIOUVFITSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOUVFITS* ObitIOUVFITSCopy  (ObitIOUVFITS *in, 
				       ObitIOUVFITS *out, ObitErr *err)
{
  const ObitIOClassInfo *ParentClass;
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
    out = newObitIOUVFITS(outName, NULL, err);
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
} /* end ObitIOUVFITSCopy */

/**
 * Initialize structures and open file.
 * The file etc. info should have been stored in the ObitInfoList.
 * The image descriptor is read if ReadOnly or ReadOnly and
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
ObitIOCode ObitIOUVFITSOpen (ObitIOUVFITS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;
  long naxes[2] = {777777701, 0};
  gint32 dim[UV_MAXDIM];
  ObitInfoType type;
  gchar tempStr[201];
  ObitUVDesc* desc;
  gchar *routine = "ObitIOUVFITSOpen";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  g_assert (in->myDesc != NULL);
  g_assert (in->mySel != NULL);
  errno = 0;  /* reset any system error */

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
  if (err->error) Obit_traceback_val (err, routine, in->FileName, retCode);
  
  /* open file by access type */
  /*------------------------ Read Only ---------------------------------*/
  if (access == OBIT_IO_ReadOnly) {
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);
    ObitTrimTrail(tempStr);  /* Trim any trailing blanks */

    /* cfitsio refuses to open a file readwrite after it has been opened
       readonly so first try opening readwrite even if requested ReadOnly.
       If that fails try readonly. */
    /* Test open readwrite */
    fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status);
    if ((status==FILE_NOT_OPENED) || (status==READONLY_FILE)) {
      /* Failed - try readonly */
      status = 0;
      fits_clear_errmsg();   /* Clear cfitsio error stack */
      if (fits_open_file(&(in->myFptr), (char*)tempStr, READONLY, &status) ) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR %d opening input FITS file %s", status, in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
    }
    /* Other errors? */
    if (status!=0) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening output FITS file %s", 
		     in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      retCode = OBIT_IO_OpenErr;
      return retCode;
    }

  /*------------------------ Read/Write ---------------------------------*/
} else if ((access == OBIT_IO_ReadWrite)  || (access == OBIT_IO_ReadCal)) {
    /* Output file may already exist test open */
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);
    ObitTrimTrail(tempStr);  /* Trim any trailing blanks */

    /* Initialize output file */
    if ( fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status) ) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR %d opening output FITS file %s", 
		     status, in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
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
    ObitTrimTrail(tempStr);  /* Trim any trailing blanks */

    /* Open read/write to see if it's there */
    fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status);

    if (status==FILE_NOT_OPENED) { /* Not there - initialize output file */
      status = 0;
      fits_clear_errmsg();   /* Clear error stack */
      if (fits_create_file(&(in->myFptr), (char*)in->FileName, &status) ) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR %d opening output FITS file %s", 
		       status, in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }

      /* Secret Header so AIPS will recognize it */
      fits_write_imghdr (in->myFptr, 8, 2, naxes, &status);
      fits_modify_comment (in->myFptr, "NAXIS1", "AIPS secret sign", &status);
      if (status!=0) { /* error */
	Obit_log_error(err, OBIT_Error, 
		       "ERROR %d opening output FITS file %s", 
		       status, in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
    } else if (status!=0) { /* error */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR %d opening output FITS file %s", 
		     status, in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
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
  desc->firstVis   = 0;
  desc->numVisBuff = 0;
  
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOUVFITSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOUVFITSClose (ObitIOUVFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */
  /* don't bother if it's not open */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) 
    return OBIT_IO_OK;

  fits_close_file (in->myFptr, &status);
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR %d closing FITS file", status);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_CloseErr;
    return retCode;
  }

  /* Delete any compression buffers */
  if (in->compBuff) in->compBuff = ObitMemFree (in->compBuff);  
  in->compBuffSize = 0;
  if (in->decompVis) in->decompVis = ObitMemFree (in->decompVis);  

  in->myStatus = OBIT_Inactive;
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOUVFITSClose */

/**
 * initialize I/O, position to beginning of uvdata.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVFITSSet (ObitIOUVFITS *in, ObitInfoList *info, 
			       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  int status = 0;

  /* reset vis pointers */
  ((ObitUVDesc*)(in->myDesc))->firstVis   = 0;
  ((ObitUVDesc*)(in->myDesc))->numVisBuff = 0;

  /* Position to "AIPS UV" table version 1 */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, "AIPS UV", 1, &status);
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR %d positioning FITS file", status);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    if ((in->access==OBIT_IO_ReadOnly) ||
	(in->access==OBIT_IO_ReadWrite)) retCode = OBIT_IO_ReadErr;
    else  retCode = OBIT_IO_WriteErr;
    return retCode;
  }

  return retCode;
} /* end ObitIOUVFITSSet */

/**
 * Read data from disk.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities attempted is mySel->numVisRead; 
 * actual value saved as myDesc->numVisBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVFITSRead (ObitIOUVFITS *in, ofloat *data, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  long size;
  int status = 0;
  olong len, i, ip, op, need;
  gboolean done, compressed;
  ofloat wtscl[2], *IOBuff = data;
  gchar *routine = "ObitIOUVFITSRead";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);
  errno = 0;  /* reset any system error */

  desc = in->myDesc; /* UV descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* what next ? */
  done = ObitIOUVFITSNext (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* check if done - all visibilities read */
  if (done) {
    /* ObitIOUVFITSClose (in, err); Close */
    desc->numVisBuff = 0; /* no data in buffer */
    return OBIT_IO_EOF;
  }

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;
  if (compressed) {
    if (in->compBuff==NULL) {
      in->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
      in->compBuff     = ObitMemAllocName (in->compBuffSize, "UVFITS comp. buff");
      
      /* check that weight and scale are available */
      if (desc->ilocws < 0) {
	Obit_log_error(err, OBIT_Error, 
		       "UV data does not have weight and scale %s", 
		       in->name);
	return retCode;
      }
    }
    IOBuff = in->compBuff; /* Use for I/O */
    /* Make sure compression buffer large enough (sel->nVisPIO) */
    need = desc->lrec*sel->numVisRead*sizeof(ofloat); 
    if (need > in->compBuffSize) {
      Obit_log_error(err, OBIT_Error, 
		     "Decompression buffer (%d) too small, need %d for %s", 
		     in->compBuffSize, need, in->name);
      return retCode;
    }
  } /* end compression setup */
  
  len = desc->lrec; /* How big is a visibility */

  /* read block of sel->numVisRead visibilities at a time  */
  /* transfer size in bytes */
  size = (long)sel->numVisRead * len * sizeof(ofloat); 

  /* Read */
  fits_read_tblbytes (in->myFptr, (long)desc->firstVis, 1, size, (guchar*)IOBuff, 
		      &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR %d reading FITS uv data for %s", status, in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* if compressed data uncompress to output buffer */
  if (compressed ) {
    ip = op = 0; /* array pointers */
    for (i=0; i<sel->numVisRead; i++) {
      /* Copy random parameters */
      ObitIOUVFITSfF2H (sel->nrparmUC, &IOBuff[ip], &data[op]);

      /* uncompress/byte swap data */
      ObitIOUVFITSfF2H (2, &IOBuff[ip+desc->ilocws], wtscl);
      ObitIOUVFITSsF2H (2*desc->ncorr, 
			(gshort*)&IOBuff[ip+desc->nrparm],
			(gshort*)&IOBuff[ip+desc->nrparm]);
      ObitIOUVFITSUncompress (desc->ncorr, &IOBuff[ip+desc->nrparm], 
			      wtscl, &data[op+sel->nrparmUC]);
      ip += desc->lrec;   /* index in i/O array */
      op += sel->lrecUC;  /* index in output array */
    } /* end decompression loop */

  } else { /* no decompression, just byte swap */
    ip = op = 0; /* array pointers */
    for (i=0; i<sel->numVisRead; i++) {
      /* Copy random parameters */
      ObitIOUVFITSfF2H (sel->nrparmUC, &IOBuff[ip], &data[op]);
      /* Copy visibility array */
      ObitIOUVFITSfF2H (3*desc->ncorr, &IOBuff[ip+desc->nrparm],
		       &data[op+sel->nrparmUC]);
      ip += desc->lrec;   /* index in i/O array */
      op += sel->lrecUC;  /* index in output array */
    }
  } /* end decompression/byte swap */

  /* how many read */
  desc->numVisBuff = sel->numVisRead;
  
  return  OBIT_IO_OK;
} /* end ObitIOUVFITSRead */

/**
 * Read data from disk applying selection and any calibration.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff (which
 * may be zero); number attempted is mySel->numVisRead.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVFITSReadSelect (ObitIOUVFITS *in, ofloat *data, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  olong size;
  int status = 0;
  olong len, i, k, ip, op, need, numVisBuff;
  gboolean done, compressed, OK;
  ofloat *workVis, *wtscl, *IOBuff = data;
  gchar *routine = "ObitIOUVFITSReadSelect";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */
  /* Make sure access was set correctly */
  if (in->access!=OBIT_IO_ReadCal) {
    Obit_log_error(err, OBIT_Error, 
		   "ObitIOUVFITSReadSelect: access not ReadCal for %s", 
		   in->name);
    return retCode;
  }
  g_assert (ObitUVCalIsA((ObitUVCal*)in->myCal));
  g_assert (data != NULL);

  desc = in->myDesc; /* UV descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  desc->numVisBuff = 0; /* no data in buffer yet */
  
  /* what next ? */
  done = ObitIOUVFITSNext (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* check if done - all visibilities read */
  if (done) {
    /* ObitIOUVFITSClose (in, err); Close */
    return OBIT_IO_EOF;
  }

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;

  /* Always use compression buffer as the output vis may have a 
     different size from the input */
  if (in->compBuff==NULL) {
    in->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
    in->compBuff = ObitMemAllocName (in->compBuffSize, "UVFITS comp. buff");
    
    /* check that weight and scale are available */
    if (compressed && (desc->ilocws < 0)) {
      Obit_log_error(err, OBIT_Error, 
		     "UV data does not have weight and scale %s", 
		     in->name);
      return retCode;
    }
  }

  IOBuff = in->compBuff; /* Use for I/O */
  /* Make sure compression buffer large enough (sel->nVisPIO) */
  need = desc->lrec*sel->numVisRead*sizeof(ofloat);
  if (need > in->compBuffSize) {
    Obit_log_error(err, OBIT_Error, 
		   "Decompression buffer (%d) too small, need %d for %s", 
		   in->compBuffSize, need, in->name);
    return retCode;
  }

  /* Visibility decompression buffer */
  if (compressed && (in->decompVis==NULL)) {
    /* conservative guess of the size of the record */
    need = 100 + desc->nrparm + (desc->lrec-desc->nrparm)*3;
    in->decompVis =  ObitMemAllocName (need*sizeof(ofloat), "UVFITS comp. vis");
  }
  len = desc->lrec; /* How big is a visibility */

  /* read block of sel->numVisRead visibilities at a time  */
  /* transfer size in bytes */
  size = (long)(sel->numVisRead * len * sizeof(ofloat)); 

  /* Read */
  fits_read_tblbytes (in->myFptr, (long)desc->firstVis, 1, size, (guchar*)IOBuff, 
		      &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR %d reading FITS uv data for %s", status, in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* uncompress/calibrate/edit/select transform... to output */
  ip = op = 0;          /* array pointers */
  numVisBuff = 0; /* How many valid visibilities */
  for (i=0; i<sel->numVisRead; i++) {

    /* Do byteswap if necessary */
      /* Random parameters */
      ObitIOUVFITSfF2H (desc->nrparm, &IOBuff[ip], &IOBuff[ip]);
      /* Visibility array */
      if (compressed ) {
	ObitIOUVFITSsF2H (2*desc->ncorr, 
			  (gshort*)&IOBuff[ip+desc->nrparm],
			  (gshort*)&IOBuff[ip+desc->nrparm]);
      } else { /* uncompressed */
	ObitIOUVFITSfF2H (3*desc->ncorr, &IOBuff[ip+desc->nrparm],
			 &IOBuff[ip+desc->nrparm]);
      }

    /* Decompress */
    if (compressed) {
      /* copy random parameters to visibility work buffer */
      for (k=0; k<sel->nrparmUC; k++) in->decompVis[k] = IOBuff[ip+k];
      
      /* uncompress data */
      wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
      ObitIOUVFITSUncompress (desc->ncorr, &IOBuff[ip+desc->nrparm], 
			      wtscl, &in->decompVis[sel->nrparmUC]);
      workVis = in->decompVis; /* working visibility pointer */

    } else {
      /* Data not compressed - work out of I/O buffer */
      workVis = &IOBuff[ip];
    }
    
    /* Calibrate and transform */
    OK = ObitUVCalApply ((ObitUVCal*)in->myCal, workVis, &data[op], err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    ip += desc->lrec;   /* index in i/O array */
    if (OK) { /* at least some of the data unflagged */
      op += sel->lrecUC;  /* index in output array */
      numVisBuff++;       /* count number */
    }
  } /* end calibration loop */

  desc->numVisBuff =  numVisBuff; /* How many good */
  return  OBIT_IO_OK;
} /* end ObitIOUVFITSReadSelect */

/**
 * Write information to disk.
 * The data in the buffer will be written starting at visibility
 * myDesc->firstVis and the number written will be myDesc->numVisBuff
 * which should not exceed mySel->nVisPIO if the internal buffer is used.
 * myDesc->firstVis will be maintained and need not be changed for
 * sequential writing.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVFITSWrite (ObitIOUVFITS *in, ofloat *data, 
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  int status = 0;
  olong size, len, i, j, ip, op, need;
  gboolean compressed;
  ofloat *wtscl, *IOBuff = data;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);
  errno = 0;  /* reset any system error */

  desc = in->myDesc; /* UV descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;
  if (sel->Compress) {
    /* buffer if necessary create */
    if (in->compBuff==NULL) {
      in->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
      in->compBuff = ObitMemAllocName (in->compBuffSize, "UVFITS comp. buff");

      /* check that weight and scale are available */
      if (desc->ilocws < 0) {
       Obit_log_error(err, OBIT_Error, 
		     "UV data does not have weight and scale %s", 
		      in->name);
       return retCode;
     }
    }
    IOBuff = in->compBuff; /* Use for I/O */
    /* Make sure compression buffer large enough (sel->nVisPIO) */
    need = desc->lrec*desc->numVisBuff*sizeof(ofloat);
    if (need > in->compBuffSize) {
      Obit_log_error(err, OBIT_Error, 
		     "Compression buffer (%d) too small, need %d for %s", 
		     in->compBuffSize, need, in->name);
      return retCode;
    }
  } /* end of compressed data set up */

  len = desc->lrec; /* How big is a visibility */
  size = desc->numVisBuff * len * sizeof(ofloat); /* transfer size in bytes */

  /* write block of sel->nVisPIO visibilities at a time  */

  /* if output compressed data */
  if (compressed) {
    ip = op = 0; /* array pointers */
    for (i=0; i<desc->numVisBuff; i++) {
       /* Copy random parameters */
      for (j=0; j<sel->nrparmUC; j++) IOBuff[ip+j] = data[op+j];
     
      /* compress/byteswap data */
      wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
      ObitIOUVFITSCompress (desc->ncorr, &data[op+sel->nrparmUC], 
			    wtscl, &IOBuff[ip+desc->nrparm]);
      /* Byteswap Random parameters */
      ObitIOUVFITSfH2F (desc->nrparm, &IOBuff[ip], &IOBuff[ip]);
      ip += desc->lrec;   /* index in i/O array (compressed) */
      op += sel->lrecUC;  /* index in input array */
    } /* end compression loop */

  } else { /* no compression, just byte swap */
    ip = op = 0; /* array pointers */
    for (i=0; i<desc->numVisBuff; i++) {
      /* byteswap whole visibility */
      ObitIOUVFITSfH2F (desc->lrec, &data[op], &IOBuff[ip]);
      ip += desc->lrec;   /* index in i/O array */
      op += sel->lrecUC;  /* index in output array */
    } /* end byteswap loop */
  } /* end compression / byteswap */

  /* Write*/
  fits_write_tblbytes (in->myFptr, (long)desc->firstVis, 1, size, (guchar*)IOBuff, 
		       &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR %d writing FITS uv data for %s", status, in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_WriteErr;
    return retCode;
  }

  /* keep track of number of visibilities */
  desc->nvis = MAX (desc->nvis, desc->firstVis+desc->numVisBuff-1);

  /* where will the next write start */
  desc->firstVis += desc->numVisBuff;

  in->myStatus = OBIT_Modified; /* file has been modified */

  return  OBIT_IO_OK;
} /* end ObitIOUVFITSWrite */

/**
 * Read image Descriptor data from disk.
 * \param in Pointer to object with ObitUVDesc to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVFITSReadDescriptor (ObitIOUVFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar commnt[FLEN_COMMENT], keyword[FLEN_KEYWORD+1];
  gchar cdata[FLEN_CARD], tunit[FLEN_CARD];
  gchar typechar[4], *today=NULL;
  int i, viscol, ncol, naxis, nhdu, hdutype, status = 0, xstatus = 0;
  double scale, zero;
  long nvis, repeat, extver, inaxes[10];
  olong temp;
  long ltemp;
  float ftemp;
  double dtemp;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  ObitTableList* tableList;
  gchar *routine = "ObitIOUVFITSReadDescriptor";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */

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
	fits_read_key_str (in->myFptr, "EXTNAME", (char*)cdata, (char*)commnt, &status);
	/* version number default to 0 */
	extver = 0;
	fits_read_key_lng (in->myFptr, "EXTVER", &extver, (char*)commnt, &xstatus);
	if (status==0) { /* Add to TableList unless it's the uv data */
	  if (strcmp (cdata, "AIPS UV")) {
	    temp = (olong)extver;
	    ObitTableListPut (tableList, cdata, &temp, NULL, err);
	    if (err->error)
	      Obit_traceback_val (err, routine, tableList->name, OBIT_IO_OpenErr);
	  }
	}
      }
    } /* end loop indexing file */
  } /* end update Table List */

  /* Position to "AIPS UV" table version 1 */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, "AIPS UV", 1, &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR %d finding UV table - may not be UV data", status);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* How many rows in table? */
  fits_get_num_rows (in->myFptr, &nvis, &status);
  desc->nvis = (olong)nvis;

  /* How many columns in table? */
  fits_get_num_cols (in->myFptr, &ncol, &status);

  /* Get random parameter names, "VISIBILITIES" should be the last
     column and contain the visibility data. */
  desc->nrparm = (olong)(ncol-1);
  for (i=1; i<=desc->nrparm; i++) {
    /* Read column info */
    fits_get_bcolparms (in->myFptr, i, 
			(char*)cdata, (char*)tunit, (char*)typechar, &repeat, &scale, &zero,
			NULL, NULL, &status);
    /* Get column label */
    if (status==0) {
      /* Copy string fixing bug in cfitsio */
      ObitIOUVFITSFixBug(desc->ptype[i-1], cdata, UVLEN_KEYWORD);
    }
    /* Assure that data is a scalar float */
    if ((typechar[0]!='E') && (repeat!=1)) {
      retCode = OBIT_IO_SpecErr;
      Obit_log_error(err, OBIT_Error, 
		     "UV table column %d not correct type/dim ", i);
    }

    /* If this is the TIME/TIME1/DATE parameter get offset */
    if ((!strncmp (desc->ptype[i-1], "DATE",  4)) ||
	(!strncmp (desc->ptype[i-1], "TIME",  4)) ||
	(!strncmp (desc->ptype[i-1], "TIME1", 5))) {
      desc->JDObs = zero;
      /* Make sure ptype is "TIME1" */
      strcpy (desc->ptype[i-1], "TIME1");
    }
    /* U,V,W scaled to wavelengths */
    if ((!strncmp (desc->ptype[i-1], "UU-", 3)) ||
	(!strncmp (desc->ptype[i-1], "VV-", 3)) ||
	(!strncmp (desc->ptype[i-1], "WW-", 3))) {
      desc->ptype[i-1][3] = 'L';
    }
  } /* end loop reading random parameters */

  /* get visibility dimensionality array */
  viscol = (int)(desc->nrparm+1);
  fits_read_tdim (in->myFptr, viscol, UV_MAXDIM, &naxis,
		  inaxes, &status);
  desc->naxis = (olong)naxis; /* type mismatch for call argument */
  for (i=0; i<naxis; i++) desc->inaxes[i] = (olong)inaxes[i];

  /* Be sure last column is "VISIBILITIES" */
  fits_get_bcolparms (in->myFptr, viscol,
		      (char*)cdata, NULL, (char*)typechar, &repeat, &scale, &zero,
		      NULL, NULL, &status);
  /* fix string bug in cfitsio */
  ObitIOUVFITSFixBug(cdata, cdata, FLEN_CARD);
  if (strncmp (cdata, "VISIBILITIES", 12)) {
    retCode = OBIT_IO_SpecErr;
    Obit_log_error(err, OBIT_Error, 
		   "Last UV table column NOT VISIBILITIES ");
  }
  
  /* visibility data must be either TFLOAT with first axis dimensioned 3
     or TSHORT with first axis dimensioned 2. */
  if (!(((typechar[0]=='E') && (desc->inaxes[0]==3)) ||
      ((typechar[0]=='I') && (desc->inaxes[0]==2)))) {
    retCode = OBIT_IO_SpecErr;
    Obit_log_error(err, OBIT_Error, 
		   "Illegal visibility type/dimension %c %d ",
		   typechar[0], desc->inaxes[0]);
  }

  /* Data array definition */

  /* Axis labels */
  for (i=0; i<UV_MAXDIM; i++) strncpy (desc->ctype[i], "        ", UVLEN_KEYWORD-1);
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCTYP%d", i+1, viscol);
    fits_read_key_str (in->myFptr, (char*)keyword, (char*)cdata, (char*)commnt, &status);
    if (status==0) strncpy (desc->ctype[i], cdata, UVLEN_KEYWORD-1);
  }

  /* Axis increments */
  for (i=0; i<UV_MAXDIM; i++) desc->cdelt[i] = 0.0; /* defaults */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCDLT%d", i+1, viscol);
    ftemp = 0.0;
    fits_read_key_flt (in->myFptr, (char*)keyword, &ftemp, (char*)commnt, 
		       &status);
    if (status==KEY_NO_EXIST) status = 0;
    desc->cdelt[i] = (ofloat)ftemp;
  }
  
  /* Axis reference pixel */
  for (i=0; i<UV_MAXDIM; i++) desc->crpix[i] = 1.0;/* defaults */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCRPX%d", i+1, viscol);
    ftemp = 0.0;
    fits_read_key_flt (in->myFptr, (char*)keyword, &ftemp, (char*)commnt, 
		       &status);
    if (status==KEY_NO_EXIST) status = 0;
    desc->crpix[i] = (ofloat)ftemp;
  }

  /* Axis rotation */
  for (i=0; i<UV_MAXDIM; i++) desc->crota[i] = 0.0;/* defaults */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCROT%d", i+1, viscol);
    ftemp = 0.0;
    fits_read_key_flt (in->myFptr, (char*)keyword, &ftemp, (char*)commnt, 
		       &status);
    if (status==KEY_NO_EXIST) status = 0;
    desc->crota[i] = (ofloat)ftemp;
  }

  /* Axis coordinate value at reference pixel */
  for (i=0; i<UV_MAXDIM; i++) desc->crval[i] = 0.0;
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD, "%dCRVL%d", i+1, viscol);
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
  
  strncpy (desc->teles, "        ", UVLEN_VALUE-1);
  fits_read_key_str (in->myFptr, "TELESCOP", (char*)cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->teles, cdata, UVLEN_VALUE-1);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->instrument, "        ", UVLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "INSTRUME", (char*)cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->instrument, cdata, UVLEN_VALUE-1);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->observer, "        ", UVLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "OBSERVER", (char*)cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->observer, cdata, UVLEN_VALUE-1);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->object, "        ", UVLEN_VALUE-1);
  fits_read_key_str (in->myFptr, "OBJECT", (char*)cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->object, cdata, UVLEN_VALUE-1);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->bunit, "        ", UVLEN_VALUE-1);
  fits_read_key_str (in->myFptr, "BUNIT", (char*)cdata, (char*)commnt, &status);
  if (status==0) strncpy (desc->bunit, cdata, UVLEN_VALUE);
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->obsdat, "        ", UVLEN_VALUE-1);
  fits_read_key_str (in->myFptr, "DATE-OBS", (char*)cdata, (char*)commnt, &status);
  if (status==0)  strncpy (desc->obsdat, cdata, UVLEN_VALUE); 
  if (status==KEY_NO_EXIST) status = 0;

  today = ObitToday();
  strncpy (desc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  fits_read_key_str (in->myFptr, "DATE-MAP", (char*)cdata, (char*)commnt, &status);
  if (status==0)  strncpy (desc->date, cdata, UVLEN_VALUE-1);  desc->date[UVLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->origin, "        ", UVLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "ORIGIN", (char*)cdata, (char*)commnt, &status);
  if (status==0)  strncpy (desc->origin, cdata, UVLEN_VALUE-1); desc->origin[UVLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "OBSRA", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->obsra = (odouble)dtemp;

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "OBSDEC", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->obsdec = (odouble)dtemp;

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "ALTRVAL", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->altRef = (odouble)dtemp;

  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "ALTRPIX", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->altCrpix = (ofloat)ftemp;

  ltemp = 0;
  fits_read_key_lng (in->myFptr, "VELREF", &ltemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->VelDef = ltemp / 256;
  desc->VelReference = ltemp - 256*desc->VelDef;

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "RESTFREQ", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->restFreq = (odouble)dtemp;

  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "XSHIFT", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->xshift = (ofloat)ftemp;

  ftemp = 0.0;
  fits_read_key_flt (in->myFptr, "YSHIFT", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->yshift = (ofloat)ftemp;

  /* sort order - written in ancient AIPSish */
  ObitIOUVFITSSortRead (in, &status);

  /* Look for anything else and add it to the InfoList on desc */
  ObitIOUVKeysOtherRead(in, &status, err);
  if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);

  /* was there an error? */
  if (err->error) return retCode;
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR %d reading input FITS file header", status);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* enforce defaults */
  ObitUVSelDefault(desc, sel);

  return retCode;
} /* end ObitIOUVFITSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in Pointer to object with ObitUVDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitIOUVFITSWriteDescriptor (ObitIOUVFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar keyword[FLEN_KEYWORD], commnt[FLEN_COMMENT+1];
  gchar *ttype[50], tttype[50][14], *tform[50], *tunit[50];
  gchar *unitSeconds="SECONDS ", *unitDays="DAYS    ", *blank="        ";
  gchar *formRP="1E      ", *VISIBILITIES = "VISIBILITIES";
  gchar dtype, keyName[FLEN_KEYWORD+1], *keyNameP;
  gchar cdata[FLEN_CARD];
  int i, k, tfield, ndata, viscol, status = 0;
  long naxes[UV_MAXDIM], extver, ltemp, nrows;
  gboolean doFill;
  odouble dtemp;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  ObitInfoType keyType;
  gint32 nkey, dim[MAXINFOELEMDIM];
  union blobEquiv {
    gchar    s[21];
    double   d;
    float    f;
    gboolean b;
    olong    i;
  } blob;
 gchar *routine = "ObitIOUVFITSWriteDescriptor";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */

  retCode = OBIT_IO_OK; /* until proven otherwise */

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* enforce defaults */
  ObitUVSelDefault(desc, sel);

  /* First put AIPSish nonsense in main file header */
  retCode = WriteAIPSUVHeader (in, err);
  if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);

  /* Position to "AIPS UV" table version 1 if it exists */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, "AIPS UV", 1, &status);
  /* Create if it doesn't already exist */
  if (status == BAD_HDU_NUM) {
    status = 0;
    /* fill descriptive arrays */
    tfield = (int)(desc->nrparm+1);  /* number of columns */
    for (i=0; i<desc->nrparm; i++) {
      tform[i] = formRP;  /* Scalar float */
      strncpy (&tttype[i][0], desc->ptype[i], 9);
      tttype[i][8] = 0;
      ttype[i] = &tttype[i][0];
      ObitTrimTrail(tttype[i]);
      if ((i==desc->ilocu) || (i==desc->ilocv) || (i==desc->ilocw)) {
	/* Indicate time units */
	tttype[i][3] = '-';
	tunit[i] = unitSeconds;
      } else if (i==desc->iloct) {
	/* Indicate time units now JD */
	strncpy (&tttype[i][0], "DATE    ", 9); 
	tunit[i] = unitDays;
      } else
	tunit[i] = blank;
    } /* end converting random parameters to columns */

    /* Data column */
    tunit[tfield-1] = desc->bunit;
    strncpy (&tttype[tfield-1][0], VISIBILITIES, 13);
    ttype[tfield-1] = &tttype[tfield-1][0];
    /* how many data values? */
    if (sel->Compress) { /* Compressed */
      ndata = 2 * desc->ncorr;
      dtype = 'I';
    } else { /* uncompressed */
      ndata = 3 * desc->ncorr;
      dtype = 'E';
    }
    /* length and type*/
    g_snprintf (cdata, FLEN_CARD-1, "%d%c", ndata, dtype);
    tform[tfield-1] = cdata;
    
    /* create table */
    fits_create_tbl (in->myFptr, BINARY_TBL, (long)desc->nvis,
		     tfield, ttype, tform, tunit,
		     "AIPS UV", &status);
    if (status!=0) {
      Obit_cfitsio_error(err); 
      ObitFileErrMsg(err);     /* system error message*/
    }

    /* Add version number 1 */
    extver = 1;
    strncpy (commnt, "Table version number", FLEN_COMMENT);
    fits_update_key_lng (in->myFptr, "EXTVER", extver, (char*)commnt, &status);

    /* set visibility dimensionality array as TDIM */
    viscol = (int)(desc->nrparm+1);
    /* set dimensionality array */
    for (i=0; i<desc->naxis; i++) naxes[i] = (long)desc->inaxes[i];
    if (sel->Compress) naxes[0] = 2;
    else naxes[0] = 3;
    fits_write_tdim (in->myFptr, viscol, (long)desc->naxis, naxes, &status);
    if (status!=0) {
      Obit_cfitsio_error(err); 
      ObitFileErrMsg(err);     /* system error message*/
    }

    /* Random parameter scaling */
    /* U,V,W to time units */
    dtemp = 1.0 / (MAX (1.0, desc->freq));
    if (desc->ilocu>=0) {
      strncpy (commnt, "Scaling for U to time units", FLEN_COMMENT);
      g_snprintf (keyword, FLEN_KEYWORD-1, "TSCAL%d", desc->ilocu+1);
      fits_update_key_dbl (in->myFptr, (char*)keyword, dtemp, 12, (char*)commnt, 
			   &status);
   }
    if (desc->ilocv>=0) {
      strncpy (commnt, "Scaling for V to time units", FLEN_COMMENT);
      g_snprintf (keyword, FLEN_KEYWORD-1, "TSCAL%d", desc->ilocv+1);
      fits_update_key_dbl (in->myFptr, (char*)keyword, dtemp, 12, (char*)commnt, 
			   &status);
   }
    if (desc->ilocw>=0) {
      strncpy (commnt, "Scaling for W to time units", FLEN_COMMENT);
      g_snprintf (keyword, FLEN_KEYWORD-1, "TSCAL%d", desc->ilocw+1);
      fits_update_key_dbl (in->myFptr,(char*) keyword, dtemp, 12, (char*)commnt, 
			   &status);
   }
    /* Offset time to JD */
    if (desc->iloct>=0) {
      strncpy (commnt, "Offset of Date from JD", FLEN_COMMENT);
      g_snprintf (keyword, FLEN_KEYWORD-1, "TZERO%d", desc->iloct+1);
      fits_update_key_dbl (in->myFptr, (char*)keyword, (double)desc->JDObs, 12, 
			   (char*)commnt,  &status);
   }

  } /* end initialize new table */
  /* Data array definition */
    
  /* Axis labels */
  strncpy (commnt, "        ", 9);
  viscol = desc->nrparm+1;
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCTYP%d", i+1, viscol);
    strncpy (cdata, desc->ctype[i], 9); cdata[8] = 0;
    fits_update_key_str (in->myFptr, (char*)keyword, (char*)cdata, (char*)commnt, &status);
  }
  
  /* Axis increments */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCDLT%d", i+1, viscol);
    fits_update_key_flt (in->myFptr, (char*)keyword, (float)desc->cdelt[i], 6, (char*)commnt, 
			 &status);
  }
  
  /* Axis reference pixel */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCRPX%d", i+1, viscol);
    fits_update_key_flt (in->myFptr, (char*)keyword, (float)desc->crpix[i], 6, (char*)commnt, 
			 &status);
  }
  
  /* Axis rotation */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCROT%d", i+1, viscol);
    fits_update_key_flt (in->myFptr, (char*)keyword, (float)desc->crota[i], 6, (char*)commnt, 
			 &status);
  }
  
  /* Axis coordinate value at reference pixel */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "%dCRVL%d", i+1, viscol);
    fits_update_key_dbl (in->myFptr, (char*)keyword, (double)desc->crval[i], 12, (char*)commnt, 
			 &status);
  }
  
  /* descriptive information */
  /* Write keyword values */
  strncpy (commnt, "Name of object", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "OBJECT", (char*)desc->object,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Telescope used", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "TELESCOP", (char*)desc->teles,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Instrument used", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "INSTRUME", (char*)desc->instrument,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Observer/project", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "OBSERVER", (char*)desc->observer,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Date (yyyy-mm-dd) of observation", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "DATE-OBS", (char*)desc->obsdat, (char*)commnt, 
		       &status);
  strncpy (commnt, "Date (yyyy-mm-dd) created ", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "DATE-MAP", (char*)desc->date, (char*)commnt, 
		       &status);
  strncpy (commnt, "Software last writing file", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "ORIGIN", (char*)desc->origin, (char*)commnt, 
		       &status);
  strncpy (commnt, "Celestial coordiate equinox", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "EPOCH", (float)desc->epoch, 6,  (char*)commnt, 
		       &status);

  strncpy (commnt, "Visibility units", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "BUNIT", (char*)desc->bunit, (char*)commnt, &status);

  strncpy (commnt, "Observed Right Ascension", FLEN_COMMENT);
  fits_update_key_dbl (in->myFptr, "OBSRA", (double)desc->obsra, 12, (char*)commnt, 
		       &status);
  strncpy (commnt, "Observed declination ", FLEN_COMMENT);
  fits_update_key_dbl (in->myFptr, "OBSDEC", (double)desc->obsdec, 12, (char*)commnt, 
		       &status);

  if (desc->altRef != 0.0 ) {
    strncpy (commnt, "Alternate reference value", FLEN_COMMENT);
    fits_update_key_dbl (in->myFptr, "ALTRVAL", (double)desc->altRef, 12, (char*)commnt, 
		       &status);
  }

  if (desc->altCrpix != 0.0) {
    strncpy (commnt, "Alternate reference pixel", FLEN_COMMENT);
    fits_update_key_flt (in->myFptr, "ALTRPIX", (float)desc->altCrpix, 6, (char*)commnt, 
		       &status);
  }

  ltemp = (long)(desc->VelReference + 256*desc->VelDef);
  if (ltemp!=0) {
     strncpy (commnt, ">256 radio, 1 LSR, 2 Hel, 3 Obs", FLEN_COMMENT);
    fits_update_key_lng (in->myFptr, "VELREF", ltemp, (char*)commnt, 
		       &status);
  }

  if (desc->restFreq != 0.0) {
    strncpy (commnt, "Line rest frequency (Hz)", FLEN_COMMENT);
    fits_update_key_dbl (in->myFptr, "RESTFREQ", (double)desc->restFreq, 12, (char*)commnt, 
		       &status);
  }

  if (desc->xshift != 0.0) {
     strncpy (commnt, "Net shift of Phase center in x", FLEN_COMMENT);
     fits_update_key_flt (in->myFptr, "XSHIFT", (float)desc->xshift, 6, (char*)commnt, 
			  &status);
  }

  if (desc->yshift != 0.0) {
     strncpy (commnt, "Net shift of Phase center in y", FLEN_COMMENT);
     fits_update_key_flt (in->myFptr, "YSHIFT", (float)desc->yshift, 6, (char*)commnt, 
			  &status);
  }

  /* sort order - to be done in ancient AIPSish */
  ObitIOUVFITSSortWrite (in, err);
  if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);

 /* Number of vis NAXIS2 - truncate if too many */
  nrows = (long)desc->nvis;
  fits_get_num_rows (in->myFptr, &nrows, &status);
  if (nrows > desc->nvis) {  /* Truncate */
    nrows = (long)(nrows - desc->nvis);
    fits_delete_rows (in->myFptr, (long)(desc->nvis+1), nrows, &status);
  }

  /* Write other keywords from descriptor */
  if (desc->info) nkey = desc->info->number; /* How many keywords? */
  else nkey = 0;
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
    Obit_log_error(err, OBIT_Error, "ERROR %d updating FITS file header",status);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOUVFITSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVFITSFlush (ObitIOUVFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* cfitsio does the buffer flushing on close */

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOUVFITSFlush */

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
ObitIOUVFITSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOUVFITS *in, ObitInfoList *info, 
			     ObitErr *err)
{
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  errno = 0;  /* reset any system error */
  Obit_return_if_fail(((in->myStatus==OBIT_Modified) ||
		       (in->myStatus==OBIT_Active)), 
		      err,
		      "Cannot define buffer, I/O not currently active");

  /* get size */
  *size = ObitUVSelBufferSize(in->myDesc, in->mySel);

  /* (re)allocate */
  if (*data) *data = ObitMemRealloc (*data, (*size)*sizeof(ofloat));
  else *data = ObitMemAlloc0Name((*size)*sizeof(ofloat), "UVBuffer");

} /* end ObitIOUVFITSCreateBuffer */

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
newObitIOUVFITSTable (ObitIOUVFITS *in, ObitIOAccess access, 
			 gchar *tabType, olong *tabVer, ObitErr *err)
{
  ObitTable *out;
  olong version;
  gboolean gotIt;
  gchar *outName, tabName[51];
  gchar *routine = "newObitIOUVFITSTable";

  /* error check */
  if (err->error) return NULL;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  g_assert(tabType!=NULL);
  g_assert(tabVer!=NULL);
  errno = 0;  /* reset any system error */

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
} /* end newObitIOUVFITSTable */

/**
 * Update any disk resident structures about the current tables.
 * Nothing is needed for FITS files.
 * \param in   Pointer to object to be updated.
 * \param info ObitInfoList of parent object (not used here).
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOUVFITSUpdateTables (ObitIOUVFITS *in, ObitInfoList *info,
				     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  return retCode;
} /* end ObitIOUVFITSUpdateTables */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOUVFITSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOUVFITSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOUVFITSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOUVFITSClassInfoDefFn (gpointer inClass)
{
  ObitIOUVFITSClassInfo *theClass = (ObitIOUVFITSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOUVFITSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOUVFITSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOUVFITSGetClass;
  theClass->newObit    = NULL;
  theClass->newObitIO  = (newObitIOFP)newObitIOUVFITS;
  theClass->ObitIOSame = (ObitIOSameFP)ObitIOUVFITSSame;
  theClass->ObitIORename  = (ObitIORenameFP)ObitIOUVFITSRename;
  theClass->ObitIOZap  = (ObitIOZapFP)ObitIOUVFITSZap;
  theClass->ObitCopy   = (ObitCopyFP)ObitIOUVFITSCopy;
  theClass->ObitClone  = NULL;
  theClass->ObitClear  = (ObitClearFP)ObitIOUVFITSClear;
  theClass->ObitInit   = (ObitInitFP)ObitIOUVFITSInit;
  theClass->ObitIOOpen = (ObitIOOpenFP)ObitIOUVFITSOpen;
  theClass->ObitIOClose= (ObitIOCloseFP)ObitIOUVFITSClose;
  theClass->ObitIOSet  = (ObitIOSetFP)ObitIOUVFITSSet;
  theClass->ObitIORead = (ObitIOReadFP)ObitIOUVFITSRead;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOUVFITSReadSelect;
  theClass->ObitIOWrite   = 
    (ObitIOWriteFP)ObitIOUVFITSWrite;
  theClass->ObitIOFlush   = 
    (ObitIOFlushFP)ObitIOUVFITSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOUVFITSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOUVFITSWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOUVFITSCreateBuffer;
  theClass->ObitIOFreeBuffer   = 
    (ObitIOFreeBufferFP)ObitIOFreeBuffer;
  theClass->newObitIOTable = 
    (newObitIOTableFP)newObitIOUVFITSTable; 
  theClass->ObitIOUpdateTables   = 
    (ObitIOUpdateTablesFP)ObitIOUVFITSUpdateTables;

} /* end ObitIOUVFITSClassDefFn */

/*--------------- Private functions --------------------------*/

/**
 * Creates empty member objects.
 * for each parent class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOUVFITSInit  (gpointer inn)
{
  const ObitClassInfo *ParentClass;
  ObitIOUVFITS *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && (ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->FileName     = NULL;
  in->decompVis    = NULL;
  in->compBuff     = NULL;
  in->compBuffSize = 0;

} /* end ObitIOUVFITSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOUVFITSClear (gpointer inn)
{
  ObitIOUVFITS *in = inn;
  const ObitClassInfo *ParentClass;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->myStatus==OBIT_Active) ||(in->myStatus==OBIT_Modified)) {
    err = newObitErr();
    ObitIOUVFITSClose (in, err); 
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  if (in->FileName) g_free(in->FileName);
  in->FileName = NULL;
  if (in->compBuff)  in->compBuff  = ObitMemFree (in->compBuff); 
  if (in->decompVis) in->decompVis = ObitMemFree (in->decompVis);

 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && (ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOUVFITSClear */

/**
 * Uses selector member to decide which visibilities to read next.
 * Leaves values in myDesc->firstVis and mySel->numVisRead
 * \param  in   Pointer to the object.
 * \param  err  ObitErr for reporting errors.
 */
static gboolean ObitIOUVFITSNext (ObitIOUVFITS *in, ObitErr *err)
{
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gboolean done = FALSE;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitIsA(in, &myClassInfo));
 
  sel  = in->mySel;  /* selector pointer */
  desc = in->myDesc; /* UV descriptor pointer */
 
  /* Let Selector decide */
  done = ObitUVSelNext (sel, desc, err);

  return done;  
} /* end ObitIOUVFITSNext */

/**
 * Compresses UV into scaled shorts
 * Compressed data stores a common weigh and scaling factors as 
 * random parameters and the real and imaginary parts as scaled shorts.
 * For flagged data, the values of the short pair is -32767
 * Input data assumed inhost order and will be swapped to FITS order 
 * (bigendian) if different.
 * \param  ncorr  Number of weighted complex numbers
 * \param  visin  Expanded visibility array
 * \param  wtscl  (out) Weight and Scale needed to uncompress.
 * \param  visout (out) Compressed visibility array.
 */
static void 
ObitIOUVFITSCompress (gint ncorr, const ofloat *visin, ofloat *wtscl, 
		      ofloat *visout)
{
  olong i;
  ofloat maxwt, maxvis, scl;
  gshort *packed = (gshort*)visout;
  union sequiv inu1, inu2, outu1, outu2;

  /* error tests */
  if (ncorr <1) return;
  g_assert (visin != NULL);
  g_assert (wtscl != NULL);
  g_assert (visout != NULL);

  /* find maximum weight and visibility component */
  maxwt = maxvis = 0.0;
  for (i=0; i<ncorr; i++) {
    if (visin[i*3+2] > 0.0) { /* Valid? */
      maxvis = MAX (maxvis, fabs(visin[i*3]));
      maxvis = MAX (maxvis, fabs(visin[i*3+1]));
      maxwt  = MAX (maxwt, visin[i*3+2]);
    }
  }

  /* output weighting and scaling */
  wtscl[0] = maxwt;
  wtscl[1] = maxvis / 32760.;
  scl = 1.0;
  if (wtscl[1] > 1.0e-10) scl = 1.0 / wtscl[1];

  /* loop over visibilities packing them in. */
  for (i=0; i<ncorr; i++) { 
    /* blanked or unblanked */
    if (visin[i*3+2] > 0.0) { /* OK - round values */
      if (visin[i*3] > 0.0)
	inu1.full = (gshort)(scl*visin[i*3] + 0.5);
      else
	inu1.full = (gshort)(scl*visin[i*3] - 0.5);
      if (visin[i*3+1] > 0.0)
	inu2.full = (gshort)(scl*visin[i*3+1] + 0.5);
      else
	inu2.full = (gshort)(scl*visin[i*3+1] - 0.5);
    } else { /* flag */
      inu1.full = -32767;
      inu2.full = -32767;
    }

    /* copy to output, byte swap as needed */
#if G_BYTE_ORDER==G_BIG_ENDIAN  /* no byte swap needed */
    packed[i*2]   = outu1.full;
    packed[i*2+1] = outu2.full ;
    
#elif G_BYTE_ORDER==G_LITTLE_ENDIAN /* byte swap */
    outu1.parts[0] = inu1.parts[1]; 
    outu1.parts[1] = inu1.parts[0]; 
    outu2.parts[0] = inu2.parts[1]; 
    outu2.parts[1] = inu2.parts[0]; 
    packed[i*2]   = outu1.full;
    packed[i*2+1] = outu2.full ;
    
#else /* unknown */
    g_error("ObitIOUVFITSH2F: Unsupported host byte order");
#endif
  } /* end loop over correlations */
} /* end ObitIOUVFITSCompress */

/**
 * Uncompresses UV from scaled shorts.
 * Compressed data stores a common weigh and scaling factors as 
 * random parameters and the real and imaginary parts as scaled shorts.
 * Values of -32767 in both of a pair of shorts indicate a flagged value.
 * \param  ncorr  Number of weighted complex numbers
 * \param  visin  Compressed visibility array.
 * \param  wtscl  Weight and Scale needed to uncompress.
 * \param  visout (out) Expanded visibility array.
 */
static void 
ObitIOUVFITSUncompress (gint ncorr, const ofloat *visin, 
			const ofloat *wtscl, ofloat *visout)
{
  olong i;
  ofloat wt, scl;
  gshort *packed = (gshort*)visin;

  /* error tests */
  if (ncorr <1) return;
  g_assert (visin != NULL);
  g_assert (wtscl != NULL);
  g_assert (visout != NULL);

  /* weighting and scaling */
  wt  = wtscl[0];
  scl = wtscl[1];

  /* loop over visibilities */
  for (i=0; i<ncorr; i++) { 

    /* blanked or unblanked */
    if ( packed[i*2] == -32767) { /* Flagged */
      visout[i*3]   = 0.0;
      visout[i*3+1] = 0.0;
      visout[i*3+2] = 0.0;
    } else { /* OK */
      visout[i*3]   = scl * packed[i*2];
      visout[i*3+1] = scl * packed[i*2+1];
      visout[i*3+2] = wt;
    }
  }
} /* end ObitIOUVFITSUncompress */

/**
 * Swaps byte order in floats if host order differs from FITS.
 * Input data assumed in host order and will be swapped to FITS order 
 * (bigendian) if different.
 * Will work in place.
 * \param  n    Number of floats
 * \param  in   Array of input floats in host byte order.
 * \param  out  Array of output floats in FITS order
 */
static void ObitIOUVFITSfH2F (gint n, ofloat *in, ofloat *out)
{
  olong i;
  union fequiv inu, outu;

#if G_BYTE_ORDER==G_BIG_ENDIAN  /* no byte swap needed */
  /* if the input and output point to the same place - just return */
  if (in==out) return;
  for (i=0; i<n; i++) out[i] = in[i]; /* simple copy */

#elif G_BYTE_ORDER==G_LITTLE_ENDIAN   /* byte swap */
  for (i=0; i<n; i++) {
    inu.full = in[i];
    outu.parts[0] = inu.parts[3]; 
    outu.parts[1] = inu.parts[2]; 
    outu.parts[2] = inu.parts[1]; 
    outu.parts[3] = inu.parts[0]; 
    out[i] = outu.full;
  }

#else /* unknown */
  g_error("ObitIOUVFITSfH2F: Unsupported host byte order");
#endif
} /* end ObitIOUVFITSfH2F */

/**
 * Swaps byte order in floats if host order differs from FITS.
 * Input data assumed in FITS byte order (bigendian) and will be
 * swapped to host order if different.
 * Will work in place.
 * \param  n    Number of floats
 * \param  in   Array of input floats in FITS order
 * \param  out  Array of output floats in host byte order.
 */
static void ObitIOUVFITSfF2H (gint n, ofloat *in, ofloat *out)
{
  olong i;
  union fequiv inu, outu;

#if G_BYTE_ORDER==G_BIG_ENDIAN  /* no byte swap needed */
  /* if the input and output point to the same place - just return */
  if (in==out) return;
  for (i=0; i<n; i++) out[i] = in[i]; /* simple copy */

#elif G_BYTE_ORDER==G_LITTLE_ENDIAN   /* byte swap */
  for (i=0; i<n; i++) {
    inu.full = in[i];
    outu.parts[0] = inu.parts[3]; 
    outu.parts[1] = inu.parts[2]; 
    outu.parts[2] = inu.parts[1]; 
    outu.parts[3] = inu.parts[0]; 
    out[i] = outu.full;
  }

#else /* unknown */
  g_error("ObitIOUVFITSfF2H: Unsupported host byte order");
#endif
} /* end ObitIOUVFITSfF2H */

/**
 * Swaps byte order in shorts if host order differs from FITS.
 * Input data assumed in host order and will be swapped to FITS order 
 * (bigendian) if different.
 * Will work in place.
 * \param  n    Number of shorts
 * \param  in   Array of input shorts in host byte order.
 * \param  out  Array of output shorts in FITS order
 */
static void ObitIOUVFITSsH2F (gint n, gshort *in, gshort *out)
{
  olong i;
  union sequiv inu, outu;

#if G_BYTE_ORDER==G_BIG_ENDIAN  /* no byte swap needed */
  /* if the input and output point to the same place - just return */
  if (in==out) return;
  for (i=0; i<n; i++) out[i] = in[i]; /* simple copy */

#elif G_BYTE_ORDER==G_LITTLE_ENDIAN   /* byte swap */
  for (i=0; i<n; i++) {
    inu.full = in[i];
    outu.parts[0] = inu.parts[1]; 
    outu.parts[1] = inu.parts[0]; 
    out[i] = outu.full;
  }

#else /* unknown */
  g_error("ObitIOUVFITSsH2F: Unsupported host byte order");
#endif
} /* end ObitIOUVFITSH2F */

/**
 * Swaps byte order in shorts if host order differs from FITS.
 * Input data assumed in FITS byte order (bigendian) and will be
 * swapped to host order if different.
 * Will work in place.
 * \param  n    Number of shorts
 * \param  in   Array of input shorts in FITS order
 * \param  out  Array of output shorts in host byte order.
 */
static void ObitIOUVFITSsF2H (gint n, gshort *in, gshort *out)
{
  olong i;
  union sequiv inu, outu;

#if G_BYTE_ORDER==G_BIG_ENDIAN  /* no byte swap needed */
  /* if the input and output point to the same place - just return */
  if (in==out) return;
  for (i=0; i<n; i++) out[i] = in[i]; /* simple copy */

#elif G_BYTE_ORDER==G_LITTLE_ENDIAN   /* byte swap */
  for (i=0; i<n; i++) {
    inu.full = in[i];
    outu.parts[0] = inu.parts[1]; 
    outu.parts[1] = inu.parts[0]; 
    out[i] = outu.full;
  }

#else /* unknown */
  g_error("ObitIOUVFITSsF2H: Unsupported host byte order");
#endif
} /* end ObitIOUVFITSsF2H */

/**
 * Look for rational keyword (SORTORD) for the Sort Order and
 * failing this, look in AIPS history keyword.
 * Descriptor value isort.
 * \param in      Pointer to ObitIOUVFITS.
 * \param status (Output) cfitsio status.
 * \return return code, 0=> OK
 */
static void  ObitIOUVFITSSortRead(ObitIOUVFITS *in, olong *lstatus)
{
  gchar commnt[FLEN_COMMENT], card[FLEN_COMMENT], cdata[FLEN_COMMENT];
  int k, keys, morekeys, status=(int)*lstatus;
  ObitUVDesc *desc;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  desc = in->myDesc; /* set descriptor */

  /* init */
  desc->isort[0] = ' '; desc->isort[1] = ' '; desc->isort[2] = 0; 

  /* Attempt rational keywords */
  fits_read_key_str (in->myFptr, "SORTORD", (char*)cdata, (char*)commnt, &status);
  if (status==0) {
    desc->isort[0]=cdata[0]; 
    desc->isort[1]=cdata[1]; 
    return;
  }
  if (status==KEY_NO_EXIST) status = 0;

  /* Oh Well, parse all the header cards looking for: 
          1         2         3  
012345678901234567890123456789012
HISTORY AIPS   SORT ORDER = 'TB'
  */

  /* how many keywords to look at? */
  fits_get_hdrspace (in->myFptr, &keys, &morekeys, &status);
  for (k=1; k<=keys; k++) {
    fits_read_record (in->myFptr, k, (char*)card, &status);
    if (status==0) {
      if (!strncmp ("HISTORY AIPS   SORT ORDER", card, 25)) {
	/* Parse card */
	desc->isort[0]=card[29]; 
	desc->isort[1]=card[30]; 
      }
    }
  } /* end loop over header cards */
  *lstatus = (olong)status;
} /* end ObitIOUVFITSSortRead */

/**
 * Write both rational keyword and HISTORY AIPS card
 * Descriptor value isort
 * \param in Pointer to ObitIOUVFITS.
 * \param err    ObitErr stack.
 * \return return code, 0=> OK
 */
static void  ObitIOUVFITSSortWrite (ObitIOUVFITS *in, ObitErr *err)
{
  gchar so[3], commnt[FLEN_COMMENT+1];
  int status = 0;
  ObitUVDesc *desc;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (err->error) return;

  desc = in->myDesc; /* set descriptor */
  /* Don't bother if not given */
  if (((desc->isort[0]==' ') && (desc->isort[1]==' ')) ||
      (desc->isort[0]==0)) return;
 
  /* Rational keyword */
  strncpy (commnt, "Sort Order code (in AIPSish)", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "SORTORD", (char*)desc->isort, (char*)commnt, 
		       &status);
  if (status!=0) { /* error */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR %d Writing Sort Order to file %s", 
		   status, in->FileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    return;
  }

  /* Now Hide it where AIPS can find it */
  /*  g_snprintf (card, FLEN_COMMENT,
      "AIPS   SORT ORDER = '%c%c'",
      desc->isort[0], desc->isort[1]);
      fits_write_history (in->myFptr, card, status);*/
  so[0] = desc->isort[0]; so[1] =  desc->isort[1]; so[2] = 0;
  ObitFileFITSWriteHisKeyStr (in->myFptr, "AIPS    SORT ORDER", so, 
			      commnt, err);
  
} /* end ObitIOUVFITSSortWrite */

/**
 * Look for additional descriptive keywords, any that are not 
 * on the exclusion list are copied to the descriptor InfoList.
 * \param in       Pointer to ObitIOUVFITS.
 * \param lstatus (Output) cfitsio status.
 * \param err      ObitErr stack.
 * \return return code, 0=> OK
 */
void  ObitIOUVKeysOtherRead(ObitIOUVFITS *in, olong *lstatus, 
			       ObitErr *err)
{
  gchar keywrd[FLEN_KEYWORD], value[FLEN_VALUE], commnt[FLEN_COMMENT+1];
  gchar *first, *last, *anF, *aT, dtype, svalue[FLEN_VALUE];
  int i, j, k, l, keys, morekeys, status=(int)*lstatus;
  olong ivalue;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  double dvalue;
  ObitUVDesc *desc;
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
  olong number, *len=NULL;
  gboolean bvalue, bad=FALSE;
  gchar *routine = "ObitIOUVKeysOtherRead";

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete old InfoList and restart */
  ((ObitUVDesc*)in->myDesc)->info = ObitInfoListUnref (((ObitUVDesc*)in->myDesc)->info);
  ((ObitUVDesc*)in->myDesc)->info = (gpointer)newObitInfoList ();
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
    fits_read_keyn (in->myFptr, k, (char*)keywrd, (char*)value, (char*)commnt, &status);
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
	fits_get_keytype (value, (char*)&dtype, &status);
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
	  ObitInfoListPut(desc->info, keywrd, OBIT_long, dim, 
			  (gconstpointer)&ivalue, err);
	  break;
	case 'F':  /* Float - use double */
	  /* AIPS uses 'D' for double exponent */
	  for (l=0; l<strlen(value); l++) if (value[l]=='D') value[l]='e';
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

  if (len) g_free(len);
  *lstatus = (olong)status;
} /* end ObitIOUVKeysOtherRead */

/**
 * Work around for cfitsio bug, trailing blanks in keywords are dropped.
 * This routine blank fills out and null terminates at out[maxn-1].
 * \param out    Output string 
 * \param in     Input string from cfitsio keyword read
 * \param maxn   length of out
 */
void ObitIOUVFITSFixBug (gchar *out, gchar *in, olong maxn)
{
  olong i, len;

  len = strlen(in);
  for (i=0; i<len; i++)    out[i] = in[i];
  for (i=len; i<maxn; i++) out[i] = ' ';
  out[maxn-1] = 0;
  
} /* end ObitIOUVFITSFixBug */

/**
 * Write Descriptor information for AIPS in main file header.
 * This stuff is needed because of the hack made in AIPS to 
 * read/write uv data as binary tables rather than random groups.
 * \param  in Pointer to object with ObitUVDesc to be written.
 * \param  err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
static ObitIOCode WriteAIPSUVHeader (ObitIOUVFITS *in, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  gchar keyword[FLEN_KEYWORD], strtemp[9];
  int i, hdutype, status = 0;
  olong velref;
  float zero, scale;
  ObitUVDesc* desc;
  gchar *routine = "WriteAIPSUVHeader";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);

  retCode = OBIT_IO_OK; /* until proven otherwise */

  desc = in->myDesc; /* descriptor pointer */

  /* Just return if descriptor not filled in */
  if (desc->nrparm<=0) return retCode;

  /* Position to main file header */
  fits_movabs_hdu (in->myFptr, 1, &hdutype, &status);
  if (status!=0) {
    Obit_cfitsio_error(err); 
    ObitFileErrMsg(err);     /* system error message*/
  }

  /* Amount of data */
  ObitFileFITSWriteHisKeyLng (in->myFptr, "AIPS    GCOUNT  ", (long)desc->nvis,
			      NULL, err);
  ObitFileFITSWriteHisKeyLng (in->myFptr, "AIPS    PCOUNT  ", (long)desc->nrparm,
			      NULL, err);

  /* Sort order */
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    SORT ORDER ");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->isort, NULL, err);

  /* Loop over random parameters */
  for (i=0; i<desc->nrparm; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    PTYPE%d ", i+1);
    /* Indicate time units for u,v,w */
    strncpy (strtemp, desc->ptype[i], 9); strtemp[8] = 0;
    if ((desc->ilocu==i) || (desc->ilocv==i) || (desc->ilocw==i)) strtemp[3] = '-';
    if (desc->iloct==i) strncpy (strtemp, "DATE    ", 9); strtemp[8] = 0;
    ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, strtemp, NULL, err);
    scale = 1.0 / (MAX (1.0, desc->freq));
    if ((desc->ilocu==i) || (desc->ilocv==i) || (desc->ilocw==i)) {
      g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    PSCAL%d ", i+1);
      ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, scale, NULL, err);
    }
    if (desc->iloct==i) {
      zero = (float)desc->JDObs;
      g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    PZERO%d ", i+1);
      ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, zero, NULL, err);
    }
  } /* end loop over random parameters */

  /* Regular axes */
  ObitFileFITSWriteHisKeyLng (in->myFptr, "AIPS    NAXIS   ", (long)desc->naxis,
			      NULL, err);
  /* Loop over regular axes */
  for (i=0; i<desc->naxis; i++) {
    g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    NAXIS%d ", i+1);
    ObitFileFITSWriteHisKeyLng (in->myFptr, (char*)keyword, (long)desc->inaxes[i], NULL, err);
    g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    CTYPE%d ", i+1);
    ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->ctype[i], NULL, err);
    g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    CRVAL%d ", i+1);
    ObitFileFITSWriteHisKeyDbl (in->myFptr, (char*)keyword, (double)desc->crval[i], NULL, err);
    g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    CDELT%d ", i+1);
    ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, (float)desc->cdelt[i], NULL, err);
    g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    CRPIX%d ", i+1);
    ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, (float)desc->crpix[i], NULL, err);
    g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    CROTA%d ", i+1);
    ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, (float)desc->crota[i], NULL, err);
  } /* end loop over Regular axes */

  /* cats and dogs */

  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    OBJECT  ");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->object, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    TELESCOP");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->teles, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    INSTRUME");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->instrument, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    OBSERVER");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->observer, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    DATE-OBS");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->obsdat, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    DATE-MAP");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->date, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    BUNIT   ");
  ObitFileFITSWriteHisKeyStr (in->myFptr, (char*)keyword, (char*)desc->bunit, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    BSCALE  ");
  ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, (float)1.0, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    BZERO   ");
  ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, (float)0.0, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    ALTRPIX ");
  ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, (float)desc->altRef, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    OBSRA   ");
  ObitFileFITSWriteHisKeyDbl (in->myFptr, (char*)keyword, (double)desc->obsra, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    OBSDEC  ");
  ObitFileFITSWriteHisKeyDbl (in->myFptr, (char*)keyword, (double)desc->obsdec, NULL, err);
  
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    BLANK   ");
  ObitFileFITSWriteHisKeyFlt (in->myFptr, (char*)keyword, (float)-1.0, NULL, err);
  
  velref = desc->VelReference + 256*desc->VelDef;
  g_snprintf (keyword, FLEN_KEYWORD-1, "AIPS    VELREF  ");
  ObitFileFITSWriteHisKeyLng (in->myFptr, (char*)keyword, (long)velref, NULL, err);
  
  /* was there an error? */
  if (err->error) Obit_traceback_val (err, routine, in->name, OBIT_IO_WriteErr);
  
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s: ERROR %d updating AIPSish FITS file header", 
		   routine, status);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_WriteErr;
  }
  
  return OBIT_IO_OK;
} /* end  WriteAIPSUVHeader */
