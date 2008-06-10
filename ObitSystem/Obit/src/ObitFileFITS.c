/* $Id: ObitFileFITS.c,v 1.9 2008/03/03 21:14:00 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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

#include "ObitFileFITS.h"
#include "ObitFITS.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFileFITS.c
 * ObitFileFITS class function definitions.
 *
 * This implementation uses cfitsio.
 * This class allows generic access to FITS files.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitFileFITS";

/** Function to obtain parent ClassInfo - ObitFile */
static ObitGetClassFP ObitParentGetClass = ObitFileGetClass;

/**
 * ClassInfo global structure ObitFileFITSClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFileFITSClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFileFITSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFileFITSClear (gpointer in);

/** Private: Find header card. */
static olong FindHeaderCard(fitsfile *inFptr, gchar *Name, gchar *card,
			   ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitFileFITSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitFileFITS* newObitFileFITS (gchar* name)
{
  ObitFileFITS* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFileFITSClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFileFITS));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFileFITSInit((gpointer)out);

  return out;
} /* end newObitFileFITS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitFileFITSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFileFITSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * The contents of any files are not modified.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitFileFITS* ObitFileFITSCopy  (ObitFileFITS *in, ObitFileFITS *out, ObitErr *err)
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
    out = newObitFileFITS(outName);
    g_free(outName);
  }

  /* Copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->myFptr    = in->myFptr;

  return out;
} /* end ObitFileFITSCopy */

/**
 * Determine if a given file name exists.
 * \param fileName  Name of file to test.
 * \param err       ObitErr for reporting errors.
 * \return TRUE if exists, else FALSE.
 */
gboolean ObitFileFITSExist (gchar *fileName, ObitErr *err)
{

  /* Use parent class function */
  return ObitFileExist (fileName, err);
} /* end ObitFileFITSExist */

/**
 * Delete the file.
 * Objust must have been fully instantiated (opened) first
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return NULL as value of pointer to in, on failure returns in.
 */
ObitFileFITS* ObitFileFITSZap (ObitFileFITS *in, ObitErr *err)
{
  int status = 0;
  gchar *routine = "ObitFileFITSZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return in;

  /* Zap file */
  fits_delete_file (in->myFptr, &status);
    if (status) {             /* it went wrong */
      Obit_log_error(err, OBIT_Error, 
		     "%s: ERROR deleting file %s", routine, in->name);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      return in;
    }
    in->myFptr = NULL;
 
  /* delete the rest of the structure */
  in = ObitFileFITSUnref(in); 
  return in;
} /* end ObitFileFITSZap */

/**
 * Delete the current HDU.
 * Object must have been fully instantiated (opened) first
 * Refuses for the first (main) HDU
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitFileFITSZapHDU (ObitFileFITS *in, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int hdunum, hdutype, status = 0;
  gchar *routine = "ObitFileFITSZapHDU";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;

  /* Make sure not main HDU */
  fits_get_hdu_num (in->myFptr, &hdunum);

  /* Check */
  if (hdunum==1) { /* NOPE */
    Obit_log_error(err, OBIT_Error, "%s: I WILL NOT delete the main HDU in %s", 
		   routine, in->name);
    return retCode;
  }

  /* Zap current HDU */
  fits_delete_hdu (in->myFptr, &hdutype, &status);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, "%s: ERROR deleting HDU %d %s", 
		   routine, hdunum, in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return retCode;
  }
 
  return OBIT_IO_OK;
} /* end ObitFileFITSZap */

/**
 * Initialize structures and open file.
 * The file will be positioned at the beginning.
 * \param in        Pointer to object to be opened.
 * \param fileName  Name of file to open. Prepend '!' to overwite
 * \param type      Obit FITS disk number
 * \param access    access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param err       ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSOpen (ObitFileFITS *in, gchar *fileName, olong disk,
		  ObitIOAccess access, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar *tempStr;
  gboolean exist;
  gchar *routine = "ObitFileFITSOpen";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;
  g_assert (fileName!=NULL);

  /* Save call arguments */
  in->access    = access;

  in->status = OBIT_ErrorExist; /* in case something goes wrong */

  /* form file name for file */
  in->disk = disk;
  if (in->fileName) g_free(in->fileName); /* release old */
  in->fileName = ObitFITSFilename (in->disk, fileName, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
   
  /* open file by access type */
  /*------------------------ Read Only ---------------------------------*/
  if (access == OBIT_IO_ReadOnly) {
    /* must strip any leading "!" for read/write */
    tempStr = in->fileName;
    if (in->fileName[0]=='!') tempStr = in->fileName+1;
    if ( fits_open_file(&(in->myFptr), tempStr, READONLY, &status) ) {
      Obit_log_error(err, OBIT_Error, 
		  "ERROR opening input FITS file %s", in->fileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */

      retCode = OBIT_IO_OpenErr;
      return retCode;
    }

  /*------------------------ Read/Write ---------------------------------*/
  } else if (access == OBIT_IO_ReadWrite) {
    /* must strip any leading "!" for read/write */
    tempStr = in->fileName;
    if (in->fileName[0]=='!') tempStr = in->fileName+1;
    if ( fits_open_file(&(in->myFptr), tempStr, READWRITE, &status) ) {
      Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", 
		     in->fileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      retCode = OBIT_IO_OpenErr;
      return retCode;
    }

  /*------------------------ Write Only ---------------------------------*/
  } else if (access == OBIT_IO_WriteOnly) {
    /* Initialize output file */
    /* Output file may already exist test open */
    /* must strip any leading "!" for read/write */
     tempStr = in->fileName;
    if (in->fileName[0]=='!') tempStr = in->fileName+1;

    /* Does it already exist? If not must create rathern than open */
    exist = ObitFileFITSExist (in->fileName, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    if (exist) { /* Open */
      if (fits_open_file(&(in->myFptr), in->fileName, READWRITE, &status)) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS file %s", in->fileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	return OBIT_IO_OpenErr;
      }
      
    } else { /* must create */
      if (fits_create_file(&(in->myFptr), tempStr, &status)) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS file %s", in->fileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	return OBIT_IO_OpenErr;
      }
    } /* end open or create */

  } else {  /* Unknown access */
    /* should never get here */
    g_assert_not_reached(); 
  }

  in->status = OBIT_Active;  /* seems to be OK */
  in->LastKeyword = 0;
  return OBIT_IO_OK;
} /* end ObitFileFITSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitFileFITSClose (ObitFileFITS *in, ObitErr *err)
{
  int status = 0;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitFileFITSClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* don't bother if it's not open (shutdown even if in error) */
  if ((in->status!=OBIT_Modified) && (in->status!=OBIT_Active) 
      && (in->status!=OBIT_ErrorExist)) 
    return OBIT_IO_OK;

  /* close file */
  fits_close_file(in->myFptr, &status);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing FITS file %s", 
		   routine, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_CloseErr;
  }

  /* mark file structure */
  in->myFptr = NULL;
  in->status = OBIT_Inactive;

  return OBIT_IO_OK;
} /* end ObitFileFITSClose */

/**
 * Position to given extention number
 * \param in      Pointer to object to be positioned
 * \param hdunum  Desired (1-rel) HDU
 * \param hdutype cfitsio code for type, NULL => not wanted
 *                has values (IMAGE_HDU, ASCII_TBL, BINARY_TBL)
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitFileFITSPosNum (ObitFileFITS *in, olong hdunum, olong *hdutype, 
			    ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0, ffhdunum=(int)hdunum, ffhdutype;
  gchar *routine = "ObitFileFITSPosNum";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* reposition */
  fits_movabs_hdu(in->myFptr, ffhdunum, &ffhdutype, &status);
  *hdutype = (olong)ffhdutype;
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, "%s: ERROR positioning FITS file %s", 
		   routine, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_CloseErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSPosNum */

/**
 * Position to given extention nam2
 * \param in      Pointer to object to be positioned
 * \param hdutype (IMAGE_HDU, ASCII_TBL, BINARY_TBL, ANY_HDU)
 * \param extname  extension (table) name
 * \param extver   extension (table) version, 0=>next
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitFileFITSPosName (ObitFileFITS *in, olong hdutype, gchar *extname, 
			    olong extver, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0, ffextver=(int)extver, ffhdutype=(int)hdutype;
  gchar *routine = "ObitFileFITSPosName";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* reposition */
  fits_movnam_hdu(in->myFptr, ffhdutype, (char*)extname, ffextver, &status);
  if (status) {             /* it went wrong */
    Obit_log_error(err, OBIT_Error, "%s: ERROR positioning FITS file %s", 
		   routine, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_CloseErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSPosName */

/**
 * Read string keyword from current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param Value   [out] Keyword value
 * \param Comment [out] Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_NotFoundErr => not found
 */
ObitIOCode 
ObitFileFITSReadKeyStr (ObitFileFITS *in, gchar *Name, 
			gchar *Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar *routine = "ObitFileFITSReadKeyStr";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;
 
  fits_read_key_str (in->myFptr, (char*)Name, (char*)Value, (char*)Comment, &status);
  /* Not found is OK */
  if (status==KEY_NO_EXIST) return OBIT_IO_NotFoundErr;

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR reading keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_ReadErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSReadKeyStr */

/**
 * Read float keyword from current HDU
 * \param in      Pointer to object to be readwritten
 * \param Name    Keyword name
 * \param Value   [out] Keyword value
 * \param Comment [out] Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_NotFoundErr => not found
 */
ObitIOCode 
ObitFileFITSReadKeyFlt (ObitFileFITS *in, gchar *Name, 
			ofloat *Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0;
  float ffValue;
  gchar *routine = "ObitFileFITSReadKeyFlt";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;
 
  fits_read_key_flt (in->myFptr, (char*)Name, &ffValue, (char*)Comment, &status);
  *Value = (ofloat)ffValue;
  /* Not found is OK */
  if (status==KEY_NO_EXIST) return OBIT_IO_NotFoundErr;

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR reading keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_ReadErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSReadKeyFlt */

/**
 * Read double keyword from current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param Value   [out] Keyword value
 * \param Comment [out] Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_NotFoundErr => not found
 */
ObitIOCode 
ObitFileFITSReadKeyDbl (ObitFileFITS *in, gchar *Name, 
			odouble *Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0;
  double ffValue;
  gchar *routine = "ObitFileFITSReadKeyDbl";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;
 
  fits_read_key_dbl (in->myFptr, (char*)Name, &ffValue, (char*)Comment, &status);
  *Value = (odouble)ffValue;
  /* Not found is OK */
  if (status==KEY_NO_EXIST) return OBIT_IO_NotFoundErr;

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR reading keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_ReadErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSReadKeyDbl */

/**
 * Read long keyword from current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param Value   [out] Keyword value
 * \param Comment [out] Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_NotFoundErr => not found
 */
ObitIOCode 
ObitFileFITSReadKeyLng (ObitFileFITS *in, gchar *Name, 
			olong *Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  long ffivalue;
  int status=0;
  gchar *routine = "ObitFileFITSReadKeyLng";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;
 
  fits_read_key_lng (in->myFptr, Name, &ffivalue, Comment, &status);
  /* Not found is OK */
  if (status==KEY_NO_EXIST) return OBIT_IO_NotFoundErr;
  *Value = (olong)ffivalue;

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR reading keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_ReadErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSReadKeyLng */

/**
 * Read next HISTORY card from current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param hiCard  [out] HISTORY card (80 characters beginning with "HISTORY ")
 *                memory must be externally allocated.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_NotFoundErr => no more.
 */
ObitIOCode 
ObitFileFITSReadHistory (ObitFileFITS *in, gchar *hiCard, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int i, keysexist, status = 0;
  gchar tempStr[81];
  gchar *routine = "ObitFileFITSReadHistory";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* How many keys total */
  fits_get_hdrspace (in->myFptr, &keysexist, NULL, &status);

  /* cfitsio easy loses it's attention, remind it where it's at */
  /* Loop through keys looking for HISTORY */
  for (i=in->LastKeyword+1; i<=keysexist; i++) {
    fits_read_record (in->myFptr, i, tempStr, &status);
    if (!strncmp (tempStr, "HISTORY ", 8)) {  /* Found one?*/
      strncpy (hiCard, tempStr, 80);
      in->LastKeyword = i;  /* remember */
      return OBIT_IO_OK;
    }
  }
    
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR reading HISTORY from FITS file %s", 
		   routine, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_ReadErr;
  }

   in->LastKeyword = i;  /* remember */
  /* If you get here there are no more */
  return OBIT_IO_NotFoundErr;
} /* ObitFileFITSReadHistory */

/**
 * Write string keyword to current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param update  If true use "Update" rather than write (usually what you want)
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteKeyStr (ObitFileFITS *in, gchar *Name, gboolean update,
			gchar *Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar *routine = "ObitFileFITSWriteKeyStr";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* Write or update? */
  if (update) {
    fits_update_key_str (in->myFptr, (char*)Name, (char*)Value, (char*)Comment, &status);
  } else {  /* write */
    fits_write_key_str (in->myFptr, (char*)Name, (char*)Value, (char*)Comment, &status);
  }
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;
} /* ObitFileFITSWriteKeySt */

/**
 * Write float keyword to current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param update  If TRUE use "Update" rather than "write" (usually what you want)
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteKeyFlt (ObitFileFITS *in, gchar *Name, gboolean update,
			ofloat Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0, decimals = -6;
  gchar *routine = "ObitFileFITSWriteKeyFlt";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* Write or update? */
  if (update) {
    fits_update_key_flt (in->myFptr, (char*)Name, (float)Value, decimals, (char*)Comment, &status);
  } else {  /* write */
    fits_write_key_flt (in->myFptr, (char*)Name, (float)Value, decimals, (char*)Comment, &status);
  }
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;
} /* ObitFileFITSWriteKeyFlt */

/**
 * Write double keyword to current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param update  If true use "Update" rather than write (usually what you want)
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteKeyDbl (ObitFileFITS *in, gchar *Name, gboolean update,
			odouble Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0, decimals = -12;
  gchar *routine = "ObitFileFITSWriteKeyDbl";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* Write or update? */
  if (update) {
    fits_update_key_dbl (in->myFptr, (char*)Name, (double)Value, decimals, (char*)Comment, &status);
  } else {  /* write */
    fits_write_key_dbl (in->myFptr, (char*)Name, (double)Value, decimals, (char*)Comment, &status);
  }
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;
} /* ObitFileFITSWriteKeyDbl */

/**
 * Write long keyword to current HDU
 * \param in      Pointer to object to be read
 * \param Name    Keyword name
 * \param update  If true use "Update" rather than write (usually what you want)
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteKeyLng (ObitFileFITS *in, gchar *Name, gboolean update,
			olong Value, gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0;
  long ffValue = (long)Value;
  gchar *routine = "ObitFileFITSWriteKeyLng";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  /* Write or update? */
  if (update) {
    fits_update_key_lng (in->myFptr, (char*)Name, ffValue, (char*)Comment, &status);
  } else {  /* write */
    fits_write_key_lng (in->myFptr, (char*)Name, ffValue, (char*)Comment, &status);
  }
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file %s", 
		   routine, Name, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;
} /* ObitFileFITSWriteKeyLng */

/**
 * Write next HISTORY card to current HDU
 * \param in      Pointer to object to be written
 * \param Name    Keyword name
 * \param hiCard  HISTORY card if longer than 70 characters, it will be split.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK, OBIT_IO_NotFoundErr => not found
 */
ObitIOCode 
ObitFileFITSWriteHistory (ObitFileFITS *in, gchar *hiCard, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar *routine = "ObitFileFITSWriteHistory";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;
 
  fits_write_history (in->myFptr, (char*)hiCard, &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR %d writing HISTORY to FITS file %s", 
		   routine, status, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSWriteHistory */

/**
 * Increase the number of keywords available in current HDU
 * \param in       Pointer to object to be modified
 * \param morekeys Number of keyword entries to be added.
 * \param err      ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSAddKeys (ObitFileFITS *in, olong morekeys, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0, ffmorekeys=(int)morekeys;
  gchar *routine = "ObitFileFITSAddKeys";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;
  if (morekeys<=0) return OBIT_IO_OK;  /* anything to do? */
 
  fits_set_hdrsize (in->myFptr, ffmorekeys, &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR %d expanding header of FITS file %s", 
		   routine, status, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

  return OBIT_IO_OK;
} /* ObitFileFITSWriteAddKeys */

/**
 * Write Date keyword to current HDU
 * \param in      Pointer to object to be written
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteDate (ObitFileFITS *in, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar *routine = "ObitFileFITSWriteDate";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->status==OBIT_ErrorExist) return retCode;

  fits_write_date (in->myFptr, &status);
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing date in FITS file %s", 
		   routine, in->fileName);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }
  return OBIT_IO_OK;
} /* ObitFileFITSWriteDate */

/**
 * Write/update string HISTORY keyword to current HDU
 * Looks for HISTORY cards in current HDU with "HISTORY "
 * immediately followed by Name (e.g. "AIPS   NAXIS   ", exact match)
 * The card will be updated or created.
 * \param inFptr  Pointer to object to be read
 * \param Name    Keyword name
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteHisKeyStr (fitsfile *inFptr, gchar *Name, gchar *Value,
			    gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int keynum, lenCard, status = 0;
  gchar hiCard[FLEN_CARD], *chkName;
  gchar *routine = "ObitFileFITSWriteHisKeyStr";

  /* error checks */
  if (err->error) return retCode;

  /* Look it up */
  chkName = g_strconcat ("HISTORY ", Name, NULL);
  keynum = (int)FindHeaderCard (inFptr, chkName, hiCard, err);
  g_free (chkName);
  if (err->error) Obit_traceback_val (err, routine, Name, retCode);

  /* Fill in hiCard */
  g_snprintf (hiCard, FLEN_CARD-1, "HISTORY %s = '%s'", Name, Value);
  lenCard = strlen(hiCard);
  if (Comment)
      g_snprintf (&hiCard[lenCard], FLEN_CARD-1-lenCard, " / %s", Comment);

  /* Replace if exists (keynum>0) or write new */
  if (keynum>0) {
    fits_modify_record (inFptr, keynum, (char*)hiCard, &status);
  } else {
    fits_write_record (inFptr, (char*)hiCard, &status);
  }

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file", 
		   routine, Name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;
} /* end ObitFileFITSWriteHisKeyStr */

/**
 * Write/update float HISTORY keyword to current HDU
 * Looks for HISTORY cards in current HDU with "HISTORY "
 * immediately followed by Name (e.g. "AIPS   NAXIS   ", exact match)
 * The card will be updated or created.
 * \param inFptr  Pointer to object to be read
 * \param Name    Keyword name
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteHisKeyFlt (fitsfile *inFptr, gchar *Name, ofloat Value, 
			    gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int keynum, lenCard, status = 0;
  gchar hiCard[FLEN_CARD], *chkName;
  gchar *routine = "ObitFileFITSWriteHisKeyFlt";

  /* error checks */
  if (err->error) return retCode;

  /* Look it up */
  chkName = g_strconcat ("HISTORY ", Name, NULL);
  keynum = (int)FindHeaderCard (inFptr, chkName, hiCard, err);
  g_free (chkName);
  if (err->error) Obit_traceback_val (err, routine, Name, retCode);

  /* Fill in Card */
  g_snprintf (hiCard, FLEN_CARD-1, "HISTORY %s = %12.5E", Name, (float)Value);
  lenCard = strlen(hiCard);
  if (Comment)
      g_snprintf (&hiCard[lenCard], FLEN_CARD-1-lenCard, " / %s", Comment);

  /* Replace if exists (keynum>0) or write new */
  if (keynum>0) {
    fits_modify_record (inFptr, keynum, (char*)hiCard, &status);
  } else {
    fits_write_record (inFptr, (char*)hiCard, &status);
  }

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file", 
		   routine, Name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;

} /* end ObitFileFITSWriteHisKeyFlt */

/**
 * Write/update double HISTORY keyword to current HDU
 * Looks for HISTORY cards in current HDU with "HISTORY "
 * immediately followed by Name (e.g. "AIPS   NAXIS   ", exact match)
 * The card will be updated or created.
 * \param inFptr  Pointer to object to be read
 * \param Name    Keyword name
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteHisKeyDbl (fitsfile *inFptr, gchar *Name, odouble Value, 
			    gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int keynum, lenCard, status = 0;
  gchar hiCard[FLEN_CARD], *chkName;
  gchar *routine = "ObitFileFITSWriteHisKeyDbl";

  /* error checks */
  if (err->error) return retCode;

  /* Look it up */
  chkName = g_strconcat ("HISTORY ", Name, NULL);
  keynum = (int)FindHeaderCard (inFptr, chkName, hiCard, err);
  g_free (chkName);
  if (err->error) Obit_traceback_val (err, routine, Name, retCode);

  /* Fill in hiCard */
  g_snprintf (hiCard, FLEN_CARD-1, "HISTORY %s = %20.13lE", Name, (double)Value);
  lenCard = strlen(hiCard);
  if (Comment)
      g_snprintf (&hiCard[lenCard], FLEN_CARD-1-lenCard, " / %s", Comment);

  /* Replace if exists (keynum>0) or write new */
  if (keynum>0) {
    fits_modify_record (inFptr, keynum, (char*)hiCard, &status);
  } else {
    fits_write_record (inFptr, (char*)hiCard, &status);
  }

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file", 
		   routine, Name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;
} /* end ObitFileFITSWriteHisKeyDbl */

/**
 * Write/update long HISTORY keyword to current HDU
 * Looks for HISTORY cards in current HDU with "HISTORY "
 * immediately followed by Name (e.g. "AIPS   NAXIS   ", exact match)
 * The card will be updated or created.
 * \param inFptr  Pointer to object to be read
 * \param Name    Keyword name
 * \param Value   Keyword value
 * \param Comment Comment value, NULL=> not wanted
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode 
ObitFileFITSWriteHisKeyLng (fitsfile *inFptr, gchar *Name, olong Value, 
			    gchar *Comment, ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int keynum, lenCard, status = 0;
  gchar hiCard[FLEN_CARD], *chkName;
  gchar *routine = "ObitFileFITSWriteHisKeyLng";

  /* error checks */
  if (err->error) return retCode;

  /* Look it up */
  chkName = g_strconcat ("HISTORY ", Name, NULL);
  keynum = (int)FindHeaderCard (inFptr, chkName, hiCard, err);
  g_free (chkName);
  if (err->error) Obit_traceback_val (err, routine, Name, retCode);

  /* Fill in hiCard */
  g_snprintf (hiCard, FLEN_CARD-1, "HISTORY %s =  %ld", Name, (long)Value);
  lenCard = strlen(hiCard);
  if (Comment)
      g_snprintf (&hiCard[lenCard], FLEN_CARD-1-lenCard, " / %s", Comment);

  /* Replace if exists (keynum>0) or write new */
  if (keynum>0) {
    fits_modify_record (inFptr, keynum, (char*)hiCard, &status);
  } else {
    fits_write_record (inFptr, (char*)hiCard, &status);
  }

  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR writing keyword %s FITS file", 
		   routine, Name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return OBIT_IO_WriteErr;
  }

 return OBIT_IO_OK;
} /* end ObitFileFITSWriteHisKeyLng */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFileFITSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFileFITSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFileFITSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFileFITSClassInfoDefFn (gpointer inClass)
{
  ObitFileFITSClassInfo *theClass = (ObitFileFITSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFileFITSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFileFITSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFileFITSGetClass;
  theClass->newObit       = (newObitFP)newObitFileFITS;
  theClass->ObitFileZap   = (ObitFileZapFP)ObitFileFITSZap;
  theClass->ObitCopy      = (ObitCopyFP)ObitFileFITSCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitRef       = (ObitRefFP)ObitRef;
  theClass->ObitUnref     = (ObitUnrefFP)ObitUnref;
  theClass->ObitIsA       = (ObitIsAFP)ObitIsA;
  theClass->ObitClear     = (ObitClearFP)ObitFileFITSClear;
  theClass->ObitInit      = (ObitInitFP)ObitFileFITSInit;
  theClass->ObitFileExist = (ObitFileExistFP)ObitFileFITSExist;
  theClass->ObitFileOpen  = NULL; /* different call sequence */
  theClass->ObitFileClose = (ObitFileCloseFP)ObitFileFITSClose;
  theClass->ObitFileEnd       = NULL;
  theClass->ObitFileRead      = NULL;
  theClass->ObitFileReadLine  = NULL;
  theClass->ObitFileWrite     = NULL;
  theClass->ObitFileWriteLine = NULL;
  theClass->ObitFilePad       = NULL;
  theClass->ObitFilePadFile   = NULL;
  theClass->ObitFileFlush     = NULL;

} /* end ObitFileFITSClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitFileFITSInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFileFITS *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myFptr = NULL;
} /* end ObitFileFITSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitFileFITSClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  int status = 0;
  ObitFileFITS *in = inn;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->status==OBIT_Active) || (in->status==OBIT_Modified)) {
    err = newObitErr();
    ObitFileFITSClose (in, err); 
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* free this class members */
  if (in->myFptr) fits_close_file(in->myFptr, &status); in->myFptr = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitFileFITSClear */


/**
 * Search cards in current current HDU for one beginning with Name
 * \param inFptr  Pointer FITS file
 * \param Name    Character string to compare (exact)
 * \param Card    [out] Contents of card, if found, at least 81 char
 * \param err     ObitErr for reporting errors.
 * \return card number in current HDU (1-rel) or 0 if not found
 */
static olong FindHeaderCard(fitsfile *inFptr, gchar *Name, gchar *Card,
			   ObitErr *err)
{
  ObitIOCode  retCode = OBIT_IO_SpecErr;
  int i, keysexist, nameLen, status = 0;
  gchar *routine = "ObitFileFITSReadHistory";

  /* error checks */
  if (err->error) return retCode;

  /* How many keys total */
  keysexist = 0;
  fits_get_hdrspace (inFptr, &keysexist, NULL, &status);

  nameLen = strlen (Name); /* How many characters to search */

  /* Loop through header cards */
  for (i=1; i<=keysexist; i++) {
    fits_read_record (inFptr, i, (char*)Card, &status);
    if (!strncmp (Name, Card, nameLen)) {  /* Found it? */
      Card[80] = 0;
      return i;
    }
  }
    
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR reading HISTORY from FITS file", 
		   routine);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return 0;
  }

  return 0;  /* Must not have found it */
} /* end FindHeaderCard */
