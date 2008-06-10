/* $Id: ObitIOImageFITS.c,v 1.40 2008/02/29 02:03:24 bcotton Exp $    */
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
#include "Obit.h"
#include "ObitIOImageFITS.h"
#include "ObitImageSel.h"
#include "ObitTableList.h"
#include "ObitFile.h"
#include "ObitFITS.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOImageFITS.c
 * ObitIOImageFITS class function definitions.
 * This class is derived from the ObitIO class.
 *
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOImageFITS";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOImageFITSClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOImageFITSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOImageFITSClear (gpointer in);

/** Private: Read AIPS CLEAN parameters. */
void  ObitIOImageAIPSCLEANRead(ObitIOImageFITS *in, olong *status);

/** Private: Write AIPS CLEAN parameters. */
void  ObitIOImageAIPSCLEANWrite (ObitIOImageFITS *in, olong *status);

/** Private: Copy other header keywords. */
void  ObitIOImageKeysOtherRead(ObitIOImageFITS *in, olong *status, 
			       ObitErr *err);

/** Private: Purge HISTORY AIPS keywords. */
static void PurgeAIPSHistory(ObitIOImageFITS *in, olong *status);

/** Private: Set Class function pointers. */
static void ObitIOImageFITSClassInfoDefFn (gpointer inClass);

/** Private:  Set up write for possible quantization*/
static gpointer WriteQuantInit (ObitIOImageFITS *in, olong size, 
				int *outType, gpointer outBlank);

/** Private:  Shuffle data for possible write quantization */
static void WriteQuantCopy (ObitIOImageFITS *in, olong size, gpointer *wbuff, 
			    ofloat *data, gpointer outBlank);

/** Private:  Cleanup from quantization write processing */
static gpointer WriteQuantCleanup (ObitIOImageFITS *in, gpointer wbuff);

/** Private: Fixup DSS coordinates */
static void fixDSS (ObitImageDesc *desc, ObitErr *err);

/** Private: Compute DSS plate coordinates */
static void dsscrd (odouble dssx[13], odouble dssy[13], odouble x, odouble y, 
		    odouble *xi, odouble *eta, odouble *dxidx, odouble *dxidy,
		    odouble *detadx, odouble *detady, olong *ierr);

/** Private: Compute DSS plate coordinates. */
static void eqstd (odouble ra0, odouble dec0, odouble ra, odouble dec, 
		   odouble *xi, odouble *eta, olong *ierr);

/** Private: DSS standard coordinates. */
static void dsseq (odouble xpixsz, odouble ypixsz, odouble ra0, odouble dec0, 
		   odouble xoff0, odouble yoff0, odouble dssx[13], odouble dssy[13], 
		   ofloat xpix, ofloat ypix, odouble *ra, odouble *dec, olong *ierr);

/* get Private: DSS celestial position from pixel value */
static void dsseq (odouble xpixsz, odouble ypixsz, odouble ra0, odouble dec0, 
		   odouble xoff0, odouble yoff0, odouble dssx[13], odouble dssy[13], 
		   ofloat xpix, ofloat ypix, odouble *ra, odouble *dec, olong *ierr);

/* get Private: DSS cpixel from celestial position */
static void dsspix (odouble scale, odouble xpixsz, odouble ypixsz, odouble ra0, 
		    odouble dec0, odouble xoff0, odouble yoff0, odouble dssx[13], 
		    odouble dssy[13], odouble ra, odouble dec, 
		    ofloat *xpix, ofloat *ypix, olong *ierr);

/** Private: Fixup IRAF coordinates */
static void fixIRAF (ObitImageDesc *desc, ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class on the first call.
 * \param name An optional name for the object.
 * \param info if non-NULL it is used to initialize the new object.
 * \param err  ObitErr for error messages.
 * \return the new object.
 */
ObitIOImageFITS* newObitIOImageFITS (gchar *name, ObitInfoList *info,
				     ObitErr *err)
{
  ObitIOImageFITS* out;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar tempStr[201];
  gchar *routine = "newObitIOImageFITS";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitIOImageFITSClassInit();

  /* allocate structure */
  out = ObitMemAlloc0Name(sizeof(ObitIOImageFITS), "ObitIOImageFITS");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOImageFITSInit((gpointer)out);

  /* Get any info from info input */
  if (info!=NULL) {
    if(!ObitInfoListGet(info, "Disk", &type, dim, 
			&out->disk, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "FileName", &type, dim, 
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
} /* end newObitIOImageFITS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOImageFITSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOImageFITSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIOImageFITSGetClass */

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
gboolean ObitIOImageFITSSame (ObitIO *in, ObitInfoList *in1, 
				ObitInfoList *in2, ObitErr *err)
{
  olong disk1, disk2;
  gchar *filename1, *filename2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOImageFITSSame";

  /* error checks */
  if (err->error) return same;
  errno = 0;  /* reset any system error */

  /* get file from info */
  if(!ObitInfoListGet(in1, "Disk", &type, dim, &disk1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in1, "FileName", &type, dim, 
		       (gpointer)&filename1)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileName not in InfoList Object %s",
		   routine, in->name);
    return same;
  }

  if(!ObitInfoListGet(in2, "Disk", &type, dim, &disk2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in2, "FileName", &type, dim, 
		       (gpointer)&filename2)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileName not in InfoList Object %s",
		   routine, in->name);
  }

  /* Compare */
  same = (disk1==disk2) && 
    !strncmp (filename1,filename2, 200);
 
  return same;
} /* end ObitIOImageFITSSame */

/**
 * Rename underlying files.
 * New name information is given on the info member:
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 * \param in   Pointer to object to be renamed
 * \param info Associated ObitInfoList
 * \param err   ObitErr for reporting errors.
 */
void ObitIOImageFITSRename (ObitIO *in, ObitInfoList *info, 
			    ObitErr *err)
{
  gchar *routine = "ObitIOImageFITSRename";

  /* Don't bother if NULL */
  if (!in) return;

  /* error checks */
  if (err->error) return;
  errno = 0;  /* reset any system error */

  /* Rename */
  ObitFITSRename (in, info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitIOImageFITSRename */


/**
 * Delete underlying files.
 * Delete the whole FITS file.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOImageFITSZap (ObitIOImageFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;
  gchar tempStr[201];
  gchar *routine = "ObitIOImageFITSZap";

   /* error check */
  if (err->error) return;
  errno = 0;  /* reset any system error */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Close if still open */
  if ((in->myStatus==OBIT_Modified) || (in->myStatus==OBIT_Active)) {
    retCode = ObitIOImageFITSClose (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Destroy using cfitsio */
  if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
  else strncpy (tempStr, in->FileName, 200);
  fits_open_file(&(in->myFptr), tempStr, READWRITE, &status);
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
} /* end ObitIOImageFITSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOImageFITS* ObitIOImageFITSCopy  (ObitIOImageFITS *in, 
				       ObitIOImageFITS *out, ObitErr *err)
{
  const ObitIOClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  if (err->error) return out;
  errno = 0;  /* reset any system error */
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitIOImageFITS(outName, NULL, err);
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
} /* end ObitIOImageFITSCopy */

/**
 * Initialize structures and open file.
 * The file etc. info should have been stored in the ObitInfoList.
 * The image descriptor is read if ReadOnly or ReadWrite and
 * written to disk if opened WriteOnly.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk"     OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) FITS file name.
 * \param in Pointer to object to be opened.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageFITSOpen (ObitIOImageFITS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int itemp1, itemp2, i, status = 0, simple, extend;
  long ltemp3[10];
  gint32 dim[IM_MAXDIM];
  long pcount, gcount;
  ObitInfoType type;
  gchar tempStr[201];
  ObitImageDesc* desc;
  ObitImageSel* sel;
  gchar *routine = "ObitIOImageFITSOpen";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  g_assert (in->myDesc != NULL);
  g_assert (in->mySel != NULL);
  errno = 0;  /* reset any system error */

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */


  /* get instructions from info */
  if(!ObitInfoListGet(info, "Disk", &type, dim, &in->disk, err))
    Obit_traceback_val (err, routine, in->name, retCode);
 
  /* set defaults */
  desc->IOsize = OBIT_IO_byRow;
  ObitInfoListGetTest(info, "IOBy", &type, dim, &desc->IOsize);
  if(!ObitInfoListGet(info, "FileName", &type, dim, 
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
    fits_open_file(&(in->myFptr), tempStr, READWRITE, &status);
    if ((status==FILE_NOT_OPENED) || (status==READONLY_FILE)) { 
      /* Failed - try readonly */
      status = 0;
      fits_clear_errmsg();   /* Clear cfitsio error stack */
      if (fits_open_file(&(in->myFptr), tempStr, READONLY, &status) ) {
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

    if (fits_read_imghdr (in->myFptr, IM_MAXDIM, &simple, &itemp1, 
			  &itemp2, ltemp3,
			  &pcount, &gcount, &extend, &status)) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR reading input required keywords in FITS file %s", 
		     in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      retCode = OBIT_IO_OpenErr;
      return retCode;
    }
    desc->bitpix = (olong)itemp1;
    desc->naxis  = (olong)itemp2;
    for (i=0; i<desc->naxis; i++) desc->inaxes[i] = (olong)ltemp3[i];

  /*------------------------ Read/Write ---------------------------------*/
  } else if (access == OBIT_IO_ReadWrite) {
    /* Initialize output file */
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);
    ObitTrimTrail(tempStr);  /* Trim any trailing blanks */

    /* If output doesn't exist open Write only */
    if (ObitFileExist(tempStr, err)) {  /* File already exists */
      if ( fits_open_file(&(in->myFptr), tempStr, READWRITE, &status) ) {
	Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", 
		       in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
      /* get header required keywords */
      if (fits_read_imghdr (in->myFptr, IM_MAXDIM, &simple, &itemp1, 
			    &itemp2, ltemp3,
			    &pcount, &gcount, &extend, &status)) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR readinging input required keywords in FITS file %s", 
		       in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
      desc->bitpix = (olong)itemp1;
      desc->naxis  = (olong)itemp2;
      for (i=0; i<desc->naxis; i++) desc->inaxes[i] = (olong)ltemp3[i];
    } else { /* File not previously extant */
      /* Create file */
      if (fits_create_file(&(in->myFptr), tempStr, &status)) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS file %s", in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
      
      /* create image */
      for (i=0; i<desc->naxis; i++) ltemp3[i] = (long)desc->inaxes[i];
      if (fits_create_img (in->myFptr, (int)desc->bitpix, (int)desc->naxis, 
			   ltemp3, &status) ) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS image %s", in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      } /* End create image */
    } /* end file does not previously exist */

  /*------------------------ Write Only ---------------------------------*/
  } else if (access == OBIT_IO_WriteOnly) {
    /* Initialize output file */
    /* Output file may already exist test open */
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);
    ObitTrimTrail(tempStr);  /* Trim any trailing blanks */

    /* CURSE YOU CFITSIO!!! */
    /* Open read/write to see if it's there */
    fits_open_file(&(in->myFptr), tempStr, READWRITE, &status);
    if (status==0) { /* If OK force required keywords */
      for (i=0; i<desc->naxis; i++) ltemp3[i] = (long)desc->inaxes[i];
      if (fits_resize_img (in->myFptr, (int)desc->bitpix, (int)desc->naxis, 
			   ltemp3, &status)) {
 	Obit_log_error(err, OBIT_Error, 
		       "ERROR resizing output FITS file %s", in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
     }
    }

    /* If it doesn't exist then create */
    if ((status==FILE_NOT_OPENED) || (status==END_OF_FILE)) {
      /* If the file exists but has zero size then use name with '!' */
      if (in->FileName[0]!='!') { /* Add a '!' */
	tempStr[0]='!';
	strncpy (tempStr+1, in->FileName, 200);
      } else {
	strncpy (tempStr, in->FileName, 200);
      }

      status = 0;
      fits_clear_errmsg();   /* Clear error stack */
      if (fits_create_file(&(in->myFptr), tempStr, &status)) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS file %s", in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
	
      } else if (status!=0) { /* error */
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS file %s", in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
      
      /* create image */
      for (i=0; i<desc->naxis; i++) ltemp3[i] = (long)desc->inaxes[i];
      if (fits_create_img (in->myFptr, (int)desc->bitpix, (int)desc->naxis, 
			   ltemp3, &status) ) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS image %s", in->FileName);
	Obit_cfitsio_error(err); /* copy cfitsio error stack */
	ObitFileErrMsg(err);     /* system error message*/
	retCode = OBIT_IO_OpenErr;
	return retCode;
      }
    } else if (status!=0) { /* other errors */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening output FITS file %s", in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      retCode = OBIT_IO_OpenErr;
      return retCode;
    } /* End create file/image */

    /* If it got all the way to here and has an '!' at the start of the 
       file name then remove it as it may cause trouble later */
    if(!ObitInfoListGet(info, "FileName", &type, dim, tempStr, err))
      Obit_traceback_val (err, routine, in->name, retCode);
    tempStr[dim[0]] = 0;  /* null terminate */
    if (tempStr[0]=='!') {
      for (i=1; i<dim[0]; i++) tempStr[i-1] = tempStr[i];
      tempStr[i-1] = 0;
    }
    dim[0] = strlen(tempStr); dim[1] = 1;
    ObitInfoListAlwaysPut(info, "FileName", OBIT_string, dim, tempStr);

  } else {
    /* should never get here */
    g_assert_not_reached(); 
  }

  /* save information */
  in->access = access;
  in->myStatus = OBIT_Active;
  /* ???  desc->areBlanks = FALSE;*/
  /* ???  desc->maxval = -1.0e20;*/
  /* ???  desc->minval =  1.0e20;*/

  /* initialize location in image */
  desc->row   = 0;
  desc->plane = 0;

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOImageFITSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOImageFITSClose (ObitIOImageFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;

  /* error checks */
  if (err->error) return retCode;
  errno = 0;  /* reset any system error */
  g_assert (ObitIsA(in, &myClassInfo));
  /* don't bother if it's not open */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) 
    return OBIT_IO_OK;

  fits_close_file (in->myFptr, &status);
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR closing FITS file");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_CloseErr;
    return retCode;
  }

  in->myStatus = OBIT_Inactive;
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOImageFITSClose */

/**
 * initialize I/O - position to beginning of image.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageFITSSet (ObitIOImageFITS *in, ObitInfoList *info, 
			       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  int hdutype, status = 0;

  /* Position to HDU 1 */
  fits_movabs_hdu (in->myFptr, 1, &hdutype, &status);
  if ((status !=0) || (hdutype!=IMAGE_HDU)) {
    Obit_log_error(err, OBIT_Error, "ERROR positioning FITS file");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    if ((in->access==OBIT_IO_ReadOnly) ||
	(in->access==OBIT_IO_ReadWrite)) retCode = OBIT_IO_ReadErr;
    else  retCode = OBIT_IO_WriteErr;
    return retCode;
  }

  /* Reset values in descriptor */
  ((ObitImageDesc*)in->myDesc)->plane = 0;
  ((ObitImageDesc*)in->myDesc)->row   = 0;

  return retCode;
} /* end ObitIOImageFITSSet */

/**
 * Read image data from disk.
 * Reads row in->myDesc->row + 1; plane in->myDesc->plane + 1
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK,
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOImageFITSRead (ObitIOImageFITS *in, ofloat *data, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong row, plane;
  long bblc[IM_MAXDIM], ttrc[IM_MAXDIM], incs[IM_MAXDIM]={1,1,1,1,1,1,1};
  long inaxes[10];
  int group=0, i, anyf, status = 0;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  ofloat fblank = ObitMagicF();

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */

  Obit_retval_if_fail (((in->myStatus==OBIT_Active) || (in->myStatus==OBIT_Modified)) , 
		       err, retCode,
		       "Cannot read, I/O not currently active");
  
  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  row   = MAX (desc->row, sel->blc[1]-1);
  plane = MAX (desc->plane, sel->blc[2]-1);
  /* set current request by desc->IOsize */
  if (desc->IOsize==OBIT_IO_byRow) {

    plane = MAX (1, plane);
   /* increment row */
    row++;
    if (row>sel->trc[1]) { /* next plane */
      row = sel->blc[1];
      plane++;
    }

    /* Set window */
    bblc[0] = sel->blc[0];
    bblc[1] = row;
    ttrc[0] = sel->trc[0];
    ttrc[1] = row;
  } else if (desc->IOsize==OBIT_IO_byPlane) {
    row = sel->blc[1];
    /* increment plane */
    plane++;
    /* Set window */
    bblc[0] = sel->blc[0];
    bblc[1] = sel->blc[1];
    ttrc[0] = sel->trc[0];
    ttrc[1] = sel->trc[1];
  }
  desc->row   = row;
  desc->plane = plane;

  /* check if done - starting on the plane past the highest. */
  if (plane > sel->trc[2]) {
    /* ObitIOImageFITSClose (in, err); Close */
    return OBIT_IO_EOF;
  }

  /* set plane */
  bblc[2] = plane;
  ttrc[2] = plane;
  for (i=3;i<desc->naxis;i++) {ttrc[i] = 1; bblc[i] = 1;}
  for (i=0; i<desc->naxis; i++) inaxes[i] = (long)desc->inaxes[i];

  /*  Read selected portion of input */
  if (fits_read_subset_flt (in->myFptr, group, (int)desc->naxis, 
			    inaxes, 
			    bblc, ttrc, incs, (float)fblank, (float*)data, 
			    &anyf, &status)) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR reading input FITS file %s plane %d row  %d", 
		   in->FileName, plane, row);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* keep track of blanking */
  desc->areBlanks = desc->areBlanks || anyf;

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOImageFITSRead */

/**
 * Write information to disk.
 * Writes row in->myDesc->row + 1; plane in->myDesc->plane + 1
 * Writing partial images is only supported in row at a time mode.
 * When OBIT_IO_EOF is returned the image has been written,
 * data in data is ignored and the I/O is closed.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOImageFITSWrite (ObitIOImageFITS *in, ofloat *data, 
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  long size, fpixel[IM_MAXDIM]={1,1,1,1,1,1,1};
  olong i, offset, len=0, iRow, nRows=0, row, plane;
  int status = 0;
  ofloat val, fblank = ObitMagicF();
  gboolean windowed;
  gpointer wbuff=NULL;
  int outType;
  gchar outBlank[12];  /* Blob of unspecified type for type dependent blanking value */
  gchar *routine = "ObitIOImageFITSWrite";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);
  g_assert (data != NULL);
  errno = 0;  /* reset any system error */
  /* is I/O active */
  Obit_retval_if_fail(((in->myStatus==OBIT_Modified) ||
		       (in->myStatus==OBIT_Active)), 
		      err, retCode, 
		      "Cannot write, I/O not currently active");

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* If windowing, writes of a plane at a time have to be done by row */
  windowed = (sel->blc[0]>1) || (sel->blc[2]>1) ||
    (sel->trc[0]!=desc->inaxes[0]) || 
    (sel->trc[1]!=desc->inaxes[1]);

  /* set cfitsio request parameters */
  row = MAX (desc->row, sel->blc[1]-1);
  plane = MAX (desc->plane, sel->blc[2]-1);

  /* set current request by desc->IOsize */
  if (desc->IOsize==OBIT_IO_byRow) {
    plane = MAX (1, plane);
    row++; /* increment row */
    nRows = 1;
    if (row>sel->trc[1]) { /* next plane */
      row = sel->blc[1];
      plane++;
    }

    len = sel->trc[0] - sel->blc[0]+1;   /* size of a transfer (row) */

  } else if (desc->IOsize==OBIT_IO_byPlane) {

    /* increment plane */
    plane++;
    row = sel->blc[1]; /* set row */

     if (windowed) {
      nRows = sel->trc[1] - sel->blc[1] + 1;
      len = sel->trc[0] - sel->blc[0]+1;   /* size of a transfer (row) */
    } else {   /* all at once */
      nRows = 1;
      len = desc->inaxes[0] * desc->inaxes[1]; /* whole plane */
    }
   }

  /* Set first pixel, size */
  fpixel[0] = sel->blc[0];
  fpixel[1] = row;
  fpixel[2] = plane;

  desc->row   = row;
  desc->plane = plane;

   /* check if done - starting on the plane past the highest. */
  if (plane > sel->trc[2]) {
    /* ObitIOImageFITSClose (in, err); Close */
    Obit_log_error(err, OBIT_Error, 
		   "%s: Attempt to write past end of %s", 
		   routine, in->FileName);
    return OBIT_IO_EOF;
  }

  size = len;           /* transfer size in floats */

  offset = 0; /* offset in input buffer */

  /* Set up for possible quantization */
  wbuff = WriteQuantInit (in, size, &outType, (gpointer)outBlank);

  /* write file one row/plane at a time - loop for windowed planes */
  for (iRow=0; iRow<nRows; iRow++) {

    /* Shuffle data for possible quantization */
    WriteQuantCopy (in, size, &wbuff, &data[offset], outBlank);

    /* write image data */
    if (fits_write_pixnull (in->myFptr, outType, fpixel, size,
			    (void*)wbuff, (void*)outBlank, &status)) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR %d writing output FITS file %s plane %d row  %d", 
		     status, in->FileName, plane, row);
      Obit_cfitsio_error(err);  /* copy cfitsio error stack */
      ObitFileErrMsg(err);      /* system error message*/
      if (wbuff) g_free(wbuff); /* cleanup */
      return OBIT_IO_WriteErr;
    }

    /* keep track on max/min/blanking */
    for (i=0; i<size; i++) {
      val = data[offset+i];
      if (val==fblank) {
	desc->areBlanks = TRUE;
      } else { /* OK */
	desc->maxval = MAX (desc->maxval, val);
	desc->minval = MIN (desc->minval, val);
      }
    }

    row ++;            /* next row */
    fpixel[1] = row;
    offset   += len;   /* offset in data buffer */
  } /* end loop writing */

  /* Cleanup from quantization processing */
  wbuff = WriteQuantCleanup (in, wbuff);

  in->myStatus = OBIT_Modified;
  in->dataMod  = TRUE;
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOImageFITSWrite */

/**
 * Read image Descriptor data from disk.
 * \param in Pointer to object with ObitImageDesc to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageFITSReadDescriptor (ObitIOImageFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar commnt[FLEN_COMMENT], *today=NULL;
  gchar cdata[IM_MAXDIM][FLEN_CARD], sdata[FLEN_CARD], *cdum[IM_MAXDIM];
  gchar wcsname[80], pltlabel[80];
  odouble cdtest;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  int i, nfound, nhdu, hdutype, status = 0, xstatus = 0;
  int temp=0;
  long extver, ltemp;
  float ftemp, farr[10];
  double dtemp, darr[10];
  olong otemp;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  ObitTableList* tableList;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */
  tableList = (ObitTableList*)in->tableList;

  /* Index tables in file and update TableList if not already done*/
  if (tableList->number <= 0) {
    fits_get_num_hdus (in->myFptr, &nhdu, &status);
    for (i=1; i<=nhdu; i++) {
      fits_movabs_hdu (in->myFptr, i, &hdutype, &status);
      if (hdutype==BINARY_TBL) { /* If it's a table enter it in the list */
	/* table name */
	fits_read_key_str (in->myFptr, "EXTNAME", (char*)sdata, (char*)commnt, &status);
	/* version number default to 0 */
	extver = 0;
	fits_read_key_lng (in->myFptr, "EXTVER", &extver, commnt, &xstatus);
	if (status==0) { /* Add to TableList */
	  otemp = (olong)extver;
	  ObitTableListPut (tableList, sdata, &otemp, NULL, err);
	  if (err->error)
	    Obit_traceback_val (err, "ObitIOImageFITSReadDescriptor", 
				tableList->name, OBIT_IO_OpenErr);
	}
      }
    } /* end loop indexing file */
  } /* end update Table List */

  /* Position to HDU 1, the image 1 */
  fits_movabs_hdu (in->myFptr, 1, &hdutype, &status);

  /* Read keyword values, use default where possible */
  ftemp = (float)desc->maxval;
  fits_read_key_flt (in->myFptr, "DATAMAX", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0; 
  desc->maxval = (ofloat)ftemp;

  ftemp = (float)desc->minval;
  fits_read_key_flt (in->myFptr, "DATAMIN", &ftemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->minval = (ofloat)ftemp;

  strncpy (desc->teles, "        ", IMLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "TELESCOP", (char*)cdata[0], (char*)commnt, &status);
  if (status==0) strncpy (desc->teles, cdata[0], IMLEN_VALUE); desc->teles[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->instrument, "        ", IMLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "INSTRUME", (char*)cdata[0], (char*)commnt, &status);
  if (status==0) strncpy (desc->instrument, cdata[0], IMLEN_VALUE); desc->instrument[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->observer, "        ", IMLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "OBSERVER", (char*)cdata[0], (char*)commnt, &status);
  if (status==0) strncpy (desc->observer, cdata[0], IMLEN_VALUE);desc->observer[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->object, "        ", IMLEN_VALUE-1);
  fits_read_key_str (in->myFptr, "OBJECT", (char*)cdata[0], (char*)commnt, &status);
  if (status==0) strncpy (desc->object, cdata[0], IMLEN_VALUE); desc->object[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->bunit, "        ", IMLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "BUNIT", (char*)cdata[0], (char*)commnt, &status);
  if (status==0) strncpy (desc->bunit, cdata[0], IMLEN_VALUE); desc->bunit[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->obsdat, "        ", IMLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "DATE-OBS", (char*)cdata[0], (char*)commnt, &status);
  if (status==0)  strncpy (desc->obsdat, cdata[0], IMLEN_VALUE); desc->obsdat[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  today = ObitToday();
  strncpy (desc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
  fits_read_key_str (in->myFptr, "DATE-MAP", (char*)cdata[0], (char*)commnt, &status);
  if (status==0)  strncpy (desc->date, cdata[0], IMLEN_VALUE); desc->date[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  strncpy (desc->origin, "        ", IMLEN_VALUE-1); 
  fits_read_key_str (in->myFptr, "ORIGIN", (char*)cdata[0], (char*)commnt, &status);
  if (status==0)  strncpy (desc->origin, cdata[0], IMLEN_VALUE); desc->origin[IMLEN_VALUE-1] = 0;
  if (status==KEY_NO_EXIST) status = 0;

  for (i=0; i<IM_MAXDIM; i++) strncpy (desc->ctype[i], "        ", IMLEN_KEYWORD-1);
  for (i=0; i<IM_MAXDIM; i++) cdum[i] = &cdata[i][0];
  fits_read_keys_str (in->myFptr, "CTYPE", 1, IM_MAXDIM, cdum, &nfound, 
		      &status);
    if (status==0) {
      for (i=0; i<nfound; i++) strncpy (desc->ctype[i], cdata[i], IMLEN_KEYWORD-1); desc->ctype[i][IMLEN_VALUE-1] = 0;
    }
  if (status==KEY_NO_EXIST) status = 0;

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
  
  for (i=0; i<IM_MAXDIM; i++) farr[i] = 0.0;
  fits_read_keys_flt (in->myFptr, "CDELT", 1, IM_MAXDIM, farr,
		      &nfound, &status);
  if (status==KEY_NO_EXIST) status = 0;
  for (i=0; i<IM_MAXDIM; i++) desc->cdelt[i] = (gfloat)farr[i];

  for (i=0; i<IM_MAXDIM; i++) farr[i] = 0.0;
  fits_read_keys_flt (in->myFptr, "CRPIX", 1, IM_MAXDIM, farr, 
		      &nfound, &status);
  if (status==KEY_NO_EXIST) status = 0;
  for (i=0; i<IM_MAXDIM; i++) desc->crpix[i] = (ofloat)farr[i];
     
  for (i=0; i<IM_MAXDIM; i++) farr[i] = 0.0;
  fits_read_keys_flt (in->myFptr, "CROTA", 1, IM_MAXDIM, farr,
		      &nfound, &status);
  if (status==KEY_NO_EXIST) status = 0;
  for (i=0; i<IM_MAXDIM; i++) desc->crota[i] = (ofloat)farr[i];

  for (i=0; i<IM_MAXDIM; i++) darr[i] = 0.0;
  fits_read_keys_dbl (in->myFptr, "CRVAL", 1, IM_MAXDIM, darr, 
		      &nfound, &status);
  if (status==KEY_NO_EXIST) status = 0;
  for (i=0; i<IM_MAXDIM; i++) desc->crval[i] = (odouble)darr[i];

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "OBSRA", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->obsra = (odouble)dtemp;

  dtemp = 0.0;
  fits_read_key_dbl (in->myFptr, "OBSDEC", &dtemp, (char*)commnt, &status);
  if (status==KEY_NO_EXIST) status = 0;
  desc->obsdec = (odouble)dtemp;
 
  /*----------------new --------------------------*/
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
  desc->VelReference = temp - 256*desc->VelDef;

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

  /* AIPS Clean parameters */
  ObitIOImageAIPSCLEANRead (in, &status);

  /* Look for anything else and add it to the InfoList on desc */
  ObitIOImageKeysOtherRead(in, &status, err);
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, "ObitIOImageFITSReadDescriptor", 
			in->name, retCode);

  /* Trap DSS images  WCSNAME, PLTLABEL */
  strcpy (wcsname, "     ");
  ObitInfoListGetTest(desc->info, "WCSNAME", &type, dim, wcsname);
  strcpy (pltlabel, "     ");
  ObitInfoListGetTest(desc->info, "PLTLABEL", &type, dim, pltlabel);
  if (!strncmp(wcsname,"DSS   ",6) && (strncmp(pltlabel,"      ",6))) 
    fixDSS (desc, err);

  /* Trap IRAF images - CD1_1 */
  if (ObitInfoListGetTest(desc->info, "CD1_1", &type, dim, &cdtest))
    fixIRAF (desc, err);

  /* was there an error? */
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file header");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* enforce defaults */
  ObitImageSelDefault(desc, sel);

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOImageFITSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in Pointer to object with ObitImageDesc to be written.
 *           If infoList member contains "Quant" entry and it is
 *           > 0.0 and an integer output (Bitpix 16, 32) is specified 
 *           then the output will be quantized at this level.  
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitIOImageFITSWriteDescriptor (ObitIOImageFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar keywrd[FLEN_KEYWORD], commnt[FLEN_COMMENT+1];
  int i, status = 0;
  long ltemp;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  ofloat quant = 0.0;
  gchar keyName[FLEN_KEYWORD+1], *keyNameP;
  ObitInfoType keyType;
  gint32 nkey, dim[ MAXINFOELEMDIM];
  union blobEquiv {
    gchar    s[201];
    double   d;
    float    f;
    gboolean b;
    olong    i;
  } blob;
  gchar *routine = "ObitIOImageFITSWriteDescriptor";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* enforce defaults */
  ObitImageSelDefault(desc, sel);

  /* if output is an integer type, data has been modified and meaningful max/min - 
     set scaling, blanking */
  if ((desc->bitpix>0) && in->dataMod && (desc->maxval > desc->minval)) {

    /* See if Quantization desired */
    ObitInfoListGetTest(desc->info, "Quant", &keyType, dim, (gpointer*)&quant);
    if (quant>0.0) {  /* Tell about it */
      Obit_log_error(err, OBIT_InfoErr, "%s: quantizing output at %f",routine,quant);
    }
    
    /* Set image scaling */
    ObitIOImageFITSUpdateScale (in, quant, err);
    if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);
    
    /* end special handling for integers */
  }

  /*  Write keyword values */
  for (i=0; i<desc->naxis; i++) {
    strncpy (commnt, "Axis type ", FLEN_COMMENT);
    sprintf (keywrd, "CTYPE%d",i+1);
    fits_update_key_str (in->myFptr, (char*)keywrd, (char*)desc->ctype[i], (char*)commnt, 
			 &status);
    strncpy (commnt, "Axis coordinate increment", FLEN_COMMENT);
    sprintf (keywrd, "CDELT%d",i+1);
    fits_update_key_flt (in->myFptr, (char*)keywrd, (float)desc->cdelt[i], 6, (char*)commnt, 
			 &status);
    strncpy (commnt, "Axis coordinate reference pixel", FLEN_COMMENT);
    sprintf (keywrd, "CRPIX%d",i+1);
    fits_update_key_flt (in->myFptr, (char*)keywrd, (float)desc->crpix[i], 6, (char*)commnt, 
			 &status);
    strncpy (commnt, "Axis coordinate rotation", FLEN_COMMENT);
    sprintf (keywrd, "CROTA%d",i+1);
    fits_update_key_flt (in->myFptr, (char*)keywrd, (float)desc->crota[i], 6, (char*)commnt, 
			 &status);
    strncpy (commnt, "Axis coordinate value at CRPIX", FLEN_COMMENT);
    sprintf (keywrd, "CRVAL%d",i+1);
    fits_update_key_dbl (in->myFptr, (char*)keywrd, (double)desc->crval[i], 12, (char*)commnt, 
			 &status);
  }
  strncpy (commnt, "Observed Right Ascension", FLEN_COMMENT);
  fits_update_key_dbl (in->myFptr, "OBSRA", (double)desc->obsra, 12, (char*)commnt, 
		       &status);
  strncpy (commnt, "Observed declination ", FLEN_COMMENT);
  fits_update_key_dbl (in->myFptr, "OBSDEC", (double)desc->obsdec, 12, (char*)commnt, 
		       &status);
  strncpy (commnt, "Name of object", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "OBJECT", (char*)desc->object,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Telescope used", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "TELESCOP", (char*)desc->teles,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Instrument used", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "INSTRUME", (char*)desc->instrument,  commnt, 
		       &status);
  strncpy (commnt, "Observer/project", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "OBSERVER", (char*)desc->observer,  commnt, 
		       &status);
  strncpy (commnt, "Date (yyyy-mm-dd) of observ(char*)ation", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "DATE-OBS", desc->obsdat, (char*)commnt, 
		       &status);
  strncpy (commnt, "Date (yyyy-mm-dd) created ", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "DATE-MAP", (char*)desc->date, (char*)commnt, 
		       &status);
  strncpy (commnt, "Software last writing file", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "ORIGIN", (char*)desc->origin, (char*)commnt, 
		       &status);
  strncpy (commnt, "Celestial coordiate epoch", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "EPOCH", (float)desc->epoch, 6,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Celestial coordiate equinox", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "EQUINOX", (float)desc->epoch, 6,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Maximum in data array", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "DATAMAX", (float)desc->maxval, 8,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Minimum in data array", FLEN_COMMENT);
  fits_update_key_flt (in->myFptr, "DATAMIN", (float)desc->minval, 8,  (char*)commnt, 
		       &status);
  strncpy (commnt, "Image pixel units", FLEN_COMMENT);
  fits_update_key_str (in->myFptr, "BUNIT", (char*)desc->bunit, (char*)commnt, &status);
  
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

  /* AIPS Clean parameters */
  ObitIOImageAIPSCLEANWrite (in, &status);

  /* Write other keywords from descriptor */
  if (desc->info) nkey = desc->info->number; /* How many keywords? */
  else nkey = 0;
  retCode = OBIT_IO_WriteErr;
  strncpy (commnt, "             ", FLEN_COMMENT);
  for (i=0; i<nkey; i++) {
    /* Read from ObitInfoList */
    ObitInfoListGetNumber(desc->info, i, &keyNameP, &keyType, dim, 
			  blob.s, err);
    /* Copy, possibly truncating name */
    strncpy (keyName, keyNameP, FLEN_KEYWORD); keyName[FLEN_KEYWORD-1] = 0;
    if (err->error)  Obit_traceback_val (err, routine, in->name, retCode);
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
    } else if (keyType==OBIT_long) { 
      fits_update_key_lng (in->myFptr, (char*)keyName, (long)blob.i, (char*)commnt, 
			   &status);
    } else if (keyType==OBIT_oint) { 
      fits_update_key_lng (in->myFptr, (char*)keyName, (long)blob.i, (char*)commnt, 
			   &status);
    } else if (keyType==OBIT_bool) { 
      fits_update_key_log (in->myFptr, (char*)keyName, (int)blob.b, (char*)commnt, 
			   &status);
    }
  } /* end loop writing additional keywords */

  /* an error? */
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR writing output FITS file header");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_WriteErr;
    return retCode;
  }

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOImageFITSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageFITSFlush (ObitIOImageFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* cfitsio does the buffer flushing on close */

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOImageFITSFlush */

/**
 * Create buffer approptiate for I/O request.
 * Not actually used for Images.
 * Should be called after ObitIO is opened.
 * \param data (output) pointer to data array
 * \param size (output) size of data array in floats.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 */
void 
ObitIOImageFITSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOImageFITS *in, ObitInfoList *info, 
			     ObitErr *err)
{
  /* just return - never called */
} /* end ObitIOImageFITSCreateBuffer */

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
newObitIOImageFITSTable (ObitIOImageFITS *in, ObitIOAccess access, 
			 gchar *tabType, olong *tabVer, ObitErr *err)
{
  ObitTable *out;
  olong version;
  gboolean gotIt;
  gchar *outName, tabName[51];

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
  if (err->error)
    Obit_traceback_val (err, "newObitIOImageFITSTable", in->name, NULL);

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
 if (err->error)
   Obit_traceback_val (err, "newObitIOImageFITSTable", in->name, NULL);
 
 
 /* register it in the TableList */
 ObitTableListPut ((ObitTableList*)in->tableList, tabType, &version, 
		   out, err);
 if (err->error)
   Obit_traceback_val (err, "newObitIOImageFITSTable", in->name, NULL);

 return (Obit*)out;
} /* end newObitIOImageFITSTable */

/**
 * Update any disk resident structures about the current tables.
 * Nothing is needed for FITS files.
 * \param in   Pointer to object to be updated.
 * \param info ObitInfoList of parent object (not used here).
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOImageFITSUpdateTables (ObitIOImageFITS *in, ObitInfoList *info, 
					ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  return retCode;
} /* end ObitIOImageFITSUpdateTables */

/**
 * Update Header scaling parameters
 * If integer image then the scaling and any blanking values are set.
 * \param in    Pointer to object to be updated.
 * \param quant Quantization desired, 0.0 use pixel value range to set scaling
 * \param err   ObitErr for reporting errors.
 */
void ObitIOImageFITSUpdateScale (ObitIOImageFITS *in, ofloat quant,
				 ObitErr *err)
{
  long magic=0;
  float bscale = -1.0, bzero = 0.0;
  ofloat rx, rn;
  /*  odouble scale, zero;*/
  int status = 0;
  ObitImageDesc* desc;
  gchar commnt[FLEN_COMMENT+1];
  gchar *routine = "ObitIOImageFITSUpdateScale";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* Init scaling values */
  in->bscale = bscale;
  in->bzero  = bzero;
  in->magic  = magic;

  desc = in->myDesc; /* descriptor pointer */

  /* See if anything useful to be done, if integer image, set scaling, blanking */
  if ((desc->bitpix<=0) || (!in->dataMod) || (desc->maxval<=desc->minval)) return;

  /* scaling for integer types */
  switch (desc->bitpix) {
  case 8:
    bscale = (desc->maxval-desc->minval) /  254.0;
    bzero = desc->minval;
    magic = 255;
    break;
  case 16:
    bscale = (desc->maxval-desc->minval) / 32760.0;
    magic = -32768;
    bzero = desc->minval;
    /* quantizing? */
    if (quant>0.0) { 
      bzero = 0.0;
      bscale = quant;
      rx = fabs(desc->maxval) / quant;
      rn = fabs(desc->minval) / quant;
      rx = MAX(rx, rn);
      if (rx > 32760.0) { /* can't fit in 16 bits */
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: quantization %f too fine for 16 bits",
		       routine, quant);
	Obit_log_error(err, OBIT_InfoErr, "%s: max %f min %f",
		       routine, rx, rn);
	return;
      }
    }
    break;
  case 32:
    bscale =  (desc->maxval-desc->minval) / 2147483600.0;
    bzero = desc->minval;
    magic = -2147483648L; 
    /* quantizing? */
    if (quant>0.0) { 
      bzero = 0.0;
      bscale = quant;
      rx = fabs(desc->maxval) / quant;
      rn = fabs(desc->minval) / quant;
      rx = MAX(rx, rn);
      if (rx > 2147483600.0) { /* can't fit in 32 bits */
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: quantization %f too fine for 32 bits",
		       routine, quant);
	return;
      }
    }
    break;
    magic = 0;
  default:
    break;
  } /* end switch on bitpix */

  /* integer blanking? */
  /* Note this function is misspelled in the cfitsio documentation */
  if (magic!=0) fits_set_imgnull (in->myFptr, magic, &status);
  
  strncpy (commnt, "Blanking value", FLEN_COMMENT);
  fits_update_key (in->myFptr, TLONG, "BLANK", (void*)&magic, (char*)commnt, 
		   &status);

  /* set data scaling */
  strncpy (commnt, "Data scaling value", FLEN_COMMENT);
  fits_update_key (in->myFptr, TFLOAT, "BSCALE", (void*)&bscale, (char*)commnt, 
		   &status);
  strncpy (commnt, "Data offset value", FLEN_COMMENT);
  fits_update_key (in->myFptr, TFLOAT, "BZERO", (void*)&bzero, (char*)commnt, 
		   &status);
  /* Save quantization/blanking values */
  in->bscale = (ofloat)bscale;
  in->bzero  = (ofloat)bzero;
  in->magic  = (olong)magic;
  /* Scaling now done in Write */
  /* Make it REALLY do it
     scale = bscale;
     zero = bzero;
     fits_set_bscale(in->myFptr, scale, zero, &status); */

  /* an error? */
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR setting FITS scaling");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
  }
} /* end ObitIOImageFITSUpdateScale */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOImageFITSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOImageFITSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOImageFITSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOImageFITSClassInfoDefFn (gpointer inClass)
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
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOImageFITSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOImageFITSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOImageFITSGetClass;
  theClass->newObitIO     = (newObitIOFP)newObitIOImageFITS;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOImageFITSSame;
  theClass->ObitIORename  = (ObitIORenameFP)ObitIOImageFITSRename;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOImageFITSZap;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOImageFITSCopy;
  theClass->ObitClear     = (ObitClearFP)ObitIOImageFITSClear;
  theClass->ObitInit      = (ObitInitFP)ObitIOImageFITSInit;
  theClass->ObitIOOpen    = (ObitIOOpenFP)ObitIOImageFITSOpen;
  theClass->ObitIOClose   = (ObitIOCloseFP)ObitIOImageFITSClose;
  theClass->ObitIOSet     = (ObitIOSetFP)ObitIOImageFITSSet;
  theClass->ObitIORead    = (ObitIOReadFP)ObitIOImageFITSRead;
  theClass->ObitIOReadSelect = (ObitIOReadSelectFP)ObitIOImageFITSRead;
  theClass->ObitIOWrite   = (ObitIOWriteFP)ObitIOImageFITSWrite;
  theClass->ObitIOFlush   = (ObitIOFlushFP)ObitIOImageFITSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOImageFITSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOImageFITSWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOImageFITSCreateBuffer;
  theClass->newObitIOTable   = 
    (newObitIOTableFP)newObitIOImageFITSTable; 
  theClass->ObitIOUpdateTables   = 
    (ObitIOUpdateTablesFP)ObitIOImageFITSUpdateTables;
} /* end ObitIOImageFITSClassDefFn */

/*--------------- Private functions --------------------------*/
/**
 * Creates empty member objects.
 * The GType constructors will call the corresponding routines
 * for each parent class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOImageFITSInit  (gpointer inn)
{
  ObitIOImageFITS *in = inn;
  const ObitClassInfo *ParentClass;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && (ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->FileName = NULL;
  in->dataMod = FALSE;
  /* Init scaling values */
  in->bscale = 1.0;
  in->bzero  = 0.0;
  in->magic  = -32765;

} /* end ObitIOImageFITSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOImageFITSClear (gpointer inn)
{
  ObitIOImageFITS *in = inn;
  const ObitClassInfo *ParentClass;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->myStatus==OBIT_Active) ||(in->myStatus==OBIT_Modified)) {
    err = newObitErr();
    ObitIOImageFITSClose (in, err);
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  if (in->FileName) g_free(in->FileName);
  in->FileName = NULL;

 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && (ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOImageFITSClear */


/**
 * Look for rational keywords for the CLEAN parameters and
 * failing this, look in AIPS history keywords.
 * Descriptor values niter, beamMaj, beamMin, beamPA
 * \param in Pointer to ObitIOImageFITS.
 * \param lstatus (Output) cfitsio status.
 * \return return code, 0=> OK
 */
void  ObitIOImageAIPSCLEANRead(ObitIOImageFITS *in, olong *lstatus)
{
  gchar commnt[FLEN_COMMENT], card[FLEN_COMMENT], temp[FLEN_COMMENT];
  int i, j, k, keys, morekeys, status=0;
  long ltemp;
  float ftemp;
  gboolean gotBeam=FALSE, gotNiter=FALSE;
  ObitImageDesc *desc;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (*lstatus) return;  /* existing cfitsio error? */
  status = (int)(*lstatus);

  desc = in->myDesc; /* set descriptor */

  /* Attempt rational keywords */
  ftemp = -1.0;
  fits_read_key_flt (in->myFptr, "CLEANBMJ", &ftemp, (char*)commnt, &status);
  gotBeam = gotBeam || (status != KEY_NO_EXIST);
  if (status==KEY_NO_EXIST) status = 0;
  desc->beamMaj = (ofloat)ftemp;

  ftemp = -1.0;
  fits_read_key_flt (in->myFptr, "CLEANBMN", &ftemp, (char*)commnt, &status);
  gotBeam = gotBeam || (status != KEY_NO_EXIST);
  if (status==KEY_NO_EXIST) status = 0;
  desc->beamMin = (ofloat)ftemp;

  ftemp = -1.0;
  fits_read_key_flt (in->myFptr, "CLEANBPA", &ftemp, (char*)commnt, &status);
  gotBeam = gotBeam || (status != KEY_NO_EXIST);
  if (status==KEY_NO_EXIST) status = 0;
  desc->beamPA = (ofloat)ftemp;

  ltemp = -1;
  fits_read_key_lng (in->myFptr, "CLEANNIT", &ltemp, (char*)commnt, &status);
  gotNiter = gotNiter || (status != KEY_NO_EXIST);
  if (status==KEY_NO_EXIST) status = 0;
  desc->niter = (olong)ltemp;

  /* If this worked, we're done */
  if (gotBeam && gotNiter) return;

  /* Oh Well, parse all the header cards looking for: 
          1         2         3         4         5         6
0123456789012345678901234567890123456789012345678901234567890123456789
HISTORY AIPS   CLEAN BMAJ=  1.3432E-07 BMIN=  4.3621E-08 BPA= -43.11  
HISTORY AIPS   CLEAN NITER=     1000 PRODUCT=1   / NORMAL   
  */

  /* how many keywords to look at? */
  fits_get_hdrspace (in->myFptr, &keys, &morekeys, &status);
  for (k=1; k<=keys; k++) {
    fits_read_record (in->myFptr, k, card, &status);
    if (status==0) {
      if (!strncmp ("HISTORY AIPS   CLEAN BMAJ", card, 25)) {
	/* Parse card */
	for (j=0,i=26; i<38; i++) temp[j++] = card[i]; temp[j] = 0;
	sscanf (temp, "%f", &ftemp);
	desc->beamMaj = (ofloat)ftemp;
	for (j=0,i=44; i<56; i++) temp[j++] = card[i]; temp[j] = 0;
	sscanf (temp, "%f", &ftemp);
	desc->beamMin = (ofloat)ftemp;
	for (j=0,i=61; i<68; i++) temp[j++] = card[i]; temp[j] = 0;
	sscanf (temp, "%f", &ftemp);
	desc->beamPA = (ofloat)ftemp;
	gotBeam = TRUE;
      } else if (!strncmp ("HISTORY AIPS   CLEAN NITER", card, 26)) {
	for (j=0,i=27; i<36; i++) temp[j++] = card[i]; temp[j] = 0;
	sscanf (temp, "%ld", &ltemp);
	desc->niter = (olong)ltemp;
      }
      /* Are we there yet? */
      if (gotBeam && gotNiter) return;
    }
  }
  *lstatus = (olong)status;  /* return status */

} /* end ObitIOImageAIPSCLEANRead */

/**
 * Write both rational keywords and AIPS HISTORY cards
 * Descriptor values niter, beamMaj, beamMin, beamPA
 * \param in Pointer to ObitIOImageFITS.
 * \param status (Output) cfitsio status.
 * \return return code, 0=> OK
 */
void  ObitIOImageAIPSCLEANWrite (ObitIOImageFITS *in, olong *lstatus)
{
  gchar commnt[FLEN_COMMENT+1], card[FLEN_CARD+1];
  int  status=0;
  long ltemp;
  ObitImageDesc *desc;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (*lstatus) return;  /* existing error? */
  status = (int)(*lstatus);

  desc = in->myDesc; /* set descriptor */

  /* Rational keywords */
  if (desc->beamMaj > 0.0) {
     strncpy (commnt, "Convolving Gaussian major axis FWHM (deg)", 
	      FLEN_COMMENT);
     fits_update_key_flt (in->myFptr, "CLEANBMJ", (float)desc->beamMaj, 
			  6, (char*)commnt,  &status);
     strncpy (commnt, "Convolving Gaussian minor axis FWHM (deg)", 
	      FLEN_COMMENT);
     fits_update_key_flt (in->myFptr, "CLEANBMN", (float)desc->beamMin,  
			  6, (char*)commnt, &status);
     strncpy (commnt, "Convolving Gaussian position angle (deg)", 
	      FLEN_COMMENT);
     fits_update_key_flt (in->myFptr, "CLEANBPA", (float)desc->beamPA,  
			  6, (char*)commnt, &status);
  }
  if (desc->niter > 0) {
    ltemp = (long)desc->niter;
    strncpy (commnt, "Number of Clean iterations", FLEN_COMMENT);
    fits_update_key_lng (in->myFptr, "CLEANNIT", ltemp,
			 (char*)commnt, &status);
  }

  /* Purge previous version */
  PurgeAIPSHistory (in, &status);

  /* Hide 'em where AIPS can find them */
  if (desc->beamMaj > 0.0) {
    g_snprintf (card, FLEN_COMMENT,
		"AIPS   CLEAN BMAJ=%12.4e BMIN=%12.4e BPA=%7.2f",
		desc->beamMaj, desc->beamMin, desc->beamPA);
    fits_write_history (in->myFptr, (char*)card, &status);
  }

   if (desc->niter > 0) {
    g_snprintf (card, FLEN_COMMENT,
		"AIPS   CLEAN NITER=%9d", desc->niter);
    fits_write_history (in->myFptr, (char*)card, &status);
  }
} /* end ObitIOImageAIPSCLEANWrite */

/**
 * Look for additional descriptive keywords, any that are not 
 * on the exclusion list are copied to the descriptor InfoList.
 * \param in      Pointer to ObitIOImageFITS.
 * \param status (Output) cfitsio status.
 * \param err    ObitErr stack.
 * \return return code, 0=> OK
 */
void  ObitIOImageKeysOtherRead(ObitIOImageFITS *in, olong *lstatus, 
			       ObitErr *err)
{
  char keywrd[FLEN_KEYWORD], value[FLEN_VALUE], commnt[FLEN_COMMENT+1];
  char *first, *last, *anF, *aT, dtype, svalue[FLEN_VALUE];
  olong i, j, l;
  olong ivalue;
  int  k, status=0, keys, morekeys;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  double dvalue;
  ObitImageDesc *desc;
  char *exclude[] = 
  {"SIMPLE", "BITPIX", "EXTEND", "HISTORY", "COMMENT", "BLANK", "        ",
   "BSCALE", "BZERO", "NAXIS", "BLOCKED", 
   "CTYPE", "CDELT", "CRPIX", "CROTA", "CRVAL", "OBSRA", "OBSDEC", 
   "OBJECT", "TELESCOP", "DATE", "EPOCH", "DATAMAX", "DATAMIN", "BUNIT", 
   "ALTRVAL", "ALTRPIX", "VELREF", "RESTFREQ", "XSHIFT", "YSHIFT", 
   "CLEAN", NULL};
  olong number, *len;
  gboolean bvalue, bad=FALSE;
  gchar *routine = "ObitIOImageKeysOtherRead";

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (*lstatus) return;  /* existing error? */
  status = (int)(*lstatus);

  /* delete old InfoList and restart */
  ((ObitImageDesc*)in->myDesc)->info = ObitInfoListUnref (((ObitImageDesc*)in->myDesc)->info);
  ((ObitImageDesc*)in->myDesc)->info = (gpointer)newObitInfoList ();
  desc = in->myDesc; /* set descriptor */

  /* get number and length of exclusion strings */
  number = 0;
  i = 0;
  while (exclude[i]!=NULL) {
    number++;
    i++;
  }
  len = ObitMemAlloc0Name(number*sizeof(olong), routine);
  for (i=0; i<number; i++) len[i] = strlen(exclude[i]);

  /* how many keywords to look at? */
  fits_get_hdrspace (in->myFptr, &keys, &morekeys, &status);
  for (k=1; k<=keys; k++) {
    fits_read_keyn (in->myFptr, k, (char*)keywrd, value, (char*)commnt, &status);
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
	  ObitInfoListPut(desc->info, (char*)keywrd, OBIT_string, dim, 
			  (gconstpointer)svalue, err);
	  
	  break;
	case 'L':  /* logical 'T', 'F' */
	  anF   = index (value,'F'); /* Logical */
	  aT    = index (value,'T'); /* Logical */
	  bvalue = FALSE;
	  if (aT!=NULL) bvalue = TRUE;
	  /* add to InfoList */
	  dim[0] = 1;
	  ObitInfoListPut(desc->info, (char*)keywrd, OBIT_bool, dim, 
			  (gconstpointer)&bvalue, err);
	  break;
	case 'I':  /* Integer */
	  ivalue = strtol(value, NULL, 10);
	  /* add to InfoList */
	  dim[0] = 1;
	  ObitInfoListPut(desc->info, (char*)keywrd, OBIT_long, dim, 
			  (gconstpointer)&ivalue, err);
	  break;
	case 'F':  /* Float - use double */
	  /* AIPS uses 'D' for double exponent */
	  for (l=0; l<strlen(value); l++) if (value[l]=='D') value[l]='e';
	  dvalue = strtod(value, &last);
	  /* add to InfoList */
	  dim[0] = 1;
	  ObitInfoListPut(desc->info, (char*)keywrd, OBIT_double, dim, 
			  (gconstpointer)&dvalue, err);
	  break;
	case 'X':  /* Complex - can't handle */
	default:
	  g_assert_not_reached(); /* unknown, barf */
	}; /* end switch on type */
	
	/* error check */
	if (err->error)  {/* add trace and return on error */
	  ObitMemFree(len);
	  Obit_traceback_msg (err, routine, in->name);
	}
      }
    }
  } /* end loop over keywords */

  /* cleanup */
  ObitMemFree(len);

} /* end ObitIOImageKeysOtherRead */

/**
 * Delete selected "HISTORY AIPS" keywords from history
 * These include restoring beam and 
 * very dangerous traps AIPS lays for itself
 * \param in      Pointer to ObitIOImageFITS.
 * \param lstatus (Output) cfitsio status.
 * \return return code, 0=> OK
 */
static void PurgeAIPSHistory(ObitIOImageFITS *in, olong *lstatus)
{
  gchar card[FLEN_CARD+1];
  int k, keys, morekeys, status=(int)*lstatus;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (status) return;  /* existing cfitsio error? */

  /* Oh Well, parse all the header cards looking for: 
          1         2         3         4         5         6
0123456789012345678901234567890123456789012345678901234567890123456789
HISTORY AIPS   CLEAN BMAJ=  1.3432E-07 BMIN=  4.3621E-08 BPA= -43.11  
HISTORY AIPS   CLEAN NITER=     1000 PRODUCT=1   / NORMAL   
  */

  /* how many keywords to look at? */
  fits_get_hdrspace (in->myFptr, &keys, &morekeys, &status);
  for (k=keys; k>=1; k--) {
    fits_read_record (in->myFptr, k, (char*)card, &status);
    if (status==0) {
      if (!strncmp ("HISTORY AIPS   CLEAN BMAJ", card, 25)) {
	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   CLEAN NITER", (char*)card, 26)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   CTYPE", (char*)card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   CRVAL", (char*)card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   CRPIX", (char*)card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   CROTA", (char*)card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   NAXIS", (char*)card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   BSCALE", card, 21)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   BZERO", card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   BUNIT", card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   BLANK", card, 20)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY AIPS   DATE", card, 19)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY 1C", (char*)card, 10)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY 3C", (char*)card, 10)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY 4C", (char*)card, 10)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY 5C", (char*)card, 10)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY 6C", (char*)card, 10)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      } else if (!strncmp ("HISTORY 7C", (char*)card, 10)) {
 	/* Delete card */
	fits_delete_record (in->myFptr, k, &status);
      }
    }
  }
} /* end  PurgeAIPSHistory */

/**
 * Initialize possible output image pixel quantization.
 * If floats are being written, set up to pass them through.
 * \param in   Pointer to object to be written.
 * \param size Size in pixels to be written
 * \param outType  cfitsio code for output data type
 * \param outBlank [out] untyped blob containing value for output blanking
 * \return pointer to quantization buffer if created, else NULL
 */
static gpointer WriteQuantInit (ObitIOImageFITS *in, olong size, 
				int *outType, gpointer outBlank)
{
  gpointer out = NULL;
  ofloat fblank = ObitMagicF();
  double scale, zero;
  int status = 0;
  ObitImageDesc *desc = (ObitImageDesc*)in->myDesc;

  /* By output type */
  if (desc->bitpix<0) {  /* Float */
    memcpy (outBlank, &fblank, sizeof(ofloat));
    *outType = TFLOAT;
  } else if (desc->bitpix==8)  {  /* 8-bit integer */
    memcpy (outBlank, &in->magic, sizeof(olong));
    *outType = TBYTE;
    out = g_malloc0(size*sizeof(gchar));
  } else if (desc->bitpix==16) {  /* 16-bit integer */
    memcpy (outBlank, &in->magic, sizeof(olong));
    *outType = TSHORT;
    out = g_malloc0(size*sizeof(gshort));
  } else if (desc->bitpix==32) {  /* 32-bit integer */
    memcpy (outBlank, &in->magic, sizeof(olong));
    *outType = TLONG;
    out = g_malloc0(size*sizeof(olong));
  } else  {  /* Shouldn't happen - pretend float */
    memcpy (outBlank, &fblank, sizeof(ofloat));
    *outType = TFLOAT;
  }

  /* Turn off cfitsio scaling */
  if (desc->bitpix>0) {  /* integer */
    scale = 1.0;
    zero  = 0.0;
    fits_set_bscale(in->myFptr, scale, zero, &status); 
  }

  return out;
} /* end WriteQuantInit */

/**
 * If quantizing, quantize data to wbuff, else set it to data
 * \param in     Pointer to object to be written.
 * \param size   Size in pixels to be written
 * \param wbuff  [in/out] pointer of data buffer to pass to cfitsio
 * \param data   Pixel values as floats to be written.
 */
static void WriteQuantCopy (ObitIOImageFITS *in, olong size, gpointer *wbuff, 
			    ofloat *data, gpointer outBlank)
{
  olong i;
  ofloat val, zero, iscale=1.0, fblank = ObitMagicF();
  gchar  cblank, *cbuff = (gchar*)*wbuff;
  gshort sblank, *sbuff = (gshort*)*wbuff;
  olong  lblank, *lbuff = (olong*)*wbuff;
  ObitImageDesc *desc = (ObitImageDesc*)in->myDesc;

  if (in->bscale!=0) iscale = 1.0 / in->bscale;
  zero = in->bzero;
  /* By output type */
  if (desc->bitpix<0) {  /* Float */
    /* Pass through */
    *wbuff = data;
  } else if (desc->bitpix==8)  {  /* unsigned 8-bit integer */
    cblank = (gchar)in->magic;
    for (i=0; i<size; i++) {
      val = 0.5 + (data[i]-zero) * iscale;
      if (val==fblank) {
	cbuff[i] = cblank;
      } else {
	cbuff[i] = (gchar)(MAX (0.0, MIN (val,255.0)));
      }
    }
  } else if (desc->bitpix==16) {  /* signed 16-bit integer */
    sblank = (gshort)in->magic;
    for (i=0; i<size; i++) {
      val = (data[i]-zero) * iscale;
      if (val==fblank) {
	sbuff[i] = sblank;
      } else {
	if (val>0.0) sbuff[i] = (gshort)(val + 0.5);
	else sbuff[i] = (gshort)(val - 0.5);
      }
    }
  } else if (desc->bitpix==32) {  /* 32-bit integer */
    lblank = (olong)in->magic;
    for (i=0; i<size; i++) {
      val = (data[i]-zero) * iscale;
      if (val==fblank) {
	lbuff[i] = lblank;
      } else {
	if (val>0.0) lbuff[i] = (olong)(val + 0.5);
	else lbuff[i] = (olong)(val - 0.5);
      }
    }
  } else  {  /* Shouldn't happen - pretend float */
    *wbuff = data;
  }
} /* end WriteQuantCopy */


/**
 * Deallocate quantization buffer if one was created.
 * \param in    Pointer to object to be written.
 * \param wbuff Pointer to quantization buffer
 * \return NULL
 */
static gpointer WriteQuantCleanup (ObitIOImageFITS *in, gpointer wbuff)
{
  double scale, zero;
  int status = 0;
  ObitImageDesc *desc = (ObitImageDesc*)in->myDesc;

  /* free wbuff if output integer */
  if (desc->bitpix>0) {  /* integer */
    if (wbuff!=NULL) g_free(wbuff);
    /* Turn cfitsio scaling back on */
    scale = in->bscale;
    zero  = in->bzero;
    fits_set_bscale(in->myFptr, scale, zero, &status); 
  }
  return NULL;
} /* end WriteQuantCleanup */


/**
 * Convert DSS plate coordinates to WCS in descriptor
 *   Fix up coordinates in header (2 dim); give standard WCS positions.
 *   Compute position at center using DSS formula and axis increments
 *   from finite differences.  The rotation is determined by obtaining
 *   the pixel location of a position 10 sec north of the field center.
 * \param desc   Image descriptor, DSS information expected on info member
 * \param err    ObitErr stack.
 */
static void fixDSS (ObitImageDesc *desc, ObitErr *err)
{
  odouble xpixsz, ypixsz, amdx[20], amdy[20], ppo[6];
  odouble scale, equinox, ras, decs, objctx, objcty, ra0, dec0;
  olong   cnpix[2], i, rah, ram, decd, decm, ierr=0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat xcen, ycen, xdss, ydss, rot, xt, yt;
  odouble rat, dect, rae, dece, ran, decn;
  gchar keyword[32], pltlabel[80], decsign[10], objctra[80], objctdec[80] ;
  gchar *routine = "fixDSS";

  if (err->error) return;  /* previous error */

  /* Get parameters from header */
  ObitInfoListGet(desc->info, "PLTLABEL", &type, dim, pltlabel, err);
  ObitInfoListGet(desc->info, "CNPIX1  ", &type, dim, &cnpix[0], err);
  ObitInfoListGet(desc->info, "CNPIX2  ", &type, dim, &cnpix[1], err);
  ObitInfoListGet(desc->info, "XPIXELSZ", &type, dim, &xpixsz, err);
  ObitInfoListGet(desc->info, "YPIXELSZ", &type, dim, &ypixsz, err);
  ObitInfoListGet(desc->info, "PLTSCALE", &type, dim, &scale, err);
  ObitInfoListGet(desc->info, "EQUINOX ", &type, dim, &equinox, err);
  ObitInfoListGet(desc->info, "PLTRAH  ", &type, dim, &rah, err);
  ObitInfoListGet(desc->info, "PLTRAM  ", &type, dim, &ram, err);
  ObitInfoListGet(desc->info, "PLTRAS  ", &type, dim, &ras, err);
  ObitInfoListGet(desc->info, "PLTDECSN", &type, dim, decsign, err);
  ObitInfoListGet(desc->info, "PLTDECD ", &type, dim, &decd, err);
  ObitInfoListGet(desc->info, "PLTDECM ", &type, dim, &decm, err);
  ObitInfoListGet(desc->info, "PLTDECS ", &type, dim, &decs, err);
  ObitInfoListGet(desc->info, "OBJCTX  ", &type, dim, &objctx, err);
  ObitInfoListGet(desc->info, "OBJCTY  ", &type, dim, &objcty, err);
  ObitInfoListGet(desc->info, "OBJCTRA ", &type, dim, objctra, err);
  ObitInfoListGet(desc->info, "OBJCTDEC", &type, dim, objctdec, err);
  for (i=1; i<=6; i++) {
    sprintf (keyword, "PPO%d",i);
    ObitInfoListGet(desc->info, keyword, &type, dim, &ppo[i-1], err);
  }
  for (i=1; i<=20; i++) {
    sprintf (keyword, "AMDX%d",i);
    ObitInfoListGet(desc->info, keyword, &type, dim, &amdx[i-1], err);
    sprintf (keyword, "AMDY%d",i);
    ObitInfoListGet(desc->info, keyword, &type, dim, &amdy[i-1], err);
  }
  if (err->error) Obit_traceback_msg (err, routine, desc->name);

  /* Correct descriptor */
  /* Coordinate types */
  strcpy (desc->ctype[desc->jlocr], "RA---ARC");
  strcpy (desc->ctype[desc->jlocd], "DEC--ARC");

  /* Get center RA, Dec */
  xcen = desc->inaxes[desc->jlocr] / 2.0;
  ycen = desc->inaxes[desc->jlocd] / 2.0;
  ra0  = 15.0*((odouble)rah + ram/60.0 + ras/3600.0);
  dec0 = ((odouble)decd + decm/60.0 + decs/3600.0);
  if (decsign[0]=='-') dec0 = -dec0;

  /* Initial values */
  desc->crpix[desc->jlocr] = xcen;
  desc->crpix[desc->jlocd] = ycen;
  desc->crota[desc->jlocr] = 0.0;
  desc->crota[desc->jlocd] = 0.0;
  desc->cdelt[desc->jlocr] = -1.7/3600.0;
  desc->cdelt[desc->jlocd] =  1.7/3600.0;
  desc->crval[desc->jlocr] = ra0;
  desc->crval[desc->jlocd] = dec0;

  /*  Pixel numbers in original */
  xdss = xcen + cnpix[0] - 0.5;
  ydss = ycen + cnpix[1] - 0.5;

  /* Get center position */
  dsseq (xpixsz, ypixsz, ra0, dec0, ppo[2], ppo[5], amdx, amdy,
	 xdss, ydss, &desc->crval[desc->jlocr], &desc->crval[desc->jlocd], &ierr);
  if (ierr!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR determining DSS coordinates");
    return;
  }

  /* Get rotation - go 10" north */
  rat  = desc->crval[desc->jlocr];
  dect = desc->crval[desc->jlocd] + 10.0 / 3600.0;
  dsspix (scale, xpixsz, ypixsz, ra0, dec0, ppo[2], ppo[5],
          amdx, amdy, rat, dect, &xt, &yt, &ierr);
  if (ierr!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR determining DSS coordinates");
    return;
  }
  
  /* Rotation on sky */
  rot = atan2 (xt-xdss, yt-ydss);
  desc->crota[desc->jlocd] = -rot * 57.2957795;

  /*  1 pixel east */
  xdss -=  1.0;
  dsseq (xpixsz, ypixsz, ra0, dec0, ppo[2], ppo[5], amdx, amdy,
         xdss, ydss, &rae, &dece, &ierr);
  if (ierr!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR determining DSS coordinates");
    return;
  }
  /* Ra increment */
  desc->cdelt[desc->jlocr] = -(rae -  desc->crval[desc->jlocr]) * 
    cos (desc->crval[desc->jlocd]*1.74533e-2) / cos (rot);

  /* 1 pixel north */
  xdss = xdss + 1.0;
  ydss = ydss + 1.0;
  dsseq (xpixsz, ypixsz, ra0, dec0, ppo[2], ppo[5], amdx, amdy,
         xdss, ydss, &ran, &decn, &ierr);
  if (ierr!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR determining DSS coordinates");
    return;
  }

  /* Declination increment */
  desc->cdelt[desc->jlocd] = (decn - desc->crval[desc->jlocd]) / cos (rot);

  /* Observed position */
  desc->obsra  = ra0;
  desc->obsdec = dec0;
  desc->equinox = equinox;
} /* end of fixDSS */

/**
 *     dsseq computes the J2000.0 equatorial coordinates of the specified
 *     Digitized Sky Survey pixel coordinates.
 *     From the AIPS program SKYVE by Mark Calebretta 
 *
 *     Called:
 *          none
 *
 *     Algorithm:
 *          The equations for converting DSS pixel coordinates to
 *          offsets from the plate centre are given on page 10 of the
 *          booklet supplied with the Digitized Sky Survey CD set.
 *
 *          The equations for computing standard DSS plate coordinates
 *          from plate offsets are given on page 11 of the booklet
 *          supplied with the Digitized Sky Survey CD:
 *
 *             xi  = a1*x + a2*y + a3 + a4*x*x + a5*x*y + a6*y*y +
 *                   a7*(x*x+y*y) + a8*x*x*x + a9*x*x*y + a10*x*y*y +
 *                   a11*y*y*y + a12*x*(x*x+y*y) + a13*x*(x*x+y*y)**2
 *
 *             eta = b1*y + b2*x + b3 + b4*y*y + b5*x*y + b6*x*x +
 *                   b7*(x*x+y*y) + b8*y*y*y + b9*x*y*y + b10*x*x*y +
 *                   b11*x*x*x + b12*y*(x*x+y*y) + b13*y*(x*x+y*y)**2
 *
 *          The equations for computing J2000.0 right ascension and
 *          declination from the standard coordinates are given on page
 *          11 of the booklet supplied with the Digitized Sky Survey CD.
 *          However, note that there is a misprint in the equation for
 *          declination, the COS(RA0) term should be COS(RA-RA0).
 *
 *     Notes:
 *       1)
 *
 *     Author:
 *          Mark Calabretta, Australia Telescope.
 *          Origin; 1994/08/03  Code last modified; 1994/08/04
 *          Translated to c by W. D. Cotton, NRAO
 *
 * \param xpixsz  Plate pixel size in X, micron
 * \param ypixsz  Plate pixel size in Y, micron
 * \param ra0     Plate center J2000.0 right ascension (deg)
 * \param dec0    Plate center J2000.0 eclination (deg)
 * \param xoff0   Plate center x offset, micron.
 * \param yoff0   Plate center y offset, micron.
 * \param dssx    The x coefficients for the plate solution for
 *                computing standard plate coordinates
 * \param dssy    The y coefficients for the plate solution for
 *                computing standard plate coordinates
 * \param xpix   DSS x pixel coordinate
 * \param ypix   DSS y pixel coordinate
 * \param ra     [out]  Required J2000.0 right ascension (deg)
 * \param dec    [out]  Required J2000.0 declinatioin (deg)
 * \param ierr   [out]  error code, 0=>OK
 */
static void dsseq (odouble xpixsz, odouble ypixsz, odouble ra0, odouble dec0, 
		   odouble xoff0, odouble yoff0, odouble dssx[13], odouble dssy[13], 
		   ofloat xpix, ofloat ypix, odouble *ra, odouble *dec, olong *ierr)
{
  odouble cdec, eta, f, tdec, x, xi, xx, xxyy, xy, y, yy;
  odouble pi = 3.141592653589793238462643; /*     Pi. */
  odouble d2r = pi/180.0;  /* Factor to convert degrees to radians. */
  odouble as2r = d2r/3600; /* Factor to convert arcsec to radians. */

  *ierr = 0;

  /*  Compute offsets. */
  x = (xoff0 - xpixsz*xpix)/1000.0;
  y = (ypixsz*ypix - yoff0)/1000.0;

  /* Compute temporaries. */
  xx = x*x;
  yy = y*y;
  xy = x*y;
  xxyy = xx + yy;

  /* Compute standard coordinates. */
  xi =     dssx[2] +
    x*(dssx[0] + x*(dssx[3] + x*dssx[7]))  +
    y*(dssx[1] + y*(dssx[5] + y*dssx[10])) +
    xy*(dssx[4] + x*dssx[8]  + y*dssx[9])  +
    xxyy*(dssx[6] + x*dssx[11] + dssx[12]*x*xxyy);
  eta =    dssy[2] +
    y*(dssy[0] + y*(dssy[3] + y*dssy[7]))  +
    x*(dssy[1] + x*(dssy[5] + x*dssy[10])) +
    xy*(dssy[4] + y*dssy[8]  + x*dssy[9])  +
    xxyy*(dssy[6] + y*dssy[11] + dssy[12]*y*xxyy);

  /* Convert to radians. */
  xi  *= as2r;
  eta *= as2r;

  /*  Compute J2000.0 coordinates. */
  cdec = cos(dec0*d2r);
  tdec = tan(dec0*d2r);
  f = 1.0 - eta*tdec;
  *ra = atan((xi/cdec)/f)/d2r + ra0;
  *dec = atan(((eta+tdec)*cos((*ra-ra0)*d2r))/f)/d2r;
  
  if (fabs (dec0-*dec) > 90.0) {
    /* Wrong solution */
    *dec = -*dec;
    if (*ra>180.0) {
      *ra = *ra - 180.0; }
    else {
      *ra = *ra + 180.0; }
  } /* end of check for wrong solution */
} /* end dsseq */

/**
 *     DSSPIX computes the pixel coordinates in a Digitized Sky Survey
 *     plate corresponding to the specified J2000.0 equatorial
 *     coordinate.  This requires inversion of the plate solution and
 *     this it does by iteration.
 *     From the AIPS program SKYVE by Mark Calebretta 
 *
 *     Given:
 *          scale       D     Approximate plate scale, arcsec/mm.
 *          xpixsz      D     Plate pixel size in X and Y, micron.
 *      and ypixsz      D
 *          ra0,dec0    D     Plate centre J2000.0 right ascension and
 *                            declination, in degrees.
 *          xoff0,yoff0 D     Plate centre offsets, micron.
 *          dssx(13)    D     The coefficients for the plate solution for
 *          dssy(13)          computing standard plate coordinates
 *                            (XI,ETA) in arcsec from plate offsets
 *                            (X,Y), in mm.
 *          ra,dec      D     Required J2000.0 right ascension and
 *                            declination, in degrees.
 *
 *     Returned:
 *          xpix,ypix   R     DSS pixel coordinates.
 *          ierr        I     Error status, 0 means success.
 *
 *     Called:
 *          {DSSCRD, EQSTD}
 *
 *     Algorithm:
 *          The iteration formula is obtained by simultaneously solving
 *          for DX and DY from the equations for the total differentials:
 *
 *              Dx = DX*dx/dX + Dy*dx/dY
 *              Dy = DX*dy/dX + Dy*dy/dY
 *
 *          where
 *
 *                       x -> xi
 *                       y -> eta
 *                  DX, DY -> total differential of X and Y
 *              d/dX, d/dY -> partial derivative with respect to X and Y
 *
 *          The equations for converting DSS pixel coordinates to
 *          offsets from the plate centre are given on page 10 of the
 *          booklet supplied with the Digitized Sky Survey CD set.
 *
 *     Notes:
 *       1)
 *
 *     Author:
 *          Mark Calabretta, Australia Telescope.
 *          Origin; 1994/07/27  Code last modified; 1994/08/04
 *          Translated to c by W. D. Cotton
 *
 * \param scale   Approximate plate scale, arcsec/mm.
 * \param xpixsz  Plate pixel size in X, micron
 * \param ypixsz  Plate pixel size in Y, micron
 * \param ra0     Plate center J2000.0 right ascension (deg)
 * \param dec0    Plate center J2000.0 eclination (deg)
 * \param xoff0   Plate center x offset, micron.
 * \param yoff0   Plate center y offset, micron.
 * \param dssx    The x coefficients for the plate solution for
 *                computing standard plate coordinates
 * \param dssy    The y coefficients for the plate solution for
 *                computing standard plate coordinates
 * \param ra     Required J2000.0 right ascension (deg)
 * \param dec    Required J2000.0 declinatioin (deg)
 * \param xpix   [out] DSS x pixel coordinate
 * \param ypix   [out] DSS y pixel coordinate
 * \param ierr   [out]  error code, 0=>OK
 */
static void dsspix (odouble scale, odouble xpixsz, odouble ypixsz, odouble ra0, 
	     odouble dec0, odouble xoff0, odouble yoff0, odouble dssx[13], 
	     odouble dssy[13], odouble ra, odouble dec, 
	     ofloat *xpix, ofloat *ypix, olong *ierr)
{
  olong   iter, niter;
  odouble    deta, detadx, detady, dx, dxi, dxidx, dxidy, dy, eta, eta0;
  odouble    tol, xoff, xi, xi0, yoff, z;

  *ierr = 0;

  /* Initialize. */
  niter = 50;
  tol = (MIN(xpixsz, ypixsz)/100.0)/1000.0;

  /*  Convert to standard coordinates. */
  eqstd (ra0, dec0, ra, dec, &xi0, &eta0, ierr);

  /* Initial guess for plate offset. */
  xoff =  xi0/scale;
  yoff = eta0/scale;

  /*  Iterate. */
  for (iter=0; iter<=niter; iter++)  {
    /*  Compute standard coordinates and their derivatives. */
      dsscrd (dssx, dssy, xoff, yoff, &xi, &eta, &dxidx, &dxidy,
	      &detadx, &detady, ierr);

      /* Error terms. */
      dxi   =  xi0 - xi;
      deta  = eta0 - eta;
      /* Compute correction. */
      z  = dxidx*detady-dxidy*detadx;
      dx = (dxi*detady - deta*dxidy)/z;
      dy = (deta*dxidx - dxi*detadx)/z;

      /* Apply correction. */
      xoff = xoff + dx;
      yoff = yoff + dy;
      /* Test for convergence. */

      if ((fabs(dx)<tol) && (fabs(dy)<tol)) break;
    } /* end of iteration loop */

  /* Convert offsets to pixels. */
  *xpix = (xoff0 - 1000.0*xoff)/xpixsz;
  *ypix = (yoff0 + 1000.0*yoff)/ypixsz;
} /* end dsspix */

/**
 *     DSSCRD computes the standard DSS plate coordinates and their
 *     partial derivatives for use in inverting the plate solution
 *     equations.
 *     From the AIPS program SKYVE by Mark Calebretta 
 *
 *     Called:
 *          none
 *
 *     Algorithm:
 *          The equations for computing standard DSS plate coordinates
 *          from plate offsets are given on page 11 of the booklet
 *          supplied with the Digitized Sky Survey CD:
 *
 *             xi  = a1*x + a2*y + a3 + a4*x*x + a5*x*y + a6*y*y +
 *                   a7*(x*x+y*y) + a8*x*x*x + a9*x*x*y + a10*x*y*y +
 *                   a11*y*y*y + a12*x*(x*x+y*y) + a13*x*(x*x+y*y)**2
 *
 *             eta = b1*y + b2*x + b3 + b4*y*y + b5*x*y + b6*x*x +
 *                   b7*(x*x+y*y) + b8*y*y*y + b9*x*y*y + b10*x*x*y +
 *                   b11*x*x*x + b12*y*(x*x+y*y) + b13*y*(x*x+y*y)**2
 *
 *     Notes:
 *       1) Adapted from the C function pltmodel() in the "getimage"
 *          library supplied with the Digitized Sky Survey.  Note that
 *          this routine has a bug in the computation of the derivative
 *          of ETA with respect to Y.
 *
 *     Author:
 *          Mark Calabretta, Australia Telescope.
 *          Origin; 1994/07/26  Code last modified; 1994/08/04
 *
 * \param dssx    The x coefficients for the plate solution for
 *                computing standard plate coordinates
 * \param dssy    The y coefficients for the plate solution for
 *                computing standard plate coordinates
 * \param x       x plate offset, in mm.
 * \param y       y plate offset, in mm.
 * \param xi     [out] Standard plate xi coordinate, in arcsec
 * \param eta    [out] Standard plate eta coordinate, in arcsec
 * \param dxidx  [out] Derivative of xi  with respect to x.
 * \param dxidy  [out] Derivative of xi  with respect to y.
 * \param detadx [out] Derivative of eta with respect to x.
 * \param detady [out] Derivative of eta with respect to y.
 * \param ierr   [out]  error code, 0=>OK
 */
static void dsscrd (odouble dssx[13], odouble dssy[13], odouble x, odouble y, 
		    odouble *xi, odouble *eta, odouble *dxidx, odouble *dxidy,
		    odouble *detadx, odouble *detady, olong *ierr)
{
  odouble xx, xxyy, xy, yy;

  *ierr = 0;

  /* Compute temporaries. */
  xx = x*x;
  yy = y*y;
  xy = x*y;
  xxyy = xx + yy;

  /* Compute XI. */
  *xi =     dssx[2] +
    x*(dssx[0] + x*(dssx[3] + x*dssx[7]))  +
    y*(dssx[1] + y*(dssx[5] + y*dssx[10])) +
    xy*(dssx[4] + x*dssx[8]  + y*dssx[9])  +
    xxyy*(dssx[6] + x*dssx[11] + dssx[12]*x*xxyy);
  
  /* Derivative of XI wrt X. */
  *dxidx =  dssx[0] +
    x*(2.0*(dssx[3] + dssx[6]) + 3.0*x*(dssx[7] + dssx[11])) +
    y*(dssx[4] + y*(dssx[9] + dssx[11])) +
    2.0*xy*dssx[8] + xxyy*(4.0*xx + xxyy)*dssx[12];
  
  /* Derivative of XI wrt Y. */
  *dxidy =  dssx[1] +
    y*(2.0*(dssx[5] + dssx[6]) + 3.0*y*dssx[10]) +
    x*(dssx[4] + x*dssx[8]) +
    2.0*xy*(dssx[9] + dssx[11]) + 4.0*xy*xxyy*dssx[12];
  
  /* Compute ETA. */
  *eta =    dssy[2] +
    y*(dssy[0] + y*(dssy[3] + y*dssy[7]))  +
    x*(dssy[1] + x*(dssy[5] + x*dssy[10])) +
    xy*(dssy[4] + y*dssy[8]  + x*dssy[9])  +
    xxyy*(dssy[6] + y*dssy[11] + dssy[12]*y*xxyy);
  
  /* Derivative of ETA wrt X. */
  *detadx = dssy[1] +
    x*(2.0*(dssy[5] + dssy[6]) + 3.0*x*dssy[10]) +
    y*(dssy[4] + y*dssy[8]) +
    2.0*xy*(dssy[9] + dssy[11]) + 4.0*xy*xxyy*dssy[12];
  
  /* Derivative of ETA wrt Y. */
  *detady = dssy[0] +
    y*(2.0*(dssy[3] + dssy[6]) + 3.0*y*(dssy[7] + dssy[11])) +
    x*(dssy[4] + x*(dssy[9] + dssy[11])) +
    2.0*xy*dssy[8] + xxyy*(4.0*yy + xxyy)*dssy[12];
} /* end dsscrd */

/**
 *     eqstd converts J2000.0 equatorial coordinates to standard
 *     coordinates on a Digitized Sky Survey plate.
 *     From the AIPS program SKYVE by Mark Calebretta 
 *
 *     Given:
 *          ra0,dec0    D     Plate centre J2000.0 right ascension and
 *                            declination, in degrees.
 *          ra,dec      D     Required J2000.0 right ascension and
 *                            declination, in degrees.
 *
 *     Returned:
 *          xi,eta      D     Standard plate coordinates, in arcsec.
 *          ierr        I     Error status, 0 means success.
 *
 *     Called:
 *          none
 *
 *     Algorithm:
 *          The equations for computing J2000.0 right ascension and
 *          declination from the standard coordinates are given on page
 *          11 of the booklet supplied with the Digitized Sky Survey CD
 *          and these are readily invertible.  However, note that there
 *          is a misprint in the equation for declination, the COS(RA0)
 *          term should be COS(RA-RA0).
 *
 *     Notes:
 *       1) Adapted from the C function transeqstd() in the "getimage"
 *          library supplied with the Digitized Sky Survey.
 *
 *     Author:
 *          Mark Calabretta, Australia Telescope.
 *          Origin; 1994/07/26  Code last modified; 1994/08/05
 *          Translated to c by W. D. Cotton, NRAO
 *
 * \param ra0    Plate center J2000.0 right ascension (deg)
 * \param dec0   Plate center J2000.0 declination (deg)
 * \param ra     Required J2000.0 right ascension  (deg)
 * \param dec    Required J2000.0 declination (deg)
 * \param xi     [out] Standard plate xi coordinate, in arcsec
 * \param eta    [out] Standard plate eta coordinate, in arcsec
 * \param ierr   [out]  error code, 0=>OK
 */
static void eqstd (odouble ra0, odouble dec0, odouble ra, odouble dec, 
		   odouble *xi, odouble *eta, olong *ierr)
{
  odouble cdec, cdec0, cdra, f, sdec, sdec0, sdra, z;
  odouble pi = 3.141592653589793238462643; /*     Pi. */
  odouble d2r = pi/180.0;         /* Factor to convert degrees to radians. */
  odouble r2as = 180.0*3600.0/pi; /* Factor to convert radians to arcsec */

  *ierr = 0;

  /* Cache trigonometric evaluations. */
  z = dec*d2r;
  cdec  = cos(z);
  sdec  = sin(z);
  z = dec0*d2r;
  cdec0 = cos(z);
  sdec0 = sin(z);
  z = (ra-ra0)*d2r;
  cdra = cos(z);
  sdra = sin(z);

  /* compute common factor. */
  f = r2as/(sdec*sdec0 + cdec*cdec0*cdra);

  /* Compute standard coordinates. */
  *xi  = cdec*sdra*f;
  *eta = (sdec*cdec0 - cdec*sdec0*cdra)*f;
} /* end eqstd*/

/**
 * Convert IRAF CD matrix coordinates to WCS in descriptor.
 * Only increment and rotation need be calculated here.
 * Method is from Hanish and Wells 1988 WCS draft memo (never adopted). 
 * \param desc   Image descriptor, CD matrix information expected on info member
 * \param err    ObitErr stack.
 */
static void fixIRAF (ObitImageDesc *desc, ObitErr *err)
{
  ofloat rot1, rot2, det, sdet;
  odouble cd1[2], cd2[2];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "fixIRAF";

  if (err->error) return;  /* previous error */

  /* Get parameters from header */
  ObitInfoListGet(desc->info, "CD1_1", &type, dim, &cd1[0], err);
  ObitInfoListGet(desc->info, "CD1_2", &type, dim, &cd1[1], err);
  ObitInfoListGet(desc->info, "CD2_1", &type, dim, &cd2[0], err);
  ObitInfoListGet(desc->info, "CD2_2", &type, dim, &cd2[1], err);
  if (err->error) Obit_traceback_msg (err, routine, desc->name);

  /* coordinate increments */
  desc->cdelt[desc->jlocr] = sqrt (cd1[0]*cd1[0] + cd2[0]*cd2[0]);
  desc->cdelt[desc->jlocd] = sqrt (cd1[1]*cd1[1] + cd2[1]*cd2[1]);

  /* Work out signs*/
  det = cd1[0]*cd2[1] - cd1[1]*cd2[0];
  if (det>0.0) sdet = 1.0; /* sign function */
  else {
    sdet = -1.0;
    /*   if negative, it must be RA*/
    desc->cdelt[desc->jlocr] = -desc->cdelt[desc->jlocr];
  }

  /*  rotation, average over skew */
  rot1 = 57.296 * atan2 ( sdet*cd1[1], cd2[1]);
  rot2 = 57.296 * atan2 (-sdet*cd2[0], cd1[0]);

  /* coordinate rotation */
  desc->crota[desc->jlocr] = 0.0;
  desc->crota[desc->jlocd] = -0.5 * (rot1+rot2);
} /* end of fixIRAF */

