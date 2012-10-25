/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2012                                          */
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

#include "ObitIOTableFITS.h"
#include "ObitFITS.h"
#include "ObitMem.h"
#include "ObitFile.h"
#include <errno.h>
#include "ObitUVDesc.h"
#include "ObitImageDesc.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOTableFITS.c
 * ObitIOTableFITS class function definitions.
 * This class is derived from the ObitIO class.
 *
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOTableFITS";

/** Function to obtain parent ClassInfo - ObitIO */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOTableFITSClassInfo myClassInfo = {FALSE};

/** Lookup table for converting ObitInfoList types to FITS data types 
 * should match ObitInfoType definition in ObitTypes.h */
static const olong DataTypeLookup[] =
{TBYTE, TSHORT, TINT, TINT, TLONG, TBYTE, TUSHORT, TUINT, TULONG, 
 TFLOAT, TDOUBLE, TCOMPLEX, TDBLCOMPLEX, TSTRING, TLOGICAL, 
 TBIT};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOTableFITSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOTableFITSClear (gpointer in);

/** Private: Determine next row to read */
static gboolean ObitIOTableFITSNext (ObitIOTableFITS *in, olong rowno);

/** Private: Copy other header keywords. */
void  ObitIOTableKeysOtherRead(ObitIOTableFITS *in, olong *status, 
			       ObitErr *err);

/** Private: Fix bug in cfitsio keyword parsing */
void ObitIOTableFITSFixBug(gchar *out, gchar *in, olong maxn);

/** Private: Convert an array of bits packed in an oint to a byte array */
static void ObitIOTableFITSoint2byte (olong count, oint* in, gchar* out);

/** Private: Convert an array of bits packed in a byte array to an oint */
static void ObitIOTableFITSbyte2oint (olong count, gchar *in, oint *out);

/** Private: Set Class function pointers. */
static void ObitIOTableFITSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class on the first call.
 * \param name An optional name for the object.
 * \param info if non-NULL it is used to initialize the new object.
 * \param err  ObitErr for error messages.
 * \return the new object.
 */
ObitIOTableFITS* newObitIOTableFITS (gchar *name, ObitInfoList *info,
				     ObitErr *err)
{
  ObitIOTableFITS* out;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar tempStr[201];
  gchar *routine = "newObitIOTableFITS";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitIOTableFITSClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitIOTableFITS));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOTableFITSInit((gpointer)out);

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

    /* table name */
    if(!ObitInfoListGet(info, "TabName", &type, (gint32*)dim, tempStr, err))
      Obit_traceback_val (err, routine, name, out);
    tempStr[dim[0]] = 0;  /* null terminate */
    /* fetch table name from temporary buffer, null terminate. */
    if (out->tabName) g_free(out->tabName); /* release old */
    out->tabName = g_strdup(tempStr);
    
    /* table version number */
    if(!ObitInfoListGet(info, "Ver", &type, (gint32*)dim, 
			&out->tabVer, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
  } /* end of initialize from info */

  return out;
} /* end newObitIOTableFITS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOTableFITSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOTableFITSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIOTableFITSGetClass */

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
gboolean ObitIOTableFITSSame (ObitIO *in, ObitInfoList *in1, 
			   ObitInfoList *in2, ObitErr *err)
{
  olong disk1, ver1, disk2, ver2;
  gchar *filename1, *filename2;
  gchar *type1, *type2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOTableFITSSame";

  /* error checks */
  if (err->error) return same;
  errno = 0;  /* reset any system error */

  /* get file from info */
  if(!ObitInfoListGet(in1, "Disk", &type, dim, &disk1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in1, "Ver", &type, dim, &ver1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in1, "FileName", &type, dim, 
		       (gpointer)&filename1)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }

  if(!ObitInfoListGetP(in1, "TabName", &type, dim, (gpointer)&type1)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry TabName not in InfoList Object %s",
		   routine, in->name);
  }

  if(!ObitInfoListGet(in2, "Disk", &type, dim, &disk2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "Ver", &type, dim, &ver2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in2, "FileName", &type, dim, (gpointer)&filename2)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }

  if(!ObitInfoListGetP(in2, "TabName", &type, dim, (gpointer)&type2)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry TabName not in InfoList Object %s",
		   routine, in->name);
  }


  /* Compare */
  same = (disk1==disk2) && (ver1==ver2) && 
    !strncmp (filename1,filename2, 200) &&
    !strncmp (type1, type2, 50);

  return same;
} /* end ObitIOTableFITSSame */

/**
 * Delete underlying files.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOTableFITSZap (ObitIOTableFITS *in, ObitErr *err)
{
  int status = 0;
  gboolean doClose = FALSE;
  gchar tempStr[201];

 /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if (ObitFileErrMsg(err)) return; /* Existing system error? */
  errno = 0;  /* reset any system error */

  /* Open file if needed */
  if (in->myStatus == OBIT_Inactive) {
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);
    ObitTrimTrail(tempStr);  /* Trim any trailing blanks */

    if ( fits_open_file(&(in->myFptr), tempStr, READWRITE, &status) ) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening output FITS file %s", 
		     in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      return;
    }
    doClose = TRUE;
  }
    
  /* Position to table if it exists */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, (char*)in->tabName, (int)in->tabVer, &status);
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR positioning to FITS table for %s",
		   in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    return;
  }

  /* delete the HDU */
  fits_delete_hdu (in->myFptr, NULL, &status );
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR deleting FITS table for %s",
		   in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    return;
  }

  /* Close if needed */
  if (doClose) {
    fits_close_file (in->myFptr, &status);
    if (status != 0) {
      Obit_log_error(err, OBIT_Error, "ERROR closing FITS file %s"
		     ,in->name);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      return;
    }
  }

 return;
} /* end ObitIOTableFITSZap */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGetIOTableFITSClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOTableFITSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOTableFITSClass */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOTableFITS* ObitIOTableFITSCopy  (ObitIOTableFITS *in, 
				       ObitIOTableFITS *out, ObitErr *err)
{
  const ObitIOClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitIOTableFITS(outName, NULL, err);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy this class */
  if (out->FileName!=NULL) g_free(out->FileName);
  out->FileName = g_strdup(in->FileName);
  if (out->tabName!=NULL) g_free(out->tabName);
  out->tabName = g_strdup(in->tabName);
  out->tabVer = in->tabVer;

  return out;
} /* end ObitIOTableFITSCopy */

/**
 * Initialize structures and open file.
 * The file etc. info should have been stored in the ObitInfoList.
 * The table descriptor is read if ReadOnly or ReadOnly and
 * written to disk if opened WriteOnly.
 * For accessing FITS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk"     OBIT_int (1,1,1) FITS "disk" number.
 * \li "FileName" OBIT_string (?,1,1) FITS file name.
 * \li "TabName"  OBIT_string (?,1,1) Table name (e.g. "AIPS CC").
 * \li "Ver"      OBIT_int    (1,1,1) Table version number
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableFITSOpen (ObitIOTableFITS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar tempStr[201];
  ObitTableDesc* desc;
  gchar *routine = "ObitIOTableFITSOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  g_assert (in->mySel != NULL);
  errno = 0;  /* reset any system error */

  desc = in->myDesc;

  if(!ObitInfoListGet(info, "Disk", &type, (gint32*)dim, 
		      &in->disk, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
 
 /* Get control info from InfoList */
  if(!ObitInfoListGet(info, "FileName", &type, (gint32*)dim, 
		      tempStr, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* form file name for file */
  /* fetch file name from temporary buffer, null terminate. */
  tempStr[dim[0]] = 0;
  ObitTrimTrail(tempStr);  /* Trim any trailing blanks */
  if (in->FileName) g_free(in->FileName); /* release old */
  in->FileName = ObitFITSFilename (in->disk, tempStr, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  /* table name */
  if(!ObitInfoListGet(info, "TabName", &type, (gint32*)dim, tempStr, err))
    Obit_traceback_val (err, routine, in->name, retCode);
  tempStr[dim[0]] = 0;  /* null terminate */
  /* fetch table name from temporary buffer, null terminate. */
  if (in->tabName) g_free(in->tabName); /* release old */
  in->tabName = g_strdup(tempStr); 


  /* table version number */
  if(!ObitInfoListGet(info, "Ver", &type, (gint32*)dim, 
		      &in->tabVer, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

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
    fits_open_file(&(in->myFptr), tempStr, READWRITE, &status);
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
  } else if (access == OBIT_IO_ReadWrite) {
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);
    /* Initialize output file */
    if ( fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status) ) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening output FITS file %s", 
		     in->FileName);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      retCode = OBIT_IO_OpenErr;
      return retCode;
    }

  /*------------------------ Write Only ---------------------------------*/
  } else if (access == OBIT_IO_WriteOnly) {
    /* Output file may already exist test open */
    /* must strip any leading "!" for read/write */
    if (in->FileName[0]=='!') strncpy (tempStr, (gchar*)&in->FileName[1], 200);
    else strncpy (tempStr, in->FileName, 200);

    /* Open read/write to see if it's there */
    fits_open_file(&(in->myFptr), (char*)tempStr, READWRITE, &status);
    /*fits_open_table(&(in->myFptr), tempStr, READWRITE, &status);*/

    if (status==FILE_NOT_OPENED) { /* Not there - initialize output file */
      status = 0;
      if (fits_create_file(&(in->myFptr), (char*)in->FileName, &status) ) {
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening output FITS file %s", 
		       in->FileName);
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

  /* Input table positioned in ReadDescriptor function */
  /* Output table created in WriteDescriptor function */

  } else {
    /* should never get here */
    g_assert_not_reached(); 
  }

  /* save information */
  in->access = access;
  in->myStatus = OBIT_Active;
  if (desc) desc->firstRow = 0;  

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableFITSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOTableFITSClose (ObitIOTableFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int status = 0;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  errno = 0;  /* reset any system error */

  /* don't bother if it's not open */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) 
    return OBIT_IO_OK;

  fits_close_file (in->myFptr, &status);
  if (status !=0) {
    Obit_log_error(err, OBIT_Error, "ERROR closing FITS file %s"
		   ,in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);
    retCode = OBIT_IO_CloseErr;
    return retCode;
  }

  in->myStatus = OBIT_Inactive;
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableFITSClose */

/**
 * initialize I/O - not used for FITS.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableFITSSet (ObitIOTableFITS *in, ObitInfoList *info, 
			       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  /* just return - nothing to do */
  return retCode;
} /* end ObitIOTableFITSSet */

/**
 * Read table data from disk and specifying start row.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in buffer) and the I/O has been closed.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param rowno Starting row number (1-rel) -1=> next.
 * \param data pointer to buffer to receive results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK,
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOTableFITSReadRow (ObitIOTableFITS *in, olong rowno,
				   ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int anynull, iCol, ftype, izero = 0, status = 0;
  long iRow, ndo, iElem;
  olong len, ilen, offset, nleft, i, j, ioff, nRows, row;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  gboolean done, rewrite, *outBool, last;
  gchar bitarr[32], *grumble[1], *tstr, nulstr[81], *cdata = (gchar*)data;
  ofloat fblank = ObitMagicF();
  olong  *idata = (olong*)data;
  gchar *routine = "ObitIOTableFITSReadRow";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */
  Obit_retval_if_fail ((in->myStatus==OBIT_Active)||(in->myStatus==OBIT_Modified), 
		       err, retCode,
		       "Cannot read, I/O not currently active");

  desc = in->myDesc; /* Table descriptor pointer */
  sel  = in->mySel;  /* selector pointer */
  ilen = desc->lrow / sizeof(olong); /* Size of row in gints */

  /* Has any data been modified? If so rewrite this buffer first */
  rewrite = FALSE;
  if ((desc->firstRow>0) && (desc->numRowBuff>0) && 
      (desc->statusOff>=0)) {
    ioff = 0;
    for (iRow=0; iRow<desc->numRowBuff; iRow++) {
      rewrite = rewrite || (idata[ioff+desc->statusOff]>0);
      ioff  += ilen;       /* offset in data buffer */
    }
    /* don't attempt rewrite if opened ReadOnly */
    if (desc->access==OBIT_IO_ReadOnly) rewrite = FALSE;
  }

  /* Rewrite previous buffer if needed */
  if (rewrite && (sel->nRowPIO>1)) {
    /* NB: This probably needs a better way to determine where 
       (which rows) the current buffer gets written */
    retCode = ObitIOTableFITSWrite (in, data, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  } /* end of rewrite modified data */

  /* what next ? */
  done = ObitIOTableFITSNext (in, rowno);

  /* check if done */
  if (done) {
    return OBIT_IO_EOF;
  }

  row    = desc->firstRow;            /* which row to start */
  nRows  = desc->numRowBuff;          /* How many rows to do */
  len    = desc->lrow;                /* Size of row in bytes */
  offset = 0;                         /* offset in buffer */

  sprintf (nulstr, "    "); /* Initialize */

 /* read file one row at a time */
  retCode = OBIT_IO_ReadErr; /* in case something goes wrong */
  for (iRow=row; iRow<=row+nRows-1; iRow++) {

    /* Two loops over column to keep cfitsio from clobbering other data types 
       with strings, do strings first loop */
    /* Loop over columns reading into array*/
    for (iCol = 1; iCol<desc->nfield; iCol++) {
      ftype = DataTypeLookup[desc->type[iCol-1]]; /* data type */
      /* Strings,bit arrays need to be handled differently */
      if (ftype==TSTRING) { /* string */
	grumble[0] = &cdata[offset+desc->byteOffset[iCol-1]];
	fits_read_col_str(in->myFptr, iCol, iRow, 1, 1,
			  nulstr, grumble, &anynull, &status);
	/* replace any NULLs with blank */
	tstr = grumble[0];
	last = FALSE;
	for (j=0; j<desc->repeat[iCol-1]; j++) 
	  {if (tstr[j]==0) last=TRUE; if (last) tstr[j] = ' ';}
      } 
    } /* end loop over string column */

    /* Loop over non string columns reading into array */
    for (iCol = 1; iCol<desc->nfield; iCol++) {
      ftype = DataTypeLookup[desc->type[iCol-1]]; /* data type */
      /* Strings,bit arrays need to be handled differently */
      if (ftype==TSTRING) { /* strings done last loop */
      } else if (ftype==TFLOAT) { /* float allow conversion from Nan to fblank */
	fits_read_col_flt(in->myFptr, 
			  iCol, iRow, 1, (long)desc->repeat[iCol-1], fblank, 
			  (void*)&cdata[offset+desc->byteOffset[iCol-1]], 
			  &anynull, &status); 
      } else if (ftype==TDOUBLE) { /* double */
	fits_read_col_dbl(in->myFptr,
			  iCol, iRow, 1, (long)desc->repeat[iCol-1], 0.0,
			  (void*)&cdata[offset+desc->byteOffset[iCol-1]], &anynull, &status);
      } else if (ftype==TBYTE) { /* byte */
	fits_read_col_byt(in->myFptr, 
			  iCol, iRow, 1, (long)desc->repeat[iCol-1], 0,
			  (void*)&cdata[offset+desc->byteOffset[iCol-1]], &anynull, &status);
      } else if (ftype==TSHORT) { /* short integer */
	fits_read_col_sht(in->myFptr,
			  iCol, iRow, 1, (long)desc->repeat[iCol-1], 0,
			  (void*)&cdata[offset+desc->byteOffset[iCol-1]], &anynull, &status);
      } else if (ftype==TINT) { /* int */
	fits_read_col_int(in->myFptr,
			  iCol, iRow, 1, (long)desc->repeat[iCol-1], 0,
			  (void*)&cdata[offset+desc->byteOffset[iCol-1]], &anynull, &status);
      } else if (ftype==TLONG) { /* long */
	fits_read_col_lng(in->myFptr,
			  iCol, iRow, 1, (long)desc->repeat[iCol-1], 0,
			  (void*)&cdata[offset+desc->byteOffset[iCol-1]], &anynull, &status);
      } else if (ftype==TLOGICAL) { /* boolean - stored as char in FITS */
	 /* use bitarray buffer */
	 outBool = (gboolean*)&cdata[offset+desc->byteOffset[iCol-1]];
	 nleft = desc->repeat[iCol-1];
	 iElem = 1;
	 while (nleft>0) {
	   ndo = MIN (32, nleft);
	   fits_read_col_log(in->myFptr,
			 iCol, iRow, iElem, ndo, 0, bitarr, &anynull, &status);
	   for (i=0; i<ndo; i++) { /* Copy to output as gboolean */
	     if (bitarr[i]) *outBool++ = TRUE;
	     else *outBool++ = FALSE;
	   }
	   iElem += ndo;
	   nleft -= ndo;
	 } /* end loop reading */
     } else if (ftype==TBIT) { /* bit arrays */
	fits_read_col(in->myFptr, ftype, 
		      iCol, iRow, 1, (long)desc->repeat[iCol-1], NULL, 
		      bitarr, &anynull, &status);
	/* For compatability with AIPS, bit arrays are passed as packed into 
	   oints.  cfitsio expects them as the first bit in an array of bytes.
	   The following should work for up to 32 bits */
	ObitIOTableFITSbyte2oint (MIN (desc->repeat[iCol-1], 32), bitarr,
				  (oint*)&cdata[offset+desc->byteOffset[iCol-1]]) ;
     } else { /* other */
	fits_read_col(in->myFptr, ftype, 
		      iCol, iRow, 1, (long)desc->repeat[iCol-1], NULL, 
		      &cdata[offset+desc->byteOffset[iCol-1]], 
		      &anynull, &status);
     }
		      
    } /* end loop over non string columns */
    
    if (status!=0) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR reading FITS table data for %s",
		     in->name);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      retCode = OBIT_IO_ReadErr;
      return retCode;
    }

    /* init '_status' column to 0 */
    g_memmove (&cdata[offset+desc->byteOffset[desc->nfield-1]], 
	       (gchar*)&izero, sizeof(olong));

    offset  += len;       /* offset in data buffer in bytes */
  } /* end loop reading rows */

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableFITSReadRow */

/**
 * Read table data from disk.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in buffer) and the I/O has been closed.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to receive results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK,
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOTableFITSRead (ObitIOTableFITS *in, ofloat *data, 
				ObitErr *err)
{
  /* Call bitIOTableFITSReadRow */
  return ObitIOTableFITSReadRow (in, -1, data, err);
} /* end ObitIOTableFITSRead */

/**
 * Read table data from disk applying selection and specifying start row.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param rowno Starting row number (1-rel), -1=> next.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK,
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOTableFITSReadRowSelect (ObitIOTableFITS *in, olong rowno, 
					 ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  int anynull, ftype, izero = 0, status = 0;
  olong len, ilen, offset, nleft, ndo, i, ioff, iRow, iCol, iElem, nRows, row;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  ofloat fblank = ObitMagicF();
  gboolean done, rewrite, *outBool;
  gchar bitarr[32], *grumble[1], nulstr[81], *cdata = (gchar*)data;
  olong  *idata = (olong*)data;
  gchar *routine = "ObitIOTableFITSReadRowSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */
  Obit_retval_if_fail ((in->myStatus==OBIT_Active)||(in->myStatus==OBIT_Modified), 
		       err, retCode,
		       "Cannot read, I/O not currently active");

  desc = in->myDesc; /* Table descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  ilen = desc->lrow / sizeof(olong); /* Size of row in gints */

  /* Has any data been modified? If so rewrite this buffer first */
  rewrite = FALSE;
  if ((desc->firstRow>0) && (desc->numRowBuff>0) && 
      (desc->statusOff>=0)) {
    ioff = 0;
    for (iRow=0; iRow<desc->numRowBuff; iRow++) {
      rewrite = rewrite || (idata[ioff+desc->statusOff]>0);
      ioff  += ilen;       /* offset in data buffer */
    }
    /* don't attempt rewrite if opened ReadOnly */
    if (desc->access==OBIT_IO_ReadOnly) rewrite = FALSE;
  }

  /* Rewrite previous buffer if needed */
  if (rewrite && (sel->nRowPIO>1)) {
    /* NB: This probably needs a better way to determine where 
       (which rows) the current buffer gets written */
    retCode = ObitIOTableFITSWrite (in, data, err);
    if (err->error)  /* add trace and return on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  } /* end of rewrite modified data */

  /* what next ? */
  done = ObitIOTableFITSNext (in, rowno);

  /* check if done */
  if (done) {
    return OBIT_IO_EOF;
  }

  row    = desc->firstRow;            /* which row to start */
  nRows  = desc->numRowBuff;          /* How many rows to do */
  len    = desc->lrow;                /* Size of row in bytes */
  offset = 0;                         /* offset in buffer */

 /* read file one row at a time */
  retCode = OBIT_IO_ReadErr; /* in case something goes wrong */
  for (iRow=row; iRow<=row+nRows-1; iRow++) {
    /* Loop over columns reading into array*/
    for (iCol = 1; iCol<desc->nfield; iCol++) {
      ftype = DataTypeLookup[desc->type[iCol-1]]; /* data type */
       /* Strings,bit arrays need to be handled differently */
      if (ftype==TSTRING) { /* string */
	grumble[0] = &cdata[offset+desc->byteOffset[iCol-1]];
	fits_read_col_str(in->myFptr, iCol, iRow, 1, 1,
			  nulstr, grumble, &anynull, &status);
       } else if (ftype==TFLOAT) { /* float allow conversion from Nan to fblank */
	fits_read_col(in->myFptr, ftype, 
		      iCol, iRow, 1, (long)desc->repeat[iCol-1], &fblank, 
		      &cdata[offset+desc->byteOffset[iCol-1]], &anynull, &status);
       } else if (ftype==TLOGICAL) { /* boolean - stored as char in FITS */
	 /* use bitarray buffer */
	 outBool = (gboolean*)&cdata[offset+desc->byteOffset[iCol-1]];
	 iElem = 1;
	 nleft = desc->repeat[iCol-1];
	 while (nleft>0) {
	   ndo = MIN (32, nleft);
	   fits_read_col(in->myFptr, ftype, 
			 iCol, iRow, iElem, ndo, NULL, bitarr, &anynull, &status);
	   for (i=0; i<ndo; i++) { /* Copy to output as gboolean */
	     if (bitarr[i]) *outBool++ = TRUE;
	     else *outBool++ = FALSE;
	   }
	   iElem += ndo;
	   nleft -= ndo;
	 } /* end loop reading */
       } else if (ftype==TBIT) { /* bit arrays */
	fits_read_col(in->myFptr, ftype, 
		      iCol, iRow, 1, (long)desc->repeat[iCol-1], NULL, 
		      bitarr, &anynull, &status);
	/* For compatability with AIPS, bit arrays are passed as packed into 
	   oints.  cfitsio expects them as the first bit in an array of bytes.
	   The following should work for up to 32 bits */
	ObitIOTableFITSbyte2oint (MIN (desc->repeat[iCol-1], 32), bitarr,
				  (oint*)&cdata[offset+desc->byteOffset[iCol-1]]) ;
       } else { /* other */
	 fits_read_col(in->myFptr, ftype, 
		       iCol, iRow, 1, (long)desc->repeat[iCol-1], NULL, 
		       &cdata[offset+desc->byteOffset[iCol-1]], 
		       &anynull, &status);
       }
    } /* end loop over column */
    
    if (status!=0) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR reading FITS table data for %s",
		     in->name);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      retCode = OBIT_IO_ReadErr;
      return retCode;
    }

    /* init '_status' column to 0 */
    g_memmove (&cdata[offset+desc->byteOffset[desc->nfield-1]], 
	       (gchar*)&izero, sizeof(olong));

    offset  += len;       /* offset in data buffer in bytes */
  } /* end loop reading rows */

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableFITSReadRowSelect */

/**
 * Read table data from disk applying selection.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK,
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOTableFITSReadSelect (ObitIOTableFITS *in, ofloat *data, 
				   ObitErr *err)
{
  /* Call ObitIOTableFITSSReadRowSelect */
  return ObitIOTableFITSReadRowSelect (in, -1, data, err);
} /* end  ObitIOTableFITSReadSelect */

/**
 * Write information to disk specifying starting row.
 * When OBIT_IO_EOF is returned the data has been written,
 * data in data is ignored and the I/O is closed.
 * \param in Pointer to object to be written.
 * \param rowno Row number (1-rel) to write
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOTableFITSWriteRow (ObitIOTableFITS *in, olong rowno, 
				    ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  ofloat fblank = ObitMagicF();
  gboolean  *inBool;
  int ftype, iCol, status = 0;
  long iRow, iElem, ndo;
  olong len, offset, nRows, nleft, i, row;
  gchar bitarr[32], *grumble[1], *cdata = (gchar*)data;

  /* error checks */
  g_assert (ObitErrIsA(err));
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

  desc = in->myDesc; /* Table descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* which row to start , specified or highest? */
  if (rowno>0) desc->firstRow = rowno;
  else desc->firstRow = MAX (1, desc->nrow+1);
  row    = desc->firstRow;  
  nRows  = desc->numRowBuff;          /* How many rows to do */
  len    = desc->lrow;                /* Size of row in bytes */
  offset = 0;                         /* offset in buffer */

 /* write file one row at a time */
  retCode = OBIT_IO_ReadErr; /* in case something goes wrong */
  for (iRow=row; iRow<=row+nRows-1; iRow++) {
    /* Loop over columns reading into array*/
    for (iCol = 1; iCol<desc->nfield; iCol++) {
      if (desc->repeat[iCol-1]>0) { /* cfitsio doesn't handle no data */
	ftype = DataTypeLookup[desc->type[iCol-1]]; /* data type */
	if (ftype==TSTRING) { /* string */
	  grumble[0] = &cdata[offset+desc->byteOffset[iCol-1]];
	  fits_write_col_str(in->myFptr, iCol, iRow, 1, 1,
			     grumble, &status);
	} else if (ftype==TBYTE) { /* bytes */
	  fits_write_col_byt(in->myFptr, iCol, iRow, 1,
			 (long)desc->repeat[iCol-1], (void*)&cdata[offset+desc->byteOffset[iCol-1]],
			  &status);
	} else if (ftype==TINT) { /* int */
	  fits_write_col_int(in->myFptr, iCol, iRow, 1,
			 (long)desc->repeat[iCol-1], (void*)&cdata[offset+desc->byteOffset[iCol-1]],
			  &status);
	} else if (ftype==TSHORT) { /* short integer */
	  fits_write_col_sht(in->myFptr, iCol, iRow, 1,
			 (long)desc->repeat[iCol-1], (void*)&cdata[offset+desc->byteOffset[iCol-1]],
			  &status);
	} else if (ftype==TLONG) { /* long integer */
	  fits_write_col_lng(in->myFptr, iCol, iRow, 1,
			 (long)desc->repeat[iCol-1], (void*)&cdata[offset+desc->byteOffset[iCol-1]],
			  &status);
	} else if (ftype==TFLOAT) { /* float */
	  /* Allow conversion from fblank to NaN */
	  fits_write_colnull_flt(in->myFptr, iCol, iRow, 1,
				 (long)desc->repeat[iCol-1], (void*)&cdata[offset+desc->byteOffset[iCol-1]],
				 fblank, &status);
	} else if (ftype==TDOUBLE) { /* double */
	  fits_write_col_dbl(in->myFptr, iCol, iRow, 1,
			     (long)desc->repeat[iCol-1], (void*)&cdata[offset+desc->byteOffset[iCol-1]],
			     &status);
       } else if (ftype==TLOGICAL) { /* boolean - stored as char in FITS */
	 /* use bitarray buffer */
	 inBool = (gboolean*)&cdata[offset+desc->byteOffset[iCol-1]];
	 nleft = desc->repeat[iCol-1];
	 iElem = 1;
	 while (nleft>0) {
	   ndo = MIN (32, nleft);
	   for (i=0; i<ndo; i++) { /* Copy to bitarray as char */
	     if (inBool[i]) bitarr[i] = 1;
	     else bitarr[i] = 0;
	   }
	   fits_write_col_log(in->myFptr, iCol, iRow, iElem, ndo, bitarr, &status);
	   iElem += ndo;
	   nleft -= ndo;
	 } /* end loop writing */
	} else if (ftype==TBIT) { /* bit arrays */
	  /* For compatability with AIPS, bit arrays are passed as packed into 
	     oints.  cfitsio expects them as the first bit in an array of bytes.
	     The following should work for up to 32 bits */
	  ObitIOTableFITSoint2byte (MIN (desc->repeat[iCol-1], 32),
				    (oint*)&cdata[offset+desc->byteOffset[iCol-1]], 
				    bitarr);
	  fits_write_col(in->myFptr, ftype, iCol, iRow, 1,
			 (long)desc->repeat[iCol-1], bitarr,
			 &status);
	} else { /* other */
	  fits_write_col(in->myFptr, ftype, iCol, iRow, 1,
			 (long)desc->repeat[iCol-1],
			 &cdata[offset+desc->byteOffset[iCol-1]],
			 &status);
	}
      }
    } /* end loop over column */
    
    if (status!=0) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR %d writing FITS table data for %s",
		     status, in->name);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      ObitFileErrMsg(err);     /* system error message*/
      retCode = OBIT_IO_ReadErr;
      return retCode;
    }

    offset  += len;       /* offset in data buffer in bytes */
  } /* end loop reading rows */

  /* where will the next write start */
   desc->firstRow += desc->numRowBuff;
   if (rowno>0) desc->firstRow = rowno; /* Override sequence */

  in->myStatus = OBIT_Modified; /* file has been modified */

  return OBIT_IO_OK;
} /* end ObitIOTableFITSWriteRow */

/**
 * Write information to disk.
 * When OBIT_IO_EOF is returned the data has been written,
 * data in data is ignored and the I/O is closed.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0(OBIT_IO_OK)=> OK
 *          OBIT_IO_EOF => image finished.
 */
ObitIOCode ObitIOTableFITSWrite (ObitIOTableFITS *in, ofloat *data, 
				 ObitErr *err)
{
  /* Call ObitIOTableFITSWriteRow */
  return ObitIOTableFITSWriteRow (in, -1, data, err);
} /* end ObitIOTableFITSWrite */

/**
 * Read table Descriptor data from disk.
 * If the table version number is 0, then the highest numbered 
 * table of the same name is used.
 * \param in Pointer to object with ObitTableDesc to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableFITSReadDescriptor (ObitIOTableFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar commnt[FLEN_COMMENT];
  gchar typechar[4], cdata[FLEN_CARD], unit[FLEN_CARD];
  int i, j, nhdu, hdutype, naxis, ncol, status = 0;
  olong nfield;
  long  extver, ltemp, repeat, sorto, naxes[10];
  double scale, zero;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  gchar *routine = "ObitIOTableFITSReadDescriptor";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);
  errno = 0;  /* reset any system error */

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* Position to table */
  /* if the version number is 0 then find the highest numbered one */
  if (in->tabVer<=0) {
    fits_get_num_hdus (in->myFptr, &nhdu, &status); /* how many? */
    for (i=1; i<=nhdu; i++) {
      fits_movabs_hdu (in->myFptr, i, &hdutype, &status);
      if (hdutype==BINARY_TBL) { /* If it's a table enter it in the list */
	/* table name */
	fits_read_key_str (in->myFptr, "EXTNAME", (char*)cdata, (char*)commnt, &status);
	/* version number */
	if (fits_read_key_lng (in->myFptr, "EXTVER", &extver, commnt, &status)) {
	  /* Default to 1 */
	  extver = 1;
	  /* Clear cfitsio */
	  fits_clear_errmsg();
	  status = 0;
	}
	/* if this is a match get the version if it's highest */
	if (!strncmp (cdata, in->tabName, FLEN_CARD-1)) 
	  in->tabVer = MAX ((olong)extver, in->tabVer);
      }
    } /* end loop over extensions */
  } /* End find highest version */

  /* Find it? */
  if (in->tabVer<=0) {
    /* If read/write just return and how the write header fixes it */
    if (in->access==OBIT_IO_ReadWrite) return OBIT_IO_OK;
    Obit_log_error(err, OBIT_Error, 
		   "Found NO tables of type %s for %s",
		   in->tabName, in->name);
    return retCode;
  }

  /* Position to table */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, in->tabName, (int)in->tabVer, &status);
  if (status!=0) {
   /* If read/write just return and how the write header fixes it */
    if ((in->access==OBIT_IO_ReadWrite) && (status==BAD_HDU_NUM)) return OBIT_IO_OK;
    Obit_log_error(err, OBIT_Error, 
		   "Found NO tables of type %s ver %d for %s",
		   in->tabName, in->tabVer, in->name);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    return OBIT_IO_OpenErr;
  }

  fits_get_num_rows (in->myFptr, &ltemp, &status);
  desc->nrow = (olong)ltemp;

  /* How many columns in table? */
  fits_get_num_cols (in->myFptr, &ncol, &status);
  nfield = ncol+1; /* Will add "_status" column */
  
  /* Table name */
  if (desc->TableName) g_free (desc->TableName);
  desc->TableName = g_strdup(in->tabName);

  /* If number of columns has changed, reallocate */
  if (nfield != desc->nfield) ObitTableDescRealloc (desc, nfield);

  /* Read keyword values, where possible */
  /* loop over columns */
  for (i=0; i<ncol; i++) {
    /* Get column basic information */
    fits_get_bcolparms (in->myFptr, i+1, 
			(char*)cdata, (char*)unit, (char*)typechar, &repeat, &scale, &zero,
			NULL, NULL, &status);
    /* fix string bug in cfitsio */
    ObitIOTableFITSFixBug(cdata, cdata, FLEN_CARD);
    ObitIOTableFITSFixBug(unit,  unit, FLEN_CARD);
    if (desc->FieldName[i]) g_free(desc->FieldName[i]);
    desc->FieldName[i] = g_strdup (cdata); /* save name */
    if (desc->FieldUnit[i]) g_free(desc->FieldUnit[i]);
    desc->FieldUnit[i] = g_strdup (unit); 
    desc->firstRow = 0;  /* where reading in table */

    /* Get dimensionality array */
    fits_read_tdim (in->myFptr, i+1, MAXINFOELEMDIM, &naxis,
		    naxes, &status);
    for (j=0; j<MAXINFOELEMDIM; j++) desc->dim[i][j] = 1;
    for (j=0; j<naxis; j++) desc->dim[i][j] = (olong)naxes[j];

    /* convert cfitsio type to InfoList type */
    if (typechar[0]=='B') {
      desc->type[i] = OBIT_ubyte;
    } else if (typechar[0]=='I') {
      desc->type[i] = OBIT_short;
    } else if (typechar[0]=='J') {
      desc->type[i] = OBIT_oint;
    } else if (typechar[0]=='A') {
      desc->type[i] = OBIT_string;
    } else if (typechar[0]=='E') {
      desc->type[i] = OBIT_float;
    } else if (typechar[0]=='D') {
      desc->type[i] = OBIT_double;
    } else if (typechar[0]=='C') {
      desc->type[i] = OBIT_complex;
    } else if (typechar[0]=='L') {
      desc->type[i] = OBIT_bool;
    } else if (typechar[0]=='X') {
      desc->type[i] = OBIT_bits;
    }
  } /* end loop reading field info */

  /* Add status column */
  if (desc->FieldName[desc->nfield-1]) g_free(desc->FieldName[desc->nfield-1]);
  desc->FieldName[desc->nfield-1] = g_strdup("_status");
  if (desc->FieldUnit[desc->nfield-1]) g_free(desc->FieldUnit[desc->nfield-1]);
  desc->FieldUnit[desc->nfield-1] = g_strdup("       ");
  desc->type[desc->nfield-1] = OBIT_long;
  desc->dim[desc->nfield-1][0] = 1; desc->dim[desc->nfield-1][1] = 1;
  desc->dim[desc->nfield-1][2] = 1; desc->dim[desc->nfield-1][3] = 1;

  /* sort order  */
   fits_read_key_lng (in->myFptr, "ISORTORD", &sorto, (char*)commnt, &status);
   if (status ==KEY_NO_EXIST) {
     desc->sort[0] = 0;
     desc->sort[1] = 0;
     status = 0;
   } else {
     desc->sort[0] = (olong)sorto;
     desc->sort[1] = 0;
   }

  /* Look for anything else and add it to the InfoList on desc */
  ObitIOTableKeysOtherRead(in, &status, err);
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* was there an error? */
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR reading input FITS file header");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_ReadErr;
    return retCode;
  }

  /* enforce defaults */
  ObitTableSelDefault(desc, sel);

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableFITSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * If the table version number is 0, then the highest numbered 
 * table of the same name +1 is used.
 * \param in Pointer to object with ObitTableDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitIOTableFITSWriteDescriptor (ObitIOTableFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar dtype=' ', **tform, cdata[FLEN_CARD], commnt[FLEN_COMMENT+1];
  int i,  nhdu, hdutype, tfield, ndata, nkey, status = 0;
  long extver, nrows, naxes[MAXINFOELEMDIM];
  olong j, k;
  gchar keyName[FLEN_KEYWORD+1], *keyNameP;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean doFill;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  ObitInfoType keyType;
  union blobEquiv {
    gchar    s[201];
    odouble  d;
    ofloat   f;
    gboolean b;
    oint     o;
    olong    i;
  } blob;
  gchar *routine = "ObitIOTableFITSWriteDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  if (desc->nfield<=0) {
    Obit_log_error(err, OBIT_Error, "%s: Table  %s has no fields defined",
		   routine, in->name);
    return retCode;
  }

  /* enforce defaults */
  ObitTableSelDefault(desc, sel);

  /* Position to table */
  /* if the version number is 0 then find the highest numbered one */
  if (in->tabVer<=0) {
    fits_get_num_hdus (in->myFptr, &nhdu, &status); /* how many? */
    for (i=1; i<=nhdu; i++) {
      fits_movabs_hdu (in->myFptr, i, &hdutype, &status);
      if (hdutype==BINARY_TBL) { /* If it's a table enter it in the list */
	/* table name */
	fits_read_key_str (in->myFptr, "EXTNAME", (char*)cdata, (char*)commnt, &status);
	/* version number */
	fits_read_key_lng (in->myFptr, "EXTVER", &extver, commnt, &status);
	/* if this is a match get the version if it's highest */
	if (!strncmp (cdata, in->tabName, FLEN_CARD-1)) 
	  in->tabVer = MAX ((olong)extver, in->tabVer);
      }
    } /* end loop over extensions */

    /* Want to create a new one */
    in->tabVer++;
  } /* End find highest version */

  /* Position to table if it exists */
  fits_movnam_hdu (in->myFptr, BINARY_TBL, in->tabName, (int)in->tabVer, &status);
  /* Create if it doesn't already exist */
  if (status == BAD_HDU_NUM) {
    status = 0;

    /* fill descriptive arrays */
    tfield = desc->nfield-1;  /* number of columns - drop _status (last) */
    /* allocate arrays */
    tform = g_malloc0(tfield*sizeof(gchar*));
    for (i=0; i<tfield; i++) {
      /* count elements */
      ndata = desc->dim[i][0];  /* First may be zero */
      for (j=1; j<MAXINFOELEMDIM; j++) ndata *= MAX (1, desc->dim[i][j]);
      /* data type is as ObitInfoType */
      if (desc->type[i]==OBIT_ubyte) {
	dtype = 'B';
      } else if (desc->type[i]==OBIT_short) {
 	dtype = 'I';
      } else if (desc->type[i]==OBIT_oint) {
	dtype = 'J';
      } else if (desc->type[i]==OBIT_string) {
 	dtype = 'A';
      } else if (desc->type[i]==OBIT_float) {
	dtype = 'E';
      } else if (desc->type[i]==OBIT_double) {
	dtype = 'D';
      } else if (desc->type[i]==OBIT_complex) {
	dtype = 'C';
      } else if (desc->type[i]==OBIT_bool) {
	dtype = 'L';
      } else if (desc->type[i]==OBIT_bits) {
	dtype = 'X';
	 /* ??? ndata = 1+ (ndata-1) / 8;packed into bytes */
      }
      g_snprintf (cdata, FLEN_CARD-1, "%d%c", ndata, dtype);
      tform[i] = g_strdup(cdata);
    } /* end defining columns */

    /* create table */
    fits_create_tbl (in->myFptr, BINARY_TBL, (long)desc->nrow,
		     tfield, desc->FieldName, tform, desc->FieldUnit,
		     (char*)in->tabName, &status);
    if (status!=0) {
      Obit_cfitsio_error(err); 
      ObitFileErrMsg(err);     /* system error message*/
      /* delete format array */
      for (i=0; i<tfield; i++) if (tform[i]) g_free(tform[i]);
      g_free(tform); tform = NULL;

      return OBIT_IO_WriteErr;
    }

    /* what version am I? */
    extver = (long)in->tabVer;
    strncpy (commnt, "Table version number", FLEN_COMMENT);
    fits_update_key_lng (in->myFptr, "EXTVER", extver, commnt, &status);

    /* delete format array */
    for (i=0; i<tfield; i++) if (tform[i]) g_free(tform[i]);
    g_free(tform); tform = NULL;

    /* set table dimensionality arrays as TDIM */
    for (i=0; i<tfield; i++) {
      if (desc->dim[i][1] > 1) {/* multi dimensional array? */
	/* count dimensions */
	ndata = 0;
	for (j=0; j<MAXINFOELEMDIM; j++) {
	  if (desc->dim[i][j]>1) ndata = j+1;
	  naxes[j] = (long)desc->dim[i][j];
	}
	fits_write_tdim (in->myFptr, i+1, ndata, naxes, &status);
	if (status!=0) {
	  Obit_cfitsio_error(err); 
	  ObitFileErrMsg(err);     /* system error message*/
	}
      }
    }
  }

  /* sort order  */
  strncpy (commnt, "AIPSish primary sort key column code", FLEN_COMMENT);
  fits_update_key_lng (in->myFptr, "ISORTORD", (long)desc->sort[0], 
		      (char*)commnt, &status);

  /* Number of rows NAXIS2  */
  nrows = desc->nrow;
  fits_get_num_rows (in->myFptr, &nrows, &status);
  if (nrows > desc->nrow) {  /* Truncate */
    nrows = nrows - desc->nrow;
    fits_delete_rows (in->myFptr, desc->nrow+1, nrows, &status);
  }


  /* Write other keywords from descriptor */
  if (desc->info) nkey = desc->info->number+0; /* How many keywords? */
  else nkey = 0;
  retCode = OBIT_IO_WriteErr;
  strncpy (commnt, "             ", FLEN_COMMENT);
  for (i=1; i<=nkey; i++) {
    /* Copy from ObitInfoList */
    ObitInfoListGetNumber(desc->info, (olong)i, &keyNameP, &keyType, dim, 
			  blob.s, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
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

  /* an error? */
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, "ERROR writing FITS file header");
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    ObitFileErrMsg(err);     /* system error message*/
    retCode = OBIT_IO_WriteErr;
    return retCode;
  }

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableFITSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableFITSFlush (ObitIOTableFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* cfitsio does the buffer flushing on close */

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableFITSFlush */

/**
 * Create buffer approptiate for I/O request.
 * Should be called after ObitIO is opened.
 * \param data (output) pointer to data array
 * \param size (output) size of data array in bytes.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 */
void 
ObitIOTableFITSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOTableFITS *in, ObitInfoList *info, 
			     ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  Obit_return_if_fail(((in->myStatus==OBIT_Modified) ||
		       (in->myStatus==OBIT_Active)), 
		      err,
		      "Cannot define buffer, I/O not currently active");

  /* get size */
  *size = ObitTableSelBufferSize(in->myDesc, in->mySel);

  /* (re)allocate */
  if (*data) *data = ObitMemRealloc (*data, (*size)*sizeof(ofloat));
  else *data = ObitMemAlloc0Name((*size)*sizeof(ofloat), "TableBuffer");

} /* end ObitIOTableFITSCreateBuffer */

/**
 * Get underlying file information in entries to an ObitInfoList
 * Following entries for AIPS files ("xxx" = prefix):
 * \param in      Object of interest.
 * \param myInfo  InfoList on basic object with selection
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 * \param outList InfoList to write entries into
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFileName OBIT_string  FITS file name
 * \li xxxDisk     OBIT_oint    FITS file disk number
 * \li xxxDir      OBIT_string  Directory name for xxxDisk
 *
 * For all File types types:
 * \li xxxDataType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxFileType OBIT_oint File type as "FITS"
 *    
 * For xxxDataType = "Table"
 * \li xxxTableParent OBIT_string  Table parent type (e.g. "MA", "UV")
 * \li xxxTab   OBIT_string  (Tables only) Table type (e.g. "AIPS CC")
 * \li xxxVer   OBIT_oint    (Tables Only) Table version number
 *    
 * \param err     ObitErr for reporting errors.
 */
void ObitIOTableFITSGetFileInfo (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
				 ObitInfoList *outList, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *FileType="FITS", *dirname;
  gchar *DataType[]={"Table","UV","MA"};
  gchar tempStr[201], *here="./";
  gpointer listPnt;
  olong disk, ver, i, ptype;
  gchar *parm[] = {NULL};
  gchar *routine = "ObitIOTableFITSGetFileInfo";

  if (err->error) return;

  /* Set basic information */
  /* FITS */
  if (prefix) keyword =  g_strconcat (prefix, "FileType", NULL);
  else keyword =  g_strdup ("FileType");
  dim[0] = strlen(FileType);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, FileType);
  g_free(keyword);
  
   /* type - Table */
  if (prefix) keyword =  g_strconcat (prefix, "DataType", NULL);
  else keyword =  g_strdup ("DataType");
  dim[0] = strlen(DataType[0]);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, DataType[0]);
  g_free(keyword);
  
  /* parent type */
  if (ObitUVDescIsA(in->myDesc))    ptype = 1;
  if (ObitImageDescIsA(in->myDesc)) ptype = 2;
  if (prefix) keyword =  g_strconcat (prefix, "TableParent", NULL);
  else keyword =  g_strdup ("TableParent");
  dim[0] = strlen(DataType[ptype]);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, DataType[ptype]);
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

  /* Table type */
  if (!ObitInfoListGet(myInfo, "TableType", &type, dim, tempStr, err)) 
      Obit_traceback_msg (err, routine, in->name);
  if (prefix) keyword =  g_strconcat (prefix, "Tab", NULL);
  else keyword =  g_strdup ("Tab");
  dim[0] = strlen (tempStr);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, tempStr);
  g_free(keyword);

  /* Table Version */
  if (!ObitInfoListGet(myInfo, "Ver", &type, dim, &ver, err)) 
      Obit_traceback_msg (err, routine, in->name);
  if (prefix) keyword =  g_strconcat (prefix, "Ver", NULL);
  else keyword =  g_strdup ("Ver");
  dim[0] = 1;
  ObitInfoListAlwaysPut (outList, keyword, OBIT_long, dim, &ver);
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
} /* end ObitIOTableFITSGetFileInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOTableFITSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOTableFITSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOTableFITSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOTableFITSClassInfoDefFn (gpointer inClass)
{
  ObitIOTableFITSClassInfo *theClass = (ObitIOTableFITSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOTableFITSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOTableFITSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOTableFITSGetClass;
  theClass->newObit       = NULL;
  theClass->newObitIO     = (newObitIOFP)newObitIOTableFITS;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOTableFITSSame;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOTableFITSZap;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOTableFITSCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = 
    (ObitClearFP)ObitIOTableFITSClear;
  theClass->ObitInit      = 
    (ObitInitFP)ObitIOTableFITSInit;
  theClass->ObitIOOpen    = 
    (ObitIOOpenFP)ObitIOTableFITSOpen;
  theClass->ObitIOClose   = 
    (ObitIOCloseFP)ObitIOTableFITSClose;
  theClass->ObitIOSet     = 
    (ObitIOSetFP)ObitIOTableFITSSet;
  theClass->ObitIORead    = 
    (ObitIOReadFP)ObitIOTableFITSRead;
  theClass->ObitIOReadRow = 
    (ObitIOReadRowFP)ObitIOTableFITSReadRow;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOTableFITSReadSelect;
  theClass->ObitIOReadRowSelect = 
    (ObitIOReadRowSelectFP)ObitIOTableFITSReadRowSelect;
  theClass->ObitIOWrite   = 
    (ObitIOWriteFP)ObitIOTableFITSWrite;
  theClass->ObitIOWriteRow   = 
    (ObitIOWriteRowFP)ObitIOTableFITSWriteRow;
  theClass->ObitIOFlush   = 
    (ObitIOFlushFP)ObitIOTableFITSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOTableFITSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOTableFITSWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOTableFITSCreateBuffer;
  theClass->ObitIOFreeBuffer   = 
    (ObitIOFreeBufferFP)ObitIOFreeBuffer;
  theClass->ObitIOGetFileInfo   =
    (ObitIOGetFileInfoFP)ObitIOTableFITSGetFileInfo;

} /* end ObitIOTableFITSClassDefFn */


/*--------------- Private functions --------------------------*/

/**
 * Creates empty member objects.
 * for each parent class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOTableFITSInit  (gpointer inn)
{
  const ObitClassInfo *ParentClass;
  ObitIOTableFITS *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && (ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->FileName     = NULL;
  in->tabName      = NULL;
  in->tabVer       = -1;

} /* end ObitIOTableFITSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOTableFITSClear (gpointer inn)
{
  ObitIOTableFITS *in = inn;
  const ObitClassInfo *ParentClass;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->myStatus==OBIT_Active) ||(in->myStatus==OBIT_Modified)) {
    err = newObitErr();
    ObitIOTableFITSClose (in, err); 
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  if (in->FileName) g_free(in->FileName);
  in->FileName = NULL;
  if (in->tabName) g_free(in->tabName);
  in->tabName = NULL;

 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && (ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOTableFITSClear */

/**
 * Uses selector member to decide which visibilities to
 * read next.
 * Leaves values in myDesc as firstVis and numVisBuff.
 * \param  in Pointer to the object.
  * \param rowno Starting row number (1-rel), -1=> next;
* \return TRUE if all data read, else FALSE
 */
static gboolean ObitIOTableFITSNext (ObitIOTableFITS *in, olong rowno)
{
  ObitTableDesc* desc;
  ObitTableSel* sel;
  olong nleft;
  gboolean done = FALSE;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
 
  desc = in->myDesc; /* Table descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* Update beginning - is this the first  pass? */
  if (desc->firstRow<=0) {
    desc->firstRow = 1;
  } else {
    desc->firstRow += desc->numRowBuff; /* These should have been done */
  }
  desc->numRowBuff = sel->nRowPIO; /* number to do next time */

   /* Check is particular row specified */
  if (rowno>0) desc->firstRow = rowno;

 /* but not more than all */
  nleft = desc->nrow - desc->firstRow + 1;
  desc->numRowBuff = MAX (1,  desc->numRowBuff);
  desc->numRowBuff = MIN (nleft,  desc->numRowBuff);
  done = (nleft<=0);
  
  return done;
} /* end ObitIOTableFITSNext */

/**
 * Look for additional descriptive keywords, any that are not 
 * on the exclusion list are copied to the descriptor InfoList.
 * \param in      Pointer to ObitIOTableFITS.
 * \param lstatus (Output) cfitsio status.
 * \param err    ObitErr stack.
 * \return return code, 0=> OK
 */
void  ObitIOTableKeysOtherRead(ObitIOTableFITS *in, olong *lstatus, 
			       ObitErr *err)
{
  gchar keywrd[FLEN_KEYWORD], value[FLEN_VALUE], commnt[FLEN_COMMENT+1];
  gchar *first, *last, *anF, *aT, dtype, svalue[FLEN_VALUE];
  int i, j, k, l, keys, morekeys, status = (int)*lstatus;
  olong ivalue;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  double dvalue;
  ObitTableDesc *desc;
  gchar *exclude[] = 
  {"SIMPLE", "BITPIX", "EXTEND", "HISTORY", "COMMENT", "BLANK", "        ",
   "XTENSION", "PCOUNT", "GCOUNT", "EXTNAME", "EXTVER", "ISORTORD",
   "TFORM", "TTYPE", "TUNIT", "TSCAL", "TZERO", "TDIM", "TDISP",
   "1CTYP", "2CTYP", "3CTYP", "4CTYP", "5CTYP", "6CTYP", "7CTYP", 
   "1CRVL", "2CRVL", "3CRVL", "4CRVL", "5CRVL", "6CRVL", "7CRVL", 
   "1CDLT", "2CDLT", "3CDLT", "4CDLT", "5CDLT", "6CDLT", "7CDLT", 
   "1CRPX", "2CRPX", "3CRPX", "4CRPX", "5CRPX", "6CRPX", "7CRPX", 
   "1CROT", "2CROT", "3CROT", "4CROT", "5CROT", "6CROT", "7CROT", 
   "BSCALE", "BZERO", "NAXIS", "TFIELDS", NULL};
  olong number, *len=NULL;
  gboolean bvalue, bad=FALSE;
  gchar *routine = "ObitIOTableKeysOtherRead";

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete old InfoList and restart */
  ((ObitTableDesc*)in->myDesc)->info = ObitInfoListUnref (((ObitTableDesc*)in->myDesc)->info);
  ((ObitTableDesc*)in->myDesc)->info = (gpointer)newObitInfoList ();
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
	if (err->error)  /* add trace and return on error */
	  Obit_traceback_msg (err, routine, in->name);
      }
    }
  } /* end loop over keywords */

  /* Cleanup */
  if (len) g_free(len);
  *lstatus = (olong)status;

} /* end ObitIOTableKeysOtherRead */

/**
 * Work around for cfitsio bug, trailing blanks in keywords are dropped.
 * This routine blank fills out and null terminates at out[maxn-1].
 * \param out    Output string 
 * \param in     Input string from cfitsio keyword read
 * \param maxn   length of out
 */
void ObitIOTableFITSFixBug (gchar *out, gchar *in, olong maxn)
{
  olong i, len;

  len = strlen(in);
  for (i=0; i<len; i++)    out[i] = in[i];
  for (i=len; i<maxn; i++) out[i] = ' ';
  out[maxn-1] = 0;
  
} /* end ObitIOTableFITSFixBug */

/**
 * Convert bit array packed into an oint into the first bits of a byte array.
 * \param count  Number of bits
 * \param in     Input bit array packed into an oint
 * \param out    Output array with one bit per byte
 */
static void ObitIOTableFITSoint2byte (olong count, oint* in, gchar *out)
{
  olong i, j, bit;
  oint  current, mask;

  /* Get first oint to convert */
  j = 0;
  bit = 0;
  current = in[j++];

  /* Loop over array */
  for (i=0; i<count; i++) {
    
    /* Time for another oint? */
    if (bit>=32) {
      bit = 0;
      current = in[j++];    
    }
    
    mask = 1<<bit; /* Mask for bit */
    out[i] = current & mask;
    bit++;  /* next bit */
  } /* end loop over array */
} /* end ObitIOTableFITSoint2byte */

/**
 * Convert bit array packed into the first bits of a byte array to bits in an oint.
 * \param count  Number of bits
 * \param in     Input bit array with one bit per byte
 * \param out    Output array packed into an oint
 */
static void ObitIOTableFITSbyte2oint (olong count, gchar *in, oint  *out)
{
  olong i, j, bit;
  oint  *current, mask;

  /* Get first oint to convert */
  j = 0;
  bit = 0;
  current = &out[j++];
  *current = 0;

  /* Loop over array */
  for (i=0; i<count; i++) {
    
    /* Time for another oint? */
    if (bit>=32) {
      bit = 0;
      current = &out[j++];    
   *current = 0;
   }
    
    if (in[i]) mask = 1<<bit; /* Mask for bit */
    else mask = 0;
    *current = *current | mask; /* turn on bit if selected */
    bit++;  /* next bit */
  } /* end loop over array */
} /* end ObitIOTableFITSbyte2oint */
