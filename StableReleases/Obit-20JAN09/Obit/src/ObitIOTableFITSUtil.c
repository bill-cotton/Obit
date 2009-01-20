/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006,2008                                          */
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

#include "ObitIOTableFITSUtil.h"
#include "ObitFITS.h"
#include "fitsio.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOTableFITSUtil.c
 * ObitIOTableFITSUtil -  FITS table utilities
 */

/*----------------------Public functions---------------------------*/

/**
 * Read Table list from an FITS catalog entry
 * \param in  Pointer to object to be read.
 * \param err ObitErr for reporting errors.
 */
void ObitIOTableFITSUtilReadTableList(ObitData *in, ObitErr *err)
{
  ObitTableList* tableList;
  int i, disk, nhdu, hdutype, status = 0;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  gchar commnt[FLEN_COMMENT], tempStr[201];
  gchar cdata[FLEN_CARD], *FileName=NULL;
  /** cfitsio file pointer */
  fitsfile *myFptr;
  long extver;
  olong otemp;
  gchar *routine = "ObitIOTableFITSUtilReadTableList";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitDataIsA(in));

  /* File info */
  if (!ObitInfoListGet(in->info, "Disk", &type, dim,  &disk, err)) 
    Obit_traceback_msg (err, routine, in->name);
 
  if (!ObitInfoListGet(in->info, "FileName", &type, dim, tempStr, err))
    Obit_traceback_msg (err, routine, in->name);

  /* form file name for file - fetch file name from temporary buffer, 
     null terminate. */ 
  tempStr[dim[0]] = 0;
  FileName = ObitFITSFilename (disk, tempStr, err);
  if (err->error) Obit_traceback_msg (err, routine, FileName);

  /* must strip any leading "!" for read/write */
  if (FileName[0]=='!') strncpy (tempStr, (gchar*)&FileName[1], 200);
  else strncpy (tempStr, FileName, 200);
  if (FileName) g_free(FileName); 
  ObitTrimTrail(tempStr);  /* Trim any trailing blanks */

  /* cfitsio refuses to open a file readwrite after it has been opened
     readonly so first try opening readwrite
     If that fails try readonly. */
  /* Test open readwrite */
  fits_open_file(&(myFptr), tempStr, READWRITE, &status);
  if ((status==FILE_NOT_OPENED) || (status==READONLY_FILE)) {
    /* Failed - try readonly */
    status = 0;
    fits_clear_errmsg();   /* Clear cfitsio error stack */
    if (fits_open_file(&(myFptr), tempStr, READONLY, &status) ) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR %d opening input FITS file %s", status, tempStr);
      Obit_cfitsio_error(err); /* copy cfitsio error stack */
      return;
    }
  }
  /* Other errors? */
  if (status!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR opening output FITS file %s", tempStr);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    return;
  }
  
  tableList = (ObitTableList*)in->tableList;
  
  /* Index tables in file and update TableList if not already done*/
  if (tableList->number <= 0) {
    fits_get_num_hdus (myFptr, &nhdu, &status); /* how many? */
    for (i=1; i<=nhdu; i++) {
      fits_movabs_hdu (myFptr, i, &hdutype, &status);
      if (hdutype==BINARY_TBL) { /* If it's a table enter it in the list */
	/* table name */
	fits_read_key_str (myFptr, "EXTNAME", (char*)cdata, (char*)commnt, &status);
	/* version number */
	fits_read_key_lng (myFptr, "EXTVER", &extver, (char*)commnt, &status);
	if (status==0) { /* Add to TableList unless it's uv or OTF data */
	  if (strcmp (cdata, "AIPS UV") && strcmp (cdata, "OTF")) {
	    otemp = (olong)extver;
	    ObitTableListPut (tableList, cdata, &otemp, NULL, err);
	    if (err->error) goto cleanup;
	  }
	}
      }
    } /* end loop indexing file */
  } /* end update Table List */

  /* Cleanup */
 cleanup:
  fits_close_file (myFptr, &status);
  /* cfitsio error? */
  if (status) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: ERROR %d reading table list in %s", 
		   routine, status, tempStr);
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
  }
  if (err->error) Obit_traceback_msg (err, routine, tableList->name);

} /* End of ObitIOTableFITSUtilReadTableList */

/**
 * Write Table list to an FITS catalog entry
 * Not really needed for FITS files
 * \param in  Pointer to object to be updated
 * \param err ObitErr for reporting errors.
 */
void ObitIOTableFITSUtilWriteTableList(ObitData *in, ObitErr *err)
{
  /* If (in->myStatus==OBIT_Modified) need to update tables */
  return;
} /* End of ObitIOTableFITSUtilWriteTableList */

/**
 * Get an associated table
 * \param in  Pointer to object with tables
 * \param err ObitErr for reporting errors.
 * \return Table or NULL on error
 */
ObitTable* ObitIOTableFITSUtilGetTable(ObitData *in,ObitIOAccess access, 
				       gchar *tabType, olong *tabVer, 
				       ObitErr *err)
{
  ObitTable *out;
  olong disk;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong version;
  gboolean gotIt;
  gchar *outName, FileName[256], tabName[51];
  gchar *routine = "ObitIOTableFITSUtilGetTable";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitDataIsA(in));
  g_assert(tabType!=NULL);
  g_assert(tabVer!=NULL);

  /* the Tablelist object must be present */
  if (in->tableList==NULL) {
      Obit_log_error(err, OBIT_Error, 
		     "my tableList member is NULL, open %s first", 
		     in->name);
      return NULL;
  }

  /* File info */
  if (!ObitInfoListGet(in->info, "Disk", &type, dim,  &disk, err)) 
    Obit_traceback_val (err, routine, in->name, NULL);
  if (!ObitInfoListGet(in->info, "FileName", &type, dim, FileName, err))
    Obit_traceback_val (err, routine, in->name, NULL);

  /* Do we already have this one? */
  version = *tabVer;
  gotIt = ObitTableListGet (in->tableList, tabType, &version, &out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

  /* Check if we're forcing a new table */
  if ((access==OBIT_IO_WriteOnly) && (*tabVer <= 0)) {
    version++;
    out = ObitTableUnref(out);
  }

  /* Set output table version */
  *tabVer = version;

  if (gotIt && (out!=NULL)) return out; /* that was easy */

  /* If it doesn't exist and request is read only - return NULL */
  if ((!gotIt) && (access==OBIT_IO_ReadOnly)) return NULL;

  /* Create one - make descriptive name */
  g_snprintf (tabName, 50, "%s table %d for ",tabType, *tabVer);
  outName =  g_strconcat (tabName, in->name, NULL);
  out = newObitTable (outName);
  g_free(outName);
  
  /* Setup info needed for access */
  ObitTableSetFITS(out, disk, FileName, tabType, version, 1, err);
 if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
 
 /* register it in the TableList */
 ObitTableListPut (in->tableList, tabType, &version, out, err);
 if (err->error)   Obit_traceback_val (err, routine, in->name, NULL);

 return out;
} /* end ObitIOTableFITSUtilGetTable */

