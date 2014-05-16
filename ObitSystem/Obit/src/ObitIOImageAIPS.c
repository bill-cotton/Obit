/* $Id$   */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include <stdio.h>
#include "Obit.h"
#include "ObitAIPS.h"
#include "ObitIOImageAIPS.h"
#include "ObitAIPSCat.h"
#include "ObitImageSel.h"
#include "ObitTableList.h"
#include "ObitMem.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOImageAIPS.c
 * ObitIOImageAIPS class function definitions.
 */

/*------------------- file globals ------------------------------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOImageAIPS";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOImageAIPSClassInfo myClassInfo = {FALSE}; 

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOImageAIPSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOImageAIPSClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitIOImageAIPSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class on the first call.
 * \param name An optional name for the object.
 * \param info if non-NULL it is used to initialize the new object.
 * \param err  ObitErr for error messages.
 * \return the new object.
 */
ObitIOImageAIPS* newObitIOImageAIPS (gchar *name, ObitInfoList *info,
				     ObitErr *err)
{
  ObitIOImageAIPS* out;
  gint32 i, dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "newObitIOImageAIPS";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOImageAIPSClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitIOImageAIPS), "ObitIOImageAIPS");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOImageAIPSInit((gpointer)out);

  /* Get any info from info input */
  if (info!=NULL) {
    type = OBIT_oint; for (i=0; i<MAXINFOELEMDIM; i++) dim[i] = 1;
    if(!ObitInfoListGet(info, "Disk", &type, (gint32*)dim, 
			&out->disk, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "User", &type, (gint32*)dim, 
			&out->UserId, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "CNO", &type, (gint32*)dim, 
			&out->CNO, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
  }  /* end of initialize from info */
  
  return out;
} /* end newObitIOImageAIPS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOImageAIPSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOImageAIPSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIOImageAIPSGetClass */

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
gboolean ObitIOImageAIPSSame (ObitIO *in, ObitInfoList *in1, 
				ObitInfoList *in2, ObitErr *err)
{
  olong CNO1, UserId1, disk1, CNO2, UserId2, disk2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOImageAIPSSame";

  /* error checks */
  if (err->error) return same;

  /* get instructions from info */
  if(!ObitInfoListGet(in1, "Disk", &type, dim, &disk1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in1, "User", &type, dim, &UserId1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in1, "CNO", &type, dim, &CNO1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "Disk", &type, dim, &disk2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "User", &type, dim, &UserId2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "CNO", &type, dim, &CNO2, err))
    Obit_traceback_val (err, routine, in->name, same);

  /* Compare */
  same = (disk1==disk2) && (CNO1==CNO2) && (UserId1==UserId2);

  return same;
} /* end ObitIOImageAIPSSame */

/**
 * Rename underlying files.
 * New name information is given on the info member:
 * \li "newName" OBIT_string (12,1,1) New AIPS Name 
 *      absent or Blank = don't change
 * \li "newClass" OBIT_string (6,1,1) New AIPS Class
 *      absent or Blank = don't changeO
 * \li "newSeq" OBIT_int (1,1,1) New AIPS Sequence
 *      0 => unique value
 * \param in Pointer to object to be zapped.
 * \param info Associated ObitInfoList
 * \li "Disk" OBIT_int (1,1,1)           Disk number
 * \li "CNO" OBIT_int (1,1,1)            Catalog slot number
 * \li "newName" OBIT_string (12,1,1) New AIPS Name 
 *      absent or Blank = don't change
 * \li "newClass" OBIT_string (6,1,1) New AIPS Class
 *      absent or Blank = don't changeO
 * \li "newSeq" OBIT_int (1,1,1) New AIPS Sequence
 *      0 => unique value
 * \param err ObitErr for reporting errors.
 */
void ObitIOImageAIPSRename (ObitIO *in, ObitInfoList *info, 
			    ObitErr *err)
{
  gchar *routine = "ObitIOImageAIPSRename";

  /* Don't bother if NULL */
  if (!in) return;

  /* error checks */
  if (err->error) return;

  /* Rename */
  ObitAIPSRename (in, info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitIOImageAIPSRename */

/**
 * Delete underlying files.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOImageAIPSZap (ObitIOImageAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableList *tableList=NULL;
  ObitTable *table=NULL;
  ObitFile *myFile=NULL;
  gchar *tabType=NULL;
  olong i, tabVer;
  gchar *routine = "ObitIOImageAIPSZap";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

   /* Close if still open */
  if ((in->myStatus==OBIT_Modified) || (in->myStatus==OBIT_Active)) {
    retCode = ObitIOImageAIPSClose (in, err);
    if (err->error) /* add traceback on error */
      Obit_traceback_msg (err, routine, in->name);
  }

  /* Clear entry from catalog - this may fail if entry too busy */
  ObitAIPSDirRemoveEntry(in->disk, in->UserId, in->CNO, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Delete any tables on the TableList */
  tableList = (ObitTableList*)in->tableList;
  if (tableList) {
    for (i=tableList->number; i>=1; i--) {
      /* Get info */
      ObitTableListGetNumber (tableList, i, &tabType, &tabVer, 
			      &table, err);
      
      /* setup input table if not instantiated */
      if (table==NULL) {
	table = (ObitTable*)newObitIOImageAIPSTable (in, OBIT_IO_ReadOnly, 
						     tabType, &tabVer, err);
      }
      
      /* Remove table from list */
      ObitTableListRemove (tableList, tabType, tabVer);
      
      /* destroy the table */
      table = ObitTableZap (table, err);
      table = ObitTableUnref(table);
      if (tabType) g_free(tabType);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    } /* End loop deleting tables */
  } /* end if tableList exists */
  while (tableList) tableList = ObitTableUnref(tableList);  /* Get table list */
  in->tableList = NULL;

  /* Get MA file */
  myFile = newObitFile("Files to Zap");
  myFile->fileName = ObitAIPSFilename (OBIT_AIPS_Image, in->disk, in->CNO, 
				 in->UserId, NULL, 0, err);
  myFile = ObitFileZap(myFile, err);
  myFile = ObitFileUnref(myFile);
  /* May be called a scratch 'SC' file */
  myFile = newObitFile("Files to Zap");
  myFile->fileName = ObitAIPSFilename (OBIT_AIPS_Scratch, in->disk, in->CNO, 
				 in->UserId, NULL, 0, err);
  myFile = ObitFileZap(myFile, err);
  myFile = ObitFileUnref(myFile);

  /* Get CB file */
  myFile = newObitFile("Files to Zap");
  myFile->fileName = ObitAIPSFilename (OBIT_AIPS_Header, in->disk, in->CNO, 
				 in->UserId, NULL, 0, err);
  myFile = ObitFileZap(myFile, err);
  myFile = ObitFileUnref(myFile);

  /* Trace back if error detected */
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  return;
} /* end ObitIOImageAIPSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOImageAIPS* ObitIOImageAIPSCopy  (ObitIOImageAIPS *in, 
				       ObitIOImageAIPS *out, ObitErr *err)
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
    out = newObitIOImageAIPS(outName, NULL, err);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy this class */
  out->disk   = in->disk;
  out->UserId = in->UserId;
  out->CNO    = in->CNO;
  if (out->AIPSFileName!=NULL) g_free(out->AIPSFileName);
  out->AIPSFileName = g_strdup(in->AIPSFileName);

  return out;
} /* end ObitIOImageAIPSCopy */

/**
 * Initialize structures and open file.
 * The file etc. info should have been stored in the ObitInfoList.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageAIPSOpen (ObitIOImageAIPS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  gchar *routine = "ObitIOImageAIPSOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));
  g_assert (in->myDesc != NULL);
  g_assert (in->mySel != NULL);

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* set defaults */
  desc->IOsize = OBIT_IO_byRow;

  /* get instructions from info */
  ObitInfoListGetTest(info, "IOBy", &type, (gint32*)dim, 
		      &desc->IOsize);

  if(!ObitInfoListGet(info, "Disk", &type, (gint32*)dim, 
		      &in->disk, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  if(!ObitInfoListGet(info, "User", &type, (gint32*)dim, 
		      &in->UserId, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  if(!ObitInfoListGet(info, "CNO", &type, (gint32*)dim, 
		      &in->CNO, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* form file name for image file */
  if (in->AIPSFileName) g_free(in->AIPSFileName); /* free old */
  in->AIPSFileName = 
    ObitAIPSFilename (OBIT_AIPS_Image, in->disk, in->CNO, 
		      in->UserId, NULL, 0, err);

  /* If in->myFile still connected unreference */
  if (in->myFile) in->myFile = ObitFileUnref (in->myFile);
  in->myFile = newObitFile(in->name);  /* new one */

  /* open file */
  retCode = OBIT_IO_OpenErr; /* in case something goes wrong */
  if (ObitFileOpen (in->myFile, in->AIPSFileName, access,  OBIT_IO_Binary,
		     0L, err) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  /* If it was just created, write header file */
  if (!in->myFile->exist) {
    if (ObitIOImageAIPSWriteDescriptor(in, err)|| (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* save some information */
  in->access = access;
  in->myStatus = OBIT_Active;
  desc->areBlanks = FALSE;
  /* ???  desc->maxval = -1.0e20;*/
  /* ???  desc->minval =  1.0e20;*/
    
  /* initialize location in image */
  desc->row    = 0;
  desc->plane  = 0;
  desc->plane4 = 0;
  desc->plane5 = 0;
  desc->plane6 = 0;
  desc->plane7 = 0;
  in->filePos  = 0;
  
  return OBIT_IO_OK;
} /* end ObitIOImageAIPSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOImageAIPSClose (ObitIOImageAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOImageAIPSClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* don't bother if it's not open */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) 
    return OBIT_IO_OK;

  if (ObitFileClose (in->myFile, err) || (err->error)) 
    /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  in->myStatus = OBIT_Inactive;
  return OBIT_IO_OK;
} /* end ObitIOImageAIPSClose */

/**
 * initialize I/O
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageAIPSSet (ObitIOImageAIPS *in, ObitInfoList *info, 
			       ObitErr *err)
{
  /* Reset values in descriptor */
  ((ObitImageDesc*)in->myDesc)->plane  = 0;
  ((ObitImageDesc*)in->myDesc)->plane4 = 0;
  ((ObitImageDesc*)in->myDesc)->plane5 = 0;
  ((ObitImageDesc*)in->myDesc)->plane6 = 0;
  ((ObitImageDesc*)in->myDesc)->plane7 = 0;
  ((ObitImageDesc*)in->myDesc)->row    = 0;

  return OBIT_IO_OK;
} /* end ObitIOImageAIPSSet */

/**
 * Read image data from disk.
 * Reads row in->myDesc->row + 1; plane in->myDesc->plane + 1
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageAIPSRead (ObitIOImageAIPS *in, ofloat *data, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  gsize size;
  olong offset, len=0, lRow, iRow, nRows=0;
  olong row, plane, plane4, plane5, plane6, plane7;
  ObitFilePos wantPos;
  olong  i, ipos[IM_MAXDIM];
  gchar *routine = "ObitIOImageAIPSRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  /* Flush if modified */
  if (in->myStatus==OBIT_Modified) {
    retCode =  ObitIOImageAIPSFlush(in, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) { /* add traceback on error */
  	Obit_log_error(err, OBIT_Error, 
		       "%s: Read flushing buffer in %s", routine, in->name);
    }
  }

  desc = in->myDesc; /* Image descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* where are we in image */
  row    = MAX (desc->row,    sel->blc[1]-1);
  plane  = MAX (desc->plane,  sel->blc[2]-1);
  plane4 = MAX (desc->plane4, sel->blc[3]);
  plane5 = MAX (desc->plane5, sel->blc[4]);
  plane6 = MAX (desc->plane6, sel->blc[5]);
  plane7 = MAX (desc->plane7, sel->blc[6]);
 
  /* set current request by desc->IOsize */
  if (desc->IOsize==OBIT_IO_byRow) {
    row++; /* increment row */
    if (row>sel->trc[1]) { /* next plane */
      row = sel->blc[1];
      plane++;
    }
    nRows = 1;
    len = sel->trc[0] - sel->blc[0]+1;   /* size of a transfer (row) */

  } else if (desc->IOsize==OBIT_IO_byPlane) {
    
    row = sel->blc[1]; /* set row */
    plane++; /* increment plane */
    /* must read a row at a time */
    nRows = sel->trc[1] - sel->blc[1] + 1;
    len = sel->trc[0] - sel->blc[0] + 1;   /* size of a transfer (row) */
  }
  if (plane>sel->trc[2]) { /* next plane4 */
    plane4++;
    plane = sel->blc[2];
  }
  if (plane4>sel->trc[3]) { /* next plane5 */
    plane5++;
    plane4 = sel->blc[3];
  }
  if (plane5>sel->trc[4]) { /* next plane6 */
    plane6++;
    plane5 = sel->blc[4];
  }
  if (plane6>sel->trc[5]) { /* next plane7 */
    plane7++;
    plane6 = sel->blc[5];
  }
  desc->row    = row;
  desc->plane  = plane;
  desc->plane4 = plane4;
  desc->plane5 = plane5;
  desc->plane6 = plane6;
  desc->plane7 = plane7;

  /* check if done - starting on the plane past the highest. */
  if ((plane>sel->trc[2]) || (plane4>sel->trc[3]) || (plane5>sel->trc[4]) ||
      (plane5>sel->trc[5]) || (plane7>sel->trc[6])) {
    /* ObitIOImageAIPSClose (in, err); Close */
    return OBIT_IO_EOF;
  }

  /* position of first pixel to access */
  for (i=0; i<IM_MAXDIM; i++) ipos[i] = sel->blc[i];
  ipos[1] = row;
  ipos[2] = MAX(1,plane);
  ipos[3] = MAX(1,plane4);
  ipos[4] = MAX(1,plane5);
  ipos[5] = MAX(1,plane6);
  ipos[6] = MAX(1,plane7);

  size = len * sizeof(ofloat);           /* transfer size in bytes */
  lRow = desc->inaxes[0]*sizeof(ofloat); /* length of transfer in bytes */

  offset = 0; /* offset in buffer */

  /* read file one row/plane at a time - loop for windowed planes */
  retCode = OBIT_IO_ReadErr; /* in case something goes wrong */
  for (iRow=0; iRow<nRows; iRow++) {

    /* get file position offset (AIPS images have funny rules) */
    wantPos = ObitAIPSImageFileOffset(desc->naxis, (olong*)desc->inaxes, ipos);

    /* Read */
    retCode = ObitFileRead (in->myFile, wantPos, size, 
			    (gchar*)&data[offset], err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) { /* add traceback on error */
      if (retCode==OBIT_IO_EOF)
	Obit_log_error(err, OBIT_Error, 
		       "%s: Hit EOF in %s", routine, in->name);
      else
	Obit_log_error(err, OBIT_Error, 
		       "%s: Read error in %s", routine, in->name);
      return retCode;
    } /* end error trap */
    in->filePos = in->myFile->filePos; /* remember current file position */

    ipos[1]++;            /* next row */
    offset  += len;       /* offset in data buffer */
  }  /* end loop reading rows */
  
  return  OBIT_IO_OK;
} /* end ObitIOImageAIPSRead */

/**
 * Write information to disk.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageAIPSWrite (ObitIOImageAIPS *in, ofloat *data, 
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  gsize size;
  ofloat val, fblank = ObitMagicF();
  olong i, offset, len=0, lRow, iRow, nRows=0;
  olong row, plane, plane4, plane5, plane6, plane7;
  ObitFilePos wantPos;
  olong  ipos[IM_MAXDIM], inaxes[IM_MAXDIM];
  gchar *routine = "ObitIOImageAIPSWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitFileIsA(in->myFile));
  g_assert (data != NULL);

  desc = in->myDesc; /* Image descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* where are we in image */
  row    = MAX (desc->row,    sel->blc[1]-1);
  plane  = MAX (desc->plane,  sel->blc[2]-1);
  plane4 = MAX (desc->plane4, sel->blc[3]);
  plane5 = MAX (desc->plane5, sel->blc[4]);
  plane6 = MAX (desc->plane6, sel->blc[5]);
  plane7 = MAX (desc->plane7, sel->blc[6]);

  /* set current request by desc->IOsize */
  if (desc->IOsize==OBIT_IO_byRow) {
    
    plane = MAX (sel->blc[2], plane);
    row++;   /* increment row */
    if (row>sel->trc[1]) { /* next plane */
      row = sel->blc[1];
      plane++;
    }
    nRows = 1;
    len = sel->trc[0] - sel->blc[0]+1;   /* size of a transfer (row) */

  } else if (desc->IOsize==OBIT_IO_byPlane) {
    
    row = sel->blc[1]; /* set row */
    plane++; /* increment plane */
     /* how many rows to write  */
    nRows = sel->trc[1] - sel->blc[1] + 1;
    len = sel->trc[0] - sel->blc[0]+1;   /* size of a transfer (row) */
  }

  if (plane>sel->trc[2]) { /* next plane4 */
    plane4++;
    plane = sel->blc[2];
  }
  if (plane4>sel->trc[3]) { /* next plane5 */
    plane5++;
    plane4 = sel->blc[3];
  }
  if (plane5>sel->trc[4]) { /* next plane6 */
    plane6++;
    plane5 = sel->blc[4];
  }
  if (plane6>sel->trc[5]) { /* next plane7 */
    plane7++;
    plane6 = sel->blc[5];
  }
  desc->row    = row;
  desc->plane  = plane;
  desc->plane4 = plane4;
  desc->plane5 = plane5;
  desc->plane6 = plane6;
  desc->plane7 = plane7;

  /* check if done - starting on the plane past the highest. */
  if (plane > sel->trc[2]) {
    /* ObitIOImageAIPSClose (in, err); Close */
    return OBIT_IO_EOF;
  }

  for (i=0; i<IM_MAXDIM; i++) inaxes[i] = desc->inaxes[i];
  /* position of first pixel to access */
  for (i=0; i<IM_MAXDIM; i++) ipos[i] = sel->blc[i];
  ipos[1]  = row;
  ipos[2]  = MAX(1,plane);
  ipos[3]  = MAX(1,plane4);
  ipos[4]  = MAX(1,plane5);
  ipos[5]  = MAX(1,plane6);
  ipos[6]  = MAX(1,plane7);

  size = len * sizeof(ofloat);           /* transfer size in bytes */
  lRow = desc->inaxes[0]*sizeof(ofloat); /* length of a row in bytes */

  offset = 0; /* offset in output buffer */

  /* get initial file position offset (AIPS images have funny rules) */
  wantPos = ObitAIPSImageFileOffset(desc->naxis, inaxes, ipos);

  /* write file one row/plane at a time - loop for windowed planes */
  for (iRow=0; iRow<nRows; iRow++) {

    /* keep track on max/min/blanking */
    for (i=0; i<len; i++) {
      val = data[offset+i];
      if (val==fblank) {
	desc->areBlanks = TRUE;
      } else { /* OK */
	desc->maxval = MAX (desc->maxval, val);
	desc->minval = MIN (desc->minval, val);
      }
    }
    /* Write */
    retCode = ObitFileWrite (in->myFile, wantPos, size, 
			    (gchar*)&data[offset], err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */

    ipos[1]++;            /* next row */
    offset  += len;       /* offset on data buffer */
    wantPos += lRow;      /* where we want to start next read */  

    /* next file position offset (AIPS images have funny rules) */
    wantPos = ObitAIPSImageFileOffset(desc->naxis, (olong*)desc->inaxes, ipos);

    /* Do we need to pad out a block just finished? */
    if ((in->myFile->access == OBIT_IO_WriteOnly) && 
	(iRow>0) && (wantPos>in->myFile->filePos)) {
      retCode = ObitFilePad (in->myFile, wantPos, (gchar*)data, size, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      in->filePos = in->myFile->filePos; /* remember current file position */
    }

  }  /* end loop writing rows */

  /* AIPS' wonky image file format requires that the last block be written */
  /* if a plane was just finished - fill to the beginning of the next plane */
  /* where is start of the next plane */
  if (ipos[1] >= sel->trc[1]) {
    ipos[0] = 1;
    ipos[1] = 1;
    ipos[2]++;
    wantPos = ObitAIPSImageFileOffset(desc->naxis, (olong*)desc->inaxes, ipos);
  
    /* Do we need to pad out a block just finished? */
    if (((in->myFile->access == OBIT_IO_WriteOnly) || 
	 (in->myFile->access == OBIT_IO_ReadWrite)) &&
	(wantPos > in->myFile->filePos)) {
      retCode = ObitFilePad (in->myFile, wantPos, (gchar*)data, size, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      in->filePos = in->myFile->filePos; /* remember current file position */
    }
  } /* end of padding section */
    
 
  in->myStatus = OBIT_Modified; /* file has been modified */

  return  OBIT_IO_OK;
} /* end ObitIOImageAIPSWrite */

/**
 * Read image Descriptor data from disk.
 * \param in Pointer to object with ObitImageDesc to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitIOImageAIPSReadDescriptor (ObitIOImageAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  gchar *HeaderFile, keyName[9], blob[256];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  AIPSint buffer[260];
  olong i, j, ip, ndo, nkey, ikey, nrec;
  gsize size;
  ObitFilePos wantPos;
  ObitFile *myFile=NULL;
  ObitInfoType keyType;
  gchar *routine = "ObitIOImageAIPSReadDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* debug */
  for (i=0; i<256; i++) buffer[i] = 0;

  desc = in->myDesc; /* Image descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* Set file name */
  HeaderFile = ObitAIPSFilename (OBIT_AIPS_Header, in->disk, in->CNO, 
				 in->UserId, NULL, 0, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* swallow file - open */
  myFile = newObitFile(in->name);
  size = 260 * sizeof(AIPSint);
  if (ObitFileOpen (myFile, HeaderFile, OBIT_IO_ReadOnly, 
		     OBIT_IO_Binary, size, err) ||
      (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  g_free(HeaderFile); HeaderFile = NULL;  /* cleanup */

  /* read */
  wantPos = 0;
  retCode = ObitFileRead (myFile, wantPos, size, (gchar*)buffer, err);
  if (retCode==OBIT_IO_EOF) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Hit EOF in AIPS Header for %s", routine, in->name);
    return retCode;
  }
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Convert to internal structure */
  ObitAIPSCatImageGetDesc (desc, (gchar*)buffer, err);
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Get table Info to Table List */
  ObitAIPSCatGetTable ((ObitTableList*)in->tableList, (gchar*)buffer,
		       in->UserId, in->disk, in->CNO, err);
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /*++++++++++++ Read keyword/values ++++++++++++++++++*/
  /* the number of keywords is in the beginning of the 2nd 256 AIPSint
     block and has been already read into buffer */
  nrec = buffer[256]-1; /* How many blocks? */
  nkey = buffer[257];   /* How many keywords? */

  /* delete old InfoList and restart */
  ((ObitImageDesc*)in->myDesc)->info = ObitInfoListUnref (((ObitImageDesc*)in->myDesc)->info);
  ((ObitImageDesc*)in->myDesc)->info = (gpointer)newObitInfoList ();
  desc = in->myDesc; /* Table descriptor pointer */
    
  if (nkey>0) {
    size    = 256 * sizeof(AIPSint);
    wantPos = 256 * sizeof(AIPSint); /* Initial File location */

    /* loop reading and parsing */
    ikey = 0;
    for (i=0; i<nrec; i++) {
      /* read block */
      retCode = ObitFileRead (myFile, wantPos, size, (gchar*)buffer, err);
      if (retCode==OBIT_IO_EOF) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: Hit EOF in AIPS Header for %s", routine, in->name);
	return retCode;
      }
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      wantPos = -1L; /* now sequential */

      /* Parse - how many? */
      ndo = 256/5;
      ndo = MIN (ndo, (nkey+1-ikey)); /* Not more than number of keywords */
      ip = 0;
      for (j=0; j<ndo; j++) {
	/* ignore the first in the first block - this is something eles */
	if (ikey==0) {ip +=6; ikey++; continue;}
	g_memmove (&keyName[0], (gchar*)&buffer[ip],   4);
	g_memmove (&keyName[4], (gchar*)&buffer[ip+1], 4); keyName[8] = 0;
	/* Save 8 bytes of data */
	g_memmove (blob, (gchar*)&buffer[ip+2], 8); blob[8] = 0;
	/* type as ObitInfoType */
	keyType = OBIT_oint;
	if (buffer[ip+4]==1) keyType = OBIT_double;
	else if (buffer[ip+4]==2) keyType = OBIT_float;
	else if (buffer[ip+4]==3) keyType = OBIT_string;
	else if (buffer[ip+4]==4) keyType = OBIT_oint;
	else if (buffer[ip+4]==5) keyType = OBIT_bool;
	
	/* Save on ObitInfoList */
	dim[0] = 1;
	if (keyType == OBIT_string) dim[0] = 8;
	ObitInfoListPut(desc->info, keyName, keyType, dim, 
			(gconstpointer)blob, err);
	if (err->error)  /* add trace and return on error */
	  Obit_traceback_val (err, routine,  in->name, retCode);
	ip +=5; /* next unit in record */
      }
    } /* end loop reading keyword/values */
  } /* end of read keywords section */

  /* close file */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* delete */
  myFile = ObitFileUnref(myFile);

  return OBIT_IO_OK;
} /* end ObitIOImageAIPSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in Pointer to object with ObitImageDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOImageAIPSWriteDescriptor (ObitIOImageAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitImageDesc* desc;
  ObitImageSel* sel;
  gchar *HeaderFile, keyName[FLEN_KEYWORD+1], *keyNameP, blob[256], *ctemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, ii, j, k, ip, ndo, nkey, ikey, nrec;
  oint oitemp[4]={0,0,0,0};
  gboolean doFill, doInit;
  gsize size;
  ObitFilePos wantPos;
  AIPSint buffer[256];
  ObitFile *myFile=NULL;
  ObitAIPSDirCatEntry *dirEntry = NULL;
  ObitInfoType keyType;
  gchar *routine = "ObitIOImageAIPSWriteDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* enforce descriptor defaults */
  desc = in->myDesc; /* Image descriptor pointer */
  sel  = in->mySel;  /* selector pointer */
  ObitImageSelDefault(desc, sel);

  for (i=0; i<256; i++) buffer[i] = 0;

  /* Set file name */
  HeaderFile = ObitAIPSFilename (OBIT_AIPS_Header, in->disk, in->CNO, 
				 in->UserId, NULL, 0, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* open Header file */
  myFile = newObitFile(in->name);
  size = 256 * sizeof(AIPSint);
  if (ObitFileOpen (myFile, HeaderFile, OBIT_IO_ReadWrite, 
		     OBIT_IO_Binary, size, err) ||
      (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  g_free(HeaderFile); HeaderFile = NULL;  /* cleanup */

  /* if it exists read old and update */
  if (myFile->exist) {
    /* read */
    wantPos = 0;
    retCode = ObitFileRead (myFile, wantPos, size, (gchar*)buffer, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  } /* end of read old file section */
  
  /* convert descriptor to header */
  retCode = OBIT_IO_ReadErr;
  /* Get catalog descriptor */
  dirEntry = ObitAIPSDirGetEntry(in->disk, in->UserId, in->CNO, err);
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* do conversion */
  doInit = (!myFile->exist) || (in->access==OBIT_IO_WriteOnly);
  ObitAIPSCatImageSetDesc (in->myDesc, (gchar*)buffer, doInit, 
			   dirEntry, err);
  if (dirEntry) g_free(dirEntry); /* free catalog directory entry */
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Add associated table information */
  ObitAIPSCatSetTable ((ObitTableList*)in->tableList, (gchar*)buffer, err);
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Now (re)write it */
  wantPos = 0; /* File location */
  /* write it  */
  retCode = ObitFileWrite (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /*++++++++++++ Write keyword/values ++++++++++++++++++*/
  /* the number of keywords is in the beginning of the 2nd 256 AIPSint
     block and has been already read into buffer */

  nkey = 0; /* How many keywords+description? */
  if (desc->info!=NULL)  nkey = desc->info->number+1;
    
  nrec = 1 + (nkey/(256/5));   /* How many blocks? */
  /* Must write at least the first record */
    
  wantPos = 256 * sizeof(AIPSint); /* File location */
  /* loop writing */
  ikey = 0;
  for (i=0; i<nrec; i++) {
    for (ii=0; ii<256; ii++) buffer[ii] = 0; /* init buffer */
    /* Parse - how many? */
    ndo = 256/5;
    ndo = MIN (ndo, (nkey-ikey)); /* Not more than number of keywords */
    ip = 0;
    for (j=0; j<ndo; j++) {
      /* First entry of first block is something else - write it */
      if (ikey==0) {
	buffer[0] = (AIPSint)(nrec+1);
	buffer[1] = (AIPSint)(nkey-1);
	ikey++;
	ip +=6; /* next unit in buffer */
      } else { /* subsequent */
	
	/* Read from ObitInfoList */
	if (desc->info!=NULL) { /* Only if infoList exists */
	  ikey++;
	  ObitInfoListGetNumber(desc->info, ikey-1, &keyNameP, &keyType, 
				dim, blob, err);
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
	  /* Copy to buffer */
	  g_memmove ((gchar*)&buffer[ip],  &keyName[0], 4);
	  g_memmove ((gchar*)&buffer[ip+1],&keyName[4], 4); 
	  /* Save 8 bytes of data */
	  g_memmove ((gchar*)&buffer[ip+2], blob, 8); blob[8] = 0;
	  /* Convert type to AIPSish */
	  buffer[ip+4] = 4; /* default int */
	  if (keyType==OBIT_double)      buffer[ip+4] = 1;
	  else if (keyType==OBIT_float)  buffer[ip+4] = 2;
	  else if (keyType==OBIT_string) {
	    /* Replace any non printing characters with blanks */
	    ctemp = (gchar*)&buffer[ip+2];
	    for (k=0; k<8; k++) if (!g_ascii_isprint(ctemp[k])) ctemp[k]=' ';
	    buffer[ip+4] = 3;
	  }
	  else if (keyType==OBIT_oint)   buffer[ip+4] = 4;
	  else if (keyType==OBIT_long) { /* May have to convert long->oint */
	    buffer[ip+4] = 4;
	    oitemp[0] = (oint)*(olong*)blob;
	    g_memmove ((gchar*)&buffer[ip+2], oitemp, 8); blob[8] = 0;
	  }
	  else if (keyType==OBIT_bool)   buffer[ip+4] = 5;
	  ip +=5; /* next unit in buffer */
	} /* end if exists */
      }
    } /* end loop filling block */
    
    /* write block */
    retCode = ObitFileWrite (myFile, wantPos, size, (gchar*)buffer, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    wantPos = -1L; /* now sequential */
  }  /* end write keywords section */

  /* flush/close file */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* delete */
  myFile = ObitFileUnref(myFile);

  return OBIT_IO_OK;
} /* end ObitIOImageAIPSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOImageAIPSFlush (ObitIOImageAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  retCode = ObitFileFlush(in->myFile, err);

  if (in->myStatus==OBIT_Modified) in->myStatus = OBIT_Active;
  return retCode;
} /* end ObitIOImageAIPSFlush */

/**
 * Create buffer appropriate for I/O request
 * Not actually used for Images.
 * \param data (output) pointer to data array
 * \param size (output) size of data array in floats.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 */
void 
ObitIOImageAIPSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOImageAIPS *in, ObitInfoList *info, 
			     ObitErr *err)
{
  /* just return - never called */
} /* end ObitIOImageAIPSCreateBuffer */

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
newObitIOImageAIPSTable (ObitIOImageAIPS *in, ObitIOAccess access, 
			 gchar *tabType, olong *tabVer, ObitErr *err)
{
  ObitTable *out;
  olong version;
  gboolean gotIt;
  gchar ttype[3], *outName, tabName[51];
  gchar *routine = "newObitIOImageAIPSTable";

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
  ttype[0] = tabType[5]; ttype[1] = tabType[6]; ttype[2] = 0;
  ObitTableSetAIPS(out, in->disk, in->CNO, ttype, version, 
		  in->UserId, 25, err);
  if (err->error)
    Obit_traceback_val (err, routine, in->name, NULL);
 
 
 /* register it in the TableList */
  ObitTableListPut ((ObitTableList*)in->tableList, tabType, &version, 
		    out, err);
  if (err->error)
    Obit_traceback_val (err, routine, in->name, NULL);
  
  /* Force Write to disk
  ObitIOImageAIPSWriteDescriptor (in, err); */
  if (err->error)
    Obit_traceback_val (err, routine, in->name, NULL);
  
  return (Obit*)out;
} /* end newObitIOImageAIPSTable */

/**
 * Update any disk resident structures about the current tables.
 * \param in   Pointer to object to be updated.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOImageAIPSUpdateTables (ObitIOImageAIPS *in, 
					ObitInfoList *info, 
					ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOImageAIPSUpdateTables";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  retCode = ObitIOImageAIPSWriteDescriptor(in, err);
  if ((retCode!= OBIT_IO_OK) || err->error)
    Obit_traceback_val (err, routine, in->name, retCode);
    
  return retCode;
} /* end ObitIOImageAIPSUpdateTables */

/**
 * Get underlying file information in entries to an ObitInfoList
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
 * For all File types types:
 * \li xxxDataType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxFileType OBIT_oint File type as ObitIOType, OBIT_IO_FITS, OBIT_IO_AIPS
 *    
 * For xxxDataType = "MA"
 * \li xxxBLC   OBIT_oint[7] (Images only) 1-rel bottom-left corner pixel
 * \li xxxTRC   OBIT_oint[7] (Images Only) 1-rel top-right corner pixel
 * \param err     ObitErr for reporting errors.
 */
void ObitIOImageAIPSGetFileInfo (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
				 ObitInfoList *outList, ObitErr *err)
{
  ObitInfoType type;
  ObitAIPSDirCatEntry *entry=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *FileType="Image", *DataType="AIPS", *dirname;
  olong disk, user, cno, seq, i;
  gpointer     listPnt;
  gchar *parm[] = {"BLC", "TRC", NULL};
  gchar *routine = "ObitIOImageAIPSGetFileInfo";

  if (err->error) return;

  /* Set basic information */
  if (prefix) keyword =  g_strconcat (prefix, "FileType", NULL);
  else keyword =  g_strdup ("FileType");
  dim[0] = strlen(FileType);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, FileType);
  g_free(keyword);
  
   /* AIPS */
  if (prefix) keyword =  g_strconcat (prefix, "DataType", NULL);
  else keyword =  g_strdup ("DataType");
  dim[0] = strlen(DataType);
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, DataType);
  g_free(keyword);
  
  /* Disk number */
  if (!ObitInfoListGet(myInfo, "Disk", &type, dim, &disk, err)) 
      Obit_traceback_msg (err, routine, in->name);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (prefix) keyword =  g_strconcat (prefix, "Disk", NULL);
  else keyword =  g_strdup ("Disk");
  ObitInfoListAlwaysPut (outList, keyword, type, dim, &disk);
  g_free(keyword);
 
  /* AIPS user number */
  if (!ObitInfoListGet(myInfo, "User", &type, dim, &user, err)) 
      Obit_traceback_msg (err, routine, in->name);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (prefix) keyword =  g_strconcat (prefix, "User", NULL);
  else keyword =  g_strdup ("User");
  ObitInfoListAlwaysPut (outList, keyword, type, dim, &user);
  ObitInfoListAlwaysPut (outList, "AIPSuser", type, dim, &user);
  g_free(keyword);
 
  /* Catalog slot number */
  if (!ObitInfoListGet(myInfo, "CNO", &type, dim, &cno, err)) 
      Obit_traceback_msg (err, routine, in->name);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (prefix) keyword =  g_strconcat (prefix, "CNO", NULL);
  else keyword =  g_strdup ("CNO");
  ObitInfoListAlwaysPut (outList, keyword, type, dim, &disk);
  g_free(keyword);

  /* Look up info */
  entry = ObitAIPSDirGetEntry (disk, user, cno, err);

  /*  AIPS name */
  if (prefix) keyword =  g_strconcat (prefix, "Name", NULL);
  else keyword =  g_strdup ("Name");
  dim[0] = 12;
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, entry->name);
  g_free(keyword);

  /*  AIPS class */
  if (prefix) keyword =  g_strconcat (prefix, "Class", NULL);
  else keyword =  g_strdup ("Class");
  dim[0] = 6;
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, entry->class);
  g_free(keyword);

  /*  AIPS sequence */
  seq = (olong)entry->seq;
  if (prefix) keyword =  g_strconcat (prefix, "Seq", NULL);
  else keyword =  g_strdup ("Seq");
  dim[0] = 1;
  ObitInfoListAlwaysPut (outList, keyword, OBIT_long, dim, &seq);
  g_free(keyword);
  if (entry) g_free(entry);
 
  /* Disk directory  */
  dirname = ObitAIPSDirname (disk, err);
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
} /* end ObitIOImageAIPSGetFileInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOImageAIPSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOImageAIPSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOImageAIPSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOImageAIPSClassInfoDefFn (gpointer inClass)
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
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOImageAIPSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOImageAIPSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOImageAIPSGetClass;
  theClass->newObitIO     = (newObitIOFP)newObitIOImageAIPS;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOImageAIPSSame;
  theClass->ObitIORename  = (ObitIORenameFP)ObitIOImageAIPSRename;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOImageAIPSZap;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOImageAIPSCopy;
  theClass->ObitClear     = (ObitClearFP)ObitIOImageAIPSClear;
  theClass->ObitInit      = (ObitInitFP)ObitIOImageAIPSInit;
  theClass->ObitIOOpen    = (ObitIOOpenFP)ObitIOImageAIPSOpen;
  theClass->ObitIOClose   = (ObitIOCloseFP)ObitIOImageAIPSClose;
  theClass->ObitIOSet     = 
    (ObitIOSetFP)ObitIOImageAIPSSet;
  theClass->ObitIORead    = 
    (ObitIOReadFP)ObitIOImageAIPSRead;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOImageAIPSRead;
  theClass->ObitIOWrite   = 
    (ObitIOWriteFP)ObitIOImageAIPSWrite;
  theClass->ObitIOFlush   = 
    (ObitIOFlushFP)ObitIOImageAIPSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOImageAIPSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOImageAIPSWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOImageAIPSCreateBuffer;
  theClass->ObitIOFreeBuffer   = 
    (ObitIOFreeBufferFP)ObitIOFreeBuffer;
  theClass->newObitIOTable   = 
    (newObitIOTableFP)newObitIOImageAIPSTable;
  theClass->ObitIOUpdateTables   = 
    (ObitIOUpdateTablesFP)ObitIOImageAIPSUpdateTables;
  theClass->ObitIOGetFileInfo   =
    (ObitIOGetFileInfoFP)ObitIOImageAIPSGetFileInfo;
} /* end ObitIOImageAIPSClassDefFn */

/*--------------- Private functions --------------------------*/

/**
 * Creates empty member objects.
 * for each parent class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOImageAIPSInit  (gpointer inn)
{
  const ObitClassInfo *ParentClass;
  ObitIOImageAIPS *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && (ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->AIPSFileName = NULL;
  in->disk         = 0;
  in->UserId       = 0;
  in->myFile       = NULL;
  in->filePos      = 0;

} /* end ObitIOImageAIPSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOImageAIPSClear (gpointer inn)
{
  ObitIOImageAIPS *in = inn;
  const ObitClassInfo *ParentClass;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->myStatus==OBIT_Active) ||(in->myStatus==OBIT_Modified)) {
    err = newObitErr();
    ObitIOImageAIPSClose (in, err); 
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  if (in->AIPSFileName) g_free(in->AIPSFileName); 
  in->AIPSFileName = NULL;
  if (in->myFile) ObitFileUnref (in->myFile);

 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && (ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOImageAIPSClear */

