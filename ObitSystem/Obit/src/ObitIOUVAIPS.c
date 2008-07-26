/* $Id$      */
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
#include <stdio.h>
#include "ObitIOUVAIPS.h"
#include "ObitAIPSCat.h"
#include "ObitAIPS.h"
#include "ObitTableList.h"
#include "ObitFile.h"
#include "ObitMem.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOUVAIPS.c
 * ObitIOUVAIPS class function definitions.
 */

/*------------------- file globals ------------------------------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOUVAIPS";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOUVAIPSClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOUVAIPSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOUVAIPSClear (gpointer in);

/** Private: Determine next visibility to read */
static gboolean ObitIOUVAIPSNext (ObitIOUVAIPS *in, ObitErr *err);

/** Private: Compress visibilities. */
static void 
ObitIOUVAIPSCompress (gint ncorr, const ofloat *visin, ofloat *wtscl, 
		      ofloat *visout);

/** Private: Uncompress visibilities. */
static void 
ObitIOUVAIPSUncompress (gint ncorr, const ofloat *visin, 
			const ofloat *wtscl, ofloat *visout);

/** Private: Set Class function pointers. */
static void ObitIOUVAIPSClassInfoDefFn (gpointer inClass);

/** Private: Check file validity */
static void check (ObitIOUVAIPS *in, ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class on the first call.
 * \param name An optional name for the object.
 * \param info if non-NULL it is used to initialize the new object.
 * \param err  ObitErr for error messages.
 * \return the new object.
 */
ObitIOUVAIPS* newObitIOUVAIPS (gchar *name, ObitInfoList *info,
			       ObitErr *err)
{
  ObitIOUVAIPS* out;
  gint32 i, dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "newObitIOUVAIPS";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOUVAIPSClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitIOUVAIPS));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOUVAIPSInit((gpointer)out);

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
} /* end newObitIOUVAIPS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOUVAIPSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOUVAIPSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIOUVAIPSGetClass */

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
gboolean ObitIOUVAIPSSame (ObitIO *in, ObitInfoList *in1, 
			   ObitInfoList *in2, ObitErr *err)
{
  olong CNO1, UserId1, disk1, CNO2, UserId2, disk2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOUVAIPSSame";

  /* error checks */
  if (err->error) return same;
  g_assert (ObitIOUVAIPSIsA(in));

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
} /* end ObitIOUVAIPSSame */

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
void ObitIOUVAIPSRename (ObitIO *in, ObitInfoList *info, 
			    ObitErr *err)
{
  gchar *routine = "ObitIOUVAIPSRename";

  /* Don't bother if NULL */
  if (!in) return;

  /* error checks */
  if (err->error) return;

  /* Rename */
  ObitAIPSRename (in, info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitIOUVAIPSRename */

/**
 * Delete underlying files.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOUVAIPSZap (ObitIOUVAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableList *tableList=NULL;
  ObitTable *table=NULL;
  ObitFile *myFile=NULL;
  gchar *tabType=NULL;
  olong i, tabVer;
  gchar *routine = "ObitIOUVAIPSZap";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Close if still open */
  if ((in->myStatus==OBIT_Modified) || (in->myStatus==OBIT_Active)) {
    retCode = ObitIOUVAIPSClose (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Clear entry from catalog - this may fail if entry too busy */
  ObitAIPSDirRemoveEntry(in->disk, in->UserId, in->CNO, err);
  if (err->error)
    Obit_traceback_msg (err, routine, in->name);

  /* Delete any tables on the TableList */
  tableList = (ObitTableList*)in->tableList;
  for (i=tableList->number; i>=1; i--) {
    /* Get info */
    ObitTableListGetNumber (tableList, i, &tabType, &tabVer, 
			    &table, err);

    /* setup input table if not instantiated */
    if (table==NULL) {
      table = (ObitTable*)newObitIOUVAIPSTable (in, OBIT_IO_ReadOnly, 
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
  while (tableList) tableList = ObitTableUnref(tableList);  /* Get table list */
  in->tableList = NULL;

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Get UV file */
  myFile = newObitFile("Files to Zap");
  myFile->fileName = ObitAIPSFilename (OBIT_AIPS_UVdata, in->disk, in->CNO, 
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
} /* end ObitIOUVAIPSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOUVAIPS* ObitIOUVAIPSCopy  (ObitIOUVAIPS *in, 
				       ObitIOUVAIPS *out, ObitErr *err)
{
  const ObitIOClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  gchar *routine = "ObitIOUVAIPSCopy";

  /* error checks */
  g_assert (err!=NULL);
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* more complete validity tests */
  check (in, err);
  check (out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

   /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitIOUVAIPS(outName, NULL, err);
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
} /* end ObitIOUVAIPSCopy */

/**
 * Initialize structures and open file.
 * The file etc. info should have been stored in the ObitInfoList.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite, 
 *               OBIT_IO_ReadCal).
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVAIPSOpen (ObitIOUVAIPS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gchar *routine = "ObitIOUVAIPSOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (in->myDesc != NULL);
  g_assert (in->mySel  != NULL);

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  desc = in->myDesc; /* descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* get instructions from info */
  if(!ObitInfoListGet(info, "Disk", &type, (gint32*)dim, 
		      &in->disk, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  if(!ObitInfoListGet(info, "User", &type, (gint32*)dim, 
		      &in->UserId, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  if(!ObitInfoListGet(info, "CNO", &type, (gint32*)dim, 
		      &in->CNO, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* form file name for file */
  if (in->AIPSFileName) g_free(in->AIPSFileName); /* free old */
  in->AIPSFileName = 
    ObitAIPSFilename (OBIT_AIPS_UVdata, in->disk, in->CNO, 
		      in->UserId, NULL, 0, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* If in->myFile still connected unreference */
  if (in->myFile) ObitFileUnref (in->myFile);
  in->myFile = newObitFile(in->name);  /* new one */

  /* open file */
  retCode = OBIT_IO_OpenErr; /* in case something goes wrong */
  if (ObitFileOpen (in->myFile, in->AIPSFileName, access,  OBIT_IO_Binary,
		     0L, err) || (err->error)) 
    Obit_traceback_val (err, routine, in->name, retCode);

   /* If it was just created, write header file */
  if (!in->myFile->exist) {
    if (ObitIOUVAIPSWriteDescriptor(in, err)|| (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);
  }

 /* save some information */
  in->access = access;
  in->myStatus = OBIT_Active;
  if ((access == OBIT_IO_ReadWrite) || (access == OBIT_IO_ReadOnly) ||
      (access == OBIT_IO_ReadCal)) {
  }
  
  /* initialize location in data */
  desc->firstVis   = 0;
  desc->numVisBuff = 0;
  in->filePos = 0;
  
  /* Validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOUVAIPSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOUVAIPSClose (ObitIOUVAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOUVAIPSClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* don't bother if it's not open */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) 
    return OBIT_IO_OK;

  if (ObitFileClose (in->myFile, err) || (err->error)) 
    /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Delete any compression buffers */
  if (in->compBuff)  in->compBuff  = ObitMemFree (in->compBuff);  
  in->compBuffSize = 0;
  if (in->decompVis) in->decompVis = ObitMemFree (in->decompVis);  

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  in->myStatus = OBIT_Inactive;
  return OBIT_IO_OK;
} /* end ObitIOUVAIPSClose */

/**
 * initialize I/O
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVAIPSSet (ObitIOUVAIPS *in, ObitInfoList *info, 
			       ObitErr *err)
{
  gchar *routine = "ObitIOUVAIPSSet";

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, OBIT_IO_OK);

  /* reset vis pointers */
  ((ObitUVDesc*)(in->myDesc))->firstVis   = 0;
  ((ObitUVDesc*)(in->myDesc))->numVisBuff = 0;
  return OBIT_IO_OK;
} /* end ObitIOUVAIPSSet */

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
ObitIOCode ObitIOUVAIPSRead (ObitIOUVAIPS *in, ofloat *data, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gsize size;
  olong len, i, j, ip, op, need;
  ObitFilePos wantPos;
  gboolean done, compressed;
  ofloat *wtscl, *IOBuff = data;
  gchar *routine = "ObitIOUVAIPSRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  desc = in->myDesc; /* UV descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* what next ? */
  done = ObitIOUVAIPSNext (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* check if done - all visibilities read */
  if (done) {
    /* ObitIOUVAIPSClose (in, err);  Close */
    desc->numVisBuff = 0; /* no data in buffer*/
    return OBIT_IO_EOF;
  }

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;
  if (compressed) {
    /* buffer if necessary create */
    if (in->compBuff==NULL) {
      in->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
      in->compBuff = ObitMemAllocName (in->compBuffSize, "UVAIPS comp. buf");

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
		     "Decompression buffer ( %d) too small, need %d for %s", 
		     in->compBuffSize, need, in->name);
      return retCode;
    }
  } /* end of compressed data set up */

  len = desc->lrec; /* How big is a visibility */

  /* read block of sel->numVisRead visibilities at a time  */
  /* get file position offset */
  wantPos = (ObitFilePos)(desc->firstVis-1) * desc->lrec * sizeof(ofloat);
  /* transfer size in bytes */
  size = sel->numVisRead * len * sizeof(ofloat); 

  /* Read */
  retCode = ObitFileRead (in->myFile, wantPos, size, 
			  (gchar*)IOBuff, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  in->filePos = in->myFile->filePos; /* remember current file position */

  /* if compressed data uncompress to output buffer */
  if (compressed) {
    ip = op = 0; /* array pointers */
    for (i=0; i<sel->numVisRead; i++) {
      /* Copy random parameters */
      for (j=0; j<sel->nrparmUC; j++) data[op+j] = IOBuff[ip+j];
      
      /* uncompress data */
      wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
      ObitIOUVAIPSUncompress (desc->ncorr, &IOBuff[ip+desc->nrparm], 
			      wtscl, &data[op+sel->nrparmUC]);
      ip += desc->lrec;   /* index in i/O array */
      op += sel->lrecUC;  /* index in output array */
    } /* end decompression loop */
  } /* end compression */
  
  desc->numVisBuff = sel->numVisRead; /* How many read */
  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSRead */

/**
 * Read data from disk applying selection and any calibration.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff (which
 * may be zero).
 * The number attempted in a read is mySel->numVisRead.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVAIPSReadSelect (ObitIOUVAIPS *in, ofloat *data, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gsize size;
  olong len, i, k, ip, op, need, numVisBuff;
  ObitFilePos wantPos;
  gboolean done, compressed, OK;
  ofloat *workVis, *wtscl, *IOBuff = data;
  gchar *routine = "ObitIOUVAIPSReadSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Make sure access was set correctly */
  if (in->access!=OBIT_IO_ReadCal) {
    Obit_log_error(err, OBIT_Error, 
		   "ObitIOUVAIPSReadSelect: access not ReadCal for %s", 
		   in->name);
    return retCode;
  }
  g_assert (ObitUVCalIsA((ObitUVCal*)in->myCal));
  g_assert (data != NULL);

  desc = in->myDesc; /* UV descriptor pointer */
  sel  = in->mySel;  /* selector pointer */
  desc->numVisBuff = 0; /* no data in buffer yet */

  /* what next ? */
  done = ObitIOUVAIPSNext (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* check if done - all visibilities read */
  if (done) {
    /* ObitIOUVAIPSClose (in, err); Close */
    return OBIT_IO_EOF;
  }

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;

  /* Always use compression buffer as the output vis may have a 
     different size from the input */
  if (in->compBuff==NULL) {
    in->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
    in->compBuff = ObitMemAllocName (in->compBuffSize, "UVAIPS comp. buf");
    
    /* check that weight and scale are available if decompressing */
    if (compressed && (desc->ilocws < 0)) {
      Obit_log_error(err, OBIT_Error, 
		     "UV data does not have weight and scale %s", 
		     in->name);
      return retCode;
    }
  }

  IOBuff = in->compBuff; /* Use for I/O */
  /* Make sure buffer large enough (sel->nVisPIO) */
  need = desc->lrec*sel->numVisRead*sizeof(ofloat);
  if (need > in->compBuffSize) {
    Obit_log_error(err, OBIT_Error, 
		   "Decompression buffer ( %d) too small, need %d for %s", 
		   in->compBuffSize, need, in->name);
    return retCode;
  }
  
  /* Visibility decompression buffer */
  if (compressed && (in->decompVis==NULL)) 
    /* Add some extra in buffer - sometimes overrun */
    in->decompVis = ObitMemAllocName ((sel->nrparmUC+5+3*desc->lrec)*sizeof(ofloat), "UVAIPS comp. vis");


  len = desc->lrec; /* How big is a visibility */

  /* read block of sel->numVisRead visibilities at a time  */
  /* get file position offset */
  wantPos = (ObitFilePos)(desc->firstVis-1) * desc->lrec * sizeof(ofloat);
  /* transfer size in bytes */
  size = sel->numVisRead * len * sizeof(ofloat); 

  /* Read */
  retCode = ObitFileRead (in->myFile, wantPos, size, (gchar*)IOBuff, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  in->filePos = in->myFile->filePos; /* remember current file position */

  /* uncompress/calibrate/edit/select transform... to output */
  ip = op = 0;          /* array pointers */
  numVisBuff = 0;       /* How many valid visibilities */
  for (i=0; i<sel->numVisRead; i++) {
 
    /* Decompress */
    if (compressed) {
      /* copy random parameters to visibility work buffer */
      for (k=0; k<sel->nrparmUC; k++) in->decompVis[k] = IOBuff[ip+k];
      
      /* uncompress data */
      wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
      ObitIOUVAIPSUncompress (desc->ncorr, &IOBuff[ip+desc->nrparm], 
			      wtscl, &in->decompVis[sel->nrparmUC]);
      workVis = in->decompVis; /* working visibility pointer */

    } else {
      /* Data not compressed - work out of I/O buffer */
      workVis = &IOBuff[ip];
    }
    
    /* Calibrate and transform */
    OK = ObitUVCalApply ((ObitUVCal*)in->myCal, workVis, &data[op], err);
    if (err->error) /* add traceback,return on error */
      Obit_traceback_val (err, routine, in->name, retCode);

    ip += desc->lrec;   /* index in i/O array */
    if (OK) { /* at least some of the data unflagged - increment output */
      op += sel->lrecUC;  /* index in output array */
      numVisBuff++;       /* count number */
    }
  } /* end compression */

  desc->numVisBuff =  numVisBuff; /* How many good */
  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSReadSelect */

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
ObitIOCode ObitIOUVAIPSWrite (ObitIOUVAIPS *in, ofloat *data, 
				 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gsize size;
  olong len, i, j, ip, op, need;
  ObitFilePos wantPos;
  gboolean compressed;
  ofloat *wtscl, *IOBuff = data;
  gchar *routine = "ObitIOUVAIPSWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);

  desc = in->myDesc; /* UV descriptor pointer */
  sel  = in->mySel;  /* selector pointer */

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;
  if (compressed) {
    /* buffer if necessary create */
    if (in->compBuff==NULL) {
      in->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
      in->compBuff = ObitMemAllocName (in->compBuffSize, "UVAIPS comp. buf");

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
		     "Compression buffer ( %d) too small, need %d for %s", 
		     in->compBuffSize, need, in->name);
      return retCode;
    }
  } /* end of compressed data set up */

  len = desc->lrec; /* How big is a visibility */
  size = desc->numVisBuff * len * sizeof(ofloat); /* transfer size in bytes */

  /* write block of sel->nVisPIO visibilities at a time  */
  /* get file position offset  */
  wantPos = (ObitFilePos)(desc->firstVis-1) * desc->lrec * sizeof(ofloat);

  /* if compressed data uncompress to output buffer */
  if (compressed) {
    ip = op = 0; /* array pointers */
    for (i=0; i<desc->numVisBuff; i++) {
      /* Copy random parameters */
      for (j=0; j<sel->nrparmUC; j++) IOBuff[ip+j] = data[op+j];
      
      /* compress data */
      wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
      ObitIOUVAIPSCompress (desc->ncorr, &data[op+sel->nrparmUC], 
			    wtscl, &IOBuff[ip+desc->nrparm]);
      ip += desc->lrec;   /* index in i/O array (compressed) */
      op += sel->lrecUC;  /* index in input array */
    } /* end compression loop */
  }

  /* Write */
  retCode = ObitFileWrite (in->myFile, wantPos, size, (gchar*)IOBuff, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  in->filePos = in->myFile->filePos; /* remember current file position */

  /* keep track of number of visibilities */
  desc->nvis = MAX (desc->nvis, desc->firstVis+desc->numVisBuff-1);

  /* where will the next write start */
  desc->firstVis += desc->numVisBuff;

  in->myStatus = OBIT_Modified; /* file has been modified */

  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSWrite */

/**
 * Read data Descriptor data from disk.
 * \param in Pointer to object with ObitUVDesc to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitIOUVAIPSReadDescriptor (ObitIOUVAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  gchar *HeaderFile, keyName[9], blob[256];
  AIPSint buffer[260];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, j, ip, ndo, nkey, ikey, nrec;
  gsize size;
  ObitFilePos wantPos;
  ObitFile *myFile=NULL;
  ObitInfoType keyType;
  gchar *routine = "ObitIOUVAIPSReadDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* debug */
  for (i=0; i<256; i++) buffer[i] = 0;

  desc = in->myDesc; /* UV descriptor pointer */

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
  if ((retCode!=OBIT_IO_OK) || (err->error)) {
    if (retCode==OBIT_IO_EOF) 
      Obit_log_error(err, OBIT_Error, 
		     "Empty AIPS header file for %s", in->name);
    Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Convert to internal structure */
  ObitAIPSCatUVGetDesc (in->myDesc, (gchar*)buffer, err);
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
  ((ObitUVDesc*)in->myDesc)->info = ObitInfoListUnref (((ObitUVDesc*)in->myDesc)->info);
  ((ObitUVDesc*)in->myDesc)->info = (gpointer)newObitInfoList ();
  desc = in->myDesc; /* Table descriptor pointer */
    
  if (nkey>0) {
    wantPos = 256 * sizeof(AIPSint); /* File location */
    size    = 256 * sizeof(AIPSint);

    /* loop reading and parsing */
    ikey = 0;
    for (i=0; i<nrec; i++) {
      /* read block */
      retCode = ObitFileRead (myFile, wantPos, size, (gchar*)buffer, err);
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
	  Obit_traceback_val (err, routine, in->name, retCode);
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
} /* end ObitIOUVAIPSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in Pointer to object with ObitUVDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVAIPSWriteDescriptor (ObitIOUVAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gchar *HeaderFile, keyName[FLEN_KEYWORD+1], *keyNameP, blob[256], *ctemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, ii, j, k, ip, ndo, nkey, ikey, nrec;
  oint oitemp[4]={0,0,0,0};
  ObitFilePos wantPos;
  gsize size;
  gboolean doFill;
  AIPSint buffer[256];
  ObitFile *myFile=NULL;
  ObitAIPSDirCatEntry *dirEntry = NULL;
  ObitInfoType keyType;
  gchar *routine = "ObitIOUVAIPSWriteDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* enforce descriptor defaults */
  desc = in->myDesc; /* UV descriptor pointer */
  sel = in->mySel;   /* UV selector pointer */

  g_assert (ObitIsA(desc,  ObitUVDescGetClass()));
  g_assert (ObitInfoListIsA (desc->info));

  ObitUVSelDefault(desc, sel);

  /* debug */
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
  } /* end of read old section */
    
  /* convert descriptor to header */
  retCode = OBIT_IO_ReadErr;
  /* Get catalog descriptor */
  dirEntry = ObitAIPSDirGetEntry(in->disk, in->UserId, in->CNO, err);
  if (err->error)  /* add trace and return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* do conversion */
  ObitAIPSCatUVSetDesc (in->myDesc, (gchar*)buffer, !myFile->exist, 
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
  if (desc->info) nkey = desc->info->number+1; /* How many keywords? */
  else nkey = 1;
  nrec = 1 + (nkey/(256/5)); /* How many blocks? */
  /* Must write at least the first record */

  wantPos = 256 * sizeof(AIPSint); /* File location */
  size    = 256 * sizeof(AIPSint);

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
	ikey++;
	ObitInfoListGetNumber(desc->info, ikey-1, &keyNameP, &keyType, dim, blob, err);
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
      }
    } /* end loop filling block */

    /* write block */
    retCode = ObitFileWrite (myFile, wantPos, size, (gchar*)buffer, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    wantPos = -1L; /* now sequential */
  }
  /* end write keywords section */

  /* flush/close file */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* delete */
  myFile = ObitFileUnref(myFile);
  
  return OBIT_IO_OK;
} /* end ObitIOUVAIPSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * File padded out to integral number of AIPS blocks to keep Moma AIPS happy.
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOUVAIPSFlush (ObitIOUVAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitFilePos wantPos;
  olong size;
  ofloat *padd;
  gchar *routine="ObitIOUVAIPSFlush";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Wonky padding at end of file if needed */
  wantPos = ObitAIPSUVWonkyPad (in->myFile->filePos);
  if (wantPos > in->myFile->filePos) {
    size = 1024;
    padd = g_malloc0(size*sizeof(ofloat));
    retCode = ObitFilePad (in->myFile, wantPos, (gchar*)padd, size, err);
    g_free(padd);
    if ((retCode!=OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */
  }

  retCode = ObitFileFlush(in->myFile, err);

  return retCode;
} /* end ObitIOUVAIPSFlush */

/**
 * Create buffer appropriate for I/O request
 * \param data (output) pointer to data array
 * \param size (output) size of data array in floats.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 */
void 
ObitIOUVAIPSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOUVAIPS *in, ObitInfoList *info, 
			     ObitErr *err)
{
  gchar *name;
  
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
  *size = ObitUVSelBufferSize(in->myDesc, in->mySel);

  /* (re)allocate */
  name =  g_strconcat ("UVBuffer:", in->name, NULL);
  if (*data) *data = ObitMemRealloc (*data, (*size)*sizeof(ofloat));
  else *data = ObitMemAlloc0Name((*size)*sizeof(ofloat), name);
  g_free(name);

} /* end ObitIOUVAIPSCreateBuffer */

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
newObitIOUVAIPSTable (ObitIOUVAIPS *in, ObitIOAccess access, 
			 gchar *tabType, olong *tabVer, ObitErr *err)
{
  ObitTable *out = NULL;
  olong version;
  gboolean gotIt;
  gchar ttype[3], *outName, tabName[51];
  gchar *routine = "newObitIOUVAIPSTable";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  g_assert(tabType!=NULL);
  g_assert(tabVer!=NULL);

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

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
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  /* register it in the TableList */
  ObitTableListPut ((ObitTableList*)in->tableList, tabType, &version, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  /* Force Write to disk 
  ObitIOUVAIPSWriteDescriptor (in, err);*/
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  return (Obit*)out;
} /* end newObitIOUVAIPSTable */

/**
 * Update any disk resident structures about the current tables.
 * \param in   Pointer to object to be updated.
 * \param info ObitInfoList of parent object (not used here).
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOUVAIPSUpdateTables (ObitIOUVAIPS *in, ObitInfoList *info, 
				     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOUVAIPSUpdateTables";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* more complete validity test */
  check (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  retCode = ObitIOUVAIPSWriteDescriptor(in, err);
  if ((retCode!= OBIT_IO_OK) || err->error)
    Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitIOUVAIPSUpdateTables */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOUVAIPSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOUVAIPSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOUVAIPSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOUVAIPSClassInfoDefFn (gpointer inClass)
{
  ObitIOUVAIPSClassInfo *theClass = (ObitIOUVAIPSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOUVAIPSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOUVAIPSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOUVAIPSGetClass;
  theClass->newObit    = NULL;
  theClass->newObitIO  = (newObitIOFP)newObitIOUVAIPS;
  theClass->ObitIOSame = (ObitIOSameFP)ObitIOUVAIPSSame;
  theClass->ObitIORename  = (ObitIORenameFP)ObitIOUVAIPSRename;
  theClass->ObitIOZap  = (ObitIOZapFP)ObitIOUVAIPSZap;
  theClass->ObitCopy   = (ObitCopyFP)ObitIOUVAIPSCopy;
  theClass->ObitClone  = NULL;
  theClass->ObitClear  = (ObitClearFP)ObitIOUVAIPSClear;
  theClass->ObitInit   = (ObitInitFP)ObitIOUVAIPSInit;
  theClass->ObitIOOpen = (ObitIOOpenFP)ObitIOUVAIPSOpen;
  theClass->ObitIOClose= (ObitIOCloseFP)ObitIOUVAIPSClose;
  theClass->ObitIOSet  = (ObitIOSetFP)ObitIOUVAIPSSet;
  theClass->ObitIORead = (ObitIOReadFP)ObitIOUVAIPSRead;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOUVAIPSReadSelect;
  theClass->ObitIOWrite   = 
    (ObitIOWriteFP)ObitIOUVAIPSWrite;
  theClass->ObitIOFlush   = 
    (ObitIOFlushFP)ObitIOUVAIPSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOUVAIPSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOUVAIPSWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOUVAIPSCreateBuffer;
  theClass->ObitIOFreeBuffer   = 
    (ObitIOFreeBufferFP)ObitIOFreeBuffer;
  theClass->newObitIOTable   = 
    (newObitIOTableFP)newObitIOUVAIPSTable;
  theClass->ObitIOUpdateTables   = 
    (ObitIOUpdateTablesFP)ObitIOUVAIPSUpdateTables;

} /* end ObitIOUVAIPSClassDefFn */

/*--------------- Private functions --------------------------*/

/**
 * Creates empty member objects.
 * for each parent class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOUVAIPSInit  (gpointer inn)
{
  const ObitClassInfo *ParentClass;
  ObitIOUVAIPS *in = inn;

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
  in->decompVis    = NULL;
  in->compBuff     = NULL;
  in->compBuffSize = 0;
} /* end ObitIOUVAIPSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOUVAIPSClear (gpointer inn)
{
  ObitIOUVAIPS *in = inn;
  const ObitClassInfo *ParentClass;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->myStatus==OBIT_Active) ||(in->myStatus==OBIT_Modified)) {
    err = newObitErr();
    ObitIOUVAIPSClose (in, err); 
    if (err->error) ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  if (in->AIPSFileName) g_free(in->AIPSFileName); 
  in->AIPSFileName = NULL;
  if (in->myFile) ObitFileUnref (in->myFile);
  if (in->compBuff)  in->compBuff  = ObitMemFree (in->compBuff); 
  if (in->decompVis) in->decompVis = ObitMemFree (in->decompVis);

 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && (ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOUVAIPSClear */

/**
 * Uses selector member to decide which visibilities to
 * read next.
 * Leaves values in myDesc->firstVis and mySel->numVisRead.
 * \param  in Pointer to the object.
 * \return TRUE is finished, else FALSE
 */
static gboolean ObitIOUVAIPSNext (ObitIOUVAIPS *in, ObitErr *err)
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
} /* end ObitIOUVAIPSNext */

/**
 * Compresses UV into scaled shorts
 * Compressed data stores a common weigh and scaling factors as 
 * random parameters and the real and imaginary parts as scaled shorts.
 * If the first short of a pair is -32767 the value is considered invalid.
 * \param  ncorr  Number of weighted complex numbers
 * \param  visin  Expanded visibility array
 * \param  wtscl  (out) Weight and Scale needed to uncompress.
 * \param  visout (out) Compressed visibility array.
 */
static void 
ObitIOUVAIPSCompress (gint ncorr, const ofloat *visin, ofloat *wtscl, 
		      ofloat *visout)
{
  olong i;
  ofloat maxwt, maxvis, scl;
  gshort *packed = (gshort*)visout;

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
	packed[i*2] = (gshort)(scl*visin[i*3] + 0.5);
      else
	packed[i*2] = (gshort)(scl*visin[i*3] - 0.5);
      if (visin[i*3+1] > 0.0)
	packed[i*2+1] = (gshort)(scl*visin[i*3+1] + 0.5);
      else
	packed[i*2+1] = (gshort)(scl*visin[i*3+1] - 0.5);
    } else { /* flag */
      packed[i*2]   = -32767;
      packed[i*2+1] = -32767;
    }
  }
} /* end ObitIOUVAIPSCompress */

/**
 * Uncompresses UV from scaled shorts.
 * Compressed data stores a common weigh and scaling factors as 
 * random parameters and the real and imaginary parts as scaled shorts.
 * If the first short of a pair is -32767 the value is considered invalid.
 * \param  ncorr  Number of weighted complex numbers
 * \param  visin  Compressed visibility array.
 * \param  wtscl  Weight and Scale needed to uncompress.
 * \param  visout (out) Expanded visibility array.
 */
static void 
ObitIOUVAIPSUncompress (gint ncorr, const ofloat *visin, 
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
    if (packed[i*2] == -32767) { /* Flagged */
      visout[i*3]   = 0.0;
      visout[i*3+1] = 0.0;
      visout[i*3+2] = 0.0;
    } else { /* OK */
      visout[i*3]   = scl * packed[i*2];
      visout[i*3+1] = scl * packed[i*2+1];
      visout[i*3+2] = wt;
    }
  }
} /* end ObitIOUVAIPSUncompress */

/**
 * Check validity of object
 * \param in  ObitIO for test
 * \param err ObitErr for reporting errors.
 */
static void check (ObitIOUVAIPS *in, ObitErr *err)
{
  if (in->myFile) {
    if (!ObitFileIsA(in->myFile)) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR myFile member corrupted on %s", in->name);
    }
  }
  /* Check Table List */
  if (in->tableList) {
    ObitTableListCheck ((ObitTableList*)in->tableList, err);
  }
} /* end check */

