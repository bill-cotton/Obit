/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
ObitIOUVAIPSCompress (olong ncorr, const ofloat *visin, ofloat *wtscl, 
		      ofloat *visout);

/** Private: Uncompress visibilities. */
static void 
ObitIOUVAIPSUncompress (olong ncorr, const ofloat *visin, 
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
 * \return return code, OBIT_IO_OK => OK
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
		   "%s: access not ReadCal for %s", 
		   routine, in->name);
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
  } /* end loop over visibilities */

  desc->numVisBuff =  numVisBuff; /* How many good */
  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSReadSelect */

/**
 * Read data from disk to multiple buffers.
 * All buffers must be the same size and the underlying dataset the same.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities attempted is mySel->numVisRead; 
 * actual value saved as myDesc->numVisBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data)
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitIOUVAIPSReadMulti (olong nBuff, ObitIOUVAIPS **in, ofloat **data, 
		       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gsize size;
  olong len, i, j, ip, op, need, ib;
  ObitFilePos wantPos;
  gboolean done, compressed;
  ofloat *wtscl, *IOBuff = data[0];
  gchar *routine = "ObitIOUVAIPSReadMulti";

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  desc = in[0]->myDesc; /* UV descriptor pointer */
  sel  = in[0]->mySel;  /* selector pointer */

  /* what next ? */
  done = ObitIOUVAIPSNext (in[0], err);
  if (err->error) Obit_traceback_val (err, routine, in[0]->name, retCode);

  /* check if done - all visibilities read */
  if (done) {
    /* ObitIOUVAIPSClose (in[0], err);  Close */
    desc->numVisBuff = 0; /* no data in buffer*/
    return OBIT_IO_EOF;
  }

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;
  if (compressed) {
    /* buffer if necessary create */
    if (in[0]->compBuff==NULL) {
      in[0]->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
      in[0]->compBuff = ObitMemAllocName (in[0]->compBuffSize, "UVAIPS comp. buf");

      /* check that weight and scale are available */
      if (desc->ilocws < 0) {
       Obit_log_error(err, OBIT_Error, 
		     "UV data does not have weight and scale %s", 
		      in[0]->name);
       return retCode;
     }
    }
    IOBuff = in[0]->compBuff; /* Use for I/O */
    /* Make sure compression buffer large enough (sel->nVisPIO) */
    need = desc->lrec*sel->numVisRead*sizeof(ofloat); 
    if (need > in[0]->compBuffSize) {
      Obit_log_error(err, OBIT_Error, 
		     "Decompression buffer ( %d) too small, need %d for %s", 
		     in[0]->compBuffSize, need, in[0]->name);
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
  retCode = ObitFileRead (in[0]->myFile, wantPos, size, 
			  (gchar*)IOBuff, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in[0]->name, retCode);
  in[0]->filePos = in[0]->myFile->filePos; /* remember current file position */
  
  /* if compressed data uncompress to output buffer */
  if (compressed) {
    ip = op = 0; /* array pointers */
    for (i=0; i<sel->numVisRead; i++) {
      /* Copy random parameters */
      for (j=0; j<sel->nrparmUC; j++) data[0][op+j] = IOBuff[ip+j];
      
      /* uncompress data */
      wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
      ObitIOUVAIPSUncompress (desc->ncorr, &IOBuff[ip+desc->nrparm], 
			      wtscl, &data[0][op+sel->nrparmUC]);
      ip += desc->lrec;   /* index in i/O array */
      op += sel->lrecUC;  /* index in output array */
    } /* end loop decompressing vis */
  }  /* end decompression */
  
  /* transfer size in bytes */
  size = sel->numVisRead * sel->lrecUC * sizeof(ofloat); 

  /* Loop copying buffers */
  for (ib=1; ib<nBuff; ib++) {
    memcpy (data[ib], data[0], size);
  } /* end loop over buffer */
  
  desc->numVisBuff = sel->numVisRead; /* How many read */
  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSReadMulti */

/**
 * Reread data from disk to multiple buffers.
 * Retreives data read in a previous call to ObitIOUVAIPSReadMulti
 * NOTE: this depends on retreiving the data from the first element in 
 * in which should be the same as in the call to ObitIOUVAIPSReadMulti
 * All buffers must be the same size and the underlying dataset the same.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities attempted is mySel->numVisRead; 
 * actual value saved as myDesc->numVisBuff.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data)
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitIOUVAIPSReReadMulti (olong nBuff, ObitIOUVAIPS **in, ofloat **data, 
			 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gsize size;
  olong ib;
  /* gchar *routine = "ObitIOUVAIPSReReadMulti";*/

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  desc = in[0]->myDesc; /* UV descriptor pointer */
  sel  = in[0]->mySel;  /* selector pointer */

  /* transfer size in bytes */
  size = desc->numVisBuff * sel->lrecUC * sizeof(ofloat); 

  /* Loop copying buffers */
  for (ib=1; ib<nBuff; ib++) {
    memcpy (data[ib], data[0], size);
  } /* end loop over buffer */
  
  desc->numVisBuff = sel->numVisRead; /* How many read */
  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSReReadMulti */

/**
 * Read data from disk applying selection to multiple buffers 
 * and any (output dependent) calibration.
 * If amp/phase calibration being applied, it is done independently
 * for each buffer using the myCal in the associated in,
 * otherwise, the first buffer is processed and copied to the others.
 * All selected buffer sizes etc must be the same.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff (which
 * may be zero).
 * The number attempted in a read is mySel->numVisRead.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitIOUVAIPSReadMultiSelect (olong nBuff, ObitIOUVAIPS **in, ofloat **data, 
			     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  gsize size;
  olong i, k, ip, op, ib,  need, numVisBuff, nFull;
  ObitFilePos wantPos;
  gboolean done, compressed, OK=FALSE;
  ofloat *workVis, *wtscl, *IOBuff = NULL;
  gchar *routine = "ObitIOUVAIPSReadMultiSelect";

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  /* Make sure access was set correctly */
  if (in[0]->access!=OBIT_IO_ReadCal) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: access not ReadCal for %s", 
		   routine, in[0]->name);
    return retCode;
  }
  desc = in[0]->myDesc; /* UV descriptor pointer */
  sel  = in[0]->mySel;  /* selector pointer */
  desc->numVisBuff = 0; /* no data in buffer yet */

  /* what next ? */
  done = ObitIOUVAIPSNext (in[0], err);
  if (err->error) Obit_traceback_val (err, routine, in[0]->name, retCode);

  /* check if done - all visibilities read */
  if (done) {
    /* ObitIOUVAIPSClose (in[0], err); Close */
    return OBIT_IO_EOF;
  }

  /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;

  /* Always use compression buffer as the output vis may have a 
     different size from the input */
  if (in[0]->compBuff==NULL) {
    in[0]->compBuffSize = desc->lrec*sel->nVisPIO*sizeof(ofloat);
    in[0]->compBuff = ObitMemAllocName (in[0]->compBuffSize, "UVAIPS comp. buf");
    
    /* check that weight and scale are available if decompressing */
    if (compressed && (desc->ilocws < 0)) {
      Obit_log_error(err, OBIT_Error, 
		     "UV data does not have weight and scale %s", 
		     in[0]->name);
      return retCode;
    }
  }

  IOBuff = in[0]->compBuff; /* Use for I/O */
  /* Make sure buffer large enough (sel->nVisPIO) */
  need = desc->lrec*sel->numVisRead*sizeof(ofloat);
  if (need > in[0]->compBuffSize) {
    Obit_log_error(err, OBIT_Error, 
		   "Decompression buffer ( %d) too small, need %d for %s", 
		   in[0]->compBuffSize, need, in[0]->name);
    return retCode;
  }
  
  /* Visibility decompression buffer */
  if (compressed && (in[0]->decompVis==NULL)) 
    /* Add some extra in buffer - sometimes overrun */
    in[0]->decompVis = 
      ObitMemAllocName ((sel->nrparmUC+5+3*desc->lrec)*sizeof(ofloat), 
			"UVAIPS comp. vis");

  /* read block of sel->numVisRead visibilities at a time  */
  /* get file position offset */
  wantPos = (ObitFilePos)(desc->firstVis-1) * desc->lrec * sizeof(ofloat);
  /* transfer size in bytes */
  size = sel->numVisRead * desc->lrec * sizeof(ofloat); 

  /* Read use work buffers on first in */
  retCode = ObitFileRead (in[0]->myFile, wantPos, size, (gchar*)IOBuff, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in[0]->name, retCode);
  in[0]->filePos = in[0]->myFile->filePos; /* remember current file position */

  /* How many buffers need the full treatment?
     none unless sel->doCal - just copy if all the same */
  if (sel->doCal) nFull = nBuff;
  else  nFull = 1;

  /* uncompress/calibrate/edit/select transform... to output */
  ip = op = 0;          /* array pointers */
  numVisBuff = 0;       /* How many valid visibilities */
  for (i=0; i<sel->numVisRead; i++) {

    /* Loop over buffers being calibrated independently */
    for (ib=0; ib<nFull; ib++) {
      
      /* Decompress */
      if (compressed) {
	/* copy random parameters to visibility work buffer */
	for (k=0; k<sel->nrparmUC; k++) in[0]->decompVis[k] = IOBuff[ip+k];
	
	/* uncompress data */
	wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
	ObitIOUVAIPSUncompress (desc->ncorr, &IOBuff[ip+desc->nrparm], 
				wtscl, &in[0]->decompVis[sel->nrparmUC]);
	workVis = in[0]->decompVis; /* working visibility pointer */
	
      } else {
	/* Data not compressed - work out of I/O buffer */
	workVis = &IOBuff[ip];
      }
      
      /* Calibrate and transform to individual buffer */
      OK = ObitUVCalApply ((ObitUVCal*)in[ib]->myCal, workVis, &data[ib][op], err);
      if (err->error) /* add traceback,return on error */
	Obit_traceback_val (err, routine, in[ib]->name, retCode);
      
    } /* end loop over buffer */
    ip += desc->lrec;   /* index in i/O array */
    if (OK) { /* at least some of the data unflagged - increment output */
      op += sel->lrecUC;  /* index in output array */
      numVisBuff++;       /* count number */
    }
  } /* end loop over visibilities */
  
  /* transfer size in bytes */
  size = numVisBuff * sel->lrecUC * sizeof(ofloat); 

  /* May need to copy the rest */
  for (ib=nFull; ib<nBuff; ib++) {
    memcpy (data[ib], data[0], size);
  }
  
  desc->numVisBuff =  numVisBuff; /* How many good */
  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSReadMultiSelect */

/**
 * Reread data from disk applying selection to multiple buffers 
 * and any (output dependent) calibration.
 * Retreives data read in a previous call to ObitIOUVAIPSReadMultiSelect
 * possibly applying new calibration.
 * If amp/phase calibration being applied, it is done independently
 * for each buffer using the myCal in the associated in,
 * otherwise, data from the first buffer copied to the others.
 * NOTE: this depends on retrieving the data from the first element in 
 * in which should be the same as in the call to ObitIOUVAIPSReadMultiSelect
 * All selected buffer sizes etc must be the same.
 * The number read will be mySel->nVisPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstVis
 * and the number of visibilities is myDesc->numVisBuff (which
 * may be zero).
 * The number attempted in a read is mySel->numVisRead.
 * When OBIT_IO_EOF is returned all data has been read (then is no new
 * data in data) and the I/O has been closed.
 * \param nBuff Number of buffers to be filled
 * \param in    Array of pointers to to object to be read; 
 *              must all be to same underlying data set but with 
 *              independent calibration
 * \param data  array of pointers to buffers to write results.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitIOUVAIPSReReadMultiSelect (olong nBuff, ObitIOUVAIPS **in, ofloat **data, 
			       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitUVDesc* desc;
  ObitUVSel* sel;
  olong i, k, ip, op, ib, numVisBuff, nFull;
  gboolean compressed, OK=FALSE;
  ofloat *workVis, *wtscl, *IOBuff = NULL;
  gsize size;
  gchar *routine = "ObitIOUVAIPSReReadMultiSelect";

  /* error checks */
  if (err->error) return retCode;
  if (nBuff<=0)   return retCode;

  /* Make sure access was set correctly */
  if (in[0]->access!=OBIT_IO_ReadCal) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: access not ReadCal for %s", 
		   routine, in[0]->name);
    return retCode;
  }
  desc = in[0]->myDesc; /* UV descriptor pointer */
  sel  = in[0]->mySel;  /* selector pointer */

  IOBuff = in[0]->compBuff; /* Use for I/O */

  /* How many buffers need the full treatment?
     none unless sel->doCal - just copy if all the same */
  if (sel->doCal) nFull = nBuff;
  else  nFull = 0;

   /* Is the data compressed? If so need to uncompress */
  compressed = sel->Compress;

  /* uncompress/calibrate/edit/select transform... to output */
  ip = op = 0;          /* array pointers */
  if (nFull>0) numVisBuff = 0;       /* How many valid visibilities */
  else numVisBuff         = desc->numVisBuff;

  if (nFull>0) { /* calibrate something? */
    for (i=0; i<sel->numVisRead; i++) {
      
      /* Loop over buffers being calibrated independently */
      for (ib=0; ib<nFull; ib++) {
	
	/* Decompress */
	if (compressed) {
	  /* copy random parameters to visibility work buffer */
	  for (k=0; k<sel->nrparmUC; k++) in[0]->decompVis[k] = IOBuff[ip+k];
	  
	  /* uncompress data */
	  wtscl = &IOBuff[ip+desc->ilocws]; /* weight and scale array */
	  ObitIOUVAIPSUncompress (desc->ncorr, &IOBuff[ip+desc->nrparm], 
				  wtscl, &in[0]->decompVis[sel->nrparmUC]);
	  workVis = in[0]->decompVis; /* working visibility pointer */
	  
	} else {
	  /* Data not compressed - work out of I/O buffer */
	  workVis = &IOBuff[ip];
	}
	
	/* Calibrate and transform to individual buffer */
	OK = ObitUVCalApply ((ObitUVCal*)in[ib]->myCal, workVis, &data[ib][op], err);
	if (err->error) /* add traceback,return on error */
	  Obit_traceback_val (err, routine, in[ib]->name, retCode);
	
      } /* end loop over buffer */
      ip += desc->lrec;   /* index in i/O array */
      if (OK) { /* at least some of the data unflagged - increment output */
	op += sel->lrecUC;  /* index in output array */
	numVisBuff++;       /* count number */
      }
    } /* end loop over visibilities */
  } /* end if calibrate something */

  /* transfer size in bytes */
  size = numVisBuff * sel->lrecUC * sizeof(ofloat); 

  /* May need to copy the rest */
  nFull = MAX (1, nFull); /* Don't copy to itself */
  for (ib=nFull; ib<nBuff; ib++) {
    memcpy (data[ib], data[0], size);
  }
  
  desc->numVisBuff =  numVisBuff; /* How many good */
  return  OBIT_IO_OK;
} /* end ObitIOUVAIPSReReadMultiSelect */

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
 * \li xxxUser  OBIT_oint    AIPS User number
 * \li xxxCNO   OBIT_oint    AIPS Catalog slot number
 * \li xxxDir   OBIT_string  Directory name for xxxDisk
 *
 * For all File types types:
 * \li xxxDataType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxFileType OBIT_oint File type as ObitIOType, OBIT_IO_FITS, OBIT_IO_AIPS
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
void ObitIOUVAIPSGetFileInfo (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
			      ObitInfoList *outList, ObitErr *err)
{
  ObitInfoType type;
  ObitAIPSDirCatEntry *entry=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *FileType="UV", *DataType="AIPS", *dirname;
  olong disk, user, cno, seq, i;
  gpointer listPnt;
  gchar *parm[] = {"nVisPIO", "DoCalSelect", "Stokes", "BChan", "EChan", "BIF", "EIF",
		   "doPol", "doCalib", "gainUse", "flagVer", "BLVer", "BPVer",
		   "Subarray", "dropSubA", "FreqID", "timeRange", "UVRange",
		   "InputAvgTime", "Sources", "souCode", "Qual", "Antennas",
		   "corrType", "passAll", "doBand", "Smooth", "SubScanTime",
		   NULL};
  gchar *routine = "ObitIOUVAIPSGetFileInfo";

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
} /* end ObitIOUVAIPSGetFileInfo */

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
  theClass->ObitIOReReadMulti = (ObitIOReReadMultiFP)ObitIOUVAIPSReReadMulti;
  theClass->ObitIOReadMulti = (ObitIOReadMultiFP)ObitIOUVAIPSReadMulti;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOUVAIPSReadSelect;
  theClass->ObitIOReadMultiSelect = 
    (ObitIOReadMultiSelectFP)ObitIOUVAIPSReadMultiSelect;
  theClass->ObitIOReReadMultiSelect = 
    (ObitIOReReadMultiSelectFP)ObitIOUVAIPSReReadMultiSelect;
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
  theClass->ObitIOGetFileInfo   =
    (ObitIOGetFileInfoFP)ObitIOUVAIPSGetFileInfo;

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
ObitIOUVAIPSCompress (olong ncorr, const ofloat *visin, ofloat *wtscl, 
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
ObitIOUVAIPSUncompress (olong ncorr, const ofloat *visin, 
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

