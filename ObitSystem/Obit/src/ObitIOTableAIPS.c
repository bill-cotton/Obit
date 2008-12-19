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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include <stdio.h>
#include "ObitIOTableAIPS.h"
#include "ObitAIPSCat.h"
#include "ObitAIPS.h"
#include "ObitMem.h"
#include "ObitUVDesc.h"
#include "ObitImageDesc.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOTableAIPS.c
 * ObitIOTableAIPS class function definitions.
 */

/*------------------- file globals ------------------------------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOTableAIPS";

/** Function to obtain parent ClassInfo - ObitIO */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/** List of AIPS "extensions which are not Tables */
gchar *NotAIPSTable[] = {"HI", "PL", "SL", NULL};

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOTableAIPSClassInfo myClassInfo = {FALSE}; 

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOTableAIPSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOTableAIPSClear (gpointer in);

/** Private: Determine next row to read */
static gboolean ObitIOTableAIPSNext (ObitIOTableAIPS *in, olong rowno);


/** Private: Is this an AIPS table? */
static gboolean ObitIOTableAIPSIsTable (ObitIOTableAIPS *in);

/** Private: Set Class function pointers. */
static void ObitIOTableAIPSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class on the first call.
 * \param name An optional name for the object.
 * \param info if non-NULL it is used to initialize the new object.
 * \param err  ObitErr for error messages.
 * \return the new object, NULL on failure.
 */
ObitIOTableAIPS* newObitIOTableAIPS (gchar *name, ObitInfoList *info,
				     ObitErr *err)
{
  ObitIOTableAIPS* out;
  gint32 i, dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "newObitIOTableAIPS";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOTableAIPSClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitIOTableAIPS));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOTableAIPSInit((gpointer)out);

  /* Get any info from info input */
  if (info!=NULL) {
    type = OBIT_long; for (i=0; i<MAXINFOELEMDIM; i++) dim[i] = 1;
    if(!ObitInfoListGet(info, "Disk", &type, (gint32*)dim, 
			&out->disk, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "User", &type, (gint32*)dim, 
			&out->UserId, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "CNO", &type, (gint32*)dim, 
			&out->CNO, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "Ver", &type, (gint32*)dim, 
			&out->tabVer, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    
    if(!ObitInfoListGet(info, "TableType", &type, (gint32*)dim, 
			&out->tabType, err)) /* add traceback on error */
      Obit_traceback_val (err, routine, name, out);
    out->tabType[2] = 0;

  }  /* end of initialize from info */
  
  return out;
} /* end newObitIOTableAIPS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOTableAIPSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOTableAIPSClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIOTableAIPSGetClass */

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
gboolean ObitIOTableAIPSSame (ObitIO *in, ObitInfoList *in1, 
				ObitInfoList *in2, ObitErr *err)
{
  olong CNO1, UserId1, disk1, ver1, CNO2, UserId2, disk2, ver2;
  gchar type1[3], type2[3];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOTableAIPSSame";

  /* error checks */
  if (err->error) return same;
  g_assert (ObitIOTableAIPSIsA(in));

  /* get instructions from info */
  if(!ObitInfoListGet(in1, "Disk", &type, dim, &disk1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in1, "User", &type, dim, &UserId1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in1, "CNO", &type, dim, &CNO1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in1, "TableType", &type, dim, type1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in1, "Ver", &type, dim, &ver1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "Disk", &type, dim, &disk2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "User", &type, dim, &UserId2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "CNO", &type, dim, &CNO2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "TableType", &type, dim, type2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if(!ObitInfoListGet(in2, "Ver", &type, dim, &ver2, err))
    Obit_traceback_val (err, routine, in->name, same);

  /* Compare */
  same = (disk1==disk2) && (CNO1==CNO2) && (UserId1==UserId2) &&
    (ver1==ver2) && (type1[0]==type2[0]) && (type1[1]==type2[1]);

  return same;
} /* end ObitIOTableAIPSSame */

/**
 * Delete underlying files.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitIOTableAIPSZap (ObitIOTableAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitFile *myFile=NULL;
  gchar *routine = "ObitIOTableAIPSZap";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

   /* Close if still open */
  if ((in->myStatus==OBIT_Modified) || (in->myStatus==OBIT_Active)) {
    retCode = ObitIOTableAIPSClose (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
  /* Get Table file */
  myFile = newObitFile("Table to Zap");
  myFile->fileName = ObitAIPSFilename (OBIT_AIPS_Table, in->disk, in->CNO, 
				       in->UserId, in->tabType, in->tabVer, err);
  myFile = ObitFileZap(myFile, err);
  myFile = ObitFileUnref(myFile);

  /* Trace back if error detected */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

 return;
} /* end ObitIOTableAIPSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIOTableAIPS* ObitIOTableAIPSCopy  (ObitIOTableAIPS *in, 
				       ObitIOTableAIPS *out, ObitErr *err)
{
  const ObitIOClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  g_assert (err!=NULL);
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) return out;

   /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitIOTableAIPS(outName, NULL, err);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy this class */
  out->disk       = in->disk;
  out->UserId     = in->UserId;
  out->CNO        = in->CNO;
  out->tabVer     = in->tabVer;
  out->tabType[0] = in->tabType[0];
  out->tabType[1] = in->tabType[1];
  if (out->AIPSFileName!=NULL) g_free(out->AIPSFileName);
  out->AIPSFileName = g_strdup(in->AIPSFileName);

  return out;
} /* end ObitIOTableAIPSCopy */

/**
 * Initialize structures and open file.
 * The file etc. info should have been stored in the ObitInfoList.
 * For accessing AIPS files the following entries in the ObitInfoList 
 * are used:
 * \li "Disk" OBIT_int (1,1,1) AIPS "disk" number.
 * \li "User" OBIT_int (1,1,1) user number.
 * \li "CNO"  OBIT_int (1,1,1) AIPS catalog slot number.
 * \li "TableType" OBIT_string (2,1,1) AIPS Table type
 * \li "Ver"  OBIT_int    (1,1,1) AIPS table version number.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite, 
 *               OBIT_IO_WriteOnly).
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSOpen (ObitIOTableAIPS *in, ObitIOAccess access, 
				ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  gsize bsize;
  gchar *routine = "ObitIOTableAIPSOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) 
    {in->myStatus = OBIT_Active; in->access = access; return OBIT_IO_OK;}

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

  if(!ObitInfoListGet(info, "Ver", &type, (gint32*)dim, 
		      &in->tabVer, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  if(!ObitInfoListGet(info, "TableType", &type, (gint32*)dim, 
		      &in->tabType, err)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  in->tabType[2] = 0;

  /* form file name for file */
  if (in->AIPSFileName) g_free(in->AIPSFileName); /* free old */
  in->AIPSFileName = 
    ObitAIPSFilename (OBIT_AIPS_Table, in->disk, in->CNO, 
		      in->UserId, in->tabType, in->tabVer, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* If in->myFile still connected unreference */
  if (in->myFile) ObitFileUnref (in->myFile);
  in->myFile = newObitFile(in->name);  /* new one */

  /* set buffer size */
  bsize = 256*sizeof(AIPSint);

  /* open file */
  retCode = OBIT_IO_OpenErr; /* in case something goes wrong */
  if (ObitFileOpen (in->myFile, in->AIPSFileName, access,  OBIT_IO_Binary,
		     bsize, err) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save some information */
  in->access = access;
  in->myStatus = OBIT_Active;
  
  /* initialize location in table */
  in->filePos = 0;
  if (desc) desc->firstRow = 0;  

  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitIOTableAIPSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSClose (ObitIOTableAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong blksize;
  gchar *routine = "ObitIOTableAIPSClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) 
    {in->myStatus = OBIT_Inactive; return OBIT_IO_OK;}

  /* don't bother if it's not open */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) 
    return OBIT_IO_OK;

  /* Curse you Mama AIPS!!!!!*/
  /* Because of the archaic AIPS file structure, the file must be 
     an exact multiple of 256*sizeof(AIPSint) in size.  On the other
     hand, you don't want to write over existing data so you can't
     just pad out the current block.  Pad the end of the file. */
  blksize = 256*sizeof(AIPSint);
  if ((in->access==OBIT_IO_ReadWrite) || 
      (in->access==OBIT_IO_WriteOnly)) 
    retCode = ObitFilePadFile(in->myFile, blksize, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  if (ObitFileClose (in->myFile, err) || (err->error)) 
    /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  in->myStatus = OBIT_Inactive;
  return OBIT_IO_OK;
} /* end ObitIOTableAIPSClose */

/**
 * initialize I/O
 * Not used - returns.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSSet (ObitIOTableAIPS *in, ObitInfoList *info, 
			       ObitErr *err)
{
  return OBIT_IO_OK;
} /* end ObitIOTableAIPSSet */

/**
 * Read table data from disk and specifying start row.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param rowno Starting row number (1-rel) -1=> next.
 * \param data pointer to buffer to receive results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSReadRow (ObitIOTableAIPS *in, olong rowno, 
				   ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  gsize size;
  olong len, ilen, offset, ioff, iRow, nRows, row;
  ObitFilePos wantPos;
  olong  *idata = (olong*)data;
  gboolean done, rewrite;
  gchar *cdata = (gchar*)data;
  gchar *routine = "ObitIOTableAIPSReadRow";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) return OBIT_IO_OK;
  g_assert (ObitFileIsA(in->myFile));
  g_assert (data != NULL);

  desc = in->myDesc;                /* Table descriptor pointer */
  sel  = in->mySel;                 /* selector pointer */
  ilen = desc->lrow / sizeof(olong); /* Size of row in gints */
  len   = desc->lrow;               /* Size of row in bytes */
  size  = desc->lrow;               /* Size of transfer in bytes */

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
    retCode = ObitIOTableAIPSWrite (in, data, err);
    if (err->error)  /* add trace and return on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  } /* end of rewrite modified data */

  /* what next ? */
  done = ObitIOTableAIPSNext (in, rowno);

  /* check if done - all visibilities read */
  if (done) {
    return OBIT_IO_EOF;
  }

  row   = desc->firstRow;    /* which row to start */
  nRows = desc->numRowBuff;  /* How many rows to do */

  offset = 0; /* offset in buffer */
  ioff   = 0; /* offset in buffer in gints */

  /* read file one row at a time */
  retCode = OBIT_IO_ReadErr; /* in case something goes wrong */
  for (iRow=0; iRow<nRows; iRow++) {

    /* get file position offset (AIPS tables have funny rules) */
    wantPos = ObitAIPSTableFileOffset(desc->startData, desc->lrowIO, row);

    /* Read it */
    retCode = ObitFileRead (in->myFile, wantPos, size, 
			    (gchar*)&cdata[offset], err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) { /* add traceback on error */
      if (retCode==OBIT_IO_EOF)
	Obit_log_error(err, OBIT_Error,
		       "Read EOF at row %d in table %s",
		       row,in->name);
      Obit_traceback_val (err, routine, in->name, retCode);
    }
    in->filePos = in->myFile->filePos; /* remember current file position */

    /* Conversion of status from the ancient AIPSish */
    if ((desc->statusOff>=0) &&(idata[ioff+desc->statusOff]>0))
      idata[ioff+desc->statusOff] = 0;

    row++;                /* next row */
    offset  += len;       /* offset in data buffer in bytes */
    ioff    += ilen;      /* offset in data buffer in gints */
  }  /* end loop reading rows */
  
  return  OBIT_IO_OK;
} /* end ObitIOTableAIPSReadRow */

/**
 * Read table data from disk.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to receive results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSRead (ObitIOTableAIPS *in, ofloat *data, 
				ObitErr *err)
{
  /* Call bitIOTableAIPSReadRow */
  return ObitIOTableAIPSReadRow (in, -1, data, err);
} /* end ObitIOTableAIPSRead */

/**
 * Read table data from disk applying selection and specifying start row.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param rowno Starting row number (1-rel), -1=> next.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSReadRowSelect (ObitIOTableAIPS *in, olong rowno, 
					 ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  gsize size;
  olong len, ilen, offset, ioff, iRow, nRows, row;
  ObitFilePos wantPos;
  olong  *idata = (olong*)data;
  gboolean done, rewrite;
  gchar *cdata = (gchar*)data;
  gchar *routine = "ObitIOTableAIPSReadRowSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) return OBIT_IO_OK;
  g_assert (ObitFileIsA(in->myFile));
  g_assert (data != NULL);

  desc = in->myDesc;                /* Table descriptor pointer */
  sel  = in->mySel;                 /* selector pointer */
  ilen = desc->lrow / sizeof(olong); /* Size of row in gints */
  len  = desc->lrow;                /* Size of row in bytes */
  size = desc->lrow;                /* Size of transfer in bytes */

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
    retCode = ObitIOTableAIPSWrite (in, data, err);
    if (err->error)  /* add trace and return on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  } /* end of rewrite modified data */

  /* what next ? */
  done = ObitIOTableAIPSNext (in, rowno);

  /* check if done - all visibilities read */
  if (done) {
    return OBIT_IO_EOF;
  }

  row   = desc->firstRow;    /* which row to start */
  nRows = desc->numRowBuff;  /* How many rows to do */

  offset = 0; /* offset in buffer in bytes */
  ioff   = 0; /* offset in buffer in gints */

  /* read file one row at a time */
  for (iRow=0; iRow<nRows; iRow++) {

    /* get file position offset (AIPS tables have funny rules) */
    wantPos = ObitAIPSTableFileOffset(desc->startData, desc->lrowIO, row);

    /* Read it */
    retCode = ObitFileRead (in->myFile, wantPos, size, 
			    (gchar*)&cdata[offset], err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */

    /* Conversion of status from the ancient AIPSish */
    if ((desc->statusOff>=0) &&(idata[ioff+desc->statusOff]>0))
      idata[ioff+desc->statusOff] = 0;

    row++;                /* next row */
    offset  += len;       /* offset in data buffer in bytes */
    ioff    += ilen;      /* offset in data buffer in gints */
  }  /* end loop reading rows */
  
  return  OBIT_IO_OK;
} /* end ObitIOTableAIPSReadSelect */

/**
 * Read table data from disk applying selection.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param data pointer to buffer to write results.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSReadSelect (ObitIOTableAIPS *in, ofloat *data, 
				ObitErr *err)
{
  /* Call ObitIOTableAIPSReadRowSelect */
  return ObitIOTableAIPSReadRowSelect (in, -1, data, err);
} /* end ObitIOTableAIPSReadSelect */

/**
 * Write information to disk specifying starting row.
 * Write the desc->numRowBuff rows beginning with desc->firstRow
 * with the data in data.
 * Rows with negative values in the "_status" column are not written.
 * \param in Pointer to object to be written.
 * \param rowno Starting row number (1-rel) -1=> next.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSWriteRow (ObitIOTableAIPS *in, olong rowno, 
				    ofloat *data, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  gsize size;
  olong len, ilen, offset, ioff, iRow, nRows, row;
  ObitFilePos wantPos;
  olong  *idata = (olong*)data;
  gchar *cdata = (gchar*)data;
  gchar *routine = "ObitIOTableAIPSWriteRow";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) return OBIT_IO_OK;
  g_assert (ObitFileIsA(in->myFile));
  g_assert (data != NULL);

  desc = in->myDesc;                /* Table descriptor pointer */
  sel  = in->mySel;                 /* selector pointer */
  ilen = desc->lrow / sizeof(olong); /* Size of row in gints */
  len   = desc->lrow;               /* Size of row in bytes */
  size  = desc->lrow;               /* Size of transfer in bytes */

  /* which row to start , specified or highest? */
  /* Which row, specified or highest? */
  if (rowno>0) desc->firstRow = rowno;
  else desc->firstRow = MAX (1, desc->nrow+1);
  row   = desc->firstRow;    
  nRows = desc->numRowBuff;  /* How many rows to do */

  offset = 0; /* offset in buffer */
  ioff   = 0; /* offset in data array */

  /* get initial file position offset (AIPS tables have funny rules) */
  wantPos = ObitAIPSTableFileOffset(desc->startData, desc->lrowIO, row);

  /* write file one row at a time  */
  for (iRow=0; iRow<nRows; iRow++) {

    /* Conversion of status to AIPSish */
    if ((desc->statusOff>=0) &&(idata[ioff+desc->statusOff]==0))
      idata[ioff+desc->statusOff] = 1;

    /* Write */
    retCode = ObitFileWrite (in->myFile, wantPos, size, 
			    (gchar*)&cdata[offset], err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */

    offset  += len;       /* offset in data buffer in bytes */
    ioff    += ilen;      /* offset in data buffer in gints */
    row++;                /* next row */

    /* get next file position offset (AIPS tables have funny rules) */
    wantPos = ObitAIPSTableFileOffset(desc->startData, desc->lrowIO, row);

    /* Do we need to pad out a block just finished? */
    if ((in->myFile->access == OBIT_IO_WriteOnly) && 
	(wantPos>in->myFile->filePos)) {
      retCode = ObitFilePad (in->myFile, wantPos, (gchar*)data, size, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      in->filePos = in->myFile->filePos; /* remember current file position */
    }
  }  /* end loop writing rows */

  desc->firstRow = row;     /* which row to start next write */
  /* keep track of maximum row number written  */
  desc->nrow = MAX (row-2, desc->nrow); 

  in->myStatus = OBIT_Modified; /* file has been modified */
  
  return  OBIT_IO_OK;
} /* end ObitIOTableAIPSWriteRow */

/**
 * Write information to disk.
 * Write the desc->numRowBuff rows beginning with desc->firstRow
 * with the data in data.
 * Rows with negative values in the "_status" column are not written.
 * \param in Pointer to object to be written.
 * \param data pointer to buffer containing input data.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSWrite (ObitIOTableAIPS *in, ofloat *data, 
				 ObitErr *err)
{
  /* Call bitIOTableAIPSWriteRow */
  return ObitIOTableAIPSWriteRow (in, -1, data, err);
} /* end ObitIOTableAIPSWrite */

/**
 * Read Table Descriptor data from disk.
 * AIPS Table headers are REALLY messy!
 * This routine needs to be coordinated with 
 * #ObitAIPSCat::ObitAIPSCatTableGetDesc.
 * \param in Pointer to object with ObitTableDesc to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode 
ObitIOTableAIPSReadDescriptor (ObitIOTableAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableDesc* desc;
  gchar keyName[9], blob[256], temp[12], *ctemp;
  gint32 cdim, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  AIPSint controlBlock[256], record[256];
  olong i, j, ip, ncol, icol, kkol, ndo, nkey, ikey, titleRec, unitRec, keyRec, nrec;
  olong damn;
  ObitInfoType keyType;
  ObitFilePos wantPos;
  gsize size;
  gchar *routine = "ObitIOTableAIPSReadDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) return OBIT_IO_OK;

  desc = in->myDesc; /* Table descriptor pointer */
  size = 256 * sizeof(AIPSint);

  /* read Control block */
  wantPos = 0; /* File location */
  retCode = ObitFileRead(in->myFile, wantPos, size, 
			 (gchar*)controlBlock, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) {/* add traceback on error */
    if (retCode==OBIT_IO_EOF)
      Obit_log_error(err, OBIT_Error,
		     "AIPS Table is not properly initialized for %s",
		     in->name);
    Obit_traceback_val (err, routine, in->name, retCode);
  }
  in->filePos = in->myFile->filePos; /* remember current file position */
  
  /* Check Recognition string */
  if (strncmp("*AIPS TABLE*", (gchar*)&controlBlock[53], 12)) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR file NOT AIPS Table for %s", in->name);
    return OBIT_IO_ReadErr;
  }
  
  /* Get critical control info */
  ncol     = controlBlock[9];  /* number of data columns in table */
  nkey     = controlBlock[52]; /* Number of keywordvalue pairs */
  titleRec = controlBlock[46]; /* First record (1024 byte) for col. titles */
  unitRec  = controlBlock[47]; /* First record for col. units */
  keyRec   = controlBlock[48]; /* First record for keyword/value pairs */
    
  /* read funky info block */
  retCode = ObitFileRead(in->myFile, -1L, size, (gchar*)record, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  in->filePos = in->myFile->filePos; /* remember current file position */

  /* Convert Control block to internal structure */
  ObitAIPSCatTableGetDesc (in->myDesc, in->tabType, in->tabVer,
			     controlBlock, record, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /*++++++++++++ Read column titles ++++++++++++*/
  /* Order is Physical rather than logical */
  /* Start File location */
  wantPos = (ObitFilePos)(titleRec-1) * 256 * sizeof(AIPSint); 
  
  /* Initialize Title array */
  if (desc->FieldName) { /* get rid of old version first */
    for (i=0; i<ncol+1; i++) 
      if (desc->FieldName[i]) g_free(desc->FieldName[i]);
    g_free(desc->FieldName);
  }
  desc->FieldName = g_malloc0((ncol+1)*sizeof(gchar*));
  desc->firstRow = 0;  /* where reading in table */

  /* loop reading and parsing */
  nrec = 1 + (ncol-1) / (256/6);
  icol = 0;
  for (i=0; i<nrec; i++) {
 
   /* read block */
    retCode = ObitFileRead(in->myFile, wantPos, size, (gchar*)record, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */
    wantPos = -1; /* continue sequential */

    /* Parse - how many? */
    ndo = 42;
    ndo = MIN (ndo, (ncol-icol)); /* Not more than number of columns */
    ip = 0;
    for (j=0; j<ndo; j++) {
      /* Inverse lookup to get logical column */
      kkol = -1;
      for (damn=0; damn<desc->nfield; damn++) 
	if (desc->order[damn]==(icol+1)) kkol = damn;
      if (desc->FieldName[kkol]) g_free(desc->FieldName[kkol]);
      desc->FieldName[kkol] = g_strndup((gchar*)&record[ip], 24);
      icol++;
      ip +=6; /* next title in record */
    }
  } /* end loop reading titles */

  /* Add "_status" columnNot needed? */
  if (desc->FieldName[icol]) g_free(desc->FieldName[icol]);
  desc->FieldName[icol] = g_strndup("_status", 24); 

  /*++++++++++++ Read column units ++++++++++++++++++*/
  /* Order is Physical rather than logical */
  /* Start File location */
  wantPos = (ObitFilePos)(unitRec-1) * 256 * sizeof(AIPSint); 

  /* Initialize Units array */
  if (desc->FieldUnit) { /* get rid of old version first */
    for (i=0; i<ncol+1; i++) 
      if (desc->FieldUnit[i]) g_free(desc->FieldUnit[i]);
    g_free(desc->FieldUnit);
  }
  desc->FieldUnit = g_malloc0((ncol+1)*sizeof(gchar*));

  /* loop reading and parsing */
  nrec = 1 + (ncol-1) / (256/2);
  icol = 0;
  for (i=0; i<nrec; i++) {

   /* read block */
    retCode = ObitFileRead(in->myFile, wantPos, size, (gchar*)record, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */
    wantPos = -1; /* continue sequential */

    /* Parse - how many? */
    ndo = 128;
    ndo = MIN (ndo, (ncol-icol)); /* Not more than number of columns */
    ip = 0;
    for (j=0; j<ndo; j++) {
     /* Inverse lookup to get logical column */
      kkol = -1;
      for (damn=0; damn<desc->nfield; damn++) 
	if (desc->order[damn]==(icol+1)) kkol = damn;
      if (desc->FieldUnit[kkol]) g_free(desc->FieldUnit[kkol]);
      desc->FieldUnit[kkol] = g_strndup((gchar*)&record[ip], 8);
      icol++;
      ip +=2; /* next unit in record */
    }
  } /* end loop reading units */

  /* Add "_status" column */
  if (desc->FieldUnit[icol]) g_free(desc->FieldUnit[icol]);
  desc->FieldUnit[icol] = g_strndup("        ", 24);

  /*++++++++++++ Read keyword/values ++++++++++++++++++*/
  /* Start File location */
  wantPos = (ObitFilePos)(keyRec-1) * 256 * sizeof(AIPSint);

  /* delete old InfoList and restart */
  ((ObitTableDesc*)in->myDesc)->info = ObitInfoListUnref (((ObitTableDesc*)in->myDesc)->info);
  ((ObitTableDesc*)in->myDesc)->info = (gpointer)newObitInfoList ();
  desc = in->myDesc; /* Table descriptor pointer */

  /* loop reading and parsing */
  nrec = 1 + (nkey-1) / (256/5);
  ikey = 0;
  for (i=0; i<nrec; i++) {

   /* read block */
    retCode = ObitFileRead(in->myFile, wantPos, size, (gchar*)record, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */
    wantPos = -1; /* continue sequential */

    /* Parse - how many? */
    ndo = 256/5;
    ndo = MIN (ndo, (nkey-ikey)); /* Not more than number of columns */
    ip = 0;
    for (j=0; j<ndo; j++) {
      g_memmove (&keyName[0], (gchar*)&record[ip],   4);
      g_memmove (&keyName[4], (gchar*)&record[ip+1], 4); keyName[8] = 0;
       /* Special trap for dates - AIPS cannot store the 10 char form 
	 yyyy-mm-dd ,  AIPS is YYYYMMDD */
      if (!strncmp(keyName, "RDATE", 5)) {
	/* AN  and others table RDATE */
	ctemp = (gchar*)&record[ip+2];
	temp[0] = ctemp[0]; temp[1] = ctemp[1]; temp[2] = ctemp[2]; temp[3] = ctemp[3];
	temp[4]='-'; temp[5] = ctemp[4]; temp[6] = ctemp[5]; temp[7] = '-';
	temp[8] = ctemp[6]; temp[9] = ctemp[7];
	g_memmove (blob, temp, 10); blob[10] = 0;
	cdim = 10;
      } else {  /* nothing special */
	/* Save 8 bytes of data */
	cdim = 8;
	g_memmove (blob, (gchar*)&record[ip+2], 8); blob[8] = 0;
      }
      /* type as ObitInfoType */
      keyType = OBIT_oint;
      if (record[ip+4]==1) keyType = OBIT_double;
      else if (record[ip+4]==2) keyType = OBIT_float;
      else if (record[ip+4]==3) keyType = OBIT_string;
      else if (record[ip+4]==4) keyType = OBIT_oint;
      else if (record[ip+4]==5) keyType = OBIT_bool;
 
     /* Save on ObitInfoList */
      dim[0] = 1;
      if (keyType == OBIT_string) dim[0] = cdim;
      ObitInfoListPut(desc->info, keyName, keyType, dim, 
		      (gconstpointer)blob, err);
      if (err->error)  /* add trace and return on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      desc->nkey = desc->info->number;

      ip +=5; /* next unit in record */
    }
  } /* end loop reading keyword/values */

  return OBIT_IO_OK;
} /* end ObitIOTableAIPSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in Pointer to object with ObitTableDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSWriteDescriptor (ObitIOTableAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableDesc* desc;
  ObitTableSel* sel;
  gchar keyName[22], *keyNameP, blob[256], *ctemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  AIPSint controlBlock[256], record[256];
  olong i, j, k, ii, ip, ncol, icol, ndo, maxkey, ikey;
  oint oitemp[4]={0,0,0,0};
  olong damn, kkol;
  olong  titleRec, unitRec, keyRec, nrec;
  ObitInfoType keyType;
  ObitFilePos wantPos;
  gsize size;
  ObitFile *readFile;
  gchar FieldName[25], temp[9];
  gboolean init, doFill;
  gchar *routine = "ObitIOTableAIPSWriteDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) return OBIT_IO_OK;

  /* enforce descriptor defaults */
  desc = in->myDesc; /* Table descriptor pointer */
  sel = in->mySel;   /* Table selector pointer */
  ObitTableSelDefault(desc, sel);
  size = 256 * sizeof(AIPSint);

  /* if it exists read old and update, ignore if WriteOnly */
  if ((in->myFile->exist) && (in->access!=OBIT_IO_WriteOnly)) {

    /* read Control block */
    readFile = newObitFile("Read Table header");
    ObitFileOpen (readFile, in->AIPSFileName, OBIT_IO_ReadOnly, 
		  OBIT_IO_Binary, 0,err);
    if (err->error)
      Obit_traceback_val (err, routine, in->name, retCode);

    wantPos = 0; /* File location */
    retCode = ObitFileRead(readFile, wantPos, size, 
			   (gchar*)controlBlock, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    /*in->filePos = in->myFile->filePos;  remember current file position */

    /* Check Recognition string */
    if (strncmp("*AIPS TABLE*", (gchar*)&controlBlock[53], 12)) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR file NOT AIPS Table for %s", in->name);
      return OBIT_IO_ReadErr;
    }
  
    /* read funky info block */
    retCode = ObitFileRead(readFile, -1L, size, (gchar*)record, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    /* in->filePos = in->myFile->filePos; remember current file position */

    /* close up */
    retCode = ObitFileClose(readFile, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);

    readFile = ObitFileUnref(readFile); /* delete */

  } /* end of read old if it exists */

  /* Set control block to internal structure */
  retCode = OBIT_IO_ReadErr;
  init = (!in->myFile->exist) || (in->access==OBIT_IO_WriteOnly);
  ObitAIPSCatTableSetDesc (in->myDesc, init,
			   in->tabType, in->tabVer,
			   controlBlock, record, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Convert Control block to internal structure to be sure */
  ObitAIPSCatTableGetDesc (in->myDesc, in->tabType, in->tabVer,
			     controlBlock, record, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Get critical control info */
  ncol     = controlBlock[9];  /* number of data columns in table */
  maxkey   = controlBlock[51]; /* Number of keyword/value pairs */
  titleRec = controlBlock[46]; /* First record (1024 byte) for col. titles */
  unitRec  = controlBlock[47]; /* First record for col. units */
  keyRec   = controlBlock[48]; /* First record for keyword/value pairs */
  
  /* If the old file existed only the Control Record and the keyword 
     sections need be written - they should be done under any 
     circumstances */

  /* Write Control block */
  wantPos = 0; /* beginning of file */
  retCode = ObitFileWrite(in->myFile, wantPos, size, 
			 (gchar*)controlBlock, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  in->filePos = in->myFile->filePos; /* remember current file position */
  wantPos = -1L; /* sequential from here */

  /* write pointer/type block if new file */
  if ((!in->myFile->exist) || (in->access==OBIT_IO_WriteOnly)) {
    retCode = ObitFileWrite(in->myFile, wantPos, size, (gchar*)record, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name,  retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */
    wantPos = -1L; /* sequential from here */
  }

  /*++++++++++++ Write keyword/values ++++++++++++++++++*/
  /* Make sure they will fit */
  if (desc->info->number > maxkey) {
    Obit_log_error(err, OBIT_Error, 
      "%s: too many keywords %d %d for %s", 
       routine, desc->info->number, maxkey, in->name);
    return OBIT_IO_WriteErr;
 }

  /* Write keywords */
  wantPos = (ObitFilePos)(keyRec-1) * 256 * sizeof(AIPSint);

  /* loop writing */
  nrec = 1 + (desc->info->number-1) / (256/5);
  ikey = 0;
  for (i=0; i<nrec; i++) {
    /* Parse - how many? */
    ndo = 256/5;
    /* Not more than number of keywords */
    ndo = MIN (ndo, (desc->info->number-ikey)); 
    ip = 0;
    for (j=0; j<ndo; j++) {
     /* Read from ObitInfoList */
      ikey++;
      ObitInfoListGetNumber(desc->info, ikey, &keyNameP, &keyType, dim, 
		      blob, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      /* Copy, possibly truncating name */
      strncpy (keyName, keyNameP, 21); keyName[21] = 0;
      /* blankfill after NULL */
      doFill = FALSE;
      for (k=0; k<8; k++) {
	if (keyName[k]==0) doFill = TRUE;
	if (doFill) keyName[k] = ' ';
      }
      /* Replace any non printing characters with blanks */
      for (k=0; k<8; k++) if (!g_ascii_isprint(keyName[k])) keyName[k]=' ';
      /* Copy to record */
      g_memmove ((gchar*)&record[ip],  &keyName[0], 4);
      g_memmove ((gchar*)&record[ip+1],&keyName[4], 4); 
      /* Special trap for dates - AIPS cannot store the 10 char form 
	 yyyy-mm-dd ,  AIPS is YYYYMMDD */
      if (!strncmp(keyName, "RDATE", 5)) {
	/* AN and others table RDATE */
	temp[0] = blob[0]; temp[1] = blob[1]; temp[2] = blob[2]; temp[3] = blob[3];
	temp[4] = blob[5]; temp[5] = blob[6]; temp[6] = blob[8]; temp[7] = blob[9];
	g_memmove ((gchar*)&record[ip+2], temp, 8); 
      } else {  /* nothing special */
	/* Save 8 bytes of data */
	g_memmove ((gchar*)&record[ip+2], blob, 8); 
      }
      /* Convert type to AIPSish */
      record[ip+4] = 4; /* default int */
      if (keyType==OBIT_double)      record[ip+4] = 1;
      else if (keyType==OBIT_float)  record[ip+4] = 2;
      else if (keyType==OBIT_string) {
	/* Replace any non printing characters with blanks */
	ctemp = (gchar*)&record[ip+2];
	for (k=0; k<8; k++) if (!g_ascii_isprint(ctemp[k])) ctemp[k]=' ';
	record[ip+4] = 3;
      }
      else if (keyType==OBIT_long) { /* May have to convert long->oint */
	record[ip+4] = 4;
	oitemp[0] = (oint)*(olong*)blob;
	g_memmove ((gchar*)&record[ip+2], oitemp, 8); blob[8] = 0;
      }
      else if (keyType==OBIT_long)   record[ip+4] = 4;
      else if (keyType==OBIT_long)   record[ip+4] = 4;
      else if (keyType==OBIT_oint)   record[ip+4] = 4;
      else if (keyType==OBIT_bool)   record[ip+4] = 5;
      ip +=5; /* next unit in record */
    } /* end loop filling block */

    /* write block */
    retCode = ObitFileWrite(in->myFile, wantPos, size, (gchar*)record, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
    in->filePos = in->myFile->filePos; /* remember current file position */
    wantPos = -1L; /* sequential from here */
  } /* end loop reading keyword/values */

  /* If the old version didn't exist, the selection strings, titles and
     Units sections must be written */
  if ((!in->myFile->exist) || (in->access==OBIT_IO_WriteOnly)){

  /*++++++++++++ Write Selection Strings ++++++++++++++++++*/
    /* Selection strings - blank fill two records */
    for (i=0; i<size; i++) ((gchar*)record)[i] = ' ';
    wantPos = 2 * 256 * sizeof(AIPSint); /* File location */

    for (i=0; i<2; i++) {
      retCode = ObitFileWrite(in->myFile, wantPos, size, (gchar*)record, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      in->filePos = in->myFile->filePos; /* remember current file position */
      wantPos = -1L; /* sequential from here */
    } /* end of blanking selection strings */

    /*++++++++++++ Write column titles ++++++++++++*/
    wantPos = (ObitFilePos)(titleRec-1) * 256 * sizeof(AIPSint);
    /* loop filling a block and writing */
    nrec = 1 + (ncol-1) / (256/6);
    icol = 0;
    for (i=0; i<nrec; i++) {
      /* copy - how many? */
      ndo = 42;
      ndo = MIN (ndo, (ncol-icol)); /* Not more than number of columns */
      ip = 0;
      for (j=0; j<ndo; j++) {
	/* Inverse lookup to get logical column */
	kkol = -1;
	for (damn=0; damn<desc->nfield; damn++) 
	  if (desc->order[damn]==(icol+1)) kkol = damn;
	/* Pad field name to 24 char */
	for (k=0; k<24; k++) FieldName[k] = ' ';
	for (k=0; k<MIN (24, strlen(desc->FieldName[kkol])); k++) 
	  FieldName[k] = desc->FieldName[kkol][k];
	g_memmove ((gchar*)&record[ip], FieldName, 24);
	icol++;
	ip +=6; /* next title in record */
      }

      /* write block */
      retCode = ObitFileWrite(in->myFile, wantPos, size, (gchar*)record, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      in->filePos = in->myFile->filePos; /* remember current file position */
      wantPos = -1L; /* sequential from here */
    } /* end loop writing titles */
    
    /*++++++++++++ Write column units ++++++++++++++++++*/
    wantPos = (ObitFilePos)(unitRec-1) * 256 * sizeof(AIPSint);
    /* loop over blocks writing */
    nrec = 1 + (ncol-1) / (256/2);
    icol = 0;
    for (i=0; i<nrec; i++) {
      /* Copy - how many? */
      ndo = 128;
      ndo = MIN (ndo, (ncol-icol)); /* Not more than number of columns */
      ip = 0;
      for (j=0; j<ndo; j++) {
      /* Inverse lookup to get logical column */
	kkol = -1;
	for (damn=0; damn<desc->nfield; damn++) 
	  if (desc->order[damn]==(icol+1)) kkol = damn;
	/* blank fill units */
	for (ii=0; ii<8; ii++) temp[i] = ' ';
	strncpy(temp, desc->FieldUnit[kkol], 8);
	for (ii=0; ii<8; ii++) if (temp[i]==0) temp[i] = ' ';
	g_memmove ((gchar*)&record[ip], temp, 8);
	icol++;
	ip +=2; /* next unit in record */
      }

      /* write block */
      retCode = ObitFileWrite(in->myFile, wantPos, size, (gchar*)record, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
	Obit_traceback_val (err, routine, in->name, retCode);
      in->filePos = in->myFile->filePos; /* remember current file position */
      wantPos = -1L; /* sequential from here */
    } /* end loop writing units */
    
  } /* end of file didn't previously exist section */

  /* Flush buffer to be sure everything is actually written */
  retCode = ObitFileFlush(in->myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);
  

  return OBIT_IO_OK;
} /* end ObitIOTableAIPSWriteDescriptor */

/**
 * Flush I/O buffer if necessary, padding is added to keep
 * Mama AIPS happy.
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, 0=> OK
 */
ObitIOCode ObitIOTableAIPSFlush (ObitIOTableAIPS *in, ObitErr *err)
{
  olong blksize;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOTableAIPSFlush";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  /* Ignore if "special" AIPS extension */
  if (!ObitIOTableAIPSIsTable(in)) return OBIT_IO_OK;

  /* Flush buffer to be sure everything is actually written */
  retCode = ObitFileFlush(in->myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Curse you Mama AIPS!!!!!*/
  /* Because of the archaic AIPS file structure, the file must be 
     an exact multiple of 256*sizeof(AIPSint) in size.  On the other
     hand, you don't want to write over existing data so you can't
     just pad out the current block.  Pad the end of the file. */
  blksize = 256*sizeof(AIPSint);
  if ((in->access==OBIT_IO_ReadWrite) || 
      (in->access==OBIT_IO_WriteOnly)) 
    retCode = ObitFilePadFile(in->myFile, blksize, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Flush buffer (again) to be sure everything is actually written */
  retCode = ObitFileFlush(in->myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitIOTableAIPSFlush */

/**
 * Create buffer appropriate for I/O request
 * \param data (output) pointer to data array
 * \param size (output) size of data array in bytes.
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 */
void 
ObitIOTableAIPSCreateBuffer (ofloat **data, olong *size, 
			     ObitIOTableAIPS *in, ObitInfoList *info, 
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

} /* end ObitIOTableAIPSCreateBuffer */

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
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFileName OBIT_string  FITS file name
 * \li xxxDisk     OBIT_oint    FITS file disk number
 * \li xxxDir      OBIT_string  Directory name for xxxDisk
 *
 * For all File types types:
 * \li xxxDataType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxFileType OBIT_oint File type as "AIPS
 *    
 * For xxxDataType = "Table"
 * \li xxxTableParent OBIT_string  Table parent type (e.g. "MA", "UV")
 * \li xxxType  OBIT_string  (Tables only) Table type (e.g. "AIPS CC")
 * \li xxxVer   OBIT_oint    (Tables Only) Table version number
 *    
 * \param err     ObitErr for reporting errors.
 */
void ObitIOTableAIPSGetFileInfo (ObitIO *in, ObitInfoList *myInfo, gchar *prefix, 
				 ObitInfoList *outList, ObitErr *err)
{
  ObitInfoType type;
  ObitAIPSDirCatEntry *entry=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *FileType="AIPS", *dirname;
  gchar *DataType[]={"Table","UV","MA"};
  gchar tempStr[201];
  olong disk, user, cno, seq, ver, ptype=1;
  gchar *routine = "ObitIOTableAIPSGetFileInfo";

  if (err->error) return;

  /* Set basic information */
  /* AIPS */
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
  dirname =  ObitAIPSDirname (disk, err);
  dim[0] = strlen(dirname);
  if (prefix) keyword =  g_strconcat (prefix, "Dir", NULL);
  else keyword =  g_strdup ("Dir");
  ObitInfoListAlwaysPut (outList, keyword, OBIT_string, dim, dirname);
  g_free(keyword);

} /* end ObitIOTableAIPSGetFileInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOTableAIPSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOTableAIPSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOTableAIPSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOTableAIPSClassInfoDefFn (gpointer inClass)
{
  ObitIOTableAIPSClassInfo *theClass = (ObitIOTableAIPSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOTableAIPSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOTableAIPSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOTableAIPSGetClass;
  theClass->newObit       = NULL;
  theClass->newObitIO     = (newObitIOFP)newObitIOTableAIPS;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOTableAIPSSame;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOTableAIPSZap;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOTableAIPSCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = 
    (ObitClearFP)ObitIOTableAIPSClear;
  theClass->ObitInit      = 
    (ObitInitFP)ObitIOTableAIPSInit;
  theClass->ObitIOOpen    = 
    (ObitIOOpenFP)ObitIOTableAIPSOpen;
  theClass->ObitIOClose   = 
    (ObitIOCloseFP)ObitIOTableAIPSClose;
  theClass->ObitIOSet     = 
    (ObitIOSetFP)ObitIOTableAIPSSet;
  theClass->ObitIORead    = 
    (ObitIOReadFP)ObitIOTableAIPSRead;
  theClass->ObitIOReadRow = 
    (ObitIOReadRowFP)ObitIOTableAIPSReadRow;
  theClass->ObitIOReadRowSelect = 
    (ObitIOReadRowSelectFP)ObitIOTableAIPSReadRowSelect;
  theClass->ObitIOReadSelect = 
    (ObitIOReadSelectFP)ObitIOTableAIPSReadSelect;
  theClass->ObitIOWrite   = 
    (ObitIOWriteFP)ObitIOTableAIPSWrite;
  theClass->ObitIOWriteRow   = 
    (ObitIOWriteRowFP)ObitIOTableAIPSWriteRow;
  theClass->ObitIOFlush   = 
    (ObitIOFlushFP)ObitIOTableAIPSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOTableAIPSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOTableAIPSWriteDescriptor;
  theClass->ObitIOCreateBuffer = 
    (ObitIOCreateBufferFP)ObitIOTableAIPSCreateBuffer;
  theClass->ObitIOFreeBuffer   = 
    (ObitIOFreeBufferFP)ObitIOFreeBuffer;
  theClass->ObitIOGetFileInfo   =
    (ObitIOGetFileInfoFP)ObitIOTableAIPSGetFileInfo;

} /* end ObitIOTableAIPSClassDefFn */

/*--------------- Private functions --------------------------*/

/**
 * Creates empty member objects.
 * for each parent class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOTableAIPSInit  (gpointer inn)
{
  const ObitClassInfo *ParentClass;
  ObitIOTableAIPS *in = inn;

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
} /* end ObitIOTableAIPSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOTableAIPSClear (gpointer inn)
{
  ObitIOTableAIPS *in = inn;
  const ObitClassInfo *ParentClass;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* close I/O if still active */
  if ((in->myStatus==OBIT_Active) ||(in->myStatus==OBIT_Modified)) {
    err = newObitErr();
    ObitIOTableAIPSClose (in, err); 
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

} /* end ObitIOTableAIPSClear */

/**
 * Uses selector member to decide which visibilities to
 * read next.
 * Leaves values in myDesc as firstVis and numVisBuff.
 * \param  in Pointer to the object.
 * \param rowno Starting row number (1-rel), -1=> next;
 * \return TRUE is finished, else FALSE
 */
static gboolean ObitIOTableAIPSNext (ObitIOTableAIPS *in, olong rowno)
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
} /* end ObitIOTableAIPSNext */

/** Private: Is this an AIPS table? */
static gboolean ObitIOTableAIPSIsTable (ObitIOTableAIPS *in)
{
  gboolean isTable = TRUE;
  olong i;

  i = 0;
  while (NotAIPSTable[i]) {
    if (!strncmp (NotAIPSTable[i],in->tabType, 2))
      isTable = FALSE;
    i++;
  }
  return isTable;
} /* end ObitIOTableAIPSIsTable */
