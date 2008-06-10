/* $Id: ObitIOHistoryFITS.c,v 1.9 2007/08/23 14:50:48 bcotton Exp $ */
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

#include "fitsio.h"
#include "ObitIOHistoryFITS.h"
#include "ObitTableHistory.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOHistoryFITS.c
 * ObitIOHistoryFITS class function definitions.
 *
 * This class allows history access to FITS files - via HISTORY table
 * This class is derived from ObitIOHistory
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOHistoryFITS";

/** Function to obtain parent ClassInfo - ObitIOHistory */
static ObitGetClassFP ObitParentGetClass = ObitIOHistoryGetClass;

/**
 * ClassInfo global structure ObitIOHistoryFITSClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOHistoryFITSClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOHistoryFITSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOHistoryFITSClear (gpointer in);

/** Private: Create ObitHistory table descriptor. */
static void  CreateHiDescriptor (ObitTableDesc *desc);

/** Private: Set Class function pointers. */
static void ObitIOHistoryFITSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name Name [optional] for object
 * \param info InfoList defining file
 * \param err ObitErr for reporting errors. 
 * \return the new object.
 */
ObitIOHistoryFITS* newObitIOHistoryFITS (gchar* name, ObitInfoList *info,
				 ObitErr *err)
{
  ObitIOHistoryFITS* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOHistoryFITSClassInit();

 /* allocate/init structure */
  out = g_malloc0(sizeof(ObitIOHistoryFITS));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitIOHistoryFITSInit((gpointer)out);

  /* Copy ObitInfoList to get relevant information */
  out->info = ObitInfoListUnref(out->info); /* out with the old - empty */
  out->info = ObitInfoListCopy(info);       /* in with the new - copy */

  return out;
} /* end newObitIOHistoryFITS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOHistoryFITSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOHistoryFITSClassInit();
  
  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

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
gboolean ObitIOHistoryFITSSame (ObitIO *in, ObitInfoList *in1, 
				ObitInfoList *in2, ObitErr *err)
{
  olong disk1, disk2;
  gchar *filename1, *filename2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOHistoryFITSSame";

  /* error checks */
  if (err->error) return same;

  /* get file from info */
  if(!ObitInfoListGet(in1, "Disk", &type, dim, &disk1, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in1, "FileName", &type, dim, 
		       (gpointer)&filename1)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList Object %s",
		   routine, in->name);
  }

  if(!ObitInfoListGet(in2, "Disk", &type, dim, &disk2, err))
    Obit_traceback_val (err, routine, in->name, same);

  if (!ObitInfoListGetP(in2, "FileName", &type, dim, 
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
} /* end ObitIOHistoryFITSSame */

/**
 * Delete underlying files.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors. 
 */
void ObitIOHistoryFITSZap (ObitIOHistoryFITS *in, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *routine = "ObitIOHistoryFITSZap";

  /* Don't bother if NULL */
  if (!in) return;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Open and close if table not yet instantiated */
  if (in->table) {
    retCode = ObitIOHistoryFITSOpen (in, OBIT_IO_ReadWrite, in->info, err);
    retCode = ObitIOHistoryFITSClose(in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Zap it */
  in->table = (ObitTableHistory*)ObitTableZap ((ObitTable*)in->table, err);

  /* Don't need row either */
  in->row = ObitUnref(in->row);

  return;
} /* end ObitIOHistoryFITSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object. 
 * \return pointer to the new object.
 */
ObitIOHistoryFITS* ObitIOHistoryFITSCopy  (ObitIOHistoryFITS *in, ObitIOHistoryFITS *out, ObitErr *err)
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
    out = newObitIOHistoryFITS(outName, NULL, err);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->table  = ObitUnref(out->table);
  out->table  = ObitRef(in->table);
  out->row    = ObitUnref(out->row);
  out->row    = ObitRef(in->row);    /* seems risky */

  return out;
} /* end ObitIOHistoryFITSCopy */

/**
 * Initialize structures and open file.
 * The file and selection info member should have been stored in the ObitInfoList
 * prior to calling.  See derived classes for details.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite)
 * \param info ObitInfoList with instructions for opening
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSOpen (ObitIOHistoryFITS *in, ObitIOAccess access, 
				  ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar FileName[201], *tab="History";
  olong disk, ver = 1, nrow = 1;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "ObitIOHistoryFITSOpen";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));

  /* Create table if it doesn't exist */
  if (in->table==NULL) {
    
    /* get FITS file info */
    if(!ObitInfoListGet(in->info, "Disk", &type, dim, &disk, err))
      Obit_traceback_val (err, routine, in->name, OBIT_IO_ReadErr);
    
    if (!ObitInfoListGet(in->info, "FileName", &type, dim, FileName, err))
      Obit_traceback_val (err, routine, in->name, OBIT_IO_ReadErr);
    FileName[dim[0]] = 0;

    /* Create */
    in->table = newObitTableHistory("History");
    ObitTableSetFITS(in->table,disk,FileName,tab,ver,nrow,err);
    if (err->error) Obit_traceback_val (err, routine, in->name, OBIT_IO_ReadErr);

    /* Create descriptor */
    CreateHiDescriptor(in->table->myDesc);
  }

  /* Open table */
  retCode = ObitTableHistoryOpen (in->table, access, err);
  if (err->error) {
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Create row if it doesn't exist */
  if (in->row==NULL) {
    in->row = newObitTableHistoryRow(in->table);
  }

   /* Keep track of where and how many */
  in->MaxNumber     = ((ObitTableDesc*)in->table->myIO->myDesc)->nrow;
  in->CurrentNumber = in->MaxNumber;
  in->Last          = 0;

  /* If writing attach row to table buffer */
  if ((access==OBIT_IO_ReadWrite) || (access==OBIT_IO_WriteOnly)) {
    ObitTableHistorySetRow(in->table, in->row, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, OBIT_IO_ReadErr);
    in->Last = in->MaxNumber;
  }

  in->access = access; /* just in case */

  return retCode;
} /* end ObitIOHistoryFITSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSClose (ObitIOHistoryFITS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOHistoryFITSClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Close table */
  retCode = ObitTableHistoryClose (in->table, err);
  if (err->error) {
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    Obit_traceback_val (err, routine, in->name, retCode);
  }

  return retCode;
} /* end ObitIOHistoryFITSClose */

/**
 * initialize I/O
 * Not needed for FITS
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSSet (ObitIOHistoryFITS *in, ObitInfoList *info, ObitErr *err)
{

  /* not needed */
  return OBIT_IO_OK;
} /* end ObitIOHistoryFITSSet */

/**
 * Read specified History record
 * \param in     Pointer to object to be read.
 * \param recno  record number (1-rel) -1=> next.
 * \param hiCard output history record (70 char)
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSReadRec (ObitIOHistoryFITS *in, olong recno, 
				     gchar *hiCard, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong irow;
  gchar *routine = "ObitIOHistoryFITSReadRec";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (hiCard != NULL);

  /* which record? */
  irow = recno;
  if (irow<0) irow = in->Last+1;

  /* read */
  retCode = ObitTableHistoryReadRow (in->table, irow, in->row, err);
  if (err->error) {
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Copy to output */
  g_snprintf (hiCard, 70, "%s", in->row->entry);

  /* Keep track */
  in->Last = irow;

  return retCode;
} /* end ObitIOHistoryFITSReadRec */

/**
 * Write specified History record
 * \param in     Pointer to object to be written.
 * \param recno  Record number (1-rel) -1=> next, overwrites any existing
 * \param hiCard input history record (70 char)
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSWriteRec (ObitIOHistoryFITS *in, olong recno, 
				      gchar *hiCard, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong irow;
  gchar *ctemp;
  gchar *routine = "ObitIOHistoryFITSWriteRec";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (hiCard != NULL);

  /* which record? Default is at end */
  irow = recno;
  if (irow<0) irow = in->CurrentNumber+1;

  /* Write */
  /*g_snprintf (in->row->entry, 70, "%s", hiCard);*/
  ctemp = in->row->entry;
  in->row->entry = hiCard;
  retCode = ObitTableHistoryWriteRow (in->table, irow, in->row, err);
  in->row->entry = ctemp;
  if (err->error) {
    Obit_cfitsio_error(err); /* copy cfitsio error stack */
    Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Keep track */
  in->Last = irow;
  in->CurrentNumber = MAX (in->CurrentNumber, irow);

  return retCode;
} /* end ObitIOHistoryFITSWriteRec */

/**
 * Tell number of History records
 * \param in     Pointer to open object to be tested
 * \return  number of records, <0 => problem
 */
olong ObitIOHistoryFITSNumRec (ObitIOHistory *in)
{
  olong out;

  out = ((ObitIOHistoryFITS*)in)->table->myDesc->nrow;

  return out;
} /* end ObitIOHistoryFITSNumRec */

/**
 * Read image Descriptor data from disk.
 * The Table system does this automagically for FITS files
 * \param in Pointer to object  be read.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSReadDescriptor (ObitIOHistoryFITS *in,  ObitErr *err)
{
  /* not needed */
  return OBIT_IO_OK;

} /* end ObitIOHistoryFITSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * The Table system does this automagically for FITS files
 * \param in Pointer to object  to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSWriteDescriptor (ObitIOHistoryFITS *in, ObitErr *err)
{
  /* not needed */
  return OBIT_IO_OK;

} /* end ObitIOHistoryFITSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * This is handled by cfitsio
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFITSFlush (ObitIOHistoryFITS *in, ObitErr *err)
{

  /* not needed */
  return OBIT_IO_OK;

} /* end ObitIOHistoryFITSFlush */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOHistoryFITSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOHistoryFITSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOHistoryFITSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOHistoryFITSClassInfoDefFn (gpointer inClass)
{
  ObitIOHistoryFITSClassInfo *theClass = (ObitIOHistoryFITSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOHistoryFITSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOHistoryFITSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOHistoryFITSGetClass;
  theClass->newObit       = NULL;
  theClass->newObitIO     = (newObitIOFP)newObitIOHistoryFITS;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOHistoryFITSCopy;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOHistoryFITSSame;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOHistoryFITSZap;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitIOHistoryFITSClear;
  theClass->ObitInit      = (ObitInitFP)ObitIOHistoryFITSInit;
  theClass->ObitIOOpen    = (ObitIOOpenFP)ObitIOHistoryFITSOpen;
  theClass->ObitIOClose   = (ObitIOCloseFP)ObitIOHistoryFITSClose;
  theClass->ObitIOSet     = (ObitIOSetFP)ObitIOHistoryFITSSet;
  theClass->ObitIORead          = NULL;
  theClass->ObitIOReadSelect    = NULL;
  theClass->ObitIOReadRowSelect = NULL;
  theClass->ObitIOWriteRow      = NULL;
  theClass->ObitIOWrite         = NULL;
  theClass->ObitIOFlush         = (ObitIOFlushFP)ObitIOHistoryFITSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOHistoryFITSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOHistoryFITSWriteDescriptor;
  theClass->ObitIOCreateBuffer = NULL;
  theClass->ObitIOFreeBuffer   = NULL;
  theClass->newObitIOTable     = NULL;
  theClass->ObitIOUpdateTables = NULL;
  theClass->ObitIOHistoryReadRec  = 
    (ObitIOHistoryReadRecFP)ObitIOHistoryFITSReadRec;
  theClass->ObitIOHistoryWriteRec = 
    (ObitIOHistoryWriteRecFP)ObitIOHistoryFITSWriteRec;
  theClass->ObitIOHistoryNumRec  = 
    (ObitIOHistoryNumRecFP)ObitIOHistoryFITSNumRec;

} /* end ObitIOHistoryFITSClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOHistoryFITSInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIOHistoryFITS *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->row   = NULL;
  in->table = NULL;
} /* end ObitIOHistoryFITSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOHistoryFITSClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIOHistoryFITS *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  in->row   = ObitUnref(in->row);
  in->table = ObitUnref(in->table);
  
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOHistoryFITSClear */

/** Private: Create ObitHistory table descriptor. */
static void  CreateHiDescriptor (ObitTableDesc *desc)
{
  olong colNo, i, ncol;

  /* How many columns actually in table? */
  ncol = 1 ;
  desc->FieldName = g_malloc0((ncol+1)*sizeof(gchar*));
  desc->FieldUnit = g_malloc0((ncol+1)*sizeof(gchar*));
  desc->type      = g_malloc0((ncol+1)*sizeof(ObitInfoType));
  desc->dim       = g_malloc0((ncol+1)*sizeof(gint32*));
  for (i=0; i<ncol+1; i++) 
    desc->dim[i] = g_malloc0(MAXINFOELEMDIM*sizeof(gint32));

  desc->TableName = g_strdup("History");
  desc->sort[0] = 0;
  desc->sort[1] = 0;
  colNo = 0;

  /* Define Columns */
  desc->FieldName[colNo] = g_strdup("ENTRY  ");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_string;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  desc->dim[colNo][0] = 70;
  colNo++;
  /* Add _status column at end */
  desc->FieldName[colNo] = g_strdup("_status");
  desc->FieldUnit[colNo] = g_strdup("        ");
  desc->type[colNo] = OBIT_long;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  
  /* number of fields */
  desc->nfield = colNo + 1;

  /* index table descriptor */
  ObitTableDescIndex (desc);

}
