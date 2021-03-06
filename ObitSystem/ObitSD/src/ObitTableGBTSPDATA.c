/* $Id$  */
/* DO NOT EDIT - file generated by ObitSDTables.pl                    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2009                                              */
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
/*;         Internet email: bcotton@nrao.edu.                        */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "ObitTableGBTSPDATA.h"
#include "ObitTableList.h"
#include "ObitData.h"

/*----------------Obit: Merx mollis mortibus nuperS ------------------*/
/**
 * \file ObitTableGBTSPDATA.c
 * ObitTableGBTSPDATA class function definitions.
 * This class is derived from the ObitTable class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitTableGBTSPDATA";

/**  Function to obtain parent Table ClassInfo - ObitTable */
static ObitGetClassFP ObitParentGetClass = ObitTableGetClass;

/** name of the Row class defined in this file */
static gchar *myRowClassName = "ObitTableGBTSPDATARow";

/**  Function to obtain parent TableRow ClassInfo */
static ObitGetClassFP ObitParentGetRowClass = ObitTableRowGetClass;

/*--------------- File Global Variables  ----------------*/
/*----------------  Table Row  ----------------------*/
/**
 * ClassInfo structure ObitTableClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableGBTSPDATARowClassInfo myRowClassInfo = {FALSE};

/*------------------  Table  ------------------------*/
/**
 * ClassInfo structure ObitTableGBTSPDATAClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableGBTSPDATAClassInfo myClassInfo = {FALSE};


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated Row object. */
void  ObitTableGBTSPDATARowInit  (gpointer in);

/** Private: Deallocate Row members. */
void  ObitTableGBTSPDATARowClear (gpointer in);

/** Private: Initialize newly instantiated object. */
void  ObitTableGBTSPDATAInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTableGBTSPDATAClear (gpointer in);

/** Private: update table specific info */
static void ObitTableGBTSPDATAUpdate (ObitTableGBTSPDATA *in, ObitErr *err);

/** Private: copy table keywords to descriptor info list */
static void ObitTableGBTSPDATADumpKey (ObitTableGBTSPDATA *in, ObitErr *err);

/** Private: Set Class function pointers */
static void ObitTableGBTSPDATAClassInfoDefFn (gpointer inClass);

/** Private: Set Row Class function pointers */
static void ObitTableGBTSPDATARowClassInfoDefFn (gpointer inClass);
/*----------------------Public functions---------------------------*/

/*------------------  Table Row ------------------------*/
/**
 * Constructor.
 * If table is open and for write, the row is attached to the buffer
 * Initializes Row class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTableGBTSPDATARow* newObitTableGBTSPDATARow (ObitTableGBTSPDATA *table)
{
  ObitTableGBTSPDATARow* out;
  odouble   *dRow;
  oint      *iRow;
  gshort    *siRow;
  ofloat    *fRow;
  gchar     *cRow;
  gboolean  *lRow;
  guint8    *bRow;

  /* Class initialization if needed */
  if (!myRowClassInfo.initialized) ObitTableGBTSPDATARowClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTableGBTSPDATARow));

  /* initialize values */
  out->name = g_strdup("TableObitTableGBTSPDATA Row");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myRowClassInfo;

  /* initialize other stuff */
  ObitTableGBTSPDATARowInit((gpointer)out);
  out->myTable   = (ObitTable*)ObitTableRef((ObitTable*)table);

  /* If writing attach to buffer */
  if ((table->buffer) && (table->myDesc->access != OBIT_IO_ReadOnly) &&
      (table->myStatus != OBIT_Inactive)) {
    /* Typed pointers to row of data */  
    dRow  = (odouble*)table->buffer;
    iRow  = (oint*)table->buffer;
    siRow = (gshort*)table->buffer;
    fRow  = (ofloat*)table->buffer;
    cRow  = (gchar*)table->buffer;
    lRow  = (gboolean*)table->buffer;
    bRow  = (guint8*)table->buffer;
  
    /* Set row pointers to buffer */
    out->data = fRow + table->dataOff;
  } /* end attaching row to table buffer */

 return out;
} /* end newObitTableDATARow */

/**
 * Returns ClassInfo pointer for the Row class.
 * \return pointer to the Row class structure.
 */
gconstpointer ObitTableGBTSPDATARowGetClass (void)
{
  /* Class initialization if needed */
  if (!myRowClassInfo.initialized) ObitTableGBTSPDATARowClassInit();
  return (gconstpointer)&myRowClassInfo;
} /* end ObitTableGBTSPDATARowGetClass */

/*------------------  Table  ------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTableGBTSPDATA* newObitTableGBTSPDATA (gchar* name)
{
  ObitTableGBTSPDATA* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableGBTSPDATAClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTableGBTSPDATA));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTableGBTSPDATAInit((gpointer)out);

 return out;
} /* end newObitTableGBTSPDATA */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitTableGBTSPDATAGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableGBTSPDATAClassInit();
  return (gconstpointer)&myClassInfo;
} /* end ObitTableGBTSPDATAGetClass */

/**
 * Constructor from values.
 * Creates a new table structure and attaches to the TableList of file.
 * If the specified table already exists then it is returned.
 * Initializes class if needed on first call.
 * Forces an update of any disk resident structures (e.g. AIPS header).
 * \param name   An optional name for the object.
 * \param file   ObitData which which the table is to be associated.
 * \param ver    Table version number. 0=> add higher, value used returned
 * \param access access (OBIT_IO_ReadOnly, means do not create if it doesn't exist.
 * \param err Error stack, returns if not empty.
 * \return the new object, NULL on failure.
 */
ObitTableGBTSPDATA* newObitTableGBTSPDATAValue (gchar* name, ObitData *file, olong *ver,
 				  ObitIOAccess access,
  		    
		     ObitErr *err)
{
  ObitTableGBTSPDATA* out=NULL;
  ObitTable *testTab=NULL;
  ObitTableDesc *desc=NULL;
  ObitTableList *list=NULL;
  ObitInfoList  *info=NULL;
  gboolean exist, optional;
  olong colNo, i, ncol, highVer;
  ObitIOCode retCode;
  gchar *tabType = "DATA";
  gchar *routine = "newObitTableGBTSPDATAValue";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitDataIsA(file));

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableGBTSPDATAClassInit();

  /* Check if the table already exists */
  /* Get TableList */
  list = ((ObitData*)file)->tableList;
  info = ((ObitData*)file)->info;

  /* See if it already exists */
  exist = FALSE;
  if (*ver>0) { /* has to be specified */
    exist = ObitTableListGet(list, tabType, ver, &testTab, err);
    if (err->error) /* add traceback,return */
      Obit_traceback_val (err, routine,"", out);
  
    /* if readonly, it must exist to proceed */
    if ((access==OBIT_IO_ReadOnly) && !exist) return out;
    if (testTab!=NULL) { /* it exists, use it if is an ObitTableGBTSPDATA */
      if (ObitTableGBTSPDATAIsA(testTab)) { /* it is an ObitTableGBTSPDATA */
	out = ObitTableRef(testTab);
      } else { /* needs conversion */
 	out = ObitTableGBTSPDATAConvert(testTab);
	/* Update the TableList */
	ObitTableListPut(list, tabType, ver, (ObitTable*)out, err);
	if (err->error) Obit_traceback_val (err, routine,"", out);
      }
      testTab = ObitTableUnref(testTab); /* remove reference */
      return out; /* done */
    }
  } /* end of test for previously existing table */
  
  /* If access is ReadOnly make sure one exists */
  if (access==OBIT_IO_ReadOnly) { 
    highVer = ObitTableListGetHigh (list, "DATA");
    if (highVer<=0) return out;
  }
  
  /* create basal table */
  testTab = newObitDataTable ((ObitData*)file, access, tabType,
			       ver, err);
  if (err->error) Obit_traceback_val (err, routine,"", out);
  
  /* likely need to convert */
  if (ObitTableGBTSPDATAIsA(testTab)) { 
    out = ObitTableRef(testTab);
  } else { /* needs conversion */
    out = ObitTableGBTSPDATAConvert(testTab);
  }
  testTab = ObitTableUnref(testTab); /* remove reference */

  /* Update the TableList */
  ObitTableListPut(list, tabType, ver, (ObitTable*)out, err);
  if (err->error) /* add traceback,return */
    Obit_traceback_val (err, routine,"", out);

  /* if it previously existed merely return it */
  if (exist) return out; 

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* Set values */
  strncpy (out->ctype1, "STATE", MAXKEYCHARTABLEDATA );
  strncpy (out->ctype2, "RECEIVER", MAXKEYCHARTABLEDATA );

  /* initialize descriptor */
  desc = out->myDesc;
  /* How many columns actually in table? */
  ncol = 5 ;
  desc->FieldName = g_malloc0((ncol+1)*sizeof(gchar*));
  desc->FieldUnit = g_malloc0((ncol+1)*sizeof(gchar*));
  desc->type      = g_malloc0((ncol+1)*sizeof(ObitInfoType));
  desc->dim       = g_malloc0((ncol+1)*sizeof(gint32*));
  for (i=0; i<ncol+1; i++) 
    desc->dim[i] = g_malloc0(MAXINFOELEMDIM*sizeof(gint32));

  desc->TableName = g_strdup(tabType);
  desc->sort[0] = 0;
  desc->sort[1] = 0;
  colNo = 0;

  /* Define Columns */
  desc->FieldName[colNo] = g_strdup("SUBSCAN ");
  desc->FieldUnit[colNo] = g_strdup("CODE");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("UTDATE ");
  desc->FieldUnit[colNo] = g_strdup("DAY");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("UTCSTART ");
  desc->FieldUnit[colNo] = g_strdup("SECOND");
  desc->type[colNo] = OBIT_double;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("PSRPER ");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_double;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  optional = FALSE;
  if ((1024 > 0) || (!optional)) {
    desc->FieldName[colNo] = g_strdup("DATA ");
    desc->FieldUnit[colNo] = g_strdup("COUNT");
    desc->type[colNo] = OBIT_float;
    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
    desc->dim[colNo][0] = 1024;
    desc->dim[colNo][1] = 2;
    desc->dim[colNo][2] = 2;
    colNo++;
  }
  /* Add _status column at end */
  desc->FieldName[colNo] = g_strdup("_status");
  desc->FieldUnit[colNo] = g_strdup("        ");
  desc->type[colNo] = OBIT_long;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  
  /* number of fields */
  desc->nfield = colNo + 1;

  /* initialize descriptor keywords */
  ObitTableGBTSPDATADumpKey (out, err);
 
  /* index table descriptor */
  ObitTableDescIndex (desc);

  /* Open and Close to fully instantiate */
  retCode = ObitTableGBTSPDATAOpen(out, OBIT_IO_WriteOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out);    
  
  retCode = ObitTableGBTSPDATAClose(out, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out); 

  /* Force update of disk resident info */
  retCode = ObitIOUpdateTables (((ObitData*)file)->myIO, info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out); 
  
 return out;
} /* end newObitTableGBTSPDATAValue */

/**
 * Convert an ObitTable to an ObitTableGBTSPDATA.
 * New object will have references to members of in.
 * \param in  The object to copy, will still exist afterwards 
 *            and should be Unrefed if not needed.
 * \return pointer to the new object.
 */
ObitTableGBTSPDATA* ObitTableGBTSPDATAConvert (ObitTable* in)
{
  ObitTableGBTSPDATA *out;

  /* error check */
  g_assert(ObitTableIsA(in));

  /* create basic object */
  out = newObitTableGBTSPDATA(in->name);

  /* Delete structures on new */
  out->info   = ObitInfoListUnref(out->info);
  out->thread = ObitThreadUnref(out->thread);
  out->myDesc = ObitTableDescUnref(out->myDesc);
  out->mySel  = ObitTableSelUnref(out->mySel);
  
  /* Reference members of in */
  
  out->info   = ObitInfoListRef(in->info);
  out->thread = ObitThreadRef(in->thread);
  out->myDesc = ObitTableDescRef(in->myDesc);
  out->mySel  = ObitTableSelRef(in->mySel);

  /* Remember who I am */
 out->tabType = g_strdup(in->tabType);
 out->tabVer  = in->tabVer;
  /* Secret reference to host */ 
 out->myHost  = in->myHost;

  return out;
} /* end ObitTableGBTSPDATAConvert */


/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Objects should be closed on input and will be closed on output.
 * In order for the disk file structures to be copied, the output file
 * must be sufficiently defined that it can be written.
 * The copy will be attempted but no errors will be logged until
 * both input and output have been successfully opened.
 * ObitInfoList and ObitThread members are only copied if the output object
 * didn't previously exist.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitTableGBTSPDATA* ObitTableGBTSPDATACopy (ObitTableGBTSPDATA *in, ObitTableGBTSPDATA *out, ObitErr *err)
{
  gchar *routine = "ObitTableGBTSPDATACopy";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableGBTSPDATAClassInit();

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Use parent class to copy */
  out = (ObitTableGBTSPDATA*)ObitTableCopy ((ObitTable*)in, (ObitTable*)out, err);
  if (err->error) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);

  /* Copy this class  info */
  strncpy (out->object, in->object, MAXKEYCHARTABLEDATA );
  out->scan = in->scan;
  out->utdate = in->utdate;
  out->utcstart = in->utcstart;
  out->utcstop = in->utcstop;
  out->inttime = in->inttime;
  strncpy (out->ctype1, in->ctype1, MAXKEYCHARTABLEDATA );
  strncpy (out->ctype2, in->ctype2, MAXKEYCHARTABLEDATA );
  /* Update class specific info */
  ObitTableGBTSPDATAUpdate (out, err);
    
  return out;
} /* end ObitTableGBTSPDATACopy */

/**
 * Initialize structures and open file.
 * The image descriptor is read if OBIT_IO_ReadOnly or 
 * OBIT_IO_ReadWrite and written to disk if opened OBIT_IO_WriteOnly.
 * After the file has been opened the member, buffer is initialized
 * for reading/storing the table unless member bufferSize is <0.
 * If the requested version ("Ver" in InfoList) is 0 then the highest
 * numbered table of the same type is opened on Read or Read/Write, 
 * or a new table is created on on Write.
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "FileType" OBIT_long scalar = OBIT_IO_FITS 
 *               for file type (see class documentation for details).
 * \li "nRowPIO" OBIT_long scalar = Maximum number of table rows
 *               per transfer, this is the target size for Reads (may be 
 *               fewer) and is used to create buffers.
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *               or OBIT_IO_WriteOnly).
 *               If OBIT_IO_WriteOnly any existing data in the output file
 *               will be lost.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitTableGBTSPDATAOpen (ObitTableGBTSPDATA *in, ObitIOAccess access, 
			  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nRowPIO;
  gchar *routine = "ObitTableGBTSPDATAOpen";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableGBTSPDATAClassInit();

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

   /* Do one row at a time */
   nRowPIO = 1;
   ObitInfoListPut(in->info, "nRowPIO", OBIT_long, dim, (gconstpointer)&nRowPIO, err);
   if (err->error) /* add traceback,return */
     Obit_traceback_val (err, routine, in->name, retCode);
   
   /* use parent class open */
   retCode = ObitTableOpen ((ObitTable*)in, access, err);
   if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
     Obit_traceback_val (err, routine, in->name, retCode);
   
   /* Update class specific info */
   ObitTableGBTSPDATAUpdate (in, err);
   
   return retCode;
} /* end ObitTableGBTSPDATAOpen */

/**
 * Read a table row and return an easily digested version.
 * Scalar values are copied but for array values, pointers 
 * into the data array are returned.
 * \param in       Table to read
 * \param iDATARow   Row number, -1 -> next
 * \param row      Table Row structure to receive data
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitTableGBTSPDATAReadRow  (ObitTableGBTSPDATA *in, olong iDATARow, ObitTableGBTSPDATARow *row,
		     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  odouble   *dRow;
  oint      *iRow;
  gshort    *siRow;
  ofloat    *fRow;
  gchar     *cRow;
  gboolean  *lRow;
  guint8    *bRow;
  gchar *routine = "ObitTableGBTSPDATAReadRow";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myStatus == OBIT_Inactive) {
    Obit_log_error(err, OBIT_Error,
		   "ObitTableGBTSPDATA Table is inactive for  %s ", in->name);
    return retCode;
 }

  /* read row iDATARow */
  retCode = ObitTableRead ((ObitTable*)in, iDATARow, NULL,  err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Typed pointers to row of data */  
  dRow  = (odouble*)in->buffer;
  iRow  = (oint*)in->buffer;
  siRow = (gshort*)in->buffer;
  fRow  = (ofloat*)in->buffer;
  cRow  = (gchar*)in->buffer;
  lRow  = (gboolean*)in->buffer;
  bRow  = (guint8*)in->buffer;
  
  /* Copy scalar fields, for arrays only set pointer*/
  row->subscan = iRow[in->subscanOff];
  row->utdate = iRow[in->utdateOff];
  row->utcstart = dRow[in->utcstartOff];
  row->psrper = dRow[in->psrperOff];
  row->data = fRow + in->dataOff;
  row->status = iRow[in->myDesc->statusOff];

  return retCode;
} /*  end ObitTableGBTSPDATAReadRow */

/**
 * Attach an ObitTableRow to the buffer of an ObitTable.
 * This is only useful prior to filling a row structure in preparation .
 * for a WriteRow operation.  Array members of the Row structure are .
 * pointers to independently allocated memory, this routine allows using .
 * the table IO buffer instead of allocating yet more memory..
 * This routine need only be called once to initialize a Row structure for write..
 * \param in  Table with buffer to be written 
 * \param row Table Row structure to attach 
 * \param err ObitErr for reporting errors.
 */
void 
ObitTableGBTSPDATASetRow  (ObitTableGBTSPDATA *in, ObitTableGBTSPDATARow *row,
		     ObitErr *err)
{
  odouble   *dRow;
  oint      *iRow;
  gshort    *siRow;
  ofloat    *fRow;
  gchar     *cRow;
  gboolean  *lRow;
  guint8    *bRow;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(row, &myRowClassInfo));

  if (in->myStatus == OBIT_Inactive) {
    Obit_log_error(err, OBIT_Error,
		   "DATA Table is inactive for  %s ", in->name);
    return;
 }

  /* Typed pointers to row of data */  
  dRow  = (odouble*)in->buffer;
  iRow  = (oint*)in->buffer;
  siRow = (gshort*)in->buffer;
  fRow  = (ofloat*)in->buffer;
  cRow  = (gchar*)in->buffer;
  lRow  = (gboolean*)in->buffer;
  bRow  = (guint8*)in->buffer;
  
  /* Set row pointers to buffer */
  row->data = fRow + in->dataOff;

} /*  end ObitTableGBTSPDATASetRow */

/**
 * Write a table row.
 * Before calling this routine, the row structure needs to be initialized
 * and filled with data. The array members of the row structure are  
 * pointers to independently allocated memory.  These pointers can be set to the 
 * correct table buffer locations using ObitTableGBTSPDATASetRow  
 * \param in       Table to read
 * \param iDATARow   Row number, -1 -> next
 * \param row      Table Row structure containing data
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitTableGBTSPDATAWriteRow  (ObitTableGBTSPDATA *in, olong iDATARow, ObitTableGBTSPDATARow *row,
		      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gshort    *siRow;
  odouble   *dRow;
  oint      *iRow, i;
  ofloat    *fRow;
  gchar     *cRow;
  gboolean  *lRow;
  guint8    *bRow;
  gchar *routine = "ObitTableGBTSPDATAWriteRow";
  

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myStatus == OBIT_Inactive) {
    Obit_log_error(err, OBIT_Error,
		   "ObitTableGBTSPDATA Table is inactive for %s ", in->name);
    return retCode;
 }

  /* Typed pointers to row of data */  
  dRow  = (odouble*)in->buffer;
  siRow = (gshort*)in->buffer;
  iRow  = (oint*)in->buffer;
  fRow  = (ofloat*)in->buffer;
  cRow  = (gchar*)in->buffer;
  lRow  = (gboolean*)in->buffer;
  bRow  = (guint8*)in->buffer;
  
  /* Make full copy of all data */
  iRow[in->subscanOff] = row->subscan;
  iRow[in->utdateOff] = row->utdate;
  dRow[in->utcstartOff] = row->utcstart;
  dRow[in->psrperOff] = row->psrper;
  if (in->dataCol >= 0) { 
    for (i=0; i<in->myDesc->repeat[in->dataCol]; i++) 
      fRow[in->dataOff+i] = row->data[i];
  } 

  /* copy status */
  iRow[in->myDesc->statusOff] = row->status;
   
  /* Write one row */
  in->myDesc->numRowBuff = 1;
 
  /* Write row iDATARow */
  retCode = ObitTableWrite ((ObitTable*)in, iDATARow, NULL,  err);
  if (err->error) 
    Obit_traceback_val (err, routine,in->name, retCode);

  return retCode;
} /*  end ObitTableGBTSPDATAWriteRow */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitTableGBTSPDATAClose (ObitTableGBTSPDATA *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitTableGBTSPDATAClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;

  /* Update keywords on descriptor if not ReadOnly*/
  if (in->myDesc->access != OBIT_IO_ReadOnly) 
    ObitTableGBTSPDATADumpKey (in, err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Close */
  retCode = ObitTableClose ((ObitTable*)in, err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitTableGBTSPDATAClose */

/*---------------Private functions--------------------------*/
/*----------------  ObitTableGBTSPDATA Row  ----------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTableGBTSPDATARowInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableGBTSPDATARow *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myRowClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  /* Set array members to NULL */
  in->data = NULL;

} /* end ObitTableGBTSPDATARowInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitTableGBTSPDATARow* cast to an Obit*.
 */
void ObitTableGBTSPDATARowClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableGBTSPDATARow *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myRowClassInfo));

  /* delete this class members */
  /* Do not free data array pointers as they were not malloced */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myRowClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTableGBTSPDATARowClear */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableGBTSPDATARowClassInit (void)
{
  if (myRowClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myRowClassInfo.ClassName   = g_strdup(myRowClassName);
  myRowClassInfo.ParentClass = ObitParentGetRowClass();

  /* Set function pointers */
  ObitTableGBTSPDATARowClassInfoDefFn ((gpointer)&myRowClassInfo);
 
  myRowClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableGBTSPDATARowClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableGBTSPDATARowClassInfoDefFn (gpointer inClass)
{
  ObitTableGBTSPDATARowClassInfo *theClass = (ObitTableGBTSPDATARowClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myRowClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myRowClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableGBTSPDATARowClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableGBTSPDATARowClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableGBTSPDATARowGetClass;
  theClass->newObit         = NULL;
  theClass->newObitTableRow = (newObitTableRowFP)newObitTableGBTSPDATARow;
  theClass->ObitCopy        = NULL;
  theClass->ObitClone       = NULL;
  theClass->ObitClear       = (ObitClearFP)ObitTableGBTSPDATARowClear;
  theClass->ObitInit        = (ObitInitFP)ObitTableGBTSPDATARowInit;

} /* end ObitTableGBTSPDATARowClassDefFn */

/*------------------  ObitTableGBTSPDATA  ------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTableGBTSPDATAInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableGBTSPDATA *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */

} /* end ObitTableGBTSPDATAInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitTableGBTSPDATA* cast to an Obit*.
 */
void ObitTableGBTSPDATAClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableGBTSPDATA *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTableGBTSPDATAClear */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableGBTSPDATAClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTableGBTSPDATAClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableGBTSPDATAClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableGBTSPDATAClassInfoDefFn (gpointer inClass)
{
  ObitTableGBTSPDATAClassInfo *theClass = (ObitTableGBTSPDATAClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableGBTSPDATAClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableGBTSPDATAClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableGBTSPDATAGetClass;
  theClass->newObit       = (newObitFP)newObitTableGBTSPDATA;
  theClass->ObitCopy      = (ObitCopyFP)ObitTableGBTSPDATACopy;
  theClass->ObitClone     = (ObitCloneFP)ObitTableClone;
  theClass->ObitClear     = (ObitClearFP)ObitTableGBTSPDATAClear;
  theClass->ObitInit      = (ObitInitFP)ObitTableGBTSPDATAInit;
  theClass->ObitTableConvert = (ObitTableConvertFP)ObitTableGBTSPDATAConvert;
  theClass->ObitTableOpen    = (ObitTableOpenFP)ObitTableGBTSPDATAOpen;
  theClass->ObitTableClose   = (ObitTableCloseFP)ObitTableGBTSPDATAClose;
  theClass->ObitTableRead    = (ObitTableReadFP)ObitTableRead;
  theClass->ObitTableReadSelect = 
    (ObitTableReadSelectFP)ObitTableReadSelect;
  theClass->ObitTableWrite = (ObitTableWriteFP)ObitTableWrite;
  theClass->ObitTableReadRow = 
    (ObitTableReadRowFP)ObitTableGBTSPDATAReadRow;
  theClass->ObitTableWriteRow = 
    (ObitTableWriteRowFP)ObitTableGBTSPDATAWriteRow;

} /* end ObitTableGBTSPDATAClassDefFn */

/**
 * Get table specific information from the infolist or descriptor
 * \param info Table to update
 * \param err  ObitErr for reporting errors.
 */
static void ObitTableGBTSPDATAUpdate (ObitTableGBTSPDATA *in, ObitErr *err)
{
  olong i;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitTableDesc *desc;
   union ObitInfoListEquiv InfoReal; 

 /* error checks */
   g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Get Keywords */
   /* object */
  if (!ObitInfoListGet(in->myDesc->info, "OBJECT", &type, (gint32*)&dim, 
		       (gpointer)&in->object, err)) return;
   /* scan */
  if (!ObitInfoListGet(in->myDesc->info, "SCAN", &type, (gint32*)&dim, 
		       (gpointer)&in->scan, err)) return;
   /* utdate */
  if (!ObitInfoListGet(in->myDesc->info, "UTDATE", &type, (gint32*)&dim, 
		       (gpointer)&in->utdate, err)) return;
   /* utcstart */
  if (!ObitInfoListGet(in->myDesc->info, "UTCSTART", &type, (gint32*)&dim, 
		       (gpointer)&InfoReal, err)) return;
  if (type==OBIT_double) in->utcstart = InfoReal.dbl;
  if (type==OBIT_float)  in->utcstart = InfoReal.flt;
   /* utcstop */
  if (!ObitInfoListGet(in->myDesc->info, "UTCSTOP", &type, (gint32*)&dim, 
		       (gpointer)&InfoReal, err)) return;
  if (type==OBIT_double) in->utcstop = InfoReal.dbl;
  if (type==OBIT_float)  in->utcstop = InfoReal.flt;
   /* inttime */
  if (!ObitInfoListGet(in->myDesc->info, "INTTIME", &type, (gint32*)&dim, 
		       (gpointer)&InfoReal, err)) return;
  if (type==OBIT_double) in->inttime = InfoReal.dbl;
  if (type==OBIT_float)  in->inttime = InfoReal.flt;
   /* ctype1 */
  strncpy (in->ctype1, "STATE", MAXKEYCHARTABLEDATA); 
  ObitInfoListGetTest(in->myDesc->info, "CTYPE1", &type, (gint32*)&dim, 
		       (gpointer)&in->ctype1);
   /* ctype2 */
  strncpy (in->ctype2, "RECEIVER", MAXKEYCHARTABLEDATA); 
  ObitInfoListGetTest(in->myDesc->info, "CTYPE2", &type, (gint32*)&dim, 
		       (gpointer)&in->ctype2);

  /* initialize column numbers/offsets */
  in->subscanOff = -1;
  in->subscanCol = -1;
  in->utdateOff = -1;
  in->utdateCol = -1;
  in->utcstartOff = -1;
  in->utcstartCol = -1;
  in->psrperOff = -1;
  in->psrperCol = -1;
  in->dataOff = -1;
  in->dataCol = -1;
  /* Find columns and set offsets */
  desc = in->myDesc;
  if (desc->FieldName) {
    for (i=0; i<desc->nfield; i++) {
      if (!strncmp (desc->FieldName[i], "SUBSCAN ", 8)) {
	 in->subscanOff = desc->offset[i];
 	 in->subscanCol = i;
      }
      if (!strncmp (desc->FieldName[i], "UTDATE ", 7)) {
	 in->utdateOff = desc->offset[i];
 	 in->utdateCol = i;
      }
      if (!strncmp (desc->FieldName[i], "UTCSTART ", 9)) {
	 in->utcstartOff = desc->offset[i];
 	 in->utcstartCol = i;
      }
      if (!strncmp (desc->FieldName[i], "PSRPER ", 7)) {
	 in->psrperOff = desc->offset[i];
 	 in->psrperCol = i;
      }
      if (!strncmp (desc->FieldName[i], "DATA ", 5)) {
	 in->dataOff = desc->offset[i];
 	 in->dataCol = i;
      }
     }
  }

  Obit_return_if_fail((in->subscanOff > -1), err,
       "ObitTableDATAUpdate: Could not find column subscan");
  Obit_return_if_fail((in->utdateOff > -1), err,
       "ObitTableDATAUpdate: Could not find column utdate");
  Obit_return_if_fail((in->utcstartOff > -1), err,
       "ObitTableDATAUpdate: Could not find column utcstart");
  Obit_return_if_fail((in->psrperOff > -1), err,
       "ObitTableDATAUpdate: Could not find column psrper");
} /* end ObitTableGBTSPDATAUpdate */

/**
 * Copy table specific (keyword) information  to infolist.
 * \param info Table to update
 * \param err  ObitErr for reporting errors.
 */
static void ObitTableGBTSPDATADumpKey (ObitTableGBTSPDATA *in, ObitErr *err)
{
  ObitInfoList *info=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

 /* error checks */
   g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Set Keywords */
  if (in->myIO!=NULL) info = ((ObitTableDesc*)(in->myIO->myDesc))->info;
  else info = in->myDesc->info;
  /* object */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEDATA;
  ObitInfoListAlwaysPut(info, "OBJECT", type, (gint32*)&dim, 
		  (gpointer)&in->object);
  /* scan */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "SCAN", type, (gint32*)&dim, 
		  (gpointer)&in->scan);
  /* utdate */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "UTDATE", type, (gint32*)&dim, 
		  (gpointer)&in->utdate);
  /* utcstart */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "UTCSTART", type, (gint32*)&dim, 
		  (gpointer)&in->utcstart);
  /* utcstop */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "UTCSTOP", type, (gint32*)&dim, 
		  (gpointer)&in->utcstop);
  /* inttime */
  type  = OBIT_double;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "INTTIME", type, (gint32*)&dim, 
		  (gpointer)&in->inttime);
  /* ctype1 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEDATA;
  ObitInfoListAlwaysPut(info, "CTYPE1", type, (gint32*)&dim, 
		  (gpointer)&in->ctype1);
  /* ctype2 */
  type  = OBIT_string;
  dim[0] = MAXKEYCHARTABLEDATA;
  ObitInfoListAlwaysPut(info, "CTYPE2", type, (gint32*)&dim, 
		  (gpointer)&in->ctype2);
   
} /* end ObitTableGBTSPDATADumpKey */
