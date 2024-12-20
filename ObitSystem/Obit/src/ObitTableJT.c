/* $Id:  $   */
/* DO NOT EDIT - file generated by ObitTables.pl                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C)  2024                                              */
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
#include "ObitTableJT.h"
#include "ObitTableList.h"
#include "ObitData.h"

/*----------------Obit:  Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableJT.c
 * ObitTableJT class function definitions.
 *
 * This class is derived from the #ObitTable class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitTableJT";

/**  Function to obtain parent Table ClassInfo - ObitTable */
static ObitGetClassFP ObitParentGetClass = ObitTableGetClass;

/** name of the Row class defined in this file */
static gchar *myRowClassName = "ObitTableJTRow";

/**  Function to obtain parent TableRow ClassInfo */
static ObitGetClassFP ObitParentGetRowClass = ObitTableRowGetClass;

/*--------------- File Global Variables  ----------------*/
/*----------------  Table Row  ----------------------*/
/**
 * ClassInfo structure ObitTableClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableJTRowClassInfo myRowClassInfo = {FALSE};

/*------------------  Table  ------------------------*/
/**
 * ClassInfo structure ObitTableJTClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableJTClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated Row object. */
void  ObitTableJTRowInit  (gpointer in);

/** Private: Deallocate Row members. */
void  ObitTableJTRowClear (gpointer in);

/** Private: Initialize newly instantiated object. */
void  ObitTableJTInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTableJTClear (gpointer in);

/** Private: update table specific info */
static void ObitTableJTUpdate (ObitTableJT *in, ObitErr *err);

/** Private: copy table keywords to descriptor info list */
static void ObitTableJTDumpKey (ObitTableJT *in, ObitErr *err);

/** Private: Set Class function pointers */
static void ObitTableJTClassInfoDefFn (gpointer inClass);

/** Private: Set Row Class function pointers */
static void ObitTableJTRowClassInfoDefFn (gpointer inClass);
/*----------------------Public functions---------------------------*/

/*------------------  Table Row ------------------------*/
/**
 * Constructor.
 * If table is open and for write, the row is attached to the buffer
 * Initializes Row class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTableJTRow* newObitTableJTRow (ObitTableJT *table)
{
  ObitTableJTRow* out;
  oint      *iRow;
  ofloat    *fRow;

  /* Class initialization if needed */
  if (!myRowClassInfo.initialized) ObitTableJTRowClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTableJTRow));

  /* initialize values */
  out->name = g_strdup("TableJT Row");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myRowClassInfo;

  /* initialize other stuff */
  ObitTableJTRowInit((gpointer)out);
  out->myTable   = (ObitTable*)ObitTableRef((ObitTable*)table);

  /* If writing attach to buffer */
  if ((table->buffer) && (table->myDesc->access != OBIT_IO_ReadOnly) &&
      (table->myStatus != OBIT_Inactive)) {
    /* Typed pointers to row of data */  
    iRow  = (oint*)table->buffer;
    fRow  = (ofloat*)table->buffer;
  
    /* Set row pointers to buffer */
    out->Jones = fRow + table->JonesOff;
    out->Weight = fRow + table->WeightOff;
    out->RefAnt = iRow + table->RefAntOff;
  } /* end attaching row to table buffer */

 return out;
} /* end newObitTableJTRow */

/**
 * Returns ClassInfo pointer for the Row class.
 * \return pointer to the Row class structure.
 */
gconstpointer ObitTableJTRowGetClass (void)
{
  /* Class initialization if needed */
  if (!myRowClassInfo.initialized) ObitTableJTRowClassInit();
  return (gconstpointer)&myRowClassInfo;
} /* end ObitTableJTRowGetClass */

/*------------------  Table  ------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTableJT* newObitTableJT (gchar* name)
{
  ObitTableJT* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableJTClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTableJT));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTableJTInit((gpointer)out);

 return out;
} /* end newObitTableJT */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitTableJTGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableJTClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitJTGetClass */

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
 * \param numIF The number of IFs
 * \param numChan Number of frequency channels
 * \param err Error stack, returns if not empty.
 * \return the new object, NULL on failure.
 */
ObitTableJT* newObitTableJTValue (gchar* name, ObitData *file, olong *ver,
 	                    ObitIOAccess access,
  		     oint numIF, oint numChan,
		     ObitErr *err)
{
  ObitTableJT* out=NULL;
  ObitTable *testTab=NULL;
  ObitTableDesc *desc=NULL;
  ObitTableList *list=NULL;
  ObitInfoList  *info=NULL;
  gboolean exist, optional;
  olong colNo, i, ncol, highVer;
  ObitIOCode retCode;
  gchar *tabType = "AIPS JT";
  gchar *routine = "newObitTableJTValue";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitDataIsA(file));

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableJTClassInit();

  /* Check if the table already exists */
  /* Get TableList */
  list = ((ObitData*)file)->tableList;
  info = ((ObitData*)file)->info;

  /* Get highest version number if not specified */
  if (*ver==0) { 
    highVer = ObitTableListGetHigh (list, "AIPS JT");
    if (access==OBIT_IO_ReadOnly) *ver = highVer;
    else if (access==OBIT_IO_ReadWrite) *ver = highVer;
    else if (access==OBIT_IO_WriteOnly) *ver = highVer+1;
  }
  /* See if it already exists */
  exist = FALSE;
  if (*ver>0) { /* has to be specified */
    exist = ObitTableListGet(list, tabType, ver, &testTab, err);
    if (err->error) /* add traceback,return */
      Obit_traceback_val (err, routine,"", out);
  
    /* if readonly, it must exist to proceed */
    if ((access==OBIT_IO_ReadOnly) && !exist) return out;
    if (testTab!=NULL) { /* it exists, use it if is an ObitTableJT */
      if (ObitTableJTIsA(testTab)) { /* it is an ObitTableJT */
	out = ObitTableRef(testTab);
      } else { /* needs conversion */
 	out = ObitTableJTConvert(testTab);
	/* Update the TableList */
	ObitTableListPut(list, tabType, ver, (ObitTable*)out, err);
	if (err->error) /* add traceback,return */
	  Obit_traceback_val (err, routine,"", out);
      }
      testTab = ObitTableUnref(testTab); /* remove reference */
      return out; /* done */
    }
  } /* end of test for previously existing table */
  
  /* If access is ReadOnly make sure one exists */
  if (access==OBIT_IO_ReadOnly) { 
    highVer = ObitTableListGetHigh (list, "AIPS JT");
    if (highVer<=0) return out;
  }
  
  /* create basal table */
  testTab = newObitDataTable ((ObitData*)file, access, tabType,
			       ver, err);
  if (err->error) Obit_traceback_val (err, routine,"", out);
  
  /* likely need to convert */
  if (ObitTableJTIsA(testTab)) { 
    out = ObitTableRef(testTab);
  } else { /* needs conversion */
    out = ObitTableJTConvert(testTab);
  }
  testTab = ObitTableUnref(testTab); /* remove reference */

  /* Update the TableList */
  ObitTableListPut(list, tabType, ver, (ObitTable*)out, err);
  if (err->error) Obit_traceback_val (err, routine,"", out);

  /* if it previously existed merely return it */
  if (exist) return out; 

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* Set values */
  out->numIF = MAX (0, numIF);
  out->numChan = MAX (0, numChan);
  out->revision = 11;
  out->numAnt = 1;

  /* initialize descriptor */
  desc = out->myDesc;
  /* How many columns actually in table? */
  ncol = 9 + out->numIF*0 + out->numChan*0 ;
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
  desc->FieldName[colNo] = g_strdup("TIME    ");
  desc->FieldUnit[colNo] = g_strdup("DAYS");
  desc->type[colNo] = OBIT_double;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("TIME INTERVAL");
  desc->FieldUnit[colNo] = g_strdup("DAYS");
  desc->type[colNo] = OBIT_float;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("SOURCE ID");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("ANTENNA NO.");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("SUBARRAY");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  desc->FieldName[colNo] = g_strdup("FREQ ID");
  desc->FieldUnit[colNo] = g_strdup("");
  desc->type[colNo] = OBIT_oint;
  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
  colNo++;
  optional = FALSE;
  if ((8 > 0) || (!optional)) {
    desc->FieldName[colNo] = g_strdup("JONES");
    desc->FieldUnit[colNo] = g_strdup("");
    desc->type[colNo] = OBIT_float;
    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
    desc->dim[colNo][0] = 8;
    desc->dim[colNo][1] = numChan;
    desc->dim[colNo][2] = numIF;
    colNo++;
  }
  optional = FALSE;
  if ((numChan > 0) || (!optional)) {
    desc->FieldName[colNo] = g_strdup("WEIGHT");
    desc->FieldUnit[colNo] = g_strdup("");
    desc->type[colNo] = OBIT_float;
    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
    desc->dim[colNo][0] = numChan;
    desc->dim[colNo][1] = numIF;
    colNo++;
  }
  optional = FALSE;
  if ((numChan > 0) || (!optional)) {
    desc->FieldName[colNo] = g_strdup("REFANT");
    desc->FieldUnit[colNo] = g_strdup("");
    desc->type[colNo] = OBIT_oint;
    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;
    desc->dim[colNo][0] = numChan;
    desc->dim[colNo][1] = numIF;
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
  ObitTableJTDumpKey (out, err);
 
  /* index table descriptor */
  ObitTableDescIndex (desc);

  /* Open and Close to fully instantiate */
  retCode = ObitTableJTOpen(out, OBIT_IO_WriteOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out);    
  
  retCode = ObitTableJTClose(out, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out); 

  /* Force update of disk resident info */
  retCode = ObitIOUpdateTables (((ObitData*)file)->myIO, info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, out->name, out); 
  
 return out;
} /* end newObitTableJTValue */

/**
 * Convert an ObitTable to an ObitTableJT.
 * New object will have references to members of in.
 * \param in  The object to copy, will still exist afterwards 
 *            and should be Unrefed if not needed.
 * \return pointer to the new object.
 */
ObitTableJT* ObitTableJTConvert (ObitTable* in)
{
  ObitTableJT *out;

  /* error check */
  g_assert(ObitTableIsA(in));

  /* create basic object */
  out = newObitTableJT(in->name);

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
} /* end ObitTableJTConvert */


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
ObitTableJT* ObitTableJTCopy (ObitTableJT *in, ObitTableJT *out, ObitErr *err)
{
  gchar *routine = "ObitTableJTCopy";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableJTClassInit();

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Use parent class to copy */
  out = (ObitTableJT*)ObitTableCopy ((ObitTable*)in, (ObitTable*)out, err);
  if (err->error) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);

  /* Copy this class  info */
  out->revision = in->revision;
  out->numIF = in->numIF;
  out->numChan = in->numChan;
  out->numAnt = in->numAnt;
  /* Update class specific info */
  ObitTableJTUpdate (out, err);
    
  return out;
} /* end ObitTableJTCopy */

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
 * \li "FileType" OBIT_long scalar = OBIT_IO_FITS or OBIT_IO_AIPS 
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
ObitIOCode ObitTableJTOpen (ObitTableJT *in, ObitIOAccess access, 
			  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nRowPIO;
  gchar *routine = "ObitTableJTOpen";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableJTClassInit();

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
   ObitTableJTUpdate (in, err);
   
   return retCode;
} /* end ObitTableJTOpen */

/**
 * Read a table row and return an easily digested version.
 * Scalar values are copied but for array values, pointers 
 * into the data array are returned.
 * \param in       Table to read
 * \param iJTRow   Row number, -1 -> next
 * \param row      Table Row structure to receive data
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitTableJTReadRow  (ObitTableJT *in, olong iJTRow, ObitTableJTRow *row,
		     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  odouble   *dRow;
  oint      *iRow;
  ofloat    *fRow;
  gchar *routine = "ObitTableJTReadRow";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myStatus == OBIT_Inactive) {
    Obit_log_error(err, OBIT_Error,
		   "AIPS JT Table is inactive for  %s ", in->name);
    return retCode;
 }

  /* read row iJTRow */
  retCode = ObitTableRead ((ObitTable*)in, iJTRow, NULL,  err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Typed pointers to row of data */  
  dRow  = (odouble*)in->buffer;
  iRow  = (oint*)in->buffer;
  fRow  = (ofloat*)in->buffer;
  
  /* Copy scalar fields, for arrays only set pointer*/
  row->Time = dRow[in->TimeOff];
  row->TimeI = fRow[in->TimeIOff];
  row->SourID = iRow[in->SourIDOff];
  row->antNo = iRow[in->antNoOff];
  row->SubA = iRow[in->SubAOff];
  row->FreqID = iRow[in->FreqIDOff];
  row->Jones = fRow + in->JonesOff;
  row->Weight = fRow + in->WeightOff;
  row->RefAnt = iRow + in->RefAntOff;
  row->status = iRow[in->myDesc->statusOff];

  return retCode;
} /*  end ObitTableJTReadRow */

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
ObitTableJTSetRow  (ObitTableJT *in, ObitTableJTRow *row,
		     ObitErr *err)
{
  oint      *iRow;
  ofloat    *fRow;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(row, &myRowClassInfo));

  if (in->myStatus == OBIT_Inactive) {
    Obit_log_error(err, OBIT_Error,
		   "JT Table is inactive for  %s ", in->name);
    return;
 }

  /* Typed pointers to row of data */  
  iRow  = (oint*)in->buffer;
  fRow  = (ofloat*)in->buffer;
  
  /* Set row pointers to buffer */
  row->Jones = fRow + in->JonesOff;
  row->Weight = fRow + in->WeightOff;
  row->RefAnt = iRow + in->RefAntOff;

} /*  end ObitTableJTSetRow */

/**
 * Write a table row.
 * Before calling this routine, the row structure needs to be initialized
 * and filled with data. The array members of the row structure are  
 * pointers to independently allocated memory.  These pointers can be set to the 
 * correct table buffer locations using ObitTableJTSetRow  
 * \param in       Table to read
 * \param iJTRow   Row number, -1 -> next
 * \param row Table Row structure containing data
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode 
ObitTableJTWriteRow  (ObitTableJT *in, olong iJTRow, ObitTableJTRow *row,
		      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  odouble   *dRow;
  oint      *iRow, i;
  ofloat    *fRow;
  gchar *routine = "ObitTableJTWriteRow";
  

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myStatus == OBIT_Inactive) {
    Obit_log_error(err, OBIT_Error,
		   "AIPS JT Table is inactive for %s ", in->name);
    return retCode;
 }

  /* Typed pointers to row of data */  
  dRow  = (odouble*)in->buffer;
  iRow  = (oint*)in->buffer;
  fRow  = (ofloat*)in->buffer;
  
  /* Make full copy of all data */
  dRow[in->TimeOff] = row->Time;
  fRow[in->TimeIOff] = row->TimeI;
  iRow[in->SourIDOff] = row->SourID;
  iRow[in->antNoOff] = row->antNo;
  iRow[in->SubAOff] = row->SubA;
  iRow[in->FreqIDOff] = row->FreqID;
  if (in->JonesCol >= 0) { 
    for (i=0; i<in->myDesc->repeat[in->JonesCol]; i++) 
      fRow[in->JonesOff+i] = row->Jones[i];
  } 
  if (in->WeightCol >= 0) { 
    for (i=0; i<in->myDesc->repeat[in->WeightCol]; i++) 
      fRow[in->WeightOff+i] = row->Weight[i];
  } 
  if (in->RefAntCol >= 0) { 
    for (i=0; i<in->myDesc->repeat[in->RefAntCol]; i++) 
      iRow[in->RefAntOff+i] = row->RefAnt[i];
  } 

  /* copy status */
  iRow[in->myDesc->statusOff] = row->status;
   
  /* Write one row */
  in->myDesc->numRowBuff = 1;
 
  /* Write row iJTRow */
  retCode = ObitTableWrite ((ObitTable*)in, iJTRow, NULL,  err);
  if (err->error) 
    Obit_traceback_val (err, routine,in->name, retCode);

  return retCode;
} /*  end ObitTableJTWriteRow */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitTableJTClose (ObitTableJT *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitTableJTClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;

  /* Update keywords on descriptor if not ReadOnly*/
  if (in->myDesc->access != OBIT_IO_ReadOnly) 
    ObitTableJTDumpKey (in, err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Close */
  retCode = ObitTableClose ((ObitTable*)in, err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitTableJTClose */

/*---------------Private functions--------------------------*/
/*----------------  TableJT Row  ----------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTableJTRowInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableJTRow *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myRowClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  /* Set array members to NULL */
  in->Jones = NULL;
  in->Weight = NULL;
  in->RefAnt = NULL;

} /* end ObitTableJTRowInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitTableJTRow* cast to an Obit*.
 */
void ObitTableJTRowClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableJTRow *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myRowClassInfo));

  /* delete this class members */
  /* Do not free data array pointers as they were not malloced */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myRowClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTableJTRowClear */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableJTRowClassInit (void)
{
  if (myRowClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myRowClassInfo.ClassName   = g_strdup(myRowClassName);
  myRowClassInfo.ParentClass = ObitParentGetRowClass();

  /* Set function pointers */
  ObitTableJTRowClassInfoDefFn ((gpointer)&myRowClassInfo);
 
  myRowClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableJTRowClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableJTRowClassInfoDefFn (gpointer inClass)
{
  ObitTableJTRowClassInfo *theClass = (ObitTableJTRowClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myRowClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myRowClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableJTRowClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableJTRowClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableJTRowGetClass;
  theClass->newObit         = NULL;
  theClass->newObitTableRow = (newObitTableRowFP)newObitTableJTRow;
  theClass->ObitCopy        = NULL;
  theClass->ObitClone       = NULL;
  theClass->ObitClear       = (ObitClearFP)ObitTableJTRowClear;
  theClass->ObitInit        = (ObitInitFP)ObitTableJTRowInit;

} /* end ObitTableJTRowClassDefFn */

/*------------------  TableJT  ------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTableJTInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableJT *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */

} /* end ObitTableJTInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitTableJT* cast to an Obit*.
 */
void ObitTableJTClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableJT *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTableJTClear */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableJTClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTableJTClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableJTClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableJTClassInfoDefFn (gpointer inClass)
{
  ObitTableJTClassInfo *theClass = (ObitTableJTClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableJTClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableJTClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableJTGetClass;
  theClass->newObit       = (newObitFP)newObitTableJT;
  theClass->ObitCopy      = (ObitCopyFP)ObitTableJTCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitTableClone;
  theClass->ObitClear     = (ObitClearFP)ObitTableJTClear;
  theClass->ObitInit      = (ObitInitFP)ObitTableJTInit;
  theClass->ObitTableConvert = (ObitTableConvertFP)ObitTableJTConvert;
  theClass->ObitTableOpen    = (ObitTableOpenFP)ObitTableJTOpen;
  theClass->ObitTableClose   = (ObitTableCloseFP)ObitTableJTClose;
  theClass->ObitTableRead    = (ObitTableReadFP)ObitTableRead;
  theClass->ObitTableReadSelect = 
    (ObitTableReadSelectFP)ObitTableReadSelect;
  theClass->ObitTableWrite = (ObitTableWriteFP)ObitTableWrite;
  theClass->ObitTableReadRow = 
    (ObitTableReadRowFP)ObitTableJTReadRow;
  theClass->ObitTableWriteRow = 
    (ObitTableWriteRowFP)ObitTableJTWriteRow;

} /* end ObitTableJTClassDefFn */

/**
 * Get table specific information from the infolist or descriptor
 * \param info Table to update
 * \param err  ObitErr for reporting errors.
 */
static void ObitTableJTUpdate (ObitTableJT *in, ObitErr *err)
{
  olong i;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitTableDesc *desc;
   

 /* error checks */
   g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Get Keywords */
   /* revision */
  in->revision = 11; 
  ObitInfoListGetTest(in->myDesc->info, "REVISION", &type, dim, 
		       (gpointer)&in->revision);
   /* numIF */
  if (!ObitInfoListGet(in->myDesc->info, "NO_IF", &type, dim, 
		       (gpointer)&in->numIF, err)) return;
   /* numChan */
  if (!ObitInfoListGet(in->myDesc->info, "NO_CHAN", &type, dim, 
		       (gpointer)&in->numChan, err)) return;
   /* numAnt */
  in->numAnt = 1; 
  ObitInfoListGetTest(in->myDesc->info, "NO_ANT", &type, dim, 
		       (gpointer)&in->numAnt);

  /* initialize column numbers/offsets */
  in->TimeOff = -1;
  in->TimeCol = -1;
  in->TimeIOff = -1;
  in->TimeICol = -1;
  in->SourIDOff = -1;
  in->SourIDCol = -1;
  in->antNoOff = -1;
  in->antNoCol = -1;
  in->SubAOff = -1;
  in->SubACol = -1;
  in->FreqIDOff = -1;
  in->FreqIDCol = -1;
  in->JonesOff = -1;
  in->JonesCol = -1;
  in->WeightOff = -1;
  in->WeightCol = -1;
  in->RefAntOff = -1;
  in->RefAntCol = -1;
  /* Find columns and set offsets */
  desc = in->myDesc;
  if (desc->FieldName) {
    for (i=0; i<desc->nfield; i++) {
      if (!strncmp (desc->FieldName[i], "TIME    ", 8)) {
	 in->TimeOff = desc->offset[i];
 	 in->TimeCol = i;
      }
      if (!strncmp (desc->FieldName[i], "TIME INTERVAL", 13)) {
	 in->TimeIOff = desc->offset[i];
 	 in->TimeICol = i;
      }
      if (!strncmp (desc->FieldName[i], "SOURCE ID", 9)) {
	 in->SourIDOff = desc->offset[i];
 	 in->SourIDCol = i;
      }
      if (!strncmp (desc->FieldName[i], "ANTENNA NO.", 11)) {
	 in->antNoOff = desc->offset[i];
 	 in->antNoCol = i;
      }
      if (!strncmp (desc->FieldName[i], "SUBARRAY", 8)) {
	 in->SubAOff = desc->offset[i];
 	 in->SubACol = i;
      }
      if (!strncmp (desc->FieldName[i], "FREQ ID", 7)) {
	 in->FreqIDOff = desc->offset[i];
 	 in->FreqIDCol = i;
      }
      if (!strncmp (desc->FieldName[i], "JONES", 5)) {
	 in->JonesOff = desc->offset[i];
 	 in->JonesCol = i;
      }
      if (!strncmp (desc->FieldName[i], "WEIGHT", 6)) {
	 in->WeightOff = desc->offset[i];
 	 in->WeightCol = i;
      }
      if (!strncmp (desc->FieldName[i], "REFANT", 6)) {
	 in->RefAntOff = desc->offset[i];
 	 in->RefAntCol = i;
      }
     }
  }

  /* Check required columns */
  Obit_return_if_fail((in->TimeOff > -1), err,
       "ObitTableJTUpdate: Could not find column Time");
  Obit_return_if_fail((in->TimeIOff > -1), err,
       "ObitTableJTUpdate: Could not find column TimeI");
  Obit_return_if_fail((in->SourIDOff > -1), err,
       "ObitTableJTUpdate: Could not find column SourID");
  Obit_return_if_fail((in->antNoOff > -1), err,
       "ObitTableJTUpdate: Could not find column antNo");
  Obit_return_if_fail((in->SubAOff > -1), err,
       "ObitTableJTUpdate: Could not find column SubA");
  Obit_return_if_fail((in->FreqIDOff > -1), err,
       "ObitTableJTUpdate: Could not find column FreqID");
  Obit_return_if_fail((in->JonesOff > -1), err,
       "ObitTableJTUpdate: Could not find column Jones");
  Obit_return_if_fail((in->WeightOff > -1), err,
       "ObitTableJTUpdate: Could not find column Weight");
} /* end ObitTableJTUpdate */

/**
 * Copy table specific (keyword) information  to infolist.
 * \param info Table to update
 * \param err  ObitErr for reporting errors.
 */
static void ObitTableJTDumpKey (ObitTableJT *in, ObitErr *err)
{
  ObitInfoList *info=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

 /* error checks */
   g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Set Keywords */
  /*if (in->myIO!=NULL) info = ((ObitTableDesc*)(in->myIO->myDesc))->info;
    else info = in->myDesc->info;*/
  info = in->myDesc->info;
  /* revision */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "REVISION", type, dim, 
		  (gpointer)&in->revision);
  /* numIF */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "NO_IF", type, dim, 
		  (gpointer)&in->numIF);
  /* numChan */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "NO_CHAN", type, dim, 
		  (gpointer)&in->numChan);
  /* numAnt */
  type  = OBIT_oint;
  dim[0] = 1;
  ObitInfoListAlwaysPut(info, "NO_ANT", type, dim, 
		  (gpointer)&in->numAnt);
   
} /* end ObitTableJTDumpKey */
