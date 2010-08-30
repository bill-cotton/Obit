/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2010                                          */
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

#include "ObitTable.h"
#include "ObitIOTableFITS.h"
#include "ObitIOTableAIPS.h"
#include "ObitTableDesc.h"
#include "ObitTableSel.h"
#include "ObitMem.h"
#include "ObitData.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTable.c
 * ObitTable class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitTable";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/** name of the Row class defined in this file */
static gchar *myRowClassName = "ObitTableRow";

/** Function to obtain parent ClassInfo */
ObitGetClassFP ObitParentGetRowClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/*----------------  Table Row  ----------------------*/
/**
 * ClassInfo structure ObitTableClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableRowClassInfo myRowClassInfo = {FALSE};

/*------------------  Table  ------------------------*/
/**
 * ClassInfo structure ObitTableClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated Row object. */
void  ObitTableRowInit  (gpointer in);

/** Private: Deallocate Row members. */
void  ObitTableRowClear (gpointer in);

/** Private: Initialize newly instantiated object. */
void  ObitTableInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTableClear (gpointer in);

/** Private: Read selection parameters from ObitInfoList. */
static void ObitTableGetSelect (ObitInfoList *info, ObitTableSel *sel,
				ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitTableClassInfoDefFn (gpointer inClass);

/** Private: Set Row Class function pointers. */
static void ObitTableRowClassInfoDefFn (gpointer inClass);


/*----------------------Public functions---------------------------*/
/*------------------  Table Row ------------------------*/
/**
 * Constructor.
 * Initializes Row class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTableRow* newObitTableRow (ObitTable *table)
{
  ObitTableRow* out;
  gchar *name;

  /* Class initialization if needed */
  if (!myRowClassInfo.initialized) ObitTableRowClassInit();

  /* allocate/init structure */
  name =  g_strconcat ("TRow:", table->tabType, NULL);
  out = ObitMemAlloc0Name(sizeof(ObitTableRow), name);
  g_free(name);

  /* initialize values */
  out->name = g_strdup("Table Row");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myRowClassInfo;

  /* initialize other stuff */
  ObitTableRowInit((gpointer)out);
  out->myTable   = ObitTableRef((ObitTable*)table);

 return out;
} /* end newObitTableRow */

/**
 * Returns ClassInfo pointer for the Row class.
 * \return pointer to the Row class structure.
 */
gconstpointer ObitTableRowGetClass (void)
{
  /* Class initialization if needed */
  if (!myRowClassInfo.initialized) ObitTableRowClassInit();

  return (gconstpointer)&myRowClassInfo;
} /* end ObitTableRowGetClass */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableRowClassInit (void)
{
  if (myRowClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myRowClassInfo.ClassName   = g_strdup(myRowClassName);
  myRowClassInfo.ParentClass = ObitParentGetRowClass();

  /* Set function pointers */
  ObitTableRowClassInfoDefFn ((gpointer)&myRowClassInfo);
 
  myRowClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableRowClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableRowClassInfoDefFn (gpointer inClass)
{
  ObitTableRowClassInfo *theClass = (ObitTableRowClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myRowClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableRowClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableRowClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableRowGetClass;
  theClass->newObit         = NULL;
  theClass->newObitTableRow = (newObitTableRowFP)newObitTableRow;
  theClass->ObitCopy        = NULL;
  theClass->ObitClone       = NULL;  
  theClass->ObitClear       = (ObitClearFP)ObitTableRowClear;
  theClass->ObitInit        = (ObitInitFP)ObitTableRowInit;
} /* end ObitTableRowClassDefFn */

/*------------------  Table  ------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTable* newObitTable (gchar* name)
{
  ObitTable* out;
  gchar *lname;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableClassInit();

  /* allocate/init structure */
  lname =  g_strconcat ("Tab:", name, NULL);
  out = ObitMemAlloc0Name(sizeof(ObitTable), lname);
  g_free(lname);

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTableInit((gpointer)out);

 return out;
} /* end newObitTable */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitTableGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitTableGetClass */

/**
 * Delete underlying files and the basic object.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitTable* ObitTableZap (ObitTable *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitTableZap";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return in;
  if (in==NULL) return NULL;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));


  /* Open and close if needed */
  if (in->myIO==NULL) {
    in->bufferSize = -1;  /* Don't need buffer */
    retCode = ObitTableOpen (in, OBIT_IO_ReadWrite, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, in);
    retCode = ObitTableClose (in, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, in);
  }
  
  /* Delete the file */
  ObitIOZap (in->myIO, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);

  in->myIO = ObitIOUnref(in->myIO);       /*  delete IO */
  in->info = ObitInfoListUnref(in->info);  /* delete infoList */

  /* Get memory resident bits as well - loop until truely deleted
  Bad Idea while (in->ReferenceCount>1) ObitTableUnref(in);*/
  in = ObitTableUnref(in);

  return in;
} /* end ObitTableZap */

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
ObitTable* ObitTableCopy (ObitTable *in, ObitTable *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitIOCode iretCode, oretCode;
  gboolean oldExist;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nRowPIO;
  gchar *outName;
  gchar *routine = "ObitTableCopy";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitTable(outName);
    if (outName) g_free(outName); outName = NULL;
    if (out->tabType) g_free(out->tabType); out->tabType = NULL;
    if (in->tabType) out->tabType = g_strdup(in->tabType);
 }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes additions only if out newly created */
  if (!oldExist) {
    /* copy */
    out->myDesc = (gpointer)ObitTableDescCopy(in->myDesc, out->myDesc, err);
    /* Don't copy selector */
    if (out->mySel) out->mySel = ObitUnref (out->mySel);
    out->mySel = newObitTableSel (out->name);
    out->info = ObitInfoListUnref(out->info);
    out->info = ObitInfoListRef(in->info);
    /* don't copy ObitTableSel, ObitThread */
  }

  /* If the output object was created this call it cannot be fully
     defined so we're done */
  if (!oldExist) return out;

  /* Only works one row at a time */
  nRowPIO = 1;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(in->info, "nRowPIO", OBIT_long, dim, &nRowPIO);

  /* if input has file designated, copy data */
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitTableOpen (in, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((iretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_val (err, routine,in->name, out);
  }

  /* copy Descriptor - this time with full information */
  out->myDesc = ObitTableDescCopy(in->myDesc, out->myDesc, err);


 /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (out->buffer) ObitIOFreeBuffer(out->buffer); /* free existing */
  out->buffer     = NULL;
  out->bufferSize = -1;

  /* test open output */
  oretCode = ObitTableOpen (out, OBIT_IO_WriteOnly, err);
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    out->buffer = NULL;
    out->bufferSize = 0;
    Obit_traceback_val (err, routine,in->name, out);
  }

  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitTableRead (in, -1, in->buffer, err);
    /* How many */
    out->myDesc->numRowBuff = in->myDesc->numRowBuff;
    if (iretCode!=OBIT_IO_OK) break;
    oretCode = ObitTableWrite (out, -1, in->buffer, err);
  }
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  out->buffer = NULL;
  out->bufferSize = 0;
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);
  
  /* close files to be sure */
  iretCode = ObitTableClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,in->name, out);
  
  /* close files to be sure */
  oretCode = ObitTableClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,out->name, out);
  
  return out;
} /* end ObitTableCopy */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
ObitTable* ObitTableClone  (ObitTable *in, ObitTable *out)
{
  const ObitClassInfo *myClass, *ParentClass;
  gboolean oldExist;
  ObitErr    *err = NULL;
  gchar *outName;
  gchar *routine = "ObitTableClone";

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Clone: ",in->name,NULL);
    out = newObitTable(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* shallow copy any parent class members */
   myClass     = in->ClassInfo;
   ParentClass = myClass->ParentClass;
   if ((ParentClass!=NULL) && (ParentClass->ObitClone!=NULL)
       && (ParentClass->ObitClone!=(ObitCloneFP)ObitTableClone))
     ParentClass->ObitClone ((Obit*)in, (Obit*)out);

   if (!oldExist) { /* only copy ObitInfoList if just created */
     out->info = ObitInfoListUnref(out->info);
     out->info = ObitInfoListRef(in->info);
   }
     
   /* copy/set this classes additions */
   /* don't copy ObitTableSel, ObitThread or ObitInfoList */
   err = newObitErr();
   out->myDesc = ObitTableDescCopy(in->myDesc, out->myDesc, err);
   if (err->error) 
     Obit_log_error(err, OBIT_Error, 
		    "%s: Error copying table descriptor for %s", 
		    routine, out->name);
   ObitErrLog (err);
   err = ObitErrUnref(err);
   if (out->tabType) g_free(out->tabType); out->tabType = NULL;
   if (in->tabType) out->tabType = g_strdup(in->tabType);

  return out;
} /* end ObitTableClone */

/**
 * Copy the contents of table in to the end of table out
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
void ObitTableConcat (ObitTable *in, ObitTable *out, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  olong outRec;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nRowPIO;
  gchar *routine = "ObitTableConcat";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Only works one row at a time */
  nRowPIO = 1;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(in->info, "nRowPIO", OBIT_long, dim, &nRowPIO);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitTableOpen (in, OBIT_IO_ReadOnly, err);
  /* if it didn't work bail out */
  if ((iretCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_msg (err, routine, in->name);
  }

 /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (out->buffer) ObitIOFreeBuffer(out->buffer); /* free existing */
  out->buffer     = NULL;
  out->bufferSize = -1;

  /* test open output */
  oretCode = ObitTableOpen (out, OBIT_IO_ReadWrite, err);
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    out->buffer = NULL;
    out->bufferSize = 0;
    Obit_traceback_msg (err, routine, in->name);
  }

  /* Check that the two tables are compatable */
  if (!ObitTableDescCompatible (in->myDesc, out->myDesc)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Tables %s and %s incompatible", 
		   routine, in->name, out->name);
    return;
  }

  /* start writing at end */
  outRec = out->myDesc->nrow + 1;

  /* we're in business, copy */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitTableRead (in, -1, in->buffer, err);
    /* How many */
    out->myDesc->numRowBuff = in->myDesc->numRowBuff;
    if (iretCode!=OBIT_IO_OK) break;
    oretCode = ObitTableWrite (out, outRec, in->buffer, err);
    outRec += out->myDesc->numRowBuff;
  }
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  out->buffer = NULL;
  out->bufferSize = 0;
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  
  /* close files to be sure */
  iretCode = ObitTableClose (in, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  
  /* close files to be sure */
  oretCode = ObitTableClose (out, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, out->name);
  
  return;
} /* end ObitTableConcat */

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
 * \param in     Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *               or OBIT_IO_WriteOnly).
 *               If OBIT_IO_WriteOnly any existing data in the output file
 *               will be lost.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitTableOpen (ObitTable *in, ObitIOAccess access, 
			  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong need;
  gchar *routine = "ObitTableOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Same type of access on descriptor */
  in->myDesc->access = access;

  /* If the file is already open - close it  first */
  if ((in->myStatus==OBIT_Active) || (in->myStatus==OBIT_Modified)) {
    if (in->myIO) retCode = ObitTableClose (in, err);
    else retCode = OBIT_IO_OK;
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* get selection parameters */
  ObitTableGetSelect (in->info, in->mySel, err);
    if (err->error) /* add traceback,return on error */
      Obit_traceback_val (err, routine,in->name, retCode);

  /* create appropriate ObitIO */
  /* unlink any existing IO structure */
  in->myIO = ObitUnref (in->myIO);
  if (in->mySel->FileType==OBIT_IO_FITS) {
    in->myIO = (ObitIO*)newObitIOTableFITS(in->name, in->info, err);
    /* copy selector pointer */
    ((ObitIOTableFITS*)in->myIO)->mySel = 
      ObitUnref(((ObitIOTableFITS*)in->myIO)->mySel);
    ((ObitIOTableFITS*)in->myIO)->mySel = ObitRef(in->mySel);
    /* copy descriptor if write enabled */
    if ((access == OBIT_IO_ReadWrite) || (access == OBIT_IO_WriteOnly))
      ((ObitIOTableFITS*)in->myIO)->myDesc = 
	ObitTableDescCopy(in->myDesc, 
			  ((ObitIOTableFITS*)in->myIO)->myDesc, err);
    else
      ((ObitIOTableAIPS*)in->myIO)->myDesc = newObitTableDesc("descriptor");

  } else if (in->mySel->FileType==OBIT_IO_AIPS) {
    in->myIO = (ObitIO*)newObitIOTableAIPS(in->name, in->info, err);
    /* copy selector pointer */
    ((ObitIOTableAIPS*)in->myIO)->mySel = 
      ObitUnref(((ObitIOTableAIPS*)in->myIO)->mySel);
    ((ObitIOTableAIPS*)in->myIO)->mySel = ObitRef(in->mySel);
    /* copy descriptor if write enabled */
    if ((access == OBIT_IO_ReadWrite) || (access == OBIT_IO_WriteOnly))
      ((ObitIOTableAIPS*)in->myIO)->myDesc = 
	ObitTableDescCopy(in->myDesc, 
			  ((ObitIOTableAIPS*)in->myIO)->myDesc, err);
    else
      ((ObitIOTableAIPS*)in->myIO)->myDesc = newObitTableDesc("descriptor");
  }

  in->myIO->access = access; /* save access type */

 /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOOpen (in->myIO, access, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* read or write Headers */
  if ((access == OBIT_IO_ReadOnly) || (access == OBIT_IO_ReadWrite)) {
    /* read header info */
    retCode = ObitIOReadDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  } 

  if (((access == OBIT_IO_ReadWrite) || (access == OBIT_IO_WriteOnly)) &&
      /* Don't attempt if descriptor not fully defined */
      (in->myIO->myDesc!=NULL) && (((ObitTableDesc*)in->myIO->myDesc)->lrowIO>0)) {
    /* Write header info */
    retCode = ObitIOWriteDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Set descriptors for the output on in to reflect the selection
     by in->mySel,  the descriptors on in->myIO will still describe
     the external representation */
  ObitTableSelSetDesc ((ObitTableDesc*)in->myIO->myDesc,
    (ObitTableSel*)in->myIO->mySel, in->myDesc, err);
  if (err->error) /* add traceback,return on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Allocate buffer - resize if necessary */
  /* buffer size < 0 => no buffer desired */
  if (in->bufferSize >= 0) {
    need = ObitTableSelBufferSize(in->myDesc, in->mySel);
    /* is current one big enough? */
    if ((in->buffer!=NULL) && (need>in->bufferSize)) {
      /* no - deallocate */
      if (in->buffer) ObitIOFreeBuffer(in->buffer);
      in->buffer = NULL;
      in->bufferSize = 0;
    }
    /* new one if needed */
    if (in->buffer==NULL)  
      ObitIOCreateBuffer (&in->buffer, &in->bufferSize, in->myIO, 
			  in->info, err);
  } /* end buffer allocation */

  /* init I/O */
  retCode = ObitIOSet (in->myIO, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Set I/O to beginning of the file */
  ((ObitTableDesc*)in->myIO->myDesc)->firstRow = 0;
  ((ObitTableDesc*)in->myIO->myDesc)->numRowBuff = 0;
  /* For WriteOnly the file is truncated.*/
  if (access == OBIT_IO_WriteOnly) {
    ((ObitTableDesc*)in->myIO->myDesc)->firstRow = 1;
    ((ObitTableDesc*)in->myIO->myDesc)->nrow = 0;
    in->myDesc->nrow = 0;
  }

  /* set Status */
  in->myStatus = OBIT_Active;

  /* save current location */
  in->myDesc->firstRow   = ((ObitTableDesc*)in->myIO->myDesc)->firstRow;
  in->myDesc->numRowBuff = ((ObitTableDesc*)in->myIO->myDesc)->numRowBuff;
  in->myDesc->nrow       = ((ObitTableDesc*)in->myIO->myDesc)->nrow;
  return retCode;
} /* end ObitTableOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitTableClose (ObitTable *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitTableClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  /* Something going on? */
  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;
  if (in->myIO == NULL) return OBIT_IO_OK;

  /* flush buffer if writing */
  if (((in->myIO->access==OBIT_IO_ReadWrite) || 
       (in->myIO->access==OBIT_IO_WriteOnly)) &&
      (in->myStatus == OBIT_Modified)) {
    retCode = ObitIOFlush (in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);

    /* Update descriptor on myIO */
    ObitTableDescCopyDesc(in->myDesc, (ObitTableDesc*)in->myIO->myDesc, err);
    if (err->error)
      Obit_traceback_val (err, routine, in->name, retCode);

    /* Update header on disk if writing */
    retCode = OBIT_IO_OK;
    if ((in->myIO->myStatus != OBIT_Inactive) &&
	/* Don't attempt if descriptor not fully defined */
	(((ObitTableDesc*)in->myIO->myDesc)->lrowIO>0))    
      retCode = ObitIOWriteDescriptor(in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);    
  }

  /* Close actual file */
  retCode = ObitIOClose (in->myIO, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Delete buffer */
  if (in->buffer)  ObitIOFreeBuffer(in->buffer); in->buffer = NULL;
  in->bufferSize = 0; 

  /* set Status */
  in->myStatus = OBIT_Inactive;

  return retCode;
} /* end ObitTableClose */

/**
 * Ensures full instantiation of object - basically open to read/write header
 * and verify or create file.
 * If object has previously been opened, as demonstrated by the existance
 * of its myIO member, this operation is a no-op.
 * Virtual - calls actual class member
 * \param in     Pointer to object
 * \param exist  TRUE if object should previously exist, else FALSE
 * \param err    ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
void ObitTableFullInstantiate (ObitTable *in, gboolean exist, ObitErr *err)
{
  ObitIOAccess access;
  gchar *routine = "ObitTableFullInstantiate";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myIO) return;  /* is this needed */

  /* Open readonly if it should exist, else writeonly */
  if (exist) access = OBIT_IO_ReadOnly;
  else access = OBIT_IO_WriteOnly;

  in->bufferSize = -1;  /* Don't need to assign buffer here */

  /* Open and close */
  ObitTableOpen(in, access, err);
  ObitTableClose(in, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);
  in->bufferSize = 0;  /* May need buffer later */
} /* end ObitTableFullInstantiate */

/**
 * Remove any previously existing rows and fully instantiate.
 * \param in     Pointer to object
 * \param err    ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
void ObitTableClearRows (ObitTable *in, ObitErr *err)
{
  gchar *routine = "ObitTableClearRows";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  in->bufferSize = -1;  /* Don't need to assign buffer here */

  /* Open and close */
  ObitTableOpen(in, OBIT_IO_ReadWrite, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);

  /* reset count to zero */
  in->myDesc->nrow = 0;
  /* The one that counts is in the IO */
  ((ObitTableDesc*)(in->myIO->myDesc))->nrow = 0;
  /* Mark as changed */
  in->myStatus = OBIT_Modified;
  
  ObitTableClose(in, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);
  in->bufferSize = 0;  /* May need buffer later */
} /* end ObitTableClearRows */

/**
 * Read table data from disk.
 * The ObitTableDesc maintains the current location in the table.
 * The number read will be mySel->nRowPIO (until the end of the selected
 * range of rows in which case it will be smaller).
 * The first row number after a read is myDesc->firstRow
 * and the number of row is myDesc->numRowBuff.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param rowno Row number to start reading, -1 = next;
 * \param data pointer to buffer to write results.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitTableRead (ObitTable *in, olong rowno, ofloat *data, 
			  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitTableRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

 /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadOnly;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Read enabled? */
  Obit_retval_if_fail(((in->myIO->access==OBIT_IO_ReadWrite)||
		       (in->myIO->access==OBIT_IO_ReadOnly)), err, retCode,
		      "%s: Read not enabled for %s", routine, in->name);

  /* select internal or external buffer */
  if (myBuf==NULL) {
    myBuf = in->buffer;
    /* Check that internal buffer (defined in bytes) large enough */
    if (in->bufferSize>0) need = in->bufferSize;
    else need = ObitTableSelBufferSize (in->myDesc, in->mySel);
    /* need = ObitTableSelBufferSize (in->myDesc, in->mySel);   DEBUG */
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  g_assert (myBuf != NULL); /* check it */

  retCode = ObitIOReadRow (in->myIO, rowno, myBuf, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location */
  in->myDesc->firstRow   = ((ObitTableDesc*)in->myIO->myDesc)->firstRow;
  in->myDesc->numRowBuff = ((ObitTableDesc*)in->myIO->myDesc)->numRowBuff;

  return retCode;
} /* end ObitTableRead */

/**
 * Read data from disk applying selection.
 * The number read will be mySel->nRowPIO (until the end of the selected
 * range of visibilities in which case it will be smaller).
 * The first visibility number after a read is myDesc->firstRow
 * and the number of visibilities is myDesc->numRowBuff.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in Pointer to object to be read.
 * \param rowno Row number to start reading, -1 = next;
 * \param data pointer to buffer to write results.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitTableReadSelect (ObitTable *in, olong rowno, ofloat *data, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitTableReadSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

 /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Active) && (in->myStatus!=OBIT_Modified)) {
    access = OBIT_IO_ReadOnly;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Read enabled? */
  Obit_retval_if_fail(((in->myIO->access==OBIT_IO_ReadWrite)||
		       (in->myIO->access==OBIT_IO_ReadOnly)), err, retCode,
		      "%s: Read not enabled for %s", routine, in->name);

  /* select internal or external buffer */
  if (myBuf==NULL) {
    myBuf = in->buffer;
    /* Check that internal buffer ( defined in gfloats) large enough */
    need = ObitTableSelBufferSize (in->myDesc, in->mySel);
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  g_assert (myBuf != NULL); /* check it */

  retCode = ObitIOReadRowSelect (in->myIO, rowno, myBuf, err);
  if ((retCode > OBIT_IO_EOF) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* save current location */
  in->myDesc->firstRow   = ((ObitTableDesc*)in->myIO->myDesc)->firstRow;
  in->myDesc->numRowBuff = ((ObitTableDesc*)in->myIO->myDesc)->numRowBuff;

  return retCode;
} /* end ObitTableReadSelect */

/**
 * Write information to disk.
 * The data in the buffer will be written starting at visibility
 * myDesc->firstRow and the number written will be myDesc->numRowBuff
 * which should not exceed mySel->nRowPIO if the internal buffer is used.
 * myDesc->firstRow will be maintained and need not be changed for
 * sequential writing.
 * \param in Pointer to object to be written.
 * \param rowno Row number (1-rel) to start reading, -1 = next;
 * \param data pointer to buffer containing input data.
 *             if NULL, use the buffer member of in.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitTableWrite (ObitTable *in, olong rowno, ofloat *data, 
			   ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOAccess access;
  ofloat *myBuf = data;
  olong need;
  gchar *routine = "ObitTableWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* check and see if its open - if not attempt */
  if ((in->myStatus!=OBIT_Modified) && (in->myStatus!=OBIT_Active)) {
    access = OBIT_IO_WriteOnly;
    retCode = ObitIOOpen (in->myIO, access, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_val (err, routine, in->name, retCode);
  }

  /* Write enabled? */
  Obit_retval_if_fail(((in->myIO->access==OBIT_IO_ReadWrite)||
		       (in->myIO->access==OBIT_IO_WriteOnly)), err, retCode,
		      "%s: Write not enabled for %s", routine, in->name);

  /* select internal or external buffer */
  if (myBuf==NULL) {
    myBuf = in->buffer;
    /* Check that internal buffer ( defined in gfloats) large enough */
    need = ObitTableSelBufferSize (in->myDesc, in->mySel);
    if (need > in->bufferSize) {
      Obit_log_error(err, OBIT_Error, 
		     "IO buffer ( %d) too small, need %d for %s", 
		     in->bufferSize, need, in->name);
      return retCode;
    }
  } 
  g_assert (myBuf != NULL); /* check it */

  /* set number and location to write on myIO descriptor */
  ((ObitTableDesc*)in->myIO->myDesc)->firstRow   = in->myDesc->firstRow;
  ((ObitTableDesc*)in->myIO->myDesc)->numRowBuff = in->myDesc->numRowBuff;

  /* most of the instructions for the I/O are in the ObitInfoList */
  retCode = ObitIOWriteRow (in->myIO, rowno, myBuf, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Modified;

  /* save current location */
  if (rowno>0) {
    if (((ObitTableDesc*)in->myIO->myDesc)->nrow < rowno)
      ((ObitTableDesc*)in->myIO->myDesc)->nrow = rowno;
  } else ((ObitTableDesc*)in->myIO->myDesc)->nrow++;
  in->myDesc->firstRow   = ((ObitTableDesc*)in->myIO->myDesc)->firstRow;
  in->myDesc->numRowBuff = ((ObitTableDesc*)in->myIO->myDesc)->numRowBuff;
  in->myDesc->nrow       = ((ObitTableDesc*)in->myIO->myDesc)->nrow;

  return retCode;
} /* end ObitTableWrite */

/**
 * Convert table to a derived type
 * \param in  Pointer to object to be converted still exists after call.
 * \param err ObitErr for reporting errors.
 * \return converted table
 */
ObitTable* ObitTableConvert (ObitTable *in)
{
  ObitTable *out = NULL;
  const ObitTableClassInfo *myClass;

 /* error checks */
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* This only is defined if this is a derived object */
   myClass = in->ClassInfo;
   if ((gconstpointer)myClass != (gconstpointer)&myClassInfo) {
     
     /* Is function defined in derived class? */
     g_assert (myClass->ObitTableConvert != NULL);

     /* call actual function */
     out = myClass->ObitTableConvert (in);
   }
  
  return out;
} /* end ObitTableConvert */

/**
 * Read one row of table data from disk.
 * The ObitTableDesc maintains the current location in the table.
 * If there are existing rows in the buffer marked as modified 
 * ("_status" column value =1) the buffer is rewritten to disk before 
 * the new buffer is read.
 * \param in    Pointer to object to be read.
 * \param rowno Row number to start reading, -1 = next;
 * \param row   pointer to Table row Structure to accept data.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitTableReadRow (ObitTable *in, olong rowno, ObitTableRow *row,
			  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitTableReadRow";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  g_assert (ObitTableRowIsA(row));

  /* Only one row */
  in->mySel->nRowPIO = 1;

  /* read row rowno */
  retCode = ObitTableRead ((ObitTable*)in, rowno, NULL,  err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* set pointer to buffer */
  row->myRowData = (gpointer)in->buffer;

  return retCode;
} /* end ObitTableReadRow */

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
ObitTableSetRow  (ObitTable *in, ObitTableRow *row,
		  ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(row, &myRowClassInfo));

  if (in->myStatus == OBIT_Inactive) {
    Obit_log_error(err, OBIT_Error,
		   "Table is inactive for  %s ", in->name);
    return;
 }

  /* set pointer to buffer */
  row->myRowData = (gpointer)in->buffer;

} /*  end ObitTableSetRow */

/**
 * Write one row of table data to disk.
 * \param in    Pointer to object to be read.
 * \param rowno Row number to start reading, -1 = next;
 * \param row   pointer to Table row Structure to accept data.
 * \param err   ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitTableWriteRow (ObitTable *in, olong rowno, ObitTableRow *row,
			  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitTableWriteRow";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  g_assert (ObitTableRowIsA(row));

  /* Only one row */
  in->mySel->nRowPIO = 1;
  in->myDesc->numRowBuff = 1;

  /* write row rowno */
  retCode = ObitTableWrite ((ObitTable*)in, rowno, (ofloat*)row->myRowData,  
			    err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

   in->myStatus = OBIT_Modified;  /* Modified */

  /* Mark as unsorted */
  in->myDesc->sort[0] = 0;
  in->myDesc->sort[1] = 0;
  ((ObitTableDesc*)in->myIO->myDesc)->sort[0] = 0;
  ((ObitTableDesc*)in->myIO->myDesc)->sort[1] = 0;

  return retCode;
} /* end ObitTableWriteRow */

/**
 * Determine table type (name of type, e.g. "AIPS AN")
 * \param in    Pointer to object to be read.
 * \param err   ObitErr for reporting errors.
 * \return pointer to string
 */
gchar* ObitTableGetType (ObitTable *in, ObitErr *err)
{
  return in->tabType;
} /* end ObitTableGetType */

/**
 * Determine table version number (1-rel)
 * \param in    Pointer to object to be read.
 * \param err   ObitErr for reporting errors.
 * \return version number
 */
olong ObitTableGetVersion (ObitTable *in, ObitErr *err)
{
  return in->tabVer;
} /* end ObitTableGetVersion */

/**
 * Get underlying file information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
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
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFileName OBIT_string  FITS file name
 * \li xxxDisk     OBIT_oint    FITS file disk number
 * \li xxxDir      OBIT_string  Directory name for xxxDisk
 *
 * For all File types types:
 * \li xxxDataType OBIT_string "UV" = UV data, "MA"=>image, "Table"=Table, 
 *                "OTF"=OTF, etc
 * \li xxxFileType OBIT_string "AIPS", "FITS"
 *    
 * For xxxDataType = "Table"
 * \li xxxTab   OBIT_string  (Tables only) Table type (e.g. "AIPS CC")
 * \li xxxVer   OBIT_oint    (Tables Only) Table version number
 *    
 * \param err     ObitErr for reporting errors.
 */
void ObitTableGetFileInfo (ObitTable *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err)
{
  const ObitIOClassInfo *myIOClass;
  gchar *routine = "ObitTableGetFileInfo";

  /* Fully instantiate */
  ObitTableFullInstantiate(in, TRUE, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);

  /* use class function on myIo */
  myIOClass = in->myIO->ClassInfo;
  myIOClass->ObitIOGetFileInfo (in->myIO, in->info, prefix, outList, err);

} /* end ObitTableGetFileInfo */

/**
 * Create a data object with selection parameters set from an InfoList
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param inList InfoList to extract object information from
 * Following InfoList entries for AIPS files ("xxx" = prefix):
 * \li xxxName  OBIT_string  AIPS file name
 * \li xxxClass OBIT_string  AIPS file class
 * \li xxxDisk  OBIT_oint    AIPS file disk number
 * \li xxxSeq   OBIT_oint    AIPS file Sequence number
 * \li AIPSuser OBIT_oint    AIPS User number
 * \li xxxCNO   OBIT_oint    AIPS Catalog slot number
 * \li xxxDir   OBIT_string  Directory name for xxxDisk
 *
 * Following entries for FITS files ("xxx" = prefix):
 * \li xxxFileName OBIT_string  FITS file name
 * \li xxxDisk     OBIT_oint    FITS file disk number
 * \li xxxDir      OBIT_string  Directory name for xxxDisk
 *
 * For all File types:
 * \li xxxFileType OBIT_string "AIPS", "FITS"
 *    
 * For xxxDataType = "Table"
 * \li xxxTableParent OBIT_string  (Tables only) Table parent type (e.g. "MA")
 * \li xxxTab   OBIT_string  (Tables only) Table type (e.g. "AIPS CC")
 * \li xxxVer   OBIT_oint    (Tables Only) Table version number
 *    
 * \param err     ObitErr for reporting errors.
 * \return new data object with selection parameters set
 */
ObitTable* ObitTableFromFileInfo (gchar *prefix, ObitInfoList *inList, 
				  ObitErr *err)
{
  ObitTable *out = NULL;
  ObitData  *parent = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong tabVer;
  gchar *keyword=NULL, *tabType=NULL;
  gchar *routine = "ObitTableFromFileInfo";

  /* GetParent File */ 
  parent = ObitDataFromFileInfo(prefix, inList, err);
  if (err->error) Obit_traceback_val (err, routine, routine, out);

  /* Get basic information - table type */
  if (prefix) keyword =  g_strconcat (prefix, "Tab", NULL);
  else keyword =  g_strconcat ("Tab", NULL);
  if (!ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&tabType)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry %s not in InfoList", routine, keyword);
    g_free(keyword);
    goto finish;
  }
  g_free(keyword);

  /* Table version number */
  if (prefix) keyword =  g_strconcat (prefix, "Ver", NULL);
  else keyword =  g_strconcat ("Ver", NULL);
  tabVer = 0;
  if (!ObitInfoListGetTest(inList, keyword, &type, dim, &tabVer)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry %s not in InfoList", routine, keyword);
    g_free(keyword);
    goto finish;
  }
  g_free(keyword);

  /* Get Table */
  if (parent)
    out = newObitDataTable (parent, OBIT_IO_ReadWrite, tabType, &tabVer, err);
 finish:
  parent = ObitDataUnref(parent);
  if (err->error) Obit_traceback_val (err, routine, routine, out);

  return out;
} /* end ObitTableFromFileInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTableClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableClassInfoDefFn (gpointer inClass)
{
  ObitTableClassInfo *theClass = (ObitTableClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableGetClass;
  theClass->newObit       = (newObitFP)newObitTable;
  theClass->ObitTableFromFileInfo = (ObitTableFromFileInfoFP)ObitTableFromFileInfo;
  theClass->ObitTableZap  = (ObitTableZapFP)ObitTableZap;
  theClass->ObitCopy      = (ObitCopyFP)ObitTableCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitTableClone;
  theClass->ObitClear     = (ObitClearFP)ObitTableClear;
  theClass->ObitInit      = (ObitInitFP)ObitTableInit;
  theClass->ObitTableConvert = (ObitTableConvertFP)ObitTableConvert;
  theClass->ObitTableOpen    = (ObitTableOpenFP)ObitTableOpen;
  theClass->ObitTableClose   = (ObitTableCloseFP)ObitTableClose;
  theClass->ObitTableRead    = (ObitTableReadFP)ObitTableRead;
  theClass->ObitTableClearRows= (ObitTableClearRowsFP)ObitTableClearRows;
  theClass->ObitTableFullInstantiate = 
    (ObitTableFullInstantiateFP)ObitTableFullInstantiate;
  theClass->ObitTableReadSelect = 
    (ObitTableReadSelectFP)ObitTableReadSelect;
  theClass->ObitTableWrite = 
    (ObitTableWriteFP)ObitTableWrite;
  theClass->ObitTableReadRow = 
    (ObitTableReadRowFP)ObitTableReadRow;
  theClass->ObitTableWriteRow = 
    (ObitTableWriteRowFP)ObitTableWriteRow;
  theClass->ObitTableSetRow = 
    (ObitTableSetRowFP)ObitTableSetRow;
  theClass->ObitTableGetType = 
    (ObitTableGetTypeFP)ObitTableGetType;
  theClass->ObitTableGetVersion = 
    (ObitTableGetVersionFP)ObitTableGetVersion;
  theClass->ObitTableGetFileInfo = 
    (ObitTableGetFileInfoFP)ObitTableGetFileInfo;

} /* end ObitTableClassDefFn */

/*---------------Private functions--------------------------*/
/*----------------  Table Row  ----------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTableRowInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableRow *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myRowClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myTable   = NULL;
  in->myRowData = NULL;

} /* end ObitTableRowInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitTableRow* cast to an Obit*.
 */
void ObitTableRowClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableRow *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myRowClassInfo));

  /* delete this class members */
  in->myTable = ObitTableUnref(in->myTable);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTableRowClear */


/*------------------  Table  ------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTableInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTable *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList(); 
  in->myIO      = NULL;
  in->myDesc    = newObitTableDesc(in->name);
  in->mySel     = newObitTableSel(in->name);
  in->myStatus  = OBIT_Inactive;
  in->buffer    = NULL;
  in->bufferSize= 0;
  in->tabVer    = -1;
  in->tabType   = NULL;
  in->myHost    = NULL;

} /* end ObitTableInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitTable* cast to an Obit*.
 */
void ObitTableClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTable *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->info   = ObitInfoListUnref(in->info);
  in->thread = ObitThreadUnref(in->thread);
  in->myIO   = ObitUnref(in->myIO);
  in->myDesc = ObitTableDescUnref(in->myDesc);
  in->mySel  = ObitTableSelUnref(in->mySel);
  /* myHost a secret reference don't unreference */ 
  if (in->buffer)  ObitIOFreeBuffer(in->buffer); 
  if (in->tabType) g_free(in->tabType); in->tabType = NULL;
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTableClear */

/**
 * Get requested information from the ObitInfoList
 * \param info Pointer to InfoList
 * \param sel  pointer to uvdata selector to update.
 * \param err  ObitErr for reporting errors.
 */
static void ObitTableGetSelect (ObitInfoList *info, ObitTableSel *sel,
				ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(info));
  g_assert (ObitIsA(sel, ObitTableSelGetClass()));

  /* what type of underlying file? */
  if (!ObitInfoListGet(info, "FileType", &type, dim, 
		       (gpointer)&sel->FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		"ObitTableGetSelect: entry FileType not in InfoList Object %s",
		sel->name);
  }

  /* Maximum number of rows per read/write? [default 1] */
  sel->nRowPIO = 1;
  ObitInfoListGetTest(info, "nRowPIO", &type, (gint32*)dim, &sel->nRowPIO);
  sel->nRowPIO = MAX (1, sel->nRowPIO); /* no fewer than 1 */


} /* end ObitTableGetSelect */


