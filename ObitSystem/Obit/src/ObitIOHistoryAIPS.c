/* $Id: ObitIOHistoryAIPS.c,v 1.10 2006/06/26 16:48:22 bcotton Exp $  */
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

#include "ObitIOHistoryAIPS.h"
#include "ObitAIPS.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOHistoryAIPS.c
 * ObitIOHistoryAIPS class function definitions.
 *
 * This class allows history access to AIPS files
 * This class is derived from ObitIOHistory
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOHistoryAIPS";

/** Function to obtain parent ClassInfo - ObitIOHistory */
static ObitGetClassFP ObitParentGetClass = ObitIOHistoryGetClass;

/** Number of AIPS history records ber AIPS block */
#define AIPSHistPerBlock 14;

/**
 * ClassInfo global structure ObitIOHistoryAIPSClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOHistoryAIPSClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOHistoryAIPSInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOHistoryAIPSClear (gpointer in);

/** Private: Determine HI filename. */
static gchar* GetFileName (ObitIOHistoryAIPS *in, ObitErr *err);

/** Private: Read AIPS Block. */
static ObitIOCode ReadBlock (ObitIOHistoryAIPS *in, olong block, ObitErr *err);

/** Private: Write AIPS Block currently in buffer. */
static ObitIOCode WriteBlock (ObitIOHistoryAIPS *in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitIOHistoryAIPSClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name Name [optional] for object
 * \param info InfoList defining file, must fully define a file (disk, CNU, user)
 * \param err ObitErr for reporting errors. 
 * \return the new object.
 */
ObitIOHistoryAIPS* newObitIOHistoryAIPS (gchar* name, ObitInfoList *info,
					 ObitErr *err)
{
  ObitIOHistoryAIPS* out;
  gchar *routine = "newObitIOHistoryAIPS";
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOHistoryAIPSClassInit();
  
  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitIOHistoryAIPS));
  
  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");
  
  /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;
  
  /* initialize other stuff */
  ObitIOHistoryAIPSInit((gpointer)out);
  
  /* Copy ObitInfoList to get relevant information */
  out->info = ObitInfoListUnref(out->info); /* out with the old - empty */
  out->info = ObitInfoListCopy(info);       /* in with the new - copy */
  
  /* Basic file IO object */
  out->myFile = newObitFile(out->name);
  if (out->FileName) g_free(out->FileName); /* free old */
  out->FileName = GetFileName(out, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);
  
 return out;
} /* end newObitIOHistoryAIPS */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOHistoryAIPSGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOHistoryAIPSClassInit();
  
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
gboolean ObitIOHistoryAIPSSame (ObitIO *in, ObitInfoList *in1, 
				ObitInfoList *in2, ObitErr *err)
{
  olong CNO1, UserId1, disk1, CNO2, UserId2, disk2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean same = FALSE;
  gchar *routine = " ObitIOHistoryAIPSSame";

  /* error checks */
  if (err->error) return same;
  g_assert (ObitIOHistoryAIPSIsA(in));

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
} /* end ObitIOHistoryAIPSSame */

/**
 * Delete underlying files.
 * Note: this does not delete the object
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors. 
 */
void ObitIOHistoryAIPSZap (ObitIOHistoryAIPS *in, ObitErr *err)
{
  gchar *routine = " ObitIOHistoryAIPSZap";

  /* Don't bother if NULL */
  if (!in) return;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Zap myFile */
  in->myFile = ObitFileZap (in->myFile, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  return;
} /* end ObitIOHistoryAIPSZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object. 
 * \return pointer to the new object.
 */
ObitIOHistoryAIPS* ObitIOHistoryAIPSCopy  (ObitIOHistoryAIPS *in, ObitIOHistoryAIPS *out, ObitErr *err)
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
    out = newObitIOHistoryAIPS(outName, NULL, err);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->access = OBIT_IO_None; /* not currently reading file */
  if (out->FileName) g_free(out->FileName); /* free old */
  out->FileName = g_strdup(in->FileName);

  return out;
} /* end ObitIOHistoryAIPSCopy */

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
ObitIOCode ObitIOHistoryAIPSOpen (ObitIOHistoryAIPS *in, ObitIOAccess access, 
				  ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong size;
  gchar *routine = "ObitIOHistoryAIPSOpen";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));

  /* Get current filename */
  if (in->FileName) g_free(in->FileName); /* free old */
  in->FileName = GetFileName(in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Open */
  size = 256 * sizeof(AIPSint);  /* AIPS blocking */
  retCode = ObitFileOpen (in->myFile, in->FileName, access, OBIT_IO_Binary, 
			  size, err);
  if (retCode == OBIT_IO_EOF) return retCode;
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Save position, init other info */
  in->filePos  = 0;
  in->CurrentBlock = -1;
  in->dirty = FALSE;
  in->access = access; /* just in case */

  return retCode;
} /* end ObitIOHistoryAIPSOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryAIPSClose (ObitIOHistoryAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOHistoryAIPSClose";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Flush if needed */
  if (in->dirty) retCode = WriteBlock(in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  retCode = ObitFileClose (in->myFile, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Save position, init other info */
  in->filePos = 0;
  in->CurrentBlock = -1;

  return OBIT_IO_OK;
} /* end ObitIOHistoryAIPSClose */

/**
 * initialize I/O
 * Not needed
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryAIPSSet (ObitIOHistoryAIPS *in, ObitInfoList *info, ObitErr *err)
{
  /* Not needed */
  return OBIT_IO_OK;
} /* end ObitIOHistoryAIPSSet */

/**
 * Read specified History record
 * NB: there is something of a mismatch in the AIPS size of a HI record and Obit
 * only the first 70 characters are kept
 * \param in     Pointer to object to be read.
 * \param recno  record number (1-rel) -1=> next.
 * \param hiCard output history record (70 char+NULL)
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryAIPSReadRec (ObitIOHistoryAIPS *in, olong recno, 
				     gchar *hiCard, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, record, block, number, offset, NL;
  gchar *routine = "ObitIOHistoryAIPSReadRec";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (hiCard != NULL);

  /* which record? */
  if (recno>0) record = recno;
  else record = in->Last+1;    /* default is next */
  in->Last = record;           /* remember */

  /* Too far? */
  if (record>in->CurrentNumber) return OBIT_IO_EOF;

  /* which block? */
  NL = AIPSHistPerBlock;  /* number of lines per block */
  block = 1 + (record-1) / NL;

  /* Need new block? */
  if (block!=in->CurrentBlock) retCode = ReadBlock(in, block, err);
  if (retCode == OBIT_IO_EOF) return retCode;
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Where in record? */
  number = record - (block-1)*NL;
  offset = 4 * sizeof(AIPSint) + (number-1)*72;

  /* extract 70 characters */
  for (i=0; i<70; i++) hiCard[i] = in->buffer[offset+i]; hiCard[i] = 0;

  return OBIT_IO_OK;
} /* end ObitIOHistoryAIPSReadRec */

/**
 * Write specified History record
 * \param in     Pointer to object to be written.
 * \param recno  Record number (1-rel) -1=> next, overwrites any existing
 * \param hiCard input history record (70 char)
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryAIPSWriteRec (ObitIOHistoryAIPS *in, olong recno, gchar *hiCard, 
				      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, record, block, number, offset, NL;
  gchar *routine = "ObitIOHistoryAIPSWriteRec";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (hiCard != NULL);

  /* which record? */
  if (recno>0) record = recno;
  else record = in->CurrentNumber+1;    /* default is after last */
  in->Last = record;  /* remember */
  in->CurrentNumber = MAX (record, in->CurrentNumber);  /* keep track */

  /* which block? */
  NL = AIPSHistPerBlock;  /* number of lines per block */
  block = 1 + (record-1) / NL;

  /* Need new block? - rewrite old if modified */
  if (in->dirty && (block!=in->CurrentBlock)) retCode = WriteBlock(in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  /* Read current block if needed */
  if (block!=in->CurrentBlock) retCode = ReadBlock(in, block, err);
  in->CurrentBlock = block;
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  in->MaxNumber    = MAX (in->MaxNumber, block*NL);

  /* Where in record? */
  number = record - (block-1)*NL;
  offset = 4 * sizeof(AIPSint) + (number-1)*72;

  /* save 70 characters */
  for (i=0; i<70; i++) in->buffer[offset+i] =  hiCard[i];
  in->buffer[offset+i] = ' ';  in->buffer[offset+i+1] = ' '; /* Blank fill to AIPS size */
  in->dirty = TRUE;

  return OBIT_IO_OK;
} /* end ObitIOHistoryAIPSWriteRec */

/**
 * Tell number of History records
 * \param in     Pointer to open object to be tested
 * \return  number of records, <0 => problem
 */
olong ObitIOHistoryAIPSNumRec (ObitIOHistory *in)
{
  olong out;

  out = in->CurrentNumber;

  return out;
} /* end ObitIOHistoryAIPSNumRec */

/**
 * Read image Descriptor data from disk.
 * \param in Pointer to object to be read.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryAIPSReadDescriptor (ObitIOHistoryAIPS *in,  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong block;
  AIPSint *AI;
  gchar *routine = "ObitIOHistoryAIPSReadDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* read the first block and extract info */
  block = 1;
  if (block!=in->CurrentBlock) retCode = ReadBlock(in, block, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* extract */
  AI = (AIPSint*)in->buffer;
  in->CurrentNumber = AI[0];
  in->MaxNumber     = AI[1];

  return OBIT_IO_OK;
} /* end ObitIOHistoryAIPSReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in  Pointer to object whose header is to be updated..
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryAIPSWriteDescriptor (ObitIOHistoryAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong block;
  AIPSint *AI;
  gchar *routine = "ObitIOHistoryAIPSWriteDescriptor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* read the first block and set info */
  block = 1;
  if (block!=in->CurrentBlock) retCode = ReadBlock(in, block, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* set values */
  AI = (AIPSint*)in->buffer;
  AI[0] = in->CurrentNumber;
  AI[1] = in->MaxNumber;
  in->dirty = TRUE;

  /* rewrite block */
  in->CurrentBlock = 1;
  retCode = WriteBlock(in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return OBIT_IO_OK;
} /* end ObitIOHistoryAIPSWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryAIPSFlush (ObitIOHistoryAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitIOHistoryAIPSFlush";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Flush if needed */
  if (in->dirty) retCode = WriteBlock(in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return OBIT_IO_OK;
} /* end ObitIOHistoryAIPSFlush */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOHistoryAIPSClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOHistoryAIPSClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOHistoryAIPSClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOHistoryAIPSClassInfoDefFn (gpointer inClass)
{
  ObitIOHistoryAIPSClassInfo *theClass = (ObitIOHistoryAIPSClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOHistoryAIPSClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOHistoryAIPSClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOHistoryAIPSGetClass;
  theClass->newObit       = NULL;
  theClass->newObitIO     = (newObitIOFP)newObitIOHistoryAIPS;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOHistoryAIPSCopy;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOHistoryAIPSSame;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOHistoryAIPSZap;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitIOHistoryAIPSClear;
  theClass->ObitInit      = (ObitInitFP)ObitIOHistoryAIPSInit;
  theClass->ObitIOOpen    = (ObitIOOpenFP)ObitIOHistoryAIPSOpen;
  theClass->ObitIOClose   = (ObitIOCloseFP)ObitIOHistoryAIPSClose;
  theClass->ObitIOSet     = (ObitIOSetFP)ObitIOHistoryAIPSSet;
  theClass->ObitIORead          = NULL;
  theClass->ObitIOReadSelect    = NULL;
  theClass->ObitIOReadRowSelect = NULL;
  theClass->ObitIOWriteRow      = NULL;
  theClass->ObitIOWrite         = NULL;
  theClass->ObitIOFlush         = (ObitIOFlushFP)ObitIOHistoryAIPSFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOHistoryAIPSReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOHistoryAIPSWriteDescriptor;
  theClass->ObitIOCreateBuffer = NULL;
  theClass->ObitIOFreeBuffer   = NULL;
  theClass->newObitIOTable     = NULL;
  theClass->ObitIOUpdateTables = NULL;
  theClass->ObitIOHistoryReadRec  = 
    (ObitIOHistoryReadRecFP)ObitIOHistoryAIPSReadRec;
  theClass->ObitIOHistoryWriteRec = 
    (ObitIOHistoryWriteRecFP)ObitIOHistoryAIPSWriteRec;
  theClass->ObitIOHistoryNumRec  = 
    (ObitIOHistoryNumRecFP)ObitIOHistoryAIPSNumRec;

} /* end ObitIOHistoryAIPSClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOHistoryAIPSInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIOHistoryAIPS *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myFile       = NULL;
  in->FileName     = NULL;
  in->CurrentBlock = 0;
  in->dirty        = FALSE;

 } /* end ObitIOHistoryAIPSInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOHistoryAIPSClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIOHistoryAIPS *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  in->myFile = ObitFileUnref(in->myFile);
  if (in->FileName) g_free(in->FileName); in->FileName = NULL;
  
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOHistoryAIPSClear */

/**
 * Determine the full name of the HI file to be operated on
 * \param in       Object whose HI file name is desired
 * \param err      ObitErr for reporting errors. 
 * \return name string, deallocate when finished.
 */
static gchar* GetFileName (ObitIOHistoryAIPS *in, ObitErr *err)
{
  gchar *out = NULL;
  olong CNO, UserId, disk;
  gint32 dim[IM_MAXDIM];
  ObitInfoType type;
  gchar *routine = "GetFileName";

  /* get instructions from info */
  if(!ObitInfoListGet(in->info, "Disk", &type, dim, &disk, err))
    Obit_traceback_val (err, routine, in->name, out);

  if(!ObitInfoListGet(in->info, "User", &type, dim, &UserId, err))
    Obit_traceback_val (err, routine, in->name, out);

  if(!ObitInfoListGet(in->info, "CNO", &type, dim, &CNO, err))
    Obit_traceback_val (err, routine, in->name, out);

  /* form file name for image file */
  out = ObitAIPSFilename (OBIT_AIPS_History, disk, CNO, UserId, "HI", 1, err);
  return out;

} /* end GetFileName */

/**
 * Read specified block into buffer
 * File size extended if necessary and new blocks blank filled.
 * Previous block written if modified.
 * \param in     Pointer to object to be written.
 * \param block  block number (1-rel) of AIPS HI file.
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
static ObitIOCode ReadBlock (ObitIOHistoryAIPS *in, olong block, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, size, nblock;
  gchar *routine = "ReadBlock";

/* check is existing block is dirty, if so write firrst,
  if block is beyond current EOF, extend file and blank fill, 
  update current block
  update max number of entries */

  /* Save the old? */
  if (in->dirty) retCode = WriteBlock(in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* Is this past the end of the current file? */
  size = 256 * sizeof(AIPSint);  /* AIPS blocking */
  nblock = in->MaxNumber/AIPSHistPerBlock;
  nblock = MAX (1, nblock);
   if (block>nblock) {
    /* Better only be by 1 */
    if (block>(nblock+1)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Requesting extention of HI file by > 1 block %d for %s", 
		     routine, block-nblock,in->name);
      return OBIT_IO_SpecErr;
   }
    /* Extend file - blank buffer */
    for (i=0; i<size; i++) in->buffer[i] = ' ';

    /* new block */
    in->CurrentBlock = block;
    in->filePos = size * (in->CurrentBlock-1);
    retCode = ObitFileWrite(in->myFile, in->filePos, size, in->buffer, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    /* new maximum number of records */
    in->MaxNumber = in->CurrentBlock * AIPSHistPerBlock;
  }

  size = 256 * sizeof(AIPSint);  /* AIPS blocking */
  in->filePos = size * (block-1);
  retCode = ObitFileRead(in->myFile, in->filePos, size, in->buffer, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  in->dirty = FALSE;
  in->CurrentBlock = block;
  return retCode;
} /* end ReadBlock */

/**
 * Write current block in buffer 
 * Writes header info into block if it's no. 1
 * \param in     Pointer to object to be written.
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
static ObitIOCode WriteBlock (ObitIOHistoryAIPS *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong size;
  AIPSint *AI;
  gchar *routine = "WriteBlock";

  /* Just to be sure, write header info in first block */
  if (in->CurrentBlock==1) {
    AI = (AIPSint*)in->buffer;
    AI[0] = in->CurrentNumber;
    AI[1] = in->MaxNumber;
  }

  size = 256 * sizeof(AIPSint);  /* AIPS blocking */
  in->filePos = size * (in->CurrentBlock-1);
  retCode = ObitFileWrite(in->myFile, in->filePos, size, in->buffer, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  in->dirty = FALSE;
  return retCode;
} /* end WriteBlock */

