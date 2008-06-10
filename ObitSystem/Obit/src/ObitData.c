/* $Id$          */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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

#include "ObitData.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitIOHistoryAIPS.h"
#include "ObitImage.h"
#include "ObitIOTableAIPSUtil.h"
#include "ObitIOTableFITSUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitData.c
 * ObitData class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitData";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitDataClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDataClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDataInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDataClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDataClassInfoDefFn (gpointer inClass);

/** Private: Set Class function pointers. */
static void ObitDataClassInfoDefFn (gpointer inClass);

/** Private: Determine file type. */
static ObitIOType getIOType (ObitInfoList *list, ObitErr *err);

/** 
 * Macro to determine if an object is a generic ObitData
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGenericData(in) in->ClassInfo==ObitDataGetClass()

/*----------------------Public functions---------------------------*/
/**
 * Constructor.  
 * A generic ObitData object allows access to tables byt not the main data
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitData* newObitData (gchar* name)
{
  ObitData* out = NULL;
  /*gchar *routine = "newObitData";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDataClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitData));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDataInit((gpointer)out);

 return out;
} /* end newObitData */

/**
 * Create a scratch file suitable for accepting the data to be read from in.
 * A scratch Data is more or less the same as a normal Data except that it is
 * automatically deleted on the final unreference.
 * The output will have the underlying files of the same type as in already 
 * allocated.
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in  The object to copy
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitData* newObitDataScratch (ObitData *in, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  gchar *routine = "newObitDataScratch";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return NULL;
  } else { /* use class function */
    myClass = (ObitDataClassInfo*)in->ClassInfo;
    return myClass->newObitDataScratch (in, err);
  }

} /* end newObitDataScratch */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDataGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDataClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDataGetClass */

/**
 * Test if two ObitDatas have the same underlying structures.
 * This test is done using values entered into the #ObitInfoList
 * in case the object has not yet been opened.
 * Not supported for Generic ObitData
 * \param in1 First object to compare
 * \param in2 Second object to compare
 * \param err ObitErr for reporting errors.
 * \return TRUE if to objects have the same underlying structures
 * else FALSE
 */
gboolean ObitDataSame (ObitData *in1, ObitData *in2, ObitErr *err )
{
  gboolean same = FALSE;
  gchar *routine = "ObitDataSame";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return same;
  g_assert (ObitIsA(in1, &myClassInfo));
  g_assert (ObitIsA(in2, &myClassInfo));

  /* If pointers are the same they must be the same */
  if (in1 == in2) return TRUE;

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in1) || ObitGenericData(in2)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in1->name);
    return FALSE;
  } else { /* use class function */
    /* Get IO if needed */
    if (!in1->myIO) ObitDataSetupIO (in1, err);
    if (err->error) Obit_traceback_val (err, routine, in1->name, FALSE);
    if (!in2->myIO) ObitDataSetupIO (in2, err);
    if (err->error) Obit_traceback_val (err, routine, in2->name, FALSE);
    
    /* Ask IO */
    same = ObitIOSame(in1->myIO, in1->info, in2->info, err);
    if (err->error) Obit_traceback_val (err, routine, in1->name, FALSE);
  }

  return same;
} /* end ObitDataSame */

/**
 * Delete underlying files and the basic object.
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 * \return pointer for input object, NULL if deletion successful
 */
ObitData* ObitDataZap (ObitData *in, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  ObitHistory *inHistory=NULL;
  gchar *routine = "ObitDataZap";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Delete history */
  inHistory  = newObitDataHistory (in, OBIT_IO_WriteOnly, err);
  if (inHistory)  inHistory =  ObitHistoryZap(inHistory, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, in);
 

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return NULL;
  } else { /* use class function */
    myClass = in->ClassInfo;
    return myClass->ObitDataZap (in, err);
  }
} /* end ObitDataZap */

/**
 * Rename underlying files.
 * New name information depends on the underlying file type and is
 * given on the info member.
 * Not supported for Generic ObitData
 * For FITS files:
 * \li "newFileName" OBIT_string (?,1,1) New Name of disk file.
 *
 * For AIPS:
 * \li "newName" OBIT_string (12,1,1) New AIPS Name 
 *      Blank = don't change
 * \li "newClass" OBIT_string (6,1,1) New AIPS Class
 *      Blank = don't changeO
 * \li "newSeq" OBIT_int (1,1,1) New AIPS Sequence
 *      0 => unique value
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors.
 */
void ObitDataRename (ObitData *in, ObitErr *err)
{
  const ObitIOClassInfo *myIOClass;
  gchar *routine = "ObitDataRename";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return;
  } else { /* use class function */
    myIOClass = in->myIO->ClassInfo;
    myIOClass->ObitIORename (in->myIO, in->info, err);
  }
} /* end ObitDataRename */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Objects should be closed on input and will be closed on output.
 * In order for the disk file structures to be copied, the output file
 * must be sufficiently defined that it can be written.
 * The copy will be attempted but no errors will be logged until
 * both input and output have been successfully opened.
 * If the contents of the data are copied, all associated tables are 
 * copied first.
 * ObitInfoList and ObitThread members are only copied if the output object
 * didn't previously exist.
 * Parent class members are included but any derived class info is ignored.
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitData* ObitDataCopy (ObitData *in, ObitData *out, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataCopy";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return NULL;
  } else { /* use class function */
    myClass = in->ClassInfo;
    return myClass->ObitDataCopy (in, out, err);
  }

} /* end ObitDataCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create a data object similar to the input one.
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 */
void ObitDataClone  (ObitData *in, ObitData *out, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataClone";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return;
  } else { /* use class function */
    myClass = in->ClassInfo;
    myClass->ObitDataClone (in, out, err);
  }

} /* end ObitDataClone */

/**
 * Initialize structures and open file.
 * Virtual - calls actual class member
 * Reads table list if in generic ObitData Object
 * \param in Pointer to object to be opened.
 * \param access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *               OBIT_IO_ReadCal or OBIT_IO_WriteOnly).
 *               If OBIT_IO_WriteOnly any existing data in the output file
 *               will be lost.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitDataOpen (ObitData *in, ObitIOAccess access, ObitErr *err)
{
  ObitIOType FileType;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataOpen";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* For generic ObitData read table list*/
  if (ObitGenericData(in)) {
    FileType = getIOType (in->info, err);
    if (FileType==OBIT_IO_FITS) {
      ObitIOTableFITSUtilReadTableList(in, err);
    } else if (FileType==OBIT_IO_AIPS) {
      ObitIOTableAIPSUtilReadTableList(in, err);
    } else { /* Unknown type */
      Obit_log_error(err, OBIT_Error, 
		     "%s: unknown file type %d for %s", 
		     routine, FileType, in->name);
    }
    if (err->error) retCode = OBIT_IO_ReadErr;
    else retCode = OBIT_IO_OK;
    
  } else { /* use class function */
    myClass = in->ClassInfo;
    retCode = myClass->ObitDataOpen (in, access, err);
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Active;

  /* set aaccess */
  in->DataAccess = access;

  return retCode;
} /* end ObitDataOpen */

/**
 * Shutdown I/O.
 * Virtual - calls actual class member
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitDataClose (ObitData *in, ObitErr *err)
{
  ObitIOType FileType;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataClose";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* For generic ObitData update table list if write enabled */
  if (ObitGenericData(in)) {
    if (((in->DataAccess==OBIT_IO_ReadWrite) || 
	(in->DataAccess==OBIT_IO_WriteOnly)) && 
	(in->myStatus==OBIT_Modified)) {
    FileType = getIOType (in->info, err);
    if (FileType==OBIT_IO_FITS) {
      ObitIOTableFITSUtilWriteTableList(in, err);
    } else if (FileType==OBIT_IO_AIPS) {
      ObitIOTableAIPSUtilWriteTableList(in, err);
    } else { /* Unknown type */
      Obit_log_error(err, OBIT_Error, 
		     "%s: unknown file type %d for %s", 
		     routine, FileType, in->name);
    }
    if (err->error) retCode = OBIT_IO_ReadErr;
    else retCode = OBIT_IO_OK;
    }

  } else { /* use class function */
    myClass = in->ClassInfo;
    retCode = myClass->ObitDataClose (in, err);
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* set Status */
  in->myStatus = OBIT_Inactive;

  return retCode;
} /* end ObitDataClose */

/**
 * Ensures full instantiation of object - basically open to read/write header
 * and verify or create file.
 * If object has previously been opened, as demonstrated by the existance
 * of its myIO member, this operation is a no-op.
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in     Pointer to object
 * \param exist  TRUE if object should previously exist, else FALSE
 * \param err    ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
void ObitDataFullInstantiate (ObitData *in, gboolean exist, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataFullInstantiate ";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->myIO) return;  /* is this needed? */

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return;
  } else { /* use class function */
    myClass = in->ClassInfo;
    myClass->ObitDataFullInstantiate (in, exist, err);
  }
} /* end ObitDataFullInstantiate */

/**
 * Return a ObitTable Object to a specified table associated with
 * the input ObitData.  
 * If such an object exists, a reference to it is returned,
 * else a new object is created and entered in the ObitTableList.
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
 *                 and the highest+1 for write.
 * \param err      ObitErr for reporting errors.
 * \return pointer to created ObitTable, NULL on failure.
 */
ObitTable* 
newObitDataTable (ObitData *in, ObitIOAccess access, 
		  gchar *tabType, olong *tabVer, ObitErr *err)
{
  ObitTable *out=NULL;
  ObitIOType FileType;
  gboolean doClose, isImage;
  ObitIOCode retCode;
  gchar *routine = "newObitDataTable";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  g_assert(tabType!=NULL);
  g_assert(tabVer!=NULL);

  /* Trap Generic ObitData */
  if (ObitGenericData(in)) {
    FileType = getIOType (in->info, err);
    if (FileType==OBIT_IO_FITS) {
      out = ObitIOTableFITSUtilGetTable(in, access, tabType, tabVer, err);
    } else if (FileType==OBIT_IO_AIPS) {
      out = ObitIOTableAIPSUtilGetTable(in, access, tabType, tabVer, err);
    } else { /* Unknown type */
      Obit_log_error(err, OBIT_Error, 
		     "%s: unknown file type %d for %s", 
		     routine, FileType, in->name);
    }
    if (err->error) retCode = OBIT_IO_ReadErr;
    else retCode = OBIT_IO_OK;
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    /* end Generic ObitData */

  } else {  
    /* Specific ObitData */
    /* the IO object must be present - create if necessary */
    if (in->myIO==NULL) {
      ObitDataSetupIO (in, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
      
      /* add table list reference */
      in->myIO->tableList = (Obit*)ObitUnref(in->myIO->tableList);
      in->myIO->tableList = (Obit*)ObitRef(in->tableList);
    }
    
    /* details depend on underlying file type,
       pass it down to ObitIO to call relevant routine */
    out = (ObitTable*)newObitIOTable(in->myIO, access, tabType, 
				     tabVer, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
    
    /* Find it? */
    if (out==NULL) {
      /*Obit_log_error(err, OBIT_InfoErr, "Could not find %s %s table  %d",in->name, 
	tabType, *tabVer); */
      return out;
    }
  } /* end specific ObitData */

  /* Is this an image? */
  isImage = ObitImageIsA(in);

  /* set Status unless readonly or not active */
  if (access!=OBIT_IO_ReadOnly) {
    /* Open if not */
    if (in->myStatus==OBIT_Inactive) {
      doClose = TRUE;
      /* Don't need to assign image buffer here */
      if (isImage) ((ObitImage*)in)->extBuffer = TRUE;  
      ObitDataOpen(in, OBIT_IO_ReadWrite, err);
    } else doClose = FALSE;
    if (in->myStatus != OBIT_Inactive) in->myStatus = OBIT_Modified;
    if (doClose) {
      ObitDataClose(in, err);
      if (isImage) ((ObitImage*)in)->extBuffer = FALSE;  /* May need buffer later */
    }
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

  /* Add info to table */
  if (out->tabType != tabType) { /* they could be the same */
    if (out->tabType) g_free(out->tabType);
    out->tabType = g_strdup(tabType);
  }
  out->tabVer  = *tabVer;
  /* Secret reference to host */ 
  out->myHost  = (Obit*)in;

  return out;
} /* end newObitDataTable */

/**
 * Destroy a specified table(s) associated with the input ObitData.  
 * The table is removed from the ObitTableList but it is not updated.
 * A call to ObitDataUpdateTables to update disk structures.
 * \param in       Pointer to object with associated tables.
 * \param tabType  The table type (e.g. "AIPS CC").
 * \param tabVer   Desired version number, may be zero in which case
 *                 the highest extant version is returned for read
 *                 and the highest+1 for write.
 *                 -1 => all versions of tabType
 * \param err      ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitDataZapTable (ObitData *in, gchar *tabType, olong tabVer, 
			     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTable *out=NULL;
  olong highVer, ver, iver;
  gchar *tabtype=NULL;
  gchar *routine = "ObitDataZapTable";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));
  g_assert(tabType!=NULL);

  /* Details depend on underlying file type,
     pass it down to ObitIO to call relevant routine */
  ver = tabVer;
  if (ver<0) ver = 0;
 
  /* Delete table or tables */
  retCode = OBIT_IO_DeleteErr;
  tabtype = g_strdup(tabType);  /* Make local copy of table type */
  if (tabVer>=0) {             /* fixed version number or highest */
    out = newObitDataTable(in, OBIT_IO_ReadOnly, tabtype, &ver, err);

    /* Does it exist? */
    if (out==NULL) {
      if (tabtype) g_free (tabtype); tabtype = NULL;
      return  OBIT_IO_OK;
    }

    /* Zap */
    out = ObitTableZap (out, err);
    if (err->error) {
      if (tabtype) g_free (tabtype); tabtype = NULL;
      Obit_traceback_val (err, routine, in->name, retCode);
    }
    /* remove from table list */
    ObitTableListRemove (in->tableList, tabtype, ver);  
    if (!in->isScratch)  /* No messages if a scratch object */
      /*Obit_log_error(err, OBIT_InfoErr, "Deleted %s %s table  %d",in->name, 
		     tabType, ver); */
    if (in->myStatus!=OBIT_Inactive) in->myStatus = OBIT_Modified;

  } else if (tabVer<0) {              /* all */
    highVer = ObitTableListGetHigh (in->tableList, tabtype);
    while (highVer>0) {
      iver = highVer;
      out = newObitDataTable(in, OBIT_IO_ReadWrite, tabtype, &iver, err);
      /* Zap */
      out = ObitTableZap (out, err);
      if (err->error) {
	if (tabtype) g_free (tabtype); tabtype = NULL;
	Obit_traceback_val (err, routine, in->name, retCode);
      }
      /* remove from table list */
      ObitTableListRemove (in->tableList, tabtype, iver);  
      if (!in->isScratch)  /* No messages if a scratch object */
	/*Obit_log_error(err, OBIT_InfoErr, "Deleted %s %s table  %d",in->name, 
		       tabType, iver); */
      if (in->myStatus!=OBIT_Inactive) in->myStatus = OBIT_Modified;

      /* more to do? */
      highVer = ObitTableListGetHigh (in->tableList, tabtype);
    }
  }
  
  if (tabtype) g_free (tabtype); tabtype = NULL;
  retCode = OBIT_IO_OK;
  return retCode;
} /* end ObitDataZapTable */

/**
 * Return a ObitHistory Object to the history associated with
 * the input ObitData.  
 * If such an object exists, a reference to it is returned,
 * else a new object is created and if access=OBIT_IO_WriteOnly
 * entered in the ObitTableList if appropriate (AIPS)
 * \param in       Pointer to object with associated tables.
 *                 This MUST have been opened before this call.
 * \param access   access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,
 *                 or OBIT_IO_WriteOnly).
 * \param err      ObitErr for reporting errors.
 * \return pointer to created ObitTable, NULL on failure.
 */
ObitHistory* 
newObitDataHistory (ObitData *in, ObitIOAccess access, ObitErr *err)
{
  ObitHistory *out;
  gboolean doClose;
  olong ver;
  gchar *name;
  ObitIOStatus myStatus;
  ObitIOCode retCode;
  gchar *routine = "newObitDataHistory";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA((Obit*)in, &myClassInfo));

  /* Create output */
  name = g_strconcat ("History of: ",in->name, NULL);
  out = newObitHistoryValue (name, in->info, err);
  g_free(name);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

  /* If WriteOnly, the IO object must be present - create if necessary */
  myStatus = in->myStatus;  /* May not change */
  if ((access==OBIT_IO_WriteOnly) && (out->myIO==NULL)) {

    /* Open and close to fully instantiate */
    retCode =  ObitHistoryOpen (out, access, err);
    retCode =  ObitHistoryClose (out, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  }

  /* If AIPS and WriteOnly pretend it's a table */
  if ((access==OBIT_IO_WriteOnly) && ObitIOHistoryAIPSIsA(out->myIO)) {
    ver = 1;
    ObitTableListPut(in->tableList, "AIPS HI", &ver, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    myStatus = OBIT_Modified;
  }

  /* set Status unless readonly or not active to force disk header update */
  if ((access!=OBIT_IO_ReadOnly) && ( myStatus==OBIT_Modified)) {
    /* Open if not */
    if (in->myStatus==OBIT_Inactive) {
      doClose = TRUE;
      ObitDataOpen(in, OBIT_IO_ReadWrite, err);
    } else doClose = FALSE;
    if (myStatus != OBIT_Inactive) in->myStatus = myStatus;
    if (doClose) ObitDataClose(in, err);
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
  
  return out;
} /* end newObitDataHistory */

/**
 * Copies the associated tables from one ObitData to another.
 * \param in      The ObitData with tables to copy.
 * \param out     An ObitData to copy the tables to, old ones replaced.
 * \param exclude a NULL termimated list of table types NOT to copy.
 *                If NULL, use include
 * \param include a NULL termimated list of table types to copy.
 *                ignored if exclude nonNULL.
 * \param err     ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitDataCopyTables (ObitData *in, ObitData *out, gchar **exclude,
			       gchar **include, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitIOType FileType;
  olong i, j, n, version;
  gchar *tabType = NULL;
  ObitTable *intable = NULL, *outtable=NULL;
  gboolean doit, isAIPS;
  gchar *deadly[] = {"AIPS HI", "AIPS SL", "AIPS PL", "History",
		     NULL}; /* NEVER try to copy these */
  gchar *routine = "ObitDataCopyTables";

  /* error check */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Is input or output AIPS? - special restrictions */
  isAIPS = FALSE;
  if (ObitGenericData(in)) {
    FileType = getIOType (in->info, err);
    isAIPS = FileType==OBIT_IO_AIPS;
  }
  if (ObitGenericData(out)) {
    FileType = getIOType (out->info, err);
    isAIPS = isAIPS || (FileType==OBIT_IO_AIPS);
  }

  retCode = OBIT_IO_ReadErr; /* in case of trouble */
  /* Loop through tables on in TableList */
  n = in->tableList->number;
  for (i=1; i<=n; i++) {
    /* Get info */
    ObitTableListGetNumber (in->tableList, i, &tabType, &version, 
			    &intable, err);

    /* are we looking for inclusions or exclusions? */
    if (exclude!=NULL) {
      /* Is this in the exclusion list? */
      doit = TRUE;
      j = 0;
      while (exclude[j]!=NULL) {
	doit = doit && strcmp (tabType, exclude[j]);
	j++;
      }
    } else {
      /* check the inclusion list */
      g_assert(include!=NULL);
      doit = FALSE;
     j = 0;
       while (include[j]!=NULL) {
	doit = doit || (!strcmp (tabType, include[j]));
	j++;
      }
   }

    /* Check for deadly AIPS extention files */
    j = 0;
    while (isAIPS && (deadly[j]!=NULL)) {
      doit = doit && (strcmp (tabType, deadly[j]));
      j++;
    }

    /* copy if wanted */
    if (doit) {
      /* setup input table if not instantiated */
      if (intable==NULL) {
 	intable = newObitDataTable (in, OBIT_IO_ReadOnly, tabType, 
				  &version, err);
	if (intable==NULL) {
	  Obit_log_error(err, OBIT_Error, 
			 "%s: No %s table found for %s", 
			 routine, tabType, in->name);
	  goto cleanup;
	}
	if (err ->error) goto cleanup;
      } /*  end init intable */
   
      /* setup output table */
      outtable = newObitDataTable (out, OBIT_IO_ReadWrite, tabType, 
				  &version, err);
      if (err ->error) goto cleanup;
      
      /* copy */
      outtable = ObitTableCopy (intable, outtable, err);
      if (err ->error) goto cleanup;
      
    } /* end of copy tables section */
      
    /* cleanup */
  cleanup: 
    if (tabType) g_free(tabType);
    intable  = ObitTableUnref (intable);
    outtable = ObitTableUnref (outtable);    
    if (err ->error) Obit_traceback_val (err, routine, in->name, retCode);
    if (out->myStatus != OBIT_Inactive) out->myStatus = OBIT_Modified;

  } /* end loop over tables on input */
  
  return OBIT_IO_OK;
} /* end ObitDataCopyTables */

/**
 * Update any disk resident structures about the current tables.
 * \param in   Pointer to object to be updated.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitDataUpdateTables (ObitData *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gboolean openClose;
  gchar *routine="ObitDataUpdateTables";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Trap Generic ObitData */
  if (ObitGenericData(in)) {
    /* set Status unless readonly or not active to force disk header update */
    if ((in->DataAccess!=OBIT_IO_ReadOnly) && (in->myStatus==OBIT_Modified)) {
      /* Open if not */
      if (in->myStatus==OBIT_Inactive) {
	openClose = TRUE;
	retCode = ObitDataOpen(in, OBIT_IO_ReadWrite, err);
      } else openClose = FALSE;
      if (openClose) retCode = ObitDataClose(in, err);
    }
    if (err->error) Obit_traceback_val (err, routine, in->name,retCode);
    return retCode;
  } /* end Generic ObitData */
  
  /* Specific ObitData */  
  /* Need to open and close? */
  openClose = !((in->myStatus==OBIT_Active) || 
		(in->myStatus==OBIT_Modified));

  if (openClose) {
    retCode = ObitIOOpen (in->myIO, OBIT_IO_ReadWrite, in->info, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  }
 
  /* do update */
  retCode = ObitIOUpdateTables(in->myIO, in->info, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_val (err, routine, in->name, retCode);

  if (openClose) {
    retCode = ObitIOClose (in->myIO, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
      Obit_traceback_val (err, routine, in->name, retCode);
  }
 
  return retCode;
} /* end ObitDataUpdateTables */

/**
 * Copies an associated table from one ObitData to another.
 * Any previous data in the output table will be lost.
 * \param in      The ObitData with tables to copy.
 * \param out     An ObitData to copy the tables to, old ones replaced.
 * \param tabType Table type, e.g. "AIPS CC"
 * \param inver   Input table version number, 0=>highest, actual returned
 * \param outver  Output table version number, 0=>new,  actual returned
 * \param err     ObitErr for reporting errors.
 */
void ObitDataCopyTable (ObitData *in, ObitData *out, gchar *tabType, 
			olong *inver, olong *outver, ObitErr *err)
{
  ObitTable *inTab=NULL, *outTab=NULL;
  gchar *routine="ObitDataCopyTable";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Get input table */
  inTab = newObitDataTable (in, OBIT_IO_ReadOnly, tabType, inver, err);
  if ((inTab==NULL) || (err->error)) goto cleanup;

  /* Instantiate */
  ObitTableFullInstantiate (inTab, TRUE, err);
  if (err->error) goto cleanup;

  /* Get output table */
  outTab = newObitDataTable (out, OBIT_IO_WriteOnly, tabType, outver, err);
  if ((outTab==NULL) || (err->error)) goto cleanup;

  /* Copy descriptor */
  outTab->myDesc = ObitTableDescCopy (inTab->myDesc, outTab->myDesc, err);
  if (err->error) goto cleanup;

  /* Instantiate */
  ObitTableFullInstantiate (outTab, FALSE, err);
  if (err->error) goto cleanup;

  /* Truncate output table */
  ObitTableClearRows (outTab, err);
  if (err->error) goto cleanup;

  /* Copy */
  outTab = ObitTableCopy (inTab, outTab, err);
  if (err->error) goto cleanup;

  /* Cleanup */
  cleanup: inTab = ObitTableUnref(inTab);
  outTab = ObitTableUnref(outTab);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  return;
} /* end ObitDataCopyTable */

/**
 * Reposition IO to beginning of file
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in   Pointer to object to be rewound.
 * \param err  ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitDataIOSet (ObitData *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataIOSet";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return retCode;
  } else { /* use class function */
    myClass = in->ClassInfo;
    retCode = myClass->ObitDataIOSet (in, err);
  }
  return retCode;
} /* end ObitDataIOSet */

/**
 * Create myIO object depending on value of FileType in in->info.
 * This is the principle place where the underlying file type is known.
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in   object to attach myIO
 * \param err  ObitErr for reporting errors.
 */
void ObitDataSetupIO (ObitData *in, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataSetupIO";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return;
  } else { /* use class function */
    myClass = in->ClassInfo;
    myClass->ObitDataSetupIO (in, err);
  }
} /* end ObitDataSetupIO */

/**
 * Write header keyword/value
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in   object to update, must be open during call with Write access
 * \param name The label (keyword) of the information. Max 8 char
 * \param type Data type of data element (enum defined in ObitInfoList class).
 * \param dim  Dimensionality of datum.  Only scalars and strings up to 8 char are allowed
 *             Note: for strings, the first element is the length in char.
 * \param data Pointer to the data. 
 * \param err  ObitErr for reporting errors.
 */
void ObitDataWriteKeyword (ObitData *in, 
			   gchar* name, ObitInfoType type, gint32 *dim, 
			   gconstpointer data, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataWriteKeyword";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return;
  } else { /* use class function */
    myClass = in->ClassInfo;
    /* Make sure defined */
    g_assert (myClass->ObitDataWriteKeyword!=
	      (ObitDataWriteKeywordFP)ObitDataWriteKeyword);
    myClass->ObitDataWriteKeyword (in, name, type, dim, data, err);
  }

} /* end ObitDataWriteKeyword */

/**
 * Read header keyword/value
 * Virtual - calls actual class member; not supported for Generic ObitData
 * \param in   object to update, must be fully instantiated
 * \param name [out] The label (keyword) of the information. Max 8 char
 * \param type [out] Data type of data element (enum defined in ObitInfoList class).
 * \param dim  [out] Dimensionality of datum.  Only scalars and strings up to 8 char 
 *                   are supported
 *                   Note: for strings, the first element is the length in char.
 * \param data [out] Pointer to the data. 
 * \param err  ObitErr for reporting errors.
 */
void ObitDataReadKeyword (ObitData *in, 
			  gchar* name, ObitInfoType *type, gint32 *dim, 
			  gpointer data, ObitErr *err)
{
  const ObitDataClassInfo *myClass;
  gchar *routine = "ObitDataReadKeyword";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Cannot do for generic ObitData */
  if (ObitGenericData(in)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Function not supported for generic ObitData %s", 
		   routine, in->name);
    return;
  } else { /* use class function */
    myClass = in->ClassInfo;
    /* Make sure defined */
    g_assert (myClass->ObitDataReadKeyword!=
	      (ObitDataReadKeywordFP)ObitDataReadKeyword);
    myClass->ObitDataReadKeyword (in, name, type, dim, data, err);
  }
  
} /* end ObitDataReadKeyword */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDataClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();
  myClassInfo.hasScratch  = TRUE; /* Scratch files allowed */

  /* Set function pointers */
  ObitDataClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDataClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDataClassInfoDefFn (gpointer inClass)
{
  ObitDataClassInfo *theClass = (ObitDataClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDataClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDataClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDataGetClass;
  theClass->newObit         = (newObitFP)newObitData;
  theClass->ObitGetClass    = (ObitGetClassFP)ObitDataGetClass;
  theClass->newObitDataScratch  = (newObitDataScratchFP)newObitDataScratch;
  theClass->ObitDataSame    = (ObitDataSameFP)ObitDataSame;
  theClass->ObitDataSetupIO = (ObitDataSetupIOFP)ObitDataSetupIO;
  theClass->ObitDataRename  = (ObitDataRenameFP)ObitDataRename;
  theClass->ObitDataZap     = (ObitDataZapFP)ObitDataZap;
  theClass->ObitClone       = NULL;  /* Different call */
  theClass->ObitClear       = (ObitClearFP)ObitDataClear;
  theClass->ObitInit        = (ObitInitFP)ObitDataInit;
  theClass->ObitDataCopy    = (ObitDataCopyFP)ObitDataCopy;
  theClass->ObitDataClone   = (ObitDataCloneFP)ObitDataClone;
  theClass->ObitDataOpen    = (ObitDataOpenFP)ObitDataOpen;
  theClass->ObitDataClose   = (ObitDataCloseFP)ObitDataClose;
  theClass->ObitDataIOSet   = (ObitDataIOSetFP)ObitDataIOSet;
  theClass->newObitDataTable= (newObitDataTableFP)newObitDataTable;
  theClass->ObitDataZapTable= (ObitDataZapTableFP)ObitDataZapTable;
  theClass->newObitDataHistory= 
    (newObitDataHistoryFP)newObitDataHistory;
  theClass->ObitDataFullInstantiate= 
    (ObitDataFullInstantiateFP)ObitDataFullInstantiate;
  theClass->ObitDataCopyTables= 
    (ObitDataCopyTablesFP)ObitDataCopyTables;
  theClass->ObitDataUpdateTables= 
    (ObitDataUpdateTablesFP)ObitDataUpdateTables;
  theClass->ObitDataCopyTable= 
    (ObitDataCopyTableFP)ObitDataCopyTable;
  theClass->ObitDataWriteKeyword= 
    (ObitDataWriteKeywordFP)ObitDataWriteKeyword;
  theClass->ObitDataReadKeyword= 
    (ObitDataReadKeywordFP)ObitDataReadKeyword;
} /* end ObitDataClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDataInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitData *in = inn;

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
  in->tableList = newObitTableList(in->name);
  in->isScratch = FALSE;

} /* end ObitDataInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *             Actually it should be an ObitData* cast to an Obit*.
 */
void ObitDataClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitData *in = inn;
  ObitErr *err;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Delete underlying files if isScratch */
  if (in->isScratch) {
    err = newObitErr();     /* for possible messages */
    /* Remove from ObitSystem list */
    ObitSystemFreeScratch ((Obit*)in, err);
    in->isScratch = FALSE;  /* avoid infinite recursion */
    ObitDataZap (in, err);    /* delete files */
    ObitErrLog(err);
    err = ObitErrUnref(err);
  }

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->myIO      = ObitUnref(in->myIO);
  in->tableList = ObitUnref(in->tableList);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDataClear */

/**
 * Determine File Type
 * \param  Infolist to read
 * \param  err ObitErr for reporting errors.
 * \return Io type, OBIT_IO_FITS, OBIT_IO_AIPS
 */
static ObitIOType getIOType (ObitInfoList *list, ObitErr *err)
{
  ObitIOType FileType=OBIT_IO_Binary;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitData:getIOTyp";

  if (err->error) return FileType;
  
  /* Get FileType */
  if (!ObitInfoListGet(list, "FileType", &type, dim, 
		       &FileType, err)) {
    /* couldn't find it - add message to err and return */
    Obit_log_error(err, OBIT_Error, 
		   "%s: entry FileType not in InfoList", routine);
  }

    return FileType;
} /* end  getIOType*/
