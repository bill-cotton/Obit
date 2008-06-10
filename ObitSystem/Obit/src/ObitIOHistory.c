/* $Id: ObitIOHistory.c,v 1.6 2006/06/19 14:51:27 bcotton Exp $   */
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

#include "ObitIOHistory.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIOHistory.c
 * ObitIOHistory class function definitions.
 *
 * This is a virtual base class and should never be directly instantiated.
 * Derived classes provide an I/O interface to various underlying disk
 * structures.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitIOHistory";

/** Function to obtain parent ClassInfo - ObitIO */
static ObitGetClassFP ObitParentGetClass = ObitIOGetClass;

/**
 * ClassInfo global structure ObitIOHistoryClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIOHistoryClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIOHistoryInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIOHistoryClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitIOHistoryClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name Name [optional] for object
 * \param info InfoList defining file
 * \param err ObitErr for reporting errors. 
 * \return the new object.
 */
ObitIOHistory* newObitIOHistory (gchar* name, ObitInfoList *info,
				 ObitErr *err)
{
  ObitIOHistory* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOHistoryClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitIOHistory));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitIOHistoryInit((gpointer)out);

  return out;
} /* end newObitIOHistory */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitIOHistoryGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIOHistoryClassInit();
  
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
gboolean ObitIOHistorySame (ObitIOHistory *in, ObitInfoList *in1, 
			    ObitInfoList *in2, ObitErr *err)
{
  const ObitIOHistoryClassInfo *myClass;

  /* Don't bother if NULL */
  if (!in) return FALSE;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return FALSE;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOSame != NULL);

  /* call actual function */
  return myClass->ObitIOSame ((ObitIO*)in, in1, in2, err);

} /* end ObitIOHistorySame */

/**
 * Delete underlying files.
 * \param in Pointer to object to be zapped.
 * \param err ObitErr for reporting errors. 
 */
void ObitIOHistoryZap (ObitIOHistory *in, ObitErr *err)
{
  const ObitIOHistoryClassInfo *myClass;

  /* Don't bother if NULL */
  if (!in) return;
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOZap != NULL);

  /* call actual function */
  myClass->ObitIOZap ((ObitIO*)in, err);

  return;
} /* end ObitIOHistoryZap */

/**
 * Make a copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object. 
 * \return pointer to the new object.
 */
ObitIOHistory* ObitIOHistoryCopy  (ObitIOHistory *in, ObitIOHistory *out, ObitErr *err)
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
    out = newObitIOHistory(outName, NULL, err);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->info = ObitInfoListUnref(out->info);
  out->info = ObitInfoListCopy(in->info);
  out->MaxNumber     = in->MaxNumber;
  out->CurrentNumber = in->CurrentNumber;
  out->Last          = in->Last;
  return out;
} /* end ObitIOHistoryCopy */

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
ObitIOCode ObitIOHistoryOpen (ObitIOHistory *in, ObitIOAccess access, ObitInfoList *info, 
			      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIOHistoryIsA(in));
  g_assert (ObitInfoListIsA (info));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOOpen != NULL);

  /* call actual function */
  retCode = myClass->ObitIOOpen ((ObitIO*)in, access, info, err);
  in->access = access; /* just in case */

  return retCode;
} /* end ObitIOHistoryOpen */

/**
 * Shutdown I/O.
 * \param in Pointer to object to be closed.
 * \param err ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryClose (ObitIOHistory *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOClose != NULL);

  /* call actual function */
  retCode = myClass->ObitIOClose ((ObitIO*)in, err);

  return retCode;
} /* end ObitIOHistoryClose */

/**
 * initialize I/O
 * \param in Pointer to object to be accessed.
 * \param info ObitInfoList with instructions
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistorySet (ObitIOHistory *in, ObitInfoList *info, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitInfoListIsA (info));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOSet != NULL);

  /* call actual function */
  retCode = myClass->ObitIOSet ((ObitIO*)in, info, err);

  return retCode;
} /* end ObitIOHistoryInit */

/**
 * Read specified History record
 * \param in     Pointer to object to be read.
 * \param recno  record number (1-rel) -1=> next.
 * \param hiCard output history record (70 char)
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryReadRec (ObitIOHistory *in, olong recno, gchar *hiCard, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (hiCard != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOHistoryReadRec != NULL);

  /* call actual function */
  retCode = myClass->ObitIOHistoryReadRec (in, recno, hiCard, err);

  return retCode;
} /* end ObitIOHistoryReadRec */

/**
 * Write specified History record
 * \param in     Pointer to object to be written.
 * \param recno  Record number (1-rel) -1=> next, overwrites any existing
 * \param hiCard input history record (70 char)
 * \param err    ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryWriteRec (ObitIOHistory *in, olong recno, gchar *hiCard, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (hiCard != NULL);

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOHistoryWriteRec != NULL);

  /* call actual function */
  retCode = myClass->ObitIOHistoryWriteRec (in, recno, hiCard, err);

  return retCode;
} /* end ObitIOHistoryWriteRow */

/**
 * Tell number of History records
 * \param in     Pointer to open object to be tested
 * \return  number of records, <0 => problem
 */
olong ObitIOHistoryNumRec (ObitIOHistory *in)
{
  olong out;
  const ObitIOHistoryClassInfo *myClass;

  myClass = in->ClassInfo;
  out =  myClass->ObitIOHistoryNumRec(in);

  return out;
} /* end ObitIOHistoryNumRec */

/**
 * Read image Descriptor data from disk.
 * \param in Pointer to object  with ObitImageDescto be read.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryReadDescriptor (ObitIOHistory *in,  ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOReadDescriptor !=NULL);

  /* call actual function */
  retCode = myClass->ObitIOReadDescriptor ((ObitIO*)in, err);

  return retCode;
} /* end ObitIOHistoryReadDescriptor */

/**
 * Write Descriptor information to disk.
 * \param in Pointer to object with ObitImageDesc to be written.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryWriteDescriptor (ObitIOHistory *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOWriteDescriptor != NULL);

  /* call actual function */
  retCode = myClass->ObitIOWriteDescriptor ((ObitIO*)in, err);

  return retCode;
} /* end ObitIOHistoryWriteDescriptor */

/**
 * Flush I/O buffer if necessary 
 * \param in Pointer to object to be accessed.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitIOHistoryFlush (ObitIOHistory *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitIOHistoryClassInfo *myClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));

  /* this is a virtual function, see if actual one defined */
  myClass = in->ClassInfo;
  g_assert (myClass->ObitIOFlush != NULL);

  /* call actual function */
  retCode = myClass->ObitIOFlush ((ObitIO*)in, err);

  return retCode;
} /* end ObitIOHistoryFlush */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIOHistoryClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIOHistoryClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIOHistoryClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIOHistoryClassInfoDefFn (gpointer inClass)
{
  ObitIOHistoryClassInfo *theClass = (ObitIOHistoryClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIOHistoryClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIOHistoryClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIOHistoryGetClass;
  theClass->newObit       = NULL;
  theClass->newObitIO     = (newObitIOFP)newObitIOHistory;
  theClass->ObitCopy      = (ObitCopyFP)ObitIOHistoryCopy;
  theClass->ObitIOSame    = (ObitIOSameFP)ObitIOHistorySame;
  theClass->ObitIOZap     = (ObitIOZapFP)ObitIOHistoryZap;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitIOHistoryClear;
  theClass->ObitInit      = (ObitInitFP)ObitIOHistoryInit;
  theClass->ObitIOOpen    = (ObitIOOpenFP)ObitIOHistoryOpen;
  theClass->ObitIOClose   = (ObitIOCloseFP)ObitIOHistoryClose;
  theClass->ObitIOSet     = (ObitIOSetFP)ObitIOHistorySet;
  theClass->ObitIORead          = NULL;
  theClass->ObitIOReadSelect    = NULL;
  theClass->ObitIOReadRowSelect = NULL;
  theClass->ObitIOWriteRow      = NULL;
  theClass->ObitIOWrite         = NULL;
  theClass->ObitIOFlush   = (ObitIOFlushFP)ObitIOHistoryFlush;
  theClass->ObitIOReadDescriptor  = 
    (ObitIOReadDescriptorFP)ObitIOHistoryReadDescriptor;
  theClass->ObitIOWriteDescriptor = 
    (ObitIOWriteDescriptorFP)ObitIOHistoryWriteDescriptor;
  theClass->ObitIOCreateBuffer = NULL;
  theClass->ObitIOFreeBuffer   = NULL;
  theClass->newObitIOTable     = NULL;
  theClass->ObitIOUpdateTables = NULL;
  /* added this class */
  theClass->ObitIOHistoryReadRec  = 
    (ObitIOHistoryReadRecFP)ObitIOHistoryReadRec;
  theClass->ObitIOHistoryWriteRec = 
    (ObitIOHistoryWriteRecFP)ObitIOHistoryWriteRec;
  theClass->ObitIOHistoryNumRec  = 
    (ObitIOHistoryNumRecFP)ObitIOHistoryNumRec;

} /* end ObitIOHistoryClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitIOHistoryInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIOHistory *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->info = newObitInfoList();
} /* end ObitIOHistoryInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitIOHistoryClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIOHistory *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  in->info = ObitInfoListUnref(in->info);
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitIOHistoryClear */

