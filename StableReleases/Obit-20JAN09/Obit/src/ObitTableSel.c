/* $Id$    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
#include "Obit.h"
#include "ObitTableSel.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableSel.c
 * ObitTableSel Obit uv data selector class definition.
 * This contains information about data selection and calibration.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitTableSel";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableSelClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitTableSelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTableSelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitTableSelClassInfoDefFn (gpointer inClass);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitTableSel* newObitTableSel (gchar *name)
{
  ObitTableSel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableSelClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitTableSel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTableSelInit((gpointer)out);

  return out;
} /* end newObitTableSel */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitTableSelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitTableSelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitTableSelGetClass */

/**
 * Copy constructor.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitTableSel* ObitTableSelCopy (ObitTableSel* in, ObitTableSel* out, 
			  ObitErr *err)
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
    out = newObitTableSel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* This class members */
  out->FileType = in->FileType;

  return out;
} /* end ObitTableSelCopy */

/**
 * Determines how large a buffer (in bytes) is needed
 * for data transfers as described by data members.
 * \param desc Pointer input descriptor.
 * \param sel Table selector.
 * \return size in bytes needed for I/O.
 */
olong ObitTableSelBufferSize (ObitTableDesc* desc, 
			       ObitTableSel* sel)
{
  olong size = 0;

  /* error checks */
  if (desc==NULL) return size; 
  g_assert (ObitIsA(desc, ObitTableDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));

  /* make sure defaults filled in */
  ObitTableSelDefault (desc, sel);

  /* size in bytes needed */
  size = sel->nRowPIO * desc->lrow;

  return size;
} /* end ObitTableSelBufferSize */

/**
 * Enforces any defaults in the descriptor.
 * Also indexes structure.
 * \param in Pointer to descriptor.
 * \param sel Table selector.
 */
void ObitTableSelDefault (ObitTableDesc* in, ObitTableSel* sel)
{

  /* error checks */
  g_assert (ObitIsA(in, ObitTableDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));

  /* Index as well */
  ObitTableDescIndex(in);
} /* end ObitTableSelDefault */

/**
 * Apply selection criteria to input descriptor to derive output.
 * \param in Pointer to input descriptor, this described the data
 *           as they appear on disk (possible compressed).
 * \param sel Table selector, blc, trc members changed if needed.
 * \param out Pointer to output descriptor, this describes the data 
 *            after any processing when read, or before any compression
 *            on output.
 * \param err Obit error stack
 */
void ObitTableSelSetDesc (ObitTableDesc* in, ObitTableSel* sel,
			  ObitTableDesc* out, ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, ObitTableDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));
  g_assert (ObitIsA(out, ObitTableDescGetClass()));

  /* make sure defaults filled in */
  ObitTableSelDefault (in, sel);

  /* copy most values */
  ObitTableDescCopy (in, out, err);
  if (err->error) /* add traceback, return on error */
      Obit_traceback_msg (err, "ObitTableSelSetDesc", 
			  in->name);

  /* make sure defaults, indices filled in */
  ObitTableSelDefault (in, sel);
  ObitTableSelDefault (out, sel);

} /* end ObitTableSelSetDesc */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableSelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTableSelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTableSelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTableSelClassInfoDefFn (gpointer inClass)
{
  ObitTableSelClassInfo *theClass = (ObitTableSelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTableSelClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTableSelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTableSelGetClass;
  theClass->newObit       = (newObitFP)newObitTableSel;
  theClass->ObitCopy      = (ObitCopyFP)ObitTableSelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitTableSelClear;
  theClass->ObitInit      = (ObitInitFP)ObitTableSelInit;

} /* end ObitTableSelClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitTableSelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableSel *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
} /* end ObitTableSelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitTableSelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableSel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitTableSelClear */

