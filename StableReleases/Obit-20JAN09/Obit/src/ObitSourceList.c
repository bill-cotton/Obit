/* $Id$  */
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

#include "ObitSourceList.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitSourceList.c
 * ObitSourceList class function definitions.
 *
 * This is a list of (astronomical) sources.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitSourceList";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitSourceListClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSourceListClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSourceListInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSourceListClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSourceListClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitSourceList* newObitSourceList (gchar* name)
{
  ObitSourceList* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSourceListClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSourceList));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSourceListInit((gpointer)out);

  return out;
} /* end newObitSourceList */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitSourceListGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSourceListClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Creates an ObitSourceList of a given size.
 * \param name  An optional name for the object.
 * \param nsou  Number of sources (actually, the highest Source ID).
 * \return the new object.
 */
ObitSourceList* ObitSourceListCreate (gchar* name, olong nsou)
{
  ObitSourceList* out;
  olong i;
  gchar sname[31];

  /* Create basic structure */
  out = newObitSourceList (name);

  /* create data array if wanted */
  if (nsou<0) return out;

  /* Save information */
  out->number = nsou;

  /* create array */
  out->SUlist = g_malloc0 (nsou*sizeof(ObitSource*));

  /* Create elements of ObitSource list */
  for (i=0; i<out->number; i++) {
    g_snprintf (sname, 30, "Source  %d", i+1);
    out->SUlist[i] = newObitSource (sname);
  }

  return out;
} /* end ObitSourceListCreate */

/**
 * Make a copy of a object.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSourceList* ObitSourceListCopy  (ObitSourceList *in, ObitSourceList *out, 
				   ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;
  gchar *routine = "ObitSourceListCopy";

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
    out = newObitSourceList(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* number of elements */
  out->number = in->number;

  /* Reallocate list if needed */
  out->SUlist = g_realloc (out->SUlist, out->number*sizeof(ObitSource*));

  /* loop through list copying elements */
  for (i=0; i<out->number; i++) 
    out->SUlist[i] = ObitSourceCopy (in->SUlist[i], out->SUlist[i], err);
  
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  
  return out;
} /* end ObitSourceListCopy */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSourceListClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSourceListClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSourceListClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSourceListClassInfoDefFn (gpointer inClass)
{
  ObitSourceListClassInfo *theClass = (ObitSourceListClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSourceListClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSourceListClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSourceListGetClass;
  theClass->newObit       = (newObitFP)newObitSourceList;
  theClass->ObitCopy      = (ObitCopyFP)ObitSourceListCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSourceListClear;
  theClass->ObitInit      = (ObitInitFP)ObitSourceListInit;

} /* end ObitSourceListClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitSourceListInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSourceList *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->SUlist = NULL;
  in->number = 0;
} /* end ObitSourceListInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitSourceListClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSourceList *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  /* loop through list copying elements */
  for (i=0; i<in->number; i++) 
    in->SUlist[i] = ObitSourceUnref (in->SUlist[i]);

  /* delete members */
  if (in->SUlist) g_free(in->SUlist); in->SUlist = NULL;
  
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitSourceListClear */

