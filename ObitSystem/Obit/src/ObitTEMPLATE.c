/* $Id: ObitTEMPLATE.c,v 1.7 2007/01/29 15:15:37 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007                                               */
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

#include "ObitTEMPLATE.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTEMPLATE.c
 * ObitTEMPLATE class function definitions.
 * This class is derived from the Obit base class.
 */

/* *************** CHANGE HERE *********************************  */
/** name of the class defined in this file */
static gchar *myClassName = "ObitTEMPLATE";

/* *************** CHANGE HERE *********************************  */
/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitTEMPLATEClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTEMPLATEClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitTEMPLATEInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTEMPLATEClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitTEMPLATEClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTEMPLATE* newObitTEMPLATE (gchar* name)
{
  ObitTEMPLATE* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTEMPLATEClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTEMPLATE));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTEMPLATEInit((gpointer)out);

 return out;
} /* end newObitTEMPLATE */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitTEMPLATEGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTEMPLATEClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitTEMPLATEGetClass */

/**
 * Make a deep copy of an ObitTEMPLATE.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitTEMPLATE* ObitTEMPLATECopy  (ObitTEMPLATE *in, ObitTEMPLATE *out, ObitErr *err)
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
    out = newObitTEMPLATE(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* *************** CHANGE HERE *********************************  */

  return out;
} /* end ObitTEMPLATECopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an TEMPLATE similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitTEMPLATEClone  (ObitTEMPLATE *in, ObitTEMPLATE *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* *************** CHANGE HERE *********************************  */

} /* end ObitTEMPLATEClone */

  /* *************** CHANGE HERE *********************************  */
/**
 * Creates an ObitTEMPLATE 
 * \param name  An optional name for the object.
 * \return the new object.
 */
  /* *************** CHANGE HERE *********************************  */
ObitTEMPLATE* ObitTEMPLATECreate (gchar* name)
{
  ObitTEMPLATE* out;

  /* Create basic structure */
  out = newObitTEMPLATE (name);

  /* *************** CHANGE HERE *********************************  */
  return out;
} /* end ObitTEMPLATECreate */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTEMPLATEClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTEMPLATEClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTEMPLATEClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTEMPLATEClassInfoDefFn (gpointer inClass)
{
  ObitTEMPLATEClassInfo *theClass = (ObitTEMPLATEClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTEMPLATEClassInit;
  theClass->newObit       = (newObitFP)newObitTEMPLATE;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTEMPLATEClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTEMPLATEGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitTEMPLATECopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitTEMPLATEClear;
  theClass->ObitInit      = (ObitInitFP)ObitTEMPLATEInit;
  theClass->ObitTEMPLATECreate = (ObitTEMPLATECreateFP)ObitTEMPLATECreate;

  /* *************** CHANGE HERE *********************************  */

} /* end ObitTEMPLATEClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTEMPLATEInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTEMPLATE *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  /* *************** CHANGE HERE *********************************  */

} /* end ObitTEMPLATEInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitTEMPLATE* cast to an Obit*.
 */
void ObitTEMPLATEClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTEMPLATE *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  /* *************** CHANGE HERE *********************************  */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTEMPLATEClear */

