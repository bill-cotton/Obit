/* $Id: ObitAntenna.c,v 1.5 2007/06/13 18:58:29 bcotton Exp $     */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitAntenna.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitAntenna.c
 * ObitAntenna class function definitions.
 *
 * This is antenna information.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitAntenna";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitAntennaClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitAntennaClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitAntennaInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitAntennaClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitAntennaClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitAntenna* newObitAntenna (gchar* name)
{
  ObitAntenna* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitAntennaClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitAntenna));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitAntennaInit((gpointer)out);

  return out;
} /* end newObitAntenna */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitAntennaGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitAntennaClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Make a copy of a object.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitAntenna* ObitAntennaCopy (ObitAntenna *in, ObitAntenna *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
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
    out = newObitAntenna(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  for (i=0; i<8; i++) out->AntName[i] = in->AntName[i];
  out->AntID = in->AntID;
  for (i=0; i<3; i++) out->AntXYZ[i] = in->AntXYZ[i];
  out->AntMount = in->AntMount;
  out->AntLong  = in->AntLong;
  out->AntLat   = in->AntLat;
  out->AntRad   = in->AntRad;
  out->FeedAPA  = in->FeedAPA;
  out->FeedBPA  = in->FeedBPA;
  out->FeedBType= in->FeedBType;
  out->FeedAType= in->FeedAType;
  out->numPCal  = in->numPCal;
  out->FeedAPCal =  g_realloc(out->FeedAPCal, in->numPCal*sizeof(ofloat));
  for (i=0; i<out->numPCal; i++) out->FeedAPCal[i] = in->FeedAPCal[i];
  out->FeedBPCal =  g_realloc(out->FeedBPCal, in->numPCal*sizeof(ofloat));
  for (i=0; i<out->numPCal; i++) out->FeedBPCal[i] = in->FeedBPCal[i];

  return out;
} /* end ObitAntennaCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitAntennaClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitAntennaClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitAntennaClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitAntennaClassInfoDefFn (gpointer inClass)
{
  ObitAntennaClassInfo *theClass = (ObitAntennaClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitAntennaClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitAntennaClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitAntennaGetClass;
  theClass->newObit       = (newObitFP)newObitAntenna;
  theClass->ObitCopy      = (ObitCopyFP)ObitAntennaCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitAntennaClear;

} /* end ObitAntennaClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitAntennaInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitAntenna *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->AntID     = -1;
  in->numPCal   = 0;
  in->FeedAPCal = NULL;
  in->FeedBPCal = NULL;
} /* end ObitAntennaInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitAntennaClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitAntenna *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->FeedAPCal) g_free(in->FeedAPCal); in->FeedAPCal = NULL;
  if (in->FeedBPCal) g_free(in->FeedBPCal); in->FeedBPCal = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitAntennaClear */

  
