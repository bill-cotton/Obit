/* $Id: ObitSource.c,v 1.9 2006/07/06 11:25:55 bcotton Exp $      */
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

#include "ObitSource.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitSource.c
 * ObitSource class function definitions.
 *
 * This is (astronomical) source information.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitSource";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitSourceClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSourceClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSourceInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSourceClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSourceClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitSource* newObitSource (gchar* name)
{
  ObitSource* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSourceClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSource));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSourceInit((gpointer)out);

  return out;
} /* end newObitSource */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitSourceGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSourceClassInit();

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
ObitSource* ObitSourceCopy  (ObitSource *in, ObitSource *out, 
				   ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;

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
    out = newObitSource(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->SourID = in->SourID;
  out->Qual   = in->Qual;
  out->numIF  = in->numIF;
  out->equinox= in->equinox;
  out->SourID = in->SourID;
  out->RAApp  = in->RAApp;
  out->DecApp = in->DecApp;
  out->RAMean = in->RAMean;
  out->DecMean= in->DecMean;
  out->Bandwidth = in->Bandwidth;
  for (i=0; i<16; i++) out->SourceName[i] = in->SourceName[i];
  out->SourceName[i] = 0;
  for (i=0; i<5; i++)  out->CalCode[i]    = in->CalCode[i];
  out->CalCode[i] = 0;
  if (out->IFlux)     g_free (out->IFlux);
  if (out->QFlux)     g_free (out->QFlux);
  if (out->UFlux)     g_free (out->UFlux);
  if (out->VFlux)     g_free (out->VFlux);
  if (out->FreqOff)   g_free (out->FreqOff);
  if (out->LSRVel)    g_free (out->LSRVel);
  if (out->RestFreq)  g_free (out->RestFreq);
  if (out->numIF>0) {
    out->IFlux     = g_malloc(out->numIF*sizeof(ofloat));
    out->QFlux     = g_malloc(out->numIF*sizeof(ofloat));
    out->UFlux     = g_malloc(out->numIF*sizeof(ofloat));
    out->VFlux     = g_malloc(out->numIF*sizeof(ofloat));
    out->FreqOff   = g_malloc(out->numIF*sizeof(odouble));
    out->LSRVel    = g_malloc(out->numIF*sizeof(odouble));
    out->RestFreq  = g_malloc(out->numIF*sizeof(odouble));
    for (i=0; i<out->numIF; i++)  out->IFlux[i]    = in->IFlux[i];
    for (i=0; i<out->numIF; i++)  out->QFlux[i]    = in->QFlux[i];
    for (i=0; i<out->numIF; i++)  out->UFlux[i]    = in->UFlux[i];
    for (i=0; i<out->numIF; i++)  out->VFlux[i]    = in->VFlux[i];
    for (i=0; i<out->numIF; i++)  out->FreqOff[i]  = in->FreqOff[i];
    for (i=0; i<out->numIF; i++)  out->LSRVel[i]   = in->LSRVel[i];
    for (i=0; i<out->numIF; i++)  out->RestFreq[i] = in->RestFreq[i];
  }
  return out;
} /* end ObitSourceCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitSourceClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSourceClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSourceClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSourceClassInfoDefFn (gpointer inClass)
{
  ObitSourceClassInfo *theClass = (ObitSourceClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSourceClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSourceClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSourceGetClass;
  theClass->newObit       = (newObitFP)newObitSource;
  theClass->ObitCopy      = (ObitCopyFP)ObitSourceCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSourceClear;
  theClass->ObitInit      = (ObitInitFP)ObitSourceInit;

} /* end ObitSourceClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitSourceInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSource *in = inn;
  olong i;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->SourID  = -1;
  in->Qual    = 0;
  in->numIF   = 0;
  in->equinox = 0.0;
  in->RAApp   = 0.0;
  in->DecApp  = 0.0;
  in->RAMean  = 0.0;
  in->DecMean = 0.0;
  in->Bandwidth = 0.0;
  for (i=0; i<16; i++) in->SourceName[i] = ' '; in->SourceName[i] = 0;
  for (i=0; i<4; i++)  in->CalCode[i] = ' ';    in->CalCode[i]    = 0;
  in->IFlux     = NULL;
  in->QFlux     = NULL;
  in->UFlux     = NULL;
  in->VFlux     = NULL;
  in->FreqOff   = NULL;
  in->LSRVel    = NULL;
  in->RestFreq  = NULL;

} /* end ObitSourceInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitSourceClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSource *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->IFlux)     g_free (in->IFlux);
  if (in->QFlux)     g_free (in->QFlux);
  if (in->UFlux)     g_free (in->UFlux);
  if (in->VFlux)     g_free (in->VFlux);
  if (in->FreqOff)   g_free (in->FreqOff);
  if (in->LSRVel)    g_free (in->LSRVel);
  if (in->RestFreq)  g_free (in->RestFreq);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSourceClear */

