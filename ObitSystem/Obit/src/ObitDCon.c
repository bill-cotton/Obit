/* $Id$        */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitDCon.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDCon.c
 * ObitDCon class function definitions.
 * Virtual deconvolution base class.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDCon";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitDConClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * VIRTUAL routine - should never be called - 
 * defined for convenience of derived classes 
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDCon* newObitDCon (gchar* name)
{
  ObitDCon* out;
  gchar *routine = "newObitDCon";

  /* VIRTUAL */
  g_error("%s: Virtual routine - should not be called",routine);

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDCon));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConInit((gpointer)out);

 return out;
} /* end newObitDCon */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConGetClass */

/**
 * Make a deep copy of an ObitDCon.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDCon* ObitDConCopy  (ObitDCon *in, ObitDCon *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitImageMosaicClassInfo *mosaicClass;
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
    out = newObitDCon(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* free any old mosaic */
  out->mosaic = ObitImageMosaicUnref(out->mosaic);
  /* Copy new */
  mosaicClass = (ObitImageMosaicClassInfo*)in->mosaic->ClassInfo;
  out->mosaic = (ObitImageMosaic*)mosaicClass->ObitCopy((ObitImageMosaic*)in->mosaic, NULL, err);

  return out;
} /* end ObitDConCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DCon similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConClone  (ObitDCon *in, ObitDCon *out, ObitErr *err)
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
  /* free any old mosaic */
  out->mosaic = ObitImageMosaicUnref(out->mosaic);
  /* Copy new */
  out->mosaic = ObitImageMosaicRef(in->mosaic);

} /* end ObitDConClone */

/**
 * Creates an ObitDCon 
 * VIRTUAL routine - should never be called - 
 * defined for convenience of derived classes 
 * \param name   An optional name for the object.
 * \param mosaic Image mosaic to be deconvolved.
 * \return the new object.
 */
ObitDCon* ObitDConCreate (gchar* name, ObitImageMosaic *mosaic,  
			  ObitErr *err)
{
  ObitDCon* out=NULL;
  gchar *routine = "ObitDConCreate";

  /* VIRTUAL */
  Obit_log_error(err, OBIT_Error,"%s: Virtual routine - should not be called",routine);
  return out;


  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitImageMosaicIsA(mosaic));

  /* Create basic structure */
  out = newObitDCon (name);

  /* Save Image Mosaic reference */
  out->mosaic = ObitImageMosaicRef(mosaic);
  return out;
} /* end ObitDConCreate */

/**
 * Read any base class parameters and then
 * read control parameters from the ObitInfoList member:
 * \li "Plane" OBIT_long array = Plane being processed, 1-rel indices of axes 3-?
 * \li "prtLv" OBIT_int message level  [def 2]
 *             0=none, 1=summary, 2=normal, higher numbers for diagnostics
 * \param in  The CLEAN object as base class
 * \param err Obit error stack object.
 */
void  ObitDConGetParms (ObitDCon *in, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  /* olong itemp;
     union ObitInfoListEquiv InfoReal; */
  /*gchar *routine = "ObitDConGetParms";*/

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* This is the base class, no parents */

  /* Plane to process */
  ObitInfoListGetTest(in->info, "Plane", &type, dim, &in->plane);

  /* Print level */
  in->prtLv = 2;  /* default = normal */
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &in->prtLv);
} /* end ObitDConGetParms */

/**
 * Do deconvolution, uses function on class pointer
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConDeconvolve (ObitDCon *in, ObitErr *err)
{
  const ObitDConClassInfo *inClass;
  gchar *routine = "ObitDConDeconvolve";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConIsA(in));

  inClass = (ObitDConClassInfo*)in->ClassInfo; /* class structure */

  /* Better be a derived class */
  if ((gpointer)inClass==(gpointer)&myClassInfo) {
     Obit_log_error(err, OBIT_Error,
		    "%s: Deconvolution with virtual class not allowed for %s",
		    routine, in->name);
      return;
 }
  
  /* Call derived function */
  inClass->ObitDConDeconvolve(in, err);

} /* end ObitDConDeconvolve */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConClassInfoDefFn (gpointer inClass)
{
  ObitDConClassInfo *theClass = (ObitDConClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConGetClass;
  theClass->newObit       = (newObitFP)newObitDCon;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConInit;
  theClass->ObitDConGetParms   = (ObitDConGetParmsFP)ObitDConGetParms;
  theClass->ObitDConDeconvolve = (ObitDConDeconvolveFP)ObitDConDeconvolve;

} /* end ObitDConClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitDCon *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread = newObitThread();
  in->info   = newObitInfoList(); 
  in->mosaic = NULL;
  in->prtLv = 2;  /* Normal messages */
  for (i=0; i<IM_MAXDIM-2; i++) in->plane[i] = 1;

} /* end ObitDConInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDCon* cast to an Obit*.
 */
void ObitDConClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDCon *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread = ObitThreadUnref(in->thread);
  in->info   = ObitInfoListUnref(in->info);
  in->mosaic = ObitImageMosaicUnref(in->mosaic);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConClear */


