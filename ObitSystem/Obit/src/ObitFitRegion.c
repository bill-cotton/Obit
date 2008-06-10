/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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

#include "ObitFitRegion.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFitRegion.c
 * ObitFitRegion class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitFitRegion";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitFitRegionClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFitRegionClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFitRegionInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFitRegionClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFitRegionClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitFitRegion* newObitFitRegion (gchar* name)
{
  ObitFitRegion* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFitRegionClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFitRegion));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFitRegionInit((gpointer)out);

 return out;
} /* end newObitFitRegion */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFitRegionGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFitRegionClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitFitRegionGetClass */

/**
 * Make a deep copy of an ObitFitRegion.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitFitRegion* ObitFitRegionCopy  (ObitFitRegion *in, ObitFitRegion *out, 
				   ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitFitRegion(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->corner[0] = in->corner[0];
  out->corner[1] = in->corner[1];
  out->dim[0]    = in->dim[0];
  out->dim[1]    = in->dim[1];
  out->peak      = in->peak;
  out->peakResid = in->peakResid;
  out->RMSResid  = in->RMSResid;
  out->fluxResid = in->fluxResid;
  out->nmodel    = in->nmodel;
  if (out->models) {
    for (i=0; i<out->nmodel; i++) 
      out->models[i] = ObitFitModelUnref(out->models[i]);
  }
  out->models    = g_realloc(out->models, out->nmodel*sizeof(ObitFitModel*));
  for (i=0; i<out->nmodel; i++) {
    out->models[i] = NULL;
    out->models[i] = ObitFitModelCopy(in->models[i], out->models[i], err);
  }

  return out;
} /* end ObitFitRegionCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an FitRegion similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitFitRegionClone  (ObitFitRegion *in, ObitFitRegion *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nmodel    = in->nmodel;
  out->models    = g_realloc(out->models, out->nmodel*sizeof(ObitFitModel*));

} /* end ObitFitRegionClone */

/**
 * Change number of FitReqionModels allowed
 * Existing models will be lost.
 * \param in     The object to resize
 * \param nmodel New number of models
 */
void ObitFitRegionResize  (ObitFitRegion *in, olong nmodel)
{
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Unref existing */
  if (in->models) {
    for (i=0; i<in->nmodel; i++) {
      in->models[i] = ObitFitModelUnref(in->models[i]);
    } 
  }

  /* Resize */
  in->nmodel = nmodel;
  in->models = g_realloc(in->models, nmodel*sizeof(ObitFitModel*));

} /* end ObitFitRegionResize */

/**
 * Creates an ObitFitRegion 
 * \param name       A unique name for the object.
 * \param corner     bottom left corner in selected region of image (0-rel)
 * \param dim        dimension of region
 * \param peak       peak in region
 * \param peakResid  peak in region residual after model subtraction
 * \param RMSResid   RMS residual
 * \param fluxResid  Sum of pixel values in residual
 * \param nmodel     Number of models
 * \param models     Array of Models, steals references
 * \return the new object.
 */
ObitFitRegion* 
ObitFitRegionCreate (gchar* name, olong corner[2], olong dim[2],
		     ofloat peak, ofloat peakResid, ofloat RMSResid,
		     ofloat fluxResid, olong nmodel, ObitFitModel **models)
{
  ObitFitRegion* out;
  olong i;

  /* Create basic structure */
  out = newObitFitRegion (name);

  /*  copy info */
  out->corner[0] = corner[0];
  out->corner[1] = corner[1];
  out->dim[0]    = dim[0];
  out->dim[1]    = dim[1];
  out->peak      = peak;
  out->peakResid = peakResid;
  out->RMSResid  = RMSResid;
  out->fluxResid = fluxResid;
  out->nmodel    = nmodel;
  out->models    = g_malloc0(out->nmodel*sizeof(ObitFitModel*));

  for (i=0; i<out->nmodel; i++) out->models[i] = models[i];
    return out;
} /* end ObitFitRegionCreate */

/**
 * Generate the region name from an index.
 * Each region has the name "regnnnnnn" where nnnnnn is the 1-rel
 * \param index  index number of region
 * \return name, should be g_freeed when done
 */
gchar* ObitFitRegionName (gint indx)
{
  gchar *out=NULL;
  
  out = g_malloc0 (10);
  sprintf (out, "reg%6.6d", indx);
  return out;
  
} /* end ObitFitRegionName */

/**
 * Subtract ObitFitRegionList from its image
 * \param in     region model in image.
 * \param image  Image with attached image buffer pixel array
 * \param err    Obit Error/message stack
 */
void ObitFitRegionSubtract (ObitFitRegion* reg, ObitImage *image, ObitErr *err)
{
  olong i, j, ndim = 2, pos1[2], pos2[2];
  ofloat Cen[2], GauMod[3];
  ObitFArray *work=NULL;

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(reg, &myClassInfo));
  g_assert (ObitImageIsA(image));

  /* Make work FArray the size of the region */
  pos1[0] = reg->dim[0]; pos1[1] = reg->dim[1];
  work =  ObitFArrayCreate ("Subtract", ndim, pos1);

  /* Loop adding (neg)  model to work */
  for (i=0; i<reg->nmodel; i++) {
    Cen[0] = reg->models[i]->DeltaX;
    Cen[1] = reg->models[i]->DeltaY;
    for (j=0; j<reg->models[i]->nparm; j++) GauMod[j] = reg->models[i]->parms[j];
    GauMod[0] = reg->models[i]->parms[0]; /* ObitFArray2DEGauss from X axis */
    GauMod[1] = reg->models[i]->parms[1];
    GauMod[2] = 90.0+reg->models[i]->parms[2]*RAD2DG;  /* Fooey! */
    ObitFArray2DEGauss (work, reg->models[i]->Peak, Cen, GauMod);
  }

  /* Subtract */
  pos1[0] = reg->corner[0] + reg->dim[0]/2; 
  pos1[1] = reg->corner[1] + reg->dim[1]/2; ;
  pos2[0] = 1+reg->dim[0]/2; 
  pos2[1] = 1+reg->dim[1]/2;
  ObitFArrayShiftAdd (image->image, pos1, work, pos2, -1.0, image->image);

  work = ObitFArrayUnref(work); /* cleanup */
} /* end ObitFitRegionSubtract */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFitRegionClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFitRegionClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFitRegionClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFitRegionClassInfoDefFn (gpointer inClass)
{
  ObitFitRegionClassInfo *theClass = (ObitFitRegionClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFitRegionClassInit;
  theClass->newObit       = (newObitFP)newObitFitRegion;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFitRegionClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFitRegionGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitFitRegionCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFitRegionClear;
  theClass->ObitInit      = (ObitInitFP)ObitFitRegionInit;
  theClass->ObitFitRegionCreate = (ObitFitRegionCreateFP)ObitFitRegionCreate;

} /* end ObitFitRegionClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitFitRegionInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFitRegion *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->corner[0] = 0;
  in->corner[1] = 0;
  in->dim[0]    = 0;
  in->dim[1]    = 0;
  in->peak      = 0;
  in->peakResid = 0;
  in->RMSResid  = 0;
  in->fluxResid = 0;
  in->nmodel    = 0;
  in->models    = NULL;
 
} /* end ObitFitRegionInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitFitRegion* cast to an Obit*.
 */
void ObitFitRegionClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitFitRegion *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->models) {
    for (i=0; i<in->nmodel; i++) {
      in->models[i] = ObitFitModelUnref(in->models[i]);
    } 
    g_free (in->models);
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && (ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitFitRegionClear */

