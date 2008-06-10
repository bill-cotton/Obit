/* $Id: ObitDConCleanPxHist.c,v 1.6 2007/05/10 21:18:08 bcotton Exp $ */
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

#include "ObitDConCleanPxHist.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanPxHist.c
 * ObitDConCleanPxHist class function definitions.
 * This class determines the pixel histogram of an image.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanPxHist";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitDConCleanPxHistClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanPxHistClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/
/** Number of cells in histogram. */
olong histSize = 8192;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanPxHistInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanPxHistClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConCleanPxHistClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanPxHist* newObitDConCleanPxHist (gchar* name)
{
  ObitDConCleanPxHist* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxHistClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanPxHist));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanPxHistInit((gpointer)out);

 return out;
} /* end newObitDConCleanPxHist */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanPxHistGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxHistClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanPxHistGetClass */

/**
 * Make a deep copy of an ObitDConCleanPxHist.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanPxHist* 
ObitDConCleanPxHistCopy  (ObitDConCleanPxHist *in, ObitDConCleanPxHist *out, 
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
    out = newObitDConCleanPxHist(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->ncell = in->ncell;
  out->histMax = in->histMax;
  out->histMin = in->histMin;
  if ((out->hist) && (ObitMemValid (out->hist))) 
    out->hist = ObitMemFree (out->hist);
  out->hist = ObitMemAlloc0Name(out->ncell*sizeof(olong),"Pixel histogram");
  for (i=0; i<in->ncell; i++) out->hist[i] = in->hist[i];

  return out;
} /* end ObitDConCleanPxHistCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanPxHist similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanPxHistClone  (ObitDConCleanPxHist *in, ObitDConCleanPxHist *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;

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
  out->histMax = in->histMax;
  out->histMin = in->histMin;
  out->ncell = in->ncell;
  if ((out->hist) && (ObitMemValid (out->hist))) 
    out->hist = ObitMemFree (out->hist);
  out->hist = ObitMemAlloc0Name(out->ncell*sizeof(olong),"Pixel histogram");
  for (i=0; i<in->ncell; i++) out->hist[i] = in->hist[i];

} /* end ObitDConCleanPxHistClone */

/**
 * Update histogram using given image and window
 * \param in     The Pixel histogram object 
 * \param field  Which field? (1-rel)
 * \param plane  1-rel indices on dimensions 3-?
 * \param mosaic Image Mosaic with images
 * \param window Corresponding windows in mosaic images.
 *        Only pixels inside of the CLEAN window are used.
 * \param err Obit error stack object.
 */
void ObitDConCleanPxHistUpdate (ObitDConCleanPxHist *in, olong field, 
				olong *plane, ObitImageMosaic *mosaic,
				ObitDConCleanWindow *window, 
				ObitErr *err)
{
  ObitIOCode retCode;
  ObitImage *image=NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong i, ix, iy, icell, nx, ny, pos[2];
  ofloat *data, tmax, tmin, tfact;
  gboolean *mask=NULL;
  gchar *routine = "ObitDConCleanPxHistUpdate";

  /* error checks */
  if (err->error) return;

  if ((field<=0) || (field>mosaic->numberImages)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1-%d in %s",
                   routine, field, mosaic->numberImages, mosaic->name);
      return;
  }

  /* Which image? */
  image = mosaic->images[field-1];

  /* Set output to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = 1;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = plane[i];
  ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
  trc[0] = trc[1] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = plane[i];
  ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);
 
  retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  retCode = ObitImageRead (image, image->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Create histogram */
  in->ncell =  histSize;  /* Pick some number */
  /* allocate or reallocate */
  if (in->hist) in->hist = ObitMemRealloc(in->hist, in->ncell*sizeof(olong));
  else in->hist = ObitMemAlloc0Name(in->ncell*sizeof(olong),"image histogram");
  for (ix=0; ix<in->ncell; ix++) in->hist[ix] = 0;  /* initialize */

  /* pointer to data */
  pos[0] = pos[1] = 0;
  data = ObitFArrayIndex(image->image, pos);

  /* Loop over image getting max, min values */
  nx = image->myDesc->inaxes[0];
  ny = image->myDesc->inaxes[1];
  tmax = -1.0e20;
  tmin = 1.2e20;
  for (iy=0; iy<ny; iy++) {
    /* Get window mask */
    if (ObitDConCleanWindowRow(window, field, iy+1, &mask, err)) {
      for (ix=0; ix<nx; ix++) {
	if (mask[ix]) {
	  tmax = MAX (tmax, fabs(data[ix]));
	  tmin = MIN (tmin, data[ix]);
	}
      }
    }
    data += nx;
  }

  /* Find anything? */
  if (tmax<0.0) {
    tmax = tmin = 0.0;
  }

  /* save extrema */
  in->histMax = tmax;
  in->histMin = tmin;

  /* Now compute histogram */
  if (tmax != tmin) tfact = (ofloat)(in->ncell-1) / (tmax - tmin);
  else tfact = 1.0;
  pos[0] = pos[1] = 0;
  data = ObitFArrayIndex(image->image, pos);
  for (iy=0; iy<ny; iy++) {
    /* Get window mask */
    if (ObitDConCleanWindowRow(window, field, iy+1, &mask, err)) {
      for (ix=0; ix<nx; ix++) {
	if (mask[ix]) {
	  icell = tfact * (fabs(data[ix]) - tmin) + 0.5;
	  icell = MIN (icell, in->ncell-1);
	  in->hist[icell]++;    /* count it */
	}
      }
    }
    data += nx;
  }

  retCode = ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Free Image array? */
  image->image = ObitFArrayUnref(image->image);

  /* Finish histogram, want cumulative distribution */
  for (icell=in->ncell-2; icell>=0; icell--) {
    in->hist[icell] += in->hist[icell+1];
  }

  /* Cleanup */
  if ((mask) && (ObitMemValid (mask))) mask = ObitMemFree (mask);

} /* end ObitDConCleanPxHistUpdate */

/**
 * Tell the number of pixels with abs. value larger than value.
 * Only pixels inside of the CLEAN window are used.
 * \param in    The Beam histogram object 
 * \param value The value of interest
 * \param err   Obit error stack object.
 * \return number of pixels of abs value > value
 */
olong ObitDConCleanPxHistNumber (ObitDConCleanPxHist *in, ofloat value,
				 ObitErr *err)
{
  olong out = 0;
  olong icell;
  gchar *routine = "ObitDConCleanPxHistPeak";

  /* error checks */
  if (err->error) return out;

  /* Do we have a histogram? */
  if ((in->ncell<=0) || (!in->hist)) {
    /* No histogram! */
    Obit_log_error(err, OBIT_Error,"%s: NO Beam histogram", routine);
    return out;
  }

  /* look up in table - closest value */
  icell = ((ofloat)(in->ncell-1) / (in->histMax - in->histMin)) * 
    (fabs(value) - in->histMin) + 0.5;
  icell = MIN (icell, in->ncell-1);
  icell = MAX (icell, 0);
  out = in->hist[icell];

  return out;
} /* end ObitDConCleanPxHistPeak */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanPxHistClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanPxHistClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanPxHistClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanPxHistClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanPxHistClassInfo *theClass = (ObitDConCleanPxHistClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanPxHistClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanPxHistClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanPxHistGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanPxHist;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanPxHistCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanPxHistClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanPxHistInit;

} /* end ObitDConCleanPxHistClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanPxHistInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanPxHist *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
    in->hist = NULL;

} /* end ObitDConCleanPxHistInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanPxHist* cast to an Obit*.
 */
void ObitDConCleanPxHistClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanPxHist *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if ((in->hist) && (ObitMemValid (in->hist))) in->hist = ObitMemFree (in->hist);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanPxHistClear */


