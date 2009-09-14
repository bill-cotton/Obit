/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2009                                          */
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

#include "ObitDConCleanBmHist.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanBmHist.c
 * ObitDConCleanBmHist class function definitions.
 * This class determines the maximum beam sidelobe exterior to a given 
 * x or y distance.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanBmHist";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitDConCleanBmHistClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanBmHistClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanBmHistInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanBmHistClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConCleanBmHistClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanBmHist* newObitDConCleanBmHist (gchar* name)
{
  ObitDConCleanBmHist* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanBmHistClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanBmHist));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanBmHistInit((gpointer)out);

 return out;
} /* end newObitDConCleanBmHist */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanBmHistGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanBmHistClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanBmHistGetClass */

/**
 * Make a deep copy of an ObitDConCleanBmHist.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanBmHist* 
ObitDConCleanBmHistCopy  (ObitDConCleanBmHist *in, ObitDConCleanBmHist *out, 
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
    out = newObitDConCleanBmHist(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->ncell = in->ncell;
  if ((out->hist) && (ObitMemValid (out->hist))) 
    out->hist = ObitMemFree (out->hist);
  out->hist = ObitMemAlloc0Name(out->ncell*sizeof(ofloat),"Beam histogram");
  for (i=0; i<in->ncell; i++) out->hist[i] = in->hist[i];

  return out;
} /* end ObitDConCleanBmHistCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanBmHist similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanBmHistClone  (ObitDConCleanBmHist *in, ObitDConCleanBmHist *out, ObitErr *err)
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
  out->ncell = in->ncell;
  if ((out->hist) && (ObitMemValid (out->hist))) 
    out->hist = ObitMemFree (out->hist);
  out->hist = ObitMemAlloc0Name(out->ncell*sizeof(ofloat),"Beam histogram");
  for (i=0; i<in->ncell; i++) out->hist[i] = in->hist[i];

} /* end ObitDConCleanBmHistClone */

/**
 * Creates an ObitDConCleanBmHist 
 * \param name  An optional name for the object.
 * \param Beam  from which to create object
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitDConCleanBmHist* 
ObitDConCleanBmHistCreate (gchar* name, ObitImage *Beam, ObitErr *err)
{
  ObitDConCleanBmHist* out=NULL;
  olong plane[5] = {1,1,1,1,1};
  gchar *routine = "ObitDConCleanBmHistCreate";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitImageIsA(Beam));

  /* Create basic structure */
  out = newObitDConCleanBmHist (name);

  /* Update with beam */
  ObitDConCleanBmHistUpdate (out, Beam, plane, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  return out;
} /* end ObitDConCleanBmHistCreate */

/**
 * Update histogram using given beam image
 * Distance measure is MAX (delta_x, delta_y)
 * \param in   The Beam histogram object 
 * \param plane  1-rel indices on dimensions 3-?
 * \param Beam The Beam image
 * \param err Obit error stack object.
 */
void ObitDConCleanBmHistUpdate (ObitDConCleanBmHist *in, ObitImage *Beam,
				olong *plane, ObitErr *err)
{
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong i, ix, iy, icenx, iceny, irad, nx, ny, pos[2];
  ofloat *data, beamMax, last=0.0, fblank = ObitMagicF();
  gchar *routine = "ObitDConCleanBmHistUpdate";

  /* Set output to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = blc[2] = blc[3] = blc[4] = blc[5] = 1;
  trc[0] = trc[1] = trc[2] = trc[3] = trc[4] = trc[5] = 0;
  /* multiplane? */
  if (Beam->myDesc->inaxes[2]>1) {
    for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = plane[i];
    for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = plane[i];
  }
  ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (Beam->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (Beam->info, "IOBy", OBIT_long, dim, &IOsize, err);
 
  retCode = ObitImageOpen (Beam, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageRead (Beam, Beam->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Create histogram */
  in->ncell = MAX (Beam->myDesc->inaxes[0], Beam->myDesc->inaxes[1]) / 2;
  /* allocate or reallocate */
  if (in->hist) in->hist = ObitMemRealloc(in->hist, in->ncell*sizeof(ofloat));
  else in->hist = ObitMemAlloc0Name(in->ncell*sizeof(ofloat),"Beam histogram");
  for (ix=0; ix<in->ncell; ix++) in->hist[ix] = 0.0;  /* initialize */

  /* pointer to data */
  pos[0] = 0; pos[1] = 5;
  data = ObitFArrayIndex(Beam->image, pos);

  /* Center pixel - make 0-rel */
  icenx =  Beam->myDesc->crpix[0] - Beam->myDesc->xPxOff - 0.5;
  iceny =  Beam->myDesc->crpix[1] - Beam->myDesc->yPxOff - 0.5;

  /* Loop over image - ignore outer pixels */
  nx = Beam->myDesc->inaxes[0];
  ny = Beam->myDesc->inaxes[1];
  for (iy=5; iy<ny-5; iy++) {
    for (ix=5; ix<nx-5; ix++) {
      if (data[ix]!=fblank) {
	irad = MAX (abs(ix-icenx), abs (iy-iceny));
	if (irad>=in->ncell) irad = in->ncell - 1;
	in->hist[irad] = MAX (in->hist[irad], data[ix]);
      }
    }
    data += nx;
  }

  retCode = ObitImageClose (Beam, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Free Image array? */
  Beam->image = ObitFArrayUnref(Beam->image);

  /* Finish histogram, currently have maximum value at each distance */
  beamMax = in->hist[in->ncell-1];
  for (ix=in->ncell-1; ix>=0; ix--) {
    /* Keep it real */
    if (ix>0) in->hist[ix] = MIN (0.99, in->hist[ix]);
    if (in->hist[ix]>beamMax) beamMax = in->hist[ix];
    else in->hist[ix]= beamMax;
  }
  /* Fill with last non zero value to end */
  for (ix=0; ix<in->ncell; ix++) {
    if (in->hist[ix]!=0.0) last = in->hist[ix];
    else if (in->hist[ix]==0.0) in->hist[ix] = last;
  }

} /* end ObitDConCleanBmHistUpdate */

/**
 * Determine max. abs. value of sidelobe exterior to a given radius.
 * \param in     The Beam histogram object 
 * \param radius How far from center in cells?
 *   (Not radius in the usual sense but a beam patch sense).
 * \param err    Obit error stack object.
 * \return the absolute values of the maximum exterior sidelobe.
 */
ofloat ObitDConCleanBmHistPeak (ObitDConCleanBmHist *in, olong radius,
				ObitErr *err)
{
  ofloat out = 0.0;
  olong irad;
  gchar *routine = "ObitDConCleanBmHistPeak";

  /* error checks */
  if (err->error) return out;

  /* Do we have a histogram? */
  if ((in->ncell<=0) || (!in->hist)) {
    /* No histogram! */
    Obit_log_error(err, OBIT_Error,"%s: NO Beam histogram", routine);
    return out;
  }

  /* look up in table - closest value */
  irad = radius;
  if (irad<in->ncell) out = in->hist[irad];
  else out = in->hist[in->ncell-1];

  return out;
} /* end ObitDConCleanBmHistPeak */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanBmHistClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanBmHistClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanBmHistClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanBmHistClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanBmHistClassInfo *theClass = (ObitDConCleanBmHistClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanBmHistClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanBmHistClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanBmHistGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanBmHist;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanBmHistCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitRef       = (ObitRefFP)ObitRef;
  theClass->ObitUnref     = (ObitUnrefFP)ObitUnref;
  theClass->ObitIsA       = (ObitIsAFP)ObitIsA;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanBmHistClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanBmHistInit;

} /* end ObitDConCleanBmHistClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanBmHistInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanBmHist *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
    in->hist = NULL;

} /* end ObitDConCleanBmHistInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanBmHist* cast to an Obit*.
 */
void ObitDConCleanBmHistClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanBmHist *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
    if ((in->hist) && (ObitMemValid (in->hist))) in->hist = ObitMemFree (in->hist);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanBmHistClear */


