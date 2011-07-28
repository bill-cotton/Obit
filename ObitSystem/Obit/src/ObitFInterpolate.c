/* $Id$ */
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

#include "ObitFInterpolate.h"
#include "ObitPosition.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFInterpolate.c
 * ObitFInterpolate class function definitions.
 * This class is derived from the Obit base class.
 * This class supports 1 and 2-D interpolation in ObitFArrays using 
 * Lagrange interpolation.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitFInterpolate";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitFInterpolateClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFInterpolateClassInfo myClassInfo = {FALSE};


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFInterpolateInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFInterpolateClear (gpointer in);

/** Private: Set convolution kernal */
static void SetConvKernal (ofloat Target, olong naxis, olong hwidth, 
			   ofloat *denom, olong *Start, ofloat *Kernal);

/** Private: Set Class function pointers. */
static void ObitFInterpolateClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitFInterpolate* newObitFInterpolate (gchar* name)
{
  ObitFInterpolate* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFInterpolateClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFInterpolate));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFInterpolateInit((gpointer)out);

 return out;
} /* end newObitFInterpolate */

/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name   An optional name for the object.
 * \param array  The ObitFarray to be interpolated.
 * \param desc   if nonNULL, an image descriptor to be attached to the output.
 * \param hwidth Half width of interpolation kernal (range [1,4] allowed).
 * \return the new object.
 */
ObitFInterpolate* newObitFInterpolateCreate (gchar* name, ObitFArray *array, 
					     ObitImageDesc *desc, olong hwidth)
{
  ObitFInterpolate* out;
  olong iwid, i, j;
  ofloat prod;

 /* error checks */
  g_assert (ObitFArrayIsA(array));

  /* Create/init output structure */
  out = newObitFInterpolate (name);

  /* Attach array */
  out->myArray = ObitFArrayRef(array);
  out->hwidth = MAX (1, MIN (4, hwidth));

  /* Image descriptor if given */
  if (desc) out->myDesc = ObitImageDescRef(desc);

  /* Get array info */
  out->array  = array->array;
  out->nx     = array->naxis[0];
  out->ny     = array->naxis[1];

  /* Use same thread object */
  out->thread = ObitThreadUnref(out->thread);
  out->thread = ObitThreadRef(array->thread);

  /* Init Lagrangian denominators for hwidth */
  iwid = 1 + (2*out->hwidth);
  for (j= 1; j<=iwid; j++) {
    prod = 1.0;
    for (i= 1; i<=iwid; i++) {
      if (i != j) prod = prod * (j - i);
    } 
    out->denom[j-1] = 1.0 / prod;
  } 


 return out;
} /* end newObitFInterpolateCreate */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFInterpolateGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFInterpolateClassInit();

  return (gconstpointer)&myClassInfo;
} /* end  ObitFInterpolateGetClass */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitFInterpolate* ObitFInterpolateCopy (ObitFInterpolate *in, ObitFInterpolate *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitFInterpolate(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes additions */
  out->myArray = ObitFArrayCopy(in->myArray, out->myArray, err);
  out->myDesc  = ObitImageDescCopy(in->myDesc, out->myDesc, err);
  out->array   = out->myArray->array;
  out->nx      = out->myArray->naxis[0];
  out->ny      = out->myArray->naxis[1];
  out->hwidth  = in->hwidth;
  for (i=0; i<10; i++) out->denom[i] = in->denom[i];
  
  return out;
} /* end ObitFInterpolateCopy */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
ObitFInterpolate* ObitFInterpolateClone  (ObitFInterpolate *in, ObitFInterpolate *out)
{
  const ObitClassInfo *myClass, *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Clone: ",in->name,NULL);
    out = newObitFInterpolate(outName);
    g_free(outName);
  }

  /* shallow copy any parent class members */
   myClass     = in->ClassInfo;
   ParentClass = myClass->ParentClass;
   g_assert ((ParentClass!=NULL) && (ParentClass->ObitClone!=NULL));
   ParentClass->ObitClone ((Obit*)in, (Obit*)out);

   if (!oldExist) { /* only copy ObitInfoList if just created */
     out->info = ObitInfoListUnref(out->info);
     out->info = ObitInfoListRef(in->info);
   }
     
   /* copy/set this classes additions */
   out->myArray = ObitFArrayUnref(out->myArray);
   out->myArray = ObitFArrayRef(in->myArray);
   out->myDesc  = ObitImageDescUnref(out->myDesc);
   out->myDesc  = ObitImageDescRef(in->myDesc);
   out->array = in->array;
   out->nx    = in->nx;
   out->ny    = in->ny;
   out->hwidth= in->hwidth;
   for (i=0; i<10; i++) out->denom[i] = in->denom[i];
   
  return out;
} /* end ObitFInterpolateClone */

/**
 * Replace the ObitFArray member to be interpolated.
 * \param in        The object to update
 * \param newArray  The new FArray for in
 */
void ObitFInterpolateReplace  (ObitFInterpolate *in, ObitFArray *newArray)
{
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitFArrayIsA(newArray));

  /* unreference old */
  in->myArray = ObitFArrayUnref(in->myArray);

  /* Set new */
  in->myArray = ObitFArrayRef(newArray);
  in->array   = in->myArray->array;

  /* update other data */
  in->nx = newArray->naxis[0];
  in->ny = newArray->naxis[1];
  in->xPixel = -1.0;
  in->yPixel = -1.0;
  in->xStart = -1;
  in->yStart = -1;
} /* end ObitFInterpolateReplace */

/**
 * Interpolate value at requested pixel in a plane of an n(>=2)-D array.
 * Interpolation between planes is not supported.
 * \param in    The object to interpolate
 * \param pixel Pixel location (1-rel) in planes and which plane.
 *              Should have number of dimensions equal to in.
 * \param err   Error stack if pixel not inside image.
 * \return value, magic blanked if invalid
 */
ofloat ObitFInterpolatePixel (ObitFInterpolate *in, ofloat *pixel, ObitErr *err)
{
  ofloat fblank =  ObitMagicF();
  ofloat value = fblank;
  ofloat sum, sumwt, wty, wt, prod, den, xp, yp, row[10];
  ofloat *xKernal, *yKernal, *data;
  olong i, j, k, good, xStart, yStart, iwid, indx, planeOff, iplane, iprod;
  /*gchar *routine = "ObitFInterpolatePixel";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return value;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitFArrayIsA(in->myArray));
  /* Must be inside array */
  iplane = 1;
  iprod = 1;
  for (i=0; i<in->myArray->ndim; i++) {
    if (i>1) {
      iplane *= pixel[i] * iprod;  /* How many planes deep? */
      iprod *= in->myArray->naxis[i];
    }
    if ((pixel[i]<1.0) || (pixel[i] > in->myArray->naxis[i])) {
      /* Obit_log_error(err, OBIT_Error, */
      /*	     "%s: Pixel %f, dimension %d outside of %s", */
      /*	     routine, pixel[i], i, in->name); */

      /* out of bounds - return blank */
      return fblank;
    }
  }
  if (err->error) return value;

  /* Update convolving x, y kernals as needed */
  if (pixel[0] != in->xPixel) {
    in->xPixel = pixel[0];
    SetConvKernal (in->xPixel, in->nx, in->hwidth, in->denom, &in->xStart, in->xKernal);
  }
    
  if (pixel[1] != in->yPixel) {
    in->yPixel = pixel[1];
    SetConvKernal (in->yPixel, in->ny, in->hwidth, in->denom, &in->yStart, in->yKernal);
  }

  /* Local versions of things */
  xKernal = in->xKernal;
  yKernal = in->yKernal;
  data    = in->array;
  xStart  = in->xStart;
  yStart  = in->yStart;
  iwid    = 1 + 2 * in->hwidth;

  /* Offset to start of plane */
  planeOff = (iplane-1) * in->nx * in->ny;
    
  /* Zero sums */
  sum = 0.0;
  sumwt = 0.0;
  good = 0;
 
  /* Loop over data summing values times convolving weights */
  for (j=0; j<iwid; j++) {
    wty = yKernal[j];
    for (i=0; i<iwid; i++) {
      indx = planeOff + xStart + i + ((yStart + j) * in->nx);
      wt = xKernal[i] * wty;
      if (data[indx] != fblank) {
	sumwt = sumwt + wt;
	sum = sum + data[indx] * wt;
	good = good + 1;
      } 
    }
  }

  /* normalize sum if not excessive blanking */
  if (sumwt > 0.90) {
    value = sum / sumwt;
    if (isnan(value)) g_error("INVALID VALUE, sum=%f wt=%f\n",sum, sumwt); /* DEBUG */
    return value;
  } 

  /* too much blanked data; try again  using knowledge of blanks if 
     sufficient data need at least a third of points. */
  if (good < (iwid * iwid)/3) return value;

  /* set "x" at first pixel to 1.0 */
  xp = pixel[0] - xStart;
  yp = pixel[1] - yStart;
  for (j= 1; j<=iwid; j++) { /* loop 200 */
    /* interpolate in rows */
    sum = 0.0;
    sumwt = 0.0;
    for (i=0; i<iwid; i++) { /* loop 120 */
      den = 1.0;
      prod = 1.0;
      for (k=0; k<iwid; k++) { /* loop 110 */
	indx = planeOff + xStart + k + ((yStart + j) * in->nx);
	if (data[indx] != fblank) {
	  if (i != k) {
	    den = den * (i - k);
	    prod = prod * (xp - k);
	  } 
	} 
      } /* end loop  L110: */;
 
     indx = planeOff + xStart + i + ((yStart + j) * in->nx);

      /* accumulate */
      if (data[indx] != fblank) {
	if (abs (den) > 1.0e-10) {
	  wt = prod / den;
	} else {
	  wt = 0.0;
	} 
	sumwt = sumwt + wt;
	sum = sum + wt * data[indx];
      } 
    } /* end loop  L120: */

    /* interpolate column value */
    if (sumwt > 0.5) {
      row[j-1] = sum / sumwt;
    } else {
      row[j-1] = fblank;
    } 
  } /* end loop  L200: */

  /* interpolate in column */
  sum = 0.0;
  sumwt = 0.0;
  for (i=0; i<iwid; i++) { /* loop 220 */
    den = 1.0;
    prod = 1.0;
    for (k=0; k<iwid; k++) { /* loop 210 */
      if (row[k] != fblank) {
	if (i != k) {
	  den = den * (i - k);
	  prod = prod * (yp - k);
	} 
      } 
    } /* end loop  L210: */
    
    /* accumulate */
    if (row[i-1] != fblank) {
      if (abs (den) > 1.0e-10) {
	wt = prod / den;
      } else {
	wt = 0.0;
      } 
      sumwt = sumwt + wt;
      sum = sum + wt * row[i-1];
    } 
  } /* end loop  L220: */

  /* interpolate value */
  if (sumwt > 0.5) {
    value = sum / sumwt;
  } else {
    value = fblank;
  } 

  /* this seems to prevent an optimizer related bug */
  if (isnan(value)) {
    for (i=0; i<iwid; i++) fprintf (stderr, "row %d %f\n",i, row[i]);
    g_error("INVALID VALUE2, sum=%f wt=%f\n",sum, sumwt); /* DEBUG */
  }
  return value;
} /* end ObitFInterpolatePixel */

/**
 * Interpolate value at requested pixel in 1-D array.
 * \param in    The object to interpolate
 * \param pixel Pixel location (1-rel) in array
 * \return value, blanked if invalid
 */
ofloat ObitFInterpolate1D (ObitFInterpolate *in, ofloat pixel)
{
  ofloat fblank =  ObitMagicF();
  ofloat value = fblank;
  ofloat sum, sumwt, wt;
  ofloat *xKernal, *data;
  olong i, good, xStart, iwid, indx;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitFArrayIsA(in->myArray));

  /* Must be inside 1-D array */
  if (in->myArray->ndim>2) return fblank; /* too many dimensions */
  if (pixel<0.5) return fblank; /* not in array */
  if (pixel>(in->myArray->naxis[0]+0.5)) return fblank; /* not in array */

  /* Update convolving kernal as needed */
  if (pixel != in->xPixel) {
    in->xPixel = pixel;
    SetConvKernal (in->xPixel, in->nx, in->hwidth, in->denom, &in->xStart, in->xKernal);
  }
    
  /* Local versions of things */
  xKernal = in->xKernal;
  data    = in->array;
  xStart  = in->xStart;
  iwid    = 1 + 2 * in->hwidth;

  /* Zero sums */
  sum   = 0.0;
  sumwt = 0.0;
  good  = 0;
 
  /* Loop over data summing values times convolving weights */
  for (i=0; i<iwid; i++) {
    indx = xStart + i;
    wt = xKernal[i];
    if (data[indx] != fblank) {
      sumwt = sumwt + wt;
      sum = sum + data[indx] * wt;
      good = good + 1;
    } 
  }


  /* normalize sum if not excessive blanking */
  if (sumwt > 0.20) {
    value = sum / sumwt;
    return value;
  } 

  /* No luck - return fblank */
  return fblank;
} /* end ObitFInterpolate1D */

/**
 * Interpolate value at requested coordinate in array.
 * The object must have an image descriptor to allow determing
 * pixel coordinates.
 * Interpolation between planes is not supported.
 * \param in    The object to interpolate
 * \param coord Coordinate value in plane and which plane.
 *              Should have number of dimensions equal to in.
 * \param err   Error stack is pixel not inside image.
 * \return value, blanked if invalid
 */
ofloat ObitFInterpolatePosition (ObitFInterpolate *in, odouble *coord, ObitErr *err)
{
  ofloat pixel[IM_MAXDIM], fblank = ObitMagicF();
  ObitImageDesc *desc = NULL;
  gchar *routine = "ObitFInterpolatePosition";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return fblank;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitImageDescIsA(in->myDesc));

  /* convert to pixel */
  desc = in->myDesc;
  if (ObitPositionXYpix(coord, desc, pixel)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: error converting coordinates", routine);
    return fblank;
  }
 
  /* interpolate */
  return ObitFInterpolatePixel (in, pixel, err);
} /* end ObitFInterpolatePosition */

/**
 * Interpolate value at requested offset from reference position in array. 
 * The object must have an image descriptor to allow determing
 * pixel coordinates.
 * Interpolation between planes is not supported.
 * \param in    The object to interpolate
 * \param off   Coordinate offset in plane
 * \param err   Error stack is pixel not inside image.
 * \return value, blanked if invalid
 */
ofloat ObitFInterpolateOffset (ObitFInterpolate *in, odouble *off, ObitErr *err)
{
  ofloat pixel[IM_MAXDIM], fblank = ObitMagicF();
  olong i;
  ObitImageDesc *desc = NULL;

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return fblank;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitImageDescIsA(in->myDesc));

  /* convert to pixel */
  desc = in->myDesc;
   /* Just linear transformation needed*/
  for (i=0; i<in->myArray->ndim; i++) {
    pixel[i] = desc->crpix[i] + (off[i] / desc->cdelt[i]);
  }
  /* debug
  fprintf(stderr,"pixel %f %f\n", pixel[0], pixel[1]); */

  /* interpolate */
  return ObitFInterpolatePixel (in, pixel, err);
} /* end ObitFInterpolateOffset */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFInterpolateClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFInterpolateClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFInterpolateClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFInterpolateClassInfoDefFn (gpointer inClass)
{
  ObitFInterpolateClassInfo *theClass = (ObitFInterpolateClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFInterpolateClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFInterpolateClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFInterpolateGetClass;
  theClass->newObit       = (newObitFP)newObitFInterpolate;
  theClass->ObitCopy      = (ObitCopyFP)ObitFInterpolateCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitFInterpolateClone;
  theClass->ObitClear     = (ObitClearFP)ObitFInterpolateClear;
  theClass->ObitInit      = (ObitInitFP)ObitFInterpolateInit;

} /* end ObitFInterpolateClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param in Pointer to the object to initialize.
 */
void ObitFInterpolateInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFInterpolate *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread      = newObitThread();
  in->info        = newObitInfoList(); 
  in->myArray     = NULL;
  in->array       = NULL;
  in->xPixel     = -1.0e20;
  in->yPixel     = -1.0e20;

} /* end ObitFInterpolateInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  in Pointer to the object to deallocate.
 *           Actually it should be an ObitFInterpolate* cast to an Obit*.
 */
void ObitFInterpolateClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFInterpolate *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread      = ObitThreadUnref(in->thread);
  in->info        = ObitInfoListUnref(in->info);
  in->myArray     = ObitFArrayUnref(in->myArray);
  in->myDesc      = ObitImageDescUnref(in->myDesc);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitFInterpolateClear */


/**
 * Set Lagrangian interpolation kernal taking into account ends of the grid.
 * \param  Pixel Which is the desired pixel (1-rel)?
 * \param  naxis  Number of pixels on axis 
 * \param  hwidth Half width of convolution kernal
 * \param  denom  Reciprocals of Lagrangian denominators
 * \param  Start  [out] first pixels in array for convolution
 * \param  Kernal [out] convolving kernal.
 */
static void SetConvKernal (ofloat Pixel, olong naxis, olong hwidth, 
			   ofloat *denom, olong *Start, ofloat *Kernal)
{
  ofloat prod, xx;
  olong ipos, i, j, cen, iwid;

  /* fractional pixel */
  ipos = Pixel + 0.5;
  iwid = hwidth*2 + 1;

  /* set first pixel */
  cen = ipos - hwidth;
  cen = MAX (1, MIN (cen, (naxis-iwid+1)));
  /* make 0 rel */
  cen = cen - 1;
  *Start = cen; /* returned version */

  /* set "x" at first pixel to 1.0 */
  xx = Pixel - cen;

  /* compute interpolating kernal */
  for (j= 0; j<iwid; j++) { /* loop 50 */
    prod = denom[j];
    for (i= 0; i<iwid; i++) { /* loop 30 */
      if (i != j) prod = prod * (xx - (i+1));
    } /* end loop  L30:  */;
    Kernal[j] = prod;
  } /* end loop  L50:  */;
} /* end SetConvKernal */

