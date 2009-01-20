/* $Id$ */
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

#include "ObitCInterpolate.h"
#include "ObitPosition.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCInterpolate.c
 * ObitCInterpolate class function definitions.
 * This class is derived from the Obit base class.
 * This class supports 1 and 2-D interpolation in ObitCArrays using 
 * Lagrange interpolation.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitCInterpolate";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitCInterpolateClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitCInterpolateClassInfo myClassInfo = {FALSE};


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitCInterpolateInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitCInterpolateClear (gpointer in);

/** Private: Compute convolution kernal table*/
static void InitConvKernal (ObitCInterpolate* in);

/** Private: Set convolution kernal from table*/
static void SetConvKernal (ObitCInterpolate* in, ofloat Pixel, 
			   olong n, olong *Start, olong *Nterm, 
			   ofloat **Kernal);

/** Private: Calculate convolution kernal */
static void CalcConvKernal (ObitCInterpolate* in, ofloat Pixel, 
			   olong n, olong *Start, olong *Nterm, 
			   ofloat *Kernal);

/** Private: Set Class function pointers. */
static void ObitCInterpolateClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitCInterpolate* newObitCInterpolate (gchar* name)
{
  ObitCInterpolate* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitCInterpolateClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitCInterpolate));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitCInterpolateInit((gpointer)out);

 return out;
} /* end newObitCInterpolate */

/**
 * Constructor.
 * Initializes class if needed on first call.
 * This class is intended for use interpolating complex values
 * in a UV grid which is the Fourier transform of an image.
 * This is potentially modified by OSX, OSY and herm.
 * Only 1 and 2D arrays are handled.
 * Note: The FFTW convention for halfplane complex images is different 
 * from AIPS, the "short" = n/2+1 axis is the first one rather than the second.
 * \param name   An optional name for the object.
 * \param array  The ObitCArray to be interpolated.
 * \param desc   if nonNULL, an image descriptor describing the image
 *               whose Fourier transform is to be interpolated.
 * \param OSX    Oversampling factor in X/U, use 1.0 for unpadded image
 * \param OSY    Oversampling factor in Y/V 
 * \param numConjCol Number of conjugate (neg V) columns in array, 
 *               should be at least hwidth
 * \param hwidth Half width of interpolation kernal (range [1,8] allowed).
 * \return the new object.
 */
ObitCInterpolate* newObitCInterpolateCreate (gchar* name, ObitCArray *array, 
					     ObitImageDesc *desc, ofloat OSX, ofloat OSY, 
					     olong numConjCol, olong hwidth, ObitErr *err)
{
  ObitCInterpolate* out=NULL;
  gchar *routine = "newObitCInterpolateCreate";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitCArrayIsA(array));
  if (numConjCol<hwidth) { /* check numConjCol */
    Obit_log_error(err, OBIT_Error,
		   "%s: numConjCol %d < hwidth (= %d) ", 
		   routine, numConjCol, hwidth); 
    return out;
  }

  /* Create/init output structure */
  out = newObitCInterpolate (name);

  /* Attach array */
  out->myArray = ObitCArrayRef(array);
  out->hwidth = MAX (1, MIN (16, hwidth));
  out->numConjCol = numConjCol;

  /* Get array info */
  out->array  = array->array;
  out->nx     = array->naxis[0];
  out->ny     = array->naxis[1];

  /* Use same thread object */
  out->thread = ObitThreadUnref(out->thread);
  out->thread = ObitThreadRef(array->thread);

  /* Ingest Image descriptor if given */
  if (desc) {
    out->myDesc = ObitImageDescCopy(desc, out->myDesc, err);
    if (err->error) Obit_traceback_val (err, routine, out->name, out);
    /* Replace with values of Fourier transform version, note, 
       this presumes that the image is Hermitian with numConjCol 
       columns added on the negative "U" side and NOT transposed */
    out->myDesc->crval[0] = 0.0;
    out->myDesc->crval[1] = 0.0;
    out->myDesc->crpix[0] = numConjCol;
    out->myDesc->crpix[1] = out->ny/2;
    out->myDesc->cdelt[0] =  
      1.0 / (-DG2RAD * OSY * (desc->cdelt[desc->jlocr])*desc->inaxes[desc->jlocr]);
    out->myDesc->cdelt[1] = 
      1.0 / (-DG2RAD * OSX * (desc->cdelt[desc->jlocd])*desc->inaxes[desc->jlocd]);
  }

  /* Initialize the interpolation */
  InitConvKernal (out);

 return out;
} /* end newObitCInterpolateCreate */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitCInterpolateGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitCInterpolateClassInit();

  return (gconstpointer)&myClassInfo;
} /* end  ObitCInterpolateGetClass */

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
ObitCInterpolate* ObitCInterpolateCopy (ObitCInterpolate *in, ObitCInterpolate *out, ObitErr *err)
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
    out = newObitCInterpolate(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes additions */
  out->myArray = ObitCArrayCopy(in->myArray, out->myArray, err);
  out->myKernal= ObitFArrayCopy(in->myKernal, out->myKernal, err);
  out->myDesc  = ObitImageDescCopy(in->myDesc, out->myDesc, err);
  out->array   = out->myArray->array;
  out->nx      = out->myArray->naxis[0];
  out->ny      = out->myArray->naxis[1];
  out->numConjCol  = in->numConjCol;
  out->kernalSpace = in->kernalSpace;
  out->numKTab = in->numKTab;
  out->xPixel  = in->xPixel;
  out->yPixel  = in->yPixel;
  out->yStart  = in->yStart;
  out->xNterm  = in->xNterm;
  out->yNterm  = in->yNterm;
  out->hwidth  = in->hwidth;
  for (i=0; i<40; i++) out->denom[i] = in->denom[i];
  
  return out;
} /* end ObitCInterpolateCopy */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
ObitCInterpolate* ObitCInterpolateClone  (ObitCInterpolate *in, ObitCInterpolate *out)
{
  const ObitClassInfo *myClass, *ParentClass;
  gboolean oldExist;
  olong i;
  ObitErr *err= NULL;
  gchar *outName;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Clone: ",in->name,NULL);
    out = newObitCInterpolate(outName);
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
   err = newObitErr();
   out->myArray = ObitCArrayUnref(out->myArray);
   out->myArray = ObitCArrayRef(in->myArray);
   out->myKernal= ObitFArrayUnref(out->myKernal);
   out->myKernal= ObitFArrayCopy(in->myKernal, out->myKernal, err);
   out->myDesc  = ObitImageDescUnref(out->myDesc);
   out->myDesc  = ObitImageDescCopy(in->myDesc, out->myDesc, err);
   out->array   = out->myArray->array;
   out->nx      = out->myArray->naxis[0];
   out->ny      = out->myArray->naxis[1];
   out->numConjCol  = in->numConjCol;
   out->kernalSpace = in->kernalSpace;
   out->numKTab   = in->numKTab;
   out->xPixel    = in->xPixel;
   out->yPixel    = in->yPixel;
   out->yStart    = in->yStart;
   out->xNterm    = in->xNterm;
   out->yNterm    = in->yNterm;
   out->hwidth= in->hwidth;
   for (i=0; i<40; i++) out->denom[i] = in->denom[i];
   ObitErrLog(err); /* Print any error messages */
   err = ObitErrUnref (err);
   
  return out;
} /* end ObitCInterpolateClone */

/**
 * Replace the ObitCArray member to be interpolated.
 * \param in        The object to update
 * \param newArray  The new CArray for in
 */
void ObitCInterpolateReplace  (ObitCInterpolate *in, ObitCArray *newArray)
{
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitCArrayIsA(newArray));

  /* unreference old */
  in->myArray = ObitCArrayUnref(in->myArray);

  /* Set new */
  in->myArray = ObitCArrayRef(newArray);
  in->array   = in->myArray->array;

  /* update other data */
  in->nx = newArray->naxis[0];
  in->ny = newArray->naxis[1];
  in->xPixel = -1.0;
  in->yPixel = -1.0;
  in->xStart = -1;
  in->yStart = -1;
} /* end ObitCInterpolateReplace */

/**
 * Interpolate value at requested pixel in a plane of an n(>=2)-D array.
 * Interpolation between planes is not supported.
 * \param in    The object to interpolate
 * \param pixel Pixel location (1-rel) in planes and which plane.
 *              Should have number of dimensions equal to in.
 * \param out   Complex interpolated value as (real,Imag), magic value blanked
 * \param err   Error stack if pixel not inside image.
 */
void ObitCInterpolatePixel (ObitCInterpolate *in, ofloat *pixel, ofloat out[2], 
			    ObitErr *err)
{
  ofloat fblank =  ObitMagicF();
  ofloat sumR, sumI, sumwt, wty, wt;
  ofloat *xKernal, *yKernal, *data;
  olong i, j, good, xStart, yStart, iwidX, iwidY, indx, planeOff, iplane, iprod;
  /*gchar *routine = "ObitCInterpolatePixel";*/

  /* error checks */
  out[0] = fblank; out[1] = fblank;  /* init result to blanked */
  /* g_assert(ObitErrIsA(err));*/
  if (err->error) return;
  /* g_assert (ObitIsA(in, &myClassInfo));*/
  /* g_assert (ObitCArrayIsA(in->myArray));*/

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
      return;
    }
  }
  /*if (err->error) return;*/

  /* Update convolving x, y kernals as needed */
  if (pixel[0] != in->xPixel) {
    in->xPixel = pixel[0];
    SetConvKernal (in, pixel[0], in->nx, &in->xStart, &in->xNterm, &in->xKernal);
   }
    
  if (pixel[1] != in->yPixel) {
    in->yPixel = pixel[1];
    SetConvKernal (in, pixel[1], in->ny, &in->yStart, &in->yNterm, &in->yKernal);
  }

  /* Local versions of things */
  xKernal = in->xKernal;
  yKernal = in->yKernal;
  data    = in->array;
  xStart  = in->xStart;
  yStart  = in->yStart;
  iwidX   = in->xNterm;
  iwidY   = in->yNterm;

  /* Offset to start of plane */
  planeOff = (iplane-1) * in->nx * in->ny;
    
  /* Zero sums */
  sumR = sumI = 0.0;
  sumwt = 0.0;
  good = 0;
 
  /* Loop over data summing values times convolving weights */
  for (j= 0; j<iwidY; j++) {
    wty = yKernal[j];
    for (i= 0; i<iwidX; i++) {
      indx = 2*(planeOff + xStart + i + ((yStart + j) * in->nx));
      wt = xKernal[i] * wty;
      if (data[indx] != fblank) {
	sumwt += wt;
	sumR  += data[indx]   * wt;
	sumI  += data[indx+1] * wt;
	good++;
      } 
    }
  }

  /* normalize sum if not excessive blanking */
  if (sumwt > 0.5) {
    out[0] = sumR / sumwt;
    out[1] = sumI / sumwt;
  } else {
    out[0] = fblank;
    out[1] = fblank;
  } 

  return;
} /* end ObitCInterpolatePixel */

/**
 * Interpolate value at requested pixel in 1-D array.
 * \param in    The object to interpolate
 * \param pixel Pixel location (1-rel) in array
 * \param out   Complex interpolated value as (real,Imag), magic value blanked
 */
void ObitCInterpolate1D (ObitCInterpolate *in, ofloat pixel, ofloat out[2])
{
  ofloat fblank =  ObitMagicF();
  ofloat sumR, sumI, sumwt, wt;
  ofloat *xKernal, *data;
  olong i, good, xStart, iwid, indx;

  /* error checks */
  out[0] = fblank; out[1] = fblank;  /* init result to blanked */
  /*g_assert (ObitIsA(in, &myClassInfo));*/
  /*g_assert (ObitCArrayIsA(in->myArray));*/

  /* Must be inside 1-D array */
  if (in->myArray->ndim>2) return; /* too many dimensions */
  if (pixel<0.5) return; /* not in array */
  if (pixel>(in->myArray->naxis[0]+0.5)) return; /* not in array */

  /* Update convolving kernal as needed */
  if (pixel != in->xPixel) {
    in->xPixel = pixel;
    SetConvKernal (in, pixel, in->nx, &in->xStart, &in->xNterm, &in->xKernal);
  }
    
  /* Local versions of things */
  xKernal = in->xKernal;
  data    = in->array;
  xStart  = in->xStart;
  iwid    = in->xNterm;

  /* Zero sums */
  sumR = sumI = 0.0;
  sumwt = 0.0;
  good  = 0;
 
  /* Loop over data summing values times convolving weights */
  for (i=0; i<iwid; i++) {
    indx = xStart + i;
    wt = xKernal[i];
    if (data[indx] != fblank) {
      sumwt += wt;
      sumR += data[2*(indx)] * wt;
      sumI += data[1+2*(indx)] * wt;
      good = good + 1;
    } 
  }


  /* normalize sum if not excessive blanking */
  if (sumwt > 0.20) {
    out[0] = sumR / sumwt;
    out[1] = sumI / sumwt;
    return;
  } 

  /* No luck - return fblank */
  out[0] = fblank;
  out[1] = fblank;
  return;
} /* end ObitCInterpolate1D */

/**
 * Interpolate value at requested coordinate in array.
 * The object must have an image descriptor to allow determing
 * pixel coordinates and the coordinates are assumed linear.
 * Interpolation between planes is not supported.
 * \param in    The object to interpolate
 * \param coord Coordinate value in plane and which plane.
 *              Should have number of dimensions equal to in.
 * \param out   Complex interpolated value as (real,Imag), magic value blanked
 * \param err   Error stack is pixel not inside image.
 */
void ObitCInterpolatePosition (ObitCInterpolate *in, odouble *coord, ofloat out[2], 
			       ObitErr *err)
{
  ofloat pixel[IM_MAXDIM], fblank = ObitMagicF();
  /*gchar *routine = "ObitCInterpolatePosition";*/

  /* error checks */
  out[0] = fblank; out[1] = fblank;  /* init result to blanked */
  /*g_assert(ObitErrIsA(err));*/
  if (err->error) return;
  /*g_assert (ObitIsA(in, &myClassInfo));*/
  /*g_assert (ObitImageDescIsA(in->myDesc));*/

  /* convert to pixel - assume linear */
  pixel[0] = in->myDesc->crpix[0] + (coord[0] / in->myDesc->cdelt[0]);
  pixel[1] = in->myDesc->crpix[1] + (coord[1] / in->myDesc->cdelt[1]);
 
  /* interpolate */
  ObitCInterpolatePixel (in, pixel, out, err);
} /* end ObitCInterpolatePosition */

/**
 * Interpolate value at requested offset from reference position in array. 
 * The object must have an image descriptor to allow determing
 * pixel coordinates.
 * Interpolation between planes is not supported.
 * \param in    The object to interpolate
 * \param off   Coordinate offset in plane
 * \param out   Complex interpolated value as (real,Imag), magic value blanked
 * \param err   Error stack is pixel not inside image.
 */
void ObitCInterpolateOffset (ObitCInterpolate *in, ofloat *off, ofloat out[2], 
			     ObitErr *err)
{
  ofloat pixel[IM_MAXDIM], fblank = ObitMagicF();

 /* error checks */
  out[0] = fblank; out[1] = fblank;  /* init result to blanked */
  /*g_assert(ObitErrIsA(err)); */
  if (err->error) return;
  /*g_assert (ObitIsA(in, &myClassInfo)); */
  /*g_assert (ObitImageDescIsA(in->myDesc)); */

   /* convert to pixel - assume linear */
  pixel[0] = in->myDesc->crpix[0] + (off[0] / in->myDesc->cdelt[0]);
  pixel[1] = in->myDesc->crpix[1] + (off[1] / in->myDesc->cdelt[1]);

  /* debug
  fprintf(stderr,"pixel %f %f\n", pixel[0], pixel[1]); */

  /* interpolate */
  ObitCInterpolatePixel (in, pixel, out, err);
} /* end ObitCInterpolateOffset */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitCInterpolateClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitCInterpolateClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitCInterpolateClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitCInterpolateClassInfoDefFn (gpointer inClass)
{
  ObitCInterpolateClassInfo *theClass = (ObitCInterpolateClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)theClass->ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitCInterpolateClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitCInterpolateClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitCInterpolateGetClass;
  theClass->newObit       = (newObitFP)newObitCInterpolate;
  theClass->ObitCopy      = (ObitCopyFP)ObitCInterpolateCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitCInterpolateClone;
  theClass->ObitClear     = (ObitClearFP)ObitCInterpolateClear;
  theClass->ObitInit      = (ObitInitFP)ObitCInterpolateInit;

} /* end ObitCInterpolateClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param in Pointer to the object to initialize.
 */
void ObitCInterpolateInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitCInterpolate *in = inn;

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
  in->myKernal    = NULL;
  in->array       = NULL;
  in->xPixel     = -1.0e20;
  in->yPixel     = -1.0e20;

} /* end ObitCInterpolateInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  in Pointer to the object to deallocate.
 *           Actually it should be an ObitCInterpolate* cast to an Obit*.
 */
void ObitCInterpolateClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitCInterpolate *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread      = ObitThreadUnref(in->thread);
  in->info        = ObitInfoListUnref(in->info);
  in->myArray     = ObitCArrayUnref(in->myArray);
  in->myKernal    = ObitFArrayUnref(in->myKernal);
  in->myDesc      = ObitImageDescUnref(in->myDesc);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitCInterpolateClear */



/**
 * Set Lagrangian interpolation kernal taking into account ends of the grid.,
 * \param  in     Object
 * \param  Pixel  Which is the desired pixel (1-rel)
 * \param  naxis  Number of pixels in row
 * \param  Start  [out] first pixels in array for convolution (0-rel)
 * \param  Nterm  [out] number of terms to include
 * \param  Kernal [out] pointer to convolving kernal.
 */
static void SetConvKernal (ObitCInterpolate* in, ofloat Pixel, 
			   olong naxis, olong *Start, olong *Nterm, 
			   ofloat **Kernal)
{
  ofloat off;
  olong ipos, cen, iwid, ioff, itab, ntab;

  /* fractional pixel */
  ipos = Pixel + 0.5;
  iwid = in->hwidth*2 + 1;
  *Nterm = iwid; /* value except at edge */

  /* set first pixel */
  cen = ipos - in->hwidth;
  cen = MAX (1, MIN (cen, (naxis-iwid)));
  *Start = cen-1; /* returned version*/

  /* Find offset in kernal table */
  ntab = in->numKTab - 1;
  off = Pixel - cen - in->hwidth;
  if (off>0.0) ioff = (olong)(off + 0.5);
  else ioff = (olong)(off - 0.5);
  off -= (ofloat)ioff;  /* fraction of cell [-0.5, 0.5]*/
   /* Which table */
  itab = (olong)(((0.5+off)*ntab)+0.5); 
  itab = MAX (0, MIN (itab, ntab));
  /* Set pointer into kernal table */
  *Kernal = in->myKernal->array + itab*in->myKernal->naxis[0];

  /* Corrections near edge */
  if (ioff<=-1) { /* near start */
    *Kernal -= ioff;
    *Nterm  += ioff;
  } else if (ioff>=1){ /* near end */
    *Kernal += ioff;
    *Nterm  -= ioff;
  }
} /* end SetConvKernal */

/**
 * Initialize interpolation:
 * \li Initialize Lagrangian denominators
 * \li Create myKernal structure for kernal table
 *     Table entries twice the kernal width to deal with edges.
 * \li Fill entries in  myKernal
 *     Each entry is for a given fractional pixel (-0.5 -> 0.5)
 * \param in  The object to initialize
 */
static void InitConvKernal (ObitCInterpolate* in)
{
  olong iwid, i, j, k, num;
  ofloat prod, xx, *table;
  olong ndim, naxis[2], Start, Nterm;

  /* Init Lagrangian denominators for hwidth */
  iwid = 1 + (2*in->hwidth);
  for (j= 1; j<=iwid; j++) {
    prod = 1.0;
    for (i= 1; i<=iwid; i++) {
      if (i != j) prod = prod * (j - i);
    } 
    in->denom[j-1] = 1.0 / prod;
  } 
  /* Compute convolution kernal table */
  in->kernalSpace = 0.005;  /* 1/200 */
  num = 1 + (1.0 / in->kernalSpace) + 0.5; /* number of tables */
  in->numKTab = num;

  /* Create/reallocate myKernal table */
  ndim = 2;
  naxis[0] = iwid; naxis[1] = num;
  if (in->myKernal!=NULL) 
    in->myKernal = ObitFArrayRealloc(in->myKernal, ndim, naxis);
  else in->myKernal = ObitFArrayCreate("KernalTable", ndim, naxis);
  naxis[0] = 0; naxis[1] = 0;
  table = ObitFArrayIndex(in->myKernal, naxis); /* Pointer to table */

  /* Loop over "pixel" */
  for (k=0; k<num; k++) {
    xx = in->hwidth + in->kernalSpace*(k-1) + 0.5;
    /* compute interpolating kernal for this "pixel" */
    CalcConvKernal (in, xx, 1024, &Start, &Nterm, table);
    table += iwid;  /* Next row in table */
  } /* end loop over pixels */
} /* end InitConvKernal */

/**
 * Calculate Lagrangian interpolation kernal.
 * \param  in     Object
 * \param  Pixel  Which is the desired pixel (1-rel)
 * \param  naxis  Number of pixels in row
 * \param  Start  [out] first pixels in array for convolution (0-rel)
 * \param  Nterm  [out] number of terms to include
 * \param  Kernal [out] pointer to convolving kernal.
 */
static void CalcConvKernal (ObitCInterpolate* in, ofloat Pixel, 
			   olong naxis, olong *Start, olong *Nterm, 
			   ofloat *Kernal)
{
  ofloat xx, prod;
  olong ipos, cen, iwid, i, j;

  /* fractional pixel */
  ipos = Pixel + 0.5;
  iwid = in->hwidth*2 + 1;
  *Nterm = iwid; /* value except at edge */

  /* set first pixel */
  cen = ipos - in->hwidth;
  cen = MAX (1, MIN (cen, (naxis-iwid)));
  cen--;  /*zero rel */
  *Start = cen; /* returned version*/
 
  /* set "x" at first pixel to 1.0 */
  xx = Pixel - cen + 1.0;

  /* compute interpolating kernal */
  for (j= 0; j<iwid; j++) { /* loop 50 */
    prod = in->denom[j];
    for (i= 0; i<iwid; i++) { /* loop 30 */
      if (i != j) prod = prod * (xx - (i+1));
    } /* end loop  L30:  */;
    Kernal[j] = prod;
  } /* end loop  L50:  */;
} /* end CalcConvKernal */

