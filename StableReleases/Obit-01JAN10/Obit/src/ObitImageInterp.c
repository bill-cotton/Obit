/* $Id: $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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

#include "ObitImageDesc.h"
#include "ObitImageInterp.h"
#include "ObitFInterpolate.h"
#include "ObitThread.h"
#include "ObitImageUtil.h"
#include "ObitPosition.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageInterp.c
 * ObitImageInterp class function definitions.
 * This class is derived from the Obit base class.
 *
 * ObitImageInterp Class to nterpolate  pixel values in images
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImageInterp";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitImageInterpClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageInterpClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageInterpInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageInterpClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageInterpClassInfoDefFn (gpointer inClass);

/** Private: Read Image. */
static void  ReadImage  (ObitImageInterp *in, ObitImage *image, 
		      olong iFreq, olong iIF, olong iplane, 
		      ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitImageInterp* newObitImageInterp (gchar* name)
{
  ObitImageInterp* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageInterpClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitImageInterp));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageInterpInit((gpointer)out);

 return out;
} /* end newObitImageInterp */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageInterpGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageInterpClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitImageInterpGetClass */

/**
 * Make a deep copy of an ObitImageInterp.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitImageInterp* ObitImageInterpCopy  (ObitImageInterp *in, ObitImageInterp *out, 
				     ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i, nfreq;
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
    out = newObitImageInterp(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->ImgDesc   = ObitImageDescUnref(out->ImgDesc);
  out->ImgDesc   = ObitImageDescCopy (in->ImgDesc, out->ImgDesc, err);
  out->ImgPixels = ObitFArrayUnref(out->ImgPixels);
  out->ImgPixels = ObitFArrayCopy (in->ImgPixels, out->ImgPixels, err);
  out->myInterp  = ObitFArrayUnref(out->myInterp);
  out->myInterp  = ObitFInterpolateCopy (in->myInterp, out->myInterp, err);
  out->nplanes   = in->nplanes;
  if (out->freqs) g_free(out->freqs);
  nfreq = out->ImgPixels->naxis[2];
  out->freqs = g_malloc0(nfreq*sizeof(odouble));
  for (i=0; i<nfreq; i++) out->freqs[i] = in->freqs[i];

  return out;
} /* end ObitImageInterpCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an ImageInterp similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitImageInterpClone  (ObitImageInterp *in, ObitImageInterp *out, ObitErr *err)
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
  out->ImgDesc   = ObitImageDescUnref(out->ImgDesc);
  out->ImgDesc   = ObitImageDescCopy (in->ImgDesc, out->ImgDesc, err);
  out->ImgPixels = ObitFArrayUnref(out->ImgPixels);
  ObitFArrayClone (in->ImgPixels, out->ImgPixels, err);
  out->myInterp  = ObitFArrayUnref(out->myInterp);
  out->myInterp  = ObitFInterpolateClone (in->myInterp, out->myInterp);
  out->freqs     = in->freqs;
  out->nplanes   = in->nplanes;
} /* end ObitImageInterpClone */

/**
 * Creates an ObitImageInterp, the order of planes in the output object 
 * is channels vary fastest, then IFs.
 * \param name    An optional name for the object.
 * \param image   Image from which to derive object 
 * \param hwidth  half width in pixels of interpolation, 1 or 2 usually OK
 * \return the new object.
 */
ObitImageInterp* ObitImageInterpCreate (gchar* name, ObitImage *image, 
					olong hwidth, ObitErr *err)
{
  ObitImageInterp* out=NULL;
  olong nImgFreq, nImgIF, nChIF, naxis[5], nx, ny;
  olong iplane, iFreq, iIF;
  gchar *tname;
  gchar *routine = "ObitImageInterpCreate";

  /* Error tests */
  g_assert(ObitImageIsA(image));
  if (err->error) return out;  /* Previous error */

  /* Create basic structure */
  out = newObitImageInterp (name);

  /* Ensure image fully instantiated and OK */
  ObitImageOpen(image, OBIT_IO_ReadOnly, err);
  ObitImageClose(image, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, out);

  /* How big?  */
  nImgFreq = image->myDesc->inaxes[image->myDesc->jlocf];
  if (image->myDesc->jlocif>=0) 
    nImgIF  = image->myDesc->inaxes[image->myDesc->jlocif];
  else nImgIF  = 1;

  nChIF = nImgFreq * nImgIF;
  nx = image->myDesc->inaxes[0];
  ny = image->myDesc->inaxes[1];
  naxis[0] = nx; naxis[1] = ny; naxis[2] = nChIF;
  out->nplanes = nChIF;    /* Number of planes in output */

  /* Save image header */
  out->ImgDesc = ObitImageDescCopy (image->myDesc, out->ImgDesc, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, out);

  /* Create ImgPixels array - big enought for all planes */
  if (name) tname = g_strconcat (name, ":ImgPixels", NULL);
  else tname = g_strdup("ImgPixels") ;
  out->ImgPixels = ObitFArrayCreate (tname, 3, naxis);
  g_free(tname);

  /* Frequency array - filled in in ReadImage */
  out->freqs = g_malloc0(out->nplanes*sizeof(odouble));

  /* Create interpolator */
  out->myInterp = newObitFInterpolateCreate ("Interpolator", 
					     out->ImgPixels, out->ImgDesc, 
					     hwidth);

  /* Loop over planes in order they appear in the UV data */
  iplane = 0;
  for (iIF=0; iIF<nImgIF; iIF++) {
    for (iFreq=0; iFreq<nImgFreq; iFreq++) {
      ReadImage (out, image, iFreq, iIF, iplane, err);
      iplane++;
    } /* end loop over channel */
  } /* end loop over IF */

 return out;
} /* end ObitImageInterpCreate */

/**
 * Interpolate requested beam value
 * Locks in->thread for multithreading
 * \param in      Object to interpolate
 * \param RA      Right Ascension desired (deg) @ std. equinox
 * \param Dec     Declination deg) desired @ std. equinox
 * \param Angle   (Parallactic) Angle to rotate image  (deg)
 * \param plane   Image plane 0-rel (frequency) to use
 *                Can be obtained from ObitImageInterpFindPlane
 * \param err     Obit error stack object.
 * \return interpolated beam value -  may be fblank
 */
ofloat ObitImageInterpValue (ObitImageInterp* in, 
			     odouble RA, odouble Dec, 
			     ofloat Angle, olong plane,
			     ObitErr *err)
{
  ofloat beamValue=0.0;
  ofloat pixel[3], fblank = ObitMagicF();
  odouble coord[2];
  gchar *routine = "ObitImageInterpValue";

  if (err->error) return beamValue;  /* Previous error? */


  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* Modify header for parallactic angle */
  in->ImgDesc->crota[1] = Angle; 

  /* Convert position to pixels */
  coord[0] = RA; coord[1] = Dec;
  if (ObitPositionXYpix(coord, in->ImgDesc, pixel)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: error converting coordinates", routine);
    beamValue = fblank;
    goto cleanup;
  }
  pixel[2] = MIN (plane+1, in->nplanes); 

  /* Interpolate */
  beamValue = ObitFInterpolatePixel (in->myInterp, pixel, err);

  cleanup:
  /* Unlock ObitObjects against other threads */
  ObitThreadUnlock(in->thread);

  return beamValue;
} /* end ObitImageInterpValue */

/**
 * Interpolate requested beam value given the interpolator
 * This allows non blocking interpolation in a multithreaded environment.
 * Locks in->thread for multithreading
 * Use this routine for best multi-threaded performance and in most useful
 * for interpolating beam images (x=azimuth, y=elevation)
 * \param in      Object to interpolate
 * \param interp  Interpolator
 * \param RA      X or Right Ascension desired (deg) @ std. equinox
 * \param Dec     Y or Declination deg) desired @ std. equinox
 * \param Angle   (Parallactic) Angle to rotate  (deg)
 * \param plane   Image plane 0-rel (frequency) to use
 *                Can be obtained from ObitImageInterpFindPlane
 * \param err     Obit error stack object.
 * \return interpolated beam value -  may be fblank
 */
ofloat ObitImageInterpValueInt (ObitImageInterp* in,  ObitFInterpolate* interp,
				odouble RA, odouble Dec, 
				ofloat Angle, olong plane,
				ObitErr *err)
{
  ofloat beamValue=0.0;
  ofloat pixel[3], fblank = ObitMagicF();
  odouble coord[2];
  gchar *routine = "ObitImageInterpValue";

  if (err->error) return beamValue;  /* Previous error? */
  
  
  /* Modify header for parallactic angle */
  interp->myDesc->crota[1] = Angle; 

  /* Convert position to pixels */
  coord[0] = RA; coord[1] = Dec;
  if (ObitPositionXYpix(coord, interp->myDesc, pixel)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: error converting coordinates", routine);
    beamValue = fblank;
    goto cleanup;
  }
  pixel[2] = MIN (plane+1, in->nplanes); 

  /* Interpolate */
  beamValue = ObitFInterpolatePixel (interp, pixel, err);

  cleanup:

  return beamValue;
} /* end ObitImageInterpValueInt */

/**
 * Find closest beam image plane to a given frequency
 * \param in      Object to interpolate
 * \param freq    Frequency (Ha) to lookup
 * \return closest 0-rel plane number
 */
olong ObitImageInterpFindPlane (ObitImageInterp* in, odouble freq)
{
  olong close = 0;
  olong i;
  ofloat diff, d;

  /* Initialize */
  diff = fabs (freq-in->freqs[close]);
  for (i=1; i<in->nplanes; i++) {
    d = fabs (freq-in->freqs[i]);
    if (d<diff) { close=i; diff=d;}
    if (d>diff) break;  /* Getting further away? */
  }

  return close;
} /* end ObitImageInterpFindPlane */

/**
 * Returns clone of interpolator allowing nonblocking interpolation
 * Returned interpolator has independent internals and an independent copy
 * of the image description.
 * \param in      Object to whose interpolator to clone
 * \return interpolator, Unref wnen done
 */
ObitFInterpolate*  ObitImageInterpCloneInterp (ObitImageInterp* in, ObitErr *err)
{
  ObitFInterpolate* out=NULL;

  out = ObitFInterpolateClone (in->myInterp, out);

  /* Full copy of descriptor */
  out->myDesc = ObitImageDescUnref(out->myDesc);
  out->myDesc = ObitImageDescCopy(in->myInterp->myDesc, out->myDesc, err);

  return out;
} /* end ObitImageInterpCloneInterp */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageInterpClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitImageInterpClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageInterpClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageInterpClassInfoDefFn (gpointer inClass)
{
  ObitImageInterpClassInfo *theClass = (ObitImageInterpClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageInterpClassInit;
  theClass->newObit       = (newObitFP)newObitImageInterp;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageInterpClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageInterpGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageInterpCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageInterpClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageInterpInit;
  theClass->ObitImageInterpCreate = (ObitImageInterpCreateFP)ObitImageInterpCreate;
  theClass->ObitImageInterpValue  = (ObitImageInterpValueFP)ObitImageInterpValue;
  theClass->ObitImageInterpValueInt  = 
    (ObitImageInterpValueIntFP)ObitImageInterpValueInt;
  theClass->ObitImageInterpFindPlane = 
    (ObitImageInterpFindPlaneFP)ObitImageInterpFindPlane;
  theClass->ObitImageInterpCloneInterp = 
    (ObitImageInterpCloneInterpFP)ObitImageInterpCloneInterp;
} /* end ObitImageInterpClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageInterpInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageInterp *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread     = newObitThread();
  in->info       = newObitInfoList(); 
  in->ImgDesc   = NULL;
  in->ImgPixels = NULL;
  in->myInterp   = NULL;
  in->freqs      = NULL;
  in->nplanes    = 0;

} /* end ObitImageInterpInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImageInterp* cast to an Obit*.
 */
void ObitImageInterpClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageInterp *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->info       = ObitInfoListUnref(in->info);
  in->thread     = ObitThreadUnref(in->thread);
  in->ImgDesc   = ObitImageDescUnref(in->ImgDesc);
  in->ImgPixels = ObitFArrayUnref(in->ImgPixels);
  in->myInterp   = ObitFInterpolateUnref(in->myInterp);
  if (in->freqs) g_free(in->freqs);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageInterpClear */

/**
 * Read Image
 * Reads  image and leaves in in->ImgPixels
 * \param in     Object to load image into
 * \param image  Image from which to derive object 
 * \param iFreq  0-rel freq number in image
 * \param iIF    0-rel IF number in image
 * \param iplane 0-rel plane number in in->ImgPixels
 * \param err Obit error stack object.
 */
static void  ReadImage  (ObitImageInterp *in, ObitImage *image, 
			olong iFreq, olong iIF, olong iplane, 
			ObitErr *err) 
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, ipos[3], jlocf, jlocif;
  olong blc[7] = {1,1,1,1,1,1,1}, trc[7] = {0,0,0,0,0,0,0};
  ofloat *inArr, *outArr; 
  gchar *routine = "ObitImageInterp:ReadImage";

  if (err->error) return;  /* Previous error? */

  /* Select image plane */
  if (image->myDesc->jlocf>=0)  {
    jlocf  = image->myDesc->jlocf;
    blc[jlocf] = trc[jlocf] = iFreq+1;
  }
  if (image->myDesc->jlocif>=0) {
    jlocif = image->myDesc->jlocif;
    blc[jlocif] = trc[jlocif] = iIF+1;
  }
  dim[0] = 7;
  ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc);

  /* Read input */
  ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  ObitImageRead (image, NULL, err);
  ObitImageClose (image,err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Set plane frequency */
  in->freqs[iplane] = image->myDesc->crval[image->myDesc->jlocf];

  /* Copy to in->ImgPixels */
  ipos[0] = ipos[1] = ipos[2] = 0;
  inArr = ObitFArrayIndex (image->image, ipos);
  ipos[0] = ipos[1] = 0; ipos[2] = iplane;
  outArr = ObitFArrayIndex (in->ImgPixels, ipos);
  for (i=0; i<image->image->arraySize; i++) *outArr++ = *inArr++;
} /* end ReadImage */
