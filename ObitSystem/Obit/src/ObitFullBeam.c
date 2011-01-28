/* $Id$        */
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
#include "ObitFullBeam.h"
#include "ObitFInterpolate.h"
#include "ObitThread.h"
#include "ObitImageUtil.h"
#include "ObitPosition.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFullBeam.c
 * ObitFullBeam class function definitions.
 * This class is derived from the Obit base class.
 *
 * ObitFullBeam Class to generate full beam corrections from beam images
 * An ObitFullBeam takes images or (hyper)cubes  of a primary beam in a given
 * Stokes parameter and assists in image plane corrections
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitFullBeam";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitFullBeamClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFullBeamClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFullBeamInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFullBeamClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFullBeamClassInfoDefFn (gpointer inClass);

/** Private: Read Beam. */
static void  ReadBeam  (ObitFullBeam *in, ObitImage *image, 
		      olong iFreq, olong iIF, olong iplane, 
		      ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitFullBeam* newObitFullBeam (gchar* name)
{
  ObitFullBeam* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFullBeamClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFullBeam));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFullBeamInit((gpointer)out);

 return out;
} /* end newObitFullBeam */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFullBeamGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFullBeamClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitFullBeamGetClass */

/**
 * Make a deep copy of an ObitFullBeam.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitFullBeam* ObitFullBeamCopy  (ObitFullBeam *in, ObitFullBeam *out, 
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
    out = newObitFullBeam(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->BeamDesc   = ObitImageDescUnref(out->BeamDesc);
  out->BeamDesc   = ObitImageDescCopy (in->BeamDesc, out->BeamDesc, err);
  out->BeamPixels = ObitFArrayUnref(out->BeamPixels);
  out->BeamPixels = ObitFArrayCopy (in->BeamPixels, out->BeamPixels, err);
  out->myInterp   = ObitFArrayUnref(out->myInterp);
  out->myInterp   = ObitFInterpolateCopy (in->myInterp, out->myInterp, err);
  out->nplanes    = in->nplanes;
  if (out->freqs) g_free(out->freqs);
  nfreq = out->BeamPixels->naxis[2];
  out->freqs = g_malloc0(nfreq*sizeof(odouble));
  for (i=0; i<nfreq; i++) out->freqs[i] = in->freqs[i];

  return out;
} /* end ObitFullBeamCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an FullBeam similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitFullBeamClone  (ObitFullBeam *in, ObitFullBeam *out, ObitErr *err)
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
  out->BeamDesc   = ObitImageDescUnref(out->BeamDesc);
  out->BeamDesc   = ObitImageDescCopy (in->BeamDesc, out->BeamDesc, err);
  out->BeamPixels = ObitFArrayUnref(out->BeamPixels);
  ObitFArrayClone (in->BeamPixels, out->BeamPixels, err);
  out->myInterp   = ObitFArrayUnref(out->myInterp);
  out->myInterp   = ObitFInterpolateClone (in->myInterp, out->myInterp);
  out->freqs      = in->freqs;
  out->nplanes    = in->nplanes;
} /* end ObitFullBeamClone */

/**
 * Creates an ObitFullBeam, the order of planes in the output object 
 * is channels vary fastest, then IFs.
 * \param name    An optional name for the object.
 * \param myInput InfoList with optional control parameters (may be NULL):
 * \li "antSize" OBIT_float scalar Diameter of antennas for gain. [def 24.5]
 * \li "minGain" OBIT_float scalar  Min. allowed antenna gain, [def 0.02]
 * \param image   Image from which to derive object 
 * \return the new object.
 */
ObitFullBeam* ObitFullBeamCreate (gchar* name, ObitInfoList *myInput,
				  ObitImage *image, ObitErr *err)
{
  ObitFullBeam* out=NULL;
  olong nImgFreq, nImgIF, nChIF, naxis[5], nx, ny;
  olong iplane, iFreq, iIF;
  gchar *dataParms[] = {"antSize","minGain",NULL};
  gchar *tname;
  gchar *routine = "ObitFullBeamCreate";

  /* Error tests */
  g_assert(ObitImageIsA(image));
  if (err->error) return out;  /* Previous error */

  /* Create basic structure */
  out = newObitFullBeam (name);

  /* Ensure image fully instantiated and OK */
  ObitImageOpen(image, OBIT_IO_ReadOnly, err);
  ObitImageClose(image, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, out);

  /* Save any control parameters */
  if (myInput)
    ObitInfoListCopyList (myInput, out->info, dataParms);

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

  /* Save image Beam header */
  out->BeamDesc = ObitImageDescCopy (image->myDesc, out->BeamDesc, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, out);

  /* Create BeamPixels array - big enought for all planes */
  if (name) tname = g_strconcat (name, ":BeamPixels", NULL);
  else tname = g_strdup("BeamPixels") ;
  out->BeamPixels = ObitFArrayCreate (tname, 3, naxis);
  g_free(tname);

  /* Frequency array - filled in in ReadBeam */
  out->freqs = g_malloc0(out->nplanes*sizeof(odouble));

  /* Create interpolator */
  out->myInterp = newObitFInterpolateCreate ("Interpolator", 
					     out->BeamPixels, out->BeamDesc, 2);

  /* Loop over planes in order they appear in the UV data */
  iplane = 0;
  for (iIF=0; iIF<nImgIF; iIF++) {
    for (iFreq=0; iFreq<nImgFreq; iFreq++) {
      ReadBeam (out, image, iFreq, iIF, iplane, err);
      iplane++;
    } /* end loop over channel */
  } /* end loop over IF */

 return out;
} /* end ObitFullBeamCreate */

/**
 * Interpolate requested beam value
 * Locks in->thread for multithreading
 * \param in      Object to interpolate
 * \param dRA     Right Ascension offset desired (deg) @ std. equinox
 * \param dDec    Declination offset (deg) desired @ std. equinox
 * \param PAngle  Parallactic Angle  (deg)
 * \param plane   Image plane 0-rel (frequency) to use
 *                Can be obtained from ObitFullBeamFindPlane
 * \param err     Obit error stack object.
 * \return interpolated beam value -  may be fblank
 */
ofloat ObitFullBeamValue (ObitFullBeam* in, 
			  odouble dRA, odouble dDec, 
			  ofloat PAngle, olong plane,
			  ObitErr *err)
{
  ofloat beamValue=0.0;
  ofloat pixel[3], fblank = ObitMagicF();
  odouble coord[2];
  gchar *routine = "ObitFullBeamValue";

  if (err->error) return beamValue;  /* Previous error? */


  /* Lock ObitObjects against other threads */
  ObitThreadLock(in->thread);

  /* Modify header for parallactic angle */
  in->BeamDesc->crota[1] = PAngle; 

  /* Convert position to pixels */
  coord[0] = dRA; coord[1] = dDec;
  if (ObitPositionXYpix(coord, in->BeamDesc, pixel)) {
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
} /* end ObitFullBeamValue */

/**
 * Interpolate requested beam value given the interpolator
 * This allows non blocking interpolation in a multithreaded environment.
 * Locks in->thread for multithreading
 * Use this routine for best multi-threaded performance.
 * \param in      Object to interpolate
 * \param interp  Interpolator
 * \param dRA     Right Ascension offset desired (deg) @ std. equinox
 * \param dDec    Declination offset (deg) desired @ std. equinox
 * \param PAngle  Parallactic Angle  (deg)
 * \param plane   Image plane 0-rel (frequency) to use
 *                Can be obtained from ObitFullBeamFindPlane
 * \param err     Obit error stack object.
 * \return interpolated beam value -  may be fblank
 */
ofloat ObitFullBeamValueInt (ObitFullBeam* in,  ObitFInterpolate* interp,
			     odouble dRA, odouble dDec, 
			     ofloat PAngle, olong plane,
			     ObitErr *err)
{
  ofloat beamValue=0.0;
  ofloat pixel[3], fblank = ObitMagicF();
  odouble coord[2];
  gchar *routine = "ObitFullBeamValue";

  if (err->error) return beamValue;  /* Previous error? */


  /* Modify header for parallactic angle */
  interp->myDesc->crota[1] = PAngle; 


  /* Convert position to pixels */
  coord[0] = dRA; coord[1] = dDec;
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
} /* end ObitFullBeamValueInt */

/**
 * Find closest beam image plane to a given frequency
 * \param in      Object to interpolate
 * \param freq    Frequency (Ha) to lookup
 * \return closest 0-rel plane number
 */
olong ObitFullBeamFindPlane (ObitFullBeam* in, odouble freq)
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
} /* end ObitFullBeamFindPlane */

/**
 * Returns clone of interpolator allowing nonblocking interpolation
 * Returned interpolator has independent internals and an independent copy
 * of the image description.
 * \param in      Object to whose interpolator to clone
 * \return interpolator, Unref wnen done
 */
ObitFInterpolate*  ObitFullBeamCloneInterp (ObitFullBeam* in, ObitErr *err)
{
  ObitFInterpolate* out=NULL;

  out = ObitFInterpolateClone (in->myInterp, out);

  /* Full copy of descriptor */
  out->myDesc = ObitImageDescUnref(out->myDesc);
  out->myDesc = ObitImageDescCopy(in->myInterp->myDesc, out->myDesc, err);

  return out;
} /* end ObitFullBeamCloneInterp */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFullBeamClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFullBeamClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFullBeamClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFullBeamClassInfoDefFn (gpointer inClass)
{
  ObitFullBeamClassInfo *theClass = (ObitFullBeamClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFullBeamClassInit;
  theClass->newObit       = (newObitFP)newObitFullBeam;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFullBeamClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFullBeamGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitFullBeamCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFullBeamClear;
  theClass->ObitInit      = (ObitInitFP)ObitFullBeamInit;
  theClass->ObitFullBeamCreate = (ObitFullBeamCreateFP)ObitFullBeamCreate;
  theClass->ObitFullBeamValue  = (ObitFullBeamValueFP)ObitFullBeamValue;
  theClass->ObitFullBeamValueInt  = 
    (ObitFullBeamValueIntFP)ObitFullBeamValueInt;
  theClass->ObitFullBeamFindPlane = 
    (ObitFullBeamFindPlaneFP)ObitFullBeamFindPlane;
  theClass->ObitFullBeamCloneInterp = 
    (ObitFullBeamCloneInterpFP)ObitFullBeamCloneInterp;
} /* end ObitFullBeamClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitFullBeamInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFullBeam *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread     = newObitThread();
  in->info       = newObitInfoList(); 
  in->BeamDesc   = NULL;
  in->BeamPixels = NULL;
  in->myInterp   = NULL;
  in->freqs      = NULL;
  in->nplanes    = 0;

} /* end ObitFullBeamInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitFullBeam* cast to an Obit*.
 */
void ObitFullBeamClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFullBeam *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->info       = ObitInfoListUnref(in->info);
  in->thread     = ObitThreadUnref(in->thread);
  in->BeamDesc   = ObitImageDescUnref(in->BeamDesc);
  in->BeamPixels = ObitFArrayUnref(in->BeamPixels);
  in->myInterp   = ObitFInterpolateUnref(in->myInterp);
  if (in->freqs) g_free(in->freqs);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitFullBeamClear */

/**
 * Read and pad Beam
 * Reads  image and leaves in in->BeamPixels
 * For Stokes I the image is normalized by a symmetric function.
 * \param in     Object to rotate, info member may contain:
 * \li "antSize" OBIT_float scalar Diameter of antennas for gain. [def 24.5]
 * \li "minGain" OBIT_float scalar  Min. allowed antenna gain, [def 0.02]
 * \param image  Image from which to derive object 
 * \param iFreq  0-rel freq number in image
 * \param iIF    0-rel IF number in image
 * \param iplane 0-rel plane number in in->BeamPixels
 * \param err Obit error stack object.
 */
static void  ReadBeam  (ObitFullBeam *in, ObitImage *image, 
			olong iFreq, olong iIF, olong iplane, 
			ObitErr *err) 
{
  ObitImage *scrImage=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, ipos[3], jlocf, jlocif;
  olong blc[7] = {1,1,1,1,1,1,1}, trc[7] = {0,0,0,0,0,0,0};
  olong outPlane[5] = {1,1,1,1,1};
  ofloat *inArr, *outArr, antSize=24.5, minGain = 0.02; 
  gchar *routine = "ObitFullBeam:ReadBeam";

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

  /* If Stokes I Divide by symmetric beam */
  /* DISABLED (crval==10.0 test, I=1.0) */
  if (image->myDesc->crval[image->myDesc->jlocs]==10.0) {

    /* Scratch Memory only image for symmetric beam */
    scrImage = newObitImage ("Beam");
    ObitImageCloneMem (image, scrImage, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    /* One plane */
    scrImage->myDesc->inaxes[2] = scrImage->myDesc->inaxes[3] = 1;
    
    /* Control parameters */
    ObitInfoListGetTest(in->info, "antSize", &type, dim, &antSize);
    ObitInfoListGetTest(in->info, "mingain", &type, dim, &minGain);
   
    /* Create Beam image */
    ObitImageUtilPBImage (image, scrImage, outPlane, antSize, minGain, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* Normalize */
    ObitFArrayDiv (image->image, scrImage->image, image->image);

    /* Cleanup */
    scrImage = ObitImageUnref(scrImage);
  } /* End normalize by Stokes I symmetric beam */

  /* Copy to in->BeamPixels */
  ipos[0] = ipos[1] = ipos[2] = 0;
  inArr = ObitFArrayIndex (image->image, ipos);
  ipos[0] = ipos[1] = 0; ipos[2] = iplane;
  outArr = ObitFArrayIndex (in->BeamPixels, ipos);
  for (i=0; i<image->image->arraySize; i++) *outArr++ = *inArr++;
} /* end ReadBeam */
