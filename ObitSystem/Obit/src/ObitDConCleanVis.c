/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2009                                          */
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

#include "ObitDConCleanVis.h"
#include "ObitMem.h"
#include "ObitFFT.h"
#include "ObitTableUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitSkyGeom.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanVis.c
 * ObitDConCleanVis class function definitions.
 * Image based CLEAN class.
 * This class is derived from the ObitDCon class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanVis";

/** Function to obtain parent ClassInfo - ObitDConClean */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanGetClass;

/**
 * ClassInfo structure ObitDConCleanVisClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanVisClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanVisInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanVisClear (gpointer in);

/** Private: (re)make residuals. */
static void  MakeResidual (ObitDConCleanVis *in, olong field, 
			   gboolean doBeam, ObitErr *err);

/** Private: (re)make all residuals. */
static void  MakeAllResiduals (ObitDConCleanVis *in, ObitErr *err);

/** Private: Find best residual image. */
static void  WhosBest (ObitDConCleanVis *in, olong *best, olong *second);

/** Private: Set Class function pointers. */
static void ObitDConCleanVisClassInfoDefFn (gpointer inClass);

/** Private: reset sky model. */
static gboolean ResetSkyModel (ObitDConCleanVis *in, ObitErr *err);

/** Private: Get pixel array for a given field. */
static ObitFArray* GetFieldPixArray (ObitDConCleanVis *in, olong field, 
				     ObitErr *err);

/** Private: Low accuracy subtract pixels for a given field. */
static void SubNewCCs (ObitDConCleanVis *in, olong field, olong newCC, 
		       ObitFArray *pixarray, ObitErr *err);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanVis* newObitDConCleanVis (gchar* name)
{
  ObitDConCleanVis* out;
  /*gchar *routine = "newObitDConCleanVis";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanVis));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanVisInit((gpointer)out);

 return out;
} /* end newObitDConCleanVis */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanVisGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanVisGetClass */

/**
 * Make a deep copy of an ObitDConCleanVis.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanVis* ObitDConCleanVisCopy  (ObitDConCleanVis *in, 
					     ObitDConCleanVis *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i, nfield;
  gboolean oldExist;
  gchar *outName;
  gchar *routine = "ObitDConCleanVisCopy";

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
    out = newObitDConCleanVis(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /*  copy this class */
  out->imager   = ObitUVImagerUnref(out->imager);
  out->skyModel = ObitSkyModelUnref(out->skyModel);
  out->display  = ObitDisplayUnref(out->display);
  out->imager   = ObitUVImagerRef(in->imager);
  out->skyModel = ObitSkyModelRef(in->skyModel);
  out->display  = ObitDisplayRef(in->display);

  /* Arrays */
  /* out with the old */
  out->quality = ObitMemFree (out->quality);

  /* In with the new */
  nfield       = in->nfield;
  out->quality = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  for (i=0; i<nfield; i++) {
    out->quality[i]    = in->quality[i];
  }

  return out;
} /* end ObitDConCleanVisCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanVis similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanVisClone  (ObitDConCleanVis *in, ObitDConCleanVis *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i, nfield;
  gchar *routine = "ObitDConCleanVisClone";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /*  copy this class */
  out->imager   = ObitUVImagerUnref(out->imager);
  out->skyModel = ObitSkyModelUnref(out->skyModel);
  out->display  = ObitDisplayUnref(out->display);
  out->imager   = ObitUVImagerRef(in->imager);
  out->skyModel = ObitSkyModelRef(in->skyModel);
  out->display  = ObitDisplayRef(in->display);

  /* Arrays */
  /* out with the old */
  out->quality = ObitMemFree (out->quality);

  /* In with the new */
  nfield       = in->nfield;
  out->quality = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  for (i=0; i<nfield; i++) {
    out->quality[i]    = in->quality[i];
  }
} /* end ObitDConCleanVisClone */

/**
 * Creates an ObitDConCleanVis 
 * defined for convenience of derived classes 
 * \param name   An optional name for the object.
 * \param uvdata from which to create object, should have all control
                 information defined on info member.
 * \param err    Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVis* ObitDConCleanVisCreate (gchar* name, ObitUV *uvdata,  
					  ObitErr *err)
{
  olong nfield, i;
  ObitDConCleanVis* out=NULL;
  gchar *routine = "ObitDConCleanVisCreate";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVis (name);

  /* Create UV imager and its ImageMosaic */
  out->imager = ObitUVImagerCreate("UVImager", uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Save uv Mosaic reference */
  out->mosaic = ObitUVImagerGetMosaic(out->imager, err);

  /* Create SkyModel object */
  out->skyModel = ObitSkyModelCreate ("SkyModel", out->mosaic);

   /* Copy control info to SkyModel */
  ObitInfoListCopyData(uvdata->info, out->skyModel->info);
  
  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Arrays per field - including those in parent classes */
  nfield =  out->mosaic->numberImages;
  out->nfield  = nfield;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->quality = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->maxAbsRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i] = -1.0;
    out->avgRes[i]    = -1.0;
    out->quality[i]   = -1.0;
  }

  return out;
} /* end ObitDConCleanVisCreate */

/**
 * Creates an ObitDConCleanVis from optional components 
 * defined for convenience of derived classes 
 * \param name     An optional name for the object.
 * \param uvdata   from which to create object, should have all control
                   information defined on info member.
 * \param imager   Optional ObitUVImager to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
 * \param skyModel Optional ObitSkyModel to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVis* 
ObitDConCleanVisCreate2 (gchar* name, ObitUV *uvdata,  
			 ObitUVImager *imager, ObitSkyModel *skyModel, 
			 ObitErr *err)
{
  olong nfield, i;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitDConCleanVis* out=NULL;
  ofloat ftemp;
  gchar *routine = "ObitDConCleanVisCreate";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVis (name);

  /* Use or create UV imager and create its ImageMosaic */
  if (imager==NULL) {
    out->imager = ObitUVImagerCreate("UVImager", uvdata, err);
    if (err->error) Obit_traceback_val (err, routine, name, out);
  } else out->imager = ObitUVImagerRef(imager);

  /* Save uv Mosaic reference */
  out->mosaic = ObitUVImagerGetMosaic(out->imager, err);

  /*  Use or create SkyModel object */
  if (skyModel==NULL) out->skyModel = ObitSkyModelCreate ("SkyModel", out->mosaic);
  else out->skyModel = ObitSkyModelRef(skyModel);

   /* Copy control info to SkyModel */
  ObitInfoListCopyData(uvdata->info, out->skyModel->info);
  /* Disable any value of minFlux */
  ftemp = -1.0e20;
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (out->skyModel->info, "minFlux", OBIT_float, dim, &ftemp);
  
  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Arrays per field - including those in parent classes */
  nfield =  out->mosaic->numberImages;
  out->nfield  = nfield;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->quality = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->maxAbsRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i] = -1.0;
    out->avgRes[i]    = -1.0;
    out->quality[i]   = -1.0;
  }

  return out;
} /* end ObitDConCleanVisCreate2 */

/**
 * Do deconvolution, uses function on class pointer
 * Does final flatten if FullField member of mosaic member is defined.
 * CLEAN control parameters are in the ObitInfoList member:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 50]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form
 * \li "CCVer"   OBIT_long array    = CLEAN table version for all fields
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimun flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "autoCen" OBIT_float scalar = Threshold leven for autocenter
 *               If peak exceeds this value the minFlux reset to 0.1*autoCen
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanVisDeconvolve (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanVis *in;
  ObitFArray *pixarray=NULL;
  gboolean done, fin, quit, doSub, bail, newWin, moreClean;
  olong jtemp, i, startCC, newCC, count;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  const ObitDConCleanVisClassInfo *inClass;
  const ObitUVImagerClassInfo *imagerClass;
  gchar *routine = "ObitDConCleanVisDeconvolve";

  /* DEBUG
  olong boom[2], xxcnt = 0; */
 
  /* Cast input to this type */
  in = (ObitDConCleanVis*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanVisIsA(in));

  inClass     = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* clean class structure */
  imagerClass = (ObitUVImagerClassInfo*)in->imager->ClassInfo;  /* imager class structure */

  /* Reset highest peak */
  in->peakFlux = -1000.0;
  bail = FALSE;  /* if going to recenter facets, don't need to make final residuals */

  /* Get parameters */
  inClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Recenter any autoCenter windows is needed */
  if (in->doRecenter) ObitDConCleanVisRecenter (in, in->imager->uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Visibility selection and weighting */
  if (in->doWeight) imagerClass->ObitUVImagerWeight (in->imager, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* If doing SDI Clean save copy of uvwork data */
  if (in->SDIGain>0.0) {
    in->SDIdata = newObitUVScratch (in->imager->uvwork, err);
    in->SDIdata = ObitUVCopy (in->imager->uvwork, in->SDIdata, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Create Pixel List if needed - has size changed? */
  if (in->Pixels && (in->Pixels->nfield!=in->mosaic->numberImages)) 
    in->Pixels = ObitDConCleanPxListUnref(in->Pixels);
  if (!in->Pixels) {
    in->Pixels = ObitDConCleanPxListCreate("Pixel List", in->mosaic, 
					   in->maxPixel, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* No more fields than in mosaic */
  Obit_return_if_fail ((in->nfield == in->mosaic->numberImages),
		       err, "%s: CLEAN and mosaic different number of fields %d %d",
		       routine, in->nfield, in->mosaic->numberImages);

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);

  /* Reset Sky model */
  doSub = ResetSkyModel (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  if (doSub) inClass->ObitDConCleanSub((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Reset/Init Pixel list if not done in ResetSkyModel */
  if (!doSub) ObitDConCleanPxListReset (in->Pixels, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Be sure to (re)generate residuals */
  for (i=0; i<in->nfield; i++) in->maxAbsRes[i] = -1.0;

  /* Save actual CC version if not specified */
  if (in->CCver<=0) {
    if (in->Pixels->CCver[0]>0) in->CCver = in->Pixels->CCver[0];
    jtemp = in->CCver;
    ObitInfoListAlwaysPut(in->info, "CCVer", OBIT_long, dim, &jtemp);
  }

  /* Make initial images */
  MakeAllResiduals (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);
  in->doBeam = FALSE;  /* Shouldn't need again */

  /* Loop until Deconvolution done */
  done = FALSE;
  while (!done) {
    /* Decide which field to do next, tells if finished CLEAN */
    done = ObitDConCleanVisPickNext((ObitDConCleanVis*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Does the peak flux density exceed the autocenter threshold */
    if (fabs (in->peakFlux) > in->autoCen) {
      for (i=0; i<in->nfield; i++) {
	in->minFlux[i] = MAX (0.1*in->autoCen, in->minFlux[i]);
	/* Value that counts is on the PxList */
	in->Pixels->minFlux[i] = in->minFlux[i];
      }
      in->autoCen = 1.0e20;  /* Only once */
      bail = TRUE;  /* if going to recenter facets, don't need to make final residuals */
      Obit_log_error(err, OBIT_InfoErr,"Truncating CLEAN at %g Jy for strong source processing",
		     in->minFlux[0]);
    } /* End autoCenter test*/

    /* If no cleaning really requested don't restore */
    bail = bail || (fabs (in->peakFlux) < in->minFlux[0]);

    if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
    if (done) break;

    /* Display/edit windows if enabled */
    if (in->display) {
      quit = ObitDisplayShow (in->display, (Obit*)in->mosaic, in->window, 
			      in->currentField, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      if (quit) {done=TRUE; break;}
    }

    /* If using autoWindow, iterate until clean not limited by autoWindow flux */
    startCC   = in->Pixels->iterField[in->currentField-1];  /* no. at start */
    moreClean = TRUE;
    /* Pixel array for field */
    pixarray  = GetFieldPixArray (in, in->currentField, err);  

    /* Get image/beam statistics needed for this cycle */
    newWin = inClass->ObitDConCleanPixelStats((ObitDConClean*)in, pixarray, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    count = 0;
    while (moreClean) {
      /* Check if new box added to window after first pass */
      if ((count>0) && in->autoWindow) 
      newWin = inClass->ObitDConCleanAutoWindow ((ObitDConClean*)in, in->currentField, 
						 pixarray, err);
      else {
	newWin = in->autoWindow;  /* pretend on the first pass */
      }
      if (err->error) Obit_traceback_msg (err, routine, in->name);
     

      /* Pick components for this major cycle */
      newCC = in->Pixels->iterField[in->currentField-1];  /* no. at start */
      if ((count==0) || newWin) 
	fin   = inClass->ObitDConCleanSelect((ObitDConClean*)in, pixarray, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
      else ObitErrClear(err);

      /* Did it just  stop because of autoWindow? */
      moreClean = (in->Pixels->complCode==OBIT_CompReasonAutoWin);
      /* Done do this forever */
      if ((count>0) && !newWin) moreClean = FALSE;
      if (count>10) moreClean = FALSE;
      if (!moreClean) break;

      /* Subtract these CCs */
      SubNewCCs(in, in->currentField, newCC+1, pixarray, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      count++;
    }
    pixarray = ObitFArrayUnref(pixarray);  /* Release working pixel array */

    /* Update quality list for new max value on field just CLEANed */
    if (fin) in->maxAbsRes[in->currentField-1] = 0.0;
    else in->maxAbsRes[in->currentField-1] = in->Pixels->maxResid;
    in->quality[in->currentField-1] = 
      ObitDConCleanVisQuality((ObitDConCleanVis*)in, in->currentField, err);

    /* Subtract any components from visibility data */
    if (in->Pixels->iterField[in->currentField-1]> startCC)
      inClass->ObitDConCleanSub((ObitDConClean*)in, err);
    else /* Tell if prtLv>2 */
      if (in->prtLv>2) 
	Obit_log_error(err, OBIT_InfoErr, "Field %d has no components to subtract",
		     in->currentField);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

  } /* end clean loop */
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Make final residuals */
  if ((!bail) && (in->niter>0)) MakeAllResiduals (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Final Restore */
  if (in->doRestore && !bail) {

    inClass->ObitDConCleanRestore((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
  }
  
  if (in->doXRestore && !bail) {
    /* Cross Restore if multiple overlapping fields */
    inClass->ObitDConCleanXRestore((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
  }

  /* Flatten if needed */
  if (in->doFlatten && !bail) {
    inClass->ObitDConCleanFlatten((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
    
    /* Display flattened image if enabled */
    if (in->display && in->mosaic->FullField) 
      ObitDisplayShow (in->display, (Obit*)in->mosaic->FullField, NULL, 
		       1, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

} /* end ObitDConCleanVisDeconvolve */

/**
 * Read any base class parameters and then
 * read CLEAN control parameters from the ObitInfoList member:
 * \li "Mode"      OBIT_long scalar = Model mode (ObitSkyModelMode) [def OBIT_SkyModel_Fastest]
 * \li "doRestore" OBIT_bool       = Restore image when done? [def TRUE]
 * \li "doXRestore" OBIT_bool      = Cross restore images when done? [def doRestore]
 * \li "doFlatten" OBIT_bool       = Flatten image when done? [def TRUE]
 * \li "doWeight"  OBIT_bool       = Weight UV data before imaging? [def TRUE]
 * \li "doBeam"    OBIT_bool       = Need to (re)make beam? [def TRUE]
 * \li "doRecenter" OBIT_bool      = Allow recentering autoCenter fields? [def TRUE]
 * \li "reuseFlux" OBIT_float      = Level of Components in initial 
 *                                   SkyModel to reuse, <0 -> none [def none]
 * From Parent classes:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 100]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form  (",", deg)
 * \li "CCVer"   OBIT_long array    = CLEAN table version per field
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimum flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "autoWindow" OBIT_boolean scalar = True if autoWindow feature wanted.
 * \li "autoCen" OBIT_float scalar = Threshold leven for autocenter
 *                 If peak exceeds this value the minFlux reset to 0.1*autoCen
 * \li "dispURL" OBIT_string scalar = URL of display server
 * \param in  The CLEAN object as base class
 * \param err Obit error stack object.
 */
void  ObitDConCleanVisGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanVis *in = (ObitDConCleanVis*)inn;  /* as this class */
  ObitDConClassInfo *ParentClass;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong i;
  union ObitInfoListEquiv InfoReal;
  gchar *dispURL=NULL, tname[129];
  gchar *routine = "ObitDConCleanVisGetParms";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Sky model type */
  InfoReal.itg = (olong)OBIT_SkyModel_Fastest; type = OBIT_long;
  ObitInfoListGetTest(in->info, "Mode", &type, dim, &InfoReal);
  in->modelMode = InfoReal.itg;

  /* Restore image when done? */
  ObitInfoListGetTest(in->info, "doRestore", &type, dim, &in->doRestore);

  /* Cross restore images when done? */
  in->doXRestore = in->doRestore;
  ObitInfoListGetTest(in->info, "doXRestore", &type, dim, &in->doXRestore);

  /* Flatten image when done? */
  ObitInfoListGetTest(in->info, "doFlatten", &type, dim, &in->doFlatten);

  /* Weight data? */
  ObitInfoListGetTest(in->info, "doWeight", &type, dim, &in->doWeight);

  /* Make beams? */
  in->doBeam = TRUE;
  ObitInfoListGetTest(in->info, "doBeam", &type, dim, &in->doBeam);

  /* Allow recentering? */
  in->doRecenter = TRUE;
  ObitInfoListGetTest(in->info, "doBeam", &type, dim, &in->doRecenter);

 /* Flux level to reuse in CCs */
  ObitInfoListGetTest(in->info, "reuseFlux", &type, dim, &in->reuseFlux);

  /* auto center minimum flux density */
  ObitInfoListGetTest(in->info, "autoCen", &type, dim, &in->autoCen);

  /* Image display? */
  if (!in->display) {
    ObitInfoListGetP(in->info, "dispURL", &type, dim, (gpointer)&dispURL);
    /* dispURL not NULL terminated */
    if (dispURL) {for (i=0; i<dim[0]; i++) tname[i] = dispURL[i]; tname[i]=0;}
    if (dispURL && (strncmp(tname, "None", 4))) 
      in->display = ObitDisplayCreate("Display", tname, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

} /* end ObitDConCleanVisGetParms */

/**
 * Set default CLEAN windows in mosaic
 * If mosaic member  Radius>0 then make round boxes on Fly's eye field
 * with this radius, else use rectangular box including all but outer 5 pixels
 * On outlier fields, use rectangular box of width OutlierSize.
 * If CLEANBox defined in in->info then its contents are used for field 1.
 * Assumes all images in mosaic have descriptors defined.
 * Uses base class function.
 * \param in   The CLEAN object
 * \param err Obit error stack object.
 */
void ObitDConCleanVisDefWindow(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVis *in;
  const ObitDConCleanVisClassInfo *inClass;
  gchar *routine = "ObitDConCleanVisDefWindow";

  /* Cast input to this type */
  in = (ObitDConCleanVis*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanVisIsA(in));

  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */

  /* Call actual function */
  inClass->ObitDConCleanDefWindow((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitDConCleanVisDefWindow */

/**
 * Subtract components from uv data.
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanVisSub(ObitDConCleanVis *in, ObitErr *err)
{
  olong i, field, ver;
  ObitSkyModelType modelType = OBIT_SkyModel_Comps;
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL;
  ObitIOCode retCode;
  gboolean Fl=FALSE, doCalSelect;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat ftemp;
  olong *itemp, jtemp, nfield=0;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitDConCleanVisSub";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanVisIsA(in));

  /* Setup SkyModel parameters */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(in->skyModel->info, "Mode", OBIT_long, dim, &in->modelMode);
  ObitInfoListAlwaysPut(in->skyModel->info, "ModelType", OBIT_long, dim, &modelType);
  jtemp = in->CCver;
  ObitInfoListAlwaysPut(in->skyModel->info, "CCVer", OBIT_long, dim, &jtemp);
  /* Disable any value of minFlux to suppress infinite recursion */
  ftemp = -1.0e20;
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (in->skyModel->info, "minFlux", OBIT_float, dim, &ftemp);

  /* IF this was an SDI Clean, compress CC's and restore data */
  if (in->doSDI) {
    /* restore initial data */
    in->imager->uvwork = ObitUVCopy (in->SDIdata, in->imager->uvwork, err);
    /* compress CCs */
    if (in->currentField<=0) {
      ObitSkyModelCompressCC (in->skyModel, err);  /* All */
      nfield = in->mosaic->numberImages;
      for (i=0; i<nfield; i++) in->skyModel->startComp[i] = 1;
      for (i=0; i<nfield; i++) in->skyModel->endComp[i]   = 0;
   } else { /* only in->currentField */

      field = in->currentField-1;
      /* Get CC table */
      ver = in->skyModel->CCver[field];
      tempTable = newObitImageTable (in->skyModel->mosaic->images[field],OBIT_IO_ReadOnly, 
				     tabType, &ver, err);
      if ((tempTable==NULL) || (err->error)) Obit_traceback_msg (err, routine, in->name);
      CCTable = ObitTableCCConvert(tempTable);
      tempTable = ObitTableUnref(tempTable);
      
      /* Merge */
      retCode = ObitTableCCUtilMerge (CCTable, CCTable, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->skyModel->name);
      in->Pixels->iterField[field] =  CCTable->myDesc->nrow;
      in->skyModel->startComp[field] = 1;
      for (i=0; i<nfield; i++) in->skyModel->startComp[i] = 1;
      for (i=0; i<nfield; i++) in->skyModel->endComp[i]   = 0;
      CCTable = ObitTableCCUnref (CCTable);
   } /* End Merge field */
    if (err->error) Obit_traceback_msg (err, routine, in->skyModel->name);
  } /* end SDI */

  nfield = in->mosaic->numberImages;
  itemp = ObitMemAlloc(nfield*sizeof(olong));  /* temp. array */
  dim[0] = nfield;
  for (i=0; i<nfield; i++) itemp[i] = in->skyModel->startComp[i];
  ObitInfoListAlwaysPut(in->skyModel->info, "BComp", OBIT_long, dim, itemp);
  if (in->doSDI) {
    /* must subtract everything */
    for (i=0; i<nfield; i++) itemp[i] = 1;
    ObitInfoListAlwaysPut(in->skyModel->info, "BComp", OBIT_long, dim, itemp);
    for (i=0; i<nfield; i++) itemp[i] = 0;
  } else { /* normal CLEAN */
    for (i=0; i<nfield; i++) itemp[i] = in->Pixels->iterField[i];
  }
  ObitInfoListAlwaysPut(in->skyModel->info, "EComp", OBIT_long, dim, itemp);
  itemp = ObitMemFree(itemp);  /* Deallocate */

  /* Subtract Current model */
  doCalSelect = FALSE;
  ObitInfoListGetTest (in->imager->uvwork->info, "doCalSelect", &type, dim, &doCalSelect);
  dim[0] = dim[1] = dim[2] = 1;  /* Grumble, grumble  */
  ObitInfoListAlwaysPut (in->imager->uvwork->info, "doCalSelect",OBIT_bool, dim, &Fl);
  ObitSkyModelSubUV(in->skyModel, in->imager->uvwork, in->imager->uvwork, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  dim[0] = dim[1] = dim[2] = 1;  /* Grumble, grumble  */
  ObitInfoListAlwaysPut (in->imager->uvwork->info, "doCalSelect",OBIT_bool, dim, &doCalSelect);

  /* Update CC counts */
  for (i=0; i<in->mosaic->numberImages; i++) {
    in->skyModel->startComp[i] = in->skyModel->endComp[i]+1;
  }

} /* end ObitDConCleanVisSub */

/**
 * Pick next field to clean and create residual image in mosaic member
 * The selected field is in->currentField.
 * Adopted from the AIPSish CLBSTF (QOOP:QCLEAN.FOR)
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 * \return TRUE iff reached minimum flux density or max. number  comp.
 */
gboolean ObitDConCleanVisPickNext(ObitDConCleanVis *in, ObitErr *err)
{
  olong i, best, second, loopCheck;
  gboolean *fresh, doBeam, done=TRUE;
  ofloat sumwts;
  ObitImage *theBeam=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *routine = "ObitDConCleanVisPickNext";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanVisIsA(in));

  /* Check if reached max number of components and some done */
  if ((in->Pixels->currentIter >= in->Pixels->niter) && 
      (in->Pixels->currentIter>0)) return done;

  fresh = ObitMemAlloc0(in->nfield*sizeof(gboolean));
  /* First time? */
  if (in->Pixels->currentIter<=0)
    for (i=0; i<in->nfield; i++) fresh[i] = TRUE;
  else
    for (i=0; i<in->nfield; i++) fresh[i] = FALSE;

  /* Make sure all fields initialized */
  for (i=0; i<in->nfield; i++) {
    if (in->maxAbsRes[i] < 0.0) {
      doBeam = in->doBeam;
      /* Make sure beam is made - is SUMWTS present? */
      theBeam = (ObitImage*)in->mosaic->images[i]->myBeam;
      if (!ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts))
	doBeam = TRUE;
      /* Make residual image - get statistics */
      MakeResidual(in, i+1, doBeam, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, done);
      fresh[i] = TRUE;
    }
  } /* end loop initializing fields */

  /* Check if reached max number of components */
  if (in->Pixels->currentIter >= in->Pixels->niter) return done;

  /* Ignore fields already known to be finished, see if all done */
  done = TRUE;
  for (i=0; i<in->nfield; i++) {
    if (in->maxAbsRes[i] <= in->minFlux[i]) in->quality[i] = 0.0;
    else done = FALSE;
  }
  if (done) {ObitMemFree(fresh); return done;}  /* anything left? */

  /* Find current estimates best and second best field field */
  WhosBest (in, &best, &second);

  /* Verify that it's still best which may require iteration */
  done = FALSE;
  loopCheck = 0;
  while (!done) {
    /* If not just created recalculate */
    if (!fresh[best]) {

      /* Don't loop forever */
      loopCheck++;
      if (loopCheck>2*in->nfield) break;

      /* Make residual image - get statistics */
      MakeResidual(in, best+1, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, done);
      fresh[best] = TRUE;
      /* Ignore fields with max residual less than min. */
      if (in->maxAbsRes[best] < in->minFlux[best]) in->quality[best] = 0.0;
    } /* end remake */

    if (in->nfield==1) break; /* If only one, the decision is easy */

    /* Still best? */
    if (in->quality[best] >= in->quality[second]) break;

    /* try again */
    if ((in->quality[best]>0.0) && (in->prtLv>1)) 
      Obit_log_error(err, OBIT_InfoWarn, 
		     "%s: field %d (%g) not as good as second %g",
		     routine, best+1, in->quality[best],  in->quality[second]);

    /* Find current estimates best and second best field field */
    WhosBest (in, &best, &second);

  } /* end loop finding best field */
  
  /* publish decision */
  in->currentField = best + 1;

  /* cleanup */
  fresh = ObitMemFree( fresh);

  /* highest abs value ever found in a CLEAN window */
  in->peakFlux = MAX (in->peakFlux, in->maxAbsRes[best]);

  /* We're done if best field has maxAbsRes<minFlux */
  done = in->maxAbsRes[best] <= in->minFlux[best];
  return done;
} /* end ObitDConCleanVisPickNext */

/**
 * Determine desirability of field for clean
 * Adopted from the AIPSish CLOFNB (QOOP:QCLEAN.FOR)
 * \param in   The object to deconvolve
 * \param field field number (1-rel) to test
 * \param err Obit error stack object.
 * \return quality measure, higher is more desirable
 */
ofloat ObitDConCleanVisQuality(ObitDConCleanVis *in, olong field, 
			       ObitErr *err)
{
  ofloat out = -1.0;
  gchar *routine = "ObitDConCleanVisQuality";

  /* error checks */
  if (err->error) return out;
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
    return out;
  }
  /* Get value */
  out = 0.95 * in->maxAbsRes[field-1] + 0.05*in->avgRes[field-1];

  return out;
} /* end ObitDConCleanVisQuality */

/**
 * See if an image needs to be remade because a source which exceeds
 * the flux  threshold is not centered (as determined by moments)
 * on the reference pixel (within toler pixel).
 * A new (96x96) field is added centered on the offending source and a negative
 * clean window added to the position of the source in its original window.
 * Avoid duplicates of the same source and ensure that all occurances of this 
 * source in any exant field has a negative clean window added.
 * Multiple centering sources per facet are allowed
 * A boolean entry "autoCenField" with value True is added to the info member of any
 * image members added to the mosaic member.
 * Routine originally translated from the AIPSish VLAUTIL.FOR/VLCCIN  
 * \param in  DConCleanVis to check
 *  Values on info member:
 * \li "restartFlux" OBIT_float scalar = Minimum brightness flux for CLEAN 
 *      restart, def infinite
 * \li "CCVer"      OBIT_long    scalar = CC table version to use, default 1
 * \li "toler" OBIT_float scalar = Tolerance to accept as on a pixel (cells)
 *      def 0.01
 * \param uvdata UV data being imaged.
 * \param err    Error/message stack
 * \return TRUE if image needs to be remade, generally meaning something 
 * exceeded threshold.
 */
gboolean ObitDConCleanVisReimage (ObitDConCleanVis *in, ObitUV* uvdata, 
				  ObitErr* err) 
{
  gboolean redo = FALSE;
  ObitTableCC *CCTab=NULL;
  ObitImageDesc *imDesc, *imDesc2;
  ObitImageMosaic *mosaic = in->mosaic;
  ofloat freset, tol;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  ObitDConCleanWindowType otype;
  olong   nfield, ifield, jfield, itemp, noParms, nccpos, nx, ny, nplane, nprior;
  olong  CCVer, newField, win[3], inaxes[2], *owin;
  ofloat tmax, xcenter, ycenter, xoff, yoff, radius, cells[2], pixel[2], opixel[2];
  ofloat xcen, ycen, RAShift, DecShift, deltax, deltay, delta, *farray;
  odouble pos[2], RAPnt, DecPnt;
  gboolean done, outside=FALSE, facetDone, clear, Tr=TRUE;
  gchar *routine = "ObitDConCleanVisReimage";

  /* Error checks */
  if (err->error) return redo;  /* previous error? */
  g_assert(ObitDConCleanVisIsA(in));

  /* Number of fields */
  nfield = mosaic->numberImages;
  nprior = nfield;

  /* Get cellsize */
  cells[0] =  fabs(mosaic->xCells); cells[1] = fabs(mosaic->yCells);

  /* Consider components within 2.5  cells  */
  radius = 2.5 * cells[0];

  /* Tolerance for match to pixel */
  tol = 0.01;
  ObitInfoListGetTest(mosaic->info, "toler", &type, dim, &tol);

  /* Flux level to reuse in CCs */
  freset = 1.0e20;
  ObitInfoListGetTest(mosaic->info, "restartFlux", &type, dim, &freset);

  /* CC table(s) */
  itemp = 1;
  ObitInfoListGetTest(mosaic->info, "CCVer", &type, dim, &itemp);
  CCVer = itemp;

  /* Loop over fields */
  for (ifield=0; ifield<nfield; ifield++) { /* loop 500 */

    /* Open image in case header needs update */
    ObitImageOpen (mosaic->images[ifield], OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Make temporary CC table */
    noParms = 0;
    CCTab = newObitTableCCValue ("Temp CC", (ObitData*)mosaic->images[ifield],
				 &CCVer, OBIT_IO_ReadWrite, noParms, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Loop over peaks in facet */
    facetDone = FALSE;
    clear = FALSE;  /* Clear CC table */
    while (!facetDone) {
      
      /* Check Table for peak */
      nccpos = CCTab->myDesc->nrow;
      ObitImageMosaicMaxCC (CCTab, nccpos, radius, &tmax, &xcenter, &ycenter, &xoff, &yoff, err);
      if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

      facetDone = (tmax < freset);  /* Found more? */
      if (facetDone) break;

      /* Position of bright peak */
      imDesc = mosaic->images[ifield]->myDesc;
      pixel[0] = imDesc->crpix[0] + ((xcenter+xoff) / imDesc->cdelt[0]);
      pixel[1] = imDesc->crpix[1] + ((ycenter+yoff) / imDesc->cdelt[1]);
      ObitImageDescGetPos(imDesc, pixel, pos, err);
      
      /* Check that this position has not just been added from another field. */
      done = FALSE;
      for (jfield=nprior; jfield<mosaic->numberImages; jfield++) { 
	imDesc2 = mosaic->images[jfield]->myDesc;
	deltax = (imDesc2->crval[0]-pos[0]) / imDesc->cdelt[0];
	deltay = (imDesc2->crval[1]-pos[1]) / imDesc->cdelt[1];
	delta = sqrt (deltax*deltax + deltay*deltay);
	/* If center within 2 pixel of the previous case call it the same */
	done = done || (delta<2.0);
	if (done) break;	
      } /* End loop checking for a previous occurance */
	
      if (done) goto doNext;  /* Has this one already been done? */
      
      /* See if peak is near the center of the image - if within 1/2 pixel
	 and exceeds threshold adjust field and mark as autoCenterField */ 
      if (((fabs(xcenter+xoff)<0.5*cells[0]) && 
	   (fabs(ycenter+yoff)<0.5*cells[1])) ) {
	
	/* Mark as an autoCenter Image */
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(mosaic->images[ifield]->info, 
			      "autoCenField", OBIT_bool, dim, &Tr);
	
	ObitImageDescGetPoint (imDesc, &RAPnt, &DecPnt);
	ObitSkyGeomShiftXY (RAPnt, DecPnt, ObitImageDescRotate(imDesc),
			    pos[0], pos[1], &RAShift, &DecShift);
	/* Update mosaic and shift in other places */
	mosaic->RAShift[ifield]  = RAShift;
	mosaic->DecShift[ifield] = DecShift;
	imDesc->crval[imDesc->jlocr] = pos[0];
	imDesc->crval[imDesc->jlocd] = pos[1];
	imDesc->xshift = RAShift;
	imDesc->yshift = DecShift;
	ObitInfoListGetP(uvdata->info, "xShift", &type, dim, (gpointer*)&farray);
	farray[ifield]  = RAShift;
	ObitInfoListAlwaysPut(uvdata->info, "xShift", type, dim, farray);
	ObitInfoListGetP(uvdata->info, "yShift", &type, dim, (gpointer*)&farray);
	farray[ifield]  = DecShift;
	ObitInfoListAlwaysPut(uvdata->info, "yShift", type, dim, farray);
	mosaic->images[ifield]->myStatus = OBIT_Modified; /* Update header */
	clear = TRUE;  /* Clear CC table */
	redo  = TRUE;
	
	/* Tell about it */
	if (in->prtLv>1) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Modify autoCenter field %d by %5.2f %5.3f pixels",
			 ifield+1,xoff/imDesc->cdelt[imDesc->jlocr],
			 yoff/imDesc->cdelt[imDesc->jlocd]);
	}
	done = TRUE;
      } /* end bright source at field center */
      if (done) goto doNext;  /* Has this one already been done? */
      
      /* See if shift needed */
      if (((fabs(xcenter+xoff)>tol*cells[0]) || 
	   (fabs(ycenter+yoff)>tol*cells[1])) ) {
	imDesc = mosaic->images[ifield]->myDesc;
	clear = TRUE;  /* Clear CC table */
	redo  = TRUE;
	
	/* Add field to mosaic */
	ObitImageDescGetPoint (imDesc, &RAPnt, &DecPnt);
	ObitSkyGeomShiftXY (RAPnt, DecPnt, ObitImageDescRotate(imDesc),
			    pos[0], pos[1], &RAShift, &DecShift);
	nx = ny = ObitFFTSuggestSize (96); nplane = 1;
	ObitImageMosaicAddField (in->mosaic, uvdata, nx, ny, nplane, 
				 RAShift, DecShift, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Mark as an autoCenter Image */
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(mosaic->images[mosaic->numberImages-1]->info, 
			      "autoCenField", OBIT_bool, dim, &Tr);
	
	/* Put unwindow on this position in all prior fields in which it occured */
	for (jfield=0; jfield<in->window->nfield; jfield++) { 
	  /* is pos in this field */
	  imDesc2 = mosaic->images[jfield]->myDesc;
	  if (!ObitImageDescOverlap(imDesc, imDesc2, err)) continue;
	  
	  /* Check if in outer window */
	  ObitImageDescGetPixel (imDesc2, pos, opixel, err);
	  if (ObitDConCleanWindowInfo(in->window, jfield+1, -1, &otype, &owin, err))
	    {
	      if (otype==OBIT_DConCleanWindow_rectangle) {
		outside = (opixel[0]<owin[0]) || (opixel[0]>owin[2]) ||
		  (opixel[1]<owin[1]) || (opixel[1]>owin[3]);
	      } else if (otype==OBIT_DConCleanWindow_round) {
		deltax = owin[1] - opixel[0];
		deltay = owin[2] - opixel[1];
		delta = sqrt (deltax*deltax + deltay*deltay);
		outside = delta > owin[0];
	      }
	    } else {  /* No outer window given  */
	      outside = (opixel[0]<1.0) || (opixel[0]>imDesc2->inaxes[0]) ||
		(opixel[1]<1.0) || (opixel[1]>imDesc2->inaxes[1]);
	    }
	  if  (err->error) 
	    Obit_traceback_val (err, routine, mosaic->images[jfield]->name, redo);
	  if (outside) continue;
	  
	  /* OK - add unbox */
	  win[0] = (nx/2)-15; win[1] = (olong)(opixel[0]+0.5); win[2] = (olong)(opixel[1]+0.5); 
	  /* If this is a previous autoCenter window make unbox smaller */
	  if ((jfield+1)>nprior) win[0] = 5;
	  ObitDConCleanWindowAdd (in->window, jfield+1, OBIT_DConCleanWindow_unround,
				  win, err);	
	} /* end loop adding unboxes */
	
	/* Add field to window */
	inaxes[0] = nx; inaxes[1] = ny;
	newField = ObitDConCleanWindowAddField (in->window, inaxes, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Add inner and outer windows */
	xcen = mosaic->images[newField-1]->myDesc->crpix[0];
	ycen = mosaic->images[newField-1]->myDesc->crpix[1];
	win[0] = 5; win[1] = (olong)(xcen+0.5); win[2] = (olong)(ycen+0.5); 
	ObitDConCleanWindowAdd (in->window, newField, OBIT_DConCleanWindow_round,
				win, err);
	win[0] = (nx/2)-10; win[1] = (olong)(xcen+0.5); win[2] = (olong)(ycen+0.5); 
	ObitDConCleanWindowOuter (in->window, newField, OBIT_DConCleanWindow_round,
				  win, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Add field to Clean */
	ObitDConCleanVisAddField (in, uvdata, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Add field to SkyModel */
	ObitSkyModelAddField (in->skyModel, err);
	if (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Tell about it */
	if (in->prtLv>1) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Add field due to strong source %d offset %5.2f %5.3f pixels, peak %f",
			 ifield+1,xoff/imDesc->cdelt[imDesc->jlocr],
			 yoff/imDesc->cdelt[imDesc->jlocd], tmax);
	  Obit_log_error(err, OBIT_InfoErr, 
			 "strong source position %lf %lf shift %f %f",
			 pos[0], pos[1], RAShift, DecShift);
	}
      } /* End redo this image */

      /* Zero CC entries around old peak in scratch CC table */
    doNext:
      nccpos = CCTab->myDesc->nrow;
      ObitImageMosaicFlagCC (CCTab, nccpos, radius, xcenter, ycenter, err);
      if  (err->error) Obit_traceback_val (err, routine, CCTab->name, redo);

    } /* end loop over multiple peaks */

    /* Need to clear CC table */
    if (clear) ObitTableClearRows ((ObitTable*)CCTab, err); /* Remove the entries and redo */
    if (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
    
    /* Delete temporary objects */
    CCTab = ObitTableCCUnref(CCTab);

    /* Close/update image */
    ObitImageClose(mosaic->images[ifield], err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
  } /* end loop  L500: */

  ObitErrLog(err);  /* Any messages */
  return redo;
} /* end of routine ObitDConCleanVisReimage */ 

/**
 * Add a field to object structures
 * \param in  DConCleanVis to expand
 * \param uvdata UV data being imaged.
 * \param err    Error/message stack
 */
void ObitDConCleanVisAddField (ObitDConCleanVis *in, ObitUV* uvdata, 
			       ObitErr* err) 
{
  olong newField = 0;
  olong i, oldField;
  ofloat *ftemp;
  gchar *routine = "ObitDConCleanVisAddField";

  /* Error checks */
  if (err->error) return;  /* previous error? */
  g_assert(ObitDConCleanVisIsA(in));

   /* field to add */
  oldField = in->nfield;
  newField = oldField+1;
  in->nfield = newField;

  /* Check that mosaic same size */
  Obit_return_if_fail ((in->nfield == in->mosaic->numberImages),
		       err, "%s: CLEAN and mosaic differnet number of fields %d %d",
		       routine, in->nfield, in->mosaic->numberImages);

  /* Resize/copy arrays */
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean Loop gain");
  for (i=0; i<oldField; i++) ftemp[i] = in->gain[i]; ftemp[i] = 0.0; 
  in->gain = ObitMemFree(in->gain);
  in->gain = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean minFlux");
  for (i=0; i<oldField; i++) ftemp[i] = in->minFlux[i]; ftemp[i] = 0.0; 
  in->minFlux = ObitMemFree(in->minFlux);
  in->minFlux = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean factor");
  for (i=0; i<oldField; i++) ftemp[i] = in->factor[i]; ftemp[i] = 0.0; 
  in->factor = ObitMemFree(in->factor);
  in->factor = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean quality");
  for (i=0; i<oldField; i++) ftemp[i] = in->quality[i]; ftemp[i] = 0.0; 
  in->quality = ObitMemFree(in->quality);
  in->quality = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean max res");
  for (i=0; i<oldField; i++) ftemp[i] = in->maxAbsRes[i]; ftemp[i] = 0.0; 
  in->maxAbsRes = ObitMemFree(in->maxAbsRes);
  in->maxAbsRes = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean avg res");
  for (i=0; i<oldField; i++) ftemp[i] = in->avgRes[i]; ftemp[i] = 0.0; 
  in->avgRes = ObitMemFree(in->avgRes);
  in->avgRes = ftemp;

} /* end of routine ObitDConCleanVisAddField */ 

/**
 * See if an AutoCenter field is properly centered.
 * The centroid of the CLEAN components of autoCenter fields are determined.
 * If they are not at the reference pixel to within toler, the position is adjusted.
 * Fields in the mosaic member of in are considered autoCenter fields if the image
 * info member contains a boolean entry "autoCenField" with value True.
 * Any images recentered will have the CLEAN components cleared 
 * Images should be remade if any positions are modified.
 * \param in  DConCleanVis to check
 *  Values on info member:
 * \li "CCVer"      OBIT_long    scalar = CC table version to use, default 1
 * \li "toler" OBIT_float scalar = Tolerance to accept as on a pixel (cells)
 *      def 0.001
 * \param uvdata UV data
 * \param err    Error/message stack
 * return TRUE if a position was modified
 */
gboolean ObitDConCleanVisRecenter (ObitDConCleanVis *in, ObitUV* uvdata, 
				   ObitErr* err) 
{
  gboolean redo = FALSE;
  ObitTableCC *CCTab=NULL;
  ObitImageDesc *imDesc;
  ObitImageMosaic *mosaic = in->mosaic;
  ofloat tol;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong   nfield, ifield, itemp, noParms, nccpos, nprior;
  olong  CCVer, highVer;
  ofloat tmax, xcenter, ycenter, xoff, yoff, radius, cells[2], pixel[2];
  ofloat RAShift, DecShift, *farray;
  odouble pos[2], RAPnt, DecPnt;
  gboolean autoCen, want;
  gchar *routine = "ObitDConCleanVisRecenter";

  /* Error checks */
  if (err->error)     return redo;  /* previous error? */
  g_assert(ObitDConCleanVisIsA(in));

  /* Number of fields */
  nfield = mosaic->numberImages;
  nprior = nfield;

   /* Get cellsize */
  cells[0] =  fabs(mosaic->xCells); cells[1] = fabs(mosaic->yCells);

  /* Consider components within 1.5  cells  */
  radius = 1.5 * cells[0];

  /* Tolerance for match to pixel */
  tol = 0.01;
  ObitInfoListGetTest(mosaic->info, "toler", &type, dim, &tol);

  /* CC table(s) */
  itemp = 1;
  ObitInfoListGetTest(mosaic->info, "CCVer", &type, dim, &itemp);
  CCVer = itemp;

  /* Loop over fields */
  for (ifield=0; ifield<nfield; ifield++) { /* loop 500 */

    /* Is this an autoCenter field with CLEAN components? */
    autoCen = FALSE;
    ObitInfoListGetTest(mosaic->images[ifield]->info, "autoCenField", &type, dim, &autoCen);
    if (!autoCen) continue;

    /* Open image in case header needs update */
    ObitImageOpen (mosaic->images[ifield], OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Make sure there is already a CC table */
    highVer =  ObitTableListGetHigh (mosaic->images[ifield]->tableList, "AIPS CC");
    if (highVer<=0) goto closeit;

    /* Get CC table */
    noParms = 0;
    CCTab = newObitTableCCValue ("Temp CC", (ObitData*)mosaic->images[ifield],
				 &CCVer, OBIT_IO_ReadWrite, noParms, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Is this an autoCenter field with CLEAN components? */
    autoCen = FALSE;
    ObitInfoListGetTest(mosaic->images[ifield]->info, "autoCenField", &type, dim, &autoCen);
    if (CCTab->myDesc->nrow>0) {
      
      /* Check Table */
      nccpos = CCTab->myDesc->nrow;
      ObitImageMosaicMaxCC (CCTab, nccpos, radius, &tmax, &xcenter, &ycenter, 
			    &xoff, &yoff, err);
      if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
      
      /* See if shift needed - if more then 2 cells this is probably a mistake */
      want = ((fabs(xcenter+xoff)>tol*cells[0]) || (fabs(ycenter+yoff)>tol*cells[1]));
      want = want && (fabs(xcenter)<2.0) && (fabs(ycenter)<2.0);
      if (want) {
	imDesc = mosaic->images[ifield]->myDesc;
	
	/* Position of peak */
	pixel[0] = imDesc->crpix[0] + ((xcenter+xoff) / imDesc->cdelt[0]);
	pixel[1] = imDesc->crpix[1] + ((ycenter+yoff) / imDesc->cdelt[1]);
	ObitImageDescGetPos(imDesc, pixel, pos, err);
	ObitImageDescGetPoint (imDesc, &RAPnt, &DecPnt);
	ObitSkyGeomShiftXY (RAPnt, DecPnt, ObitImageDescRotate(imDesc),
			    pos[0], pos[1], &RAShift, &DecShift);
	/* Update mosaic and shift in other places */
	mosaic->RAShift[ifield]  = RAShift;
	mosaic->DecShift[ifield] = DecShift;
	imDesc->crval[imDesc->jlocr] = pos[0];
	imDesc->crval[imDesc->jlocd] = pos[1];
	imDesc->xshift = RAShift;
	imDesc->yshift = DecShift;
	ObitInfoListGetP(uvdata->info, "xShift", &type, dim, (gpointer*)&farray);
	farray[ifield]  = RAShift;
	ObitInfoListAlwaysPut(uvdata->info, "xShift", type, dim, farray);
	ObitInfoListGetP(uvdata->info, "yShift", &type, dim, (gpointer*)&farray);
	farray[ifield]  = DecShift;
	ObitInfoListAlwaysPut(uvdata->info, "yShift", type, dim, farray);
	mosaic->images[ifield]->myStatus = OBIT_Modified; /* Update header */
	ObitTableClearRows ((ObitTable*)CCTab, err); /* Remove the entries and redo */
	redo = TRUE;

        /* Tell about it */
	if (in->prtLv>1) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Modify autoCenter field %d by %5.2f %5.3f pixels",
			 ifield+1,xoff/imDesc->cdelt[imDesc->jlocr],
			 yoff/imDesc->cdelt[imDesc->jlocd]);
	}
	
      } /* End redo this image */
    } /* end autoCenter field */
    
    /* Delete temporary table */
    CCTab = ObitTableCCUnref(CCTab);
    
    /* Close/update image */
  closeit:
    ObitImageClose(mosaic->images[ifield], err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
  } /* end loop  L500: */
  
  return redo;
} /* end of routine ObitDConCleanVisRecenter */ 

/**
 * Filter any weak, isolated components
 * \param in  DConCleanVis to check
 *  Values on info member:
 * \li "CCVer"      OBIT_long    scalar = CC table version to use, default 1
 * \param filter, filtering parameters:
 *        [0] =  radius    Radius within which to consider components. (deg) 
 *               if this is 0.0 nothing is done
 *        [1] =  minFlux   Minimum acceptable summed flux (Jy)
 * \param err    Error/message stack
 * return TRUE if any components were removed, i.e. reimaging needed
 */
gboolean ObitDConCleanVisFilter (ObitDConCleanVis *in, ofloat filter[2], 
				 ObitErr* err) 
{
  gboolean redo = FALSE;
  ObitTableCC *CCTab=NULL;
  ObitImageDesc *imDesc;
  ObitImageMosaic *mosaic = in->mosaic;
  gboolean some;
  ofloat radius;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong   nfield, ifield, itemp, noParms;
  olong  CCVer, highVer;
  gchar *routine = "ObitDConCleanVisFilter";

  /* Error checks */
  if (err->error) return redo;  /* previous error? */
  if (filter[0]==0.0) return redo; 
  g_assert(ObitDConCleanVisIsA(in));

  /* Tell about it */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Filtering comps < %f within %f pixel",
		   filter[0], filter[1]);
    ObitErrLog(err); 
  }

  /* Number of fields */
  nfield = mosaic->numberImages;

  /* CC table(s) */
  itemp = 1;
  ObitInfoListGetTest(mosaic->info, "CCVer", &type, dim, &itemp);
  CCVer = itemp;

  /* Loop over fields */
  for (ifield=0; ifield<nfield; ifield++) { /* loop over fields */

    /* Open image in case header needs update */
    ObitImageOpen (mosaic->images[ifield], OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Make sure there is already a CC table */
    highVer =  ObitTableListGetHigh (mosaic->images[ifield]->tableList, "AIPS CC");
    if (highVer<=0) goto closeit;

    /* Get CC table */
    noParms = 0;
    CCTab = newObitTableCCValue ("Temp CC", (ObitData*)mosaic->images[ifield],
				 &CCVer, OBIT_IO_ReadWrite, noParms, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Filter, convert radius to degrees */
    imDesc = mosaic->images[ifield]->myDesc;
    radius = fabs(imDesc->cdelt[imDesc->jlocr]*imDesc->cdelt[imDesc->jlocd]);
    radius = filter[1]*sqrt(radius);
    some = ObitTableCCUtilFiltCC (CCTab, radius, filter[0], err);
    redo = redo || some;
    if  (err->error) Obit_traceback_val (err, routine, CCTab->name, redo);

    /* Sort table if any removed */
    if (some) {
      ObitTableUtilAbsSort ((ObitTable*)CCTab, "FLUX    ", TRUE,  err);
      if (err->error) Obit_traceback_val (err, routine, CCTab->name, redo);
    }

    /* Delete temporary table */
    CCTab = ObitTableCCUnref(CCTab);
    
    /* Close/update image */
  closeit:
    ObitImageClose(mosaic->images[ifield], err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
  } /* end loop over fields */
  
  /* Tell if found something */
  if (in->prtLv>1) {
    if (redo) Obit_log_error(err, OBIT_InfoErr, 
			     "Some components removed, will remake residuals");
    else Obit_log_error(err, OBIT_InfoErr, 
			     "No components removed");
    ObitErrLog(err); 
  }
  return redo;
} /* end of routine ObitDConCleanVisFilter */ 

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanVisClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanVisClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanVisClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanVisClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanVisClassInfo *theClass = (ObitDConCleanVisClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanVisClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanVisClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanVisGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanVis;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanVisCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanVisClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanVisInit;
  theClass->ObitDConGetParms        = (ObitDConGetParmsFP)ObitDConCleanVisGetParms;
  theClass->ObitDConDeconvolve      = (ObitDConDeconvolveFP)ObitDConCleanVisDeconvolve;
  theClass->ObitDConCleanSub        = (ObitDConCleanSubFP)ObitDConCleanVisSub;
  theClass->ObitDConCleanVisPickNext  = (ObitDConCleanVisPickNextFP)ObitDConCleanVisPickNext;
  theClass->ObitDConCleanVisQuality   = (ObitDConCleanVisQualityFP)ObitDConCleanVisQuality;
  theClass->ObitDConCleanVisReimage   = (ObitDConCleanVisReimageFP)ObitDConCleanVisReimage;
  theClass->ObitDConCleanVisAddField  = (ObitDConCleanVisAddFieldFP)ObitDConCleanVisAddField;
  theClass->ObitDConCleanVisRecenter  = (ObitDConCleanVisRecenterFP)ObitDConCleanVisRecenter;
  theClass->ObitDConCleanVisFilter    = (ObitDConCleanVisFilterFP)ObitDConCleanVisFilter;

} /* end ObitDConCleanVisClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanVisInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVis *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->imager   = NULL;
  in->skyModel = NULL;
  in->modelMode = OBIT_SkyModel_Fastest;
  in->doRestore = TRUE;
  in->doXRestore= TRUE;
  in->doFlatten = TRUE;
  in->doWeight  = TRUE;
  in->doBeam    = TRUE;
  in->quality   = NULL;
  in->display   = NULL;
  in->SDIdata   = NULL;
  in->peakFlux  = -1000.0;
  in->reuseFlux = -1.0;
  in->autoCen   =  1.0e20;
} /* end ObitDConCleanVisInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanVis* cast to an Obit*.
 */
void ObitDConCleanVisClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVis *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->imager    = ObitUVImagerUnref(in->imager);
  in->skyModel  = ObitSkyModelUnref(in->skyModel);
  in->quality   = ObitMemFree(in->quality);
  in->display   = ObitDisplayUnref(in->display);
  in->SDIdata   = ObitUVUnref(in->SDIdata);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanVisClear */

/**
 * Make selected residual image and get statistics
 * \param in     The Clean object
 * \param field  Which (1-rel) field
 * \param doBeam If TRUE also make beam
 * \param err    Obit error stack object.
 */
static void  MakeResidual (ObitDConCleanVis *in, olong field, 
			   gboolean doBeam, ObitErr *err)
{
  gboolean doWeight, doFlatten;
  const ObitDConCleanVisClassInfo *inClass;
  ObitUVImagerClassInfo *imgClass = 
    (ObitUVImagerClassInfo*)in->imager->ClassInfo;
  gchar *routine = "MakeResidual";

  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */

  /* Are residuals fresh? */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Make residual image */
  imgClass->ObitUVImagerImage (in->imager, field, doWeight, doBeam, doFlatten, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Get statistics */
  inClass->ObitDConCleanImageStats ((ObitDConClean*)in, field, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Quality measure */
  in->quality[field-1] = ObitDConCleanVisQuality((ObitDConCleanVis*)in, field, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end MakeResidual */

/**
 * Make all residual images and get statistics
 * \param in     The Clean object
 * \param err    Obit error stack object.
 */
static void  MakeAllResiduals (ObitDConCleanVis *in, ObitErr *err)
{
  gboolean doBeam = in->doBeam;
  gboolean doWeight, doFlatten;
  olong i;
  ObitUVImagerClassInfo *imgClass = 
    (ObitUVImagerClassInfo*)in->imager->ClassInfo;
  const ObitDConCleanVisClassInfo *inClass;
  gchar *routine = "MakeAllResiduals";
  
  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */

  /* Turn off things not needed */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Parallel Image images without needing beam */
  imgClass->ObitUVImagerImage (in->imager, 0,  doWeight, doBeam, doFlatten, err);
  
  /* Loop over fields getting statistics */
  for (i=0; i<in->nfield; i++) {
    inClass->ObitDConCleanImageStats ((ObitDConClean*)in, i+1, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Quality measure */
    in->quality[i] = ObitDConCleanVisQuality((ObitDConCleanVis*)in, i+1, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end loop over field */

} /* end MakeAllResiduals */

/**
 * Find current estimates best and second field
 * \param in     The Clean object
 * \param err    Obit error stack object.
 * \param best   [out] 0-rel best field
 * \param second [out] 0-rel second best field
 */
static void WhosBest (ObitDConCleanVis *in, olong *bbest, olong *ssecond)
{
  ofloat testBest, testSecond;
  olong i, best, second;;
  
  best = second = -1; 
  testBest = testSecond = -1.0e20;
  for (i=0; i<in->nfield; i++) {
    /* First best ?*/
    if (in->quality[i]>testBest) {
      /* Move old one to second place? */
      if (testBest>testSecond) {
	testSecond = testBest;
	second = best;
      }
      testBest = in->quality[i];
      best = i;
    } else if (in->quality[i]>testSecond) {
      testSecond = in->quality[i];
      second = i;
    }
    /* Second Best */
  } /* end loop finding current best */
  
  /* results */
  *bbest   = best;
  *ssecond = second;
} /* end WhosBest */

/**
 * Reset Sky Model
 * If reuseFlux member is > 0 then any components in the sky model above 
 * this level are left and TRUE is returned if any exist.
 * Resets the number of components
 * the reuseFlux option uses the Pixels member and the list of CC Tables
 * which must be fully instantiated.
 * Sets BChan, EChan, BIF, EIF to default
 * \param in   The Clean object
 * \param err Obit error stack object.
 * \return true if components in the sky model need to be subtracted
 */
static gboolean ResetSkyModel(ObitDConCleanVis *in, ObitErr *err)
{
  gboolean doSub=FALSE;
  olong ncomp;
  ofloat sum;
  olong i, irow, it;
  gint32 dim[MAXINFOELEMDIM];
  olong *itemp=NULL, nfield;
  ObitTableCCRow *CCRow = NULL;
  gchar *routine = "ObitDConCleanVisResetSkyModel";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return doSub;
  g_assert (ObitDConCleanVisIsA(in));

  /* Reset PixelList parameters */
  ObitDConCleanPxListGetParms (in->Pixels, err);
  if (err->error) goto cleanup;
  in->Pixels->currentIter = 0;
  in->Pixels->totalFlux   = 0.0;

  /* Reset SkyModel parameters */
  nfield = in->mosaic->numberImages;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  for (i=0; i<nfield; i++) in->skyModel->startComp[i] = 1;
  for (i=0; i<nfield; i++) in->skyModel->endComp[i]   = 0;
  itemp = ObitMemAlloc(nfield*sizeof(olong));  /* temp. array */
  dim[0] = nfield;
  for (i=0; i<nfield; i++) itemp[i] = 1;
  ObitInfoListAlwaysPut(in->skyModel->info, "BComp", OBIT_long, dim, itemp);

  /* Channel/IF selection */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  it = 1;
  ObitInfoListAlwaysPut(in->skyModel->info, "BChan", OBIT_long, dim, &it);
  ObitInfoListAlwaysPut(in->skyModel->info, "BIF",   OBIT_long, dim, &it);
  it = 0;
  ObitInfoListAlwaysPut(in->skyModel->info, "EChan", OBIT_long, dim, &it);
  ObitInfoListAlwaysPut(in->skyModel->info, "EIF",   OBIT_long, dim, &it);

  for (i=0; i<nfield; i++) itemp[i] = 0;  /* initial number to subtract */

  /* Check for reuse  */
  if ((in->reuseFlux > 0.0) && (in->Pixels)) {
    in->Pixels->currentIter = 0;
    in->Pixels->totalFlux   = 0.0;
    ncomp = 0;
    sum = 0.0;
    for (i=0; i<nfield; i++) {
      in->Pixels->iterField[i] = 0;
      in->Pixels->fluxField[i] = 0.0;
      
      if (ObitTableCCIsA(in->Pixels->CCTable[i])) {
	ObitTableCCOpen (in->Pixels->CCTable[i], OBIT_IO_ReadWrite, err);
	if (err->error) goto cleanup;
	if (!CCRow) CCRow = newObitTableCCRow (in->Pixels->CCTable[i]);
	for (irow=1; irow<=in->Pixels->CCTable[i]->myDesc->nrow; irow++) {
	  ObitTableCCReadRow (in->Pixels->CCTable[i], irow, CCRow, err);
	  if (err->error) goto cleanup;
	  if (fabs(CCRow->Flux) > in->reuseFlux) {
	    itemp[i] = irow;
	    ncomp++;
	    sum += CCRow->Flux;
	    in->Pixels->iterField[i] = irow;
	    in->Pixels->fluxField[i] += CCRow->Flux;
	    /* Remember this as the brightest point is likely here */
	    in->peakFlux = MAX (in->peakFlux, CCRow->Flux);
	    doSub = TRUE;
	  } else break;
	}

	/* Reset number of rows */
	in->Pixels->CCTable[i]->myDesc->nrow = itemp[i];
	/* The one that counts is in the IO */
	((ObitTableDesc*)(in->Pixels->CCTable[i]->myIO->myDesc))->nrow = itemp[i];
	/* Mark as changed */
	in->Pixels->CCTable[i]->myStatus = OBIT_Modified;
   
	ObitTableCCClose (in->Pixels->CCTable[i], err);
	if (err->error) goto cleanup;
      }
    } /* end loop over fields */
    if ((ncomp>0) && (in->prtLv>1))
      Obit_log_error(err, OBIT_InfoErr,"Restart CLEAN with %d comps with %g Jy",
		     ncomp, sum);
    in->Pixels->currentIter = ncomp;
    in->Pixels->totalFlux   = sum;
  } /* End check for reuse */

  /* Cleanup */
 cleanup:
  ObitInfoListAlwaysPut(in->skyModel->info, "EComp", OBIT_long, dim, itemp);
  CCRow = ObitTableCCRowUnref(CCRow);  
  itemp = ObitMemFree(itemp);  /* Deallocate */
  if (err->error) Obit_traceback_val (err, routine, in->name, doSub);
   return doSub;
} /* end ResetSkyModel */

/**
 * Get pixel array for a given field
 * \param in     The Clean object
 * \param field  Field number (1-rel) in ImageMosaic
 * \param err    Obit error stack object.
 * \return true if components in the sky model need to be subtracted
 */
static ObitFArray* GetFieldPixArray (ObitDConCleanVis *in, olong field, 
				     ObitErr *err)
{
  ObitFArray *usePixels=NULL;
  ObitImage *image=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  i, blc[IM_MAXDIM], trc[IM_MAXDIM];
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gchar *routine = "GetFieldPixArray";
  
  /* error checks */
  if (err->error) return usePixels;
  g_assert (ObitDConCleanIsA(in));

  /* Set output to full image, plane at a time */
  blc[0] = blc[1] = 1;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  trc[0] = trc[1] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];

  /* Get image */
  image = in->mosaic->images[field-1];
  
  /* Set input to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc); 
  dim[0] = 1;
  ObitInfoListAlwaysPut (image->info, "IOBy", OBIT_long, dim, &IOsize);
  
  retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, usePixels);
  
  retCode = ObitImageRead (image, image->image->array, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, usePixels);
  
  retCode = ObitImageClose (image, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, usePixels);

  /* Make copy */
  usePixels    = ObitFArrayCopy(image->image, NULL, err);
  image->image = ObitFArrayUnref(image->image);  /* Free buffer */
  if (err->error) Obit_traceback_val (err, routine, image->name, usePixels);
 
  return usePixels;
} /* end GetFieldPixArray */

/**
 * Low accuracy subtract pixels from image for a given field.
 * Uses BeamPatch to subtract list of components
 * \param in       The Clean object
 * \param field    Field number (1-rel) in ImageMosaic
 * \param newCC    Start CC for subtraction
 * \param pixarray pixelarray to subtract from
 * \param err    Obit error stack object.
 */
static void SubNewCCs (ObitDConCleanVis *in, olong field, olong newCC, 
		       ObitFArray *pixarray, ObitErr *err)
{
  olong i, ver, ncc, pos1[2], pos2[2], offset[2], len;
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL;
  ObitImageDesc *imDesc;
  ObitFArray *comps=NULL;
  ofloat ftemp, parms[20], flux;
  gchar *tabType = "AIPS CC";
  gchar *routine = "SubNewCCs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Something to do? */
  ncc     =  in->Pixels->iterField[field-1];
  if (newCC>=ncc) return;

  /* Get CC table */
  ver = in->skyModel->CCver[field-1];
  tempTable = newObitImageTable (in->skyModel->mosaic->images[field-1],OBIT_IO_ReadOnly, 
				 tabType, &ver, err);
  if ((tempTable==NULL) || (err->error)) Obit_traceback_msg (err, routine, in->name);
  CCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);

  /* Get selected CCs after compression */
  comps   = ObitTableCCUtilMergeSel (CCTable, newCC, ncc, parms, err);
  CCTable = ObitTableUnref(CCTable);
  Obit_return_if_fail ((comps!=NULL),
		       err, "%s: Error merging CLEAN components for field %d",
		       routine, in->nfield);

  /* Loop over CCs - assumes all points */
  len       = comps->naxis[0];
  pos2[0]   = in->BeamPatch->naxis[0]/2; 
  pos2[1]   = in->BeamPatch->naxis[1]/2;  /* Center of beam patch */
  offset[0] = pixarray->naxis[0]/2; 
  offset[1] = pixarray->naxis[1]/2;       /* Center of image */
  imDesc  = in->mosaic->images[field-1]->myDesc;
  for (i=0; i<comps->naxis[1]; i++) {
    /* Get 0-rel pixel numbers */
    ftemp = comps->array[1+i*len]/imDesc->cdelt[0];
    if (ftemp>0.0) ftemp += 0.5;
    else           ftemp -= 0.5;
    pos1[0] = offset[0] + (olong)(ftemp); 
    ftemp = comps->array[2+i*len]/imDesc->cdelt[1];
    if (ftemp>0.0) ftemp += 0.5;
    else           ftemp -= 0.5;
    pos1[1] = offset[1] + (olong)(ftemp); 
    flux = comps->array[i*len];
    ObitFArrayShiftAdd(pixarray, pos1, in->BeamPatch, pos2, -flux, pixarray);
  }

  comps = ObitFArrayUnref(comps); /* Cleanup */
  
} /* end SubNewCCs */

