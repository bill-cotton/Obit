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

#include "ObitImageMosaic.h"
#include "ObitImageUtil.h"
#include "ObitDConClean.h"
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


/*---------------Private structures----------------*/
/* Image model subtraction threaded function argument */
typedef struct {
  /* CLEAN Object */
  ObitDConCleanVis *in;
  /* Input plane pixel data */
  ObitFArray *inData;
  /* Field number (1-rel) of data in inData  */
  olong      ofield;
  /* Input array of compressed CLEAN component arrays */
  ObitFArray **inComp;
  /* Number of fields in inComp, fields in in->currentFields  */
  olong      nfield;
  /* thread number, <0 -> no threading  */
  olong      ithread;
  /* Obit Thread object */
  ObitThread  *thread;
  /* Obit error stack object */
  ObitErr    *err;
} ImSubFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanVisInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanVisClear (gpointer in);

/** Private: (re)make residuals. */
static void  MakeResiduals (ObitDConCleanVis *in, olong *fields, 
			    gboolean doBeam, ObitErr *err);

/** Private: (re)make all residuals. */
static void  MakeAllResiduals (ObitDConCleanVis *in, ObitErr *err);

/** Private: Find best 3D residual image. */
static void  WhosBest (ObitDConCleanVis *in, olong *best, olong *second);

/** Private: Find best 2D residual image. */
static void  WhosBest2D (ObitDConCleanVis *in, ofloat autoCenFlux, olong *best);

/** Private: Priority order for making residual images. */
static void  OrderImage (ObitDConCleanVis *in, gboolean *fresh, ofloat autoCenFlux, 
			 olong *fldList);

/** Private: Priority order for CLEANing */
static void  OrderClean (ObitDConCleanVis *in, gboolean *fresh, ofloat autoCenFlux, 
			 olong *fldList);

/** Private: Set Class function pointers. */
static void ObitDConCleanVisClassInfoDefFn (gpointer inClass);

/** Private: reset sky model. */
static gboolean ResetSkyModel (ObitDConCleanVis *in, ObitErr *err);

/** Private: Get pixel array for a given field. */
static ObitFArray* GetFieldPixArray (ObitDConCleanVis *in, olong field, 
				     ObitErr *err);

/** Private: Low accuracy subtract CLEAN model. */
static void SubNewCCs (ObitDConCleanVis *in, olong *newCC, 
		       ObitFArray **pixarray, ObitErr *err);

/** Public: Pick next field(s) and get Residual image(s) */
static gboolean ObitDConCleanVisPickNext2D(ObitDConCleanVis *in, 
					   ObitErr *err);
static gboolean ObitDConCleanVisPickNext3D(ObitDConCleanVis *in, 
					   ObitErr *err);

/** Private: Set Beam patch, min. flux, decide on SDI CLEAN. */
static void ObitDConCleanVisDecide (ObitDConCleanVis* in, ObitErr *err);

/** Private: Determine cleanable fluxes for a field. */
static ofloat ObitDConCleanVisCleanable(ObitDConCleanVis *in, olong field, 
					ObitFArray *pixarray, ObitErr *err);

/** Private: Threaded Image subtractor */
static gpointer ThreadImSub (gpointer arg);

/** Private: Make Threaded Image subtraction args */
static olong MakeImSubFuncArgs (ObitThread *thread,
				ObitErr *err, ImSubFuncArg ***ThreadArgs);

/** Private: Delete Threaded Image subtraction args */
static void KillImSubFuncArgs (olong nargs, ImSubFuncArg **ThreadArgs);
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
  out->quality  = ObitMemFree (out->quality);
  out->cleanable = ObitMemFree (out->cleanable);

  /* In with the new */
  nfield       = in->nfield;
  out->quality   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  for (i=0; i<nfield; i++) {
    out->quality[i]    = in->quality[i];
    out->cleanable[i]  = in->cleanable[i];
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
  out->cleanable = ObitMemFree (out->cleanable);

  /* In with the new */
  nfield       = in->nfield;
  out->quality   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  for (i=0; i<nfield; i++) {
    out->quality[i]   = in->quality[i];
    out->cleanable[i] = in->cleanable[i];
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
  out->gain        = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->quality     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
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
  out->gain        = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->quality     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
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
 * \li "autoCen" OBIT_float scalar = Threshold level for autocenter
 *               If peak exceeds this value the minFlux reset to 0.1*autoCen
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanVisDeconvolve (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanVis *in;
  ObitFArray **pixarray=NULL;
  gboolean done, fin=TRUE, quit, doSub, bail, doMore, moreClean, notDone;
  gboolean redo;
  olong jtemp, i, *startCC=NULL, *newCC=NULL, count, ifld;
  olong redoCnt=0;
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
  in->doBeam = FALSE;                               /* Shouldn't need again */

  /* Restart CLEAN here if more CLEANable found at end */
 doRedo:
  in->Pixels->complCode = OBIT_CompReasonUnknown;   /* Clean not yet started */
  /* Create local arrays */
  startCC = g_malloc0(in->nfield*sizeof(olong));
  newCC   = g_malloc0(in->nfield*sizeof(olong));

  /* Loop until Deconvolution done */
  done = in->niter<=0; /* Cleaning requested? */
  while (!done) {
    /* Decide which fields to do next, tells if finished CLEAN 
       depends on type of imaging */
    if (in->mosaic->images[0]->myDesc->do3D)
      done = ObitDConCleanVisPickNext3D(in, err);
    else
      done = ObitDConCleanVisPickNext2D(in, err);
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
			      in->currentFields[0], err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      if (quit) {done=TRUE; break;}
    }

    /* If using autoWindow, iterate until clean not limited by autoWindow flux */
    pixarray   = g_malloc0(in->numCurrentField*sizeof(ObitFArray*));
    for (ifld=0; ifld<in->numCurrentField; ifld++) {
      if (in->currentFields[ifld]<=0) break;  /* List terminated? */
      startCC[ifld]   = in->Pixels->iterField[in->currentFields[ifld]-1];  /* no. at start */
      /* Pixel array for field */
      pixarray[ifld]  = GetFieldPixArray (in, in->currentFields[ifld], err);  
    }

    /* Get image/beam statistics needed for this cycle */
    notDone = inClass->ObitDConCleanPixelStats((ObitDConClean*)in, pixarray, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (!notDone) {
      /* Nothing to do in current list */
      for (ifld=0; ifld<in->numCurrentField; ifld++) {
	if (in->currentFields[ifld] <= 0) break;  /* List terminates? */
	in->quality[in->currentFields[ifld]-1]   = 0.0;
	in->cleanable[in->currentFields[ifld]-1] = 
	  -in->cleanable[in->currentFields[ifld]-1];
	
      }
    }
    
    /*fprintf (stderr,"DEBUG doMore %d done %d\n", doMore, done);*/

    count = 0;
    moreClean = TRUE;
    /* Middle CLEAN loop */
    while (moreClean) {
      /* Check if new box added to window after first pass */
      if ((count>0) && in->autoWindow) 
	doMore = inClass->ObitDConCleanAutoWindow ((ObitDConClean*)in, in->currentFields, 
						   pixarray, err);
      else {
	doMore = in->autoWindow;  /* pretend on the first pass */
      }
      if (err->error) Obit_traceback_msg (err, routine, in->name);
     
      /* Are we done */
      if (!doMore) break;

      /* Number of components before CLEAN */
      for (ifld=0; ifld<in->numCurrentField; ifld++) {
	newCC[ifld] = in->Pixels->iterField[in->currentFields[ifld]-1]; 
      }

      /* Pick components for this major cycle */
      if ((count==0) || doMore) 
	fin   = inClass->ObitDConCleanSelect((ObitDConClean*)in, pixarray, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
      else ObitErrClear(err);

      /* Did it just  stop because of autoWindow? */
      moreClean = (in->Pixels->complCode==OBIT_CompReasonAutoWin);
      /* Don't do this forever */
      if ((count>0) && !doMore) moreClean = FALSE;
      if (count>10) moreClean = FALSE;
      /*fprintf (stderr,"DEBUG doMore %d count %d moreClean %d Peak/RMS %g\n", 
	doMore, count, moreClean, in->imgPeakRMS[31]);*/
      if (!moreClean) break;

      /* Subtract these CCs from images in pixarray if possible */
      SubNewCCs(in, newCC, pixarray, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      count++;
    } /* End middle CLEAN loop */
    
    /* Update quality list for new max value on fields just CLEANed */
    for (ifld=0; ifld<in->numCurrentField; ifld++) {
      if (in->currentFields[ifld] <= 0) break;  /* List terminates? */
      if (fin) in->maxAbsRes[in->currentFields[ifld]-1] = 0.0;
      else in->maxAbsRes[in->currentFields[ifld]-1] = in->Pixels->maxResid;
      in->quality[in->currentFields[ifld]-1] = 
	ObitDConCleanVisQuality((ObitDConCleanVis*)in, in->currentFields[ifld], err);
      in->cleanable[in->currentFields[ifld]-1] = in->maxAbsRes[in->currentFields[ifld]-1];
    }

    /* Release working pixel arrays */
    if (pixarray) {
      for (ifld=0; ifld<in->numCurrentField; ifld++) {
	if (in->currentFields[ifld]<=0) break;  /* List terminated? */
	pixarray[ifld] = ObitFArrayUnref(pixarray[ifld]);
      }
      g_free(pixarray);  pixarray = NULL;
    }

    /* Clear BeamPatches array */
    if (in->BeamPatches) {
      for (i=0; i<in->mosaic->numberImages; i++) 
	if (in->BeamPatches[i]) ObitFArrayUnref(in->BeamPatches[i]);
      in->BeamPatches =  ObitMemFree (in->BeamPatches);
    }

    /* Subtract any components from visibility data */
    inClass->ObitDConCleanSub((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

  } /* end clean loop */
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Cleanup */
  if (startCC) g_free(startCC);
  if (newCC) g_free(newCC);

    /* Release working pixel arrays */
    if (pixarray) {
      for (ifld=0; ifld<in->numCurrentField; ifld++) {
	if (in->currentFields[ifld]<=0) break;  /* List terminated? */
	pixarray[ifld] = ObitFArrayUnref(pixarray[ifld]);
      }
      g_free(pixarray); pixarray = NULL;
    }

  /* Make final residuals */
  if ((!bail) && (in->niter>0)) {
    MakeAllResiduals (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
    /* Check if any fields now CLEANable and autoWin and more clean allowed */
    if ((in->autoWindow) && 
	(in->Pixels->currentIter<in->Pixels->niter)) {
      redo = FALSE;
      for (ifld=0; ifld<in->nfield; ifld++) {
	redo = redo || ((in->cleanable[ifld]>=in->minFlux[ifld])
			&& (in->imgPeakRMS[ifld]>=5.0));
      }
      /* Something left to CLEAN? */
      if (redo && (redoCnt<1)) {
	redoCnt++;
	Obit_log_error(err, OBIT_InfoErr,"Found more to CLEAN");
	goto doRedo;
      }
    } /* End check for more on autoWin */
  }
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

    /* If 2D imaging concatenate CC tables */
    if (!in->mosaic->images[0]->myDesc->do3D) 
      ObitImageMosaicCopyCC (in->mosaic, err);
    
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
  olong *itemp, jtemp, ifld, nfield=0;
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
    /* Loop over fields */
    for (ifld=0; ifld<in->numCurrentField; ifld++) {
      if (in->currentFields[ifld]<=0) break; /* List terminated? */
      /* compress CCs */
      if (in->currentFields[ifld]<=0) {
	ObitSkyModelCompressCC (in->skyModel, err);  /* All */
	nfield = in->mosaic->numberImages;
	for (i=0; i<nfield; i++) in->skyModel->startComp[i] = 1;
	for (i=0; i<nfield; i++) in->skyModel->endComp[i]   = 0;
      } else { /* only in->currentField */
	
	field = in->currentFields[ifld]-1;
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
    } /* end loop over fields */
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
 * 2D version - multiple parallel facets possible
 * Pick next fields to clean and create residual image in mosaic member
 * The selected fields are in->currentFields
 * Looks for "autoCenFlux" value in the first autoWindow image.
 * Adopted from the AIPSish CLBSTF (QOOP:QCLEAN.FOR)
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 * \return TRUE iff reached minimum flux density or max. number  comp.
 *         or no fields with SNR>5
 */
static gboolean ObitDConCleanVisPickNext2D(ObitDConCleanVis *in, ObitErr *err)
{
  olong i, best, lastBest=-1, loopCheck, indx, NumPar;
  olong *fldList=NULL;
  gboolean *fresh, doBeam=FALSE, done=TRUE, found;
  ofloat sumwts, autoCenFlux=0.0;
  ObitImage *theBeam=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  const ObitUVImagerClassInfo *imagerClass;
  gchar *routine = "ObitDConCleanVisPickNext2D";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanVisIsA(in));

  /* Check if reached max number of components and some done */
  if ((in->Pixels->currentIter >= in->Pixels->niter) && 
      (in->Pixels->currentIter>0)) return done;

  /* If only one field - it's the one */
  if (in->nfield==1) {
    in->currentFields[0] = 1;
    in->currentFields[1] = 0;
    in->numCurrentField  = 1;
    MakeResiduals(in, in->currentFields, FALSE, err);
    done = (in->cleanable[0] <= in->minFlux[0]) || 
      ((in->Pixels->maxResid <= in->minFlux[0]) && (in->Pixels->currentIter>1)) || 
      (in->imgPeakRMS[0]<MAX (4.0,(0.1*in->beamPeakRMS[0])));
    return done;
  }

  /* Check if there are autoCenter windows, if so get level */
  for (i=0; i<in->mosaic->numberImages; i++) {
    if (in->mosaic->isAuto[i]>0) {
      ObitInfoListGetTest(in->mosaic->images[i]->info, "autoCenFlux", &type, 
			  dim, &autoCenFlux);
      break;
    }
  }

  fresh   = ObitMemAlloc0(in->nfield*sizeof(gboolean));
  fldList = ObitMemAlloc0((in->nfield+3)*sizeof(olong));

  /* How many images in parallel? */
  imagerClass = (ObitUVImagerClassInfo*)in->imager->ClassInfo;
  NumPar = imagerClass->ObitUVImagerGetNumPar (in->imager, err);
  NumPar = MIN (NumPar, in->nfield);  /* No more than what's there */
  
  /* First time? */
  if (in->Pixels->currentIter > 0) {
    
    /* Already have done some CLEANing - need to remake residuals */
    for (i=0; i<in->nfield; i++) fresh[i] = FALSE;
    
    /* Make sure all fields initialized */
    indx = 0;
    for (i=0; i<in->nfield; i++) {
      if (in->maxAbsRes[i] < 0.0) {
	doBeam = in->doBeam;
	/* Make sure beam is made - is SUMWTS present? */
	theBeam = (ObitImage*)in->mosaic->images[i]->myBeam;
	if (!ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts))
	  doBeam = TRUE;
	/* Make residual image - get statistics */
	fldList[indx++] = i+1;
	fresh[i] = TRUE;
      }
    } /* end loop initializing fields */
    /* Make residual images */
    if (fldList[0]>0) MakeResiduals(in, fldList, doBeam, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);

    /* Check if reached max number of components */
    if (in->Pixels->currentIter >= in->Pixels->niter) return done;
    
    /* Ignore fields already known to be finished or low SNR, see if all done */
    done = TRUE;
    for (i=0; i<in->nfield; i++) {
      if ((in->maxAbsRes[i] <= in->minFlux[i]) || (in->imgPeakRMS[i]<5.0)) {
	in->quality[i]   = 0.0;
	in->cleanable[i] = -in->cleanable[i];
      }
      else done = FALSE;
    }
    
    /* anything left? */
    if (done) {ObitMemFree(fresh); ObitMemFree(fldList); return done;} 
    
  } else {
    /* First pass - all images should be OK */  
    for (i=0; i<in->nfield; i++) fresh[i] = TRUE;
  }
  
  /* Shouldn't need to make beams again */
  in->doBeam = FALSE;
    
  /* Loop remaking blocks of images until something suitable to CLEAN */
  done = FALSE;
  loopCheck = 0;
  while (!done) {
    
    /* Get ordered list */
    OrderImage (in, fresh, autoCenFlux, fldList); 
    
    /* No more than NumPar */
    fldList[NumPar] = 0;
    
    /* Count */
    in->numCurrentField = 0;
    for (i=0; i<in->nfield; i++) {
      if (fldList[i]>0) in->numCurrentField++;  /* Count */
      else break;
    }

    /* Make residual images if needed */
    if (fldList[0]>0) MakeResiduals(in, fldList, FALSE, err);
    /* If this is not the first call we should need to make residuals */
    else if (in->Pixels->currentIter>0) break;  /* Apparently nothing to do for now */
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    
    /* Which ones remade? */
    for (i=0; i<in->nfield; i++) {
      if (fldList[i]==0) break;
      /* Ignore fields with max residual less than min. */
      if (in->maxAbsRes[fldList[i]-1] < in->minFlux[fldList[i]-1]) {
	in->quality[fldList[i]-1]   = 0.0;
	in->cleanable[fldList[i]-1] = -in->cleanable[fldList[i]-1];
      }
      fresh[fldList[i]-1] = TRUE;
    }
    
    /* See if all done */
    done = TRUE;
    for (i=0; i<in->nfield; i++) {
      if (in->maxAbsRes[i] > in->minFlux[i]) done = FALSE;
    }
    if (done) {ObitMemFree(fresh); ObitMemFree(fldList); return done;} 
    
    /* Which fields ready to CLEAN? 
     if the best one is an autoCenter field only it is done */
    OrderClean (in, fresh, autoCenFlux, fldList);
    
    /* Find current estimates best field */
    WhosBest2D (in, autoCenFlux, &best);

    /* No best means we're done */
    if (best<0) return TRUE;

    /* Don't loop if best not changing - it may not be cleanable */
    if (best==lastBest) break;
    lastBest = best;
    
    /* make sure current best in fldList */
    found = FALSE;
    for (i=0; i<in->nfield; i++) {
      if (fldList[i]==0) break;
      if ((best+1)==fldList[i]) {found = TRUE; break;}
    }

    if (found && (fldList[0] > 0)) break;  /* At least one */
    /* Are these CLEANable? */
    if (in->quality[best]==0.0) break;
    
    /* Don't loop forever */
    loopCheck++;
    if (loopCheck>2*in->nfield) break;
    
    /* Message */
    if (fldList[0]>0)
      Obit_log_error(err, OBIT_InfoWarn,
		     "%s:  There may be something better - %f<%f",
		     routine, in->quality[fldList[0]-1], in->quality[best]);
    else
      Obit_log_error(err, OBIT_InfoWarn,
		     "%s: Nothing - try again ",
		     routine);

  } /* end loop reimaging */

  /* Get highest abs value in any CLEAN window */
  for (i=0; i<in->nfield; i++) {
    in->peakFlux = MAX (in->peakFlux, in->maxAbsRes[i]);
  }

  /* We're done if best field has maxAbsRes<minFlux */
  done = in->maxAbsRes[fldList[0]-1] <= in->minFlux[fldList[0]-1];

  /* Save list to clean, count */
  in->numCurrentField = 0;
  for (i=0; i<in->nfield; i++) {
    if (fldList[i]>0) in->numCurrentField++;  /* Count */
    else break;
  }
  for (i=0; i<in->nfield; i++) {
    in->currentFields[i] = fldList[i];
  }
  in->currentFields[i] = 0;
 
  /* cleanup */
  fresh   = ObitMemFree(fresh);
  fldList = ObitMemFree(fldList);

  return done;
} /* end ObitDConCleanVisPickNext2D */

/**
 * 3D version
 * Pick next field to clean and create residual image in mosaic member
 * The selected field is in->currentFields[0]
 * Adopted from the AIPSish CLBSTF (QOOP:QCLEAN.FOR)
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 * \return TRUE iff reached minimum flux density or max. number  comp.
 *         or no fields with SNR>5
 */
static gboolean ObitDConCleanVisPickNext3D(ObitDConCleanVis *in, ObitErr *err)
{
  olong i, best, second, loopCheck, indx, NumPar;
  olong *fldList;
  gboolean *fresh, doBeam=FALSE, done=TRUE;
  ofloat sumwts;
  ObitImage *theBeam=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  const ObitUVImagerClassInfo *imagerClass;
  gchar *routine = "ObitDConCleanVisPickNext3D";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanVisIsA(in));

  /* Check if reached max number of components and some done */
  if ((in->Pixels->currentIter >= in->Pixels->niter) && 
      (in->Pixels->currentIter>0)) return done;

  /* If only one field - it's the one */
  if (in->nfield==1) {
    in->currentFields[0] = 1;
    in->currentFields[1] = 0;
    in->numCurrentField  = 1;
    MakeResiduals(in, in->currentFields, FALSE, err);
    done = (in->cleanable[0] <= in->minFlux[0]) || 
      ((in->Pixels->maxResid <= in->minFlux[0]) && (in->Pixels->currentIter>1)) || 
      (in->imgPeakRMS[0]<MAX (4.0,(0.1*in->beamPeakRMS[0])));
    return done;
  }

  fresh   = ObitMemAlloc0(in->nfield*sizeof(gboolean));
  fldList = ObitMemAlloc0((in->nfield+3)*sizeof(olong));

  /* How many images in parallel? */
  imagerClass = (ObitUVImagerClassInfo*)in->imager->ClassInfo;
  NumPar = imagerClass->ObitUVImagerGetNumPar (in->imager, err);
  NumPar = MIN (NumPar, in->nfield);  /* No more than what's there */
  /* No point in doing too many */
  NumPar = MIN (NumPar, 2); 

  /* First time? */
  if (in->Pixels->currentIter > 0) {
    for (i=0; i<in->nfield; i++) fresh[i] = FALSE;
    
    /* Make sure all fields initialized */
    indx = 0;
    for (i=0; i<in->nfield; i++) {
      if (in->maxAbsRes[i] < 0.0) {
	doBeam = in->doBeam;
	/* Make sure beam is made - is SUMWTS present? */
	theBeam = (ObitImage*)in->mosaic->images[i]->myBeam;
	if (!ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts))
	  doBeam = TRUE;
	/* Make residual image - get statistics */
	fldList[indx++] = i+1;
	fresh[i] = TRUE;
      }
    } /* end loop initializing fields */
    
    /* Make images */
    if (fldList[0]>0) MakeResiduals (in, fldList, doBeam, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);

    /* Check if reached max number of components */
    if (in->Pixels->currentIter >= in->Pixels->niter) return done;
    
    /* Ignore fields already known to be finished or low SNR, see if all done */
    done = TRUE;
    for (i=0; i<in->nfield; i++) {
      if ((in->maxAbsRes[i] <= in->minFlux[i]) || (in->imgPeakRMS[i]<5.0)) {
	in->quality[i]   = 0.0;
	in->cleanable[i] = -in->cleanable[i];
      }
      else done = FALSE;
    }
    /* anything left? */
    if (done) {ObitMemFree(fresh); ObitMemFree(fldList); return done;} 
    
    /* Find current estimates best and second best field field */
    WhosBest (in, &best, &second);

    /* No best means we're done */
    if (best<0) return TRUE;
    
    /* Verify that it's still best which may require iteration */
    done = FALSE;
    loopCheck = 0;
    while (!done) {
      /* If not just created recalculate */
      if ((best>=0) && (!fresh[best])) {
	
	/* Make best two residual image - get statistics */
	fldList[0] = best+1;
	if ((second>=0) && (!fresh[second]) && (NumPar>1)) fldList[1] = second+1;
	else fldList[1] = 0;
	fldList[2] = 0;
	MakeResiduals(in, fldList, FALSE, err);
	if (err->error) Obit_traceback_val (err, routine, in->name, done);
	fresh[best]   = TRUE;
	if ((second>=0) && (fldList[1]>0)) fresh[second] = TRUE;
	/* Ignore fields with max residual less than min. */
	if (in->maxAbsRes[best] < in->minFlux[best]) {
	  in->quality[best]   = 0.0;
	  in->cleanable[best] = -in->cleanable[best];
	}
      } /* end remake */
      
	/* Don't loop forever */
      loopCheck++;
      if (loopCheck>2*in->nfield) break;
	
      if (in->nfield==1) break; /* If only one, the decision is easy */
      
      /* Still best? */
      if (best<0) break;
      if ((second<0) || (in->quality[best] >= in->quality[second])) break;
      
      /* try again */
      if ((in->quality[best]>0.0) && (in->prtLv>1)) 
	Obit_log_error(err, OBIT_InfoWarn, 
		       "%s: field %d (%g) not as good as second %g",
		       routine, best+1, in->quality[best],  in->quality[second]);
      
      /* Find current estimates best and second best field field */
      WhosBest (in, &best, &second);
      
    } /* end loop finding best field */
    
  } else {
    /* First pass - all images should be OK */  
    for (i=0; i<in->nfield; i++) fresh[i] = TRUE;
    WhosBest (in, &best, &second);
  }

  /* Shouldn't need to make beams again */
  in->doBeam = FALSE;
    
  /* publish decision */
  in->currentFields[0] = best + 1;
  in->currentFields[1] = 0;
  in->numCurrentField  = 1;

  /* Diagnostic */
  if (in->prtLv>2) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "PickNext3D: best %d quality %g cleanable %g maxRes %g SNR %g",
		   best+1, in->quality[best],  in->cleanable[best], 
		   in->maxAbsRes[best], in->imgPeakRMS[best]);
  }
  
  /* highest abs value ever found in a CLEAN window */
  if (best>=0) in->peakFlux = MAX (in->peakFlux, in->maxAbsRes[best]);

  /* We're done if best field has maxAbsRes<minFlux */
  if (best>=0) done = in->maxAbsRes[best] <= in->minFlux[best];

  /* cleanup */
  fresh   = ObitMemFree(fresh);
  fldList = ObitMemFree(fldList);

  return done;
} /* end ObitDConCleanVisPickNext3D */

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
 * A new (128x128) field is added centered on the offending source and a negative
 * clean window added to the position of the source in its original window.
 * Avoid duplicates of the same source and ensure that all occurances of this 
 * source in any exant field has a negative clean window added.
 * Multiple centering sources per facet are allowed
 * A boolean entry "autoCenField" with value True is added to the info member of any
 * image members added to the mosaic member.
 * A float entry  "autoCenField" with the value of restartFlux is added 
 * to the info member of any autoCen fields.
 * If 2D imaging is being used, for each autoCenter field, a shifted
 * field is also added.
 * Routine originally translated from the AIPSish VLAUTIL.FOR/VLCCIN  
 * \param in  DConCleanVis to check
 *  Values on info member:
 * \li "restartFlux" OBIT_float scalar = Minimum brightness flux for autoCenter
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
  odouble pos[2], pos2[2], RAPnt, DecPnt;
  gboolean done, outside=FALSE, facetDone, clear, Tr=TRUE, Fl=FALSE;
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

  /* Flux level to autoCenter */
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
	/* Get position of center */
	pixel[0] = 1.0 + imDesc2->inaxes[0]*0.5;
	pixel[1] = 1.0 + imDesc2->inaxes[1]*0.5;
	ObitImageDescGetPos(imDesc2, pixel, pos2, err);
	deltax = (pos2[0]-pos[0]) / imDesc->cdelt[0];
	deltay = (pos2[1]-pos[1]) / imDesc->cdelt[1];
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
	ObitInfoListAlwaysPut(mosaic->images[ifield]->info, 
			      "autoCenFlux", OBIT_float, dim, &freset);
	mosaic->isAuto[ifield] =  mosaic->numberImages+9999999; /* fill in later */
	
	ObitImageDescGetPoint (imDesc, &RAPnt, &DecPnt);
	ObitSkyGeomShiftXY (RAPnt, DecPnt, ObitImageDescRotate(imDesc),
			    pos[0], pos[1], &RAShift, &DecShift);
	/* Update mosaic and shift in other places */
	mosaic->RAShift[ifield]  = RAShift;
	mosaic->DecShift[ifield] = DecShift;
	if (imDesc->do3D) {
	  /* 3D shift position */
	  imDesc->crval[imDesc->jlocr] = pos[0];
	  imDesc->crval[imDesc->jlocd] = pos[1];
	} else {
	  /* Not 3D - shift reference pixel */
	  imDesc->crpix[imDesc->jlocr] -=  xoff/imDesc->cdelt[imDesc->jlocr];
	  imDesc->crpix[imDesc->jlocd] -=  yoff/imDesc->cdelt[imDesc->jlocd];
	  imDesc->xPxOff -= xoff/imDesc->cdelt[imDesc->jlocr];
	  imDesc->yPxOff -= yoff/imDesc->cdelt[imDesc->jlocd];
	}
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
	nx = ny = ObitFFTSuggestSize (128); nplane = 1;
	ObitImageMosaicAddField (in->mosaic, uvdata, nx, ny, nplane, 
				 RAShift, DecShift, TRUE, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

	/* Mark as an autoCenter Image */
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(mosaic->images[mosaic->numberImages-1]->info, 
			      "autoCenField", OBIT_bool, dim, &Tr);
	ObitInfoListAlwaysPut(mosaic->images[mosaic->numberImages-1]->info, 
			      "autoCenFlux", OBIT_float, dim, &freset);
	
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
	xcen = mosaic->images[newField-1]->myDesc->crpix[0] - 
	  mosaic->images[newField-1]->myDesc->xPxOff;
	ycen = mosaic->images[newField-1]->myDesc->crpix[1] - 
	  mosaic->images[newField-1]->myDesc->yPxOff;
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
	
	/* Allow shift wrt 2D grid */
	if (!mosaic->images[mosaic->numberImages-1]->myDesc->do3D) {
	  dim[0] = dim[1] = 1;
	  ObitInfoListAlwaysPut(mosaic->images[mosaic->numberImages-1]->myDesc->info, 
				"doGrid", OBIT_bool, dim, &Fl);
	}
	
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

  /* If 2D imaging, make any shifted images needed */
  if (!in->mosaic->images[0]->myDesc->do3D) {
   nfield = mosaic->numberImages;
   for (ifield=0; ifield<nfield; ifield++) {
      if (in->mosaic->isAuto[ifield] > nfield) {
	/* Duplicate of autoCenter field but with crpix at closest integer */
	nx = mosaic->images[ifield]->myDesc->inaxes[0];
	ny = mosaic->images[ifield]->myDesc->inaxes[1];
	nplane   = 1;
	RAShift  = mosaic->RAShift[ifield];
	DecShift = mosaic->DecShift[ifield];
	ObitImageMosaicAddField (in->mosaic, uvdata, nx, ny, nplane, 
				 RAShift, DecShift, FALSE, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

	/* Add field to window */
	inaxes[0] = nx; inaxes[1] = ny;
	newField = ObitDConCleanWindowAddField (in->window, inaxes, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

	/* Add inner and outer windows */
	xcen = mosaic->images[newField-1]->myDesc->crpix[0] - 
	  mosaic->images[newField-1]->myDesc->xPxOff;
	ycen = mosaic->images[newField-1]->myDesc->crpix[1] - 
	  mosaic->images[newField-1]->myDesc->yPxOff;
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
	
	/* Indicate association */
	in->mosaic->isAuto[ifield]  = mosaic->numberImages;
	in->mosaic->isShift[mosaic->numberImages-1] = ifield+1;
      }
    } /* end loop over fields */   
  } /* end of adding shifted 2D images */

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
  olong i, oldField, *itemp;
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
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean cleanable");
  for (i=0; i<oldField; i++) ftemp[i] = in->cleanable[i]; ftemp[i] = 0.0; 
  in->cleanable = ObitMemFree(in->cleanable);
  in->cleanable = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean max res");
  for (i=0; i<oldField; i++) ftemp[i] = in->maxAbsRes[i]; ftemp[i] = 0.0; 
  in->maxAbsRes = ObitMemFree(in->maxAbsRes);
  in->maxAbsRes = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean avg res");
  for (i=0; i<oldField; i++) ftemp[i] = in->avgRes[i]; ftemp[i] = 0.0; 
  in->avgRes = ObitMemFree(in->avgRes);
  in->avgRes = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Image Peak/RMS");
  for (i=0; i<oldField; i++) ftemp[i] = in->imgPeakRMS[i]; ftemp[i] = 0.0; 
  in->imgPeakRMS = ObitMemFree(in->imgPeakRMS);
  in->imgPeakRMS = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Beam Peak/RMS");
  for (i=0; i<oldField; i++) ftemp[i] = in->beamPeakRMS[i]; ftemp[i] = 0.0; 
  in->beamPeakRMS = ObitMemFree(in->beamPeakRMS);
  in->beamPeakRMS = ftemp;
  itemp = ObitMemAlloc0Name((newField+3)*sizeof(olong),"currentFields");
  for (i=0; i<oldField; i++) itemp[i] = in->currentFields[i]; itemp[i] = 0; 
  in->currentFields = ObitMemFree(in->currentFields);
  in->currentFields = itemp;

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
	if (imDesc->do3D) {
	  /* 3D shift position */
	  imDesc->crval[imDesc->jlocr] = pos[0];
	  imDesc->crval[imDesc->jlocd] = pos[1];
	} else {
	  /* Not 3D - shift reference pixel */
	  imDesc->crpix[imDesc->jlocr] -=  xoff/imDesc->cdelt[imDesc->jlocr];
	  imDesc->crpix[imDesc->jlocd] -=  yoff/imDesc->cdelt[imDesc->jlocd];
	  imDesc->xPxOff -= xoff/imDesc->cdelt[imDesc->jlocr];
	  imDesc->yPxOff -= yoff/imDesc->cdelt[imDesc->jlocd];
	}
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
 * Automatically set a box in window if needed
 * autoWindow feature will automatically set CLEAN windows inside 
 * a predefined outer window.  Each cycle the residuals inside the outer 
 * window are searched to the maximum value; if the peak is outside the 
 * inner window and > 4 sigma, a new round box  is added 
 * to the window.  Cleaning in each cycle will stop when the peak residual 
 * drops to the level of the highest value outside the CLEAN window.
 * For fields not currently being CLEANed, the minimum level is 0.7*peak-3 sigma.
 * Only the field with the highest quality is modified.
 * \param in       The object to restore
 * \param fields   Field numbers (1-rel) in ImageMosaic
 *                 zero terminated list, no longer than in->numCurrentField
 * \param pixarray If nonNULL, then use this array of pixels rather than in the ImageMosaic
 *                 Elements are Pixel arrays corresponfing to fields in fields
 * \param err      Obit error stack object.
 * \return TRUE if clean should be continued
 */
gboolean ObitDConCleanVisAutoWindow(ObitDConClean *inn, olong *fields, ObitFArray **pixarray, 
				    ObitErr *err)
{
  ObitDConCleanVis *in=(ObitDConCleanVis*)inn;
  gboolean newWin=FALSE, doMore=FALSE;
  ObitImage *image=NULL;
  ObitFArray *usePixels;
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  i, j, best, field, blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat PeakIn, PeakOut, RMS, bestPeak, oldAutoWinFlux, *FieldPeak=NULL;
  olong PeakInPos[2] = {0,0};
  gboolean doAbs, found;
  gchar *routine = "ObitDConCleanVisAutoWindow";

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanIsA(in));

  /* Set output to full image, plane at a time */
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];

  /* find best field */
  FieldPeak = g_malloc0(in->numCurrentField*sizeof(ofloat));
  for (field=0; field<in->numCurrentField; field++) FieldPeak[field]=0.0;

  best = 0;
  bestPeak = -1.0e20;
  for (field=0; field<in->numCurrentField; field++) {
    if (fields[field]==0) break;  /* List terminated? */

    /* Use passed pixel array or get image */
    if (pixarray==NULL) {
      image = in->mosaic->images[fields[field]-1];
      
      /* Set input to full image, plane at a time */
      dim[0] = IM_MAXDIM;
      ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
      ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc); 
      dim[0] = 1;
      ObitInfoListAlwaysPut (image->info, "IOBy", OBIT_long, dim, &IOsize);
      
      retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
      retCode = ObitImageRead (image, image->image->array, err);
      retCode = ObitImageClose (image, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
      usePixels = image->image;  /* Pointer to image buffer */
      
    } else { /* Use array passed */
      usePixels = ObitFArrayRef(pixarray[field]);
      image = in->mosaic->images[0];  /* Need one */
    }
    
    /* Allow negative for Stokes other than I */
    doAbs = fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1;
  
    /* Get statistics */
    ObitDConCleanWindowStats (in->window, fields[field], usePixels,
			      doAbs,
			      &PeakIn, &PeakInPos[0], &PeakOut, 
			      &RMS, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, newWin);

    /* Free Image array */
    if (usePixels==image->image) 
      image->image = ObitFArrayUnref(image->image);
    else
      usePixels = ObitFArrayUnref(usePixels);
    
    /* Save peak not in window */
    FieldPeak[field] = PeakOut;
    /* Diagnostic */
    if (in->prtLv>4) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "Clean field %d PeakIn %f PeakOut %f RMS %f",
		     fields[field], PeakIn, PeakOut, RMS);
    }

    /* Check for best */
    if (PeakIn>bestPeak) {
      best = field;
      bestPeak = PeakIn;
    }

    /* Reset peak/RMS */
    in->imgPeakRMS[fields[field]-1] = PeakIn/RMS;
  } /* end loop gathering statistics */

  /* Use best field */
  field = fields[best];

  /* Use passed pixel array or get image */
  if ((pixarray==NULL) || (pixarray[best]==NULL)) {   
    image = in->mosaic->images[field-1];

    /* Set input to full image, plane at a time */
    dim[0] = IM_MAXDIM;
    ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
    ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc); 
    dim[0] = 1;
    ObitInfoListAlwaysPut (image->info, "IOBy", OBIT_long, dim, &IOsize);
    
    retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
    
    retCode = ObitImageRead (image, image->image->array, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
    
    retCode = ObitImageClose (image, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
    usePixels = image->image;  /* Pointer to image buffer */

  } else { /* Use array passed */
    usePixels = ObitFArrayRef(pixarray[best]);
    image = in->mosaic->images[0];  /* Need one */
  }
    
  /* Allow negative for Stokes other than I */
  doAbs = fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1;
  
  /* Get field info - set new box if needed */
  newWin = ObitDConCleanWindowAutoWindow (in->window, field, usePixels,
					  doAbs,
					  &PeakIn, &PeakInPos[0], &PeakOut, 
					  &RMS, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, newWin);

  /* Free Image array? */
  if (usePixels==image->image) 
    image->image = ObitFArrayUnref(image->image);
  else
    usePixels = ObitFArrayUnref(usePixels);
  
  /* Determine max. cleanable pixel value not in a window */
  for (field=0; field<in->numCurrentField; field++) {
    if (fields[field]==0) break;  /* List terminated? */
    if (field!=best) PeakOut = MAX (PeakOut,  FieldPeak[field]);
  }
  g_free(FieldPeak);

  /* Set flux limit for next cycle */
  oldAutoWinFlux = in->autoWinFlux;
  if (oldAutoWinFlux<0.0) oldAutoWinFlux = 1.0e20;
  in->autoWinFlux = MAX (PeakOut-5.0*RMS, 0.5*RMS);

  /* Don't let the autoWinFlux go much below cleanable maxima in other fields */
  for (i=0; i<in->nfield; i++) {
    /* Is this one in fields? */
    found = FALSE;
    for (j=0; j<in->numCurrentField; j++) {
      if ((fields[j]-1)==i) found = TRUE;
      /* Look for autoCenter/shifted pairs */
      if (in->mosaic->isAuto[i]  == (fields[j])) found = TRUE;
      if (in->mosaic->isShift[fields[j]-1] == (i+1)) found = TRUE;
      if (found) break;
    }
    if (!found && (in->imgPeakRMS[i]>4.0) && in->quality[i]>0.0) {
      /* Only to 70% minus 3 sigma */
      in->autoWinFlux = MAX (in->autoWinFlux, 0.7*in->cleanable[i]-3.0*RMS);
      /* Diagnostic */
      if (in->prtLv>4) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Nonclean field %d cleanable %f Peak/RMS %g quality %g maxAbsRes %g",
		       i+1, in->cleanable[i], in->imgPeakRMS[i], in->quality[i],
		       in->maxAbsRes[i]);
      }
    }
  }
  /* Don't set too low */
  in->autoWinFlux = MAX (in->autoWinFlux, 0.5*RMS);

  /* Diagnostic */
  if (in->prtLv>4) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "autoWin Flux %g previous %g best field %d",
		   in->autoWinFlux, oldAutoWinFlux, fields[best]);
  }

  /* If autoWinFlux substantially lower and lower than bestPeak, continue */
  doMore = (in->autoWinFlux<0.9*oldAutoWinFlux) && (in->autoWinFlux<bestPeak);
  /* DEBUG  
  fprintf (stderr,"DEBUG autoWinFlux %f RMS %f PeakOut %f PeakIn %f\n",
	   in->autoWinFlux, RMS, PeakOut, PeakIn); */

  /* Use minFlux on best Pixels */
  in->minFlux[fields[best]-1] = 
    MIN (in->minFlux[fields[best]-1], in->Pixels->minFlux[fields[best]-1]);

  /* Don't clean too far into the noise */
  if ((in->Pixels->currentIter>1) && (in->Pixels->maxResid<0.5*RMS) && (fields[best]>0)) {
     in->minFlux[fields[best]-1] = 0.5*RMS;
     /* Tell user */
     if (in->prtLv>0) {
       Obit_log_error(err, OBIT_InfoWarn,"Cleaned into noise %g - resetting minFlux",
		      RMS);
     }
  }
  /*fprintf (stderr,"DEBUG newWin %d doMore %d\n",newWin, doMore);*/
  return newWin||doMore;
} /* end ObitDConCleanVisAutoWindow */

/**
 * Get image and beam statistics and prepare to deconvolve
 * If autoWindow option is selected then windows may be added 
 * \param in        The object to deconvolve
 * \param pixarray  If NonNULL use instead of the flux densities from the image file.
 *                  This is an array of ObitFarrays corresponding to fields in
                    in->currentFields
 * \param err       Obit error stack object.
 * \return TRUE if CLEAN should continue
 */
gboolean ObitDConCleanVisPixelStats(ObitDConClean *inn, ObitFArray **pixarray, 
				    ObitErr *err)
{
  ObitDConCleanVis *in=(ObitDConCleanVis*)inn;
  gboolean newWin = FALSE, isLast = FALSE, allDone=TRUE;
  olong field, ifld;
  ObitImage *Beam=NULL;
  const ObitDConCleanClassInfo *inClass;
  gchar *routine = "ObitDConCleanVisPixelStats";

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanIsA(in));

  /* Something in pixarray? */
  if ((pixarray) && (pixarray[0]==NULL)) return FALSE;

  inClass = (ObitDConCleanClassInfo*)in->ClassInfo; /* class structure */

  /* Adjust windows if autoWindow */
  if (in->autoWindow) 
    newWin = inClass->ObitDConCleanAutoWindow (inn, in->currentFields, 
					       pixarray, err);
  else {
    in->autoWinFlux = -1.0e20; 
    newWin = FALSE;
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, newWin);

  /* Did the previous CLEAN reach the CLEAN limit - if so don't check
     if all fields are finished; if there was something else to CLEAN
     a new Window will be added later */
  allDone = TRUE;
  
  /* Loop over current fields */
  for (ifld=0; ifld<=in->numCurrentField; ifld++) {
    
    field = in->currentFields[ifld];
    if (field<=0) break;  /* List terminates */

    /* Is this field finished? */  
    if ((in->Pixels->complCode!=OBIT_CompReasonNiter) &&
	(in->Pixels->complCode!=OBIT_CompReasonMinFlux)) {
      allDone = allDone && 
	((in->cleanable[field-1] <= in->minFlux[field-1]) || 
	 (in->imgPeakRMS[field-1]<4.0));
    }

    /* Beam image */
    Beam = (ObitImage*)(in->mosaic->images[field-1]->myBeam);
    
    /* Get Beam histogram from first field */
    if (ifld==0) {
      if (in->BeamHist->field != in->currentFields[0]) {
	ObitDConCleanBmHistUpdate(in->BeamHist, Beam, in->plane, err);
	in->BeamHist->field = in->currentFields[0];
	if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
      }
    }

    /* Get Pixel histogram */
    if (ifld==0) { /* init on first */
      /* Is this the last one? */
      isLast = (in->currentFields[ifld+1]==0) || (in->numCurrentField==1);
      ObitDConCleanPxHistUpdate (in->PixelHist, field, in->plane, in->mosaic, 
				 in->window, isLast, err);
    } else { /* Add data */
      /* Is this the last one? */
      isLast = (in->currentFields[ifld+1]==0) || (ifld >= (in->numCurrentField-1));
      ObitDConCleanPxHistAdd (in->PixelHist, field, in->plane, in->mosaic, 
			      in->window, isLast, err);
    }
    if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
    
    /* Decide beamPatchSize, minFluxLoad, SDI Clean */
    ObitDConCleanVisDecide (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
  }

  return newWin || (!allDone);
} /* end ObitDConCleanVisPixelStats */

/**
 * Determine Beam patch size and minimum flux density to consider
 * such that the minimum flux density is the maximum ignored sidelobe
 * of the brightest pixel value.
 * Output members are beamPatchSize, minFluxLoad and depend on the members
 * currentFields, window, BeamHist and PixelHist being up to date.
 * Member doSDI is set to TRUE iff SDI CLEAN is needed
 * If the histogram is too coarse (too many in top bin) then numberSkip
 * is set to decimate the residuals so that they will fit;
 * \param in   The object to deconvolve
 * \param err  Obit error stack object.
 */
static void ObitDConCleanVisDecide (ObitDConCleanVis* in, ObitErr *err)
{
  olong minPatch, maxPatch, Patch, i;
  ofloat minFlux=0.0, fract;
  gchar *routine = "ObitDConCleanVisDecide";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* initialize values */
  in->beamPatchSize = in->minPatchSize;
  in->minFluxLoad = 0.0;
  in->numberSkip = 0;
  minPatch = in->minPatchSize;
  maxPatch = in->BeamHist->ncell;

  /* Maximum beam patch size needed for Clean window */
  maxPatch = 
    MIN (ObitDConCleanWindowSize(in->window, in->currentFields[0], err), maxPatch);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  maxPatch = MAX (maxPatch, in->minPatchSize);

  /* Minimum should not exceed maximum */
  minPatch = MIN (minPatch, maxPatch);

  /* Find largest patch size/min. flux where the number of residuals fits into
     the allowed memory (maxPixel) */
  Patch = -1;
  for (i=maxPatch; i>=minPatch; i--) {
    /* Minimum flux density */
    minFlux = ObitDConCleanBmHistPeak (in->BeamHist, i, err) * 
      in->PixelHist->histMax;
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (ObitDConCleanPxHistNumber(in->PixelHist, minFlux, err) < 
	in->maxPixel) {  /* this fits */
      Patch = i;
      /* diagnostics */
      if (in->prtLv>2) {
	Obit_log_error(err, OBIT_InfoErr, "%s field %d min %f Patch %d bmPeak %f pxMax %f",
		       routine, in->currentFields[0], minFlux, Patch, 
		       ObitDConCleanBmHistPeak (in->BeamHist, i, err), in->PixelHist->histMax);
	Obit_log_error(err, OBIT_InfoErr, "   HistNumber %d maxPixel  %d",
		       ObitDConCleanPxHistNumber(in->PixelHist, minFlux, err), in->maxPixel);
      }
      in->minFluxLoad = minFlux;
      break;
    }
  }

  /* Find solution? */
  if (Patch>0) {
    in->beamPatchSize = Patch;
    goto doSDI;
  }

  /* Loop through histogram for lowest flux limit that fits */
  for (i=0; i<in->PixelHist->ncell-1;i++) {
    if (in->PixelHist->hist[i]<in->maxPixel) { /* go team! */
     in->beamPatchSize = minPatch;
     in->minFluxLoad   = in->PixelHist->histMin + ((ofloat)(i+1)) * 
       ((in->PixelHist->histMax-in->PixelHist->histMin) / in->PixelHist->ncell);
     goto doSDI;
    }
  }

  /* If you get here the pixel histogram has too many entries in the
     highest bin, this usually means that there are many nearly 
     indistinquishable pixels.  Resort to descimation */
  /* Bottom of top pixel histogram bin */
  minFlux = in->PixelHist->histMax - (in->PixelHist->histMax-in->PixelHist->histMin) / 
    in->PixelHist->ncell;
  in->minFluxLoad = minFlux;
  in->numberSkip = 1 + ObitDConCleanPxHistNumber(in->PixelHist, minFlux, err) / 
    MAX (1, in->maxPixel);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* See if SDI CLEAN called for */
 doSDI:
  if (in->SDIGain>0.0) {
    fract = (ofloat)in->PixelHist->hist[in->PixelHist->ncell/2] / MAX(1, in->PixelHist->hist[0]);
    if (in->prtLv>1)
      Obit_log_error(err, OBIT_InfoErr,"SDI fraction %f", fract);
      
    in->doSDI = fract>in->SDIGain;
    if (in->doSDI) {
      /* minFlux = in->PixelHist->histMax * MAX (0.3333, in->ccfLim) */
      in->minFluxLoad = MAX (in->autoWinFlux, minFlux);
      in->minFluxLoad = MAX (in->minFluxLoad, in->PixelHist->histMax * MAX (0.3333, in->ccfLim));
      in->minFluxLoad = MAX (in->minFluxLoad, in->minFlux[in->currentFields[0]]);
      in->numberSkip  = 0;
      /* Find min. beam patch with exterior sidelobe <0.1 */
      i = 0;
      minFlux = 1.0;
      while ((i<in->BeamHist->ncell) && (minFlux>0.1)) {
	minFlux = ObitDConCleanBmHistPeak (in->BeamHist, i, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
	i++;
      }
      in->beamPatchSize = i;
    }
  } else {
    in->doSDI = FALSE;
    /* Impose min fraction of initial residual */
    in->minFluxLoad = MAX (in->minFluxLoad, in->PixelHist->histMax*in->ccfLim);
  }

  /* Give warning if skipping */
  if ((in->numberSkip>=1) && (in->prtLv>0)) 
    Obit_log_error(err, OBIT_InfoWarn,"%s: Too many residuals, taking 1 of every  %d",
		   routine, in->numberSkip+1);

} /* end  ObitDConCleanVisDecide */

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
  theClass->ObitDConCleanVisQuality   = (ObitDConCleanVisQualityFP)ObitDConCleanVisQuality;
  theClass->ObitDConCleanVisReimage   = (ObitDConCleanVisReimageFP)ObitDConCleanVisReimage;
  theClass->ObitDConCleanVisAddField  = (ObitDConCleanVisAddFieldFP)ObitDConCleanVisAddField;
  theClass->ObitDConCleanVisRecenter  = (ObitDConCleanVisRecenterFP)ObitDConCleanVisRecenter;
  theClass->ObitDConCleanVisFilter    = (ObitDConCleanVisFilterFP)ObitDConCleanVisFilter;
  theClass->ObitDConCleanAutoWindow = 
    (ObitDConCleanAutoWindowFP)ObitDConCleanVisAutoWindow;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanVisPixelStats;

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
  in->cleanable = NULL;
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
  in->cleanable = ObitMemFree(in->cleanable);
  in->display   = ObitDisplayUnref(in->display);
  in->SDIdata   = ObitUVUnref(in->SDIdata);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanVisClear */

/**
 * Make selected residual images and get statistics
 * \param in     The Clean object
 * \param fields zero terminated list of field numbers to image
 * \param doBeam If TRUE also make beam
 * \param err    Obit error stack object.
 */
static void  MakeResiduals (ObitDConCleanVis *in, olong *fields, 
			    gboolean doBeam, ObitErr *err)
{
  gboolean doWeight, doFlatten, found;
  const ObitDConCleanVisClassInfo *inClass;
  ObitUVImagerClassInfo *imgClass = 
    (ObitUVImagerClassInfo*)in->imager->ClassInfo;
  olong i, ifld, jfld, field=1;
  gchar *routine = "MakeResidual";

  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */

  /* Are residuals fresh? */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Make residual images */
  imgClass->ObitUVImagerImage (in->imager, fields, doWeight, doBeam, doFlatten, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Statistics */
  for (ifld=0; ifld<in->numCurrentField; ifld++) {
    field = fields[ifld];
    if (field<=0) break;
    /* Get statistics  for image */
    inClass->ObitDConCleanImageStats ((ObitDConClean*)in, field, FALSE, err);
    /* Need Beam statistics? */
    if (doBeam)
      inClass->ObitDConCleanImageStats ((ObitDConClean*)in, field, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Quality measure */
    in->quality[field-1] = ObitDConCleanVisQuality(in, field, err);

    /* Max cleanable flux */
    in->cleanable[field-1] = ObitDConCleanVisCleanable(in, field, NULL, err);
  } /* end statistics loop */
  
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* For any shifted fields get statistics */
  for (ifld=0; ifld<in->mosaic->numberImages; ifld++) {
    if (in->mosaic->isShift[ifld] > 0) {

      /* Has this field reached the min flux? maxAbsRes<minFlux */
      if (in->maxAbsRes[ifld]<=in->minFlux[ifld]) continue;
      jfld = in->mosaic->isShift[ifld]-1;
      /* If the corresponding autoCenter field just remade (as well as the shifted
	 version - determine the statistics from the image */
      found = FALSE;
      for (i=0; i<in->numCurrentField; i++) {
	if (fields[i]==0) break;
	if ((jfld+1)==fields[i]) {found = TRUE; break;}
      }

      if (found) {
	/* Get statistics  for image */
	inClass->ObitDConCleanImageStats ((ObitDConClean*)in, ifld+1, FALSE, err);
	/* Quality measure */
	in->quality[ifld] = ObitDConCleanVisQuality((ObitDConCleanVis*)in, ifld+1, err);
	/* Max cleanable flux */
	in->cleanable[field-1] = ObitDConCleanVisCleanable(in, field, NULL, err);
	in->beamPeakRMS[ifld] = in->beamPeakRMS[jfld];
      } else { /* Not remade - Use residual */
	in->maxAbsRes[ifld] = in->maxAbsRes[jfld];
	in->quality[jfld]   = in->maxAbsRes[jfld];  /* Reset Quality to residual */
	in->quality[ifld]   = in->quality[jfld];
	in->cleanable[ifld]   = MIN(in->cleanable[jfld], in->maxAbsRes[jfld]);
	in->beamPeakRMS[ifld] = in->beamPeakRMS[jfld];
      }
    }
  }

  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end MakeResiduals */

/**
 * Make all residual images and get statistics
 * \param in     The Clean object
 * \param err    Obit error stack object.
 */
static void  MakeAllResiduals (ObitDConCleanVis *in, ObitErr *err)
{
  gboolean doBeam = in->doBeam;
  gboolean doWeight, doFlatten;
  olong i, ifld, jfld, fields[2]={0,0};
  ObitUVImagerClassInfo *imgClass = 
    (ObitUVImagerClassInfo*)in->imager->ClassInfo;
  const ObitDConCleanVisClassInfo *inClass;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "MakeAllResiduals";
  
  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */

  /* Turn off things not needed */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Copy prtLv to in->mosaic->info */
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (in->mosaic->info, "prtLv", OBIT_long, dim, &in->prtLv);

  /* Parallel Image images without needing beam */
  imgClass->ObitUVImagerImage (in->imager, fields,  doWeight, doBeam, doFlatten, err);
  
  /* Loop over fields getting statistics for Image and Beam */
  for (i=0; i<in->nfield; i++) {
    inClass->ObitDConCleanImageStats ((ObitDConClean*)in, i+1, FALSE, err);
    if (doBeam)
      inClass->ObitDConCleanImageStats ((ObitDConClean*)in, i+1, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Quality measure */
    in->quality[i] = ObitDConCleanVisQuality(in, i+1, err);
    /* Max cleanable flux */
    in->cleanable[i] = ObitDConCleanVisCleanable(in, i+1, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end loop over field */

  /* For any shifted fields get statistics */
  for (ifld=0; ifld<in->mosaic->numberImages; ifld++) {
    if (in->mosaic->isShift[ifld] > 0) {
      jfld = in->mosaic->isShift[ifld]-1;
      /* Get statistics  for image */
      inClass->ObitDConCleanImageStats ((ObitDConClean*)in, ifld+1, FALSE, err);
      /* Quality measure */
      in->quality[ifld] = ObitDConCleanVisQuality(in, ifld+1, err);
      /* Max cleanable flux */
      in->cleanable[ifld] = ObitDConCleanVisCleanable(in, ifld+1, NULL, err);
      in->beamPeakRMS[ifld] = in->beamPeakRMS[jfld];
    }
  }

} /* end MakeAllResiduals */

/**
 * Find current estimates best and second field  3D
 * Ignore fields with imgPeakRMS<4.0
 * \param in     The Clean object
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
    if (in->imgPeakRMS[i]<4.0) continue;  /*ignore if too low SNR */
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
 * Find current estimates best  field  2D
 * Uses autoShift only for flux densities above 0.1 x autoCenFlux
 * Ignores fields with imgPeakRMS<4.0
 * \param in          The Clean object
 * \param autoCenFlux Cutoff level for autocentering,  0=> no autoCenter
 * \param bbest       [out] 0-rel best field
 */
static void WhosBest2D (ObitDConCleanVis *in, ofloat autoCenFlux, 
			olong *bbest)
{
  ofloat testBest;
  olong i, best;
  gboolean isAuto=FALSE;

  /* Find absolute best */
  best = -1; 
  testBest = -1.0e20;
  for (i=0; i<in->nfield; i++) {
    /* Ignore if SNR too low */
    if (in->imgPeakRMS[i]<4.0) continue;
    if (in->quality[i]>testBest) {
      testBest = in->quality[i];
      best = i;
    }
  } /* end loop finding current best */

  /* Using only autoCenter? */
  for (i=0; i<in->nfield; i++) {
    isAuto = (in->mosaic->isAuto[i]>=0) && (in->imgPeakRMS[i]>=4.0) 
      && (in->maxAbsRes[i]>0.1*autoCenFlux) 
      && (in->quality[i]>=testBest);       /* Is best a strong autoCenter image? */
    if (isAuto) break;
  }
  /* Find best */
  best = -1; 
  testBest = -1.0e20;
  for (i=0; i<in->nfield; i++) {
    /* Ignore autoCenter if not isAuto */
    if (!isAuto && (in->mosaic->isAuto[i]>=0)) continue;
    /* Ignore if SNR too low */
    if (in->imgPeakRMS[i]<4.0) continue;
    if ((in->quality[i]>testBest) &&
	(!isAuto || (in->mosaic->isAuto[i]>=0))) { /* Only autoCenter if isAuto */
      testBest = in->quality[i];
      best = i;
    }
  } /* end loop finding current best */
  
  /* results */
  *bbest   = best;
} /* end WhosBest2D */

/**
 * Find priority ordering for making residual images
 * Ranking by quality, ignored "fresh" and shifted images
 * Only includes images with a quality within 30% of the best
 * Uses autoShift only for flux densities above 0.1 x autoCenFlux
 * Fields with quality=0.0 (reached CLEAN goal) are ignored.
 * If the best is an autoWindow, only other autoWindow fields are considered
 * \param in      The Clean object
 * \param fresh   List of flags indicating freshly made
 *                In order of images in mosaic member of in
 * \param autoCenFlux Cutoff level for autocentering,  0=> no autoCenter
 * \param fldList [out] Ordered list of 1-rel field numbers
 *                The size of the number of images in mosaic+1
 */
static void OrderImage (ObitDConCleanVis *in, gboolean *fresh, 
			ofloat autoCenFlux, olong *fldList)
{
  ofloat *tmpQual=NULL, maxQual, bestQual;
  olong i, n, indx, best;
  gboolean done, isAuto=FALSE;
  
  /* Make temporary array of quality factors */
  tmpQual = g_malloc0(in->mosaic->numberImages*sizeof(ofloat));

  /* Copy quality factors to temporary array */
  n = in->mosaic->numberImages;
  bestQual = -1.0e18;
  for (i=0; i<n; i++) {
    if ((!fresh[i]) && (in->mosaic->isShift[i]<0)&& (in->quality[i]>0.0))
      tmpQual[i] = in->quality[i];
    else tmpQual[i] = -1.0e20;
    if ((tmpQual[i]>bestQual) && (in->imgPeakRMS[i]>=4.0)) {
      bestQual = tmpQual[i];                 /* Best quality */
      isAuto   = (in->mosaic->isAuto[i]>=0)
	&& (in->maxAbsRes[i]>0.1*autoCenFlux); /* Is best a strong autoCenter image? */
    }
  }

  /* Loop finding and dropping best until all gone */
  indx = 0; 
  while (1) {   /* Sort loop */
    maxQual  = -1.0e18;
    best    = -1;
    done = TRUE;
    for (i=0; i<n; i++) {
      if ((tmpQual[i]>maxQual) &&                    /* Best so far */
	  (!isAuto || (in->mosaic->isAuto[i]>=0))) { /* Only autoCenter if isAuto */
	maxQual = tmpQual[i];
	best = i;
	done = FALSE;
      }
    }
    /* Find anything? */
    if (done) break;
    
    /* Save value if within 30% of best, dummy tmpQual */
    if (maxQual>=0.7*bestQual) {
      fldList[indx++] = best+1;
    }
    if ((best>=0) && (best<n)) tmpQual[best]   = -1.0e20;
  }  /* end sort loop  */
  
  fldList[indx++] = 0; /* terminate */
  
  /* Cleanup */
  if (tmpQual) g_free(tmpQual);
  
} /* end OrderImage */

/**
 * Find priority ordering for cleaning residual images
 * Ranking by quality, using only "fresh" and no shifted images
 * Only includes images with a quality within 50% of the best
 * If the best is an autoWindow (>0.1autoCenFlux), only other 
 * autoWindow fields are considered
 * Ignores windows with in->imgPeakRMS<4.0
 * Uses autoShift only for flux densities above 0.1 x autoCenFlux
 * \param in      The Clean object
 * \param fresh   List of flags indicating freshly made
 *                In order of images in mosaic member of in
 * \param autoCenFlux Cutoff level for autocentering, 0=> no autoCenter
 * \param fldList [out] Ordered list of 1-rel field numbers
 *                The size of the number of images in mosaic+1
 */
static void OrderClean (ObitDConCleanVis *in, gboolean *fresh, 
                        ofloat autoCenFlux, olong *fldList)
{
  ofloat *tmpQual=NULL, maxQual, bestQual, testBest;
  olong i, n, indx, best;
  gboolean done, isAuto=FALSE;
  
  /* Find absolute best */
  best = -1; 
  testBest = -1.0e20;
  for (i=0; i<in->nfield; i++) {
    /* Ignore if SNR too low */
    if (in->imgPeakRMS[i]<4.0) continue;
    if (in->quality[i]>testBest) {
      testBest = in->quality[i];
      best = i;
    }
  } /* end loop finding current best */

  /* Make temporary array of quality factors */
  tmpQual = g_malloc0(in->mosaic->numberImages*sizeof(ofloat));

  /* Copy quality factors to temporary array */
  n = in->mosaic->numberImages;
  bestQual = -1.0e18;
  for (i=0; i<n; i++) {
    if (fresh[i] && (in->imgPeakRMS[i]>=4.0) && (in->mosaic->isShift[i]<0))
      tmpQual[i] = in->quality[i];
    else tmpQual[i] = -1.0e20;
    if (tmpQual[i]>bestQual) {
      bestQual = tmpQual[i];                 /* Best quality */
      isAuto   = (in->mosaic->isAuto[i]>=0)
	&& (in->maxAbsRes[i]>0.1*autoCenFlux) 
        && (in->quality[i]>=testBest);        /* Is best a strong autoCenter image? */
  }
  } 

  /* Loop finding an dropping best until all gone */
  indx = 0;
  while (1) {    /* Sort loop */
    maxQual = -1.0e18;
    best    = -1;
    done = TRUE;
    for (i=0; i<n; i++) {
      if ((tmpQual[i]>maxQual) &&                    /* Best so far */
	  (!isAuto || (in->mosaic->isAuto[i]>=0))) { /* Only autoCenter if isAuto */
	maxQual = tmpQual[i];
	best = i;
	done = FALSE;
      }
    }
    
    /* Find anything? */
    if (done) break;
    
    /* Save value if within 50% of best, dummy tmpQual */
    if (maxQual>=0.5*bestQual) {
      /* If this is any autoCenter field and !isAuto, then use the shifted field */
      if (!isAuto && (in->mosaic->isAuto[best]>=0))
	fldList[indx++] = in->mosaic->isAuto[best];
      else
	fldList[indx++] = best+1;
    }
    if ((best>=0) && (best<n)) tmpQual[best]   = -1.0e20;
  } /* end sort loop */
  
  fldList[indx++] = 0; /* terminate */
  
  /* Cleanup */
  if (tmpQual) g_free(tmpQual);
  
} /* end OrderClean */

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
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
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
 * Determine max. abs value of a cleanable pixel in an image.
 * "Cleanable" means inside the outer window and not in an unbox.
 * \param in    The object to deconvolve
 * \param field Field number (1-rel) to test
 * \param pixarray If nonNUILL, use these pixels rather than reading
 * \param err   Obit error stack object.
 * \return Maximum abs pixel value inside outer window but outside unboxes.
 */
static ofloat ObitDConCleanVisCleanable(ObitDConCleanVis *in, olong field, 
					ObitFArray *pixarray, ObitErr *err)
{
  ofloat out = -1.0;
  ObitImage *image=NULL;
  ObitFArray *usePixels=NULL;  
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat PeakIn, PeakOut, RMS;
  olong i,PeakInPos[2] = {0,0};
  gboolean doAbs;
  gchar *routine = "ObitDConCleanVisCleanable";

  /* error checks */
  if (err->error) return out;
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
    return out;
  }

  /* Set output to full image, plane at a time */
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
  
  image = in->mosaic->images[field-1];
  
  /* Set input to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc); 
  dim[0] = 1;
  ObitInfoListAlwaysPut (image->info, "IOBy", OBIT_long, dim, &IOsize);
  
  /* Read */
  if (pixarray) { /* use pixarray */
    usePixels  = ObitFArrayRef(pixarray);
  } else { /* read */
    retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
    retCode = ObitImageRead (image, image->image->array, err);
    retCode = ObitImageClose (image, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, out);
    usePixels = ObitFArrayRef(image->image);  /* Pointer to image buffer */
  }
  
  /* Allow negative for Stokes other than I */
  doAbs = fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1;
    
  /* Get statistics */
  ObitDConCleanWindowStats (in->window, field, usePixels,
			    doAbs,
			    &PeakIn, &PeakInPos[0], &PeakOut, 
			    &RMS, err);
  /* Free Image array */
  usePixels = ObitFArrayUnref(usePixels);
  if (pixarray==NULL) image->image = ObitFArrayUnref(image->image);
  
  if (err->error) Obit_traceback_val (err, routine, image->name, out);
  
  return fabs(PeakOut);
} /* end ObitDConCleanVisCleanable */

/* GLOBAL DEBUG */
ObitDConCleanVis *damn;
/**
 * Low accuracy subtract pixels from image for current CLEAN fields.
 * Loops over fields in in->currentFields
 * Uses BeamPatch to subtract list of components
 * Updates cleanable member on in
 * \param in       The Clean object
 * \param newCC    Start CC -1 per in->mosaic field for subtraction
 * \param pixarray Array of arrays of pixels to subtract
 * \param err      Obit error stack object.
 * \return TRUE if attempted, FALSE if cannot do 
 */
static void SubNewCCs (ObitDConCleanVis *in, olong *newCC, ObitFArray **pixarray, 
		       ObitErr *err)
{
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL;
  ImSubFuncArg **threadArgs;
  ObitFArray **comps=NULL;
  ofloat parms[20];
  olong i, j, ncc, ver, nThreads=0, nTh, nfield, ifield, nDo, nLeft;
  gboolean OK;
  gchar *tabType = "AIPS CC";
  gchar *routine = "SubNewCCs";
  gboolean DebugGDB=FALSE;  /*  DEBUG */

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisIsA(in));

  /* DEBUG */
  damn = in;

  /* Count number of fields */
  nfield = 0;
  for (i=0; i<in->numCurrentField; i++) if (in->currentFields[i]>0) nfield++;
  if (nfield<1) return;

  /* Collect compressed CC tables */
  comps = g_malloc0(nfield*sizeof(ObitFArray*));
  for (i=0; i<nfield; i++) comps[i] = NULL;

  /* Loop over fields */
  for (i=0; i<nfield; i++) {
    ifield = in->currentFields[i];
    if (ifield<0) break;

    /* Something to do? */
    ncc = in->Pixels->iterField[ifield-1];
    if (newCC[i]>=ncc) continue;
    
    /* Get CC table */
    ver = in->skyModel->CCver[ifield-1];
    tempTable = newObitImageTable (in->skyModel->mosaic->images[ifield-1],OBIT_IO_ReadOnly, 
				   tabType, &ver, err);
    if ((tempTable==NULL) || (err->error)) Obit_traceback_msg (err, routine, in->name);
    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);

    /* Get selected CCs after compression */
    comps[i]   = ObitTableCCUtilMergeSel (CCTable, newCC[i]+1, ncc, parms, err);
    CCTable = ObitTableUnref(CCTable);
    Obit_return_if_fail ((comps[i]!=NULL),
			 err, "%s: Error merging CLEAN components for field %d",
			 routine, ifield);

  } /* end loop over fields */

  /* setup for threaded processing
     Initialize Threading */
  nThreads = MakeImSubFuncArgs (in->thread, err, &threadArgs);
  /* No more threads than work to spread them across */
  if (nfield>1) nTh = nThreads;
  else nTh = 1;

  /* Initialize Thread args */
  for (i=0; i<nTh; i++) {
    threadArgs[i]->in     = in;
    threadArgs[i]->inComp = comps;
    threadArgs[i]->nfield = nfield;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
  }

  /* Loop over fields */
  nTh = MIN (nThreads, nfield);
  nLeft = nfield;
  for (i=0; i<nfield; i+=nTh) {
    nDo = MIN (nTh, nLeft);
    /* Args for this batch */
    for (j=0; j<nDo; j++) {
      threadArgs[j]->ofield = in->currentFields[i+j];
      threadArgs[j]->inData = pixarray[i+j];
    }
    /* Dummy any remainders */
    for (j=nDo; j<nTh; j++) {
      threadArgs[j]->ofield = 1;
      threadArgs[j]->inData = NULL;
    }

    /* fprintf (stderr,"fields %d %d nDo %d nLeft %d nfield %d\n",
       threadArgs[0]->ofield, threadArgs[1]->ofield, nDo, nLeft, nfield);*/
    /* Do operation - need to stub remainder of nTh (don't know why) */
    OK = ObitThreadIterator (in->thread, nTh, 
			     (ObitThreadFunc)ThreadImSub,
			     (gpointer**)threadArgs);

    /* fprintf (stderr,"Done OK %d\n", OK);*/
    /* Check for problems */
    if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
    if (err->error) goto cleanup;

    /* DEBUG */
    if (DebugGDB) {
      ObitImageUtilArray2Image ("DbugSubNewCC.fits", 0, pixarray[i], err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
    /* END DEBUG */

    nLeft -= nDo;
    if (nLeft<=0) break;

  } /* end loop over field subtracting */

  
  /* Cleanup */
 cleanup:
  KillImSubFuncArgs (nThreads, threadArgs);
  if (comps) {
    for (i=0; i<nfield; i++) comps[i] = ObitFArrayUnref(comps[i]);
    g_free(comps);
  }
  
} /* end SubNewCCs */

/**
 * Subtract a set of lists of CLEAN components from an image in a thread
 * Determines statistics on completed image set on in.
 * Uses BeamPatches on in.
 * \param arg Pointer to ImSubFuncArg argument with elements:
 * \li in       The ObitDConCleanVis object
 * \li inData   Pixel array to be subtracted from
 * \li field    Field number (1-rel) of data in inData 
 * \li inComp   Array of compressed arrays of point pixels
 *              Correspond to fields in->currentFields
 * \li nfield   Number of elements in inComp
 * \li ithread  thread number, <0 -> no threading
 * \li err      ObitErr Obit error stack object
 * \li thread   thread Object
 * \return NULL
 */
static gpointer ThreadImSub (gpointer args)
{
  /* Get arguments from structure */
  ImSubFuncArg *largs = (ImSubFuncArg*)args;
  ObitDConCleanVis *in  = largs->in;
  ObitFArray *pixarray  = largs->inData;
  olong       ofield    = largs->ofield;
  ObitFArray **inComp   = largs->inComp;
  olong       nfield    = largs->nfield;
  ObitErr    *err       = largs->err;
  ObitThread *thread    = largs->thread;
  /* local */
  olong i, j, ifield, len, pos1[2], pos2[2], PeakInPos[2];
  ofloat PeakIn, PeakOut, RMS;
  gboolean doAbs;
  ofloat idx, idy, ftemp, flux, offset[2];
  ObitFArray *comps;
  ObitImageDesc *inDesc, *outDesc;

  if (err->error) goto finish;
  if (pixarray==NULL) goto finish;
  /* DEBUG
     fprintf (stderr,"thread %d in %p damn %p loc %p locp %p pixarray %p\n",
     largs->ithread, in, damn, largs->in, &largs->in, pixarray->array);
     g_assert (in==damn);
     g_assert (ObitIsA(in, &myClassInfo)); */


  /* Loop over fields */
  for (i=0; i<nfield; i++) {
    ifield = in->currentFields[i];
    if (ifield<=0) break;
    comps = inComp[i];          /* Component list */
    if (comps==NULL) continue;  /* Any components? */
    g_assert (ObitFArrayIsA(comps)); /* DEBUG */
    
    /* Loop over CCs - assumes all points */
    inDesc    = in->mosaic->images[ifield-1]->myDesc;
    outDesc   = in->mosaic->images[ofield-1]->myDesc;
    len       = comps->naxis[0];
    pos2[0]   = in->BeamPatches[ifield-1]->naxis[0]/2; 
    pos2[1]   = in->BeamPatches[ifield-1]->naxis[1]/2;    /* Center of beam patch */
    offset[0] = outDesc->xPxOff + pixarray->naxis[0]/2; 
    offset[1] = outDesc->yPxOff + pixarray->naxis[1]/2;    /* Center of image */
    idx = 1.0 / inDesc->cdelt[0];
    idy = 1.0 / inDesc->cdelt[1];
    for (j=0; j<comps->naxis[1]; j++) {
      /* Get 0-rel pixel numbers */
      ftemp = comps->array[1+j*len] * idx + offset[0];
      if (ftemp>0.0) ftemp += 0.5;
      else           ftemp -= 0.5;
      pos1[0] = (olong)(ftemp); 
      ftemp = comps->array[2+j*len] * idy + offset[1];
      if (ftemp>0.0) ftemp += 0.5;
      else           ftemp -= 0.5;
      pos1[1] = (olong)(ftemp); 
      flux = comps->array[j*len];
      ObitFArrayShiftAdd(pixarray, pos1, in->BeamPatches[ifield-1], pos2, -flux, pixarray);
    } /* End loop over components */
  } /* end loop over fields */
 
  /* Remeasure cleanable flux */
  /* Allow negative for Stokes other than I */
  doAbs = fabs (in->mosaic->images[ofield-1]->myDesc->crval[in->mosaic->images[ofield-1]->myDesc->jlocs]-1.0) > 0.1;

  /* Get statistics */
  ObitDConCleanWindowStats (in->window, ofield, pixarray,
                            doAbs,
                            &PeakIn, &PeakInPos[0], &PeakOut,
                            &RMS, err);
  in->cleanable[ofield-1]  = PeakOut;
  in->imgPeakRMS[ofield-1] = PeakIn/RMS;

  /* Indicate completion */
 finish:
  if (largs->ithread>=0)
    ObitThreadPoolDone (thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadImSub */

/**
 * Make arguments for Threaded ThreadImSub
 * Copy pointers to members rather than Refing them.
 * \param thread     ObitThread object to be used for function
 * \param err        Obit error stack object.
 * \param ThreadArgs [out] Created array of ImSubFuncArg, 
 *                   delete with KillImSubFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeImSubFuncArgs (ObitThread *thread, 
				ObitErr *err, ImSubFuncArg ***ThreadArgs)
{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(ImSubFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(ImSubFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->in         = NULL;
    (*ThreadArgs)[i]->inData     = NULL;
    (*ThreadArgs)[i]->ofield     = 0;
    (*ThreadArgs)[i]->inComp     = NULL;
    (*ThreadArgs)[i]->nfield     = 0;
    (*ThreadArgs)[i]->ithread    = i;
    (*ThreadArgs)[i]->thread     = thread;
    (*ThreadArgs)[i]->err        = err;
  }

  return nThreads;
} /*  end MakeImSubFuncArgs */

/**
 * Delete arguments for ThreadImageStats
 * No objects are Unreffed
 * \param nargs      number of elements in args.
 * \param ThreadArgs Array of ImSubFuncArg
 */
static void KillImSubFuncArgs (olong nargs, ImSubFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillImSubFuncArgs */
