/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2025                                          */
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
#include "ObitImageWB.h"
#include "ObitSkyModelVMBeam.h"
#include "ObitSystem.h"

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

/** Private: Zero all residuals. */
static void  ClearAllResiduals (ObitDConCleanVis *in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitDConCleanVisClassInfoDefFn (gpointer inClass);

/** Private: Get pixel array for a given field. */
static ObitFArray* GetFieldPixArray (ObitDConCleanVis *in, olong field, 
				     ObitErr *err);

/** Private: Set Beam patch, min. flux, decide on SDI CLEAN. */
static void ObitDConCleanVisDecide (ObitDConCleanVis* in, ObitErr *err);

/** Private: Threaded Image subtractor */
static gpointer ThreadImSub (gpointer arg);

/** Private: Make Threaded Image subtraction args */
static olong MakeImSubFuncArgs (ObitThread *thread,
				ObitErr *err, ImSubFuncArg ***ThreadArgs);

/** Private: Delete Threaded Image subtraction args */
static void KillImSubFuncArgs (olong nargs, ImSubFuncArg **ThreadArgs);

/** Private: Low accuracy subtract CLEAN model. */
static void SubNewCCs (ObitDConCleanVis *in, olong *newCC, 
		       ObitFArray **pixarray, ObitErr *err);

/** Private: Create/init PxList. */
static void NewPxList (ObitDConCleanVis *in, ObitErr *err);

/** Private: Create/init Pixarray. */
static ObitFArray** NewPxArray (ObitDConCleanVis *in, olong *startCC, ObitErr *err);

/** Private: Delete Pixarray. */
static ObitFArray** KillPxArray (ObitDConCleanVis *in, ObitFArray **pixarray);

/** Private: Delete BeamPatches. */
static void KillBeamPatches (ObitDConCleanVis *in);

/** Private: Find peak brightness. */
static void FindPeak (ObitDConCleanVis *in, ObitErr *err);

/** Private: Pick next field(s) and get Residual image(s) */
static gboolean PickNext2D(ObitDConCleanVis *in, ObitErr *err);
static gboolean PickNext3D(ObitDConCleanVis *in, ObitErr *err);

/** Private: Find best 3D residual image. */
static void  WhosBest (ObitDConCleanVis *in, olong *best, olong *second);

/** Private: Find best 2D residual image. */
static void WhosBest2D (ObitDConCleanVis *in, ofloat autoCenFlux, olong *best);

/** Private: Priority order for making residual images. */
static void OrderImage (ObitDConCleanVis *in, gboolean *fresh, ofloat autoCenFlux, 
			olong *fldList);

/** Private: Priority order for CLEANing */
static void OrderClean (ObitDConCleanVis *in, gboolean *fresh, ofloat autoCenFlux, 
			olong *fldList);

/** Private: Select tapering for next CLEAN */
static gboolean SelectTaper (ObitDConCleanVis *in, gboolean *fresh, ObitErr *err);

/** Private: Checks if fresh facets are "Done" */
void CheckIfDone(ObitDConCleanVis* in, ObitErr *err);

/** Private: Are any facets modified? */
static gboolean anyDirty (ObitDConCleanVis *in);

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
  out->quality   = ObitMemFree (out->quality);
  out->cleanable = ObitMemFree (out->cleanable);
  out->fresh     = ObitMemFree (out->fresh);

  /* In with the new */
  nfield       = in->nfield;
  out->quality   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  out->fresh     = ObitMemAlloc0Name(nfield*sizeof(gboolean),"Clean fresh");
  for (i=0; i<nfield; i++) {
    out->quality[i]    = in->quality[i];
    out->cleanable[i]  = in->cleanable[i];
    out->fresh[i]      = in->fresh[i];
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
  out->quality   = ObitMemFree (out->quality);
  out->cleanable = ObitMemFree (out->cleanable);
  out->fresh     = ObitMemFree (out->fresh);

  /* In with the new */
  nfield       = in->nfield;
  out->quality   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  out->fresh     = ObitMemAlloc0Name(nfield*sizeof(gboolean),"Clean fresh");
  for (i=0; i<nfield; i++) {
    out->quality[i]   = in->quality[i];
    out->cleanable[i] = in->cleanable[i];
    out->fresh[i]     = in->fresh[i];
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
  out->gain        = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Clean factor");
  out->quality     = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Clean quality");
  out->cleanable   = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Clean cleanable");
  out->fresh       = ObitMemAlloc0Name((nfield+3)*sizeof(gboolean),"Clean fresh");
  out->maxAbsRes   = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Clean avg res");
  out->imgRMS      = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name((nfield+3)*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->fresh[i]       = FALSE;
    out->imgRMS[i]      = -1.0;
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
  out->fresh       = ObitMemAlloc0Name(nfield*sizeof(gboolean),"Clean fresh");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->fresh[i]       = FALSE;
    out->imgRMS[i]      = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
 }

  return out;
} /* end ObitDConCleanVisCreate2 */

/**
 * Determine the maximum number of auto window middle loops
 * This is determined from the ratio of the data volume of the mosaic to uv data.
 * \param in   The CLEAN object, ObitInfoList member:
 * \li "maxAWLoop" OBIT_long scalar = Override for default max middle CLEAN loop
 * \param err Obit error stack object.
 */
olong autoWinLoopCount (ObitDConCleanVis *in, ObitErr *err)
{
  olong i, maxLoop=1, maxReq=0;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  odouble ratio;
  ollong uvcnt, imcnt, ll;
  ObitUVDesc *uvdesc = in->imager->uvdata->myDesc;
  ObitImageDesc *imdesc=NULL;
  ObitImageMosaic *mosaic = in->mosaic;

  if (err->error) return maxLoop;      /* Existing error? */
  if (!in->autoWindow) return maxLoop;  /* No auto windowing? */

  /* Check for user request */
  ObitInfoListGetTest(in->info, "maxAWLoop", &type, dim, &maxReq);
  if (maxReq>0) {
    if (err->prtLv>=2) {
      Obit_log_error(err, OBIT_InfoErr,"User request max autoWin loop = %d",
		     maxReq);
      ObitErrLog(err);
    }
    return maxReq;  
  }

  /* How big is the visibility data */
  uvcnt = uvdesc->nvis;
  uvcnt *= (uvdesc->nrparm + 
	    uvdesc->inaxes[0]*uvdesc->inaxes[1]*uvdesc->inaxes[2]*uvdesc->inaxes[3]);

  /* How big are all the images in the mosaic? */
  imcnt = 0;
  for (i=0; i<mosaic->numberImages; i++) {
    imdesc = (ObitImageDesc*)mosaic->images[i]->myIO->myDesc;
    ll = imdesc->inaxes[0];
    ll *= imdesc->inaxes[1] * imdesc->inaxes[2];
    imcnt += ll;
  }

  /* Ratio */
  ratio = (odouble)uvcnt / (odouble)imcnt;
  if (ratio<0.5) maxLoop = 1;
  if (ratio>1.0) maxLoop = 2;
  if (ratio>2.0) maxLoop = 3;
  if (ratio>3.0) maxLoop = 4;
  /*if (ratio>4.0) maxLoop = 5;*/

  /* diagnostics */
  if (err->prtLv>=2) {
    Obit_log_error(err, OBIT_InfoErr,"UV/im size ratio=%lf, max autoWin loop = %d",
		   ratio, maxLoop);
    ObitErrLog(err);
  }

  return maxLoop;
} /* end autoWinLoopCount */

/**
 * Do deconvolution, uses function on class pointer
 * Does final flatten if FullField member of mosaic member is defined.
 * CLEAN control parameters are in the ObitInfoList member:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 50]
 * \li "maxAWLoop" OBIT_long scalar = Override for default max middle CLEAN loop
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form
 * \li "CCVer"   OBIT_long array    = CLEAN table version for all fields
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimun flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "autoCen" OBIT_float scalar = Threshold level for autocenter
 *               If peak exceeds this value the minFlux reset to 0.09*autoCen
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "maxBeamTaper" OBIT_float scalar = max BeamTaper facets to consider [def 0.0]
 * \li "minBeamTaper" OBIT_float scalar = min BeamTaper facets to consider [def 0.0]
 * \li "MResKnob" OBIT_float array = Resolution selection controls
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanVisDeconvolve (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanVis *in;
  ObitFArray **pixarray=NULL;
  gboolean done, fin=TRUE, quit=FALSE, doSub, bail, doMore, moreClean, notDone;
  gboolean redo, isBeamCor, isDirty, doAutoCen=FALSE;
  olong jtemp, i, *startCC=NULL, *newCC=NULL, count, maxAutoWinLoop, ifld;
  olong redoCnt=0, damnCnt=0, lastIter, lastFld, NoPixelCnt;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  const ObitDConCleanVisClassInfo *inClass;
  const ObitUVImagerClassInfo *imagerClass;
  gchar *HALF="HALF", oldStokes[6];
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

  NoPixelCnt  = 0;  /* How many times minor cycles had no residuals to CLEAN */

  /* Reset highest peak */
  in->peakFlux = -1000.0;
  bail = FALSE;  /* if going to recenter facets, don't need to make final residuals */
  in->autoWinFlux = 0.0;

  /* Get parameters */
  inClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  doAutoCen = in->autoCen<1.0e10;  /* AutoCentering? */

  /* Recenter any autoCenter windows is needed */
  if (in->doRecenter) ObitDConCleanVisRecenter (in, in->imager->uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create Pixel List if needed - has size changed? */
  inClass->NewPxList (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Visibility selection and weighting */
  if (in->doWeight) imagerClass->ObitUVImagerWeight (in->imager, err);
  /* Trap no data and return */
  if (err->error==10) {
    ObitErrClear(err);  /* Change to warning */
    Obit_log_error(err, OBIT_InfoWarn, "%s: NO Data copied for %s", routine, in->name);
    ClearAllResiduals (in, err);  /* Zero Images */
    doSub = inClass->ResetSkyModel (in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    return;
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* If doing SDI Clean save copy of uvwork data */
  if (in->SDIGain>0.0) {
    /* If BeamCor imaging, make sure both parallel hands copied */
    isBeamCor = ObitSkyModelVMBeamIsA(in->skyModel);
    if (isBeamCor) {
      strncpy (oldStokes, HALF, 5);  /* Save old value */
      ObitInfoListGetTest(in->imager->uvwork->info, "Stokes", &type, dim, oldStokes);
      dim[0] = 4; dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (in->imager->uvwork->info, "Stokes", OBIT_string, dim, HALF);
      
    }
    in->SDIdata = newObitUVScratch (in->imager->uvwork, err);
    in->SDIdata = ObitUVCopy (in->imager->uvwork, in->SDIdata, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (isBeamCor) {  /* restore previous value of Stokes */
      dim[0] = 4; dim[1] = dim[2] = 1;
      ObitInfoListAlwaysPut (in->imager->uvwork->info, "Stokes", OBIT_string, dim, oldStokes);
      
    }
  }

  /* No more fields than in mosaic */
  Obit_return_if_fail ((in->nfield == in->mosaic->numberImages),
		       err, "%s: CLEAN and mosaic different number of fields %d %d",
		       routine, in->nfield, in->mosaic->numberImages);

  /* Reset Sky model */
  doSub = inClass->ResetSkyModel (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  if (doSub) inClass->ObitDConCleanSub((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Reset/Init Pixel list if not done in ResetSkyModel */
  if (!doSub) inClass->ResetPixelList (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  for (i=0; i<in->nfield; i++) in->Pixels->minFlux[i] = in->minFlux[i];

  /* Be sure to (re)generate residuals */
  for (i=0; i<in->nfield; i++) in->maxAbsRes[i] = -1.0;

  /* Save actual CC version if not specified */
  if (in->CCver<=0) {
    if (in->Pixels->CCver[0]>0) in->CCver = in->Pixels->CCver[0];
    jtemp = in->CCver;
    ObitInfoListAlwaysPut(in->info, "CCVer", OBIT_long, dim, &jtemp);
  }

  /* Make initial images */
  inClass->MakeAllResiduals (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);
  in->doBeam = FALSE;                               /* Shouldn't need again */

  /* Restart CLEAN here if more CLEANable found at end */
  damnCnt   = lastIter = lastFld = 0;
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
      done = inClass->PickNext3D(in, err);
    else
      done = inClass->PickNext2D(in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Does the IPol peak flux density exceed the autocenter threshold */
    doAutoCen = doAutoCen || (fabs (in->peakFlux)>in->autoCen);
    if ((fabs (in->peakFlux) > in->autoCen) &&
	(in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocs]<1.01)) {
      for (i=0; i<in->nfield; i++) {
	in->minFlux[i] = MAX (0.09*in->autoCen, in->minFlux[i]);
	/* Value that counts is on the PxList */
	in->Pixels->minFlux[i] = in->minFlux[i];
      }
      in->autoCen = 1.0e20;  /* Only once */
      bail = TRUE;  /* if going to recenter facets, don't need to make final residuals */
      Obit_log_error(err, OBIT_InfoErr,"Truncating CLEAN at %g Jy for strong source processing",
		     in->minFlux[0]);
    } /* End autoCenter test*/

    if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
    if (done) break;
 
    /* Running out of time? (80% of allowed time) */
    quit = ObitSystemOutOfTime (0.80, err);
    if (quit) {done=TRUE; in->outaTime=TRUE; break;}

    /* Display/edit windows if enabled */
    if (in->display) {
      quit = ObitDisplayShow (in->display, (Obit*)in->mosaic, in->window, 
			      MAX(1,in->currentFields[0]), err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      if (quit) {done=TRUE; break;}
    }

    /* If using autoWindow, iterate until clean not limited by autoWindow flux */
    pixarray = inClass->NewPxArray (in, startCC, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Get image/beam statistics needed for this cycle */
    notDone = inClass->ObitDConCleanPixelStats((ObitDConClean*)in, pixarray, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Should this CLEAN continue? */ 
    if (!notDone) {
      /* Nothing to do in current list */
      for (ifld=0; ifld<in->nfield; ifld++) {
	if (in->currentFields[ifld] <= 0) break;
	in->quality[in->currentFields[ifld]-1]   = 0.0;
	in->cleanable[in->currentFields[ifld]-1] = 
	  -in->cleanable[in->currentFields[ifld]-1];
      }
    }
    
    /*fprintf (stderr,"DEBUG doMore %d done %d\n", doMore, done);*/

    count = 0;
    maxAutoWinLoop = 1;  /* How many middle auto window loops to allow */
    maxAutoWinLoop = autoWinLoopCount(in, err);
    moreClean = TRUE;
    /* Middle CLEAN loop */
    while (moreClean) {
      /* Check if new box added to window after first pass */
      if ((count>0) && in->autoWindow) 
	doMore = inClass->ObitDConCleanAutoWindow ((ObitDConClean*)in, in->currentFields, 
						   pixarray, err);
      else {
	doMore = TRUE;
      }
      if (err->error) Obit_traceback_msg (err, routine, in->name);
     
      /* Are we done */
      if (!doMore) break;

      /* Number of components before CLEAN */
      for (ifld=0; ifld<in->numCurrentField; ifld++) {
	if (in->currentFields[ifld]>0)
	  newCC[ifld] = in->Pixels->iterField[in->currentFields[ifld]-1]; 
      }

      /* Pick components for this major cycle */
      if ((count==0) || doMore) 
	fin   = inClass->ObitDConCleanSelect((ObitDConClean*)in, pixarray, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      /* Check if any facets are done */
      CheckIfDone(in, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
      else ObitErrClear(err);

      /* Did it just  stop because of autoWindow? */
      moreClean = (in->Pixels->complCode==OBIT_CompReasonAutoWin);
      /* Don't do this forever */
      if ((count>0) && !doMore)      moreClean = FALSE;
      if ((count+1)>=maxAutoWinLoop) moreClean = FALSE;
      /*fprintf (stderr,"DEBUG doMore %d count %d moreClean %d Peak/RMS %g\n", 
	doMore, count, moreClean, in->imgPeakRMS[31]);*/
      if (!moreClean) break;

      /* Subtract these CCs from images in pixarray if possible */
      inClass->SubNewCCs(in, newCC, pixarray, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      count++;
    } /* End middle CLEAN loop */
    
    /* If CLEAN stopped due to OBIT_CompReasonMinFlux check that all fields
       are below the limit */
    if (fin && (in->Pixels->complCode==OBIT_CompReasonMinFlux)) {
      for (ifld=0; ifld<in->nfield; ifld++) {
	if (in->maxAbsRes[ifld]>in->minFlux[ifld]) {
	  in->Pixels->complCode = OBIT_CompReasonAutoWin;
	  fin = FALSE;
	}
      } 
    }    /* end check for min CLEAN */

    /* Check if there were no pixels this cycle */
    if (in->Pixels->complCode==OBIT_CompReasonNoPixel) {
      NoPixelCnt++; /* Count times */
      if (NoPixelCnt<=(2*in->nfield)) fin = FALSE;
    }
    /* Call CLEAN done if no pixels more than 2*nfield */
    moreClean = (!fin) && NoPixelCnt<=(2*in->nfield);
    if (!moreClean) break;

    /* Update quality list for new max value on fields just CLEANed */
    for (ifld=0; ifld<in->nfield; ifld++) {
      if (in->currentFields[ifld] <= 0) break;  /* List terminates? */
      /* This clean finished? */ 
      if (in->Pixels->complCode==OBIT_CompReasonMinFlux) 
	  in->maxAbsRes[in->currentFields[ifld]-1] = 0.0;
      else in->maxAbsRes[in->currentFields[ifld]-1] = in->Pixels->maxResid;
      in->quality[in->currentFields[ifld]-1] = 
	inClass->ObitDConCleanVisQuality((ObitDConCleanVis*)in, in->currentFields[ifld], err);
      in->cleanable[in->currentFields[ifld]-1] = in->maxAbsRes[in->currentFields[ifld]-1];
      in->fresh[in->currentFields[ifld]-1]     = FALSE;
    }

    /* Release working pixel arrays */
    pixarray = inClass->KillPxArray (in, pixarray);

    /* Clear BeamPatches array */
    inClass->KillBeamPatches (in);

    /* Subtract any components from visibility data */
    inClass->ObitDConCleanSub((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Check for sombitch getting stuck */
    if (((in->Pixels->currentIter-lastIter)==1) && 
	(lastFld==in->currentFields[0])) {
      damnCnt++;
    } else {
      damnCnt = 0;
    }
    lastIter = in->Pixels->currentIter;
    lastFld = in->currentFields[0];
    if (damnCnt>10) {
      Obit_log_error(err, OBIT_InfoWarn,"CLEAN in trouble, disable field %d", lastFld);
      in->quality[lastFld-1]   = -1.0;
      in->cleanable[lastFld-1] = -1.0;
      /* If this is a autoCenter 2D shifted image - get its twin */
      if (in->mosaic->isShift[lastFld-1]>0) {
	in->quality[in->mosaic->isShift[lastFld-1]-1]   = -1.0;
	in->cleanable[in->mosaic->isShift[lastFld-1]-1] = -1.0;
      }
      damnCnt = 0;
    }

  } /* end clean loop */
  if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Any facets modified? */
  isDirty = anyDirty (in);

  /* Subtract any remaining components from visibility data */
  if ((in->niter>0) && isDirty)
    inClass->ObitDConCleanSub((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  if (startCC) g_free(startCC);
  if (newCC) g_free(newCC);

  /* Release working pixel arrays */
  pixarray = inClass->KillPxArray (in, pixarray);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Any facets modified? */
  isDirty = anyDirty (in);

  /* Make final residuals */
  if ((!bail) && (in->niter>0) && isDirty) {
    inClass->ObitDConCleanResetChDone((ObitDConClean*)in, err);  /* Any resets needed */
    /* Make all stale fields */
    for (ifld=0; ifld<in->nfield; ifld++) in->currentFields[ifld] = ifld+1;
    in->currentFields[ifld] = 0;
    inClass->MakeResiduals (in, in->currentFields, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
    /* Check if any fields now CLEANable and autoWin and more clean allowed */
    if ((in->autoWindow) && 
	(in->Pixels->currentIter<in->Pixels->niter)) {
      redo = FALSE;
      for (ifld=0; ifld<in->nfield; ifld++) {
	redo = redo || ((in->cleanable[ifld]>=in->minFlux[ifld])
			&& (in->imgPeakRMS[ifld]>=5.0));
      }
      /* Something left to CLEAN? */
      if (redo && (redoCnt<1) && (!quit)) {
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
    if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
  }
  
  if (in->doXRestore && !bail) {
    /* Cross Restore if multiple overlapping fields */
    inClass->ObitDConCleanXRestore((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
  }

  /* Flatten if needed */
  if (in->doFlatten && !bail) {
    inClass->ObitDConCleanFlatten((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);

    /* If 2D imaging concatenate CC tables */
    if (!in->mosaic->images[0]->myDesc->do3D) 
      ObitImageMosaicCopyCC (in->mosaic, in->imager->uvwork, err);
    
    /* Display flattened image if enabled */
    if (in->display && in->mosaic->FullField) 
      ObitDisplayShow (in->display, (Obit*)in->mosaic->FullField, NULL, 
		       1, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Find peak brightness in mosaic */
  inClass->FindPeak (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Clear BeamPatches array - to be sure */
  inClass->KillBeamPatches (in);
  
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
 *                 If peak exceeds this value the minFlux reset to 0.09*autoCen
 * \li "dispURL" OBIT_string scalar = URL of display server
 * \li "maxBeamTaper" OBIT_float scalar = max BeamTaper facets to consider [def 0.0]
 * \li "minBeamTaper" OBIT_float scalar = min BeamTaper facets to consider [def 0.0]
 * \li "MResKnob" OBIT_float array = Resolution selection controls
 * \param in  The CLEAN object as base class
 * \param err Obit error stack object.
 */
void  ObitDConCleanVisGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanVis *in = (ObitDConCleanVis*)inn;  /* as this class */
  ObitDConClassInfo *ParentClass;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong i, n;
  ofloat *farr;
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
  ObitInfoListGetTest(in->info, "doRecenter", &type, dim, &in->doRecenter);

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

  /* get prtLv on err */
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &err->prtLv);

  /* Max BeamTaper to consider */
  ObitInfoListGetTest(in->info, "maxBeamTaper", &type, dim, &in->maxBeamTaper);

  /* Min BeamTaper to consider */
  ObitInfoListGetTest(in->info, "minBeamTaper", &type, dim, &in->minBeamTaper);

  /* Resolution selection controls */
  if (ObitInfoListGetP(in->info, "MResKnob", &type, dim, (gpointer)&farr)) {
    n = MIN (10, dim[0]);
    for (i=0; i<n; i++) in->MResKnob[i] = farr[i];
  }
} /* end ObitDConCleanVisGetParms */

/**
 * Set default CLEAN windows in mosaic
 * If mosaic member  Radius>0 then make round boxes on Fly's eye field
 * with this radius, else use rectangular box including all but outer 5 pixels
 * On outlier fields, use rectangular box of width OutlierSize.
 * If CLEANBox defined in in->info then its contents are used for field 1.
 * Assumes all images in mosaic have descriptors defined.
 * If there are multiple tapers, any window on facet 1 are copied to
 * the corresponding facets at other resolutions.
 * Uses base class function.
 * \param in   The CLEAN object
 * \param err Obit error stack object.
 */
void ObitDConCleanVisDefWindow(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVis *in;
  const ObitDConCleanVisClassInfo *inClass;
  ObitDConCleanWindowType type;
  olong ifield, i, *window;
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

  /* Multiple resolutions? Copy first facet window to corresponding tapered ones */
  if (in->mosaic->numBeamTapes>1) {
    for (ifield=1; ifield<in->mosaic->numberImages; ifield++) {
      if ((in->mosaic->FacetNo[ifield]==0) && 
	  (in->mosaic->BeamTaper[ifield]>0.0)) {
	for (i=0; i<=in->window->maxId[0]; i++ ) {
	  if (ObitDConCleanWindowInfo (in->window, 1, i, &type, &window, err))
	    ObitDConCleanWindowUpdate (in->window, ifield+1, i, type, window, err);
	  if (err->error) Obit_traceback_msg (err, routine, in->name);
	}
      }
    }
  }

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
  gboolean Fl=FALSE, doCalSelect, subbed=FALSE;
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
    for (ifld=0; ifld<in->nfield; ifld++) {
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

  /* Update CC counts - was anything actually subtracted? */
  for (i=0; i<in->mosaic->numberImages; i++) {
    subbed = subbed || in->skyModel->endComp[i]>=in->skyModel->startComp[i];
    in->skyModel->startComp[i] = in->skyModel->endComp[i]+1;
  }
  /* Update Fresh? */
  if (subbed) {
    for (i=0; i<in->mosaic->numberImages; i++) {
      in->fresh[i] = FALSE;  /* Need to remake all images? */
      in->mosaic->images[i]->fresh = FALSE; 
    }
  }

  /* Reset max residual on Pixel List */
  in->Pixels->resMax    = -1.0e20;  /* Maximum residual */

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
static gboolean PickNext2D(ObitDConCleanVis *in, ObitErr *err)
{
  olong i, j, best=-1, lastBest=-1, loopCheck, indx, NumPar;
  olong *fldList=NULL,*fldList2=NULL;
  gboolean doBeam=FALSE, done=TRUE, found, OK, doGridGPU=FALSE;
  ofloat sumwts, autoCenFlux=0.0;
  ObitImage *theBeam=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  const ObitUVImagerClassInfo *imagerClass;
  const ObitDConCleanVisClassInfo *inClass;
  gchar *routine = "PickNext2D";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanVisIsA(in));

  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* clean class structure */

  /* GPU Gridding? */
  if (ObitInfoListGetTest(in->imager->uvwork->info, "doGPUGrid", &type, dim, &doGridGPU)) 

  /* Check if reached max number of components and some done */
  if ((in->Pixels->currentIter >= in->Pixels->niter) && 
      (in->Pixels->currentIter>0)) return done;

  /* If only one field - it's the one */
  if (in->nfield==1) {
    in->currentFields[0] = 1;
    in->currentFields[1] = 0;
    in->numCurrentField  = 1;
    /* Remake if not first pass */
    if (in->Pixels->currentIter>0)
      inClass->MakeResiduals(in, in->currentFields, FALSE, err);
    /* Is this field OK to image? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[0]>1.0 ||
	in->cleanable[0]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[0]))));
    else
      OK = in->imgPeakRMS[0]>1.0;
    done = (in->cleanable[0] <= in->minFlux[0]) || 
      ((in->Pixels->maxResid <= in->minFlux[0]) && (in->Pixels->currentIter>1)) || (!OK);
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium1");
    
    in->peakFlux = MAX (in->peakFlux, in->maxAbsRes[0]);
    return done;
  } /* end single field */
  
  /* Check if there are autoCenter windows, if so get level */
  for (i=0; i<in->mosaic->numberImages; i++) {
    if (in->mosaic->isAuto[i]>0) {
      ObitInfoListGetTest(in->mosaic->images[i]->info, "autoCenFlux", &type, 
			  dim, &autoCenFlux);
      break;
    }
  } /* end autoCen loop */
 
  fldList = ObitMemAlloc0((in->nfield+3)*sizeof(olong));
  fldList2 = ObitMemAlloc0((in->nfield+3)*sizeof(olong));
  
  /* First time? */
  if (in->Pixels->currentIter > 0) {
    
    /* Select taper for next CLEAN / remake images */
    done = inClass->SelectTaper (in, in->fresh, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium2");
    if (done) return done;
    
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
      }
    } /* end loop initializing fields */
    

    /* Make residual images */
    if (fldList[0]>0) inClass->MakeResiduals(in, fldList, doBeam, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    
    /* Check if reached max number of components */
    if (in->Pixels->currentIter >= in->Pixels->niter) return done;
    
    /* Ignore fields already known to be finished or low SNR, see if all done */
    done = TRUE;
    for (i=0; i<in->nfield; i++) {
      /* Is this field OK to image? */
      if (in->autoWindow)
	OK = in->imgPeakRMS[i]>1.0 ||
	  in->cleanable[i]/in->imgRMS[i]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
      else
	OK = in->imgPeakRMS[i]>1.0;
      if ((!OK) || (in->maxAbsRes[i] <= in->minFlux[i])) {
	in->quality[i]   = 0.0;
	in->cleanable[i] = -in->cleanable[i];
      }  else done = FALSE;
    } /* end loop over fields */
    
    /* anything left? */
    if (done) {ObitMemFree(fldList); ObitMemFree(fldList2); return done;} 
    
  } else {
    /* Select taper for next CLEAN */
    done = inClass->SelectTaper (in, in->fresh, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium3");
    if (done) return done;
  }
  
  /* Shouldn't need to make beams again */
  in->doBeam = FALSE;
  
  /* How many images in parallel? */
  if (doGridGPU) 
    NumPar = ObitUVGridGetNumPar(in->imager->uvwork, (Obit*)in->mosaic, doBeam, err);  /* GPU */
  else {
    imagerClass = (ObitUVImagerClassInfo*)in->imager->ClassInfo;
    NumPar = imagerClass->ObitUVImagerGetNumPar (in->imager, in->doBeam, err);  /* CPU */
  }
  NumPar = MIN (NumPar, in->nfield);  /* No more than what's there */
  
  /* Loop remaking blocks of images until something suitable to CLEAN */
  done = FALSE;
  loopCheck = 0;
  while (!done) {
    
    /* Get ordered list */
    if (in->mosaic->numBeamTapes<=1) 
      inClass->OrderImage (in, in->fresh, autoCenFlux, fldList); 
    else {
      inClass->OrderImage (in, in->fresh, autoCenFlux, fldList2); 
      /* Only ones with current selected taper */
      i = j = 0;
      while (fldList2[i]>0) {
	if ((in->mosaic->BeamTaper[fldList2[i]-1]>=in->minBeamTaper) &&
	    (in->mosaic->BeamTaper[fldList2[i]-1]<=in->maxBeamTaper)) {
	  /*Obit_log_error(err, OBIT_InfoErr, "Current Taper f %d j %d tap %f min %f max %f",
			 fldList2[i], j, in->mosaic->BeamTaper[fldList2[i]-1], 
			 in->minBeamTaper,in->maxBeamTaper);*/
	  fldList[j++] = fldList2[i];
	}
	i++;
      }
      fldList[j] = 0;
    } /* end if tapering */
    
    /* No more than NumPar */
    fldList[NumPar] = 0;
    
    /* Count */
    in->numCurrentField = 0;
    for (i=0; i<in->nfield; i++) {
      if (fldList[i]>0) in->numCurrentField++;  /* Count */
      else break;
    }
    
    /* Make residual images if needed */
    if (fldList[0]>0) inClass->MakeResiduals(in, fldList, FALSE, err);
    /* If this is not the first call we should need to make residuals 
       unless multiple resolutions are being used */
    else if ((in->Pixels->currentIter>0) && 
	     (in->mosaic->numBeamTapes<=1)) break;  /* Apparently nothing to do for now */
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    
    /* Which ones remade? */
    for (i=0; i<in->nfield; i++) {
      if (fldList[i]==0) break;
      /* Ignore fields with max residual less than min. */
      if (in->maxAbsRes[fldList[i]-1] < in->minFlux[fldList[i]-1]) {
	in->quality[fldList[i]-1]   = 0.0;
	in->cleanable[fldList[i]-1] = -in->cleanable[fldList[i]-1];
      }
      in->fresh[fldList[i]-1] = TRUE;  /* just to be sure */
    }
    
    /* See if all done */
    done = TRUE;
    for (i=0; i<in->nfield; i++) {
      if (in->maxAbsRes[i] > in->minFlux[i]) done = FALSE;
    }
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium4");
    if (done) {ObitMemFree(fldList); ObitMemFree(fldList2); return done;} 

    /* Which fields ready to CLEAN? 
     if the best one is an autoCenter field only it is done */
    inClass->OrderClean (in, in->fresh, autoCenFlux, fldList);
    
    /* Find current estimates best field */
    inClass->WhosBest2D (in, autoCenFlux, &best);

    /* No best means we're done */
    /* Tell if stopping */
    if ((best<0) && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium5");
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
    
    /* Test BeamTaper */
    if (((in->mosaic->BeamTaper[best]<in->minBeamTaper) || 
	 (in->mosaic->BeamTaper[best]>in->maxBeamTaper))) continue;
    
    /* Don't loop forever */
    loopCheck++;
    if (loopCheck>2*in->nfield) break;
    
    /* Message */
    if (fldList[0]>0) {
      /* best needs to be significantly (>5%) better */
      if ((in->quality[best]-in->quality[fldList[0]-1])<1.05*in->quality[best]) break;
      Obit_log_error(err, OBIT_InfoWarn,
		     "%s:  There may be something better - %f (%d)<%f (%d)",
		     routine, in->quality[fldList[0]-1], fldList[0], 
		     in->quality[best], best+1);
    } else
      Obit_log_error(err, OBIT_InfoWarn,
		     "%s: Nothing - try again (best %d)",
		     routine, best);
    /* Reset in->maxAbsRes,cleanable for any fields just made (fresh) but deemed unusable */
    for (i=0; i<in->nfield; i++) {
      if (in->fresh[i]) 
	{in->maxAbsRes[i] = in->cleanable[i] = 0.0; in->fresh[i]=FALSE;}
    }

  } /* end loop reimaging */

  /* Get highest abs value in any CLEAN window */
  in->peakFlux = 0.0;
  for (i=0; i<in->nfield; i++) {
    in->peakFlux = MAX (in->peakFlux, in->maxAbsRes[i]);
  }

  /* Make sure there is a list of fields */
  if (fldList[0]<=0) inClass->OrderClean (in, in->fresh, autoCenFlux, fldList);

  /* trap premature termination if autoCen */
  if ((autoCenFlux>0.0) && (fldList[0]<=0) && (best>0) && (best<=in->nfield) &&
      (in->maxAbsRes[fldList[best]-1] <= in->minFlux[fldList[best]-1])) {
    fldList[0] = best+1;
    if (err->prtLv>=3) {
       Obit_log_error(err, OBIT_InfoWarn,
		     "%s: best %d maxAbsRes %f cleanable %f",
		      routine, best, in->maxAbsRes[fldList[best]-1], 
		      in->cleanable[fldList[best]-1]);
    }
 }

  /* We're done if best single resolution field has maxAbsRes<minFlux */
  done = fldList[0]<=0;
  done = done || ((in->mosaic->numBeamTapes<=1) && 
		  (in->maxAbsRes[fldList[0]-1] <= in->minFlux[fldList[0]-1]));

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
 
  /* Check if all below minFlux */
  OK = FALSE;
  for (i=0; i<in->numCurrentField; i++) {
    OK = OK || (in->maxAbsRes[fldList[i]-1] > in->minFlux[fldList[i]-1]);
  }
  if ((in->numCurrentField>=1) && !OK) done = FALSE;


  /* Tell if stopping */
  if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium6");
  if (done && (err->prtLv>=3)) {
      Obit_log_error(err, OBIT_InfoWarn,
		     "%s: best %d maxAbsRes %g cleanable %g",
		     routine, best, in->maxAbsRes[fldList[best]-1], 
		     in->cleanable[fldList[best]-1]);
    }

  /* cleanup */
  fldList = ObitMemFree(fldList);
  fldList2 = ObitMemFree(fldList2);

  return done;
} /* end PickNext2D */

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
static gboolean PickNext3D(ObitDConCleanVis *in, ObitErr *err)
{
  olong i, ip, best, second, nextBest, nextSecond, loopCheck, indx, NumPar;
  olong *fldList;
  gboolean OK, doBeam=FALSE, done=TRUE;
  ofloat sumwts;
  ObitImage *theBeam=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  const ObitUVImagerClassInfo *imagerClass;
  const ObitDConCleanVisClassInfo *inClass;
  gchar *routine = "PickNext3D";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanVisIsA(in));

  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* clean class structure */

  /* Check if reached max number of components and some done */
  if ((in->Pixels->currentIter >= in->Pixels->niter) && 
      (in->Pixels->currentIter>0)) return done;

  /* If only one field - it's the one */
  if (in->nfield==1) {
    in->currentFields[0] = 1;
    in->currentFields[1] = 0;
    in->numCurrentField  = 1;
    /* Remake if not first pass */
    if (in->Pixels->currentIter>0)
      inClass->MakeResiduals(in, in->currentFields, FALSE, err);
    /* Is this field OK to image? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[0]>1.0 ||
	in->cleanable[0]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[0]))));
    else
      OK = in->imgPeakRMS[0]>1.0;
    done = (in->cleanable[0] <= in->minFlux[0]) || 
      ((in->Pixels->maxResid <= in->minFlux[0]) && (in->Pixels->currentIter>1)) || (!OK);
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium7");

    in->peakFlux = MAX (in->peakFlux, in->maxAbsRes[0]);
    return done;
  }

  fldList = ObitMemAlloc0((in->nfield+3)*sizeof(olong));

  /* How many images in parallel? */
  imagerClass = (ObitUVImagerClassInfo*)in->imager->ClassInfo;
  NumPar = imagerClass->ObitUVImagerGetNumPar (in->imager, in->doBeam, err);
  NumPar = MIN (NumPar, in->nfield);  /* No more than what's there */
  /* No point in doing too many */
  NumPar = MIN (NumPar, 2); 

  /* First time? */
  if (in->Pixels->currentIter > 0) {  /* No */
    
    /* Select taper for next CLEAN / reimage */
    done = inClass->SelectTaper (in, in->fresh, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    if (done) return done;
    
    /* Make sure all fields initialized */
    indx = 0;
    for (i=0; i<in->nfield; i++) {
      if ((in->maxAbsRes[i] < 0.0) &&
	  (((in->mosaic->BeamTaper[i]>=in->minBeamTaper) && 
	    (in->mosaic->BeamTaper[i]<=in->maxBeamTaper)))) {
	doBeam = in->doBeam;
	/* Make sure beam is made - is SUMWTS present? */
	theBeam = (ObitImage*)in->mosaic->images[i]->myBeam;
	if (!ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts))
	  doBeam = TRUE;
	/* Make residual image - get statistics */
	fldList[indx++] = i+1;
      }
      /* end loop initializing fields */
    }
    /* Make images */
    if (fldList[0]>0) inClass->MakeResiduals (in, fldList, doBeam, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);

    /* Check if reached max number of components */
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium8");
    if (in->Pixels->currentIter >= in->Pixels->niter) return done;
    
    /* Ignore fields already known to be finished or low SNR, see if all done */
    done = TRUE;
    for (i=0; i<in->nfield; i++) {
      /* Is this field OK to image? */
      if (in->autoWindow)
	OK = in->imgPeakRMS[i]>1.0 ||
	  in->cleanable[i]/in->imgRMS[i]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
      else
	OK = in->imgPeakRMS[i]>1.0;
      if ((!OK) || (in->maxAbsRes[i] <= in->minFlux[i]) || 
	  (in->cleanable[i] <= in->minFlux[i])) {
	in->quality[i]   = 0.0;
	in->cleanable[i] = -in->cleanable[i];
      } else done = FALSE;
    } /* end loop over fields */
    /* anything left? */
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium9");
    if (done) {ObitMemFree(fldList); return done;} 
    
    /* Find current estimates best and second best field field */
    inClass->WhosBest (in, &best, &second);

    /* No best means we're done */
    /* Tell if stopping */
    if ((best<0) && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium10");
    if (best<0) return TRUE;
    
    /* Verify that it's still best which may require iteration */
    done = FALSE;
    loopCheck = 0;
    while (!done) {
      ip = 0;
      /* if autoWin, Add any facets with peaks higher than best */
      if (in->autoWindow) {
	for (i=0; i<in->nfield; i++) {
	  OK = (in->imgPeakRMS[i]>4.0) && (in->cleanable[i]>=in->cleanable[best]);
	  OK = OK && (i!=best) && (!in->fresh[i]);
	  /* Test BeamTaper */
	  OK = OK && ((in->mosaic->BeamTaper[i]>=in->minBeamTaper) && 
		      (in->mosaic->BeamTaper[i]<=in->maxBeamTaper));
	  if (OK) {
	    fldList[ip++] = i+1;
	    in->fresh[i]   = TRUE;
	    if (err->prtLv>2) {
	      Obit_log_error(err, OBIT_InfoErr, 
			     "PickNext3D: image %d higher peak %g than best %g",
			     i+1, in->cleanable[i], in->cleanable[best]);
	    }
	  }
	}
      } /* End if autoWin */
     
      /* If not just created recalculate */
      if ((best>=0) && (!in->fresh[best])) {
	
	/* Make best two residual image - get statistics */
	if (!in->fresh[best]) fldList[ip++] = best+1;
	if ((second>=0) && (!in->fresh[second]) && (NumPar>1)) fldList[ip++] = second+1;
	else fldList[ip++] = 0;
	fldList[ip] = 0;
	inClass->MakeResiduals(in, fldList, FALSE, err);
	if (err->error) Obit_traceback_val (err, routine, in->name, done);
	/* Ignore fields with max residual or cleanable less than min. */
	if ((in->maxAbsRes[best] < in->minFlux[best]) || 
	    (in->cleanable[best] < in->minFlux[best])) {
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
          inClass->WhosBest (in, &nextBest, &nextSecond);
      if ((second<0) || (in->quality[best] >= in->quality[nextBest])) break;
      
      /* try again */
      if ((in->quality[best]>0.0) && (err->prtLv>1)) 
	Obit_log_error(err, OBIT_InfoWarn, 
		       "%s: field %d (%g) not as good as second %g",
		       routine, best+1, in->quality[best],  in->quality[nextBest]);
      
      /* Find current estimates best and second best field field */
      inClass->WhosBest (in, &best, &second);
      
    } /* end loop finding best field */
    
  } else {
    /* First pass - all images should be OK */  

    /* Select taper for next CLEAN */
    done = inClass->SelectTaper (in, in->fresh, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    /* Tell if stopping */
    if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium11");
    if (done) return done;

    inClass->WhosBest (in, &best, &second);
  }

  /* Shouldn't need to make beams again */
  in->doBeam = FALSE;
    
  /* publish decision */
  in->currentFields[0] = best + 1;
  in->currentFields[1] = 0;
  in->numCurrentField  = 1;

  /* Diagnostic */
  if (err->prtLv>2) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "PickNext3D: best %d quality %g cleanable %g maxRes %g SNR %g",
		   best+1, in->quality[best],  in->cleanable[best], 
		   in->maxAbsRes[best], in->imgPeakRMS[best]);
  }
  
  /* highest abs value ever found in a CLEAN window */
  if (best>=0) in->peakFlux = MAX (in->peakFlux, in->maxAbsRes[best]);

  /* We're done if best field has maxAbsRes<minFlux */
  if (best>=0) done = in->maxAbsRes[best] <= in->minFlux[best];

  if (best<0) done = TRUE; /* Anything survive? */

  /* cleanup */
  fldList = ObitMemFree(fldList);

  /* Tell if stopping */
  if (done && (err->prtLv>=2)) Obit_log_error(err, OBIT_InfoErr, "Met stopping criterium12");
  return done;
} /* end PickNext3D */

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
 * A new (256x256) field is added centered on the offending source and a negative
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
#define MAXAUTOCEN 200
  gboolean redo = FALSE;
  ObitTableCC *CCTab=NULL;
  ObitImageDesc *imDesc, *imDesc2;
  ObitImageMosaic *mosaic = in->mosaic;
  ObitImageMosaicClassInfo* mosaicClass; 
  ofloat freset, tol;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  ObitDConCleanWindowType otype;
  olong   nfield, ifield, jfield, itemp, noParms, nccpos, nx=0, ny, nplane, nprior;
  olong  CCVer, newField, win[3], inaxes[2], *owin;
  ofloat tmax, xcenter, ycenter, xoff, yoff, radius, cells[2], pixel[2], opixel[2];
  ofloat xcen, ycen, RAShift, DecShift, deltax, deltay, delta, *farray, dx, dy;
  odouble pos[2], pos2[2], RAPnt, DecPnt;
  gboolean done, outside=FALSE, facetDone, clear, Tr=TRUE, Fl=FALSE;
  const ObitDConCleanVisClassInfo *inClass;
  olong ibox, autoCenNum[MAXAUTOCEN], nAutoCen=0;
  odouble autoCenPos[2*MAXAUTOCEN];
  gchar *routine = "ObitDConCleanVisReimage";

  /* Error checks */
  if (err->error) return redo;  /* previous error? */
  g_assert(ObitDConCleanVisIsA(in));

  /* Only I Pol */
  if (in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocs]>1.0) return redo;

  inClass     = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* clean class structure */

  /* Number of fields */
  nfield = mosaic->numberImages;
  nprior = nfield;

  /* Get cellsize */
  cells[0] =  fabs(mosaic->xCells); cells[1] = fabs(mosaic->yCells);

  /* Consider components within 2.5  cells or beam  */
  radius = MAX(2.5 * cells[0],in->mosaic->images[0]->myDesc->beamMaj);

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
    mosaic->images[ifield]->extBuffer = TRUE;  /* Don't need buffer here */
    mosaic->images[ifield]->image = ObitFArrayUnref(mosaic->images[ifield]->image);
    ObitImageOpen (mosaic->images[ifield], OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Make temporary CC table object including CCs from any overlapping facets */
    CCTab = ObitImageMosaicCombineCC (mosaic, ifield+1, CCVer, err);
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
      /* DEBUG HACK -0.5* */
      pixel[0] = imDesc->crpix[0] + ((xcenter-0.5*xoff) / imDesc->cdelt[0]); 
      pixel[1] = imDesc->crpix[1] + ((ycenter-0.5*yoff) / imDesc->cdelt[1]); 
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
      
      imDesc = mosaic->images[ifield]->myDesc;
      if (imDesc->do3D) {
	dx = xcenter;
	dy = ycenter;
      } else { /* 2D */
	dx = xcenter + imDesc->xPxOff*imDesc->cdelt[0];
	dy = ycenter + imDesc->yPxOff*imDesc->cdelt[1];
      }

      /* See if peak is near the center of the image - if within 1/2 pixel
	 and exceeds threshold adjust field and mark as autoCenterField */ 
      if (((fabs(dx+xoff)<0.5*cells[0]) && 
	   (fabs(dy+yoff)<0.5*cells[1])) ) {
	
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
	if (err->prtLv>1) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Modify autoCenter field %d by %5.2f %5.3f pixels",
			 ifield+1,xoff/imDesc->cdelt[imDesc->jlocr],
			 yoff/imDesc->cdelt[imDesc->jlocd]);
	}
	done = TRUE;
      } /* end bright source at field center */
      if (done) goto doNext;  /* Has this one already been done? */
      
      /* See if shift needed */
      if (((fabs(dx+xoff)>tol*cells[0]) || 
	   (fabs(dy+yoff)>tol*cells[1])) ) {
	imDesc = mosaic->images[ifield]->myDesc;
	clear = TRUE;  /* Clear CC table */
	redo  = TRUE;
	
	/* Add field to mosaic */
	/*ObitImageDescGetPoint (imDesc, &RAPnt, &DecPnt); No - use phase center */
	RAPnt  = uvdata->myDesc->crval[uvdata->myDesc->jlocr];
	DecPnt = uvdata->myDesc->crval[uvdata->myDesc->jlocd];
	ObitSkyGeomShiftXY (RAPnt, DecPnt, ObitImageDescRotate(imDesc),
			    pos[0], pos[1], &RAShift, &DecShift);
	nx = ny = ObitFFTSuggestSize (256); nplane = 1;
	mosaicClass = (ObitImageMosaicClassInfo*)in->mosaic->ClassInfo;
	mosaicClass->ObitImageMosaicAddField (in->mosaic, uvdata, nx, ny, nplane, 
					      RAShift, DecShift, TRUE, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	/* Get size actually made */
	nx = mosaic->images[mosaic->numberImages-1]->myDesc->inaxes[0];
	ny = mosaic->images[mosaic->numberImages-1]->myDesc->inaxes[1];

	/* Mark as an autoCenter Image */
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(mosaic->images[mosaic->numberImages-1]->info, 
			      "autoCenField", OBIT_bool, dim, &Tr);
	ObitInfoListAlwaysPut(mosaic->images[mosaic->numberImages-1]->info, 
			      "autoCenFlux", OBIT_float, dim, &freset);
	
	/* Add new field to window */
	inaxes[0] = nx; inaxes[1] = ny;
	newField = ObitDConCleanWindowAddField (in->window, inaxes, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Save position of unbox */
	if (nAutoCen<MAXAUTOCEN) {
	  autoCenNum[nAutoCen]     = newField-1;
	  autoCenPos[nAutoCen*2]   = pos[0];
	  autoCenPos[nAutoCen*2+1] = pos[1];
	  nAutoCen++;
	} else { /* oh bother - blew array */
	  Obit_log_error(err, OBIT_InfoWarn,"Exceeded maximum number of autoCenter fields %d", 
			 MAXAUTOCEN);
	}

	/* Add inner and outer windows */
	xcen = mosaic->images[newField-1]->myDesc->crpix[0] - 
	  mosaic->images[newField-1]->myDesc->xPxOff;
	ycen = mosaic->images[newField-1]->myDesc->crpix[1] - 
	  mosaic->images[newField-1]->myDesc->yPxOff;
	win[0] = 5; win[1] = (olong)(xcen+0.5); win[2] = (olong)(ycen+0.5); 
	ObitDConCleanWindowAdd (in->window, newField, OBIT_DConCleanWindow_round,
				win, err);
	win[0] = (nx/2)-10; win[1] = (olong)(xcen+0.5); win[2] = (olong)(ycen+0.5); 
	/* Set size for new Wideband image */
	if (ObitImageWBIsA(mosaic->images[mosaic->numberImages-1])) win[0] = 22;
	ObitDConCleanWindowOuter (in->window, newField, OBIT_DConCleanWindow_round,
				  win, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Add field to Clean */
	inClass->ObitDConCleanVisAddField (in, uvdata, err);
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
	if (err->prtLv>1) {
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
      ObitImageMosaicFlagCC (CCTab, nccpos, 10.*radius, xcenter, ycenter, err);
      if  (err->error) Obit_traceback_val (err, routine, CCTab->name, redo);

    } /* end loop over multiple peaks */

    /* Delete temporary combined CC table */
    if (CCTab && (ifield>=0)) {
      ObitDataZapTable((ObitData*)mosaic->images[ifield], CCTab->tabType, 
		       CCTab->tabVer, err);
      if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
      CCTab = ObitTableCCUnref(CCTab);
    }

    /* Make temporary CC table to clear rows */
    noParms = 0;
    CCTab = newObitTableCCValue ("Temp CC", (ObitData*)mosaic->images[ifield],
				 &CCVer, OBIT_IO_ReadWrite, noParms, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
    if (clear) ObitTableClearRows ((ObitTable*)CCTab, err); /* Remove the entries and redo */
    CCTab = ObitTableCCUnref(CCTab);
    if (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
    
    /* Close/update image */
    ObitImageClose(mosaic->images[ifield], err);
    mosaic->images[ifield]->extBuffer = FALSE;  /* May need buffer later */
    mosaic->images[ifield]->image = ObitFArrayUnref(mosaic->images[ifield]->image);
   if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
  } /* end loop  L500: */

  /* Loop over autoCen position adding unboxes */
  for (ibox=0; ibox<nAutoCen; ibox++) {
    ifield = autoCenNum[ibox];
    imDesc = mosaic->images[ifield]->myDesc;
    pos[0] = autoCenPos[ibox*2];
    pos[1] = autoCenPos[ibox*2+1];
    /* Put unwindow on this position in all prior fields in which it occured */
    for (jfield=0; jfield<in->window->nfield; jfield++) { 
      /* Not if same fieldor an autoCen field  */
      if ((jfield!=ifield) && (mosaic->isAuto[jfield]<0)) {
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
	/* if (outside) continue; NO - add one anyway */
	
	/* OK - add unbox */
	win[0] = (nx/2)-15; win[1] = (olong)(opixel[0]+0.5); win[2] = (olong)(opixel[1]+0.5); 
	if (outside) win[0] = 20;
	/* Set size for new Wideband image */
	if (ObitImageWBIsA(mosaic->images[mosaic->numberImages-1])) win[0] = 20;
	/* If this is a previous autoCenter window make unbox smaller */
	if ((jfield+1)>nprior) win[0] = 8;
	ObitDConCleanWindowAdd (in->window, jfield+1, OBIT_DConCleanWindow_unround,
				win, err);
      } /* End not same field */	
    } /* end loop adding unboxes */
  } /* end loop over autoCen positions */
  
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
	mosaicClass = (ObitImageMosaicClassInfo*)in->mosaic->ClassInfo;
	mosaicClass->ObitImageMosaicAddField (in->mosaic, uvdata, nx, ny, nplane, 
					      RAShift, DecShift, FALSE, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	/* Get size actually made */
	nx = mosaic->images[mosaic->numberImages-1]->myDesc->inaxes[0];
	ny = mosaic->images[mosaic->numberImages-1]->myDesc->inaxes[1];

	/* Add field to window */
	inaxes[0] = nx; inaxes[1] = ny;
	newField = ObitDConCleanWindowAddField (in->window, inaxes, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

	/* Just copy window */
	ObitDConCleanWindowReplaceField (in->window, ifield+1, in->window, newField, err);
	/* Make sure at least 5 pixel radius first box */
	if (ObitDConCleanWindowInfo (in->window, newField, 1, &otype, &owin, err)) {
	  owin[0] = MAX (5, owin[0]);
	  ObitDConCleanWindowUpdate(in->window, newField, 1, otype, owin, err);
	}
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

	/* Add field to Clean */
	inClass->ObitDConCleanVisAddField (in, uvdata, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Add field to SkyModel */
	ObitSkyModelAddField (in->skyModel, err);
	if (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
	
	/* Indicate association */
	in->mosaic->isAuto[ifield]  = mosaic->numberImages;
	in->mosaic->isShift[mosaic->numberImages-1] = ifield+1;

	/* Set shift on uv data */
	ObitInfoListGetP(uvdata->info, "xShift", &type, dim, (gpointer*)&farray);
	farray[newField-1]  = mosaic->images[newField-1]->myDesc->xshift;
	ObitInfoListAlwaysPut(uvdata->info, "xShift", type, dim, farray);
	ObitInfoListGetP(uvdata->info, "yShift", &type, dim, (gpointer*)&farray);
	farray[newField-1]  = mosaic->images[newField-1]->myDesc->yshift;
	ObitInfoListAlwaysPut(uvdata->info, "yShift", type, dim, farray);
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
  gboolean *btemp;
  ofloat *ftemp;
  const ObitDConCleanVisClassInfo *inClass;
  gchar *routine = "ObitDConCleanVisAddField";

  /* Error checks */
  if (err->error) return;  /* previous error? */
  g_assert(ObitDConCleanVisIsA(in));

  /* Clear BeamPatches array */
  inClass     = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* clean class structure */
  inClass->KillBeamPatches (in);

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
  for (i=0; i<oldField; i++) {ftemp[i] = in->gain[i];} ftemp[i] = 0.0; 
  in->gain = ObitMemFree(in->gain);
  in->gain = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean minFlux");
  for (i=0; i<oldField; i++) {ftemp[i] = in->minFlux[i];} ftemp[i] = 0.0; 
  in->minFlux = ObitMemFree(in->minFlux);
  in->minFlux = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean factor");
  for (i=0; i<oldField; i++) {ftemp[i] = in->factor[i];} ftemp[i] = 0.0; 
  in->factor = ObitMemFree(in->factor);
  in->factor = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean quality");
  for (i=0; i<oldField; i++) {ftemp[i] = in->quality[i];} ftemp[i] = 0.0; 
  in->quality = ObitMemFree(in->quality);
  in->quality = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean cleanable");
  for (i=0; i<oldField; i++) {ftemp[i] = in->cleanable[i];} ftemp[i] = 0.0; 
  in->cleanable = ObitMemFree(in->cleanable);
  in->cleanable = ftemp;
  btemp = ObitMemAlloc0Name(newField*sizeof(gboolean),"Clean fresh");
  for (i=0; i<oldField; i++) {btemp[i] = in->fresh[i];} btemp[i] = FALSE; 
  in->fresh = ObitMemFree(in->fresh);
  in->fresh = btemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean max res");
  for (i=0; i<oldField; i++) {ftemp[i] = in->maxAbsRes[i];} ftemp[i] = 0.0; 
  in->maxAbsRes = ObitMemFree(in->maxAbsRes);
  in->maxAbsRes = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Clean avg res");
  for (i=0; i<oldField; i++) {ftemp[i] = in->avgRes[i];} ftemp[i] = 0.0; 
  in->avgRes = ObitMemFree(in->avgRes);
  in->avgRes = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Image RMS");
  for (i=0; i<oldField; i++) {ftemp[i] = in->imgRMS[i];} ftemp[i] = 0.0; 
  in->imgRMS = ObitMemFree(in->imgRMS);
  in->imgRMS = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Image Peak/RMS");
  for (i=0; i<oldField; i++) {ftemp[i] = in->imgPeakRMS[i];} ftemp[i] = 0.0; 
  in->imgPeakRMS = ObitMemFree(in->imgPeakRMS);
  in->imgPeakRMS = ftemp;
  ftemp = ObitMemAlloc0Name(newField*sizeof(ofloat),"Beam Peak/RMS");
  for (i=0; i<oldField; i++) {ftemp[i] = in->beamPeakRMS[i];} ftemp[i] = 0.0; 
  in->beamPeakRMS = ObitMemFree(in->beamPeakRMS);
  in->beamPeakRMS = ftemp;
  itemp = ObitMemAlloc0Name((newField+3)*sizeof(olong),"currentFields");
  for (i=0; i<oldField; i++) {itemp[i] = in->currentFields[i];} itemp[i] = 0; 
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
  ofloat tol, autoCenFlux;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong   nfield, ifield, itemp, nccpos;
  olong  CCVer, highVer;
  ofloat tmax, xcenter, ycenter, xoff, yoff, radius, cells[2], pixel[2];
  ofloat RAShift, DecShift, *farray, dx, dy;
  odouble pos[2], RAPnt, DecPnt;
  gboolean autoCen, want;
  gchar *routine = "ObitDConCleanVisRecenter";

  /* Error checks */
  if (err->error)     return redo;  /* previous error? */
  g_assert(ObitDConCleanVisIsA(in));

  /* Only Stokes I */
  if (in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocs]>1.0) return redo;

  /* Number of fields */
  nfield = mosaic->numberImages;

   /* Get cellsize */
  cells[0] =  fabs(mosaic->xCells); cells[1] = fabs(mosaic->yCells);

  /* Consider components within 2.5  cells  */
  radius = 2.5 * cells[0];

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
    mosaic->images[ifield]->extBuffer = TRUE;  /* Don't need buffer here */
    mosaic->images[ifield]->image = ObitFArrayUnref(mosaic->images[ifield]->image);
    ObitImageOpen (mosaic->images[ifield], OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Make sure there is already a CC table */
    highVer =  ObitTableListGetHigh (mosaic->images[ifield]->tableList, "AIPS CC");
    if (highVer<=0) goto closeit;

    /* Make temporary CC table object including CCs from any overlapping facets */
    CCTab = ObitImageMosaicCombineCC (mosaic, ifield+1, CCVer, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);

    /* Is this an autoCenter field with CLEAN components? */
    autoCenFlux = 1.0e20;
    ObitInfoListGetTest(mosaic->images[ifield]->info, "autoCenFlux",  &type, dim, &autoCenFlux);
    if (CCTab->myDesc->nrow>0) {
      
      /* Check Table */
      nccpos = CCTab->myDesc->nrow;
      ObitImageMosaicMaxCC (CCTab, nccpos, radius, &tmax, &xcenter, &ycenter, 
			    &xoff, &yoff, err);
      if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
      
      /* See if shift needed - if more then 2 cells this is probably a mistake */
      imDesc = mosaic->images[ifield]->myDesc;
      if (imDesc->do3D) {
	dx = xcenter;
	dy = ycenter;
      } else { /* 2D */
	dx = xcenter + imDesc->xPxOff*imDesc->cdelt[0];
	dy = ycenter + imDesc->yPxOff*imDesc->cdelt[1];
      }
      want = ((fabs(dx+xoff)>tol*cells[0]) || (fabs(dy+yoff)>tol*cells[1]));
      want = want && (fabs(dx)<2.0) && (fabs(dy)<2.0);
      want = want && (tmax>=autoCenFlux*0.5);   /* Might have been peeled */
      if (want) {
	
	/* Position of peak */
	pixel[0] = imDesc->crpix[0] + ((xcenter+xoff) / imDesc->cdelt[0]);
	pixel[1] = imDesc->crpix[1] + ((ycenter+yoff) / imDesc->cdelt[1]);
	/* DEBUG HACK -0.5* */
 	pixel[0] = imDesc->crpix[0] + ((xcenter-0.5*xoff) / imDesc->cdelt[0]);
	pixel[1] = imDesc->crpix[1] + ((ycenter-0.5*yoff) / imDesc->cdelt[1]);
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
	if (err->prtLv>1) {
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
    mosaic->images[ifield]->extBuffer = FALSE;  /* May need buffer later */
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
  if (err->prtLv>1) {
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
    mosaic->images[ifield]->extBuffer = TRUE;  /* Don't need buffer here */
    mosaic->images[ifield]->image = ObitFArrayUnref(mosaic->images[ifield]->image);
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
    mosaic->images[ifield]->extBuffer = FALSE;  /* May need buffer later */
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, redo);
  } /* end loop over fields */
  
  /* Tell if found something */
  if (err->prtLv>1) {
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
 * Facets might have a window added iff 1) it is the facet with the highest 
 * in window pixel, or 2) it has a cleanable, out of window pixel within 30% of the peak.
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
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  i, j, best, field, blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat PeakIn, PeakOut, RMS, bestPeak, oldAutoWinFlux, *FieldPeak=NULL;
  olong PeakInPos[2] = {0,0};
  gboolean doAbs, isMR, found, tnewWin;
  gchar *routine = "ObitDConCleanVisAutoWindow";

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanIsA(in));

  /* Set output to full image, plane at a time */
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];

  /* Is this a multiresolution image (to allow neg regions) */
  isMR = (in->mosaic->numBeamTapes>1);

  /* find best field */
  FieldPeak = g_malloc0(in->nfield*sizeof(ofloat));
  for (field=0; field<in->nfield; field++) FieldPeak[field]=0.0;

  best = 0;
  bestPeak = -1.0e20;
  for (field=0; field<in->nfield; field++) {
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
      
      ObitImageOpen (image, OBIT_IO_ReadOnly, err);
      ObitImageRead (image, image->image->array, err);
      ObitImageClose (image, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
      usePixels = image->image;  /* Pointer to image buffer */
      
    } else { /* Use array passed */
      usePixels = ObitFArrayRef(pixarray[field]);
      image = in->mosaic->images[0];  /* Need one */
    }
    
    /* Allow negative for Stokes other than I or multiresolution */
    doAbs = isMR || (fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1);
  
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
    if (err->prtLv>4) {
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
    in->imgRMS[fields[field]-1]     = RMS;
  } /* end loop gathering statistics */

  /* Loop over potential fields */
  for (j=0; j<in->numCurrentField; j++) {
    if (fields[j]==0) break;  /* List terminated? */
    /* Use fields within 30% of best Peak in window - or best */
    if ((FieldPeak[j]<0.7*bestPeak) && (j!=best)) continue;
    field = fields[j];

    /* Use passed pixel array or get image */
    if ((pixarray==NULL) || (pixarray[j]==NULL)) {   
      image = in->mosaic->images[field-1];
      
      /* Set input to full image, plane at a time */
      dim[0] = IM_MAXDIM;
      ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
      ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc); 
      dim[0] = 1;
      ObitInfoListAlwaysPut (image->info, "IOBy", OBIT_long, dim, &IOsize);
      
      ObitImageOpen (image, OBIT_IO_ReadOnly, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
      
      ObitImageRead (image, image->image->array, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
      
      ObitImageClose (image, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
      usePixels = image->image;  /* Pointer to image buffer */
      
    } else { /* Use array passed */
      usePixels = ObitFArrayRef(pixarray[j]);
      image = in->mosaic->images[0];  /* Need one */
    }
    
    /* Allow negative for Stokes other than I */
    doAbs = isMR || (fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1);
    
    /* Get field info - set new box if needed */
    tnewWin = ObitDConCleanWindowAutoWindow (in->window, field, usePixels,
					     doAbs,
					     &PeakIn, &PeakInPos[0], &PeakOut, 
					     &RMS, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
    newWin = newWin || tnewWin;
    
    /* If window added, redo statistics */
    if (tnewWin) {
      ObitDConCleanWindowStats (in->window, field, usePixels,
				doAbs,
				&PeakIn, &PeakInPos[0], &PeakOut, 
				&RMS, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, newWin);
      FieldPeak[j] = PeakOut;
    }

    /* Free Image array? */
    if (usePixels==image->image) 
      image->image = ObitFArrayUnref(image->image);
    else
      usePixels = ObitFArrayUnref(usePixels);
  } /* end loop over fields */
  
  /* Determine max. cleanable pixel value not in a window */
  for (field=0; field<in->nfield; field++) {
    if (fields[field]==0) break;  /* List terminated? */
    if (field!=best) PeakOut = MAX (PeakOut,  FieldPeak[field]);
  }
  g_free(FieldPeak);

  /* Set flux limit for next cycle */
  oldAutoWinFlux = in->autoWinFlux;
  if (oldAutoWinFlux<0.0) oldAutoWinFlux = 1.0e20;
  in->autoWinFlux = MAX (PeakOut-5.0*RMS, 0.25*RMS);
  in->autoWinFlux = MIN (in->autoWinFlux, 0.85*PeakOut);

  /* Don't let the autoWinFlux go much below cleanable maxima in other fields */
  for (i=0; i<in->nfield; i++) {
 
     /* Test BeamTaper */
    if (((in->mosaic->BeamTaper[i]<in->minBeamTaper) || 
	 (in->mosaic->BeamTaper[i]>in->maxBeamTaper))) continue;

    /* Is this one in fields? */
    found = FALSE;
    for (j=0; j<in->nfield; j++) {
      if (fields[j]==0) break;
      if ((fields[j]-1)==i) found = TRUE;
      /* Look for autoCenter/shifted pairs */
      if (in->mosaic->isAuto[i]  == (fields[j])) found = TRUE;
      if (in->mosaic->isShift[i] == (fields[best])) found = TRUE;
      if (in->mosaic->isShift[fields[j]-1] == (i+1)) found = TRUE;
      if (found) break;
    }
    if (!found && (in->imgPeakRMS[i]>5.0) && in->quality[i]>0.0) {
      /* Only to 70% minus 3 sigma */
      in->autoWinFlux = MAX (in->autoWinFlux, 0.7*in->cleanable[i]-3.0*RMS);
      /* Diagnostic */
      if (err->prtLv>4) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Nonclean field %d cleanable %f Peak/RMS %g quality %g maxAbsRes %g",
		       i+1, in->cleanable[i], in->imgPeakRMS[i], in->quality[i],
		       in->maxAbsRes[i]);
      }
    }
  }
  /* Don't set too low */
  in->autoWinFlux = MAX (in->autoWinFlux, 0.25*RMS);

  /* Diagnostic */
  if (err->prtLv>4) {
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
    in->minFlux[fields[best]-1] = MAX (in->minFlux[fields[best]-1], 0.5*RMS);
     /* Tell user */
     if (err->prtLv>0) {
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
  gboolean newWin = FALSE, isLast = FALSE, allDone=TRUE, OK;
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
  for (ifld=0; ifld<=in->nfield; ifld++) {
    
    field = in->currentFields[ifld];
    if (field<=0) break;  /* List terminates */

    /* Is this field OK to clean? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[field-1]>1.0 ||
	in->cleanable[field-1]/in->imgRMS[field-1]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[field-1]))));
    else
      OK = in->imgPeakRMS[field-1]>1.0;

    /* Is this field finished? */  
    if ((in->Pixels->complCode!=OBIT_CompReasonNiter) &&
	(in->Pixels->complCode!=OBIT_CompReasonMinFlux)) {
      allDone = allDone && 
	((in->cleanable[field-1] <= in->minFlux[field-1]) || !OK);
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
      if (err->prtLv>2) {
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
  if (in->minFluxLoad<=0.0)
    in->minFluxLoad = MAX (in->autoWinFlux, in->minFluxLoad);
  in->numberSkip = 1 + ObitDConCleanPxHistNumber(in->PixelHist, minFlux, err) / 
    MAX (1, in->maxPixel);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* See if SDI CLEAN called for */
 doSDI:
  if (in->SDIGain>0.0) {
    fract = (ofloat)in->PixelHist->hist[in->PixelHist->ncell/2] / MAX(1, in->PixelHist->hist[0]);
    if (err->prtLv>1)
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
  }

  /* Make sure at least autoWinFlux */
  if (in->minFluxLoad<=0.0)
    in->minFluxLoad = MAX (in->autoWinFlux, in->minFluxLoad);

  /* Give warning if skipping */
  if ((in->numberSkip>=1) && (err->prtLv>0)) 
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
  theClass->ObitDConCleanVisDefWindow = (ObitDConCleanVisDefWindowFP)ObitDConCleanVisDefWindow;
  theClass->ObitDConCleanAutoWindow = 
    (ObitDConCleanAutoWindowFP)ObitDConCleanVisAutoWindow;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanVisPixelStats;
  theClass->ObitDConCleanVisCleanable = (ObitDConCleanVisCleanableFP)ObitDConCleanVisCleanable;

  /* Private functions for derived classes */
  theClass->MakeResiduals   = (MakeResidualsFP)MakeResiduals;
  theClass->MakeAllResiduals= (MakeAllResidualsFP)MakeAllResiduals;
  theClass->SubNewCCs       = (SubNewCCsFP)SubNewCCs;
  theClass->NewPxList       = (NewPxListFP)NewPxList;
  theClass->NewPxArray      = (NewPxArrayFP)NewPxArray;
  theClass->KillPxArray     = (KillPxArrayFP)KillPxArray;
  theClass->KillBeamPatches = (KillBeamPatchesFP)KillBeamPatches;
  theClass->PickNext2D      = (PickNext2DFP)PickNext2D;
  theClass->PickNext3D      = (PickNext3DFP)PickNext3D;
  theClass->WhosBest        = (WhosBestFP)WhosBest;
  theClass->WhosBest2D      = (WhosBest2DFP)WhosBest2D;
  theClass->OrderImage      = (OrderImageFP)OrderImage;
  theClass->OrderClean      = (OrderCleanFP)OrderClean;
  theClass->SelectTaper     = (SelectTaperFP)SelectTaper;
  theClass->ResetSkyModel   = (ResetSkyModelFP)VisResetSkyModel;
  theClass->ResetPixelList  = (ResetPixelListFP)VisResetPixelList;
  theClass->FindPeak        = (FindPeakFP)FindPeak;
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
  olong i;
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
  in->mosaic2  = NULL;
  in->modelMode = OBIT_SkyModel_Fastest;
  in->doRestore = TRUE;
  in->doXRestore= TRUE;
  in->doFlatten = TRUE;
  in->doWeight  = TRUE;
  in->doBeam    = TRUE;
  in->quality   = NULL;
  in->cleanable = NULL;
  in->fresh     = NULL;
  in->display   = NULL;
  in->SDIdata   = NULL;
  in->peakFlux  = -1000.0;
  in->reuseFlux = -1.0;
  in->autoCen   =  1.0e20;
  in->maxBeamTaper = 0.0;
  in->minBeamTaper = 0.0;
  for (i=0; i<10; i++) in->MResKnob[i] = 0.0;
  in->outaTime     = FALSE;
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
  in->mosaic2   = ObitImageMosaicUnref(in->mosaic2);
  in->quality   = ObitMemFree(in->quality);
  in->cleanable = ObitMemFree(in->cleanable);
  in->fresh     = ObitMemFree(in->fresh);
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
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, ifld, jfld, field=1, *stale=NULL;
  ofloat autoCenFlux=0.0;
  gchar *routine = "MakeResidual";

  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */

  /* Are residuals fresh? */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Make copy of only list of stale images */
  stale = g_malloc0((in->nfield+1)*sizeof(olong));
  ifld = jfld = 0;
  while(fields[ifld]>0) {
    if (!in->fresh[fields[ifld]-1]) stale[jfld++] = fields[ifld];
    ifld++;
  }
  stale[jfld] = 0;

  /* If all stale then add all to list - NO, this isn't right
  if (stale[0]<=0) {
   for (ifld=0; ifld<in->nfield; ifld++) {
     stale[ifld] = ifld+1;
   }
   stale[ifld] = 0;
  }*/

  /* If all fresh, bail */
  if (stale[0]<=0) {g_free(stale); return;}

  /* Make residual images for stale fields */
  imgClass->ObitUVImagerImage (in->imager, stale, doWeight, doBeam, doFlatten, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Statistics */
  for (ifld=0; ifld<in->nfield; ifld++) {
    field = stale[ifld];
    if (field<=0) break;
    /* Get statistics  for image */
    inClass->ObitDConCleanImageStats ((ObitDConClean*)in, field, FALSE, err);
    /* Need Beam statistics? */
    if (doBeam || in->beamPeakRMS[field-1]<=0.0)
      inClass->ObitDConCleanImageStats ((ObitDConClean*)in, field, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Quality measure */
    in->quality[field-1] = inClass->ObitDConCleanVisQuality(in, field, err);

    /* Max cleanable flux */
    in->cleanable[field-1] = inClass->ObitDConCleanVisCleanable(in, field, NULL, err);
    in->fresh[field-1]    = TRUE;
  } /* end statistics loop */
  
  if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* For any shifted fields get statistics */
  for (ifld=0; ifld<in->nfield; ifld++) {
    /* get autoCenFlux if autoCentering */
    if (in->mosaic->isAuto[ifld]>0)
      ObitInfoListGetTest(in->mosaic->images[ifld]->info, "autoCenFlux", &type, 
			  dim, &autoCenFlux);
    if (in->mosaic->isShift[ifld] > 0) {

      /* Has this field reached the min flux? maxAbsRes<minFlux */
      if (in->maxAbsRes[ifld]<=in->minFlux[ifld]) continue;
      jfld = in->mosaic->isShift[ifld]-1;
      /* If the corresponding autoCenter field just remade (as well as the shifted
	 version - determine the statistics from the image */
      found = FALSE;
      for (i=0; i<in->nfield; i++) {
	if (fields[i]==0) break;
	if ((jfld+1)==stale[i]) {found = TRUE; break;}
      }

      if (found) {
	/* Get statistics  for image */
	inClass->ObitDConCleanImageStats ((ObitDConClean*)in, ifld+1, FALSE, err);
	/* Quality measure */
	in->quality[ifld] = inClass->ObitDConCleanVisQuality((ObitDConCleanVis*)in, ifld+1, err);
	/* Max cleanable flux */
	in->cleanable[ifld]   = inClass->ObitDConCleanVisCleanable(in, ifld+1, NULL, err);
	in->beamPeakRMS[ifld] = in->beamPeakRMS[jfld];
	in->fresh[ifld]       = TRUE;
     } else { /* Not remade - Use residual */
	in->maxAbsRes[ifld] = in->maxAbsRes[jfld] = in->Pixels->maxResid;
	in->quality[ifld]   = in->quality[jfld]   = in->Pixels->maxResid;
	in->cleanable[ifld] = in->cleanable[jfld] = in->Pixels->maxResid;
	in->beamPeakRMS[ifld] = in->beamPeakRMS[jfld];
	in->fresh[ifld]       = in->fresh[jfld];
      }
    }
  }

  /* For autoCenter fields with peaks below the threshold, set stats to shifted */
  for (ifld=0; ifld<in->nfield; ifld++) {
    jfld = in->mosaic->isAuto[ifld]-1;
    if ((jfld>=0) && (jfld<in->mosaic->numberImages) && 
	(in->cleanable[ifld]<(0.1*autoCenFlux))) {
      in->cleanable[ifld] = in->cleanable[jfld];
      in->maxAbsRes[ifld] = in->maxAbsRes[jfld];
      in->quality[ifld]   = in->quality[jfld];
    }
  }

  if(stale) g_free(stale);
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

  if (err->error) return;   /* existing error */
  
  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */

  /* Turn off things not needed */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Copy prtLv to in->mosaic->info */
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (in->mosaic->info, "prtLv", OBIT_long, dim, &err->prtLv);

  /* Parallel Image images without needing beam */
  imgClass->ObitUVImagerImage (in->imager, fields,  doWeight, doBeam, doFlatten, err);
  
  /* Loop over fields getting statistics for Image and Beam */
  for (i=0; i<in->nfield; i++) {
    inClass->ObitDConCleanImageStats ((ObitDConClean*)in, i+1, FALSE, err);
    if (doBeam)
      inClass->ObitDConCleanImageStats ((ObitDConClean*)in, i+1, TRUE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Quality measure */
    in->quality[i] = inClass->ObitDConCleanVisQuality(in, i+1, err);
    /* Max cleanable flux */
    in->cleanable[i] = inClass->ObitDConCleanVisCleanable(in, i+1, NULL, err);
    in->fresh[i]     = TRUE;  /* Freshly made image */
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end loop over field */

  /* For any shifted fields get statistics */
  for (ifld=0; ifld<in->mosaic->numberImages; ifld++) {
    if (in->mosaic->isShift[ifld] > 0) {
      jfld = in->mosaic->isShift[ifld]-1;
      /* Get statistics  for image */
      inClass->ObitDConCleanImageStats ((ObitDConClean*)in, ifld+1, FALSE, err);
      /* Quality measure */
      in->quality[ifld] = inClass->ObitDConCleanVisQuality(in, ifld+1, err);
      /* Max cleanable flux */
      in->cleanable[ifld]   = inClass->ObitDConCleanVisCleanable(in, ifld+1, NULL, err);
      in->beamPeakRMS[ifld] = in->beamPeakRMS[jfld];
      in->fresh[ifld]       = in->fresh[jfld];
    }
  }

} /* end MakeAllResiduals */

/**
 * Zero all residual images
 * \param in     The Clean object
 * \param err    Obit error stack object.
 */
static void  ClearAllResiduals (ObitDConCleanVis *in, ObitErr *err)
{
  olong i, ip, np;
  ObitImage *image=NULL;
  olong plane[5] = {1,1,1,1,1};
  gchar *routine = "ClearAllResiduals";

  if (err->error) return;   /* existing error */
  
  /* Loop over fields */
  for (i=0; i<in->nfield; i++) {
    image = in->mosaic->images[i];
    ObitImageOpen (image, OBIT_IO_ReadWrite, err);
    /* Loop over planes */
    np = image->myDesc->inaxes[2];
    for (ip=1; ip<=np; ip++) {
      plane[0] = ip;
      ObitImageGetPlane (image, image->image->array, plane, err);
      if (err->error) Obit_traceback_msg (err, routine, image->name);
      ObitFArrayFill(image->image, 0.0);  /* Zero */
      ObitImagePutPlane (image, image->image->array,plane,  err);
      if (err->error) Obit_traceback_msg (err, routine, image->name);
    }
    ObitImageClose (image, err);
    /* Free Image array */
    image->image = ObitFArrayUnref(image->image);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
  } /* end loop over field */

} /* end ClearAllResiduals */

/**
 * Find current estimates best and second field  3D
 * Ignore fields with imgPeakRMS<4.0 if autoWindow, else imgPeakRMS<1.0
 * \param in     The Clean object
 * \param best   [out] 0-rel best field
 * \param second [out] 0-rel second best field
 */
static void WhosBest (ObitDConCleanVis *in, olong *bbest, olong *ssecond)
{
  ofloat testBest, testSecond;
  olong i, best, second;
  gboolean OK;
  
  best = second = -1; 
  testBest = testSecond = -1.0e20;
  for (i=0; i<in->nfield; i++) {

    /* Test BeamTaper */
    if (((in->mosaic->BeamTaper[i]<in->minBeamTaper) || 
	 (in->mosaic->BeamTaper[i]>in->maxBeamTaper))) continue;
    
    /* Is this field OK to CLEAN? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[i]>1.0 ||
	in->cleanable[i]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
    else
      OK = in->imgPeakRMS[i]>1.0;

    OK = OK && (in->cleanable[i]>in->minFlux[i]);

     /* Test BeamTaper */
    OK = OK && (((in->mosaic->BeamTaper[i]>=in->minBeamTaper) && 
		 (in->mosaic->BeamTaper[i]<=in->maxBeamTaper)));

    if (!OK) continue;  /*ignore if too low SNR or done */

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
 * Ignore fields with imgPeakRMS<4.0 if autoWindow, else imgPeakRMS<1.0
 * \param in          The Clean object
 * \param autoCenFlux Cutoff level for autocentering,  0=> no autoCenter
 * \param bbest       [out] 0-rel best field
 */
static void WhosBest2D (ObitDConCleanVis *in, ofloat autoCenFlux, 
			olong *bbest)
{
  ofloat testBest;
  olong i, best, isShift;
  gboolean OK, isAuto=FALSE;

  /* Find absolute best */
  best = -1; 
  testBest = -1.0e20;
  for (i=0; i<in->nfield; i++) {

    /* Ignore shift fields if can still use autoCentered fields */
    isShift = in->mosaic->isShift[i];
    if ((isShift>0) && (autoCenFlux>0.0) &&
	(in->cleanable[isShift-1]>(0.1*autoCenFlux))) continue;
    /* Redundant if ((isShift>0) && (in->maxAbsRes[isShift-1]>(0.1*autoCenFlux))) continue;*/

    /* Ignore autoCenter fields below threshold */
    if ((in->mosaic->isAuto[i]>0) && 
	(in->maxAbsRes[i]<0.1*autoCenFlux)) continue;

     /* Test BeamTaper */
    if (((in->mosaic->BeamTaper[i]<in->minBeamTaper) || 
	 (in->mosaic->BeamTaper[i]>in->maxBeamTaper))) continue;

    /* Is this field OK to CLEAN? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[i]>1.0 ||
	in->cleanable[i]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
    else
      OK = in->imgPeakRMS[i]>1.0;
    
    if (!OK) continue;  /* Ignore if too low SNR */
    
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
    
    /* Is this field OK to CLEAN? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[i]>1.0 ||
	in->cleanable[i]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
    else
      OK = in->imgPeakRMS[i]>1.0;
    
     /* Test BeamTaper */
    OK = OK && (((in->mosaic->BeamTaper[i]>=in->minBeamTaper) && 
		 (in->mosaic->BeamTaper[i]<=in->maxBeamTaper)));

   /* Ignore if SNR too low */
    if (!OK) continue;
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
 * Only includes images with a quality within 30% of the best or 
 * with cleanables at least that of the best.
 * Uses autoShift only for flux densities above 0.1 x autoCenFlux
 * Fields with quality=0.0 (reached CLEAN goal) are ignored.
 * If the best is an autoWindow, only other autoWindow fields are considered
 * Includes allowed range of BeamTaper and all tapers of selected images
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
  olong i, j, k, n, m, indx, best, myAuto, bQual;
  gboolean OK, done;
  
  /* make sure no shift field has resid > corresponding autoCen */
  if (autoCenFlux>0.0) {
    for (i=0; i<in->nfield; i++) {
      myAuto = in->mosaic->isShift[i]-1;
      if ((myAuto>=0) && (in->cleanable[myAuto]<in->cleanable[i]))
	{
	  in->cleanable[i] = in->cleanable[myAuto];
	  in->quality[i]   = in->quality[myAuto];
	  in->maxAbsRes[i] = in->maxAbsRes[myAuto];
	}
    } /* end loop */
  } /* end if autoCenter */
  
  /* Make temporary array of quality factors */
  tmpQual = g_malloc0(in->mosaic->numberImages*sizeof(ofloat));

  /* Copy quality factors to temporary array */
  n = in->mosaic->numberImages;
  bestQual = -1.0e18;
  bQual    = 0;
  done = in->mosaic->numBeamTapes<=1; /* test for one taper case */
  for (i=0; i<n; i++) {
    /* Is this field OK to image? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[i]>1.0 ||
	in->cleanable[i]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
    else
      OK = in->imgPeakRMS[i]>1.0;

    /* Ignore autoCen facets below threshold */
    if ((in->mosaic->isAuto[i]>=0) && (in->maxAbsRes[i]<0.1*autoCenFlux)) 
      OK = FALSE;
    
    /* Ignore shifted facets above threshold */
    if ((in->mosaic->isShift[i]>=0) && (in->maxAbsRes[i]>0.1*autoCenFlux)) 
      OK = FALSE;
    
    /* Test BeamTaper */
    OK = OK && ((in->mosaic->BeamTaper[i]>=in->minBeamTaper) && 
		(in->mosaic->BeamTaper[i]<=in->maxBeamTaper));

    /* Is the CLEAN in this field done? */
    if (OK) 
      done = done && (in->quality[i]<=in->minFlux[i]);

    if (!fresh[i]) tmpQual[i] = in->quality[i];
    else tmpQual[i] = -1.0e20;
    if ((tmpQual[i]>bestQual) && (OK)) {
      bestQual = tmpQual[i];       
      bQual    = i;          /* Best quality */
    }
  }

  /* Anything to do? */
  if (done) {fldList[0]=0; goto cleanup;}

  /* Loop finding and dropping best until all gone */
  indx = 0; 
  while (1) {   /* Sort loop */
    maxQual = -1.0e18;
    best    = -1;
    done = TRUE;
    for (i=0; i<n; i++) {
      /* Ignore autoCen facets below threshold 
      if ((in->mosaic->isAuto[i]>=0) && (in->maxAbsRes[i]<0.1*autoCenFlux)) 
	continue;*/
    
      /* Ignore shifted facets above threshold */
      if ((in->mosaic->isShift[i]>=0) && (in->maxAbsRes[i]>0.1*autoCenFlux)) 
	continue;
    
      /* Test BeamTaper */
      if ((in->mosaic->BeamTaper[i]<in->minBeamTaper) || 
	  (in->mosaic->BeamTaper[i]>in->maxBeamTaper)) continue;

      if (tmpQual[i]>maxQual) { /* Best so far */
	maxQual = tmpQual[i];
	best = i;
	done = FALSE;
      }
    }
    /* Find anything? */
    if (done) break;
    
    /* Test BeamTaper */
    OK = ((in->mosaic->BeamTaper[best]>=in->minBeamTaper) && 
	  (in->mosaic->BeamTaper[best]<=in->maxBeamTaper));
    
    /* Save value if within 30% of best, or larger cleanable, same taper,  dummy tmpQual */
    if (((maxQual>=0.7*bestQual) || (in->cleanable[best]>=in->cleanable[bQual])) 
	&& (!fresh[best]) && OK) {
      fldList[indx++] = best+1;
    }
    if ((best>=0) && (best<n)) tmpQual[best]   = -1.0e20;
  }  /* end sort loop  */

  /* If multiple tapers add all matching selected fields */
  if (in->mosaic->numBeamTapes>1) {
    m = indx;
    for (j=0; j<m; j++) {
      k = fldList[j] - 1;  /* Zero rel field of selected facet */
      for (i=0; i<n; i++) {  /* Loop throu mosaic */
	if ((k!=i) && /* Not same */
	    (in->mosaic->FacetNo[i]==in->mosaic->FacetNo[k]) &&
	    (!fresh[i])) {
	  fldList[indx++] = i+1;
	}
      } /* End loop over fields */
    } /* end loop over selected images */
  } /* end add other taperings of selected images */
  
  fldList[indx++] = 0; /* terminate */
  
  /* Cleanup */
 cleanup:
  if (tmpQual) g_free(tmpQual);
  
} /* end OrderImage */

/**
 * Find priority ordering for cleaning residual images
 * Ranking by quality, using only "fresh" and no shifted images
 * Only includes images with a quality within 50% of the bes,
 * or a cleanable >= the that of the best.
 * If the best is an autoWindow (>0.1autoCenFlux), only other 
 * autoWindow fields are considered
 * Ignore fields with imgPeakRMS<4.0 if autoWindow, else imgPeakRMS<1.0
 * Uses autoShift only for flux densities above 0.1 x autoCenFlux
 * Includes allowed range of BeamTaper
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
  ofloat *tmpQual=NULL, maxQual, bestQual, testBest, bestClean;
  olong i, n, indx, myAuto, best, isShift, itemp;
  gboolean OK, done, isAuto=FALSE;

  /* make sure no shift field has resid > corresponding autoCen */
  if (autoCenFlux>0.0) {
    for (i=0; i<in->nfield; i++) {
      myAuto = in->mosaic->isShift[i]-1;
      if ((myAuto>0) && (in->cleanable[myAuto]<in->cleanable[i]))
	{
	  in->cleanable[i] = in->cleanable[myAuto];
	  in->quality[i]   = in->quality[myAuto];
	  in->maxAbsRes[i] = in->maxAbsRes[myAuto];
	}
    } /* end loop */
  } /* end if autoCenter */
  
  /* Find absolute best */
  best = -1; 
  testBest  = -1.0e20;
  bestClean = 0.0;  /* Best cleanable */
  for (i=0; i<in->nfield; i++) {
    /* Is this field OK to CLEAN? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[i]>1.0 ||
	in->cleanable[i]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
    else
      OK = in->imgPeakRMS[i]>1.0;

    /* Test BeamTaper */
    OK = OK && ((in->mosaic->BeamTaper[i]>=in->minBeamTaper) && 
		(in->mosaic->BeamTaper[i]<=in->maxBeamTaper));

    /* Ignore shifted fields with cleanable > 0.1*autoCenFlux */
    isShift = in->mosaic->isShift[i];
    if ((isShift>0) && (in->cleanable[isShift-1]>(0.1*autoCenFlux))) 
      OK = FALSE; 
 
    /* Ignore autoCenter fields with cleanable < 0.1*autoCenFlux */
    if ((in->mosaic->isAuto[i]>=0) && (in->cleanable[i]<(0.1*autoCenFlux))) 
      OK = FALSE;
    
    if (!OK) continue;  /* Ignore not OK */

    if (in->quality[i]>testBest) {
      testBest  = in->quality[i];
      bestClean = in->cleanable[i];  /* Cleanable of best */
      best = i;
    }
  } /* end loop finding current best */

  /* Make temporary array of quality factors */
  tmpQual = g_malloc0(in->mosaic->numberImages*sizeof(ofloat));

  /* Copy quality factors to temporary array */
  n = in->mosaic->numberImages;
  bestQual = -1.0e18;
  for (i=0; i<n; i++) {
    /* Is this field OK to CLEAN? */
    if (in->autoWindow)
      OK = in->imgPeakRMS[i]>1.0 ||
	in->cleanable[i]/in->imgRMS[0]>MAX (4.0,(0.1*(MIN(40.0,in->beamPeakRMS[i]))));
    else
      OK = in->imgPeakRMS[i]>1.0;

    /* Ignore shifted fields with cleanable > 0.1*autoCenFlux */
    isShift = in->mosaic->isShift[i];
    if ((isShift>0) && (in->cleanable[isShift-1]>(0.1*autoCenFlux))) 
      OK = FALSE;

    /* Ignore autoCenter fields with cleanable < (0.1*autoCenFlux)  */
    if ((in->mosaic->isAuto[i]>=0) && (in->cleanable[i]<(0.1*autoCenFlux))) 
      OK = FALSE;
    
    /* Test BeamTaper */
    OK = OK && ((in->mosaic->BeamTaper[i]>=in->minBeamTaper) && 
		(in->mosaic->BeamTaper[i]<=in->maxBeamTaper));

    if (fresh[i] && OK) tmpQual[i] = in->quality[i];
    else tmpQual[i] = -1.0e20;
    if (tmpQual[i]>bestQual) {
      bestQual = tmpQual[i];                 /* Best quality */
      isAuto   = (in->mosaic->isAuto[i]>=0)
	&& (in->maxAbsRes[i]>0.1*autoCenFlux) 
        && (in->quality[i]>=testBest);        /* Is best a strong autoCenter image? */
    }
  } 

  /* Loop finding and dropping best until all gone */
  indx = 0;
  while (1) {    /* Sort loop */
    maxQual = -1.0e18;
    best    = -1;
    done = TRUE;
    for (i=0; i<n; i++) {
      /* Test BeamTaper */
      if ((in->mosaic->BeamTaper[i]<in->minBeamTaper) || 
	  (in->mosaic->BeamTaper[i]>in->maxBeamTaper)) continue;

      /* Ignore shifted fields with cleanable > 0.1*autoCenFlux */
      isShift = in->mosaic->isShift[i];
      if ((isShift>0) && (in->cleanable[isShift-1]>(0.1*autoCenFlux)))
	continue;

      /* Ignore autoCenter fields with cleanable < 0.1*autoCenFlux */
      if ((in->mosaic->isAuto[i]>=0) && (in->cleanable[i]<(0.1*autoCenFlux)))
	continue; 
    
      if ((tmpQual[i]>maxQual) &&                    /* Best so far */
	  (!isAuto || (in->mosaic->isAuto[i]>=0))) { /* Only autoCenter if isAuto */
	maxQual = tmpQual[i];
	best = i;
	done = FALSE;
      }
    }
    
    /* Find anything? */
    if (done) break;
    
    /* Save value if within 50% of best, or cleanable>=best dummy tmpQual */
    if ((maxQual>=0.5*bestQual) || (in->cleanable[best]>=bestClean)){
      /* If this is any autoCenter field and !isAuto, then use the shifted field */
      if (!isAuto && (in->mosaic->isAuto[best]>=0) && (in->mosaic->isAuto[best]<n)) {
	itemp = in->mosaic->isAuto[best];
	if (itemp>n) itemp = in->mosaic->FacetNo[best]+1;
	fldList[indx++] = itemp;
      } else
	fldList[indx++] = best+1;
    }
    if ((best>=0) && (best<n)) tmpQual[best]   = -1.0e20;
  } /* end sort loop */
  
  fldList[indx++] = 0; /* terminate */
  
  /* Cleanup */
  if (tmpQual) g_free(tmpQual);
  
} /* end OrderClean */

/**
 * If multiple tapers are being made, select next
 * Previously CLEANed facets are remade if multiple tapers are in use
 * Ranking by quality, using only "fresh" and no shifted images
 * Only includes images with a quality within 50% of the best
 * If the best is an autoWindow (>0.1 autoCenFlux), only other 
 * autoWindow fields are considered
 * Ignore fields with imgPeakRMS<4.0 if autoWindow, else imgPeakRMS<1.0
 * Uses autoShift only for flux densities above 0.1 x autoCenFlux
 * Includes allowed range of BeamTaper
 * Sets values of in->maxBeamTaper and in->minBeamTaper
 * Sets minFlux for selected taper to the ratio of the beam areas
 * to no taper case.
 * \param in     The Clean object
 * \param fresh  List of flags indicating freshly made
 *               In order of images in mosaic member of in
 * \param err    Obit error stack object.
 * \return TRUE if CLEAN finished
 */
static gboolean SelectTaper (ObitDConCleanVis *in, gboolean *fresh, ObitErr *err)
{
  olong i, j, l, n, m, best, bestTap[30], *fldList=NULL;
  ofloat test, bestTest, bestQuality, bestSNR, maxTape;
  ofloat fact1, fact2, fact3, minT, maxT, cells, minFlux=0.0, gain=0.1, beamrat=1.0;
  gboolean done=FALSE, doBeam=FALSE;
  const ObitDConCleanVisClassInfo *inClass;
  gchar *routine = "SelectTaper";

  /* Anything to do? */
  if (in->mosaic->numBeamTapes<=1) {
    in->minBeamTaper = 0.0;
    in->maxBeamTaper = 0.0;
    return done;
  }

  inClass = (ObitDConCleanVisClassInfo*)in->ClassInfo; /* class structure */
  
  /* Reimage those just CLEANed and not finished and in all resolutions */
  fldList = ObitMemAlloc0((in->nfield+3)*sizeof(olong));
  n = in->mosaic->numberImages;
  i = m = 0;
  /* Add unfinished sources in current List to fldList to reimage */
  while (in->currentFields[i]>0) {
    j = in->currentFields[i]-1;
    if (in->cleanable[j]>in->minFlux[j]) fldList[m++] = in->currentFields[i];
    i++;
  }
  i = 0;
  while (in->currentFields[i]>0) {
    l = in->currentFields[i]-1;
    if (!fresh[l]) {   /* Need to reimage?  */
      /* Find all matching resolutions */
      for (j=0; j<n; j++) {
	if ((j!=l) &&  /* Not the same */
	    (in->mosaic->FacetNo[j]==in->mosaic->FacetNo[l]) &&
	    (!fresh[j]) && (in->cleanable[j]>in->minFlux[j])) {
	  fldList[m++] = j+1;
	}
      }
    }
    i++;
  } /* end finding facets to be remade */
  
  /* Make images */
  fldList[m] = 0;  /* terminate list */
  if (fldList[0]>0) inClass->MakeResiduals (in, fldList, doBeam, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  if (fldList) g_free(fldList);  /* Cleanup */

  /* Control knobs with defaults */
  if (in->MResKnob[0]<=0.0) fact1 = 0.20;
  else                      fact1 = in->MResKnob[0];
  if (in->MResKnob[1]<=0.0) fact2 = 0.33;
  else                      fact2 = in->MResKnob[1];
  if (in->MResKnob[2]<=0.0) fact3 = 0.20;
  else                      fact3 = in->MResKnob[2];

  /* Find best facet in each taper */
  n = in->mosaic->numberImages;
  cells = fabs(in->mosaic->xCells);  /* Cell size */
  for (j=0; j<in->mosaic->numBeamTapes; j++) {
    minT = 0.95*cells*in->mosaic->BeamTapes[j];
    maxT = 1.05*cells*in->mosaic->BeamTapes[j];
    bestTest = -1.0;
    best     = -1;
    for (i=0; i<n; i++) {
      /* Desired taper? */
      if (((in->mosaic->BeamTaper[i]<minT) || 
	   (in->mosaic->BeamTaper[i]>maxT))) continue;
      test = in->quality[i];
      if (test>bestTest) {
	bestTest = test;
	best     = i;
      }
    } /* end loop over facets */
    bestTap[j] = best;
  } /* end loop over tapers */

  /* Show statistics for each resolution */
  if (err->prtLv>3) {
    Obit_log_error(err, OBIT_InfoErr,"Current taper [%f,%f]", 
		   in->minBeamTaper*3600.,  in->maxBeamTaper*3600.);
    for (i=0; i<n; i++) {
      Obit_log_error(err, OBIT_InfoErr,"Facet %d, Taper %f, No. %d, Shift %d", 
		     i+1, in->mosaic->BeamTaper[i]*3600.0, in->mosaic->FacetNo[i],
		     in->mosaic->isShift[i]);
      Obit_log_error(err, OBIT_InfoErr,"   quality %f, cleanable %f", 
		     in->quality[i], in->cleanable[i]);
      Obit_log_error(err, OBIT_InfoErr,"   maxAbsRes %f, avgRes %f", 
		     in->maxAbsRes[i], in->avgRes[i]);
      Obit_log_error(err, OBIT_InfoErr,"   imgRMS %f, imgPeakRMS %f",  
		     in->imgRMS[i], in->imgPeakRMS[i]);
      Obit_log_error(err, OBIT_InfoErr,"   beamPeakRMS %f, minFlux %f, gain %f",  
		     in->beamPeakRMS[i], in->minFlux[i], in->gain[i]);

    } /* End Loop over fields */
  } /* end print criteria */

  /* Maximum taper */
   bestTest = -1.0;
  for (i=0; i<n; i++) 
    bestTest = MAX (bestTest, in->mosaic->BeamTaper[i]);
  maxTape = MAX (0.0000001, bestTest);
 
  /* Find best quality */
  bestTest = -1.0;
  for (i=0; i<n; i++) 
    bestTest = MAX (bestTest, in->quality[i]);
  bestQuality = MAX (0.0000001, bestTest);

  /* Find best SNR */
  bestTest = -1.0;
  for (i=0; i<n; i++) 
    bestTest = MAX (bestTest, in->imgPeakRMS[i]);
  bestSNR = MAX (1.0, bestTest);

  /* Find best - decreasing bias towards lower taper */
  fact1 *= (0.1 + pow((1.0 - ((ofloat)in->Pixels->currentIter/in->Pixels->niter)), 2.0));
  bestTest = -1.0;
  best     = 0;
  done = TRUE;
  for (j=0; j<in->mosaic->numBeamTapes; j++) {
    i = bestTap[j];

    test = 
      /* Bias towards less tapering */
      fact1*((maxTape - in->mosaic->BeamTaper[i])/maxTape) +
      /* SNR */
      fact2*(in->imgPeakRMS[i]/bestSNR) + 
      /* Quality */
      fact3*(in->quality[i]/bestQuality);

    /* Is this one done? */
    if (in->cleanable[i]<in->minFlux[i]) test = 0.0;

    /* If SNR<5 degrade - careful, this measure is not always what it appears */
    if (in->imgPeakRMS[i]<1.5) test *= 0.25;
    else if (in->imgPeakRMS[i]<2.5) test *= 0.50;
    else if (in->imgPeakRMS[i]<5.0) test *= 0.75;
    if (test>bestTest) {
      bestTest = test;
      best     = j;
    }
    if (test>0.0) done = FALSE;  /* Finished CLEAN? */
    /* DEBUG */
    if ((test<=0.0) && (err->prtLv>=2)) {
      /* Explain */
      Obit_log_error(err, OBIT_InfoErr,"DEBUG best %d bestTest %f bestTap %d", 
		     best,bestTest,bestTap[j]);
      Obit_log_error(err, OBIT_InfoErr,"DEBUG cleanable %f minFlux %f quality %f", 
		     in->cleanable[i],in->minFlux[i],in->quality[i]);
      Obit_log_error(err, OBIT_InfoErr,"DEBUG in->imgPeakRMS[i] %f", 
		     in->imgPeakRMS[i]);
    }
    /* end DEBUG */
    if (err->prtLv>=2) 
      Obit_log_error(err, OBIT_InfoErr,"Resoln. %d objective fn %f", 
		     j, test);
  } /* End loop testing */

  /* Set taper */
  in->minBeamTaper = 0.95 * in->mosaic->BeamTaper[bestTap[best]];
  in->maxBeamTaper = 1.05 * in->mosaic->BeamTaper[bestTap[best]];
  if (err->prtLv>=2) 
    Obit_log_error(err, OBIT_InfoErr,"Select taper %d %f", 
		   best, in->mosaic->BeamTaper[bestTap[best]]);

  /* Set minFlux/gain for all facets with this taper by ratio of beam areas */
  if (best>0) {
    i = bestTap[0];    /* untapered */
    j = bestTap[best]; /* tapered  */
    if ((in->mosaic->images[j]->myDesc->beamMaj>0.0) && 
	(in->mosaic->images[j]->myDesc->beamMin>0.0) &&
	(in->mosaic->images[i]->myDesc->beamMaj>0.0) && 
	(in->mosaic->images[i]->myDesc->beamMin>0.0)) {
      beamrat = (in->mosaic->images[j]->myDesc->beamMaj * 
		 in->mosaic->images[j]->myDesc->beamMin) / 
	((in->mosaic->images[i]->myDesc->beamMaj * 
	  in->mosaic->images[i]->myDesc->beamMin));
      minFlux = in->minFlux[bestTap[0]] * beamrat;
      gain    = in->gain[bestTap[0]] / beamrat;
    } else {
      minFlux = in->minFlux[bestTap[0]];
      gain    = in->gain[bestTap[0]];
    }
    minT = in->minBeamTaper;
    maxT = in->maxBeamTaper;
   if (err->prtLv>=3) 
    Obit_log_error(err, OBIT_InfoErr,"MinFlux %f beam ratio %f gain %f", 
		   minFlux, beamrat, gain);
   for (i=0; i<n; i++) {
      /* Desired taper? And not native resolution */
      if (((in->mosaic->BeamTaper[i]>=minT) && 
	   (in->mosaic->BeamTaper[i]<=maxT) &&
	   (in->mosaic->BeamTaper[i]>0.0))) {
	in->Pixels->minFlux[i] = minFlux;  /* For actual cleaning */
	in->minFlux[i]         = minFlux;
	in->Pixels->gain[i]    = gain;     /* For actual cleaning */
	in->gain[i]            = gain; 
      }
    } /* end loop over facets */
  } /* end reset minFlux */
  return done;
} /* end SelectTaper */

/**
 * Checks if fresh facets are "Done"
 * if in->Pixels->complCode=OBIT_CompReasonMinFlux then fresh facets which 
 * have no more than 1 new component selected are markes as "done" in 
 * cleanable and quality.  
 * Call after Select.
 * \param in       CLEAN object
 * \param err      Obit error stack object.
 * \return the new object.
 */
void CheckIfDone(ObitDConCleanVis* in, ObitErr *err)
{
  olong i;
  gchar *routine = "ObitDConCleanVis:CheckIfDone";

  /* existing error check */
  if (err->error) return;

  /* Need to check? */
  if (in->Pixels->complCode!=OBIT_CompReasonMinFlux) return;

  /* Loop over fresh facets */
  for (i=0; i<in->nfield; i++) {
    if (in->fresh[i] && (in->mosaic->isAuto[i]<0) && 
	(in->mosaic->isShift[i]<0)) {
      /* More than 1 selected? Mark as "done"*/
      if (in->skyModel->startComp[i]>=in->Pixels->iterField[i]) {
	in->cleanable[i] = -in->cleanable[i]; 
	in->quality[i] = 0.0;
	/* Message if prtLv>=3 */
	if (err->prtLv>=3)
	  Obit_log_error(err, OBIT_InfoErr, "%s: Field %d marked done",
			 routine, i+1);
      } /* end done */
    } /* end if fresh and not shift/autoCen */
  } /* end loop */
} /* end CheckIfDone */

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
gboolean VisResetSkyModel(ObitDConCleanVis *in, ObitErr *err)
{
  gboolean doSub=FALSE;
  olong ncomp;
  ofloat sum;
  olong i, irow, it;
  gint32 dim[MAXINFOELEMDIM];
  olong *itemp=NULL, nfield;
  ObitTableCCRow *CCRow = NULL;
  const ObitDConCleanPxListClassInfo *pxListClass;
  gchar *routine = "VisResetSkyModel";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return doSub;
  g_assert (ObitDConCleanVisIsA(in));

  /* Reset PixelList parameters */
  pxListClass = (ObitDConCleanPxListClassInfo*)in->Pixels->ClassInfo; 
  pxListClass->ObitDConCleanPxListGetParms (in->Pixels, err);
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
    if ((ncomp>0) && (err->prtLv>1))
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
} /* end VisResetSkyModel */

/**
 * Reset Pixel List for beginning of a CLEAN
 * \param in   The Clean object
 * \param err Obit error stack object.
 */
void VisResetPixelList(ObitDConCleanVis *in, ObitErr *err)
{
  const ObitDConCleanPxListClassInfo *pxListClass;
  gchar *routine = "ResetPixelList";

  /* PxList class structure */
  pxListClass = (ObitDConCleanPxListClassInfo*)in->Pixels->ClassInfo; 

  pxListClass->ObitDConCleanPxListReset (in->Pixels, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end VisResetPixelList */

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
  
  ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, usePixels);
  
  ObitImageRead (image, image->image->array, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, usePixels);
  
  ObitImageClose (image, err);
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
 * \param pixarray If nonNULL, use these pixels rather than reading
 * \param err   Obit error stack object.
 * \return Maximum abs pixel value inside outer window but outside unboxes.
 */
ofloat ObitDConCleanVisCleanable(ObitDConCleanVis *in, olong field, 
				 ObitFArray *pixarray, ObitErr *err)
{
  ofloat out = -1.0;
  ObitImage *image=NULL;
  ObitFArray *usePixels=NULL;  
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat PeakIn, PeakOut, RMS;
  olong i,PeakInPos[2] = {0,0};
  gboolean doAbs, isMR;
  gchar *routine = "ObitDConCleanVisCleanable";

  /* error checks */
  if (err->error) return out;
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
    return out;
  }

  /* Is this a multiresolution image (to allow neg regions) */
  isMR = (in->mosaic->numBeamTapes>1);

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
    ObitImageOpen (image, OBIT_IO_ReadOnly, err);
    ObitImageRead (image, image->image->array, err);
    ObitImageClose (image, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, out);
    usePixels = ObitFArrayRef(image->image);  /* Pointer to image buffer */
}
  
  /* Allow negative for Stokes other than I or multiresolution */
  doAbs = isMR || (fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1);
    
  /* Get statistics */
  ObitDConCleanWindowStats (in->window, field, usePixels,
			    doAbs,
			    &PeakIn, &PeakInPos[0], &PeakOut, 
			    &RMS, err);
  /* Check for pathologies */
  if (fabs(PeakIn)>1.0e10)  PeakIn = 0.0;
  if (fabs(PeakOut)>1.0e10) PeakOut = 0.0;

  /* Free Image array */
  usePixels = ObitFArrayUnref(usePixels);
  if (pixarray==NULL) image->image = ObitFArrayUnref(image->image);
  
  if (err->error) Obit_traceback_val (err, routine, image->name, out);
  
  return MAX (fabs(PeakOut), fabs(PeakIn));
} /* end ObitDConCleanVisCleanable */

/**
 * Create Pixel list for cleaning
 * \param in       The Clean object
 * \param err      Obit error stack object.
 * \return TRUE if attempted, FALSE if cannot do 
 */
static void NewPxList (ObitDConCleanVis *in, ObitErr *err)
{
  gchar *routine = "ObitDConCleanVis:NewPxList";

  if (in->Pixels && (in->Pixels->nfield!=in->mosaic->numberImages)) 
    in->Pixels = ObitDConCleanPxListUnref(in->Pixels);
  if (!in->Pixels) {
    in->Pixels = ObitDConCleanPxListCreate("Pixel List", in->mosaic, 
					   in->maxPixel, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
  /* Reset  min. resid */
  in->Pixels->maxResid = 1.0e10;

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);

} /* end NewPxList */

/**
 * Create Pixel array for intermediate CLEAN
 * \param in       The Clean object
 * \param startCC  [out] Current number of components per field
 * \param err      Obit error stack object.
 * \return array of ObitFArrays for cleaning
 */
static ObitFArray** NewPxArray (ObitDConCleanVis *in, olong *startCC, 
				ObitErr *err)
{
  ObitFArray** pixarray=NULL;
  olong ifld;
  /*gchar *routine = "ObitDConCleanVis:NewPxArray";*/

  pixarray  = g_malloc0(in->nfield*sizeof(ObitFArray*));
  for (ifld=0; ifld<in->nfield; ifld++) {
    if (in->currentFields[ifld]<=0) break;  /* List terminated? */
    startCC[ifld]   = in->Pixels->iterField[in->currentFields[ifld]-1];  /* no. at start */
    /* Pixel array for field */
    pixarray[ifld]  = GetFieldPixArray (in, in->currentFields[ifld], err);  
  }
  
  return pixarray;
} /* end NewPxArray */

/**
 * Delete Pixel array for intermediate CLEAN
 * \param in       The Clean object
 * \param pixarray Array to delete
 * \return NULL
 */
static ObitFArray** KillPxArray (ObitDConCleanVis *in,  ObitFArray **pixarray)
{
  olong ifld;

  if (pixarray) {
    for (ifld=0; ifld<in->nfield; ifld++) {
      if (in->currentFields[ifld]<=0) break;  /* List terminated? */
      pixarray[ifld] = ObitFArrayUnref(pixarray[ifld]);
    }
    g_free(pixarray);  pixarray = NULL;
  }
  return pixarray;
} /* end KillPxArray */

/**
 * Delete Beam Patches on in
 * \param in       The Clean object
 */
static void KillBeamPatches (ObitDConCleanVis *in)
{
  olong i;
  
  if (in->BeamPatches) {
    for (i=0; i<in->mosaic->numberImages; i++) 
      if (in->BeamPatches[i]) ObitFArrayUnref(in->BeamPatches[i]);
    in->BeamPatches =  ObitMemFree (in->BeamPatches);
  }
} /* end KillBeamPatches */

/* GLOBAL DEBUG */
ObitDConCleanVis *damn;
/**
 * Low accuracy subtract pixels from image for current CLEAN fields.
 * Loops over fields in in->currentFields
 * Uses BeamPatch to subtract list of components
 * Updates cleanable member on in
 * \param in       The Clean object
 * \param newCC    Start CC -1 per in->mosaic field for subtraction
 * \param pixarray Array of arrays of pixels to subtract, in order of fields in 
 *                 in->currentFields
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
  ofloat PeakIn, PeakOut, RMS,  parms[20];
  olong i, j, ncc, ver, nThreads=0, nTh, nfield, ifield, nDo, nLeft, PeakInPos[2];
  gboolean doAbs, isMR, OK;
  gchar *tabType = "AIPS CC";
  gchar *routine = "SubNewCCs";
  gboolean DebugGDB=FALSE;  /*  DEBUG */

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisIsA(in));

  /* DEBUG */
  damn = in;

  /* Is this a multiresolution image (to allow neg regions) */
  isMR = (in->mosaic->numBeamTapes>1);

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
  /* NO, this only uses ObitFArrayShiftAdd which is threaded */
  nTh = 1;

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
  nTh = 1;   /* NO threading here */
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

  /* Image statistics */
  /* Loop over fields */
  for (i=0; i<nfield; i++) {
    ifield = in->currentFields[i]-1;
    if (ifield<0) break;
    /* Remeasure cleanable flux */
    /* Allow negative for Stokes other than I or multiresolution */
    doAbs = isMR || 
      (fabs (in->mosaic->images[ifield]->myDesc->crval[in->mosaic->images[ifield]->myDesc->jlocs]-1.0) > 0.1);
    
    /* Get statistics */
    ObitDConCleanWindowStats (in->window, ifield+1, pixarray[i],
			      doAbs,
			      &PeakIn, &PeakInPos[0], &PeakOut,
			      &RMS, err);
    in->cleanable[ifield]  = MAX(PeakOut, PeakIn);
    in->imgRMS[ifield]     = RMS;
    in->imgPeakRMS[ifield] = PeakIn/RMS;
  } /* end statistics loop over fields */
  
} /* end SubNewCCs */

/**
 * Subtract a set of lists of CLEAN components from an image in a thread
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
  olong i, j, ifield, len, pos1[2], pos2[2];
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
    /* g_assert (ObitFArrayIsA(comps)); DEBUG */
    
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
 * Delete arguments for ThreadImSub
 * No objects are Unreffed
 * \param nargs      number of elements in args.
 * \param ThreadArgs Array of ImSubFuncArg
 */
static void KillImSubFuncArgs (olong nargs, ImSubFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillImSubFuncArgs */

/**
 * Find maximum brightness using CC tables
 * Set in->peakFlux to the larger of in->peakFlux and the maximum 
 * brightness in the CC tables.
 * \param in      Object of interest.
 * \param err     ObitErr for reporting errors.
 */
static void FindPeak (ObitDConCleanVis *in, ObitErr *err)
{
  ObitImageMosaic *mosaic=in->mosaic;
  ofloat peakFlux;
  ObitTableCC *CCTab=NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  ofloat tmax, xcenter, ycenter, xoff, yoff, radius, beamrat, cells[2];
  olong   nfield, ifield, itemp, nccpos, CCVer;
  gchar *routine = "FindPeak";

  if (err->error) return;

   /* Number of fields */
  nfield = mosaic->numberImages;

  /* Get cellsize */
  cells[0] =  fabs(mosaic->xCells); cells[1] = fabs(mosaic->yCells);

  /* Consider components within 2.5  cells  */
  radius = 2.5 * cells[0];

  /* CC table(s) */
  itemp = 1;
  ObitInfoListGetTest(mosaic->info, "CCVer", &type, dim, &itemp);
  CCVer = itemp;

  /* If multiple channels */
  if (in->plane[0]>1) CCVer = in->plane[0]; 
  
  /* Loop over fields */
  peakFlux = -1.0e20;
  for (ifield=0; ifield<nfield; ifield++) { 
     /* Make CC table object including CCs from any overlapping facets */
    CCTab = ObitImageMosaicCombineCC (mosaic, ifield+1, CCVer, err);
    if  (err->error) Obit_traceback_msg (err, routine, mosaic->images[ifield]->name);

    /* Determine maximum */
    nccpos = CCTab->myDesc->nrow;
    ObitImageMosaicMaxCC (CCTab, nccpos, radius, &tmax, &xcenter, &ycenter, &xoff, &yoff, err);
    if  (err->error) Obit_traceback_msg (err, routine, mosaic->images[ifield]->name);

    /* Scale tmax by beam ratio */
    beamrat = (in->mosaic->images[0]->myDesc->beamMaj * 
	       in->mosaic->images[0]->myDesc->beamMin) / 
      ((in->mosaic->images[ifield]->myDesc->beamMaj * 
	in->mosaic->images[ifield]->myDesc->beamMin));
    tmax *= beamrat;
     
    /* get max */
    peakFlux = MAX (peakFlux, tmax);

    /* Delete temporary table */
    ObitDataZapTable((ObitData*)mosaic->images[ifield], CCTab->tabType, CCTab->tabVer, err);
    if  (err->error) Obit_traceback_msg (err, routine, mosaic->images[ifield]->name);
    CCTab = ObitTableCCUnref(CCTab);
 } /* end loop over fields */

  /* set peak */
  in->peakFlux = peakFlux;

  /* Tell */
  if (err->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Max. in mosaic chan %d %f", 
		   in->plane[0], in->peakFlux);
    ObitErrLog(err);  /* Progress Report */
  }

} /* end  FindPeak  */

/**
 * Determine if any facets are modified - have unsubtracted CCs
 * \param in     The Clean object
 * \return TRUE if some modified
 */
static gboolean anyDirty (ObitDConCleanVis *in)
{
  olong i;
  gboolean out=FALSE;
  /* Loop over facets */
  for (i=0; i<in->nfield; i++) {
    if (!in->fresh[i]) return TRUE;
  }

  return out; /* Must be OK */
} /* End anyDirty */
