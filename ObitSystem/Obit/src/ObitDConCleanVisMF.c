/* $Id: ObitDConCleanVisMF.c 149 2010-01-01 18:31:02Z bill.cotton $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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

#include "ObitImageMosaicMF.h"
#include "ObitImageUtil.h"
#include "ObitImageMF.h"
#include "ObitDConClean.h"
#include "ObitDConCleanVisMF.h"
#include "ObitMem.h"
#include "ObitFFT.h"
#include "ObitTableUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitSkyGeom.h"
#include "ObitDConCleanPxListMF.h"
#include "ObitSkyModelMF.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanVisMF.c
 * ObitDConCleanVisMF class function definitions.
 * Image based CLEAN class.
 * This class is derived from the ObitDCon class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanVisMF";

/** Function to obtain parent ClassInfo - ObitDConClean */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanVisGetClass;

/**
 * ClassInfo structure ObitDConCleanVisMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanVisMFClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private structures----------------*/
/* Image model subtraction threaded function argument */
typedef struct {
  /* CLEAN Object */
  ObitDConCleanVisMF *in;
  /* Input plane pixel data, */
  ObitFArray *inData;
  /* Input plane beam patch, */
  ObitFArray *inBeam;
  /* Field number (1-rel) of data in inData  */
  olong      ofield;
  /* Input array of compressed CLEAN component arrays */
  ObitFArray **inComp;
  /* Number of fields in inComp, fields in in->currentFields  */
  olong      nfield;
  /* Which (0-rel) channel number  */
  ofloat      iSpec;
  /* thread number, <0 -> no threading  */
  olong      ithread;
  /* Obit Thread object */
  ObitThread  *thread;
  /* Obit error stack object */
  ObitErr    *err;
} ImSubFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanVisMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanVisMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConCleanVisMFClassInfoDefFn (gpointer inClass);

/** Private: (re)make residuals. */
static void  MakeResiduals (ObitDConCleanVis *in, olong *fields, 
			    gboolean doBeam, ObitErr *err);

/** Private: (re)make all residuals. */
static void  MakeAllResiduals (ObitDConCleanVis *in, ObitErr *err);

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

/** Private: Convolve spectral CCs with a Gaussian. */
static ObitFArray* ConvlCC(ObitImage *image, olong CCVer, olong iterm, 
			   ObitErr *err);

/** Private: Cross convolve spectral CCs with a Gaussian. */
static void XConvlCC(ObitImage *in, olong CCVer, olong iterm, 
		     ObitImage *out, ObitFArray *outGrid, ObitErr *err);

/** Private: Apply Gaussian taper to uv grid. */
static void GaussTaper (ObitCArray* uvGrid,  ObitImageDesc *imDesc,
			ofloat gparm[3]);

/** Private: Read spectral plane beam patch. */
static ObitFArray* GetSpecBeamPatch (ObitDConCleanVisMF *in, ObitFArray *BP, 
				     olong ispec, ObitImageMF* image, 
				     ObitErr *err);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanVisMF* newObitDConCleanVisMF (gchar* name)
{
  ObitDConCleanVisMF* out;
  /*gchar *routine = "newObitDConCleanVisMF";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisMFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanVisMF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanVisMFInit((gpointer)out);

 return out;
} /* end newObitDConCleanVisMF */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanVisMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisMFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanVisMFGetClass */

/**
 * Make a deep copy of an ObitDConCleanVisMF.
 * \param inn  The object to copy
 * \param outt An existing object pointer for output or NULL if none exists.
 * \param err  Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanVis* ObitDConCleanVisMFCopy  (ObitDConCleanVis *inn, 
					     ObitDConCleanVis *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  ObitDConCleanVisMF *in  = (ObitDConCleanVisMF*)inn;
  ObitDConCleanVisMF *out = (ObitDConCleanVisMF*)outt;
  gchar *routine = "ObitDConCleanVisMFCopy";

  /* error checks */
  if (err->error) return (ObitDConCleanVis*)out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitDConCleanVisMF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, (ObitDConCleanVis*)out);

  /*  copy this class */
  return (ObitDConCleanVis*)out;
} /* end ObitDConCleanVisMFCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanVisMF similar to the input one.
 * \param inn  The object to copy
 * \param outt An existing object pointer for output, must be defined.
 * \param err  Obit error stack object.
 */
void ObitDConCleanVisMFClone  (ObitDConCleanVis *inn, ObitDConCleanVis *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitDConCleanVisMF *in  = (ObitDConCleanVisMF*)inn;
  ObitDConCleanVisMF *out = (ObitDConCleanVisMF*)outt;
  gchar *routine = "ObitDConCleanVisMFClone";

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
} /* end ObitDConCleanVisMFClone */

/**
 * Creates an ObitDConCleanVisMF 
 * defined for convenience of derived classes 
 * \param name   An optional name for the object.
 * \param uvdata from which to create object, should have all control
                 information defined on info member.
 * \param order  Order of the imaging, Spectral index only=1, plus curvature=2
 * \param maxFBW Max. IF center fractional bandwidth.
 * \param alpha  Spectral index correction previously applied to data.
 * \param err    Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVisMF* ObitDConCleanVisMFCreate (gchar* name, ObitUV *uvdata,  
					      olong order, ofloat maxFBW, 
					      ofloat alpha, ObitErr *err)
{
  olong nfield, i;
  ObitDConCleanVisMF* out=NULL;
  gchar *routine = "ObitDConCleanVisMFCreate";

 /* error checks */
  if (err->error) return out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVisMF (name);

  out->order = order;  /* Save order */

  /* Create UV imager and its ImageMosaic */
  out->imager = (ObitUVImager*)ObitUVImagerMFCreate("UVImagerMF", order, maxFBW, 
						    alpha, uvdata, err);
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
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->imgRMS[i]      = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
 }

  return out;
} /* end ObitDConCleanVisMFCreate */

/**
 * Creates an ObitDConCleanVisMF from optional components 
 * defined for convenience of derived classes 
 * \param name     An optional name for the object.
 * \param uvdata   from which to create object, should have all control
                   information defined on info member.
 * \param imager   Optional ObitUVImager to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
                   Should be an ObitUVImagerMF
 * \param skyModel Optional ObitSkyModel to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
 * \param order    Order of the imaging, Spectral index only=1, plus curvature=2
 * \param maxFBW   Max. IF center fractional bandwidth.
 * \param alpha  Spectral index correction previously applied to data.
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVis* 
ObitDConCleanVisMFCreate2 (gchar* name, ObitUV *uvdata,  
			   ObitUVImager *imager, ObitSkyModel *skyModel, 
			   olong order, ofloat maxFBW, ofloat alpha, 
			   ObitErr *err)
{
  olong nfield, i;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitDConCleanVisMF* out=NULL;
  ofloat ftemp;
  gchar *routine = "ObitDConCleanVisMFCreate";

 /* error checks */
  if (err->error) return (ObitDConCleanVis*)out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVisMF (name);

  out->order = order;  /* Save order */

  /* Use or create UV imager and create its ImageMosaic */
  if (imager==NULL) {
    out->imager =(ObitUVImager*) ObitUVImagerMFCreate("UVImager", order, maxFBW, 
						      alpha, uvdata, err);
    if (err->error) Obit_traceback_val (err, routine, name, (ObitDConCleanVis*)out);
  } else out->imager = ObitUVImagerRef(imager);

  /* Save uv Mosaic reference */
  out->mosaic = ObitUVImagerGetMosaic(out->imager, err);

  /*  Use or create SkyModel object */
  if (skyModel==NULL) out->skyModel = 
			(ObitSkyModel*)ObitSkyModelMFCreate ("SkyModel", out->mosaic);
  else out->skyModel = ObitSkyModelRef(skyModel);

  /* Copy control info to SkyModel */
  ObitInfoListCopyData(uvdata->info, out->skyModel->info);
  /* Disable any value of minFlux */
  ftemp = -1.0e20;
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (out->skyModel->info, "minFlux", OBIT_float, dim, &ftemp);
  
  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, (ObitDConCleanVis*)out);

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
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->imgRMS[i]      = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
 }

  return (ObitDConCleanVis*)out;
} /* end ObitDConCleanVisMFCreate2 */

/**
 * Read any base class parameters and then
 * read CLEAN control parameters from the ObitInfoList member:
 * \li "OrdFlux"    OBIT_float array  = min residual for orders 1 and 2
 * From Parent classes:
 * \li "Mode"      OBIT_long scalar = Model mode (ObitSkyModelMode) [def OBIT_SkyModel_Fastest]
 * \li "doRestore" OBIT_bool       = Restore image when done? [def TRUE]
 * \li "doXRestore" OBIT_bool      = Cross restore images when done? [def doRestore]
 * \li "doFlatten" OBIT_bool       = Flatten image when done? [def TRUE]
 * \li "doWeight"  OBIT_bool       = Weight UV data before imaging? [def TRUE]
 * \li "doBeam"    OBIT_bool       = Need to (re)make beam? [def TRUE]
 * \li "doRecenter" OBIT_bool      = Allow recentering autoCenter fields? [def TRUE]
 * \li "reuseFlux" OBIT_float      = Level of Components in initial 
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
void  ObitDConCleanVisMFGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConClassInfo *ParentClass;
  /*ObitInfoType type;*/
  /*gint32 dim[MAXINFOELEMDIM];*/
  /*olong i;*/
  /*union ObitInfoListEquiv InfoReal;*/
  ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;  /* as this class */
  gchar *routine = "ObitDConCleanVisMFGetParms";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* minimum residual flux for 1st and second order
  ObitInfoListGetTest(in->info, "OrdFlux", &type, dim, &in->OrdFlux); */

} /* end ObitDConCleanVisMFGetParms */

/** 
 * Restore components removed from the residual image(s)
 * Wideband imaging version.
 * Spectral orders higher than 0 are flux density weighted averages.
 * \param inn  The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanVisMFRestore(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  ObitFArray *convl=NULL;
  ObitImage *image=NULL;
  olong iplane, field, num, nOrd, plane[5]={1,1,1,1,1};
  gchar *routine = "ObitDConCleanVisMFRestore";

   /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisMFIsA(in));

  /* Anything to restore? */
  if (in->Pixels->currentIter<=0) return;

  /* Tell user */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Restoring components");
    ObitErrLog(err);  /* Progress Report */
  }

  /* Loop over fields */
  for (field = 0; field<in->nfield; field++) {

    /* Anything to restore? */
    if (in->Pixels->iterField[field]<=0) continue;

    /* which Image? */
    image = in->mosaic->images[field];

    /* Form combined image  */
    ObitImageMFCombine ((ObitImageMF*)image, FALSE, err);

    /* Restore Flux then individual channels */
    /* Convolve Gaussians */
    iplane = 0;
    convl = ConvlCC (image, in->CCver, iplane, err);
    /* Read image */
    plane[0] = 1;
    ObitImageGetPlane (image, NULL, plane, err);
    /* Sum */
    ObitFArrayAdd (image->image, convl, image->image);
    /* Rewrite */
    ObitImagePutPlane (image, NULL, plane, err);
    convl = ObitFArrayUnref(convl);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* individual channels */
    num  = ((ObitImageMF*)image)->nSpec;
    nOrd = ((ObitImageMF*)image)->maxOrder;
    for (iplane=1; iplane<=num; iplane++) {
      /* Convolve Gaussians */
      convl = ConvlCC (image, in->CCver, iplane, err);
      /* Read image */
      plane[0] = iplane+1 + nOrd;
      ObitImageGetPlane (image, NULL, plane, err);
      /* Sum */
      ObitFArrayAdd (image->image, convl, image->image);
      /* Rewrite */
      ObitImagePutPlane (image, NULL, plane, err);
      convl = ObitFArrayUnref(convl);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
    /* Free image memory */
    image->image = ObitFArrayUnref(image->image);
  } /* end loop over fields */

} /* end ObitDConCleanVisMFRestore */

/**
 * Restore components removed from one field but also 
 * appearing in another.  Does brute force convolution.
 * Wideband imaging version - does all spectral planes.
 * Spectral orders higher than 0 are flux density weighted averages.
 * Adopted from the AIPSish QOOP:QCLEAN.FOR(CLOVER)
 * Presumes in->mosaic and image descriptors filled in.
 * \param inn  The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanVisMFXRestore(ObitDConClean *inn, ObitErr *err)
{
 ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  ObitImage *image1=NULL, *image2=NULL;
  olong ifield, jfield, iplane, num, nOrd, plane[5]={1,1,1,1,1};
  gchar *routine = "ObitDConCleanVisMFXRestore";

   /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisMFIsA(in));

  /* Anything to restore? */
  if (in->Pixels->currentIter<=0) return;

  /* Tell user */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Cross Restoring components");
    ObitErrLog(err);  /* Progress Report */
  }

  /* Double loop over fields */
  for (jfield = 0; jfield<in->nfield; jfield++) {
    /* output image  */
    image1 = in->mosaic->images[jfield];
    num  = ((ObitImageMF*)image1)->nSpec;
    nOrd = ((ObitImageMF*)image1)->maxOrder;

    /* Restore Flux then individual channels */
    for (iplane=0; iplane<(num+1); iplane++) {
      /* Read image */
      if (iplane==0) plane[0] = 1;
      else plane[0] = 1+iplane+nOrd;
      ObitImageGetPlane (image1, NULL, plane, err);

      /* Loop over others */
      for (ifield = 0; ifield<in->nfield; ifield++) {
	/* Only cross */
	if (ifield==jfield) continue;

	/* Anything to restore? */
	if (in->Pixels->iterField[ifield]<=0) continue;

	/* which Image? */
	image2 = in->mosaic->images[ifield];
      
	/* Cross convolve Gaussians */
	XConvlCC (image2, in->CCver, iplane, image1, image1->image, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
      } /* end inner loop over fields */
      /* Rewrite */
      ObitImagePutPlane (image1, NULL, plane, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    } /* end loop over planes */
    
    image1->image = ObitFArrayUnref(image1->image);
  } /* end outer loop over fields */
} /* end ObitDConCleanVisMFXRestore */

/**
 * Flatten multiple facets if needed
 * Wideband imaging version.
 * Flattens all Spectral planes
 * Does Flatten if FullField member of mosaic member is defined.
 * \param inn  The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanVisMFFlatten(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  ObitImageMosaicClassInfo* mosaicClass; 
  gchar *routine = "ObitDConCleanVisMFFlatten";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisMFIsA(in));

  if ((in->mosaic->FullField!=NULL) && (in->mosaic->numberImages>1)) {
    /* Tell user */
    if (in->prtLv>1) {
      Obit_log_error(err, OBIT_InfoErr,"Flattening images");
      ObitErrLog(err);  /* Progress Report */
    }
    mosaicClass = (ObitImageMosaicClassInfo*)inn->mosaic->ClassInfo;
    mosaicClass->ObitImageMosaicFlatten (in->mosaic, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitDConCleanVisMFFlatten */

/**
 * Subtract components from uv data.
 * Frees any work arrays on mosaic images
 * \param inn  The object to deconvolve
 * \param err  Obit error stack object.
 */
void ObitDConCleanVisMFSub(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  const ObitDConCleanVisClassInfo *parentClass;
  olong ifld;
  gchar *routine = "ObitDConCleanVisMFSub";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisMFIsA(in));

  /* Mark All as not fresh */
  for (ifld=0; ifld<in->nfield; ifld++) {
    ((ObitImageMF*)in->mosaic->images[ifld])->fresh = FALSE;
  }
  
  /* Most work in parent class */
  parentClass = myClassInfo.ParentClass;
  parentClass->ObitDConCleanSub((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitDConCleanVisMFSub */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanVisMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanVisMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanVisMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanVisMFClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanVisMFClassInfo *theClass = (ObitDConCleanVisMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanVisMFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanVisMFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanVisMFGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanVisMF;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanVisMFCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanVisMFClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanVisMFInit;
  theClass->ObitDConGetParms     = (ObitDConGetParmsFP)ObitDConCleanVisMFGetParms;
  theClass->ObitDConCleanSub     = (ObitDConCleanSubFP)ObitDConCleanVisMFSub;
  theClass->ObitDConCleanRestore = (ObitDConCleanRestoreFP)ObitDConCleanVisMFRestore;
  theClass->ObitDConCleanFlatten = (ObitDConCleanFlattenFP)ObitDConCleanVisMFFlatten;
  theClass->ObitDConCleanXRestore= (ObitDConCleanXRestoreFP)ObitDConCleanVisMFXRestore;

  /* Private functions definitions for derived classes */
  theClass->MakeResiduals   = (MakeResidualsFP)MakeResiduals;
  theClass->MakeAllResiduals= (MakeAllResidualsFP)MakeAllResiduals;
  theClass->SubNewCCs       = (SubNewCCsFP)SubNewCCs;
  theClass->NewPxList       = (NewPxListFP)NewPxList;
} /* end ObitDConCleanVisMFClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanVisMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVisMF *in = inn;

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
} /* end ObitDConCleanVisMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanVisMF* cast to an Obit*.
 */
void ObitDConCleanVisMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVisMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanVisMFClear */

/**
 * Make selected residual images and get statistics, 
 * \param inn    The Clean object
 * \param fields zero terminated list of field numbers to image
 * \param doBeam If TRUE also make beam
 * \param err    Obit error stack object.
 */
static void  MakeResiduals (ObitDConCleanVis *inn, olong *fields, 
			    gboolean doBeam, ObitErr *err)
{
  ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  const ObitDConCleanVisClassInfo *parentClass;
  olong ifld, field;
  gchar *routine = "ObitDConCleanVisMF:MakeResiduals";
  
  if (err->error) return; /* prior error condition? */

  parentClass = myClassInfo.ParentClass;
  /* Call MakeResiduals in parent class */
  parentClass->MakeResiduals (inn, fields, doBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Mark images as fresh - make combined images */
  for (ifld=0; ifld<in->nfield; ifld++) {
    field = fields[ifld];
    if (field<=0) break;
    ((ObitImageMF*)in->mosaic->images[field-1])->fresh = TRUE;
  }
  
} /* end MakeResiduals */

/**
 * Make all residual images and get statistics
 * \param inn    The Clean object
 * \param err    Obit error stack object.
 */
static void  MakeAllResiduals (ObitDConCleanVis *inn, ObitErr *err)
{
 ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  const ObitDConCleanVisClassInfo *parentClass;
  olong ifld;
  gchar *routine = "ObitDConCleanVisMF:MakeAll Residuals";
 
  if (err->error) return; /* prior error condition? */

  parentClass = myClassInfo.ParentClass;
  /* Call MakeAllResiduals in parent class */
  parentClass->MakeAllResiduals (inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Mark images as fresh */
  for (ifld=0; ifld<in->nfield; ifld++) {
    ((ObitImageMF*)in->mosaic->images[ifld])->fresh = TRUE;
  }
  
} /* end MakeAllResiduals */

/**
 * Create Pixel list for cleaning
 * \param inn      The Clean object
 * \param err      Obit error stack object.
 */
static void NewPxList (ObitDConCleanVis *inn, ObitErr *err)
{
  ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  gchar *routine = "ObitDConCleanVis:NewPxList";

  if (in->Pixels && (in->Pixels->nfield!=in->mosaic->numberImages)) 
    in->Pixels = ObitDConCleanPxListUnref(in->Pixels);
  if (!in->Pixels) {
    in->Pixels = 
      (ObitDConCleanPxList*)ObitDConCleanPxListMFCreate("Pixel List", in->mosaic, 
							in->imager->uvwork, in->maxPixel, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
  /* Reset  min. resid */
  in->Pixels->maxResid = 1.0e10;

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);

} /* end NewPxList */

/**
 * Low accuracy subtract pixels from image for current CLEAN fields.
 * Loops over fields in in->currentFields
 * Uses pixel list to subtract list of components
 * Updates cleanable member on in
 * Threading is over spectral channels
 * \param inn      The Clean object
 * \param newCC    Start CC -1 per in->mosaic field for subtraction
 * \param pixarray Array of arrays of pixels to subtract, in->maxOrder+1 per field
 * \param err      Obit error stack object.
 * \return TRUE if attempted, FALSE if cannot do 
 */
static void SubNewCCs (ObitDConCleanVis *inn, olong *newCC, ObitFArray **pixarray, 
		       ObitErr *err)
{
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL;
  ObitImageMF *image=NULL;
  ObitImage *theBeam=NULL;
  ImSubFuncArg **threadArgs;
  ObitFArray **comps=NULL;
  ObitFArray **inFArrays, **bmFArrays;
  ObitImageDesc *outDesc;
  ofloat PeakIn, PeakOut, RMS,  parms[20];
  olong i, j, l, ip, jp, ncc, ver, nThreads=0, nTh, nfield, ifield, nDo, nLeft, PeakInPos[2];
  olong ispec, naxis[2], plane[5] = {1,1,1,1,1};
  gboolean doAbs, OK;
  gchar *tabType = "AIPS CC";
  ObitDConCleanVisMF *in = (ObitDConCleanVisMF*)inn;
  ObitImageClassInfo *imgClass;
  gchar *routine = "SubNewCCs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisMFIsA(in));

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
    comps[i] = ObitTableCCUtilMergeSelSpec (CCTable, newCC[i]+1, ncc, parms, err);
    CCTable  = ObitTableUnref(CCTable);
    Obit_return_if_fail ((comps[i]!=NULL),
			 err, "%s: Error merging CLEAN components for field %d",
			 routine, ifield);

  } /* end loop over fields collecting CCs */

  /* setup for threaded processing
     Initialize Threading */
  nThreads = MakeImSubFuncArgs (in->thread, err, &threadArgs);
  /* No more threads than work to spread them across */
  if (((ObitImageMF*)in->mosaic->images[0])->nSpec>1) nTh = nThreads;
  else nTh = 1;

  /* Initialize Thread args */
  for (i=0; i<nTh; i++) {
    threadArgs[i]->in     = in;
    threadArgs[i]->inComp = comps;
    threadArgs[i]->nfield = nfield;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
  }

  /* Create inFArrays/bmFArrays work arrays for image planes */
  inFArrays = g_malloc0(nTh*sizeof(ObitFArray*));
  bmFArrays = g_malloc0(nTh*sizeof(ObitFArray*));

 /* Loop over fields */
  for (i=0; i<nfield; i++) {
    ifield = in->currentFields[i];
    if (ifield<0) break;

    /* Which image? */
    image = (ObitImageMF*)in->mosaic->images[ifield-1];
    /* Which Beam? */
    imgClass  = (ObitImageClassInfo*)image->ClassInfo;    /* Image class */
    theBeam   = imgClass->ObitImageGetBeam((ObitImage*)image, 0, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    nTh = MIN (nThreads, image->nSpec);
    nLeft = image->nSpec;
    ip = 0;
    
    /* Create inFArrays/bmFArrays work arrays for image planes */
    naxis[0] = image->myDesc->inaxes[0];
    naxis[1] = image->myDesc->inaxes[1];
    for (l=0; l<nTh; l++) 
      inFArrays[l] = ObitFArrayCreate (NULL, 2, naxis);
    
    /* Loop over spectral planes */
    for (ispec=0; ispec<image->nSpec; ispec+=nTh) {
      nDo = MIN (nTh, nLeft);
      /* Args for this batch */
      jp = 0;
      for (j=0; j<nDo; j++) {
	threadArgs[j]->ofield = ifield;
	if (nDo>1) threadArgs[j]->ithread = j;
	else threadArgs[j]->ithread = -1;
	/* Read plane */
	plane[0] = 2+image->maxOrder+ispec+j;
	ObitImageGetPlane ((ObitImage*)image, inFArrays[jp]->array, plane, err);
	if (err->error) Obit_traceback_msg (err, routine, image->name);
	threadArgs[j]->inData = inFArrays[jp];
	threadArgs[j]->iSpec = ispec+j;
	/* Read beam patch */
	bmFArrays[jp] = GetSpecBeamPatch (in, in->BeamPatches[ifield-1], 
					  ispec+j, image, err);
	if (err->error) Obit_traceback_msg (err, routine, image->name);
	threadArgs[j]->inBeam = bmFArrays[jp];
	jp++;
      }
      /* Dummy any remainders */
      for (j=nDo; j<nTh; j++) {
	bmFArrays[jp] = ObitFArrayUnref(bmFArrays[jp]);
	threadArgs[j]->ofield = 1;
	threadArgs[j]->inData = NULL;
	threadArgs[j]->inBeam = NULL;
	jp++;
      }
      
      /* Do operation - need to stub remainder of nTh (don't know why) */
      OK = ObitThreadIterator (in->thread, nDo, 
			       (ObitThreadFunc)ThreadImSub,
			       (gpointer**)threadArgs);
      
      /* Check for problems */
      if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
      if (err->error) goto cleanup;

      /* Rewrite planes, delete spectral beam patches  */
      jp = 0;
      for (j=0; j<nDo; j++) {
	plane[0] = 2+image->maxOrder+ispec+j;
	ObitImagePutPlane ((ObitImage*)image, inFArrays[jp]->array, plane, err);
	if (err->error) Obit_traceback_msg (err, routine, image->name);
	bmFArrays[jp] = ObitFArrayUnref(bmFArrays[jp]);
 	jp++;
      }
      ip    += nDo;
      nLeft -= nDo;
      if (nLeft<=0) break;
    } /* end loop over spectral planes */

    /* Form combined image and copy to pixarray */
    ObitImageMFCombine (image, TRUE, err);
    plane[0] = 1;
    ObitImageGetPlane ((ObitImage*)image, pixarray[i]->array, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* Be sure to delete inFArrays */
    for (j=0; j<nTh; j++) {
      inFArrays[j]  = ObitFArrayUnref(inFArrays[j]);
      bmFArrays[i] = ObitFArrayUnref(bmFArrays[i]);
    }
  } /* end loop over field subtracting */

  
  /* Cleanup */
 cleanup:
  g_free(inFArrays); inFArrays = NULL;
  g_free(bmFArrays); bmFArrays = NULL;
  KillImSubFuncArgs (nThreads, threadArgs);
  if (comps) {
    for (i=0; i<nfield; i++) comps[i] = ObitFArrayUnref(comps[i]);
    g_free(comps);
  }
  
  /* Image statistics */
  /* Loop over fields */
  ip = 0;
  for (i=0; i<nfield; i++) {
    ifield = in->currentFields[i]-1;
    if (ifield<0) break;
    /* Remeasure cleanable flux */
    /* Allow negative for Stokes other than I */
    outDesc = in->mosaic->images[ifield]->myDesc;
    doAbs   = fabs (outDesc->crval[outDesc->jlocs]-1.0) > 0.1;
    
    /* Get statistics */
    ObitDConCleanWindowStats (in->window, ifield+1, pixarray[ip],
			      doAbs,
			      &PeakIn, &PeakInPos[0], &PeakOut,
			      &RMS, err);
    in->cleanable[ifield]  = MAX (PeakOut, PeakIn);
    in->imgRMS[ifield]     = RMS;
    in->imgPeakRMS[ifield] = PeakIn/RMS;
    ip++;
  } /* end statistics loop over fields */
  
} /* end SubNewCCs */

/**
 * Subtract a set of lists of CLEAN components from an image in a thread
 * Uses BeamPatches on in.
 * \param arg Pointer to ImSubFuncArg argument with elements:
 * \li in       The ObitDConCleanVis object
 * \li inData   Pixel array to be subtracted from
 * \li inBeam   Field/plane beam
 * \li field    Field number (1-rel) of data in inData 
 * \li inComp   Array of compressed arrays of point pixels
 *              Correspond to fields in->currentFields
 * \li nfield   Number of elements in inComp
 * \li iSpec    Which channel number
 * \li ithread  thread number, <0 -> no threading
 * \li err      ObitErr Obit error stack object
 * \li thread   thread Object
 * \return NULL
 */
static gpointer ThreadImSub (gpointer args)
{
  /* Get arguments from structure */
  ImSubFuncArg *largs = (ImSubFuncArg*)args;
  ObitDConCleanVisMF *in  = largs->in;
  ObitFArray *pixarray  = largs->inData;
  ObitFArray *BeamPatch = largs->inBeam;
  olong       ofield    = largs->ofield;
  ObitFArray **inComp   = largs->inComp;
  olong       iSpec     = largs->iSpec;
  olong       nfield    = largs->nfield;
  ObitErr    *err       = largs->err;
  ObitThread *thread    = largs->thread;

  /* local */
  olong i, j, ifield, len, pos1[2], pos2[2];
  ofloat idx, idy, ftemp, flux, offset[2];
  ObitFArray *comps;
  ObitImageDesc *inDesc, *outDesc;

  if (err->error)     goto finish;
  if (pixarray==NULL) goto finish;
  if (ofield<=0)      goto finish;


  /* Loop over fields with components */
  for (i=0; i<nfield; i++) {
    ifield = in->currentFields[i];
    if (ifield<=0) break;
    comps = inComp[i];          /* Component list */
    if (comps==NULL) continue;  /* Any components? */
    g_assert (ObitFArrayIsA(comps));
    
    /* Loop over CCs - assumes all points */
    inDesc    = in->mosaic->images[ifield-1]->myDesc;
    outDesc   = in->mosaic->images[ofield-1]->myDesc;
    len       = comps->naxis[0];
    pos2[0]   = BeamPatch->naxis[0]/2; 
    pos2[1]   = BeamPatch->naxis[1]/2;    /* Center of beam patch */
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
      flux = comps->array[j*len+3+iSpec];
      ObitFArrayShiftAdd(pixarray, pos1, BeamPatch, pos2, -flux, pixarray);
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
    (*ThreadArgs)[i]->inBeam     = NULL;
    (*ThreadArgs)[i]->ofield     = 0;
    (*ThreadArgs)[i]->inComp     = NULL;
    (*ThreadArgs)[i]->nfield     = 0;
    (*ThreadArgs)[i]->iSpec      = 0;
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
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillImSubFuncArgs */

/** 
 * Convolve a set of Clean components with a beam returning an image array
 * 
 * \param image  Image with CC table and defines size of image grid and
 *               with beam to convolve with.
 * \param CCVer  CC table number
 * \param iterm  Select spectral term, 0=flux, higher is a spectral channel..
 * \param err    Obit error stack object.
 * \return An array with the Clean components convolved
 */
static ObitFArray* ConvlCC(ObitImage *image, olong CCVer, olong iterm, ObitErr *err)
{
  ObitIOCode retCode;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitImageDesc *imDesc = NULL;
  olong first, last, ncomp, ver, ndim, naxis[2], ddim[2];
  gchar *tabType = "AIPS CC";
  ofloat gparm[3], bmaj, bmin, bpa;
  ObitFArray *grid = NULL;
  ObitCArray *uvGrid = NULL;
  ObitFFT *forFFT = NULL, *revFFT = NULL;
  gchar *routine = "ConvlCC";

  /* error checks */
  if (err->error) return grid;

  /* Open Image */
  retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, grid);
  
  /* Restoring beam */
  imDesc = image->myDesc;
  bmaj   = imDesc->beamMaj;
  bmin   = imDesc->beamMin;
  bpa    = imDesc->beamPA;
    
  /* Close */
  retCode = ObitImageClose (image, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, grid);

  /* Get CC table */
  ver = CCVer;
  tempTable = newObitImageTable (image, OBIT_IO_ReadWrite, tabType, &ver, err);
  if ((tempTable==NULL) || (err->error)) Obit_traceback_val (err, routine, image->name, grid);
  CCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_val (err, routine, image->name, grid);

  /* Open and close to get header */
  ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
  ObitTableCCClose (CCTable, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, grid);

  /* Grid components */
  first = 0;  /* all */
  last = 0;
  /* Need routine to use product of flux with another parameter */
  /* Return gridded flux*term */
  /* Spectral or normal */
  if ((CCTable->noParms>4) && (iterm>0)) { /* Spectral */
    retCode = ObitTableCCUtilGridSpect (CCTable, 1, iterm, &first, &last, FALSE, 
					1.0, 0.0, 1.0e20, imDesc, &grid, gparm, &ncomp, 
					err);
  } else { /* normal */
    retCode = ObitTableCCUtilGrid (CCTable, 1, &first, &last, FALSE, 
				   1.0, 0.0, 1.0e20, imDesc, &grid, gparm, &ncomp, 
				   err);
  }
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, grid);
  
  /* Free CC table */
  CCTable = ObitTableCCUnref(CCTable);
  
  /* DEBUG */
  /* ObitImageUtilArray2Image ("DbugGridCC.fits", 1, grid, err);  */
  /* if (err->error) Obit_traceback_val (err, routine, image->name, grid);*/
  /* END DEBUG */
  
  /* FFT to image plane */
  ObitFArray2DCenter (grid); /* Swaparoonie to FFT order */
  /* Make Output of FFT if needed */
  ndim = 2;
  naxis[0] = 1+grid->naxis[1]/2; naxis[1] = grid->naxis[0]; 
  if (uvGrid) uvGrid = ObitCArrayRealloc(uvGrid, ndim, naxis);
  else uvGrid = ObitCArrayCreate ("FFT output", ndim, naxis);
  /* Create Forward FFT or reuse if OK */
  ddim[0] = grid->naxis[1]; ddim[1] = grid->naxis[0];
  if ((!forFFT) || (ddim[0]!=forFFT->dim[0]) || (ddim[1]!=forFFT->dim[1])) {
    forFFT = ObitFFTUnref(forFFT);
    forFFT = newObitFFT("FFT:FTImage", OBIT_FFT_Forward, 
			OBIT_FFT_HalfComplex, 2, ddim);
  }
  /* FFT */
  ObitFFTR2C (forFFT, grid, uvGrid);
  /* Put the center at the center */
  ObitCArray2DCenter (uvGrid);
  
  /* DEBUG */
  /*tempFArray = ObitCArrayMakeF(uvGrid);*/  /* Temp FArray */
  /*ObitCArrayReal (uvGrid, tempFArray);*/   /* Get real part */
  /*tempFArray = ObitFArrayUnref(tempFArray);*/   /* delete temporary */
  /* END DEBUG */
  
  /* Taper for Gaussian - Use restoring beam or Gaussians from table if any 
     add rotation of image */
  if (gparm[0]<0.0) { /* restoring beam */
    gparm[0] = bmaj;
    gparm[1] = bmin;
    gparm[2] = bpa + image->myDesc->crota[image->myDesc->jlocd];
  } else { /* Gaussians from table */
    gparm[0] = gparm[0];
    gparm[1] = gparm[1];
    gparm[2] = gparm[2] + image->myDesc->crota[image->myDesc->jlocd];
  }
  GaussTaper (uvGrid, imDesc, gparm);
  
  /* DEBUG */
  /*tempFArray = ObitCArrayMakeF(uvGrid); */  /* Temp FArray */
  /*ObitCArrayReal (uvGrid, tempFArray); */   /* Get real part */
  /*ObitImageUtilArray2Image ("DbuguvGridAfter.fits", 1, tempFArray, err); */
  /*tempFArray = ObitFArrayUnref(tempFArray); */   /* delete temporary */
  /* END DEBUG */
  
  /* FFT back to image */
  ObitCArray2DCenter (uvGrid); /* Swaparoonie to FFT order */
  /* Create reverse FFT or reuse if OK */
  ddim[0] = grid->naxis[1]; ddim[1] = grid->naxis[0];
  if ((!revFFT) || (ddim[0]!=revFFT->dim[0]) || (ddim[1]!=revFFT->dim[1])) {
    revFFT = ObitFFTUnref(revFFT);
    revFFT = newObitFFT("FFT:FTuv", OBIT_FFT_Reverse, 
			OBIT_FFT_HalfComplex, 2, ddim);
  }
  /* FFT */
  ObitFFTC2R (revFFT, uvGrid, grid);
  /* Put the center at the center */
  ObitFArray2DCenter (grid);
  
  /* DEBUG  
  if (iterm==1)
    ObitImageUtilArray2Image ("DbugRestCC.fits", 0, grid, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, grid); */
  /* END DEBUG */
  
  /* Cleanup */
  forFFT = ObitFFTUnref(forFFT);
  revFFT = ObitFFTUnref(revFFT);
  uvGrid = ObitCArrayUnref(uvGrid);

  return grid;
} /* end ConvlCC */

/** 
 * Cross convolve a set of Clean components with a beam returning an image array
 * \param in      Image with CC table and defines size of image grid and
 *                with beam to convolve with.
 * \param CCVer   CC table number
 * \param iterm   Select spectral term, 0=flux, higher are spectral channels
 * \param out     Image defining output grid
 * \param outGrid Grid onto which to accumulate
 * \param err     Obit error stack object.
 */
static void XConvlCC(ObitImage *in, olong CCVer, olong iterm, 
		     ObitImage *out, ObitFArray *outGrid, ObitErr *err)
{
  ObitImageDesc *imDesc1=NULL, *imDesc2=NULL;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitFArray *list = NULL, *tmpArray = NULL;
  olong ver, ncomp, ndim, naxis[2];
  ofloat gparm[3], gauss[3], bmaj, bmin, bpa, sr, cr, cellx, celly;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitDConCleanVisMF:XConvlCC";

  /* error checks */
  if (err->error) return;
  
  imDesc1 = in->myDesc;
  imDesc2 = out->myDesc;
  
  /* Any overlap? */
  if (!ObitImageDescOverlap(imDesc1, imDesc2, err)) return;
  
  /* Get CC table */
  ver = CCVer;
  tempTable = newObitImageTable (in, OBIT_IO_ReadWrite, tabType, &ver, err);
  CCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Get list from jfield */
  if (iterm>0)
    list = ObitTableCCUtilCrossListSpec (CCTable, imDesc1, imDesc2, 
					 gparm, &ncomp, iterm, err);
  else
    list = ObitTableCCUtilCrossList (CCTable, imDesc1, imDesc2, 
				     gparm, &ncomp, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Free CC table */
  CCTable = ObitTableCCUnref(CCTable);
  
  /* Anything to do? */
  if (ncomp<=0) {ObitFArrayUnref(list); return;}

  /* Setup for ifield if needed */
  if (!tmpArray) {  /* Create FArray (zero filled) */
    ndim = 2; naxis[0] = imDesc2->inaxes[0]; naxis[1] = imDesc2->inaxes[1];
    tmpArray = ObitFArrayCreate ("Image for CCs", ndim, naxis);
  }
  
  /* Set Gaussian parameters  Use beam from CC table or header? */
  if (gparm[0]<0.0) {
    bmaj = imDesc1->beamMaj;
    bmin = imDesc1->beamMin;
    bpa  = imDesc1->beamPA;
  } else {
    bmaj = gparm[0];
    bmin = gparm[1];
    bpa  = gparm[2];
  }

  cellx = imDesc1->cdelt[0];
  celly = imDesc1->cdelt[1];
  cr = cos ((bpa + imDesc1->crota[imDesc1->jlocd])*DG2RAD);
  sr = sin ((bpa + imDesc1->crota[imDesc1->jlocd])*DG2RAD);
  gauss[0] = ((cr*cr)/(bmin*bmin) + (sr*sr)/(bmaj*bmaj)) *
    cellx*cellx*4.0*log(2.0);
  gauss[1] =  ((sr*sr)/(bmin*bmin) + (cr*cr)/(bmaj*bmaj)) *
    celly*celly*4.0*log(2.0);
  gauss[2] = (1.0/(bmin*bmin) - 1.0/(bmaj*bmaj)) *
    sr*cr*fabs(celly*celly)*8.0*log(2.0);

  /* Convolve list to tmpArray */
  ObitFArrayConvGaus (tmpArray, list, ncomp, gauss);
  
  /* Accumulate */
  ObitFArrayAdd (outGrid, tmpArray, outGrid);
    
  /* Cleanup */
  list = ObitFArrayUnref(list);
  tmpArray = ObitFArrayUnref(tmpArray);
  
  /* DEBUG
  ObitImageUtilArray2Image ("DbugCrossConv.fits", 0, outGrid, err); */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end XConvlCC */

/**
 * Apply Gaussian taper to a half Plane Complex grid 
 * assumed in the form from an ObitFFT.
 * NOTE: the uv grid is different in Obit (FFTW) and AIPS.
 * \param uvGrid Grid to be tapered
 * \param imDesc Image descriptor for image of which uvGrid is the FFT.
 * \param Gaussian in units of degrees, bmaj, bmin, bpa
 */
static void GaussTaper (ObitCArray* uvGrid, ObitImageDesc *imDesc,
			ofloat gparm[3])
{
  ofloat dU, dV, UU, VV, texp;
  ofloat konst, xmaj, xmin, cpa, spa, b1, b2, b3, bb2, bb3;
  ofloat taper, norm, *grid, tx, ty;
  olong i, j, nx, ny, ndim, naxis[2];

  /* Image info - descriptor should still be valid */
  nx = imDesc->inaxes[imDesc->jlocr];
  ny = imDesc->inaxes[imDesc->jlocd];
  
  /* Normalization factor */
  norm = ((ofloat)nx) * ((ofloat)ny);
  tx = MAX (1.0/sqrt(1.1331), gparm[0]/fabs(imDesc->cdelt[imDesc->jlocr]));
  ty = MAX (1.0/sqrt(1.1331), gparm[1]/fabs(imDesc->cdelt[imDesc->jlocd]));

  norm = 1.1331 * tx * ty / norm;

  /* UV cell spacing */
  dU = RAD2DG /  (nx * fabs(imDesc->cdelt[imDesc->jlocr]));
  dV = RAD2DG /  (ny * fabs(imDesc->cdelt[imDesc->jlocd]));
  
  konst = DG2RAD * G_PI * sqrt (0.5) / 1.17741022;
  xmaj = gparm[0] * konst;
  xmin = gparm[1] * konst;
  cpa = cos (DG2RAD * (90.0 + gparm[2])); /* is this right? */
  spa = sin (DG2RAD * (90.0 + gparm[2]));
  b1 = -(((cpa*xmaj)*(cpa*xmaj)) + ((spa*xmin)*(spa*xmin)));
  b2 = -(((spa*xmaj)*(spa*xmaj)) + ((cpa*xmin)*(cpa*xmin)));
  b3 = - 2.0 * spa * cpa * (xmaj*xmaj - xmin*xmin);
  
  /* pointer to complex grid */
  ndim = 2; naxis[0] = 0; naxis[1] = 0; 
  grid = ObitCArrayIndex(uvGrid, naxis);
  
  /* loop over uv array */  
  for (i=0; i<ny; i++) {
    VV = dV * (i-nx/2);
    UU = 0.0;
    bb2 = b2 * VV * VV;
    bb3 = b3 * VV;
    /* Loop down row computing, applying taper */
    for (j=0; j<1+nx/2; j++) {
      texp = b1 * UU * UU + bb2 + bb3 * UU;
      if (texp>-14.0) taper = norm * exp (texp);
      else  taper = 0.0;
      UU = UU + dU;
      grid[2*j]   *= taper;
      grid[2*j+1] *= taper;
    }
    grid += 2*uvGrid->naxis[0];
  }

} /* end GaussTaper */

/**
 * Read spectral beam patch for a spectral plane in the 
 * beam for image
 * \param in    The Clean object 
 * \param BP    Combined beam patch (used to determine size)
 * \param ispec 0-rel spectral channel number
 * \param image Image with Beam to extract beam patches
 * \param err   Obit error stack object.
 */
static ObitFArray* GetSpecBeamPatch (ObitDConCleanVisMF *in, ObitFArray *BP, 
				     olong ispec, ObitImageMF* image, 
				     ObitErr *err)
{
  ObitImage *theBeam;
  olong ablc[2], atrc[2], pos[2], plane[]={1,1,1,1,1};
  olong nx, ny, icenx, iceny, beamPatchSize;
  ofloat fmax;
  ObitImageClassInfo *imgClass;
  ObitFArray *out=NULL, *FAtemp=NULL;
  gchar *routine = "GetSpecBeamPatch";

  /* error checks */
  if (err->error) return out;

  imgClass  = (ObitImageClassInfo*)image->ClassInfo;    /* Image class */
  theBeam   = imgClass->ObitImageGetBeam((ObitImage*)image, 0, plane, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Get beam patch size same as for combined beam */
  beamPatchSize = BP->naxis[0]/2;

  /* Center pixel  */
  nx = theBeam->myDesc->inaxes[0];
  ny = theBeam->myDesc->inaxes[1];
  
  /* Read whole plane */
  plane[0] = ispec+image->maxOrder+2;
  ObitImageGetPlane (theBeam, NULL, plane, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, out);
    
  /* Setup to find peak in inner quarter - corners may be high */
  icenx = nx/2;
  iceny = ny/2;
  ablc[0] = icenx - nx/4;
  atrc[0] = icenx + nx/4;
  ablc[1] = iceny - nx/4;
  atrc[1] = iceny + nx/4;
  
  /* Find center in inner quarter */
  FAtemp = ObitFArraySubArr(theBeam->image, ablc, atrc, err);
  if (err->error) Obit_traceback_val (err, routine, theBeam->name, out);
  fmax = ObitFArrayMax (FAtemp, pos);
  FAtemp = ObitFArrayUnref(FAtemp);
  icenx = pos[0]+ablc[0];
  iceny = pos[1]+ablc[1];
  
  /* Beam patch window as 0-rel */
  ablc[0] = icenx - beamPatchSize;
  atrc[0] = icenx + beamPatchSize;
  ablc[1] = iceny - beamPatchSize;
  atrc[1] = iceny + beamPatchSize;
  
  /* Output Beam patch */
  out = ObitFArraySubArr(theBeam->image, ablc, atrc, err);
  if (err->error) Obit_traceback_val (err, routine, theBeam->name, out);
    
  /* Free Image array? */
  theBeam->image = ObitFArrayUnref(theBeam->image);
  return out;
} /* end GetSpecBeamPatch  */

