/* $Id: ObitDConCleanVisWB.c 149 2010-01-01 18:31:02Z bill.cotton $  */
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

#include "ObitImageMosaicWB.h"
#include "ObitImageUtil.h"
#include "ObitImageWB.h"
#include "ObitDConClean.h"
#include "ObitDConCleanVisWB.h"
#include "ObitMem.h"
#include "ObitFFT.h"
#include "ObitTableUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitSkyGeom.h"
#include "ObitDConCleanPxListWB.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanVisWB.c
 * ObitDConCleanVisWB class function definitions.
 * Image based CLEAN class.
 * This class is derived from the ObitDCon class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanVisWB";

/** Function to obtain parent ClassInfo - ObitDConClean */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanVisGetClass;

/**
 * ClassInfo structure ObitDConCleanVisWBClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanVisWBClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private structures----------------*/
/* Image model subtraction threaded function argument */
typedef struct {
  /* CLEAN Object */
  ObitDConCleanVisWB *in;
  /* Input plane pixel data, flux, si, curve..., nterm of these */
  ObitFArray **inData;
  /* Number of spectral terms, 0=flux, 1+=si, 2+=curve  */
  olong      nterm;
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
void  ObitDConCleanVisWBInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanVisWBClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConCleanVisWBClassInfoDefFn (gpointer inClass);

/** Private: Set Beam patch, min. flux, decide on SDI CLEAN. */
static void ObitDConCleanVisWBDecide (ObitDConCleanVis* in, ObitErr *err);

/** Private: (re)make residuals. */
static void  MakeResiduals (ObitDConCleanVis *in, olong *fields, 
			    gboolean doBeam, ObitErr *err);

/** Private: (re)make all residuals. */
static void  MakeAllResiduals (ObitDConCleanVis *in, ObitErr *err);

/** Private: Threaded Image subtractor */
static gpointer ThreadImSub (gpointer arg);

/** Private: Make Threaded Image subtraction args */
static olong MakeImSubFuncArgs (ObitThread *thread, olong nterm,
				ObitErr *err, ImSubFuncArg ***ThreadArgs);

/** Private: Delete Threaded Image subtraction args */
static void KillImSubFuncArgs (olong nargs, ImSubFuncArg **ThreadArgs);

/** Private: Low accuracy subtract CLEAN model. */
static void SubNewCCs (ObitDConCleanVis *in, olong *newCC, 
		       ObitFArray **pixarray, ObitErr *err);

/** Private: Create/init PxList. */
static void NewPxList (ObitDConCleanVis *in, ObitErr *err);

/** Private: Create/init Pixarray. */
static ObitFArray** NewPxArray (ObitDConCleanVis *in, olong *startCC, 
				ObitErr *err);

/** Private: Delete Pixarray. */
static ObitFArray** KillPxArray (ObitDConCleanVis *in, ObitFArray **pixarray);

/** Private: Delete BeamPatches. */
static void KillBeamPatches (ObitDConCleanVis *in);

/** Private: Convolve spectral CCs with a Gaussian. */
static ObitFArray* ConvlCC(ObitImage *image, olong CCVer, olong iterm, 
			   ObitErr *err);

/** Private: Cross convolve spectral CCs with a Gaussian. */
static void XConvlCC(ObitImage *in, olong CCVer, olong iterm, 
		     ObitImage *out, ObitFArray *outGrid, ObitErr *err);

/** Private: Apply Gaussian taper to uv grid. */
static void GaussTaper (ObitCArray* uvGrid,  ObitImageDesc *imDesc,
			ofloat gparm[3]);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanVisWB* newObitDConCleanVisWB (gchar* name)
{
  ObitDConCleanVisWB* out;
  /*gchar *routine = "newObitDConCleanVisWB";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisWBClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanVisWB));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanVisWBInit((gpointer)out);

 return out;
} /* end newObitDConCleanVisWB */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanVisWBGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisWBClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanVisWBGetClass */

/**
 * Make a deep copy of an ObitDConCleanVisWB.
 * \param inn  The object to copy
 * \param outt An existing object pointer for output or NULL if none exists.
 * \param err  Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanVis* ObitDConCleanVisWBCopy  (ObitDConCleanVis *inn, 
					     ObitDConCleanVis *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  ObitDConCleanVisWB *in  = (ObitDConCleanVisWB*)inn;
  ObitDConCleanVisWB *out = (ObitDConCleanVisWB*)outt;
  gchar *routine = "ObitDConCleanVisWBCopy";

  /* error checks */
  if (err->error) return (ObitDConCleanVis*)out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitDConCleanVisWB(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, (ObitDConCleanVis*)out);

  /*  copy this class */
  return (ObitDConCleanVis*)out;
} /* end ObitDConCleanVisWBCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanVisWB similar to the input one.
 * \param inn  The object to copy
 * \param outt An existing object pointer for output, must be defined.
 * \param err  Obit error stack object.
 */
void ObitDConCleanVisWBClone  (ObitDConCleanVis *inn, ObitDConCleanVis *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitDConCleanVisWB *in  = (ObitDConCleanVisWB*)inn;
  ObitDConCleanVisWB *out = (ObitDConCleanVisWB*)outt;
  gchar *routine = "ObitDConCleanVisWBClone";

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
} /* end ObitDConCleanVisWBClone */

/**
 * Creates an ObitDConCleanVisWB 
 * defined for convenience of derived classes 
 * \param name   An optional name for the object.
 * \param uvdata from which to create object, should have all control
                 information defined on info member.
 * \param order  Order of the imaging, Spectral index only=1, plus curvature=2
 * \param err    Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVisWB* ObitDConCleanVisWBCreate (gchar* name, ObitUV *uvdata,  
					      olong order, ObitErr *err)
{
  olong nfield, i;
  ObitDConCleanVisWB* out=NULL;
  gchar *routine = "ObitDConCleanVisWBCreate";

 /* error checks */
  if (err->error) return out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVisWB (name);

  out->order = order;  /* Save order */

  /* Create UV imager and its ImageMosaic */
  out->imager = (ObitUVImager*)ObitUVImagerWBCreate("UVImagerWB", order, uvdata, err);
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
} /* end ObitDConCleanVisWBCreate */

/**
 * Creates an ObitDConCleanVisWB from optional components 
 * defined for convenience of derived classes 
 * \param name     An optional name for the object.
 * \param uvdata   from which to create object, should have all control
                   information defined on info member.
 * \param imager   Optional ObitUVImager to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
                   Should be an ObitUVImagerWB
 * \param skyModel Optional ObitSkyModel to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
 * \param order    Order of the imaging, Spectral index only=1, plus curvature=2
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVis* 
ObitDConCleanVisWBCreate2 (gchar* name, ObitUV *uvdata,  
			   ObitUVImager *imager, ObitSkyModel *skyModel, 
			   olong order, ObitErr *err)
{
  olong nfield, i;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitDConCleanVisWB* out=NULL;
  ofloat ftemp;
  gchar *routine = "ObitDConCleanVisWBCreate";

 /* error checks */
  if (err->error) return (ObitDConCleanVis*)out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVisWB (name);

  out->order = order;  /* Save order */

  /* Use or create UV imager and create its ImageMosaic */
  if (imager==NULL) {
    out->imager =(ObitUVImager*) ObitUVImagerWBCreate("UVImager", order, uvdata, err);
    if (err->error) Obit_traceback_val (err, routine, name, (ObitDConCleanVis*)out);
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
} /* end ObitDConCleanVisWBCreate2 */

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
void  ObitDConCleanVisWBGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConClassInfo *ParentClass;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  /*olong i;*/
  /*union ObitInfoListEquiv InfoReal;*/
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;  /* as this class */
  gchar *routine = "ObitDConCleanVisWBGetParms";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* minimum residual flux for 1st and second order */
  ObitInfoListGetTest(in->info, "OrdFlux", &type, dim, &in->OrdFlux);

  /* Image display? */
} /* end ObitDConCleanVisWBGetParms */

/**
 * Set default CLEAN windows in mosaic
 * If mosaic member  Radius>0 then make round boxes on Fly's eye field
 * with this radius, else use rectangular box including all but outer 5 pixels
 * On outlier fields, use rectangular box of width OutlierSize.
 * If CLEANBox defined in in->info then its contents are used for field 1.
 * Assumes all images in mosaic have descriptors defined.
 * Basic setup uses parent class function.
 * Then for Wideband SW imaging, the outer 3/4 of each image is covered 
 * with an upbox. 
 * \param in   The CLEAN object
 * \param err Obit error stack object.
 */
void ObitDConCleanVisWBDefWindow(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisWB *in;
  const ObitDConCleanClassInfo *parentClass;
  ObitDConCleanWindowType type;
  olong field, *window, win[4], nx, ny;
  gchar *routine = "ObitDConCleanVisWBDefWindow";

  /* Cast input to this type */
  in = (ObitDConCleanVisWB*)inn;
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));
  
  parentClass = myClassInfo.ParentClass;
  
  /* Call actual function */
  parentClass->ObitDConCleanDefWindow((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Make unwindows on outer 3/4 - make outer window only hold inner 1/4 */
  if (in->mosaic && (in->mosaic->numberImages>0)) {
    for (field=0; field<in->mosaic->numberImages; field++) {
      nx = in->mosaic->images[field]->myDesc->inaxes[0];
      ny = in->mosaic->images[field]->myDesc->inaxes[1];
      /* Reset any outer window */
      if (ObitDConCleanWindowInfo(in->window, field+1, -1, &type, &window, err)) {
	switch (type) {
	case OBIT_DConCleanWindow_rectangle:
	case OBIT_DConCleanWindow_unrectangle:
	  window[0] = MAX((nx/4)+1,   window[0]); 
	  window[2] = MIN((3*nx/4)-1, window[2]);
	  window[1] = MAX((nx/4)+1,   window[1]); 
	  window[3] = MIN((3*nx/4)-1, window[3]);
	  break;
	case OBIT_DConCleanWindow_round:
	case OBIT_DConCleanWindow_unround:
	  window[0] = MIN (nx/4, window[0]);
	  break;
	default: /* bad juju Bwana */
	  Obit_log_error(err, OBIT_InfoWarn,"%s: Unknown window type %d",
			 routine, type);
	  return;
	}; /* end switch by window type */
	/* Update */
	ObitDConCleanWindowUpdate(in->window, field+1, -1, type, window, err);
      }
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      /* Add unboxes */
      type = OBIT_DConCleanWindow_unrectangle;
      win[0] = 1;      win[1] = 1;      win[2] = nx;   win[3] = ny/4;
      ObitDConCleanWindowAdd (in->window, field+1, type, win, err);
      win[0] = 1;      win[1] = ny/4;    win[2] = nx/4; win[3] = 3*ny/4;
      ObitDConCleanWindowAdd (in->window, field+1, type, win, err);
      win[0] = 3*nx/4; win[1] = ny/4;    win[2] = nx;    win[3]= 3*ny/4;
      ObitDConCleanWindowAdd (in->window, field+1, type, win, err);
      win[0] = 1;      win[1] = 3*ny/4; win[2] = nx;    win[3] = ny;
      ObitDConCleanWindowAdd (in->window, field+1, type, win, err);			  
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }/* end loop over fields */
  } /* end if mosaic */
} /* end ObitDConCleanVisWBDefWindow */

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
 * \param inn      The object to restore
 * \param fields   Field numbers (1-rel) in ImageMosaic
 *                 zero terminated list, no longer than in->numCurrentField
 * \param pixarray If nonNULL, then use this array of pixels rather than in the ImageMosaic
 *                 Elements are Pixel arrays corresponging to fields in fields, 
 *                 in->order+1 per field
 * \param err      Obit error stack object.
 * \return TRUE if clean should be continued
 */
gboolean ObitDConCleanVisWBAutoWindow(ObitDConClean *inn, olong *fields, ObitFArray **pixarray, 
				      ObitErr *err)
{
  ObitDConCleanVisWB *in=(ObitDConCleanVisWB*)inn;
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
  gchar *routine = "ObitDConCleanVisWBAutoWindow";

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanVisWBIsA(in));

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
    if ((pixarray==NULL) || (pixarray[field]==NULL)) {
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
    in->imgRMS[fields[field]-1]     = RMS;
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
} /* end ObitDConCleanVisWBAutoWindow */

/**
 * Get image and beam statistics and prepare to deconvolve
 * If autoWindow option is selected then windows may be added 
 * \param inn       The object to deconvolve
 * \param pixarray  If NonNULL use instead of the flux densities from the image file.
 *                  This is an array of ObitFarrays corresponding to fields in
                    in->currentFields, in->order+1 per field
 * \param err       Obit error stack object.
 * \return TRUE if CLEAN should continue
 */
gboolean ObitDConCleanVisWBPixelStats(ObitDConClean *inn, ObitFArray **pixarray, 
				      ObitErr *err)
{
  ObitDConCleanVisWB *in=(ObitDConCleanVisWB*)inn;
  gboolean newWin = FALSE, isLast = FALSE, allDone=TRUE;
  olong field, ifld;
  ObitImage *Beam=NULL;
  const ObitDConCleanClassInfo *inClass;
  gchar *routine = "ObitDConCleanVisWBPixelStats";

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanVisWBIsA(in));

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
    ObitDConCleanVisWBDecide ((ObitDConCleanVis*)in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
  }

  return newWin || (!allDone);
} /* end ObitDConCleanVisWBPixelStats */

/**
 * Determine Beam patch size and minimum flux density to consider
 * such that the minimum flux density is the maximum ignored sidelobe
 * of the brightest pixel value.
 * Wideband version, beampatch = size of beam
 * Output members are beamPatchSize, minFluxLoad and depend on the members
 * currentFields, window, BeamHist and PixelHist being up to date.
 * No SDI CLEAN this class
 * If the histogram is too coarse (too many in top bin) then numberSkip
 * is set to decimate the residuals so that they will fit;
 * \param inn   The object to deconvolve
 * \param err  Obit error stack object.
 */
static void ObitDConCleanVisWBDecide (ObitDConCleanVis* inn, ObitErr *err)
{
  olong i;
  ofloat minFlux=0.0;
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  gchar *routine = "ObitDConCleanVisWBDecide";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));

  /* initialize values */
  in->beamPatchSize = in->BeamHist->ncell;  /* Use this one */
  in->minFluxLoad = 0.0;
  in->numberSkip = 0;

  /* Loop through pixel histogram for lowest flux limit that fits */
  for (i=0; i<in->PixelHist->ncell-1;i++) {
    if (in->PixelHist->hist[i]<in->maxPixel) { /* go team! */
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

  /* SDI CLEAN makes no sense here */
 doSDI:
  in->doSDI = FALSE;
  /* Give warning if skipping */
  if ((in->numberSkip>=1) && (in->prtLv>0)) 
    Obit_log_error(err, OBIT_InfoWarn,"%s: Too many residuals, taking 1 of every  %d",
		   routine, in->numberSkip+1);

} /* end  ObitDConCleanVisWBDecide */

/**
 * Select components to be subtracted
 * Supports multiple fields
 * For ObitDConCleanVisWB, the BeamPatches member is not used, 
 * the equivalent is obtained from the Images directly with no 
 * control over the patch size.
 * \param inn        The object to deconvolve
 * \param pixarray   If NonNULL use instead of the flux densities from the image file.
 *                   Array of ObitFArrays corresponding to fields in in->currentFields,
 *                   in->order+1 per field
 * \param err        Obit error stack object.
 * \return TRUE if deconvolution is complete
 */
gboolean ObitDConCleanVisWBSelect(ObitDConClean *inn, ObitFArray **pixarray, 
				  ObitErr *err)
{
  gboolean done = FALSE;
  olong field, ifld;
  const ObitDConCleanPxListClassInfo *pxListClass;
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  gchar *routine = "ObitDConCleanSelect";

  /* error checks */
  if (err->error) return done;
  g_assert (ObitDConCleanVisWBIsA(in));

  /* Get SW arrays on images if pixarray not specified */
  if (!pixarray) {
    ifld = 0;
    field = in->currentFields[ifld];
    /* Loop over fields in in->currentFields */
    while (field>0) {
      ObitImageWBMakeWork ((ObitImageWB*)in->mosaic->images[field-1], err);
      ifld++;
      field = in->currentFields[ifld];
    } /* end loop over fields */
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
  }

  /* Read Pixel List */
  pxListClass = (ObitDConCleanPxListClassInfo*)in->Pixels->ClassInfo; 
  pxListClass->ObitDConCleanPxListUpdate (in->Pixels, in->currentFields, in->numberSkip,
					  in->minFluxLoad, in->autoWinFlux, 
					  in->window, in->BeamPatches, pixarray, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* Only BGC is sensible */
  done = pxListClass->ObitDConCleanPxListCLEAN (in->Pixels, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  return done;
} /* end ObitDConCleanVisWBSelect */

/** 
 * Decompose images into spectral components (NOT a deconvolution)
 * Wideband imaging version.
 * Writes  all spectral planes in images on mosaic replacing original residual
 * Note, this is also done in MakeResiduals and MakeAllResiduals
 * and should not be called a second time.
 * \param inn  The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanVisWBDecomp(ObitDConCleanVisWB *in, ObitErr *err)
{
  olong field;
  gchar *routine = "ObitDConCleanVisWBDecomp";

   /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));

  /* Loop over fields decomposing */
  for (field = 0; field<in->nfield; field++) {
    /* Decompose image into spectral components */
    if (ObitImageWBIsA(in->mosaic->images[field])) 
      ObitImageWBDecomp ((ObitImageWB*)in->mosaic->images[field], err);
  } /* end loop decomposing */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

}

/** 
 * Restore components removed from the residual image(s)
 * Wideband imaging version.
 * Spectral orders higher than 0 are flux density weighted averages.
 * \param inn  The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanVisWBRestore(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  ObitFArray *convl=NULL, *flux=NULL;
  ObitImage *image=NULL;
  olong iplane, field, num, plane[5]={1,1,1,1,1};
  gchar *routine = "ObitDConCleanVisWBRestore";

   /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));

  /* Anything to restore? */
  if (in->Pixels->currentIter<=0) return;

  /* Tell user */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Restoring components");
    ObitErrLog(err);  /* Progress Report */
  }

  /* Multiply spectral planes by flux plane */
  /* Loop over fields */
  for (field = 0; field<in->nfield; field++) {
    /* which Image? */
    image = in->mosaic->images[field];
    /* How many planes? */
    if (ObitImageWBIsA(image))
      num = MAX(1,(((ObitImageWB*)image)->order+1));
    else
      num = 1;
    if (num>1) {
      /* Read flux image */
      plane[0] = 1;
      ObitImageGetPlane (image, NULL, plane, err);
      flux = ObitFArrayCopy(image->image, flux, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      for (iplane=1; iplane<num; iplane++) {
	plane[0] = iplane+1;
	ObitImageGetPlane (image, NULL, plane, err);
	/* Multiply */
	ObitFArrayMul (image->image, flux, image->image);
	/*  DEBUG 
	if (iplane==1) ObitImageUtilArray2Image ("DbugRestoreSI.fits", 0, image->image, err);*/
	/* Rewrite */
	ObitImagePutPlane (image, NULL, plane, err);
      }
      flux = ObitFArrayUnref(flux);
    } /* End multiple spectral planes */
  }  /* end loop over fields */

  /* Loop over fields */
  for (field = 0; field<in->nfield; field++) {

    /* Anything to restore? */
    if (in->Pixels->iterField[field]<=0) continue;

    /* which Image? */
    image = in->mosaic->images[field];

    /* Loop over spectral planes - Do higher order by flux weighted average */
    if (ObitImageWBIsA(image))
      num = MAX(1,(((ObitImageWB*)image)->order+1));
    else
      num = 1;
    for (iplane=0; iplane<num; iplane++) {
      /* Convolve Gaussians */
      convl = ConvlCC (image, in->CCver, iplane, err);
      /* Read image */
      plane[0] = iplane+1;
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

} /* end ObitDConCleanVisWBRestore */

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
void ObitDConCleanVisWBXRestore(ObitDConClean *inn, ObitErr *err)
{
 ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  ObitImage *image1=NULL, *image2=NULL;
  olong ifield, jfield, num, iplane, plane[5]={1,1,1,1,1};
  gchar *routine = "ObitDConCleanVisWBXRestore";

   /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));

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
 
    /* Loop over spectral planes - Do higher order by flux weighted average */
    if (ObitImageWBIsA(image1))
      num = MAX(1,(((ObitImageWB*)image1)->order+1));
    else
      num = 1;
    for (iplane=0; iplane<num; iplane++) {
      /* Read image */
      plane[0] = iplane+1;
      ObitImageGetPlane (image1, NULL, plane, err);

      /* Loop others */
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
} /* end ObitDConCleanVisWBXRestore */

/** 
 * Normalize spectral channels
 * Expects spectral planes to be the convolution of flux*spectral term 
 * divides flux plane clips at sigma into spectral planes
 * Adds inn->info entry "Alpha" to spectral index image if not 0.0
 * \param inn  The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanVisWBSpecNorm(ObitDConCleanVisWB *in, ObitErr *err)
{
  olong field, num, iplane, plane[5]={1,1,1,1,1};
  ObitImage *image;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitFArray *flux=NULL;
  ofloat maxSpec[5] = {1.0e20,3.0, 1.0, 1.0, 1.0};
  ofloat sigma, alpha, fblank = ObitMagicF();
  /*gchar *routine = "ObitDConCleanVisWBSpecNorm";*/

  /* Spectral index to add? */
  alpha = 0.0;
  ObitInfoListGetTest(in->info, "Alpha", &type, dim, &alpha);

  /* Loop over fields normalizing */
  for (field = 0; field<in->nfield; field++) {
    /* which Image? */
    image = in->mosaic->images[field];
    /* Ignore non ImageWB fields */
      /* Read flux */
      plane[0] = 1;
      ObitImageGetPlane (image, NULL, plane, err);
      flux = ObitFArrayCopy(image->image, flux, err);
      /* blank below 1 sigma */
      sigma = ObitFArrayRMS (flux);
      ObitFArrayClip (flux, 3.0*sigma, 1.0e20, fblank);

      /* Loop over spectral planes */
      num = MAX(1,(((ObitImageWB*)image)->order+1));
     for (iplane=1; iplane<num; iplane++) {
       plane[0] = 1+iplane;
      ObitImageGetPlane (image, NULL, plane, err);
      /* normalize */
      ObitFArrayDiv (image->image, flux, image->image);
      /* Clip outside +/- maxSpec[iplane] */
      ObitFArrayClip (image->image, -maxSpec[iplane], maxSpec[iplane], fblank);
      /* Add alpha to si if non-zero */
      if ((iplane==1) && (alpha!=0.0)) ObitFArraySAdd (image->image, alpha);
      /* Rewrite */
      ObitImagePutPlane (image, NULL, plane, err);
     } /* end loop over planes */
  
    if (ObitImageWBIsA(image)) {
    } /* end if ImageWB */

  } /* end loop over fields */
} /* end  ObitDConCleanVisWBSpecNorm */

/**
 * Flatten multiple facets if needed
 * Wideband imaging version.
 * Flattens all Spectral planes
 * Does Flatten if FullField member of mosaic member is defined.
 * \param inn  The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanVisWBFlatten(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  ObitImageMosaicClassInfo* mosaicClass; 
  gchar *routine = "ObitDConCleanVisWBFlatten";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));

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

} /* end ObitDConCleanVisWBFlatten */

/**
 * Subtract components from uv data.
 * Frees any work arrays on mosaic images
 * \param inn  The object to deconvolve
 * \param err  Obit error stack object.
 */
void ObitDConCleanVisWBSub(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  const ObitDConCleanVisClassInfo *parentClass;
  olong ifld;
  gchar *routine = "ObitDConCleanVisWBSub";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));

  /* Clean up image spectral work arrays */
  for (ifld=0; ifld<in->nfield; ifld++) {
    ObitImageWBClearWork ((ObitImageWB*)in->mosaic->images[ifld]);
  }
  
  /* Most work in parent class */
  parentClass = myClassInfo.ParentClass;
  parentClass->ObitDConCleanSub((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitDConCleanVisWBSub */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanVisWBClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanVisWBClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanVisWBClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanVisWBClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanVisWBClassInfo *theClass = (ObitDConCleanVisWBClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanVisWBClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanVisWBClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanVisWBGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanVisWB;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanVisWBCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanVisWBClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanVisWBInit;
  theClass->ObitDConGetParms          = (ObitDConGetParmsFP)ObitDConCleanVisWBGetParms;
  theClass->ObitDConCleanSub          = (ObitDConCleanSubFP)ObitDConCleanVisWBSub;
  /*theClass->ObitDConDeconvolve        = (ObitDConDeconvolveFP)ObitDConCleanVisWBDeconvolve;*/
  /*theClass->ObitDConCleanVisQuality   = (ObitDConCleanVisQualityFP)ObitDConCleanVisWBQuality;*/
  /*theClass->ObitDConCleanVisReimage   = (ObitDConCleanVisReimageFP)ObitDConCleanVisWBReimage;*/
  /*theClass->ObitDConCleanVisAddField  = (ObitDConCleanVisAddFieldFP)ObitDConCleanVisWBAddField;*/
  /*theClass->ObitDConCleanVisRecenter  = (ObitDConCleanVisRecenterFP)ObitDConCleanVisWBRecenter;*/
  /*theClass->ObitDConCleanVisFilter    = (ObitDConCleanVisFilterFP)ObitDConCleanVisWBFilter;*/
  theClass->ObitDConCleanSelect       = (ObitDConCleanSelectFP)ObitDConCleanVisWBSelect;
  theClass->ObitDConCleanVisDefWindow = (ObitDConCleanVisDefWindowFP)ObitDConCleanVisWBDefWindow;
  theClass->ObitDConCleanAutoWindow = 
    (ObitDConCleanAutoWindowFP)ObitDConCleanVisWBAutoWindow;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanVisWBPixelStats;
  theClass->ObitDConCleanRestore = (ObitDConCleanRestoreFP)ObitDConCleanVisWBRestore;
  theClass->ObitDConCleanFlatten = (ObitDConCleanFlattenFP)ObitDConCleanVisWBFlatten;
  theClass->ObitDConCleanXRestore= (ObitDConCleanXRestoreFP)ObitDConCleanVisWBXRestore;
  theClass->ObitDConCleanVisWBDecomp= (ObitDConCleanVisWBDecompFP)ObitDConCleanVisWBDecomp;

  /* Private functions definitions for derived classes */
  theClass->MakeResiduals   = (MakeResidualsFP)MakeResiduals;
  theClass->MakeAllResiduals= (MakeAllResidualsFP)MakeAllResiduals;
  theClass->SubNewCCs       = (SubNewCCsFP)SubNewCCs;
  theClass->NewPxList       = (NewPxListFP)NewPxList;
  theClass->NewPxArray      = (NewPxArrayFP)NewPxArray;
  theClass->KillPxArray     = (KillPxArrayFP)KillPxArray;
  theClass->KillBeamPatches = (KillBeamPatchesFP)KillBeamPatches;

} /* end ObitDConCleanVisWBClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanVisWBInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVisWB *in = inn;

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
} /* end ObitDConCleanVisWBInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanVisWB* cast to an Obit*.
 */
void ObitDConCleanVisWBClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVisWB *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanVisWBClear */

/**
 * Make selected residual images and get statistics, 
 * Also rebuilds image spectral work arrays and
 * replaces residual with SW spectral estimates
 * \param inn    The Clean object
 * \param fields zero terminated list of field numbers to image
 * \param doBeam If TRUE also make beam
 * \param err    Obit error stack object.
 */
static void  MakeResiduals (ObitDConCleanVis *inn, olong *fields, 
			    gboolean doBeam, ObitErr *err)
{
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  const ObitDConCleanVisClassInfo *parentClass;
  olong ifld, field;
  gchar *routine = "ObitDConCleanVisWB:MakeResiduals";
  
  if (err->error) return; /* prior error condition? */

  /* Clean up all stale image spectral work arrays */
  for (ifld=0; ifld<in->nfield; ifld++) {
    ObitImageWBClearWork ((ObitImageWB*)in->mosaic->images[ifld]);
  }
  
  parentClass = myClassInfo.ParentClass;
  /* Call MakeResiduals in parent class */
  parentClass->MakeResiduals (inn, fields, doBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Mark images as fresh */
  for (ifld=0; ifld<in->nfield; ifld++) {
    field = fields[ifld];
    if (field<=0) break;
    ((ObitImageWB*)in->mosaic->images[field-1])->fresh = TRUE;
    /* Set imaging order */
    if (fabs(in->cleanable[field-1]) < in->OrdFlux[1]) 
      ((ObitImageWB*)in->mosaic->images[field-1])->curOrder = 1;
      if (fabs(in->cleanable[field-1]) < in->OrdFlux[0]) 
      ((ObitImageWB*)in->mosaic->images[field-1])->curOrder = 0;
  }
  
} /* end MakeResiduals */

/**
 * Make all residual images and get statistics
 * Also rebuilds image spectral work arrays and
 * replaces residual with SW spectral estimates
 * \param inn    The Clean object
 * \param err    Obit error stack object.
 */
static void  MakeAllResiduals (ObitDConCleanVis *inn, ObitErr *err)
{
 ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  const ObitDConCleanVisClassInfo *parentClass;
  olong ifld;
  gchar *routine = "ObitDConCleanVisWB:MakeAll Residuals";
 
  if (err->error) return; /* prior error condition? */

  /* Clean up all image spectral work arrays */
  for (ifld=0; ifld<in->nfield; ifld++) {
    ObitImageWBClearWork ((ObitImageWB*)in->mosaic->images[ifld]);
  }
  
  parentClass = myClassInfo.ParentClass;
  /* Call MakeAllResiduals in parent class */
  parentClass->MakeAllResiduals (inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Mark images as fresh */
  for (ifld=0; ifld<in->nfield; ifld++) {
    ((ObitImageWB*)in->mosaic->images[ifld])->fresh = TRUE;
    /* Set imaging order */
    if (fabs(in->cleanable[ifld]) < in->OrdFlux[1])
      ((ObitImageWB*)in->mosaic->images[ifld])->curOrder = 1;
    if (fabs(in->cleanable[ifld]) < in->OrdFlux[0]) 
      ((ObitImageWB*)in->mosaic->images[ifld])->curOrder = 0;
  }
  
} /* end MakeAllResiduals */

/**
 * Create Pixel list for cleaning
 * \param inn      The Clean object
 * \param err      Obit error stack object.
 */
static void NewPxList (ObitDConCleanVis *inn, ObitErr *err)
{
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  gchar *routine = "ObitDConCleanVis:NewPxList";

  if (in->Pixels && (in->Pixels->nfield!=in->mosaic->numberImages)) 
    in->Pixels = ObitDConCleanPxListUnref(in->Pixels);
  if (!in->Pixels) {
    in->Pixels = 
      (ObitDConCleanPxList*)ObitDConCleanPxListWBCreate("Pixel List", in->mosaic, 
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
 * \param inn      The Clean object
 * \param newCC    Start CC -1 per in->mosaic field for subtraction
 * \param pixarray Array of arrays of pixels to subtract, in->order+1 per field
 * \param err      Obit error stack object.
 * \return TRUE if attempted, FALSE if cannot do 
 */
static void SubNewCCs (ObitDConCleanVis *inn, olong *newCC, ObitFArray **pixarray, 
		       ObitErr *err)
{
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL;
  ImSubFuncArg **threadArgs;
  ObitFArray **comps=NULL;
  ObitImageDesc *outDesc;
  ofloat PeakIn, PeakOut, RMS,  parms[20];
  olong i, j, k, ip, ncc, ver, nThreads=0, nTh, nfield, ifield, nDo, nLeft, PeakInPos[2];
  gboolean doAbs, OK;
  gchar *tabType = "AIPS CC";
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  gchar *routine = "SubNewCCs";
  gboolean DebugGDB=FALSE;  /*  DEBUG */

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisWBIsA(in));

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
    comps[i]   = ObitTableCCUtilMergeSelSpec (CCTable, newCC[i]+1, ncc, parms, err);
    CCTable = ObitTableUnref(CCTable);
    Obit_return_if_fail ((comps[i]!=NULL),
			 err, "%s: Error merging CLEAN components for field %d",
			 routine, ifield);

  } /* end loop over fields */

  /* setup for threaded processing
     Initialize Threading */
  nThreads = MakeImSubFuncArgs (in->thread, in->order, err, &threadArgs);
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
  ip = 0;
  for (i=0; i<nfield; i+=nTh) {
    nDo = MIN (nTh, nLeft);
    /* Args for this batch */
    for (j=0; j<nDo; j++) {
      threadArgs[j]->ofield = in->currentFields[i+j];
      for (k=0; k<=in->order; k++)
	threadArgs[j]->inData[k] = pixarray[ip++];
    }
    /* Dummy any remainders */
    for (j=nDo; j<nTh; j++) {
      threadArgs[j]->ofield = 1;
      for (k=0; k<=in->order; k++)
	threadArgs[j]->inData[k] = NULL;
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
    ip += in->order+1;
  } /* end statistics loop over fields */
  
} /* end SubNewCCs */

/**
 * Create Pixel array for intermediate CLEAN
 * \param inn      The Clean object
 * \param startCC  [out] Current number of components per field
 * \param err      Obit error stack object.
 * \return array of ObitFArrays for cleaning, in->order+1 per field
 */
static ObitFArray** NewPxArray (ObitDConCleanVis *inn, olong *startCC, 
				ObitErr *err)
{
  ObitFArray** pixarray=NULL;
  ObitImageWB *image;
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  olong ifld, field, num, ip=0;
  gchar *routine = "ObitDConCleanVisWB:NewPxArray";

  if (err->error) return pixarray;

  /* 1/field for order 0, 2 for order 1, 3 for order 2 */
  num = in->numCurrentField*(1+in->order);
  pixarray  = g_malloc0(num*sizeof(ObitFArray*));
  for (ifld=0; ifld<in->numCurrentField; ifld++) {
    if (in->currentFields[ifld]<=0) break;  /* List terminated? */
    field = in->currentFields[ifld];
    image = (ObitImageWB*)in->mosaic->images[field-1];
    if (!image->ResidArr[0]) ObitImageWBMakeWork (image, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, pixarray);
    startCC[ifld]   = in->Pixels->iterField[field-1];  /* no. at start */
    /* Pixel array for field */
    pixarray[ip] = ObitFArrayRef(image->ResidArr[0]); ip++;
    if (in->order>0) {
      pixarray[ip] = ObitFArrayRef(image->ResidArr[1]); ip++;
    }
    if (in->order>1) {
      pixarray[ip] = ObitFArrayRef(image->ResidArr[2]); ip++;
    }
    /* DEBUG */
    /* ObitImageUtilArray2Image ("DbugConvlResid.fits", 0, image->ResidArr[0], err);  */
    /* ObitImageUtilArray2Image ("DbugBeam00.fits", 0, image->BeamMatx[0][0], err);  */
    /* if (err->error) Obit_traceback_val (err, routine, image->name, pixarray); */
    /* END DEBUG */

  } /* end loop over fields */
  
  return pixarray;
} /* end NewPxArray */

/**
 * Delete Pixel array for intermediate CLEAN
 * \param inn      The Clean object
 * \param pixarray Array to delete, in->order+1 per field
 * \return NULL
 */
static ObitFArray** KillPxArray (ObitDConCleanVis *inn,  ObitFArray **pixarray)
{
  ObitDConCleanVisWB *in = (ObitDConCleanVisWB*)inn;
  olong ifld, ip = 0;

  /* 1/field for order 0, 2 for order 1, 2 for order 2 */
  if (pixarray) {
    for (ifld=0; ifld<in->numCurrentField; ifld++) {
      if (in->currentFields[ifld]<=0) break;  /* List terminated? */
      
      pixarray[ip] = ObitFArrayUnref(pixarray[ip]); ip++;
      if (in->order>0) {pixarray[ip] = ObitFArrayUnref(pixarray[ip]); ip++;}
      if (in->order>1) {pixarray[ip] = ObitFArrayUnref(pixarray[ip]); ip++;}
    }
    g_free(pixarray);  pixarray = NULL;
  }
  return pixarray;
} /* end KillPxArray */

/**
 * Delete Beam Patches on in
 * NOP for ObitDConCleanVisWB, the beam patches on images are used instead
 * \param in       The Clean object
 */
static void KillBeamPatches (ObitDConCleanVis *in)
{
  /* Nothing*/
} /* end KillBeamPatches */

/**
 * Subtract a set of lists of CLEAN components from an image in a thread
 * Wideband version, use beams on PxList
 * \param arg Pointer to ImSubFuncArg argument with elements:
 * \li in       The ObitDConCleanVisWB object
 * \li inData   Pixel arrays to be subtracted from,
 *              [0]=>flux, [1]=si, [2]=curve
 * \li nterm    Number of spectral terms, 0=flux, 1+=si, 2+=curve 
 * \li ofield   Field number (1-rel) of data in inData 
 * \li inComp   Array of compressed arrays of point pixels
 *              Correspond to fields in->currentFields
 *              Spectral terms follow flux, X, Y.
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
  ObitDConCleanVisWB *in= largs->in;
  ObitFArray **pixarray = largs->inData;
  olong       nterm     = largs->nterm;
  olong       ofield    = largs->ofield;
  ObitFArray **inComp   = largs->inComp;
  olong       nfield    = largs->nfield;
  ObitErr    *err       = largs->err;
  ObitThread *thread    = largs->thread;
  /* local */
  olong i, j, ifield, len, pos1[2], pos2[2];
  ofloat idx, idy, ftemp, flux=0.0, si=0.0, curve=0.0, offset[2];
  ObitFArray *comps;
  ObitImageDesc *inDesc, *outDesc;
  ObitDConCleanPxListWB *pxList = (ObitDConCleanPxListWB*)in->Pixels;
  ObitFArray *beam00=NULL, *beam10=NULL, *beam01=NULL, *beam11=NULL;
  ObitFArray *beam02=NULL, *beam20=NULL, *beam21=NULL, *beam12=NULL, *beam22=NULL;

  if (err->error) goto finish;
  if (pixarray[0]==NULL) goto finish;

  /* Loop over fields */
  for (i=0; i<nfield; i++) {
    ifield = in->currentFields[i];
    if (ifield<=0) break;
    comps = inComp[i];          /* Component list */
    if (comps==NULL) continue;  /* Any components? */
    beam00  = pxList->BeamPatch00[ifield-1]; /* Beam patch pointer */
    if (nterm>0) { /* first order */
      beam01  = pxList->BeamPatch01[ifield-1];
      beam10  = pxList->BeamPatch10[ifield-1];
      beam11  = pxList->BeamPatch11[ifield-1];
      if (nterm>1) { /* second order */
	beam20  = pxList->BeamPatch20[ifield-1];
	beam02  = pxList->BeamPatch02[ifield-1];
	beam12  = pxList->BeamPatch12[ifield-1];
	beam21  = pxList->BeamPatch21[ifield-1];
	beam22  = pxList->BeamPatch22[ifield-1];
      } 
    } 
    
    /* Loop over CCs - assumes all points */
    inDesc    = in->mosaic->images[ifield-1]->myDesc;
    outDesc   = in->mosaic->images[ofield-1]->myDesc;
    len       = comps->naxis[0];
    pos2[0]   = pxList->BeamPatch00[ifield-1]->naxis[0]/2; 
    pos2[1]   = pxList->BeamPatch00[ifield-1]->naxis[1]/2;    /* Center of beam patch */
    offset[0] = outDesc->xPxOff + pixarray[0]->naxis[0]/2; 
    offset[1] = outDesc->yPxOff + pixarray[0]->naxis[1]/2;   /* Center of image */
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
      if (nterm>0) si    = comps->array[j*len+3];  /* Spectral index */
      if (nterm>1) curve = comps->array[j*len+3];  /* Spectral curvature */
      /* Zeroth order */
      ObitFArrayShiftAdd(pixarray[0], pos1, beam00, pos2,   -flux,  pixarray[0]);
      if (nterm>0) 
	ObitFArrayShiftAdd(pixarray[0], pos1, beam10, pos2, -si,    pixarray[0]);
      if (nterm>1) 
	ObitFArrayShiftAdd(pixarray[0], pos1, beam20, pos2, -curve, pixarray[0]);

      /* First order */
      if ((nterm>0)  && pixarray[1]) {
	ObitFArrayShiftAdd(pixarray[1], pos1, beam01, pos2, -flux,  pixarray[1]);
 	ObitFArrayShiftAdd(pixarray[1], pos1, beam11, pos2, -si,    pixarray[1]);
	if (nterm>1) 
	  ObitFArrayShiftAdd(pixarray[1], pos1, beam21, pos2, -curve, pixarray[1]);
      } /* end first order */
    
      /* Second order */
      if ((nterm>1)  && pixarray[2]) {
	ObitFArrayShiftAdd(pixarray[2], pos1, beam02, pos2, -flux,  pixarray[2]);
 	ObitFArrayShiftAdd(pixarray[2], pos1, beam12, pos2, -si,    pixarray[2]);
	ObitFArrayShiftAdd(pixarray[2], pos1, beam22, pos2, -curve, pixarray[2]);
      } /* end second order */
    
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
static olong MakeImSubFuncArgs (ObitThread *thread, olong nterm,
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
    (*ThreadArgs)[i]->inData     = g_malloc0(MAX(1,(nterm+1))*sizeof(ObitFArray*));
    (*ThreadArgs)[i]->nterm      = nterm;
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
    g_free(ThreadArgs[i]->inData);
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
 *        with beam to convolve with.
 * \param CCVer CC table number
 * \param iterm Select spectral term, 0=flux, 1=flux*si, 2=flux*curve...
 * \param err   Obit error stack object.
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
  if (CCTable->noParms>4) { /* Spectral */
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
  ObitImageUtilArray2Image ("DbugRestCC.fits", 0, grid, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, grid);  */
  /* END DEBUG */
  
  /* Cleanup */
  forFFT = ObitFFTUnref(forFFT);
  revFFT = ObitFFTUnref(revFFT);
  uvGrid = ObitCArrayUnref(uvGrid);

  return grid;
} /* end ConvlCC */

/** 
 * Cross convolve a set of Clean components with a beam returning an image array
 * 
 * \param in      Image with CC table and defines size of image grid and
 *        with beam to convolve with.
 * \param CCVer   CC table number
 * \param iterm   Select spectral term, 0=flux, 1=flux*si, 2=flux*curve...
 * \param out     Image defining output grid
 * \param outGrid Grid onto which to accumulate
 * \param err    Obit error stack object.
 */
static void XConvlCC(ObitImage *in, olong CCVer, olong iterm, 
		     ObitImage *out, ObitFArray *outGrid, ObitErr *err)
{
  ObitImageDesc *imDesc1=NULL, *imDesc2=NULL;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitFArray *list = NULL, *tmpArray = NULL;
  olong ver, ncomp, ndim, naxis[2], i=0, n, len;
  ofloat gparm[3], gauss[3], bmaj, bmin, bpa, sr, cr, cellx, celly;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitDConCleanVisWB:XConvlCC";

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
  list = ObitTableCCUtilCrossListSpec (CCTable, imDesc1, imDesc2, 
				       gparm, &ncomp, iterm, err);
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

  /* Move spectral component to flux position in array  
     weight by flux for spectral terms*/
  n   = list->naxis[1];
  len = list->naxis[0];
  if (iterm==0) { /* flux */
    for (i=0; i<n; i++) {
      list->array[i*len+2] = list->array[i*len+3];
    }
  } else { /* spectral */
      list->array[i*len+2] *= list->array[i*len+3];
  }
  
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

