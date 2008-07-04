/* $Id$     */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitDConClean.h"
#include "ObitMem.h"
#include "ObitFFT.h"
#include "ObitTableCCUtil.h"
#include "ObitImageUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConClean.c
 * ObitDConClean class function definitions.
 * Virtual CLEAN base class.
 * This class is derived from the ObitDCon class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConClean";

/** Function to obtain parent ClassInfo - ObitDCon */
static ObitGetClassFP ObitParentGetClass = ObitDConGetClass;

/**
 * ClassInfo structure ObitDConCleanClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Deallocate members. */
void  ObitDConCleanInit (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanClear (gpointer in);

/** Private: Set Beam patch, min. flux, decide on SDI CLEAN. */
void ObitDConCleanDecide (ObitDConClean* in, ObitErr *err);

/** Private: Read Beam patch. */
static void ReadBP (ObitDConClean* in, ObitErr *err);

/** Private: Apply Gaussian taper to uv grid. */
static void GaussTaper (ObitCArray* uvGrid,  ObitImageDesc *imDesc,
			ofloat gparm[3]);

/** Private: Set Class function pointers. */
static void ObitDConCleanClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Virtual routine - should never be called
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConClean* newObitDConClean (gchar* name)
{
  ObitDConClean* out;
  gchar *routine = "newObitDConClean";

  /* Virtual */
  g_error("%s: Virtual routine - should not be called",routine);

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConClean));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanInit((gpointer)out);

 return out;
} /* end newObitDConClean */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanGetClass */

/**
 * Make a deep copy of an ObitDConClean.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConClean* ObitDConCleanCopy  (ObitDConClean *in, ObitDConClean *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i, nfield;
  gchar *outName;
  gchar *routine = "ObitDConCleanCopy";

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
    out = newObitDConClean(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /*  copy this class */
  out->window = ObitDConCleanWindowUnref(out->window);
  out->window = ObitDConCleanWindowCopy(in->window,out->window, err);
  out->CCver   = in->CCver;
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

  /* Arrays */
  /* out with the old */
  out->gain    = ObitMemFree (out->gain);
  out->minFlux = ObitMemFree (out->minFlux);
  out->factor  = ObitMemFree (out->factor);
  out->maxAbsRes  = ObitMemFree (out->maxAbsRes);
  out->avgRes  = ObitMemFree (out->avgRes);

  /* In with the new */
  nfield =  in->mosaic->numberImages;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->maxAbsRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  for (i=0; i<nfield; i++) {
    out->gain[i]    = in->gain[i];
    out->minFlux[i] = in->minFlux[i];
    out->factor[i]  = in->factor[i];
    out->maxAbsRes[i]  = in->maxAbsRes[i];
    out->avgRes[i]  = in->avgRes[i];
  }

  return out;
} /* end ObitDConCleanCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConClean similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanClone  (ObitDConClean *in, ObitDConClean *out, ObitErr *err)
{
  olong i, nfield;
  const ObitClassInfo *ParentClass;
  gchar *routine = "ObitDConCleanClone";

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
  out->CCver  = in->CCver;
  out->window = ObitDConCleanWindowUnref(out->window);
  ObitDConCleanWindowClone(in->window,out->window, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Arrays */
  /* out with the old */
  out->gain    = ObitMemFree (out->gain);
  out->minFlux = ObitMemFree (out->minFlux);
  out->factor  = ObitMemFree (out->factor);
  out->maxAbsRes  = ObitMemFree (out->maxAbsRes);
  out->avgRes  = ObitMemFree (out->avgRes);

  /* In with the new */
  nfield       =  in->nfield;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->maxAbsRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  for (i=0; i<nfield; i++) {
    out->gain[i]    = in->gain[i];
    out->minFlux[i] = in->minFlux[i];
    out->factor[i]  = in->factor[i];
    out->maxAbsRes[i]  = in->maxAbsRes[i];
    out->avgRes[i]  = in->avgRes[i];
  }

} /* end ObitDConCleanClone */

/**
 * Creates an ObitDConClean 
 * VIRTUAL routine - should never be called - 
 * defined for convenience of derived classes 
 * \param name  An optional name for the object.
 * \param mosaic from which to create object
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitDConClean* ObitDConCleanCreate (gchar* name, ObitImageMosaic *mosaic,  
			  ObitErr *err)
{
  olong nfield, i;
  ObitDConClean* out=NULL;
  gchar *routine = "ObitDConCleanCreate";

  /* VIRTUAL */
  Obit_log_error(err, OBIT_Error,"%s: Virtual routine - should not be called",routine);
  return out;

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitImageMosaicIsA(mosaic));

  /* Create basic structure */
  out = newObitDConClean (name);

  /* Save Image Mosaic reference */
  out->mosaic = ObitImageMosaicRef(mosaic);

  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Arrays per field */
  nfield =  mosaic->numberImages;
  out->nfield  = nfield;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->maxAbsRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i] = -1.0;
    out->avgRes[i] = -1.0;
  }
 
  return out;
} /* end ObitDConCleanCreate */

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
 * \li "Beam"   = OBIT_float [3] = (BMAJ, BMIN, BPA) alternate form (",", deg)
 * \li "CCVer"   OBIT_long array    = CLEAN table version for all fields
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimum flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanDeconvolve (ObitDCon *inn, ObitErr *err)
{
  ObitDConClean *in;
  gboolean done;
  olong jtemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  const ObitDConCleanClassInfo *inClass;
  gchar *routine = "ObitDConCleanDeconvolve";

  /* Cast input to this type */
  in = (ObitDConClean*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  inClass = (ObitDConCleanClassInfo*)in->ClassInfo; /* class structure */

  /* Get parameters */
  inClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create Pixel List if needed */
  if (!in->Pixels) {
    in->Pixels = ObitDConCleanPxListCreate("Pixel List", in->mosaic, 
					   in->maxPixel, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);

  /* Reset/Init Pixel list*/
  ObitDConCleanPxListReset (in->Pixels, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Save actual CC version if not specified */
  if (in->CCver<=0) {
    if (in->Pixels->CCver[0]>0) in->CCver = in->Pixels->CCver[0];
    jtemp = in->CCver;
    ObitInfoListAlwaysPut(in->info, "CCVer", OBIT_long, dim, &jtemp);
  }

  /* Loop until Deconvolution done */
  done = FALSE;
  while (!done) {
    /* Get image/beam statistics needed for this cycle, decide CLEAN type */
    inClass->ObitDConCleanPixelStats(in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Pick components for this major cycle, tells if finished CLEAN */
    done = inClass->ObitDConCleanSelect(in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);

    /* Subtract components and make new residual(s) */
    inClass->ObitDConCleanSub(in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

  } /* end clean loop */
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Restore */
  inClass->ObitDConCleanRestore(in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Cross Restore if multiple overlapping fields */
  inClass->ObitDConCleanXRestore(in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

  /* Flatten if needed */
  inClass->ObitDConCleanFlatten(in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (in->prtLv>1) ObitErrLog(err);  /* Progress Report */
  else ObitErrClear(err);

} /* end ObitDConCleanDeconvolve */

/**
 * Read any base class parameters and then
 * read CLEAN control parameters from the ObitInfoList member.
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 100]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"   = OBIT_float [3] = (BMAJ, BMIN, BPA) alternate form
 * \li "CCVer"   OBIT_long array    = CLEAN table version for all fields
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 *                                   If only one given it is used for all.
 * \li "minFlux" OBIT_float array  = Minimum flux density (Jy)  per field
 *                                   If only one given it is used for all.
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "autoWindow" OBIT_boolean scalar = True if autoWindow feature wanted.
 * \li "ccfLim"  OBIT_float        = Min. fraction of residual peak to CLEAN to
 *                                   clipped to [0.0,0.9]
 * \li "SDIGain" OBIT_float        = Fraction of pixels in the upper half of the pixel
 *                                   histogram to trigger SDI mode. 
 * From Parent classes:
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 *                                   def (1,1,1,1,1)
 * \param in  The CLEAN object as base class
 * \param err Obit error stack object.
 */
void  ObitDConCleanGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConClean *in = (ObitDConClean*)inn;  /* as this class */
  ObitDConClassInfo *ParentClass;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat beam[3];
  olong i, itemp;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitDConCleanGetParms";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Niter */
  InfoReal.itg = 0;type = OBIT_oint;
  ObitInfoListGetTest(in->info, "Niter", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  in->niter  = itemp;

  /* Restoring beam */
  ObitInfoListGetTest(in->info, "BMAJ", &type, (gint32*)dim, &in->bmaj);
  ObitInfoListGetTest(in->info, "BMIN", &type, (gint32*)dim, &in->bmin);
  ObitInfoListGetTest(in->info, "BPA",  &type, (gint32*)dim, &in->bpa);
  /* Try alternate form - all in beam */
  beam[0] = in->bmaj; beam[1] = in->bmin; beam[2] = in->bpa;
  ObitInfoListGetTest(in->info, "Beam",  &type, dim, beam);
  in->bmaj = beam[0]/3600.0; in->bmin = beam[1]/3600.0; in->bpa = beam[2];

  /* Loop CC version for all fields */
  InfoReal.itg = in->CCver; type = OBIT_oint;
  ObitInfoListGetTest(in->info, "CCVer", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  in->CCver  = itemp;

  /* Loop GAIN per field */
  ObitInfoListGetTest(in->info, "Gain", &type, (gint32*)dim, in->gain);
  /* If only one, use for all */
  if ((in->nfield>1) && (dim[0]==1))
    for (i=1; i<in->nfield; i++) in->gain[i] = in->gain[0];

  /* Minimum flux density per field */
  ObitInfoListGetTest(in->info, "minFlux", &type, (gint32*)dim, in->minFlux);
  /* If only one, use for all */
  if ((in->nfield>1) && (dim[0]==1))
    for (i=1; i<in->nfield; i++) in->minFlux[i] = in->minFlux[0];

  /* CLEAN depth factor per field */
  ObitInfoListGetTest(in->info, "Factor", &type, (gint32*)dim, in->factor);
  /* If only one, use for all */
  if ((in->nfield>1) && (dim[0]==1))
    for (i=1; i<in->nfield; i++) in->factor[i] = in->factor[0];

  /* Maximum number of residual pixels in CLEAN */
  InfoReal.itg = in->maxPixel; type = OBIT_oint;
  ObitInfoListGetTest(in->info, "maxPixel", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (itemp>0) in->maxPixel = itemp;

  /* Minimum beam patch */
  InfoReal.itg = in->minPatchSize; type = OBIT_oint;
  ObitInfoListGetTest(in->info, "minPatch", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  in->minPatchSize = MAX (itemp, 50); /* Set lower bound */

  /* auto Window */
  ObitInfoListGetTest(in->info, "autoWindow", &type, dim, &in->autoWindow);
  /* Set on window object */
  in->window->autoWindow = in->autoWindow;

  /* Min fractional CLEAN */
  ObitInfoListGetTest(in->info, "ccfLim", &type, dim, &in->ccfLim);
  in->ccfLim = MAX (0.0, MIN (0.9, in->ccfLim)); /* to range */

  /* SDI Trip level */
  ObitInfoListGetTest(in->info, "SDIGain", &type, dim, &in->SDIGain);

} /* end ObitDConCleanGetParms */

/**
 * Set default CLEAN windows in mosaic, both inner and outer.
 * If mosaic member  Radius>0 then make round boxes on Fly's eye field
 * with this radius, else use rectangular box including all but outer 5 pixels
 * On outlier fields, use rectangular box of width OutlierSize.
 * If CLEANBox defined in in->info then its contents are used for field 1.
 * Sets outer windows the same as inner windows except for field 1 when 
 * CLEANBox set.
 * If autoWindow, no default inner windows are set.
 * Assumes all images in mosaic have descriptors defined.
 * Any previously existing Windows will be lost.
 * \param in   The CLEAN object, info may have CLEAN Boxes:
 * \li "CLEANBox"   OBIT_long [4,?]  = Array of Clean boxes for field 1
 *    Any entries with first element=0 are ignored.
 *
 * \param err Obit error stack object.
 */
void ObitDConCleanDefWindow(ObitDConClean *in, ObitErr *err)
{
  olong field, sfield = 1, i, j, window[4];
  ObitDConCleanWindowType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType itype;
  olong  *winArray;
  gchar *routine = "ObitDConCleanDefWindow";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Clear any existing windows */
  in->window = ObitDConCleanWindowUnref(in->window);
  in->window = ObitDConCleanWindowCreate("Clean Window", in->mosaic, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* See if auto Window given */
  ObitInfoListGetTest(in->info, "autoWindow", &itype, dim, &in->autoWindow);
  /* Set on window object */
  in->window->autoWindow = in->autoWindow;

  /* See if CLEANBox given */
  if (ObitInfoListGetP (in->info, "CLEANBox", &itype, dim, (gpointer)&winArray)) {
    if (dim[1]>0) {
      j = 0;
      for (i=0; i<dim[1]; i++) {
	/* Default if all zero */
	if ((winArray[j]!=0) && (winArray[j+1]!=0) && 
	    (winArray[j+2]!=0) && (winArray[j+2]!=0)) {
	  /* Round or rectangular? ignore if first element = 0 */
	  if (winArray[j]>0) { /* rectangular */
	    type = OBIT_DConCleanWindow_rectangle;
	    window[0] = winArray[j];
	    window[1] = winArray[j+1];
	    window[2] = winArray[j+2];
	    window[3] = winArray[j+3];
	    ObitDConCleanWindowAdd (in->window, 1, type, window, err);
	    sfield = 2;  /* Done 1 already */
	  } else if (winArray[j]<0){        /* round */
	    type = OBIT_DConCleanWindow_round;
	    window[0] = winArray[j+1];
	    window[1] = winArray[j+2];
	    window[2] = winArray[j+3];
	    ObitDConCleanWindowAdd (in->window, 1, type, window, err);
	    sfield = 2;  /* Done 1 already */
	  }
	  if (err->error) Obit_traceback_msg (err, routine, in->name);
	} /* end if window fully specified */
	j += 4;   /* loop through window array */
      }
    }
  } /* End of use CLEANBox for field 1 */

  if (in->mosaic->nFlyEye>0) {  /* Have Fly's eye use Radius if possible */
    if (in->mosaic->Radius>0.0) {  /* Use round box of radius Radius */
      type = OBIT_DConCleanWindow_round;
      window[0] = (olong)(in->mosaic->Radius);
      for (field=sfield; field<=in->mosaic->nFlyEye; field++) {
	window[1] = (olong)in->mosaic->images[field-1]->myDesc->crpix[0];
	window[2] = (olong)in->mosaic->images[field-1]->myDesc->crpix[1];
	/* Add inner if CLEANBox not specified */
	if ((field>=sfield) && (!in->autoWindow))
	  ObitDConCleanWindowAdd (in->window, field, type, window, err);
	/* Add outer window */
	ObitDConCleanWindowOuter (in->window, field, type, window, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
     }
    } else { /* Use rectangle */
      type = OBIT_DConCleanWindow_rectangle;
      window[0] = 5;
      window[1] = 5;
      for (field=sfield; field<=in->mosaic->nFlyEye; field++) {
	window[2] = in->mosaic->images[field-1]->myDesc->inaxes[0]-5;
	window[3] = in->mosaic->images[field-1]->myDesc->inaxes[1]-5;
	/* Add inner if CLEANBox not specified */
	if ((field>=sfield)  && (!in->autoWindow))
	  ObitDConCleanWindowAdd (in->window, field, type, window, err);
	/* Add outer window */
	ObitDConCleanWindowOuter (in->window, field, type, window, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
     }
    }  /* End Fly's Eye */
    /* Do outliers */
    type = OBIT_DConCleanWindow_rectangle;
    if (in->mosaic->OutlierSize) {
      for (field=in->mosaic->nFlyEye+1; field<=in->mosaic->numberImages; field++) {
	window[0] = (olong)in->mosaic->images[field-1]->myDesc->crpix[0] - 
	  in->mosaic->OutlierSize/2;
	window[1] = (olong)in->mosaic->images[field-1]->myDesc->crpix[1] - 
	  in->mosaic->OutlierSize/2;
	window[2] = (olong)in->mosaic->images[field-1]->myDesc->crpix[0] +
	  in->mosaic->OutlierSize/2;
	window[3] = (olong)in->mosaic->images[field-1]->myDesc->crpix[1] + 
	  in->mosaic->OutlierSize/2;
	/* Only use if not autoWindow */
	if (!in->autoWindow) ObitDConCleanWindowAdd (in->window, field, type, window, err);
	/* Add outer window */
	ObitDConCleanWindowOuter (in->window, field, type, window, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
     }
    } else {  /* Default for outliers */
      window[0] = 5;
      window[1] = 5;
      for (field=in->mosaic->nFlyEye+1; field<=in->mosaic->numberImages; field++) {
	window[2] = in->mosaic->images[field-1]->myDesc->inaxes[0]-5;
	window[3] = in->mosaic->images[field-1]->myDesc->inaxes[1]-5;
	if (!in->autoWindow) 
	  ObitDConCleanWindowAdd (in->window, field, type, window, err);
	/* Add outer window */
	ObitDConCleanWindowOuter (in->window, field, type, window, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
     }
   }
  } else { /* No Fly's eye - rectangular boxes for all */
    type = OBIT_DConCleanWindow_rectangle;
    window[0] = 5;
    window[1] = 5;
    for (field=sfield; field<=in->mosaic->numberImages; field++) {
      window[2] = in->mosaic->images[field-1]->myDesc->inaxes[0]-5;
      window[3] = in->mosaic->images[field-1]->myDesc->inaxes[1]-5;
      if (!in->autoWindow) 
	ObitDConCleanWindowAdd (in->window, field, type, window, err);
      /* Add outer window */
      ObitDConCleanWindowOuter (in->window, field, type, window, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
  }

} /* end ObitDConCleanDefWindow */

/**
 * Get image and beam statistics and prepare to deconvolve
 * If autoWindow option is selected then windows may be added 
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanPixelStats(ObitDConClean *in, ObitErr *err)
{
  olong field;
  ObitImage *Beam=NULL;
  const ObitDConCleanClassInfo *inClass;
  gchar *routine = "ObitDConCleanPixelStats";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  inClass = (ObitDConCleanClassInfo*)in->ClassInfo; /* class structure */

  /* DEBUG
  in->beamPatchSize = 100;
  in->minFluxLoad = 8.0;
  return; */

  field = in->currentField;
  /* Beam image */
  Beam = (ObitImage*)(in->mosaic->images[in->currentField-1]->myBeam);

  /* Get Beam histogram */
  if (in->BeamHist->field != in->currentField) 
      ObitDConCleanBmHistUpdate(in->BeamHist, Beam, in->plane, err);
  in->BeamHist->field = in->currentField;
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Adjust window if autoWindow */
  if (in->autoWindow) 
    inClass->ObitDConCleanAutoWindow (in, in->currentField, err);
  else 
    in->autoWinFlux = -1.0e20; 
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Get Pixel histogram */
  ObitDConCleanPxHistUpdate (in->PixelHist, field, in->plane, in->mosaic, 
			     in->window, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Decide beamPatchSize, minFluxLoad, SDI Clean */
  ObitDConCleanDecide (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
} /* end ObitDConCleanPixelStats */

/**
 * Get Image statistics to help decide which field is next to process
 * If autoWindow then the outer window is used to specify valid pixels, 
 * else the inner window.
 * For this version the following are calculated:
 * \li maxAbsRes Maximum absolute windowed residual value
 * \li avgRes    Average windowed residual value
 *
 * \param in    The object to deconvolve
 * \param field Which field? (1-rel) <=0 -> all;
 * \param err   Obit error stack object.
 */
void ObitDConCleanImageStats(ObitDConClean *in, olong field, ObitErr *err)
{
  ObitIOCode retCode;
  ObitImage *image=NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong lo, hi, i, ix, iy, nx, ny, count, pos[2];
  ofloat *data, tmax, sum;
  gboolean *umask=NULL, *mask=NULL, isUnbox;
  gchar *routine = "ObitDConCleanImageStats";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  if (field>in->nfield) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return;
  }

  /* Field range (0-rel) */
  if (field>0) {
    lo = field-1;
    hi = field-1;
  } else { /* all */
    lo = 0;
    hi = in->nfield-1;
  }

  /* Set output to full image, plane at a time */
  blc[0] = blc[1] = 1;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  trc[0] = trc[1] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
 
  /* Loop over selected fields */
  for (i=lo; i<=hi; i++) {

    /* get image */
    image = in->mosaic->images[i];

    /* Set output to full image, plane at a time */
    dim[0] = IM_MAXDIM;
    ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
    dim[0] = 1;
    ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);
 
    retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    retCode = ObitImageRead (image, image->image->array, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* Loop over image getting statistics in window */
    nx = image->myDesc->inaxes[0];
    ny = image->myDesc->inaxes[1];
    count = 0;
    sum = 0.0;
    tmax = -1.0e20;
    pos[0] = pos[1] = 0;
    data = ObitFArrayIndex(image->image, pos);
    for (iy=0; iy<ny; iy++) {   /* loop over rows */
      /* Use inner or outer window? */
      if (in->autoWindow) { /* autoWindow mode outer window */
	/* Get mask for unwindows */
	isUnbox = ObitDConCleanWindowUnrow(in->window, i+1, iy+1, &umask, err);
	if (ObitDConCleanWindowOuterRow(in->window, i+1, iy+1, &mask, err)) {
	  if (isUnbox) {
	    for (ix=0; ix<nx; ix++) {
	      if (mask[ix] && (!umask[ix])) {
		count++;
		sum += data[ix];
		tmax = MAX (tmax, fabs(data[ix]));
	      }
	    }
	  } else { /* no unboxes */
	    for (ix=0; ix<nx; ix++) {
	      if (mask[ix]) {
		count++;
		sum += data[ix];
		tmax = MAX (tmax, fabs(data[ix]));
	      }
	    }
	  }
	}
      } else { /* use inner window */
	/* Get window mask */
	if (ObitDConCleanWindowRow(in->window, i+1, iy+1, &mask, err)) {
	  for (ix=0; ix<nx; ix++) {
	    if (mask[ix]) {
	      count++;
	      sum += data[ix];
	      tmax = MAX (tmax, fabs(data[ix]));
	    }
	  }
	}
      } /* end branch for inner or outer window */
      data += nx;
    } /* end loop over rows */
    
    /* save statistics */
    in->maxAbsRes[i] = MAX (tmax, 0.0);
    if (count>0) in->avgRes[i] = sum/count;
    else in->avgRes[i] = 0.0;
    
    /* Save max residual on image */
    dim[0] = 1;
    ObitInfoListPut (image->info, "maxAbsResid", OBIT_long, dim, &tmax, err); 

    retCode = ObitImageClose (image, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    /* Free Image array? */
    image->image = ObitFArrayUnref(image->image);
    
    /* Cleanup */
    if (mask)  mask  = ObitMemFree (mask);
    if (umask) umask = ObitMemFree (umask);
  } /* end loop over fields */

} /* end  ObitDConCleanImageStats */

/**
 * Select components to be subtracted
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 * \return TRUE if deconvolution is complete
 */
gboolean ObitDConCleanSelect(ObitDConClean *in, ObitErr *err)
{
  gboolean done = FALSE;
  olong fields[20];
  gchar *routine = "ObitDConCleanSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanIsA(in));

  /* WHAT about multiple fields at a time?*/

  /* Read beam Patch */
  ReadBP (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* Read Pixel List */
  fields[0] = in->currentField; fields[1] = 0;
  ObitDConCleanPxListUpdate (in->Pixels, fields, in->numberSkip,
			     in->minFluxLoad, in->autoWinFlux, 
			     in->window, in->BeamPatch, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* BGC or SDI Clean? */
  if (in->doSDI) {
    done = ObitDConCleanPxListSDI (in->Pixels, err);
  } else {
    done = ObitDConCleanPxListCLEAN (in->Pixels, err);
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  return done;
} /* end ObitDConCleanSelect */

/**
 * Subtract components and generate new residual image(s).
 * Virtual routine, only defined in derived classes
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanSub(ObitDConClean *in, ObitErr *err)
{
  gchar *routine = "ObitDConCleanSub";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Virtual routine */
  g_error ("%s: Virtual routine", routine);

} /* end ObitDConCleanSub */

/**
 * Restore components removed from the residual image(s)
 * \param in   The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanRestore(ObitDConClean *in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitImage *image=NULL;
  ObitImageDesc *imDesc = NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM], ddim[2];
  olong i, field, first, last, ncomp, ver, ndim, naxis[2];
  gchar *tabType = "AIPS CC";
  ofloat gparm[3], bmaj, bmin, bpa;
  ObitFArray *grid = NULL;
  ObitCArray *uvGrid = NULL;
  ObitFFT *forFFT = NULL, *revFFT = NULL;
  gchar *routine = "ObitDConCleanRestore";

   /* DEBUG 
     ObitFArray *tempFArray = NULL;*/
  /* END DEBUG */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Tell user */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Restoring components");
    ObitErrLog(err);  /* Progress Report */
  }

  /* Loop over fields */
  for (field = 0; field<in->nfield; field++) {

    /* which Image? */
    image = in->mosaic->images[field];
    imDesc = image->myDesc;

    /* Full field, correct plane */
    dim[0] = IM_MAXDIM;
    blc[0] = blc[1] = 1;
    for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
    ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
    trc[0] = trc[1] = 0;
    for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
    ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
    dim[0] = 1;
    ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);

    /* Open Image */
    retCode = ObitImageOpen (image, OBIT_IO_ReadWrite, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);

    /* Restoring beam, use value on image if given (else bmaj==0)
       or the value in the image header */
    if (in->bmaj>0.0) { /* value on CLEAN */
      bmaj = in->bmaj;
      bmin = in->bmin;
      bpa  = in->bpa;
      image->myDesc->beamMaj = bmaj;
      image->myDesc->beamMin = bmin;
      image->myDesc->beamPA  = bpa;
    } else { /* use header value */
      bmaj = image->myDesc->beamMaj;
      bmin = image->myDesc->beamMin;
      bpa  = image->myDesc->beamPA;
    }
    
    /* Get CC table */
    ver = in->CCver;
    tempTable = newObitImageTable (image, OBIT_IO_ReadWrite, tabType, &ver, err);
    if ((tempTable==NULL) || (err->error)) Obit_traceback_msg (err, routine, in->name);
    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    in->CCver = ver;  /* save if defaulted (0) */

    /* Grid components */
    first = 0;  /* all */
    last = 0;
    retCode = ObitTableCCUtilGrid (CCTable, 1, &first, &last, FALSE, 1.0, 0.0, 1.0e20,
				   imDesc, &grid, gparm, &ncomp, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);

    /* DEBUG
    fprintf (stderr,"%s: ncomp %d  %d  %d\n",routine,ncomp,first,last); */

    /* Free CC table */
    CCTable = ObitTableCCUnref(CCTable);

    /* DEBUG */
    /* ObitImageUtilArray2Image ("DbugGridCC.fits", 1, grid, err);  */
    /* if (err->error) Obit_traceback_msg (err, routine, in->name);*/
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

    /* read residuals */
    retCode = ObitImageRead (image, image->image->array, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* DEBUG */
    /* ObitImageUtilArray2Image ("DbugRestCC.fits", 1, grid, err); */
    /* if (err->error) Obit_traceback_msg (err, routine, in->name);  */
    /* END DEBUG */

    /* Add restored components */
    ObitFArrayAdd (image->image, grid, grid);
    
    /* Close Image and reopen to reposition */
    retCode = ObitImageClose (image, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
        retCode = ObitImageOpen (image, OBIT_IO_ReadWrite, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);

    /* rewrite image */
    retCode = ObitImageWrite (image, grid->array, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* Close Image */
    retCode = ObitImageClose (image, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);

    /* Free image memory */
    image->image = ObitFArrayUnref(image->image);
  } /* end loop over fields */

    /* Cleanup */
  forFFT = ObitFFTUnref(forFFT);
  revFFT = ObitFFTUnref(revFFT);
  grid   = ObitFArrayUnref(grid);
  uvGrid = ObitCArrayUnref(uvGrid);

} /* end ObitDConCleanRestore */

/**
 * Flatten multiple facets if needed
 * Does Flatten if FullField member of mosaic member is defined.
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanFlatten(ObitDConClean *in, ObitErr *err)
{
  gchar *routine = "ObitDConCleanFlatten";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  if (in->mosaic->FullField!=NULL) {
    /* Tell user */
    if (in->prtLv>1) {
      Obit_log_error(err, OBIT_InfoErr,"Flattening images");
      ObitErrLog(err);  /* Progress Report */
    }
   
    ObitImageMosaicFlatten (in->mosaic, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitDConCleanFlatten */

/**
 * Restore components removed from one field but also 
 * appearing in another.  Does brute force convolution.
 * Adopted from the AIPSish QOOP:QCLEAN.FOR(CLOVER)
 * Presumes in->mosaic and image descriptors filled in.
 * \param in   The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanXRestore(ObitDConClean *in, ObitErr *err)
{
  olong i, ifield, jfield, ncomp, ver, ndim, naxis[2];
  ObitImage *image=NULL;
  ObitImageDesc *imDesc1=NULL, *imDesc2=NULL;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitFArray *list = NULL, *tmpArray = NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat gparm[3], gauss[3], bmaj, bmin, bpa, sr, cr, cellx, celly;
  gchar *tabType = "AIPS CC";
  gboolean gotSome;
  gchar *routine = "ObitDConCleanXRestore";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Need multiple fields for this to be of use */
  if (in->nfield<=1) return;

  /* Tell user */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Cross restoring components");
    ObitErrLog(err);  /* Progress Report */
  }

  /* Double Loop over fields */
  for (ifield = 0; ifield<in->nfield; ifield++) {
    imDesc2 = (in->mosaic->images[ifield])->myDesc;

    gotSome = FALSE; /* until proven otherwise */
    for (jfield = 0; jfield<in->nfield; jfield++) {
      if (ifield==jfield) continue; /* don't do same field */
      imDesc1 = (in->mosaic->images[jfield])->myDesc;

      /* Any overlap? */
      if (ObitImageDescOverlap(imDesc1, imDesc2, err)) {

	/* Get CC table */
	ver = in->CCver;
	tempTable = newObitImageTable (in->mosaic->images[jfield], 
				       OBIT_IO_ReadWrite, tabType, &ver, err);
	CCTable = ObitTableCCConvert(tempTable);
	tempTable = ObitTableUnref(tempTable);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
	in->CCver = ver;  /* save if defaulted (0) */

	/* Get list from jfield */
	list = ObitTableCCUtilCrossList (CCTable, imDesc1, imDesc2, 
					 gparm, &ncomp, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);

	/* Free CC table */
	CCTable = ObitTableCCUnref(CCTable);

	if ((in->prtLv>2) && (ncomp>0)) {
	  Obit_log_error(err, OBIT_InfoErr,"Restore %d components from %d to  %d",
			 ncomp, jfield+1, ifield+1);
	  ObitErrLog(err);  /* Progress Report */
	}

	/* Anything to do? */
	if (ncomp>0) {
	  gotSome = TRUE;
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

	} /* end of Anything to do? */

	/* Free list */
	list = ObitFArrayUnref(list);

      } /* end if overlap */
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }/* end inner loop over fields */
    
    /* Do we have something to add */
    if (gotSome) {

      /* Get image to update */
      image = in->mosaic->images[ifield];
      
      /* Full field, correct plane */
      dim[0] = IM_MAXDIM;
      blc[0] = blc[1] = 1;
      for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
      ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
      trc[0] = trc[1] = 0;
      for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
      ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
      dim[0] = 1;
      ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);
      
      /* Open Image */
      retCode = ObitImageOpen (image, OBIT_IO_ReadWrite, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* read image */
      retCode = ObitImageRead (image, image->image->array, err);
      if (err->error) Obit_traceback_msg (err, routine, image->name);
      
      /* Add restored components */
      ObitFArrayAdd (image->image, tmpArray, tmpArray);
      
      /* Close Image and reopen to reposition */
      retCode = ObitImageClose (image, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      retCode = ObitImageOpen (image, OBIT_IO_ReadWrite, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* rewrite image */
      retCode = ObitImageWrite (image, tmpArray->array, err);
      if (err->error) Obit_traceback_msg (err, routine, image->name);
      
      /* Close Image */
      retCode = ObitImageClose (image, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* Free image memory */
      image->image = ObitFArrayUnref(image->image);
      tmpArray = ObitFArrayUnref(tmpArray);  /* free accumulator */
      
    } /* end add to existing image */
  } /* end outer loop over fields */
} /* end ObitDConCleanXRestore */

/**
 * Automatically set a box in window if needed
 * autoWindow feature will automatically set CLEAN windows inside 
 * a predefined outer window.  Each cycle the residuals inside the outer 
 * window are searched to the maximum value; if the peak is outside the 
 * inner window and > 5 sigma, a new round box  is added 
 * to the window.  Cleaning in each cycle will stop when the peak residual 
 * drops to the level of the highest value outside the CLEAN window.
 * \param in    The object to restore
 * \param field Field number (1-rel) in ImageMosaic
 * \param err   Obit error stack object.
 */
void ObitDConCleanAutoWindow(ObitDConClean *in, olong field, ObitErr *err)
{
  ObitImage *image=NULL;
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  i, blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat PeakIn, PeakOut, RMS;
  olong PeakInPos[2] = {0,0};
  gboolean doAbs;
  gchar *routine = "ObitDConCleanAutoWindow";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Set output to full image, plane at a time */
  blc[0] = blc[1] = 1;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  trc[0] = trc[1] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];

  /* get image */
  image = in->mosaic->images[field-1];

  /* Set input to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);
  
  retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  retCode = ObitImageRead (image, image->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  retCode = ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Allow negative for Stokes other than I */
  doAbs = fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1;
  
  /* Get field info - set new box if needed */
  ObitDConCleanWindowAutoWindow (in->window, field, image->image,
				 doAbs,
				 &PeakIn, &PeakInPos[0], &PeakOut, 
				 &RMS, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Free Image array? */
  image->image = ObitFArrayUnref(image->image);
  
  /* Set flux limit for next cycle */
  in->autoWinFlux = MAX (PeakOut-5.0*RMS, 0.5*RMS);
  /* DEBUG  
  fprintf (stderr,"DEBUG autoWinFlux %f RMS %f PeakOut %f PeakIn %f\n",
	   in->autoWinFlux, RMS, PeakOut, PeakIn); */

  /* Use minFlux on Pixels */
  in->minFlux[in->currentField-1] = 
    MIN (in->minFlux[in->currentField-1], in->Pixels->minFlux[in->currentField-1]);

  /* Don't clean too far into the noise */
  if ((in->Pixels->currentIter>1) && (in->Pixels->maxResid<0.5*RMS)) {
     in->minFlux[in->currentField-1] = 0.5*RMS;
     /* Tell user */
     if (in->prtLv>0) {
       Obit_log_error(err, OBIT_InfoWarn,"Cleaned into noise %g - resetting minFlux",
		      RMS);
     }
  }

} /* end ObitDConCleanAutoWindow */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanClassInfo *theClass = (ObitDConCleanClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanGetClass;
  theClass->newObit       = (newObitFP)newObitDConClean;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanInit;
  theClass->ObitDConGetParms        = (ObitDConGetParmsFP)ObitDConCleanGetParms;
  theClass->ObitDConDeconvolve      = (ObitDConDeconvolveFP)ObitDConCleanDeconvolve;
  theClass->ObitDConCleanDefWindow  = (ObitDConCleanDefWindowFP)ObitDConCleanDefWindow;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanPixelStats;
  theClass->ObitDConCleanImageStats = (ObitDConCleanImageStatsFP)ObitDConCleanImageStats;
  theClass->ObitDConCleanSelect  = (ObitDConCleanSelectFP)ObitDConCleanSelect;
  theClass->ObitDConCleanSub     = (ObitDConCleanSubFP)ObitDConCleanSub;
  theClass->ObitDConCleanRestore = (ObitDConCleanRestoreFP)ObitDConCleanRestore;
  theClass->ObitDConCleanFlatten = (ObitDConCleanFlattenFP)ObitDConCleanFlatten;
  theClass->ObitDConCleanXRestore= (ObitDConCleanXRestoreFP)ObitDConCleanXRestore;
  theClass->ObitDConCleanAutoWindow = 
    (ObitDConCleanAutoWindowFP)ObitDConCleanAutoWindow;

  /* *************** CHANGE HERE *********************************  */

} /* end ObitDConCleanClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConClean *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->window    = NULL;
  in->BeamPatch = NULL;
  in->BeamHist  = newObitDConCleanBmHist("BeamHist") ;
  in->PixelHist = newObitDConCleanPxHist("PixelHist");
  in->gain      = NULL;
  in->minFlux   = NULL;
  in->factor    = NULL;
  in->nfield    = 0;
  in->maxAbsRes = NULL;
  in->avgRes    = NULL;
  in->CCver     = 0;
  in->bmaj      = 0.0;
  in->bmin      = 0.0;
  in->bpa       = 0.0;
  in->niter     = 0;
  in->maxPixel     = 20000;
  in->minPatchSize = 100;
  in->autoWinFlux  = -1.0e20;
  in->ccfLim       = 0.0;
  in->SDIGain      = 0.0;
  in->doSDI        = FALSE;
  in->autoWindow   = FALSE;

} /* end ObitDConCleanInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConClean* cast to an Obit*.
 */
void ObitDConCleanClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConClean *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->window    = ObitDConCleanWindowUnref(in->window);
  in->BeamPatch = ObitFArrayUnref(in->BeamPatch);
  in->BeamHist  = ObitDConCleanBmHistUnref(in->BeamHist);
  in->PixelHist = ObitDConCleanPxHistUnref(in->PixelHist);
  in->Pixels    = ObitDConCleanPxListUnref(in->Pixels);
  in->gain        = ObitMemFree (in->gain);
  in->minFlux     = ObitMemFree (in->minFlux);
  in->factor      = ObitMemFree (in->factor);
  in->maxAbsRes   = ObitMemFree(in->maxAbsRes);
  in->avgRes      = ObitMemFree(in->avgRes);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanClear */

/**
 * Determine Beam patch size and minimum flux density to consider
 * such that the minimum flux density is the maximum ignored sidelobe
 * of the brightest pixel value.
 * Output members are beamPatchSize, minFluxLoad and depend on the members
 * currentField, window, BeamHist and PixelHist being up to date.
 * Member doSDI is set to TRUE iff SDI CLEAN is needed
 * If the histogram is too coarse (too many in top bin) then numberSkip
 * is set to decimate the residuals so that they will fit;
 * \param in   The object to deconvolve
 * \param err  Obit error stack object.
 */
void ObitDConCleanDecide (ObitDConClean* in, ObitErr *err)
{
  olong minPatch, maxPatch, Patch, i;
  ofloat minFlux, fract;
  gchar *routine = "ObitDConCleanDecide";

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
    MIN (ObitDConCleanWindowSize(in->window, in->currentField, err), maxPatch);
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
		       routine, in->currentField, minFlux, Patch, 
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
      in->minFluxLoad = MAX (in->minFluxLoad, in->minFlux[in->currentField]);
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

} /* end  ObitDConCleanDecide */

/**
 * Read Beam patch into BeamPatch
 * Loads beam from image field currentField and 
 * Beam patch halfwidth of  beamPatchSize.
 * The beam patch is symmetric about the center position allowing the 
 * beam itself not to be symmetric.
 * \param in   The object to deconvolve
 * \param err  Obit error stack object.
 */
static void ReadBP (ObitDConClean* in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong  ablc[2], atrc[2], pos[2];
  olong icenx, iceny, nx, ny, mxPatch;
  ofloat fmax;
  ObitImage *Beam;
  ObitFArray *FAtemp=NULL;
  gchar *routine = "ObitDConClean:ReadBP";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));
  
  /* Beam image */
  Beam = (ObitImage*)(in->mosaic->images[in->currentField-1]->myBeam);
  
  /* Set output to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = blc[2] = blc[3] = blc[4] = blc[5] = 1;
  ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, blc, err); 
  trc[0] = trc[1] = trc[2] = trc[3] = trc[4] = trc[5] = 0;
  ObitInfoListPut (Beam->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (Beam->info, "IOBy", OBIT_long, dim, &IOsize, err);
  
  retCode = ObitImageOpen (Beam, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageRead (Beam, Beam->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Compute region of image to take */
  /* Center pixel - make 0-rel */
  nx = Beam->myDesc->inaxes[0];
  ny = Beam->myDesc->inaxes[1];

  /* Find peak in inner quarter */
  /* Window as 0-rel */
  icenx = nx/2;
  iceny = ny/2;
  ablc[0] = icenx - nx/4;
  atrc[0] = icenx + nx/4;
  ablc[1] = iceny - nx/4;
  atrc[1] = iceny + nx/4;

  /* Save Beam patch */
  FAtemp = ObitFArraySubArr(Beam->image, ablc, atrc, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* center = peak */
  fmax = ObitFArrayMax (FAtemp, pos);
  FAtemp = ObitFArrayUnref(FAtemp);
  icenx = pos[0]+ablc[0];
  iceny = pos[1]+ablc[1];

  /* How big can the patch be? */
  mxPatch = MIN (icenx-5, nx-icenx-4);
  mxPatch = MIN ( mxPatch, iceny-5);
  mxPatch = MIN ( mxPatch, ny-iceny-4);
  /* Beam patch can't be larger than this */
  in->beamPatchSize = MIN (in->beamPatchSize, mxPatch);

  /* Beam patch shouldn't be smaller than this */
  in->beamPatchSize = MAX (in->beamPatchSize, 3);

  /* Window as 0-rel */
  ablc[0] = icenx - in->beamPatchSize;
  atrc[0] = icenx + in->beamPatchSize;
  ablc[1] = iceny - in->beamPatchSize;
  atrc[1] = iceny + in->beamPatchSize;

  /* Save Beam patch */
  in->BeamPatch = ObitFArrayUnref(in->BeamPatch);
  in->BeamPatch = ObitFArraySubArr(Beam->image, ablc, atrc, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageClose (Beam, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Free Image array? */
  Beam->image = ObitFArrayUnref(Beam->image);
} /* end ReadBP */

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
  ofloat taper, norm, *grid;
  olong i, j, nx, ny, ndim, naxis[2];

  /* Image info - descriptor should still be valid */
  nx = imDesc->inaxes[imDesc->jlocr];
  ny = imDesc->inaxes[imDesc->jlocd];
  
  /* Normalization factor */
  norm = ((ofloat)nx) * ((ofloat)ny);
  norm = 1.1331 * ((gparm[0]/fabs(imDesc->cdelt[imDesc->jlocr])) * 
		   (gparm[1]/fabs(imDesc->cdelt[imDesc->jlocd]))) / norm;

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
