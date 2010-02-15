/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2010                                          */
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
#include "ObitThread.h"

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


/*---------------Private structures----------------*/
/* Image statistics threaded function argument */
typedef struct {
  /* Input plane pixel data */
  ObitFArray *inData;
  /* Input Window */
  ObitDConCleanWindow *window;
  /* Field number (1-rel) */
  olong      field;
  /* First (1-rel) row in image to process this thread */
  olong      first;
  /* Highest (1-rel) row in image to process this thread  */
  olong      last;
  /* thread number, <0 -> no threading  */
  olong      ithread;
  /* Obit Thread object */
  ObitThread  *thread;
  /* Obit error stack object */
  ObitErr    *err;
  /* Number of values in window [out] */
  olong count;
  /* Maximum value in outer window [out] */
  ofloat tmax;
  /* Maximum value in inner window [out] */
  ofloat tmax2;
  /* Sum of values [out] */
  ofloat sum;
  /* Sum of values^2 [out] */
  ofloat sum2;
} StatsFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Deallocate members. */
void  ObitDConCleanInit (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanClear (gpointer in);

/** Private: Set Beam patch, min. flux, decide on SDI CLEAN. */
void ObitDConCleanDecide (ObitDConClean* in, ObitErr *err);

/** Private: Read Beam patches. */
void ReadBP (ObitDConClean* in, ObitErr *err);

/** Private: Apply Gaussian taper to uv grid. */
static void GaussTaper (ObitCArray* uvGrid,  ObitImageDesc *imDesc,
			ofloat gparm[3]);

/** Private: Set Class function pointers. */
static void ObitDConCleanClassInfoDefFn (gpointer inClass);

/** Private: Threaded Image Statistics */
static gpointer ThreadImageStats (gpointer arg);

/** Private: Make Threaded Image statistics args */
static olong MakeStatsFuncArgs (ObitThread *thread,
				ObitDConCleanWindow *window,
				ObitErr *err, StatsFuncArg ***ThreadArgs);

/** Private: Delete Threaded Image statistics args */
static void KillStatsFuncArgs (olong nargs, StatsFuncArg **ThreadArgs);
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
  out->gain        = ObitMemFree (out->gain);
  out->minFlux     = ObitMemFree (out->minFlux);
  out->factor      = ObitMemFree (out->factor);
  out->maxAbsRes   = ObitMemFree (out->maxAbsRes);
  out->avgRes      = ObitMemFree (out->avgRes);
  out->imgRMS      = ObitMemFree (out->imgRMS);
  out->imgPeakRMS  = ObitMemFree (out->imgPeakRMS);
  out->beamPeakRMS = ObitMemFree (out->beamPeakRMS);
  out->currentFields = ObitMemFree (out->currentFields);
  if (out->BeamPatches) {
    for (i=0; i<in->mosaic->numberImages; i++) 
      out->BeamPatches[i] = ObitFArrayUnref(out->BeamPatches[i]);
    out->BeamPatches = ObitMemFree (out->BeamPatches);
  }

  /* In with the new */
  nfield =  in->mosaic->numberImages;
  out->gain        = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->numCurrentField = in->numCurrentField;
  for (i=0; i<nfield; i++) {
    out->gain[i]        = in->gain[i];
    out->minFlux[i]     = in->minFlux[i];
    out->factor[i]      = in->factor[i];
    out->maxAbsRes[i]   = in->maxAbsRes[i];
    out->avgRes[i]      = in->avgRes[i];
    out->imgRMS[i]      = in->imgRMS[i];
    out->imgPeakRMS[i]  = in->imgPeakRMS[i];
    out->beamPeakRMS[i] = in->beamPeakRMS[i];
  }
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<in->numCurrentField; i++) out->currentFields[i] = in->currentFields[i];
  for (i=in->numCurrentField; i<=nfield; i++) out->currentFields[i] = 0;

  if (in->BeamPatches) {
    out->BeamPatches = ObitMemAlloc0Name(in->mosaic->numberImages*sizeof(ObitFArray*),"Beam patches");
    for (i=0; i<in->mosaic->numberImages; i++) 
      out->BeamPatches[i] = ObitFArrayCopy(in->BeamPatches[i], NULL, err);
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
  out->gain        = ObitMemFree (out->gain);
  out->minFlux     = ObitMemFree (out->minFlux);
  out->factor      = ObitMemFree (out->factor);
  out->maxAbsRes   = ObitMemFree (out->maxAbsRes);
  out->avgRes      = ObitMemFree (out->avgRes);
  out->imgRMS      = ObitMemFree (out->imgRMS);
  out->imgPeakRMS  = ObitMemFree (out->imgPeakRMS);
  out->beamPeakRMS = ObitMemFree (out->beamPeakRMS);
  out->currentFields = ObitMemFree (out->currentFields);
  if (out->BeamPatches) {
    for (i=0; i<in->mosaic->numberImages; i++) 
      out->BeamPatches[i] = ObitFArrayUnref(out->BeamPatches[i]);
    out->BeamPatches = ObitMemFree (out->BeamPatches);
  }

  /* In with the new */
  nfield = in->nfield;
  out->gain        = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->numCurrentField = in->numCurrentField;
  for (i=0; i<nfield; i++) {
    out->gain[i]        = in->gain[i];
    out->minFlux[i]     = in->minFlux[i];
    out->factor[i]      = in->factor[i];
    out->maxAbsRes[i]   = in->maxAbsRes[i];
    out->avgRes[i]      = in->avgRes[i];
    out->imgRMS[i]      = in->imgRMS[i];
    out->imgPeakRMS[i]  = in->imgPeakRMS[i];
    out->beamPeakRMS[i] = in->beamPeakRMS[i];
  }
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<in->numCurrentField; i++) out->currentFields[i] = in->currentFields[i];
  for (i=in->numCurrentField; i<=nfield; i++) out->currentFields[i] = 0;

  if (in->BeamPatches) {
    out->BeamPatches = ObitMemAlloc0Name(in->mosaic->numberImages*sizeof(ObitFArray*),"Beam patches");
    for (i=0; i<in->numCurrentField; i++) 
      out->BeamPatches[i] = ObitFArrayCopy(in->BeamPatches[i], NULL, err);
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
  out->gain        = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->imgRMS[i]      = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
    out->currentFields[i] = 0;
  }
  out->currentFields[nfield] = 0;
  out->numCurrentField  = 1;
  out->currentFields[0] = 1;
  out->currentFields[1] = 0;
 
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
  gboolean done, newWin;
  olong jtemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  const ObitDConCleanClassInfo *inClass;
  const ObitDConCleanPxListClassInfo *pxListClass;
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
  /* Pixel list class pointer */
  pxListClass = (ObitDConCleanPxListClassInfo*)in->Pixels->ClassInfo; 

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);

  /* Reset/Init Pixel list*/
  pxListClass->ObitDConCleanPxListReset (in->Pixels, err);
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
    newWin = inClass->ObitDConCleanPixelStats(in, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Pick components for this major cycle, tells if finished CLEAN */
    done = inClass->ObitDConCleanSelect(in, NULL, err);
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
	window[1] = (olong)(in->mosaic->images[field-1]->myDesc->crpix[0]-
			    in->mosaic->images[field-1]->myDesc->xPxOff);
	window[2] = (olong)(in->mosaic->images[field-1]->myDesc->crpix[1]-
                            in->mosaic->images[field-1]->myDesc->yPxOff);
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
	window[0] = (olong)(in->mosaic->images[field-1]->myDesc->crpix[0]-
			    in->mosaic->images[field-1]->myDesc->xPxOff) - 
	  in->mosaic->OutlierSize/2;
	window[1] = (olong)(in->mosaic->images[field-1]->myDesc->crpix[1]-
                            in->mosaic->images[field-1]->myDesc->yPxOff) - 
	  in->mosaic->OutlierSize/2;
	window[2] = (olong)(in->mosaic->images[field-1]->myDesc->crpix[0]-
                            in->mosaic->images[field-1]->myDesc->xPxOff) +
	  in->mosaic->OutlierSize/2;
	window[3] = (olong)(in->mosaic->images[field-1]->myDesc->crpix[1]-
                            in->mosaic->images[field-1]->myDesc->yPxOff) + 
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
 * \param in        The object to deconvolve
 * \param pixarray  If NonNULL use instead of the flux densities from the image file.
 *                  This is an array of ObitFarrays corresponding to fields in
                    in->currentFields
 * \param err       Obit error stack object.
 * \return TRUE if a new window added
 */
gboolean ObitDConCleanPixelStats(ObitDConClean *in, ObitFArray **pixarray, 
				 ObitErr *err)
{
  gboolean newWin = FALSE;
  olong field, ifld;
  ObitImage *Beam=NULL;
  const ObitDConCleanClassInfo *inClass;
  gchar *routine = "ObitDConCleanPixelStats";

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanIsA(in));

  inClass = (ObitDConCleanClassInfo*)in->ClassInfo; /* class structure */

  /* DEBUG
  in->beamPatchSize = 100;
  in->minFluxLoad = 8.0;
  return; */

  /* Adjust windows if autoWindow */
  if (in->autoWindow) 
    newWin = inClass->ObitDConCleanAutoWindow (in, in->currentFields, pixarray, err);
  else {
    in->autoWinFlux = -1.0e20; 
    newWin = FALSE;
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, newWin);

  /* Get Beam histogram from first field */
  if (in->BeamHist->field != in->currentFields[0])
    ObitDConCleanBmHistUpdate(in->BeamHist, Beam, in->plane, err);
  in->BeamHist->field = in->currentFields[0];
  if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
    
  /* Loop over current fields */
  for (ifld=0; ifld<=in->numCurrentField; ifld++) {
    
    field = in->currentFields[ifld];
    if (field<=0) break;  /* List terminates */
    
    /* Beam image */
    Beam = (ObitImage*)(in->mosaic->images[field-1]->myBeam);
    
    /* Get Pixel histogram */
    ObitDConCleanPxHistUpdate (in->PixelHist, field, in->plane, in->mosaic, 
			       in->window, TRUE, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
    
    /* Decide beamPatchSize, minFluxLoad, SDI Clean */
    ObitDConCleanDecide (in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
  }

  return newWin;
} /* end ObitDConCleanPixelStats */

/**
 * Get Image statistics to help decide which field is next to process
 * The outer window is used to specify valid pixels, 
 * For this version the following are calculated:
 * \li maxAbsRes   Maximum absolute windowed residual value (doBEAM=FALSE)
 * \li avgRes      Average windowed residual value (doBEAM=FALSE)
 * \li imgRMS      Image RMS  (doBeam=FALSE)
 * \li imgPeakRMS  Image Peak/RMS  (doBeam=FALSE)
 * \li beamPeakRMS Beam Peak/RMS  (doBeam = TRUE)
 *
 * \param in    The object to deconvolve
 * \param field Which field? (1-rel) <=0 -> all;
 * \param doBeam If TRUE, do Beam statistics else Image
 * \param err   Obit error stack object.
 */
void ObitDConCleanImageStats(ObitDConClean *in, olong field, gboolean doBeam, 
			     ObitErr *err)
{
  ObitIOCode retCode;
  ObitImage *image=NULL;
  ObitDConCleanWindow *theWindow=NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong nThreads, i, it, nTh, nrow, nrowPerThread, hirow, lorow, count;
  olong  lo, hi;
  ofloat tmax, tmax2, sum, sum2;
  gboolean OK;
  StatsFuncArg **threadArgs;
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

  /* Initialize Threading */
  if (!doBeam) theWindow = in->window; /* Use image window */
  nThreads = MakeStatsFuncArgs (in->thread, theWindow, err, &threadArgs);

  /* Set output to full image, plane at a time */
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
 
  /* Loop over selected fields */
  for (i=lo; i<=hi; i++) {

    /* get image */
    if (doBeam)
      image = (ObitImage*)in->mosaic->images[i]->myBeam;
    else
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

    /* Divide up work */
    nrow = image->myDesc->inaxes[1];
    nrowPerThread = nrow/nThreads;
    nTh = nThreads;
    if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
    lorow = 1;
    hirow = nrowPerThread;
    hirow = MIN (hirow, nrow);
    
    /* Set up thread arguments */
    for (it=0; it<nTh; it++) {
      if (it==(nTh-1)) hirow = nrow;  /* Make sure do all */
      if (threadArgs[it]->inData)  ObitFArrayUnref(threadArgs[it]->inData);
      threadArgs[it]->inData  = ObitFArrayRef(image->image);
      threadArgs[it]->field   = i+1;
      threadArgs[it]->first   = lorow;
      threadArgs[it]->last    = hirow;
      if (nTh>1) threadArgs[it]->ithread = it;
      else threadArgs[it]->ithread = -1;
      /* Update which row */
      lorow += nrowPerThread;
      hirow += nrowPerThread;
      hirow = MIN (hirow, nrow);
    }
    
    /* Do operation */
    OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadImageStats,
			   (gpointer**)threadArgs);

    /* Check for problems */
    if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);

    /* save statistics */
    count = 0;
    tmax  = 0.0;
    tmax2 = 0.0;
    sum   = 0.0;
    sum2  = 0.0;
    for (it=0; it<nTh; it++) {
      tmax   = MAX (tmax,   threadArgs[it]->tmax);
      tmax2  = MAX (tmax2,  threadArgs[it]->tmax2);
      count += threadArgs[it]->count;
      sum   += threadArgs[it]->sum;
      sum2  += threadArgs[it]->sum2;
    }    

    /* Beam? */
    if (doBeam) {
      if (count>0) {
	in->beamPeakRMS[i] = tmax / sqrt(sum2/count);
      } else {
	in->beamPeakRMS[i] = 0.0;
      }
    } else { /* Image */
      in->maxAbsRes[i] = MAX (tmax, 0.0);
      if (count>0) {
	in->avgRes[i]     = sum/count;
 	in->imgRMS[i]     = sqrt(sum2/count);
 	in->imgPeakRMS[i] = tmax2 / in->imgRMS[i];
     } else {
	in->avgRes[i]     = 0.0;
 	in->imgRMS[i]     = 0.0;
 	in->imgPeakRMS[i] = 0.0;
     }
   }

    /* Save max residual on image */
    dim[0] = 1;
    ObitInfoListPut (image->info, "maxAbsResid", OBIT_long, dim, &tmax, err); 

    retCode = ObitImageClose (image, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    /* Free Image array? */
    image->image = ObitFArrayUnref(image->image);
    
  } /* end loop over fields */

  /* cleanup */
  KillStatsFuncArgs (nThreads, threadArgs);

} /* end  ObitDConCleanImageStats */

/**
 * Select components to be subtracted
 * Supports multiple fields
 * \param in   The object to deconvolve
 * \param pixarray   If NonNULL use instead of the flux densities from the image file.
 *                   Array of ObitFArrays corresponding to fields in in->currentFields 
 * \param err        Obit error stack object.
 * \return TRUE if deconvolution is complete
 */
gboolean ObitDConCleanSelect(ObitDConClean *in, ObitFArray **pixarray, 
			     ObitErr *err)
{
  gboolean done = FALSE;
  const ObitDConCleanClassInfo *inClass;
  const ObitDConCleanPxListClassInfo *pxListClass;
  gchar *routine = "ObitDConCleanSelect";

  /* error checks */
  if (err->error) return done;
  g_assert (ObitDConCleanIsA(in));

  /* Read beam Patch(es) */
  inClass = (ObitDConCleanClassInfo*)in->ClassInfo; /* class structure */
  inClass->ReadBP (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* Read Pixel List */
  pxListClass = (ObitDConCleanPxListClassInfo*)in->Pixels->ClassInfo; 
  pxListClass->ObitDConCleanPxListUpdate (in->Pixels, in->currentFields, in->numberSkip,
					  in->minFluxLoad, in->autoWinFlux, 
					  in->window, in->BeamPatches, pixarray, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* BGC or SDI Clean? */
  if (in->doSDI) {
    done = pxListClass->ObitDConCleanPxListSDI (in->Pixels, err);
  } else {
    done = pxListClass->ObitDConCleanPxListCLEAN (in->Pixels, err);
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
    for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
    for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
    for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
    ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
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

  if ((in->mosaic->FullField!=NULL) && (in->mosaic->numberImages>1)) {
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
      for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
      for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
      for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
      ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
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
 * Only the field with the highest peak inside the outer window is modified.
 * \param in       The object to restore
 * \param fields   Field numbers (1-rel) in ImageMosaic
 *                 zero terminated list, no longer than in->numCurrentField
 * \param pixarray If nonNULL, then use this array of pixels rather than in the ImageMosaic
 *                 Elements are Pixel arrays corresponding to fields in fields [only 1]
 * \param err      Obit error stack object.
 * \return TRUE is a new window added
 */
gboolean ObitDConCleanAutoWindow(ObitDConClean *in, olong *fields, ObitFArray **pixarray, 
				 ObitErr *err)
{
  gboolean newWin = FALSE;
  ObitImage *image=NULL;
  ObitFArray *usePixels;
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  i, field, best, blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat PeakIn, PeakOut, RMS, bestPeak, *FieldPeak=NULL;
  olong PeakInPos[2] = {0,0};
  gboolean doAbs;
  gchar *routine = "ObitDConCleanAutoWindow";

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
    image = in->mosaic->images[field];
    if (pixarray==NULL) {
      
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
    }
    
    /* Allow negative for Stokes other than I */
    doAbs = fabs (image->myDesc->crval[image->myDesc->jlocs]-1.0) > 0.1;
  
    /* Get statistics */
    ObitDConCleanWindowStats (in->window, field+1, usePixels,
			      doAbs,
			      &PeakIn, &PeakInPos[0], &PeakOut, 
			      &RMS, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, newWin);

    /* Save peak not in window */
    FieldPeak[field] = PeakOut;

    /* Check for best */
    if (PeakIn>bestPeak) {
      best = field;
      bestPeak = PeakIn;
    }
  } /* end loop gathering statistics */

  /* Use best field */
  field = fields[best];

  /* Use passed pixel array or get image */
  image = in->mosaic->images[field-1];
  if (pixarray==NULL) {
    
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
  in->autoWinFlux = MAX (PeakOut-5.0*RMS, 0.5*RMS);
  /* DEBUG  
  fprintf (stderr,"DEBUG autoWinFlux %f RMS %f PeakOut %f PeakIn %f\n",
	   in->autoWinFlux, RMS, PeakOut, PeakIn); */

  /* Use minFlux on Pixels */
  in->minFlux[in->currentFields[0]-1] = 
    MIN (in->minFlux[in->currentFields[0]-1], in->Pixels->minFlux[in->currentFields[0]-1]);

  /* Don't clean too far into the noise */
  if ((in->Pixels->currentIter>1) && (in->Pixels->maxResid<0.5*RMS)) {
     in->minFlux[in->currentFields[0]-1] = 0.5*RMS;
     /* Tell user */
     if (in->prtLv>0) {
       Obit_log_error(err, OBIT_InfoWarn,"Cleaned into noise %g - resetting minFlux",
		      RMS);
     }
  }
  return newWin;
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

  /* private functions */
  theClass->ReadBP = (ReadBPFP)ReadBP;

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
  in->window      = NULL;
  in->BeamPatches = NULL;
  in->BeamHist    = newObitDConCleanBmHist("BeamHist") ;
  in->PixelHist   = newObitDConCleanPxHist("PixelHist");
  in->gain        = NULL;
  in->minFlux     = NULL;
  in->factor      = NULL;
  in->nfield      = 0;
  in->maxAbsRes   = NULL;
  in->imgRMS      = NULL;
  in->imgPeakRMS  = NULL;
  in->beamPeakRMS = NULL;
  in->CCver       = 0;
  in->bmaj        = 0.0;
  in->bmin        = 0.0;
  in->bpa         = 0.0;
  in->niter       = 0;
  in->maxPixel     = 20000;
  in->minPatchSize = 100;
  in->autoWinFlux  = -1.0e20;
  in->ccfLim       = 0.0;
  in->SDIGain      = 0.0;
  in->doSDI        = FALSE;
  in->autoWindow   = FALSE;
  in->currentFields   = NULL;
  in->numCurrentField = 0;

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
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->BeamPatches) {
    for (i=0; i<in->mosaic->numberImages; i++) 
      if (in->BeamPatches[i]) ObitFArrayUnref(in->BeamPatches[i]);
    in->BeamPatches =  ObitMemFree (in->BeamPatches);
  }
  in->window       = ObitDConCleanWindowUnref(in->window);
  in->BeamHist     = ObitDConCleanBmHistUnref(in->BeamHist);
  in->PixelHist    = ObitDConCleanPxHistUnref(in->PixelHist);
  in->Pixels       = ObitDConCleanPxListUnref(in->Pixels);
  in->currentFields= ObitMemFree (in->currentFields);
  in->gain         = ObitMemFree (in->gain);
  in->minFlux      = ObitMemFree (in->minFlux);
  in->factor       = ObitMemFree (in->factor);
  in->maxAbsRes    = ObitMemFree (in->maxAbsRes);
  in->avgRes       = ObitMemFree (in->avgRes);
  in->imgRMS       = ObitMemFree (in->imgRMS);
  in->imgPeakRMS   = ObitMemFree (in->imgPeakRMS);
  in->beamPeakRMS  = ObitMemFree (in->beamPeakRMS);

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
 * currentFields, window, BeamHist and PixelHist being up to date.
 * Member doSDI is set to TRUE iff SDI CLEAN is needed
 * If the histogram is too coarse (too many in top bin) then numberSkip
 * is set to decimate the residuals so that they will fit;
 * \param in   The object to deconvolve
 * \param err  Obit error stack object.
 */
void ObitDConCleanDecide (ObitDConClean* in, ObitErr *err)
{
  olong minPatch, maxPatch, Patch, i;
  ofloat minFlux=0.0, fract;
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

} /* end  ObitDConCleanDecide */

/**
 * Read Beam patches into BeamPatches
 * Loads beam from image fields currentFields and 
 * Beam patch halfwidth of  beamPatchSize.
 * The beam patch is symmetric about the center position allowing the 
 * beam itself not to be symmetric.
 * \param in   The object to deconvolve
 * \param err  Obit error stack object.
 */
void ReadBP (ObitDConClean* in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong  ablc[2], atrc[2], pos[2];
  olong i, j, field, icenx, iceny, nx, ny, mxPatch;
  ofloat fmax;
  ObitImage *Beam;
  ObitFArray *FAtemp=NULL;
  gchar *routine = "ObitDConClean:ReadBP";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Clear old BeamPatches array */
  if (in->BeamPatches) {
    for (i=0; i<in->mosaic->numberImages; i++) 
      if (in->BeamPatches[i]) ObitFArrayUnref(in->BeamPatches[i]);
    in->BeamPatches =  ObitMemFree (in->BeamPatches);
  }
  
  /* Create new BeamPatches array */
  in->BeamPatches = ObitMemAlloc0Name(in->mosaic->numberImages*sizeof(ObitFArray*),"Beam patches");
  for (i=0; i<in->mosaic->numberImages; i++) in->BeamPatches[i] = NULL;

  /* Loop over beams */
  for (i=0; i<in->numCurrentField; i++) {
    field = in->currentFields[i];
    if (field<=0) break;  /* End of list? */
    
    /* Beam image */
    Beam = (ObitImage*)(in->mosaic->images[field-1]->myBeam);
    
    /* Set output to full image, plane at a time */
    dim[0] = IM_MAXDIM;
    for (j=0; j<IM_MAXDIM; j++) blc[j] = 1;
    for (j=0; j<IM_MAXDIM; j++) trc[j] = 0;
    ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, blc, err); 
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
    in->BeamPatches[field-1] = ObitFArraySubArr(Beam->image, ablc, atrc, err);
    if (err->error) Obit_traceback_msg (err, routine, Beam->name);
    
    retCode = ObitImageClose (Beam, err);
    if (err->error) Obit_traceback_msg (err, routine, Beam->name);
    
    /* Free Image array? */
    Beam->image = ObitFArrayUnref(Beam->image);
  } /* end loop over current fields */
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
 * Get image statistics for portion of an image in a thread
 * \param arg Pointer to StatsFuncArg argument with elements:
 * \li inData   ObitFArray with input plane pixel data
 * \li window   Clean Window
 * \li field    Field number (1-rel)
 * \li first    First (1-rel) row in image to process this thread
 * \li last     Highest (1-rel) row in image to process this thread
 * \li ithread  thread number, <0 -> no threading
 * \li err      ObitErr Obit error stack object
 * \li thread   thread Object
 * \return NULL
 */
static gpointer ThreadImageStats (gpointer args)
{
  /* Get arguments from structure */
  StatsFuncArg *largs = (StatsFuncArg*)args;
  ObitFArray *inData    = largs->inData;
  ObitDConCleanWindow *window = largs->window;
  olong      field      = largs->field;
  olong      loRow      = largs->first-1;
  olong      hiRow      = largs->last-1;
  ObitErr    *err       = largs->err;
  ObitThread *thread    = largs->thread;
  /* local */
  olong ix, iy, nx, ny, count, pos[2];
  ofloat *data, tmax, tmax2, sum, sum2;
  gboolean *umask=NULL, *mask=NULL, *innerMask=NULL, isUnbox;
  /*gchar *routine = "ThreadImageStats";*/

  /* Get array pointer */
  pos[0] = 0; pos[1] = loRow;
  data = ObitFArrayIndex(inData, pos);

  /* init */
  nx = inData->naxis[0];
  ny = inData->naxis[1];
  count = 0;
  sum   = 0.0;
  sum2  = 0.0;
  tmax  = -1.0e20;
  tmax2 = -1.0e20;

  for (iy = loRow; iy<=hiRow; iy++) { /* loop in y */
    /* Was a window given? */
    if (window) { /* Yes - there is a window */
      /* Use outer window - Get mask for unwindows */
      isUnbox = ObitDConCleanWindowUnrow(window, field, iy+1, &umask, err);
      if (ObitDConCleanWindowOuterRow(window, field, iy+1, &mask, err)) {
	if (isUnbox) {
	  for (ix=0; ix<nx; ix++) {
	    if (mask[ix] && (!umask[ix])) {
	      count++;
	      sum  += data[ix];
	      sum2 += data[ix]*data[ix];
	      tmax = MAX (tmax, fabs(data[ix]));
	    }
	  } /* end loop over columns - now inner window */
	  if (ObitDConCleanWindowInnerRow(window, field, iy+1, &innerMask, err)) {
	    for (ix=0; ix<nx; ix++) {
	      if (innerMask[ix] && (!umask[ix])) {
		tmax2 = MAX (tmax2, fabs(data[ix]));
	      }
	    }
	  } /* End if inner window */
	} else { /* no unboxes */
	  for (ix=0; ix<nx; ix++) {
	    if (mask[ix]) {
	      count++;
	      sum  += data[ix];
	      sum2 += data[ix]*data[ix];
	      tmax = MAX (tmax, fabs(data[ix]));
	    }
	  }
	  /* get inner max */
	  if (ObitDConCleanWindowInnerRow(window, field, iy+1, &innerMask, err)) {
	    for (ix=0; ix<nx; ix++) {
	      if (innerMask[ix]) tmax2 = MAX (tmax2, fabs(data[ix]));
	    }
	  } /* End if inner window */
	} /* end no unboxes */
      }
    } else { /* No - no window - do everything */
      for (ix=0; ix<nx; ix++) {
	count++;
	sum  += data[ix];
	sum2 += data[ix]*data[ix];
	tmax  = MAX (tmax, fabs(data[ix]));
      }
      tmax2 = tmax;
    }
    data += nx;
  } /* end loop over rows */

  /* If didn't get tmax2, set to tmax */
  if (tmax2<0.0) tmax2 = tmax;
  
  /* save statistics */
  largs->count = count;
  largs->sum   = sum;
  largs->sum2  = sum2;
  largs->tmax  = tmax;
  largs->tmax2 = tmax2;
      
  /* Cleanup */
  if (mask)  mask  = ObitMemFree (mask);
  if (innerMask)  innerMask  = ObitMemFree (innerMask);
  if (umask) umask = ObitMemFree (umask);
  
  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadImageStats */

/**
 * Make arguments for Threaded ThreadImageStats
 * \param thread     ObitThread object to be used 
 * \param inData     Input plane pixel data
 * \param window     Input Window
 * \param err        Obit error stack object.
 * \param ThreadArgs[out] Created array of StatsFuncArg, 
 *                   delete with KillStatsFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeStatsFuncArgs (ObitThread *thread, 
				ObitDConCleanWindow *window,
				ObitErr *err, StatsFuncArg ***ThreadArgs)
{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(StatsFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(StatsFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->field      = 0;
    (*ThreadArgs)[i]->inData     = NULL;
    if (window) (*ThreadArgs)[i]->window = ObitDConCleanWindowRef(window);
    else (*ThreadArgs)[i]->window = NULL;
    (*ThreadArgs)[i]->ithread    = i;
    (*ThreadArgs)[i]->thread     = thread;
    (*ThreadArgs)[i]->err        = err;
  }

  return nThreads;
} /*  end MakeStatsImageArgs */

/**
 * Delete arguments for ThreadImageStats
 * \param nargs      number of elements in args.
 * \param ThreadArgs Array of StatsFuncArg
 */
static void KillStatsFuncArgs (olong nargs, StatsFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->inData)  ObitFArrayUnref(ThreadArgs[i]->inData);
      if (ThreadArgs[i]->window)  ObitDConCleanWindowUnref(ThreadArgs[i]->window);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillStatsFuncArgs */
