/* $Id:  $ */
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

#include "ObitThread.h"
#include "ObitSkyModelVMBeam.h"
#include "ObitTableCCUtil.h"
#include "ObitFFT.h"
#include "ObitUVUtil.h"
#include "ObitImageUtil.h"
#include "ObitPBUtil.h"
#include "ObitMem.h"
#include "ObitPrecess.h"
#include "ObitTableANUtil.h"
#include "ObitSkyGeom.h"
#include "ObitSinCos.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVMBeam.c
 * ObitSkyModelVMBeam class function definitions.
 *
 * This class is derived from the #ObitSkyModelVM class
 *
 * This class represents sky models incorporating beam corrections and 
 * their Fourier transforms.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSkyModelVMBeam";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitSkyModelVMGetClass;

/**
 * ClassInfo structure ObitSkyModelVMBeamClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSkyModelVMBeamClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: FT by DFT, may be overridden in derived class */
void ObitSkyModelVMBeamFTDFT (ObitSkyModelVM *in, olong field, 
				ObitUV *uvdata, ObitErr *err);

/** Private: Initialize newly instantiated object. */
void  ObitSkyModelVMBeamInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSkyModelVMBeamClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSkyModelVMBeamClassInfoDefFn (gpointer inClass);

/** Private: Get Inputs. */
void  ObitSkyModelVMBeamGetInput (ObitSkyModel* inn, ObitErr *err);

/** Private: Threaded FTDFT */
static gpointer ThreadSkyModelVMBeamFTDFT (gpointer arg);

/** Private: Threaded FTDFT with phase correction */
static gpointer ThreadSkyModelVMBeamFTDFTPh (gpointer arg);

/** Private: get model frequency primary beam */
static ofloat getPBBeam(ObitImageDesc *desc, ofloat x, ofloat y, 
			ofloat antSize, ofloat pbmin);

/*---------------Private structures----------------*/
/* FT threaded function argument 
 Note: Derived classes MUST have the following entries at the beginning 
 of the corresponding structure */
typedef struct {
  /* type "vmbeam" in this class */
  gchar type[12];
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitSkyModel *in;
  /* Field number being processed (-1 => all) */
  olong        field;
  /* UV data set to model and subtract from current buffer */
  ObitUV       *uvdata;
  /* First (1-rel) vis in uvdata buffer to process this thread */
  olong        first;
  /* Highest (1-rel) vis in uvdata buffer to process this thread  */
  olong        last;
  /* thread number, <0 -> no threading  */
  olong        ithread;
  /* Obit error stack object */
  ObitErr      *err;
  /* Amp, phase interpolator for I pol Beam image */
  ObitFInterpolate *BeamIInterp, *BeamIPhInterp;
  /* Amp, phase interpolator for V pol Beam image */
  ObitFInterpolate *BeamVInterp, *BeamVPhInterp;
  /* Amp, phase interpolator for Q pol Beam image */
  ObitFInterpolate *BeamQInterp, *BeamQPhInterp;
  /* Amp, phase interpolator for U pol Beam image */
  ObitFInterpolate *BeamUInterp, *BeamUPhInterp;
  /* UV Interpolator for FTGrid */
  ObitCInterpolate *Interp;
  /* Start time (days) of validity of model */
  ofloat begVMModelTime;
  /* End time (days) of validity of model */
  ofloat endVMModelTime;
  /* Thread copy of Components list - not used here */
  ObitFArray *VMComps;
  /** Current uv channel number being processed.  */
  olong channel;
  /** Dimension of Rgain...  */
  olong dimGain;
  /** Array of time/spatially variable R component gain, real, imag */
  ofloat *Rgain, *Rgaini;
  /** Array of time/spatially variable L component gain, real, imag */
  ofloat *Lgain, *Lgaini;
  /** Array of time/spatially variable Q component gain, real, imag */
  ofloat *Qgain, *Qgaini;
  /** Array of time/spatially variable U component gain, real, imag */
  ofloat *Ugain, *Ugaini;
} VMBeamFTFuncArg;
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSkyModelVMBeam* newObitSkyModelVMBeam (gchar* name)
{
  ObitSkyModelVMBeam* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMBeamClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSkyModelVMBeam));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSkyModelVMBeamInit((gpointer)out);

 return out;
} /* end newObitSkyModelVMBeam */

/**
 * Initializes from ObitInfoList.
 * Initializes class if needed on first call.
 * \param out     the new object.to be initialized
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string SkyModel type, "Squint" for this class
 *      \li "xxxThreshold" ofloat Threshold flux density for doing high accuracy DFT model
 * \param err     ObitErr for reporting errors.
 */
void ObitSkyModelVMBeamFromInfo (ObitSkyModel *out, gchar *prefix, ObitInfoList *inList, 
				   ObitErr *err)
{ 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *value=NULL;
  gboolean missing;
  gchar *Type = "Squint";
  gchar *routine = "ObitSkyModelVMBeamFromInfo";
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMBeamClassInit();

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(out, &myClassInfo));

  /* check class type */
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  if ((missing) || (type!=OBIT_string) || (!strncmp(Type,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s Wrong class type %s!=%s", routine, value, Type);
    return;
  }

  /* "xxxThreshold" ofloat Threshold flux density for doing high accuracy DFT model */
  if (prefix) keyword = g_strconcat (prefix, "Threshold", NULL);
  else        keyword = g_strdup("Threshold");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->maxGrid);
  g_free(keyword);

} /* end ObitSkyModelVMBeamFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSkyModelVMBeamGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMBeamClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSkyModelVMBeamGetClass */

/**
 * Make a deep copy of an ObitSkyModelVMBeam.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSkyModelVMBeam* 
ObitSkyModelVMBeamCopy  (ObitSkyModelVMBeam *in, 
			   ObitSkyModelVMBeam *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  /*gchar *routine = "ObitSkyModelVMCopy";*/

  /* Copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL) && 
	    /* Don't call yourself */
	    (ParentClass!=(const ObitClassInfo*)&myClassInfo));
  out = ParentClass->ObitCopy (in, out, err);

  /* This class */
  out->Threshold  = in->Threshold;
  out->maxResid   = in->maxResid;
  out->doCrossPol = in->doCrossPol;
  /* Component arrays */
  out->dimGain    = in->dimGain;
  if (in->dimGain>0) {
    if (out->Rgain) g_free(out->Rgain); out->Rgain = NULL;
    out->Rgain = g_malloc0(in->dimGain*sizeof(ofloat));
    for (i=0; i<in->dimGain; i++) out->Rgain[i] = in->Rgain[i];
    if (out->Rgaini) g_free(out->Rgaini); out->Rgaini = NULL;
    if (in->Rgaini) {
      out->Rgaini = g_malloc0(in->dimGain*sizeof(ofloat));
      for (i=0; i<in->dimGain; i++) out->Rgaini[i] = in->Rgaini[i];
    }

    if (out->Lgain) g_free(out->Lgain); out->Lgain = NULL;
    out->Lgain = g_malloc0(in->dimGain*sizeof(ofloat));
    for (i=0; i<in->dimGain; i++) out->Lgain[i] = in->Lgain[i];
    if (out->Lgaini) g_free(out->Lgaini); out->Lgaini = NULL;
    if (in->Lgaini) {
      out->Lgaini = g_malloc0(in->dimGain*sizeof(ofloat));
      for (i=0; i<in->dimGain; i++) out->Lgaini[i] = in->Lgaini[i];
    }
    if (in->doCrossPol) {
      if (out->Qgain) g_free(out->Qgain); out->Qgain = NULL;
      out->Qgain = g_malloc0(in->dimGain*sizeof(ofloat));
      for (i=0; i<in->dimGain; i++) out->Qgain[i] = in->Qgain[i];
      if (in->Qgaini) {
	if (out->Qgaini) g_free(out->Qgaini); out->Qgaini = NULL;
	out->Qgaini = g_malloc0(in->dimGain*sizeof(ofloat));
	for (i=0; i<in->dimGain; i++) out->Qgaini[i] = in->Qgaini[i];
      }

      if (out->Ugain) g_free(out->Ugain); out->Ugain = NULL;
      out->Ugain = g_malloc0(in->dimGain*sizeof(ofloat));
      for (i=0; i<in->dimGain; i++) out->Ugain[i] = in->Ugain[i];
      if (in->Ugaini) {
	if (out->Ugaini) g_free(out->Ugaini); out->Ugaini = NULL;
	out->Ugaini = g_malloc0(in->dimGain*sizeof(ofloat));
	for (i=0; i<in->dimGain; i++) out->Ugaini[i] = in->Ugaini[i];
      }
    }
  }
  return out;
} /* end ObitSkyModelVMBeamCopy */

/**
 * Creates an ObitSkyModelVMBeam 
 * \param name     An optional name for the object.
 * \param mosaic   ObitImageMosaic giving one or more images/CC tables
 * \param uvData   UV data to be operated on
 * \param IBeam    I Beam normalized image
 * \param VBeam    V Beam fractional image
 * \param QBeam    Q Beam fractional image if nonNULL
 * \param UBeam    U Beam fractional image if nonNULL
 * \param IBeamPh  I Beam phase image if nonNULL
 * \param VBeamPh  V Beam phase image if nonNULL
 * \param QBeamPh  Q Beam phase image if nonNULL
 * \param UBeamPh  U Beam phase image if nonNULL
 * \return the new object.
 */
ObitSkyModelVMBeam* 
ObitSkyModelVMBeamCreate (gchar* name, ObitImageMosaic* mosaic,
			  ObitUV *uvData,
			  ObitImage *IBeam,  ObitImage *VBeam, 
			  ObitImage *QBeam,  ObitImage *UBeam, 
			  ObitImage *IBeamPh,  ObitImage *VBeamPh, 
			  ObitImage *QBeamPh,  ObitImage *UBeamPh, 
			  ObitErr *err)
{
  ObitSkyModelVMBeam* out=NULL;
  olong number, i, nchan, nif;
  gchar *routine = "ObitSkyModelVMBeamCreate";

  /* Error tests */
  if (err->error) return out;  /* Previous error */

  /* Create basic structure */
  out = newObitSkyModelVMBeam (name);

  /* Modify for input mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);
  if ((out->mosaic) && (out->mosaic->numberImages>0)) {
    number = out->mosaic->numberImages;
    out->CCver = ObitMemAlloc0 (sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = 0;
    out->startComp = ObitMemAlloc0 (sizeof(olong)*number);
    out->endComp   = ObitMemAlloc0 (sizeof(olong)*number);
  }

  /* Ensure uvData fully instantiated and OK */
  ObitUVFullInstantiate (uvData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);

  /* Swallow Beam images */
  out->doCrossPol = TRUE;
  out->IBeam = ObitImageInterpCreate("IBeam", IBeam, 2, err);
  out->VBeam = ObitImageInterpCreate("VBeam", VBeam, 2, err);
  if (QBeam)
    out->QBeam = ObitImageInterpCreate("QBeam", QBeam, 2, err);
  else
    out->doCrossPol = FALSE;
  if (UBeam)
    out->UBeam = ObitImageInterpCreate("UBeam", UBeam, 2, err);
  else
    out->doCrossPol = FALSE;
  if (err->error) Obit_traceback_val (err, routine, name, out);
  out->numPlane = out->IBeam->nplanes;

  /* Phase beams */
  if (IBeamPh)
    out->IBeamPh = ObitImageInterpCreate("IBeamPh", IBeamPh, 2, err);
  if (VBeamPh)
    out->VBeamPh = ObitImageInterpCreate("VBeamPh", VBeamPh, 2, err);
  if (QBeamPh)
    out->QBeamPh= ObitImageInterpCreate("QBeamPh", QBeamPh, 2, err);
  if (UBeamPh)
    out->UBeamPh = ObitImageInterpCreate("UBeamPh", UBeamPh, 2, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Make sure they are all consistent */
  Obit_retval_if_fail ((ObitFArrayIsCompatable(out->IBeam->ImgPixels, 
					       out->VBeam->ImgPixels)), err, out,
		       "%s: Incompatable I, L beam arrays", routine);
  if (out->doCrossPol) {
    Obit_retval_if_fail ((ObitFArrayIsCompatable(out->IBeam->ImgPixels, 
						 out->QBeam->ImgPixels)), err, out,
			 "%s: Incompatable I, Q beam arrays", routine);
    Obit_retval_if_fail ((ObitFArrayIsCompatable(out->IBeam->ImgPixels, 
						 out->UBeam->ImgPixels)), err, out,
		       "%s: Incompatable I, U beam arrays", routine);
  }
  if (out->IBeamPh) {
    Obit_retval_if_fail ((ObitFArrayIsCompatable(out->IBeam->ImgPixels, 
						 out->IBeamPh->ImgPixels)), err, out,
			 "%s: Incompatable I amp, phase beam arrays", routine);
  }
  if (out->VBeamPh) {
    Obit_retval_if_fail ((ObitFArrayIsCompatable(out->VBeam->ImgPixels, 
						 out->VBeamPh->ImgPixels)), err, out,
			 "%s: Incompatable V amp, phase beam arrays", routine);
  }
  if (out->doCrossPol && out->QBeamPh) {
    Obit_retval_if_fail ((ObitFArrayIsCompatable(out->QBeam->ImgPixels, 
						 out->QBeamPh->ImgPixels)), err, out,
			 "%s: Incompatable Q amp, phase beam arrays", routine);
  }
  if (out->doCrossPol && out->UBeamPh) {
    Obit_retval_if_fail ((ObitFArrayIsCompatable(out->UBeam->ImgPixels, 
						 out->UBeamPh->ImgPixels)), err, out,
			 "%s: Incompatable U amp, phase beam arrays", routine);
  }

 /* Get list of planes per channel */
  nchan = uvData->myDesc->inaxes[uvData->myDesc->jlocf];
  if (uvData->myDesc->jlocf>=0) 
    nif = uvData->myDesc->inaxes[uvData->myDesc->jlocif];
  else nif = 1;
  out->numUVChann = nchan*nif;
  out->FreqPlane  = g_malloc0(out->numUVChann*sizeof(olong));
  for (i=0; i<out->numUVChann; i++) 
    out->FreqPlane[i] = MAX(0, MIN (out->numPlane-1, 
				    ObitImageInterpFindPlane(out->IBeam, uvData->myDesc->freqArr[i])));
  /* Release beam buffers */
  if (IBeam && (IBeam->image)) IBeam->image = ObitImageUnref(IBeam->image);
  if (QBeam && (QBeam->image)) QBeam->image = ObitImageUnref(QBeam->image);
  if (UBeam && (UBeam->image)) UBeam->image = ObitImageUnref(UBeam->image);
  if (VBeam && (VBeam->image)) VBeam->image = ObitImageUnref(VBeam->image);
  if (IBeamPh && (IBeamPh->image)) IBeamPh->image = ObitImageUnref(IBeamPh->image);
  if (QBeamPh && (QBeamPh->image)) QBeamPh->image = ObitImageUnref(QBeamPh->image);
  if (UBeamPh && (UBeamPh->image)) UBeamPh->image = ObitImageUnref(UBeamPh->image);
  if (VBeamPh && (VBeamPh->image)) VBeamPh->image = ObitImageUnref(VBeamPh->image);

  return out;
} /* end ObitSkyModelVMBeamCreate */

/**
 * Initializes Sky Model
 * Checks that data contain RR, LL , save calibration/selection request
 * and set uv data for no selection/calibration
 * \param in      SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamInitMod (ObitSkyModel* inn, ObitUV *uvdata, 
				ObitErr *err)
{
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong iver;
  gboolean btemp;
  olong itemp, numAntList;
  ObitTableList *list=NULL;
  ObitTableAN *TableAN=NULL;
  ObitUVDesc *uvDesc;
  ofloat phase=0.5, cp, sp;
  /*gchar *blank="    ";*/
  olong i;
  VMBeamFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMBeamInitMod";

  if (err->error) return;

  /* Save/reset calibration state */
  /* It's not clear why it was doing this - not always wanted*/
  in->saveDoCalSelect = FALSE;
  ObitInfoListGetTest(uvdata->info, "doCalSelect", &type, dim, &in->saveDoCalSelect);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  btemp = FALSE;
  /*ObitInfoListAlwaysPut (uvdata->info, "doCalSelect", OBIT_bool, dim, &btemp);*/
  in->saveDoCalib = 0;
  ObitInfoListGetTest(uvdata->info, "doCalib", &type, dim, &in->saveDoCalib);
  dim[0] = dim[1] = dim[2] = 1;
  itemp = -1;
  /*ObitInfoListAlwaysPut (uvdata->info, "doCalib", OBIT_long, dim, &itemp);*/
  strncpy (in->saveStokes, "    ", 4);
  ObitInfoListGetTest(uvdata->info, "Stokes", &type, dim, in->saveStokes);
  dim[0] = 4; dim[1] = dim[2] = 1;
  /* ObitInfoListAlwaysPut (uvdata->info, "Stokes", OBIT_string, dim, blank);*/

  /* Open/close to reset */
  ObitUVOpen (uvdata, OBIT_IO_ReadOnly, err);
  ObitUVClose (uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Setup for threading - delete existing threadArgs */
  if (in->threadArgs) {
    for (i=0; i<in->nThreads; i++) {
      g_free(in->threadArgs[i]);
    }
    g_free(in->threadArgs);
    in->threadArgs = NULL;
  } 

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(VMBeamFTFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(VMBeamFTFuncArg)); 
  } /* end initialize */
  
  for (i=0; i<in->nThreads; i++) {
    args = (VMBeamFTFuncArg*)in->threadArgs[i];
    strcpy (args->type, "vmbeam");  /* Enter type as first entry */
    args->in     = inn;
    args->uvdata = uvdata;
    args->ithread = i;
    args->err    = err;
    if (in->IBeam) args->BeamIInterp = ObitImageInterpCloneInterp(in->IBeam,err);
    if (in->VBeam) args->BeamVInterp = ObitImageInterpCloneInterp(in->VBeam,err);
    if (in->QBeam) args->BeamQInterp = ObitImageInterpCloneInterp(in->QBeam,err);
    if (in->UBeam) args->BeamUInterp = ObitImageInterpCloneInterp(in->UBeam,err);
    if (in->IBeamPh) args->BeamIPhInterp = ObitImageInterpCloneInterp(in->IBeamPh,err);
    if (in->VBeamPh) args->BeamVPhInterp = ObitImageInterpCloneInterp(in->VBeamPh,err);
    if (in->QBeamPh) args->BeamQPhInterp = ObitImageInterpCloneInterp(in->QBeamPh,err);
    if (in->UBeamPh) args->BeamUPhInterp = ObitImageInterpCloneInterp(in->UBeamPh,err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    args->begVMModelTime = -1.0e20;
    args->endVMModelTime = -1.0e20;
    args->VMComps= NULL;
    args->dimGain= 0;
    args->Rgain  = NULL;
    args->Rgaini = NULL;
    args->Lgain  = NULL;
    args->Lgaini = NULL;
    args->Qgain  = NULL;
    args->Qgaini = NULL;
    args->Ugain  = NULL;
    args->Ugaini = NULL;
  }

  /* Call parent initializer */
  ObitSkyModelVMInitMod(inn, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fourier transform routines - DFT only */
  /* Are phases given? */
  if (in->IBeamPh) 
    in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamFTDFTPh;
  else /* No phase */
    in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamFTDFT;

  /* Check requested Stokes
  Obit_return_if_fail((!strncmp(in->stokes,"    ",4)), err,
		      "%s: Unsupported Stokes %s", routine, in->stokes); */

  /* Check that data contains RR, LL */
  uvDesc = uvdata->myDesc;
  Obit_return_if_fail(((uvDesc->crval[uvDesc->jlocs]==-1.0) && 
		       (uvDesc->inaxes[uvDesc->jlocs]>=2)), err,
		      "%s: RR, LL not in UV data", routine);

  /* Antenna Lists */
  /* Get TableList */
  list = uvdata->tableList;
  numAntList = ObitTableListGetHigh (list, "AIPS AN");  /* How many subarrays? */
  if (numAntList!=in->numAntList) { /* Rebuild Antenna Lists if needed */
    for (iver=1; iver<=in->numAntList; iver++) { 
      in->AntList[iver-1] = ObitAntennaListUnref(in->AntList[iver-1]);
    }
    if (in->AntList) g_free(in->AntList); in->AntList = NULL;
    in->AntList = g_malloc0((numAntList)*sizeof(ObitAntennaList*));
  }
  in->numAntList = numAntList;
  for (iver=1; iver<=numAntList; iver++) { 
    TableAN = newObitTableANValue ("AN", (ObitData*)uvdata, &iver, 
				   OBIT_IO_ReadOnly, 0, 0, err);
    in->AntList[iver-1] = ObitAntennaListUnref(in->AntList[iver-1]);
    in->AntList[iver-1] = ObitTableANGetList(TableAN, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    TableAN = ObitTableANUnref(TableAN);
  }
  
  /* Source */
  if (!in->curSource) in->curSource = newObitSource ("Source");
  /* Get mean position */
  ObitUVGetRADec  (uvdata, &in->curSource->RAMean, &in->curSource->DecMean, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Precess to get Apparent position */
  ObitPrecessUVJPrecessApp (uvdata->myDesc, in->curSource);

  /* Init Sine/Cosine calculator - just to be sure about threading */
  ObitSinCosCalc(phase, &sp, &cp);

} /* end ObitSkyModelVMBeamInitMod */

/**
 * Any shutdown operations needed for a model
 * Restore calibration/selection state
 * \param in  SkyModel to shutdown
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamShutDownMod (ObitSkyModel* inn, ObitUV *uvdata,
				    ObitErr *err)
{
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i;
  VMBeamFTFuncArg *args;

  /* Call parent shutdown */
  ObitSkyModelVMShutDownMod(inn, uvdata, err);

  if (in->threadArgs) {
    /* Check type - only handle "vmbeam" */
    if (!strncmp((gchar*)in->threadArgs[0], "vmbeam", 6)) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMBeamFTFuncArg*)in->threadArgs[i];
	args->BeamIInterp = ObitFInterpolateUnref(args->BeamIInterp);
	args->BeamVInterp = ObitFInterpolateUnref(args->BeamVInterp);
	args->BeamQInterp = ObitFInterpolateUnref(args->BeamQInterp);
	args->BeamUInterp = ObitFInterpolateUnref(args->BeamUInterp);
	args->BeamIPhInterp = ObitFInterpolateUnref(args->BeamIPhInterp);
	args->BeamVPhInterp = ObitFInterpolateUnref(args->BeamVPhInterp);
	args->BeamQPhInterp = ObitFInterpolateUnref(args->BeamQPhInterp);
	args->BeamUPhInterp = ObitFInterpolateUnref(args->BeamUPhInterp);
	if (args->Rgain)  g_free(args->Rgain);
	if (args->Lgain)  g_free(args->Lgain);
	if (args->Qgain)  g_free(args->Qgain);
	if (args->Ugain)  g_free(args->Ugain);
	if (args->Rgaini) g_free(args->Rgaini);
	if (args->Lgaini) g_free(args->Lgaini);
	if (args->Qgaini) g_free(args->Qgaini);
	if (args->Ugaini) g_free(args->Ugaini);
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
    } /* end if this a "vmbeam" threadArg */
  }

  /* Restore calibration state */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "doCalSelect", OBIT_bool, dim, &in->saveDoCalSelect);
  ObitInfoListAlwaysPut (uvdata->info, "doCalib", OBIT_long, dim, &in->saveDoCalib);
  dim[0] = 4; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "Stokes", OBIT_string, dim, in->saveStokes);

} /* end ObitSkyModelVMBeamShutDownMod */

/**
 * Initializes an ObitSkyModel for a pass through data in time order.
 * Resets current times, converts field offsets of components to pointing offsets
 * \param in  SkyModel to initialize
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamInitModel (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  ObitImageDesc *inDesc;
  olong i, npos[2], lcomp, ncomp, ifield;
  ofloat *ccData;
  odouble RAPnt, DecPnt, ra, dec;
  VMBeamFTFuncArg *args;

  /*  Reset time of current model */
  in->begVMModelTime = -1.0e20;
  in->endVMModelTime = -1.0e20;
  in->curVMModelTime = -1.0e20;
  in->modelMode = OBIT_SkyModel_DFT;  /* Only can do DFT */

  /* Threading */
  if (in->threadArgs) {
    /* Check type - only handle "vmbeam" */
    if (!strncmp((gchar*)in->threadArgs[0], "vmbeam", 6)) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMBeamFTFuncArg*)in->threadArgs[i];
	args->begVMModelTime = -1.0e20;
	args->endVMModelTime = -1.0e20;
      }
    } /* end if this a "vmbeam" threadArg */
  }

  /* Get pointing position */
  ObitImageDescGetPoint(in->mosaic->images[0]->myDesc, &RAPnt, &DecPnt);

 /* Convert position offsets in in->comps to offsets from pointing */
  npos[0] = 0; npos[1] = 0; 
  ccData = ObitFArrayIndex(in->comps, npos);
  lcomp = in->comps->naxis[0];  /* Length of row in comp table */
  ncomp = in->numComp;          /* number of components */
  for (i=0; i<ncomp; i++) {
    ifield = ccData[i*lcomp+0]+0.5;
    if (ifield<0) continue;
    inDesc = in->mosaic->images[ifield]->myDesc;
    /* Convert field offset to position */
    ObitSkyGeomXYShift (inDesc->crval[inDesc->jlocr], inDesc->crval[inDesc->jlocd],
			ccData[i*lcomp+1], ccData[i*lcomp+2], ObitImageDescRotate(inDesc),
			&ra, &dec);
    /* Get position offset from pointing */
    ObitSkyGeomShiftXY (RAPnt, DecPnt, ObitImageDescRotate(inDesc), ra, dec,
			&ccData[i*lcomp+1], &ccData[i*lcomp+2]);
  } /* End loop over components */
} /* end ObitSkyModelVMBeamInitModel */

/**
 * Update VM model with time and spatial modifications to model
 * Update model in Rgain, Lgain, Qgain, Ugain etc, members for comps member 
 * Sets begVMModelTime to current time
 * Sets endVMModelTime for when parallactic angle differs by 1 degree.
 * \param in      SkyModelVM 
 *                Parameters on info:(?)
 *                  PBmin = min. beam gain correction
 * \param time    current time (d)
 * \param suba    0-rel subarray number
 * \param uvdata  uv data being modeled.
 * \param ithread which thread (0-rel) , <0-> no threads
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamUpdateModel (ObitSkyModelVM *inn, 
				    ofloat time, olong suba,
				    ObitUV *uvdata, olong ithread, 
				    ObitErr *err)
{
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  olong npos[2], lcomp, ncomp, i, ifield, lithread, plane, kincif, kincf, kindex;
  ofloat *Rgain=NULL,  *Lgain=NULL,  *Qgain=NULL,  *Ugain=NULL, *ccData=NULL;
  ofloat *Rgaini=NULL, *Lgaini=NULL, *Qgaini=NULL, *Ugaini=NULL;
  ofloat curPA, tPA, tTime, bTime, fscale, PBCor, xx, yy;
  ofloat xr, xi, tr, ti, tri, tii;
  ofloat Ipol, Vpol, IpolPh=0.0, VpolPh=0.0, QpolPh=0.0, UpolPh=0.0;
  ofloat fblank = ObitMagicF();
  odouble x, y;
  VMBeamFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMBeamUpdateModel";

  if (err->error) return;

  /* Thread to update */
  lithread = MAX (0, ithread);
  args = (VMBeamFTFuncArg*)in->threadArgs[lithread];
  g_assert (!strncmp(args->type,"vmbeam",6));  /* Test arg */

   /* Check subarray */
  Obit_return_if_fail(((suba>=0) && (suba<in->numAntList)), err,
		      "%s: bad subarray number %d", routine, suba+1);

  /* Check gain array dimension */
  Obit_return_if_fail((args->dimGain>=in->numComp), err,
    "%s:More components %d than allowed in gains  %d", 
    routine, in->numComp, args->dimGain);

  /* Need Parallactic angle */
  tTime = time - 5.0*suba;  /* Correct time for subarray offset  from DBCON */
  bTime = tTime;
  curPA = ObitAntennaListParAng (in->AntList[suba], 1, tTime, in->curSource);

  /* Beginning of model validity */
  args->begVMModelTime = tTime + 5.0*suba;

  /* Find end time of validity - need Parallactic angle */
  tTime += 1.0/1440.0;
  tPA = ObitAntennaListParAng (in->AntList[suba], 1, tTime, in->curSource);
  /* Step by a min until the parallactic angle changes by 1 deg */
  while (fabs(tPA-curPA) < 1.0*DG2RAD) {
    tTime += 1.0/1440.0;
    tPA = ObitAntennaListParAng (in->AntList[suba], 1, tTime, in->curSource);
    /* But not forever */
    if (tTime-time>0.25) break;
  }

  /* Time for next update */
  args->endVMModelTime = tTime + 5.0*suba;

  /* Get parallactic angle half way between now and end */
  bTime = time + (tTime-bTime) / 2.0;
  curPA = ObitAntennaListParAng (in->AntList[suba], 1, bTime, in->curSource);
  curPA *= RAD2DG;  /* In deg */
  
  /* Which antennas are EVLA ? */
  if (in->isEVLA==NULL) {
    ObitThreadLock(in->thread);  /* Lock against other threads */
    if (in->isEVLA==NULL)   /* Still NULL? */
      in->isEVLA = g_malloc((100+in->AntList[suba]->number)*sizeof(gboolean));
    for (i=0; i<in->AntList[suba]->number; i++) {
      in->isEVLA[i] = 
	(in->AntList[suba]->ANlist[i]->AntName[0]=='E') &&
	(in->AntList[suba]->ANlist[i]->AntName[1]=='V') &&
	(in->AntList[suba]->ANlist[i]->AntName[2]=='L') &&
	(in->AntList[suba]->ANlist[i]->AntName[3]=='A');
    }
    ObitThreadUnlock(in->thread); 
  }

  npos[0] = 0; npos[1] = 0; 
  ccData = ObitFArrayIndex(in->comps, npos);
  Rgain  = args->Rgain;
  Lgain  = args->Lgain;
  Rgaini = args->Rgaini;
  Lgaini = args->Lgaini;
  if (in->doCrossPol) {
    Qgain  = args->Qgain;
    Ugain  = args->Ugain;
    Qgaini = args->Qgaini;
    Ugaini = args->Ugaini;
  } else {
    Qgain  = NULL;
    Ugain  = NULL;
    Qgaini = NULL;
    Ugaini = NULL;
  }
  lcomp = in->comps->naxis[0];  /* Length of row in comp table */
  ncomp = in->numComp;          /* number of components */

  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF beforefreq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  /* Scale by ratio of frequency to beam image ref. frequency */
  plane = in->FreqPlane[args->channel];
  kindex = (in->startChannelPB+(in->numberChannelPB/2)-1)*kincf +
    (in->startIFPB+(in->numberIFPB/2)-1)*kincif;
  fscale = uvdata->myDesc->freqArr[kindex] / in->IBeam->freqs[plane];
  
  /* Compute antenna gains and put in to Rgain, Lgain, Qgain, Ugain */
  for (i=0; i<ncomp; i++) {
    Lgain[i]  = 1.0;
    Rgain[i]  = 1.0;
    Lgaini[i] = 0.0;
    Rgaini[i] = 0.0;
    if (in->doCrossPol) {
      Qgain[i]  = 1.0;
      Ugain[i]  = 1.0;
      Qgaini[i] = 0.0;
      Ugaini[i] = 0.0;
    }
    /* Where in the beam? */
    xx = -ccData[i*lcomp+1];  /* AZ opposite of RA; offsets in beam images are of source */
    yy =  ccData[i*lcomp+2];
    /* Scale by ratio of frequency to beam image ref. frequency */
    x = (odouble)xx * fscale;
    y = (odouble)yy * fscale;
    ifield = ccData[i*lcomp+0]+0.5;
    if (ifield<0) continue;

    /* Interpolate gains -RR and LL as voltage gains */
    plane = in->FreqPlane[args->channel];
    Ipol = ObitImageInterpValueInt (in->IBeam, args->BeamIInterp, x, y, -curPA, plane, err);
    /* Get primary beam correction for component */
    PBCor = Ipol / getPBBeam(in->mosaic->images[ifield]->myDesc, xx, yy, in->antSize, 0.01);
    if (in->IBeamPh)
      IpolPh = DG2RAD*ObitImageInterpValueInt (in->IBeamPh, args->BeamIPhInterp, x, y, -curPA, plane, err);
    Vpol = ObitImageInterpValueInt (in->VBeam, args->BeamVInterp, x, y, -curPA, plane, err);
    if (in->VBeamPh)
      VpolPh = DG2RAD*ObitImageInterpValueInt (in->VBeamPh, args->BeamVPhInterp, x, y, -curPA, plane, err);
    if (Ipol!=fblank) {
      /* Phases given? */
      if (in->IBeamPh && in->VBeamPh) {
	Lgain[i]  = (1.0*cos(IpolPh) - Vpol*cos(VpolPh));
	Lgaini[i] = (1.0*sin(IpolPh) - Vpol*sin(VpolPh));
	Rgain[i]  = (1.0*cos(IpolPh) + Vpol*cos(VpolPh));
	Rgaini[i] = (1.0*sin(IpolPh) + Vpol*sin(VpolPh));
      } else { /* no phase */
	Lgain[i] = (PBCor - Vpol);
	Lgaini[i] = 0.0;
	Rgain[i] = (PBCor + Vpol);
	Rgaini[i] = 0.0;
      }

     } else {
      Lgain[i] = 0.001;
      Rgain[i] = 0.001;
    }

    if (in->doCrossPol) {
      Qgain[i] = ObitImageInterpValueInt (in->QBeam, args->BeamQInterp, x, y, -curPA, plane, err);
      if (Qgain[i]==fblank) Qgain[i] = 0.0;
      if (in->QBeamPh)
	QpolPh = DG2RAD*ObitImageInterpValueInt (in->QBeamPh, args->BeamQPhInterp, x, y, -curPA, plane, err);
      Ugain[i] = ObitImageInterpValueInt (in->UBeam, args->BeamUInterp, x, y, -curPA, plane, err);
      if (Ugain[i]==fblank) Ugain[i] = 0.0;
      if (in->UBeamPh)
	UpolPh = DG2RAD*ObitImageInterpValueInt (in->UBeamPh, args->BeamUPhInterp, x, y, -curPA, plane, err);
      /* Rotate Q+jU by parallactic angle */
      /* Phases given? */
      if (in->QBeamPh && in->UBeamPh) {
	xr  = cos (2.0*curPA*DG2RAD);
	xi  = sin (2.0*curPA*DG2RAD);
	tr  = Qgain[i]*cos(QpolPh);
	ti  = Ugain[i]*cos(UpolPh);
	tri = Qgain[i]*sin(QpolPh);
	tii = Ugain[i]*sin(UpolPh);
	Qgain[i] = tr*xr - xi*ti;
	Ugain[i] = tr*xi + xr*ti;
	Qgaini[i] = tri*xr - xi*tii;
	Ugaini[i] = tri*xi + xr*tii;
      } else { /* no phase */
	xr = cos (2.0*curPA*DG2RAD);
	xi = sin (2.0*curPA*DG2RAD);
	tr = Qgain[i];
	ti = Ugain[i];
	Qgain[i]  = tr*xr - xi*ti;
	Qgaini[i] = 0.0;
	Ugain[i]  = tr*xi + xr*ti;
	Ugaini[i] = 0.0;
      }
    }
    
  } /* end loop over components */

} /* end ObitSkyModelVMBeamUpdateModel */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string SkyModel type, "Squint" for this class
 *      \li "xxxThreshold" ofloat Threshold flux density for doing high accuracy DFT model
 * \param err     ObitErr for reporting errors.
 */
void ObitSkyModelVMBeamGetInfo (ObitSkyModel *inn, gchar *prefix, 
				  ObitInfoList *outList, ObitErr *err)
{ 
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *Type="Squint";
  gchar *routine = "ObitSkyModelVMBeamGetInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Use Base class */
  ObitSkyModelGetInfo(inn, prefix, outList, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* set Class type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  dim[0] = strlen(Type);
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, Type);

  /* "xxxThreshold" ofloat Threshold flux density for doing high accuracy DFT model */
  if (prefix) keyword = g_strconcat (prefix, "Threshold", NULL);
  else        keyword = g_strdup("Threshold");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->Threshold);
  g_free(keyword);

} /* end ObitSkyModelVMBeamGetInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSkyModelVMBeamClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSkyModelVMBeamClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSkyModelVMBeamClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSkyModelVMBeamClassInfoDefFn (gpointer inClass)
{
  ObitSkyModelVMBeamClassInfo *theClass = (ObitSkyModelVMBeamClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSkyModelVMBeamClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSkyModelVMBeamClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSkyModelVMBeamGetClass;
  theClass->newObit       = (newObitFP)newObitSkyModelVMBeam;
  theClass->ObitCopy      = (ObitCopyFP)ObitSkyModelVMBeamCopy;
  theClass->ObitClear     = (ObitClearFP)ObitSkyModelVMBeamClear;
  theClass->ObitInit      = (ObitInitFP)ObitSkyModelVMBeamInit;
  theClass->ObitSkyModelCreate       = (ObitSkyModelCreateFP)ObitSkyModelVMBeamCreate;
  theClass->ObitSkyModelInitMod      = (ObitSkyModelInitModFP)ObitSkyModelVMBeamInitMod;
  theClass->ObitSkyModelShutDownMod  = (ObitSkyModelShutDownModFP)ObitSkyModelVMBeamShutDownMod;
  theClass->ObitSkyModelInitModel    = (ObitSkyModelInitModelFP)ObitSkyModelVMBeamInitModel;
  theClass->ObitSkyModelFTDFT        = (ObitSkyModelFTDFTFP)ObitSkyModelVMBeamFTDFT;
  theClass->ObitSkyModelGetInput     = (ObitSkyModelGetInputFP)ObitSkyModelVMBeamGetInput;
  theClass->ObitSkyModelChose        = (ObitSkyModelChoseFP)ObitSkyModelVMBeamChose;
  theClass->ObitSkyModelVMUpdateModel=
    (ObitSkyModelVMUpdateModelFP)ObitSkyModelVMBeamUpdateModel;
  theClass->ObitSkyModelGetInfo= (ObitSkyModelGetInfoFP)ObitSkyModelVMBeamGetInfo;

} /* end ObitSkyModelVMBeamClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSkyModelVMBeamInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModelVMBeam *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->isEVLA       = NULL;
  in->Rgain        = NULL;
  in->Lgain        = NULL;
  in->Qgain        = NULL;
  in->Ugain        = NULL;
  in->Rgaini       = NULL;
  in->Lgaini       = NULL;
  in->Qgaini       = NULL;
  in->Ugaini       = NULL;
  in->IBeam        = NULL;
  in->QBeam        = NULL;
  in->UBeam        = NULL;
  in->VBeam        = NULL;
  in->IBeamPh      = NULL;
  in->QBeamPh      = NULL;
  in->UBeamPh      = NULL;
  in->VBeamPh      = NULL;
  in->AntList      = NULL;
  in->curSource    = NULL;
  in->numAntList   = 0;
  in->Threshold    = 0.0;
  in->maxResid     = 0.0;
} /* end ObitSkyModelVMBeamInit */


/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSkyModelVMBeam* cast to an Obit*.
 */
void ObitSkyModelVMBeamClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  VMBeamFTFuncArg *args;
  ObitSkyModelVMBeam *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->isEVLA) g_free(in->isEVLA); in->isEVLA = NULL;
  if (in->Rgain)  g_free(in->Rgain);  in->Rgain  = NULL;
  if (in->Lgain)  g_free(in->Lgain);  in->Lgain  = NULL;
  if (in->Qgain)  g_free(in->Qgain);  in->Qgain  = NULL;
  if (in->Ugain)  g_free(in->Ugain);  in->Ugain  = NULL;
  if (in->Rgaini) g_free(in->Rgaini); in->Rgaini = NULL;
  if (in->Lgaini) g_free(in->Lgaini); in->Lgaini = NULL;
  if (in->Qgaini) g_free(in->Qgaini); in->Qgaini = NULL;
  if (in->Ugaini) g_free(in->Ugaini); in->Ugaini = NULL;
  in->IBeam     = ObitImageInterpUnref(in->IBeam);
  in->QBeam     = ObitImageInterpUnref(in->QBeam);
  in->UBeam     = ObitImageInterpUnref(in->UBeam);
  in->VBeam     = ObitImageInterpUnref(in->VBeam);
  in->IBeamPh   = ObitImageInterpUnref(in->IBeamPh);
  in->QBeamPh   = ObitImageInterpUnref(in->QBeamPh);
  in->UBeamPh   = ObitImageInterpUnref(in->UBeamPh);
  in->VBeamPh   = ObitImageInterpUnref(in->VBeamPh);
  in->curSource = ObitSourceUnref(in->curSource);
  if (in->AntList)  {
    for (i=0; i<in->numAntList; i++) { 
      in->AntList[i] = ObitAntennaListUnref(in->AntList[i]);
    }
    g_free(in->AntList); in->AntList = NULL;
  }
    
  /* Thread stuff */
  if (in->threadArgs) {
    /* Check type - only handle "vmbeam" */
    if (!strncmp((gchar*)in->threadArgs[0], "vmbeam", 6)) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMBeamFTFuncArg*)in->threadArgs[i];
	if (args->Rgain)  g_free(args->Rgain);
	if (args->Lgain)  g_free(args->Lgain);
	if (args->Qgain)  g_free(args->Qgain);
	if (args->Ugain)  g_free(args->Ugain);
	if (args->Rgaini) g_free(args->Rgaini);
	if (args->Lgaini) g_free(args->Lgaini);
	if (args->Qgaini) g_free(args->Qgaini);
	if (args->Ugaini) g_free(args->Ugaini);
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
    } /* end if this a "vmbeam" threadArg */
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSkyModelVMBeamClear */

/**
 * Get input parameters from info member
 * If maxResid value not given or <0 then for it is determined
 * by examining each image in the Image mosaic.  
 * If any image has an info item "maxAbsResid" the the maximum of any
 * of these is used, else the MAX (fabs(minval), fabs(maxval)) is used.
 * \param in Pointer to the ObitSkyModelVMBeam .
 * \param err Obit error stack object.
 */
void  ObitSkyModelVMBeamGetInput (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitImageDesc *desc;
  gboolean lookup, gotit;
  ofloat maxv, maxAbsResid;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitSkyModelVMBeamGetInput";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelVMBeamIsA(in));
  if (!ObitInfoListIsA(in->info)) return;
  InfoReal.itg = 0;type = OBIT_oint;

  /* Call base class version */
  ObitSkyModelVMGetInput (inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Threshold for high accuracy model */
  InfoReal.flt = in->Threshold; type = OBIT_float;
  ObitInfoListGetTest(in->info, "Threshold", &type, (gint32*)dim, &InfoReal);
  in->Threshold = InfoReal.flt;

  /* Current maximum abs residual flux density */
  InfoReal.flt = -1.0; type = OBIT_float;
  lookup = !ObitInfoListGetTest(in->info, "maxResid", &type, (gint32*)dim, &InfoReal);
  if (lookup || (InfoReal.flt<0.0)) {
    /* Need to lookup from attached Mosaic */
    maxv = 0.0;
    for (i=0; i<in->mosaic->numberImages; i++) {
      gotit = ObitInfoListGetTest (in->mosaic->images[i]->info, "maxAbsResid", 
				   &type, dim, &maxAbsResid); 
      if (gotit) maxv = MAX (maxv, maxAbsResid);
    } /* end loop over mosaic */

    /* If still nothing use max/min in image headers */
    if (maxv==0.0) {
      for (i=0; i<in->mosaic->numberImages; i++) {
	desc = in->mosaic->images[i]->myDesc;
	maxv = MAX (maxv, fabs(desc->maxval));
	maxv = MAX (maxv, fabs(desc->minval));
      }/* end loop over mosaic */
    } /* end if still nothing */	  
    in->maxResid = maxv;	
  } else  in->maxResid = InfoReal.flt;

  /* Make sure doPBCor turned off, will always use beam model instead */
  in->doPBCor = FALSE;
  
} /* end ObitSkyModelVMBeamGetInput */

/**
 * Decide which method is the most appropriate to calculate the FT of a model
 * Sets currentMode member function
 * Only DFT supported
 * \param in     Pointer to the ObitSkyModel .
 * \param uvdata UV data set
 */
void  ObitSkyModelVMBeamChose (ObitSkyModel* inn, ObitUV* uvdata) 
{
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  if (in->maxResid <= 0.0) {/* May be mixed */
    if (in->Threshold>0.0) in->currentMode = OBIT_SkyModel_Mixed;
    else in->currentMode = OBIT_SkyModel_DFT;
    in->maxGrid = in->minDFT = in->Threshold;
    return;
  } else if (in->maxResid >= in->Threshold) {
    /* Want accurate model for everything */
    in->currentMode = OBIT_SkyModel_DFT;
    in->maxGrid = 0.0;
    in->minDFT  = 0.0;
    return;
  } else { /* Nothing special - use base selector */
    ObitSkyModelChose (inn, uvdata);
    in->maxGrid = 1.0e20;
    in->minDFT  = 0.0;
    return;
  }

} /* end ObitSkyModelVMBeamChose */


/**
 * Do Fourier transform using a DFT for a buffer of data.
 * If threading has been enabled by a call to ObitThreadAllowThreads 
 * this routine will divide the buffer up amount the number of processors
 * returned by ObitThreadNumProc.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * After the AIPSish QXXPTS, QPTDIV and friends
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * \param in     SkyModel with model components loaded (ObitSkyModelLoad)
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to model and subtract from current buffer
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamFTDFT (ObitSkyModelVM *inn, olong field, ObitUV *uvdata, ObitErr *err)
{
  olong i, mcomp, iComp, pos[2], nvis, lovis, hivis, nvisPerThread, nThreads;
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  VMBeamFTFuncArg *args;
  ofloat *ddata;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelVMBeamFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* Check */
  args = (VMBeamFTFuncArg*)in->threadArgs[0];
  if (strncmp (args->type, "vmbeam", 6)) {
    Obit_log_error(err, OBIT_Error,"%s: Wrong type FuncArg %s", routine,args->type);
    return;
  }

  /* Count number of actual components */
  mcomp = 0;
  pos[0] = 0; pos[1] = 0;
  ddata = ObitFArrayIndex(in->comps, pos);
  for (iComp=0; iComp<in->comps->naxis[1]; iComp++) {
    if (ddata[3]!=0.0) mcomp = iComp+1;
    ddata += in->comps->naxis[0];  /* update pointer */
  } /* end loop over components */
  in->numComp = mcomp;  /* Number of actual components */

  /* Divide up work */
  nvis = uvdata->myDesc->numVisBuff;
  if (nvis<1000) nThreads = 1;
  else nThreads = in->nThreads;
  nvisPerThread = nvis/nThreads;
  lovis = 1;
  hivis = nvisPerThread;
  hivis = MIN (hivis, nvis);

  /* Set up thread arguments */
  for (i=0; i<nThreads; i++) {
    if (i==(nThreads-1)) hivis = nvis;  /* Make sure do all */
    args = (VMBeamFTFuncArg*)in->threadArgs[i];
    args->in     = (ObitSkyModel*)inn;
    args->field  = field;
    args->uvdata = uvdata;
    args->first  = lovis;
    args->last   = hivis;
    if (nThreads>1) args->ithread= i;
    else args->ithread = -1;
    args->err    = err;
    if (args->dimGain!=in->numComp) {
      if (args->Rgain)  g_free(args->Rgain);
      if (args->Lgain)  g_free(args->Lgain);
      if (args->Qgain)  g_free(args->Qgain); args->Qgain = NULL;
      if (args->Ugain)  g_free(args->Ugain); args->Ugain = NULL;
      if (args->Rgaini) g_free(args->Rgaini);
      if (args->Lgaini) g_free(args->Lgaini);
      if (args->Qgaini) g_free(args->Qgaini); args->Qgaini = NULL;
      if (args->Ugaini) g_free(args->Ugaini); args->Ugaini = NULL;
      args->dimGain = in->numComp;
      args->Rgain  = g_malloc0(args->dimGain*sizeof(ofloat));
      args->Rgaini = g_malloc0(args->dimGain*sizeof(ofloat));
      args->Lgain  = g_malloc0(args->dimGain*sizeof(ofloat));
      args->Lgaini = g_malloc0(args->dimGain*sizeof(ofloat));
      if (in->doCrossPol) {
	args->Qgain  = g_malloc0(args->dimGain*sizeof(ofloat));
	args->Qgaini = g_malloc0(args->dimGain*sizeof(ofloat));
	args->Ugain  = g_malloc0(args->dimGain*sizeof(ofloat));
	args->Ugaini = g_malloc0(args->dimGain*sizeof(ofloat));
      }
    }
    /* Update which vis */
    lovis += nvisPerThread;
    hivis += nvisPerThread;
    hivis = MIN (hivis, nvis);
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nThreads, in->DFTFunc, in->threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
}  /* end ObitSkyModelVMBeamFTDFT */

/**
 * Do Fourier transform using a DFT for a buffer of data.
 * Version for time/spatial dependent effects.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * If doFlip member is true the Fourier transform is multiplied by sqrt(-1)
 * (for Stokes RL and LR)
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Method assumes same correction to all antennas.
 * After the AIPSish QXXPTS, QPTDIV and friends
 * Arguments are given in the VMBeamFTFuncArg structure passed as arg starting 
 * with the following:
 * \li type   String identifying structure
 * \li in     SkyModelVM with model components loaded (ObitSkyModelLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number, <0-> no threads
 * \li err    Obit error stack object.
 * \li begVMModelTime Start time (days) of validity of model
 * \li endVMModelTime End time (days) of validity of model
 * \li VMComps Thread copy of Components list - not used here
 * \li dimGain Dimension of Rgain
 * \li Rgain   Float array of time/spatially variable R component gain
 * \li Lgain   Float array of time/spatially variable L component gain
 * \li Qgain   Float array of time/spatially variable Q component gain
 * \li Ugain   Float array of time/spatially variable U component gain
 * \li channel Current UV channel being processed (used in model update ).
 * \return NULL
 */
static gpointer ThreadSkyModelVMBeamFTDFT (gpointer args)
{
  /* Get arguments from structure */
  VMBeamFTFuncArg *largs = (VMBeamFTFuncArg*)args;
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)largs->in;
  /*olong field      = largs->field;*/
  ObitUV *uvdata   = largs->uvdata;
  olong loVis      = largs->first-1;
  olong hiVis      = largs->last;
  olong ithread    = MAX (0, largs->ithread);
  ObitErr *err     = largs->err;
  /*olong dimGain    = largs->dimGain;*/
  ofloat *Rgain    = largs->Rgain;
  ofloat *Lgain    = largs->Lgain;
  ofloat *Qgain    = largs->Qgain;
  ofloat *Ugain    = largs->Ugain;

  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, ncomp;
  olong lrec, nrparm, naxis[2], channel, plane;
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, ilocb, suba, itemp, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealQ,  sumImagQ,  modRealQ=0.0,  modImagQ=0.0;
  ofloat sumRealU,  sumImagU,  modRealU=0.0,  modImagU=0.0;
  ofloat cbase, *rgain1, *lgain1, *rgain2, *lgain2, tx, ty, tz, ll, lll;
  ofloat *qgain1, *ugain1, *qgain2, *ugain2, re, im, specFact;
  ofloat amp, ampr, ampl, sp=0.0, cp=1.0, arg, freq2=0.0, freqFact, wtRR=0.0, wtLL=0.0, temp;
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArrR[FazArrSize], AmpArrL[FazArrSize], FazArr[FazArrSize];
  ofloat CosArr[FazArrSize], SinArr[FazArrSize];
  olong it, jt, itcnt;
  gboolean doCrossPol;
  odouble *freqArr;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMBeamFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;
  ilocb =  uvdata->myDesc->ilocb;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF       = in->startIFPB-1;
  numberIF      = MAX (1, in->numberIFPB);
  jincif        = uvdata->myDesc->incif;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  nUVchan       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocf];
  jincf         = uvdata->myDesc->incf;
  startPoln     = in->startPoln-1;
  numberPoln    = in->numberPoln;
  nUVpoln       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocs];
  jincs         = uvdata->myDesc->incs;  /* increment in real array */
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF beforefreq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  /* Cross or only parallel pol? */
  doCrossPol = (nUVpoln > 2) && in->doCrossPol;
  /* Only parallel for divide */
  if (in->doDivide) doCrossPol = FALSE;

  /* Get pointer for components */
  naxis[0] = 0; naxis[1] = 0; 
  Data = ObitFArrayIndex(in->comps, naxis);
  lcomp = in->comps->naxis[0];   /* Length of row in comp table */
  ncomp = in->comps->naxis[1];   /* Number of components */
  mcomp = in->numComp;           /* Actual number */

  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

  /* Current channel (0-rel) */
  channel = 0;
  plane   = in->FreqPlane[largs->channel];  /* Which plane in correction cube */

  /* Outer loop over blocks of channels */
  /* Starting parameters this pass */
  lstartIF       = startIF;
  lstartChannel  = startChannel;
  while (channel<(numberIF*numberChannel)) {
    
    /* Loop over vis in buffer */
    lrec    = uvdata->myDesc->lrec;         /* Length of record */
    visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
    nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */
    
    for (iVis=loVis; iVis<hiVis; iVis++) {
      
      /* Is current model still valid? */
      if ((visData[iloct] > largs->endVMModelTime) || (visData[iloct] < largs->begVMModelTime)) {
	/* Current channel - which plane in correction to apply? */
	largs->channel = channel;
	/* Subarray 0-rel */
	itemp = (olong)visData[ilocb];
	suba = 100.0 * (visData[ilocb]-itemp) + 0.5; 
	/* Update */
	myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba, uvdata, ithread, err);
	if (err->error) {
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
			 routine);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}
      }
      
      /* Need antennas numbers */
      cbase = visData[ilocb]; /* Baseline */
      ant1 = (cbase / 256.0) + 0.001;
      ant2 = (cbase - ant1 * 256) + 0.001;
      ant1--;    /* 0 rel */
      ant2--;    /* 0 rel */
      
      /* Loop over IFs */
      channel = lstartIF* nUVchan + lstartChannel; /* UV Channel */
      for (iIF=lstartIF; iIF<startIF+numberIF; iIF++) {
	offsetIF = nrparm + iIF*jincif; 
	for (iChannel=lstartChannel; iChannel<startChannel+numberChannel; iChannel++) {
	  offsetChannel = offsetIF + iChannel*jincf; 
	  freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */
	  
	  /* Sum over components */
	  /* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	  sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	  sumRealQ  = sumImagQ  = sumRealU  = sumImagU  = 0.0;
	  /* Set component gain lists by antenna and type */
	  ddata = Data;
	  rgain1 = Rgain;
	  rgain2 = Rgain;
	  lgain1 = Lgain;
	  lgain2 = Lgain;
	  qgain1 = Qgain;
	  qgain2 = Qgain;
	  ugain1 = Ugain;
	  ugain2 = Ugain;
	  
	  /* Sum by model type - assume phase same for RR, LL */
	  switch (in->modType) {
	  case OBIT_SkyModel_PointMod:     /* Point */
	    /* From the AIPSish QXXPTS.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArrR[itcnt] = ddata[3] * rgain1[iComp];
		AmpArrL[itcnt] = ddata[3] * lgain1[iComp];
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrR[jt]*CosArr[jt];
		sumImagRR += AmpArrR[jt]*SinArr[jt];
		sumRealLL += AmpArrL[jt]*CosArr[jt];
		sumImagLL += AmpArrL[jt]*SinArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrR[jt]+AmpArrL[jt])*CosArr[jt];
		  im = 0.5*(AmpArrR[jt]+AmpArrL[jt])*SinArr[jt];
		  sumRealQ += re * qgain1[jt];
		  sumImagQ += im * qgain1[jt];
		  sumRealU += re * ugain1[jt];
		  sumImagU += im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_PointModSpec:     /* Point + spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[6+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  FazArr[itcnt]  = freqFact * (tx + ty + tz);
		  AmpArrR[itcnt] = specFact * ddata[3] * rgain1[iComp];
		  AmpArrL[itcnt] = specFact * ddata[3] * lgain1[iComp];
		}  /* end if valid */
		ddata += lcomp;  /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);

	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		/* Parallel  pol */
		sumRealRR += AmpArrR[jt]*CosArr[jt];
		sumImagRR += AmpArrR[jt]*SinArr[jt];
		sumRealLL += AmpArrL[jt]*CosArr[jt];
		sumImagLL += AmpArrL[jt]*SinArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrR[jt]+AmpArrL[jt])*CosArr[jt];
		  im = 0.5*(AmpArrR[jt]+AmpArrL[jt])*SinArr[jt];
		  sumRealQ += re * qgain1[jt];
		  sumImagQ += im * qgain1[jt];
		  sumRealU += re * ugain1[jt];
		  sumImagU += im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	    /* From the AIPSish QGASUB.FOR  */
	    freq2 = freqArr[iIF*kincif + iChannel*kincf];    /* Frequency squared */
	    freq2 *= freq2;
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			       ddata[8]*visData[ilocv]*visData[ilocv] +
			       ddata[9]*visData[ilocu]*visData[ilocv]);
		ampr = ddata[3] * rgain1[iComp];
		if (arg>1.0e-5) ampr *= exp (arg);
		sumRealRR += ampr*cp;
		sumImagRR += ampr*sp;
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			       ddata[8]*visData[ilocv]*visData[ilocv] +
			       ddata[9]*visData[ilocu]*visData[ilocv]);
		ampl = ddata[3] * lgain1[iComp];
		if (arg>1.0e-5) ampl *= exp (arg);
		AmpArrR[itcnt] = ampr;
		AmpArrL[itcnt] = ampl;
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }  /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrR[jt]*CosArr[jt];
		sumImagRR += AmpArrR[jt]*SinArr[jt];
		sumRealLL += AmpArrL[jt]*CosArr[jt];
		sumImagLL += AmpArrL[jt]*SinArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrR[it]+AmpArrL[jt])*CosArr[jt];
		  im = 0.5*(AmpArrR[it]+AmpArrL[jt])*SinArr[jt];
		  sumRealQ += re * qgain1[jt];
		  sumImagQ += im * qgain1[jt];
		  sumRealU += re * ugain1[jt];
		  sumImagU += im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	    freq2 = freqFact*freqFact;    /* Frequency factor squared */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[9+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
				 ddata[8]*visData[ilocv]*visData[ilocv] +
				 ddata[9]*visData[ilocu]*visData[ilocv]);
		  if (arg<-1.0e-5) amp = specFact * ddata[3] * exp (arg);
		  else amp = specFact * ddata[3];
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  FazArr[itcnt] = freqFact * (tx + ty + tz);
		  AmpArrR[itcnt] = amp * rgain1[iComp];
		  AmpArrL[itcnt] = amp * lgain1[iComp];
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }  /* end inner loop over components */

	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrR[jt]*CosArr[jt];
		sumImagRR += AmpArrR[jt]*SinArr[jt];
		sumRealLL += AmpArrL[jt]*CosArr[jt];
		sumImagLL += AmpArrL[jt]*SinArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrR[jt]+AmpArrL[jt])*CosArr[jt];
		  im = 0.5*(AmpArrR[jt]+AmpArrL[jt])*SinArr[jt];
		  sumRealQ += re * qgain1[jt];
		  sumImagQ += im * qgain1[jt];
		  sumRealU += re * ugain1[jt];
		  sumImagU += im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	    /* From the AIPSish QSPSUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ddata[7];
		arg = MAX (arg, 0.1);
		arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		AmpArrR[itcnt] = ddata[3] * rgain1[iComp] * arg;
		AmpArrL[itcnt] = ddata[3] * lgain1[iComp] * arg;
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrR[jt]*CosArr[jt];
		sumImagRR += AmpArrR[jt]*SinArr[jt];
		sumRealLL += AmpArrL[jt]*CosArr[jt];
		sumImagLL += AmpArrL[jt]*SinArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrR[jt]+AmpArrL[jt])*CosArr[jt];
		  im = 0.5*(AmpArrR[jt]+AmpArrL[jt])*SinArr[jt];
		  sumRealQ += re * qgain1[jt];
		  sumImagQ += im * qgain1[jt];
		  sumRealU += re * ugain1[jt];
		  sumImagU += im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereModSpec:    /* Uniform sphere + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[8+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
					visData[ilocv]*visData[ilocv]) * ddata[7];
		  arg = MAX (arg, 0.1);
		  amp = specFact * ddata[3] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  FazArr[itcnt] = freqFact * (tx + ty + tz);
		  AmpArrR[itcnt] = amp * rgain1[iComp];
		  AmpArrL[itcnt] = amp * lgain1[iComp];
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }  /* end inner loop over components */

	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrR[jt]*CosArr[jt];
		sumImagRR += AmpArrR[jt]*SinArr[jt];
		sumRealLL += AmpArrL[jt]*CosArr[jt];
		sumImagLL += AmpArrL[jt]*SinArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrR[jt]+AmpArrL[jt])*CosArr[jt];
		  im = 0.5*(AmpArrR[jt]+AmpArrL[jt])*SinArr[jt];
		  sumRealQ += re * qgain1[jt];
		  sumImagQ += im * qgain1[jt];
		  sumRealU += re * ugain1[jt];
		  sumImagU += im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  default:
	    ObitThreadLock(in->thread);  /* Lock against other threads */
	    Obit_log_error(err, OBIT_Error,"%s Unknown Comp model type %d in %s",
			   routine, in->modType, in->name);
	    ObitThreadUnlock(in->thread); 
	    goto finish;
	  }; /* end switch by model type */
	  
	  modRealRR = sumRealRR;
	  modImagRR = sumImagRR;
	  modRealLL = sumRealLL;
	  modImagLL = sumImagLL;
	  if (doCrossPol) {
	    modRealQ = sumRealQ;
	    modImagQ = sumImagQ;
	    modRealU = sumRealU;
	    modImagU = sumImagU;
	  }
	  
	  /* Dividing? */
	  if (in->doDivide) {
	    /* Divide model - also correct weight */
	    wtRR = modRealRR * modRealRR + modImagRR * modImagRR;
	    modRealRR /= wtRR;
	    modImagRR /= wtRR;
	    wtRR = sqrt (wtRR);
	    wtLL = modRealLL * modRealLL + modImagLL * modImagLL;
	    modRealLL /= wtLL;
	    modImagLL /= wtLL;
	    wtLL = sqrt (wtLL);
	  }
	  
	  /* RR */
	  iStoke = 0;
	  offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	  
	  /* Ignore blanked data unless replacing the data */
	  if ((visData[offset+2]>0.0) || in->doReplace) {
	    /* Apply model to data */
	    if (in->doDivide) {
	      temp = modRealRR * visData[offset] + modImagRR * visData[offset+1];
	      visData[offset+1] = modRealRR * visData[offset+1] - modImagRR * visData[offset];
	      visData[offset]   = temp;
	      visData[offset+2] *= wtRR;  /* correct weight */
	    } else if (in->doReplace) {  /* replace data with model */
	      visData[offset]   = modRealRR;
	      visData[offset+1] = modImagRR;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealRR;
	      visData[offset+1] -= modImagRR;
	    }
	  } /* end RR not blanked */
	  
	  /* LL */
	  offset += jincs;
	  /* Ignore blanked data unless replacing the data */
	  if ((visData[offset+2]>0.0) || in->doReplace) {
	    /* Apply model to data */
	    if (in->doDivide) {
	      temp = modRealLL * visData[offset] + modImagLL * visData[offset+1];
	      visData[offset+1] = modRealLL * visData[offset+1] - modImagLL * visData[offset];
	      visData[offset]   = temp;
	      visData[offset+2] *= wtLL;  /* correct weight */
	    } else if (in->doReplace) {  /* replace data with model */
	      visData[offset]   = modRealLL;
	      visData[offset+1] = modImagLL;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealLL;
	      visData[offset+1] -= modImagLL;
	    }
	  } /* end LL not blanked */
	  
	  if (doCrossPol) {
	    /* RL = Q+iU */
	    iStoke = 2;
	    offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	    
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
		/* TMS visData[offset]   = modRealU + modImagQ;
		   visData[offset+1] = modImagU - modRealQ;*/
                visData[offset]   = modRealQ - modImagU;
                visData[offset+1] = modImagQ + modRealU;
	      } else {
		/* Subtract model */
                visData[offset]   -= (modRealQ - modImagU);
                visData[offset+1] -= (modImagQ + modRealU);
	      }
	    } /* end RL not blanked */
	    
	    /* LR  = Q-iU */
	    offset += jincs;
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
		/* TMS visData[offset]   = (-modRealU + modImagQ);
		   visData[offset+1] = (-modImagU - modRealQ);*/
		visData[offset]   = modRealQ + modImagU;
                visData[offset+1] = modImagQ - modRealU;
	      } else {
		/* Subtract model */
		visData[offset]   -= (modRealQ + modImagU);
                visData[offset+1] -= (modImagQ - modRealU);
	      }
	    } /* end LR not blanked */
	  } /* end crosspol */
	  
	  offsetChannel += jincf;
	  channel++; /* Finished another channel */
	  /* Have we finished this plane in the correction cubes? */
	  if (plane!=in->FreqPlane[MIN(channel, (in->numUVChann-1))]) {
	    /* Reset gains & channels if this the last vis */
	    if (iVis>=(hiVis-1)) {
	      largs->endVMModelTime = -1.0e20;  
	      lstartChannel = iChannel;
	      lstartIF      = iIF;
	      if (iChannel==(startChannel+numberChannel-1)) {  /* end of channel loop */
		lstartChannel = 0;
		lstartIF++;
	      }
	    }
	    plane   = in->FreqPlane[MIN(channel, (in->numUVChann-1))];  /* Which plane in correction cube */
	    goto newPlane;
	  } /* end if new channel */
	} /* end loop over Channel */
	offsetIF += jincif;
      } /* end loop over IF */
      
      /* Come to here is finished with correction plane */
    newPlane:
      visData += lrec; /* Update vis pointer */
    } /* end loop over visibilities */
    iStoke = 0;  /* DEBUG */
  } /* end outer frequency loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadSkyModelVMBeamFTDFT */

/**
 * Do Fourier transform using a DFT for a buffer of data.
 * Version for time/spatial dependent effects with phase corrections .
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * If doFlip member is true the Fourier transform is multiplied by sqrt(-1)
 * (for Stokes RL and LR)
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Method assumes same correction to all antennas.
 * After the AIPSish QXXPTS, QPTDIV and friends
 * Arguments are given in the VMBeamFTFuncArg structure passed as arg starting 
 * with the following:
 * \li type   String identifying structure
 * \li in     SkyModelVM with model components loaded (ObitSkyModelLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number, <0-> no threads
 * \li err    Obit error stack object.
 * \li begVMModelTime Start time (days) of validity of model
 * \li endVMModelTime End time (days) of validity of model
 * \li VMComps Thread copy of Components list - not used here
 * \li dimGain Dimension of Rgain
 * \li Rgain   Float array of time/spatially variable R component gain
 * \li Rgaini  Float array of time/spatially variable R imaginary component gain
 * \li Lgain   Float array of time/spatially variable L component gain
 * \li Lgaini  Float array of time/spatially variable L imaginary component gain
 * \li Qgain   Float array of time/spatially variable Q component gain
 * \li Qgaini  Float array of time/spatially variable Q imaginary component gain
 * \li Ugain   Float array of time/spatially variable U component gain
 * \li Ugaini  Float array of time/spatially variable U imaginary component gain
 * \li channel Current UV channel being processed (used in model update ).
 * \return NULL
 */
static gpointer ThreadSkyModelVMBeamFTDFTPh (gpointer args)
{
  /* Get arguments from structure */
  VMBeamFTFuncArg *largs = (VMBeamFTFuncArg*)args;
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)largs->in;
  /*olong field      = largs->field;*/
  ObitUV *uvdata   = largs->uvdata;
  olong loVis      = largs->first-1;
  olong hiVis      = largs->last;
  olong ithread    = MAX (0, largs->ithread);
  ObitErr *err     = largs->err;
  /*olong dimGain    = largs->dimGain;*/
  ofloat *Rgain    = largs->Rgain;
  ofloat *Lgain    = largs->Lgain;
  ofloat *Qgain    = largs->Qgain;
  ofloat *Ugain    = largs->Ugain;
  ofloat *Rgaini   = largs->Rgaini;
  ofloat *Lgaini   = largs->Lgaini;
  ofloat *Qgaini   = largs->Qgaini;
  ofloat *Ugaini   = largs->Ugaini;

  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, ncomp;
  olong lrec, nrparm, naxis[2], channel, plane;
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, ilocb, suba, itemp, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealQ,  sumImagQ,  modRealQ=0.0,  modImagQ=0.0;
  ofloat sumRealU,  sumImagU,  modRealU=0.0,  modImagU=0.0;
  ofloat cbase, tx, ty, tz, ll, lll;
  ofloat re, im, specFact;
  ofloat *rgain1, *lgain1, *rgain2, *lgain2, *rgain1i, *lgain1i, *rgain2i, *lgain2i;
  ofloat *qgain1, *ugain1, *qgain2, *ugain2, *qgain1i, *ugain1i, *qgain2i, *ugain2i;
  ofloat amp, ampr, ampl, arg, freq2=0.0, freqFact, wtRR=0.0, wtLL=0.0, temp;
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArrRr[FazArrSize], AmpArrLr[FazArrSize], AmpArrRi[FazArrSize], AmpArrLi[FazArrSize];
  ofloat FazArr[FazArrSize];
  ofloat CosArr[FazArrSize], SinArr[FazArrSize];
  olong it, jt, itcnt;
  gboolean doCrossPol;
  odouble *freqArr;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMBeamFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;
  ilocb =  uvdata->myDesc->ilocb;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF       = in->startIFPB-1;
  numberIF      = MAX (1, in->numberIFPB);
  jincif        = uvdata->myDesc->incif;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  nUVchan       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocf];
  jincf         = uvdata->myDesc->incf;
  startPoln     = in->startPoln-1;
  numberPoln    = in->numberPoln;
  nUVpoln       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocs];
  jincs         = uvdata->myDesc->incs;  /* increment in real array */
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF beforefreq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  /* Cross or only parallel pol? */
  doCrossPol = (nUVpoln > 2) && in->doCrossPol;
  /* Only parallel for divide */
  if (in->doDivide) doCrossPol = FALSE;

  /* Get pointer for components */
  naxis[0] = 0; naxis[1] = 0; 
  Data = ObitFArrayIndex(in->comps, naxis);
  lcomp = in->comps->naxis[0];   /* Length of row in comp table */
  ncomp = in->comps->naxis[1];   /* Number of components */
  mcomp = in->numComp;           /* Actual number */

  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

  /* Current channel (0-rel) */
  channel = 0;
  plane   = in->FreqPlane[largs->channel];  /* Which plane in correction cube */

  /* Outer loop over blocks of channels */
  /* Starting parameters this pass */
  lstartIF       = startIF;
  lstartChannel  = startChannel;
  while (channel<(numberIF*numberChannel)) {
    
    /* Loop over vis in buffer */
    lrec    = uvdata->myDesc->lrec;         /* Length of record */
    visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
    nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */
    
    for (iVis=loVis; iVis<hiVis; iVis++) {
      
      /* Is current model still valid? */
      if ((visData[iloct] > largs->endVMModelTime) || (visData[iloct] < largs->begVMModelTime)) {
	/* Current channel - which plane in correction to apply? */
	largs->channel = channel;
	/* Subarray 0-rel */
	itemp = (olong)visData[ilocb];
	suba = 100.0 * (visData[ilocb]-itemp) + 0.5; 
	/* Update */
	myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba, uvdata, ithread, err);
	if (err->error) {
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
			 routine);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}
      }
      
      /* Need antennas numbers */
      cbase = visData[ilocb]; /* Baseline */
      ant1 = (cbase / 256.0) + 0.001;
      ant2 = (cbase - ant1 * 256) + 0.001;
      ant1--;    /* 0 rel */
      ant2--;    /* 0 rel */
      
      /* Loop over IFs */
      channel = lstartIF* nUVchan + lstartChannel; /* UV Channel */
      for (iIF=lstartIF; iIF<startIF+numberIF; iIF++) {
	offsetIF = nrparm + iIF*jincif; 
	for (iChannel=lstartChannel; iChannel<startChannel+numberChannel; iChannel++) {
	  offsetChannel = offsetIF + iChannel*jincf; 
	  freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */
	  
	  /* Sum over components */
	  /* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	  sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	  sumRealQ  = sumImagQ  = sumRealU  = sumImagU  = 0.0;
	  /* Set component gain lists by antenna and type */
	  ddata = Data;
	  rgain1 = Rgain;
	  rgain2 = Rgain;
	  lgain1 = Lgain;
	  lgain2 = Lgain;
	  qgain1 = Qgain;
	  qgain2 = Qgain;
	  ugain1 = Ugain;
	  ugain2 = Ugain;
	  /* Imaginary parts */
	  rgain1i = Rgaini;
	  rgain2i = Rgaini;
	  lgain1i = Lgaini;
	  lgain2i = Lgaini;
	  qgain1i = Qgaini;
	  qgain2i = Qgaini;
	  ugain1i = Ugaini;
	  ugain2i = Ugaini;

	  /* Sum by model type - assume phase same for RR, LL */
	  switch (in->modType) {
	  case OBIT_SkyModel_PointMod:     /* Point */
	    /* From the AIPSish QXXPTS.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArrRr[itcnt] = ddata[3] * rgain1[iComp];
		AmpArrLr[itcnt] = ddata[3] * lgain1[iComp];
		AmpArrRi[itcnt] = ddata[3] * rgain1i[iComp];
		AmpArrLi[itcnt] = ddata[3] * lgain1i[iComp];
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt] - AmpArrRi[jt]*SinArr[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt] + AmpArrRi[jt]*CosArr[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt] - AmpArrLi[jt]*SinArr[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt] + AmpArrLi[jt]*CosArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_PointModSpec:     /* Point + spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[6+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  FazArr[itcnt]  = freqFact * (tx + ty + tz);
		  AmpArrRr[itcnt] = ddata[3] * rgain1[iComp];
		  AmpArrLr[itcnt] = ddata[3] * lgain1[iComp];
		  AmpArrRi[itcnt] = ddata[3] * rgain1i[iComp];
		  AmpArrLi[itcnt] = ddata[3] * lgain1i[iComp];
		}  /* end if valid */
		ddata += lcomp;  /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);

	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		/* Parallel  pol */
		sumRealRR += AmpArrRr[jt]*CosArr[jt] - AmpArrRi[jt]*SinArr[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt] + AmpArrRi[jt]*CosArr[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt] - AmpArrLi[jt]*SinArr[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt] + AmpArrLi[jt]*CosArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	    /* From the AIPSish QGASUB.FOR  */
	    freq2 = freqArr[iIF*kincif + iChannel*kincf];    /* Frequency squared */
	    freq2 *= freq2;
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			       ddata[8]*visData[ilocv]*visData[ilocv] +
			       ddata[9]*visData[ilocu]*visData[ilocv]);
		if (arg>1.0e-5) ampr *= exp (arg);
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			       ddata[8]*visData[ilocv]*visData[ilocv] +
			       ddata[9]*visData[ilocu]*visData[ilocv]);
		ampr = ddata[3] * rgain1[iComp];
		ampl = ddata[3] * lgain1[iComp];
		if (arg>1.0e-5) {ampl *= exp (arg); ampr *= exp (arg);}
		AmpArrRr[itcnt] = ampr;
		AmpArrLr[itcnt] = ampl;
		ampl = ddata[3] * lgain1i[iComp];
		ampr = ddata[3] * rgain1i[iComp];
		if (arg>1.0e-5) {ampl *= exp (arg); ampr *= exp (arg);}
		AmpArrRi[itcnt] = ampr;
		AmpArrLi[itcnt] = ampl;
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }  /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt] - AmpArrRi[jt]*SinArr[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt] + AmpArrRi[jt]*CosArr[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt] - AmpArrLi[jt]*SinArr[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt] + AmpArrLi[jt]*CosArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	    freq2 = freqFact*freqFact;    /* Frequency factor squared */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[9+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
				 ddata[8]*visData[ilocv]*visData[ilocv] +
				 ddata[9]*visData[ilocu]*visData[ilocv]);
		  if (arg<-1.0e-5) amp = specFact * ddata[3] * exp (arg);
		  else amp = specFact * ddata[3];
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  FazArr[itcnt] = freqFact * (tx + ty + tz);
		  AmpArrRr[itcnt] = ddata[3] * rgain1[iComp];
		  AmpArrLr[itcnt] = ddata[3] * lgain1[iComp];
		  AmpArrRi[itcnt] = ddata[3] * rgain1i[iComp];
		  AmpArrLi[itcnt] = ddata[3] * lgain1i[iComp];
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }  /* end inner loop over components */

	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt] - AmpArrRi[jt]*SinArr[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt] + AmpArrRi[jt]*CosArr[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt] - AmpArrLi[jt]*SinArr[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt] + AmpArrLi[jt]*CosArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	    /* From the AIPSish QSPSUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ddata[7];
		arg = MAX (arg, 0.1);
		arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		AmpArrRr[itcnt] = ddata[3] * rgain1[iComp];
		AmpArrLr[itcnt] = ddata[3] * lgain1[iComp];
		AmpArrRi[itcnt] = ddata[3] * rgain1i[iComp];
		AmpArrLi[itcnt] = ddata[3] * lgain1i[iComp];
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt] - AmpArrRi[jt]*SinArr[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt] + AmpArrRi[jt]*CosArr[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt] - AmpArrLi[jt]*SinArr[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt] + AmpArrLi[jt]*CosArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereModSpec:    /* Uniform sphere + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      for (iComp=it; iComp<mcomp; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[8+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
					visData[ilocv]*visData[ilocv]) * ddata[7];
		  arg = MAX (arg, 0.1);
		  amp = specFact * ddata[3] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  FazArr[itcnt] = freqFact * (tx + ty + tz);
		  AmpArrRr[itcnt] = ddata[3] * rgain1[iComp];
		  AmpArrLr[itcnt] = ddata[3] * lgain1[iComp];
		  AmpArrRi[itcnt] = ddata[3] * rgain1i[iComp];
		  AmpArrLi[itcnt] = ddata[3] * lgain1i[iComp];
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
		if (itcnt>=FazArrSize) break;
	      }  /* end inner loop over components */

	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt] - AmpArrRi[jt]*SinArr[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt] + AmpArrRi[jt]*CosArr[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt] - AmpArrLi[jt]*SinArr[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt] + AmpArrLi[jt]*CosArr[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  default:
	    ObitThreadLock(in->thread);  /* Lock against other threads */
	    Obit_log_error(err, OBIT_Error,"%s Unknown Comp model type %d in %s",
			   routine, in->modType, in->name);
	    ObitThreadUnlock(in->thread); 
	    goto finish;
	  }; /* end switch by model type */
	  
	  /* Models from sums */
	  modRealRR = sumRealRR;
	  modImagRR = sumImagRR;
	  modRealLL = sumRealLL;
	  modImagLL = sumImagLL;
	  if (doCrossPol) {
	    modRealQ = sumRealQ;
	    modImagQ = sumImagQ;
	    modRealU = sumRealU;
	    modImagU = sumImagU;
	  }
	  
	  /* Dividing? */
	  if (in->doDivide) {
	    /* Divide model - also correct weight */
	    wtRR = modRealRR * modRealRR + modImagRR * modImagRR;
	    modRealRR /= wtRR;
	    modImagRR /= wtRR;
	    wtRR = sqrt (wtRR);
	    wtLL = modRealLL * modRealLL + modImagLL * modImagLL;
	    modRealLL /= wtLL;
	    modImagLL /= wtLL;
	    wtLL = sqrt (wtLL);
	  }
	  
	  /* RR */
	  iStoke = 0;
	  offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	  
	  /* Ignore blanked data unless replacing the data */
	  if ((visData[offset+2]>0.0) || in->doReplace) {
	    /* Apply model to data */
	    if (in->doDivide) {
	      temp = modRealRR * visData[offset] + modImagRR * visData[offset+1];
	      visData[offset+1] = modRealRR * visData[offset+1] - modImagRR * visData[offset];
	      visData[offset]   = temp;
	      visData[offset+2] *= wtRR;  /* correct weight */
	    } else if (in->doReplace) {  /* replace data with model */
	      visData[offset]   = modRealRR;
	      visData[offset+1] = modImagRR;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealRR;
	      visData[offset+1] -= modImagRR;
	    }
	  } /* end RR not blanked */
	  
	  /* LL */
	  offset += jincs;
	  /* Ignore blanked data unless replacing the data */
	  if ((visData[offset+2]>0.0) || in->doReplace) {
	    /* Apply model to data */
	    if (in->doDivide) {
	      temp = modRealLL * visData[offset] + modImagLL * visData[offset+1];
	      visData[offset+1] = modRealLL * visData[offset+1] - modImagLL * visData[offset];
	      visData[offset]   = temp;
	      visData[offset+2] *= wtLL;  /* correct weight */
	    } else if (in->doReplace) {  /* replace data with model */
	      visData[offset]   = modRealLL;
	      visData[offset+1] = modImagLL;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealLL;
	      visData[offset+1] -= modImagLL;
	    }
	  } /* end LL not blanked */
	  
	  if (doCrossPol) {
	    /* RL = Q+iU */
	    iStoke = 2;
	    offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	    
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
		/* TMS visData[offset]   = modRealU + modImagQ;
		   visData[offset+1] = modImagU - modRealQ;*/
                visData[offset]   = modRealQ - modImagU;
                visData[offset+1] = modImagQ + modRealU;
	      } else {
		/* Subtract model */
                visData[offset]   -= (modRealQ - modImagU);
                visData[offset+1] -= (modImagQ + modRealU);
	      }
	    } /* end RL not blanked */
	    
	    /* LR  = Q-iU */
	    offset += jincs;
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
		/* TMS visData[offset]   = (-modRealU + modImagQ);
		   visData[offset+1] = (-modImagU - modRealQ);*/
		visData[offset]   = modRealQ + modImagU;
                visData[offset+1] = modImagQ - modRealU;
	      } else {
		/* Subtract model */
		visData[offset]   -= (modRealQ + modImagU);
                visData[offset+1] -= (modImagQ - modRealU);
	      }
	    } /* end LR not blanked */
	  } /* end crosspol */
	  
	  offsetChannel += jincf;
	  channel++; /* Finished another channel */
	  /* Have we finished this plane in the correction cubes? */
	  if (plane!=in->FreqPlane[MIN(channel, (in->numUVChann-1))]) {
	    /* Reset gains & channels if this the last vis */
	    if (iVis>=(hiVis-1)) {
	      largs->endVMModelTime = -1.0e20;  
	      lstartChannel = iChannel;
	      lstartIF      = iIF;
	      if (iChannel==(startChannel+numberChannel-1)) {  /* end of channel loop */
		lstartChannel = 0;
		lstartIF++;
	      }
	    }
	    plane   = in->FreqPlane[MIN(channel, (in->numUVChann-1))];  /* Which plane in correction cube */
	    goto newPlane;
	  } /* end if new channel */
	} /* end loop over Channel */
	offsetIF += jincif;
      } /* end loop over IF */
      
      /* Come to here is finished with correction plane */
    newPlane:
      visData += lrec; /* Update vis pointer */
    } /* end loop over visibilities */
    iStoke = 0;  /* DEBUG */
  } /* end outer frequency loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadSkyModelVMBeamFTDFTPh */

/** 
 *  Get model frequency primary beam 
 * \param desc    model image header
 * \param Angle   Angle from the pointing position (deg)
 * \param x       x offset from center (deg)
 * \param y       y offset from center (deg)
 * \param antSize Antenna diameter in meters. (defaults to 25.0)
 * \param pbmin   Minimum antenna gain Jinc 0=>0.05, poly 0=> 0.01
 * \return Relative gain at freq refFreq wrt average of Freq.
*/
static ofloat getPBBeam(ObitImageDesc *desc, ofloat x, ofloat y, 
			ofloat antSize, ofloat pbmin)
{
  ofloat PBBeam = 1.0;
  odouble freq;
  gboolean doJinc;
  odouble Angle, RAPnt, DecPnt, ra, dec, xx, yy, zz;

  if (antSize <=0.0) antSize = 25.0;  /* defailt antenna size */

  /* Get pointing position */
  ObitImageDescGetPoint(desc, &RAPnt, &DecPnt);
  RAPnt  *= DG2RAD;
  DecPnt *= DG2RAD;

  /* Convert offset to position */
  ObitSkyGeomXYShift (desc->crval[desc->jlocr], desc->crval[desc->jlocd],
		      x, y, ObitImageDescRotate(desc), &ra, &dec);

  /* Angle from center */
  xx = DG2RAD * ra;
  yy = DG2RAD * dec;
  zz = sin (yy) * sin (DecPnt) + cos (yy) * cos (DecPnt) * cos (xx-RAPnt);
  zz = MIN (zz, 1.000);
  Angle = acos (zz) * RAD2DG;

  /* Which beam shape function to use? */
  freq = desc->crval[desc->jlocf];
  doJinc = (freq >= 1.0e9);
  
  /* Gain at freq */
  if (doJinc) PBBeam = ObitPBUtilJinc(Angle, freq, antSize, pbmin);
  else        PBBeam = ObitPBUtilPoly(Angle, freq, pbmin);
 
  return PBBeam;
} /* end  getPBBeam*/

