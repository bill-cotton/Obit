/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2020                                          */
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

#include "ObitUVDesc.h"
#include "ObitThread.h"
#include "ObitBeamShape.h"
#include "ObitSkyModelVMBeam.h"
#include "ObitTableCCUtil.h"
#include "ObitTableAN.h"
#include "ObitFFT.h"
#include "ObitUVUtil.h"
#include "ObitImageUtil.h"
#include "ObitPBUtil.h"
#include "ObitMem.h"
#include "ObitPrecess.h"
#include "ObitTableANUtil.h"
#include "ObitSkyGeom.h"
#include "ObitSinCos.h"
#include "ObitExp.h"
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
static ofloat getPBBeam(ObitBeamShape *beamShape, ObitImageDesc *desc, 
			ofloat x, ofloat y, 
			ofloat antSize, odouble freq, ofloat pbmin);

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
  /* UV Interpolator for FTGrid */
  ObitCInterpolate *Interp;
  /* VM class entries */
  /* Start time (days) of validity of model */
  ofloat begVMModelTime;
  /* End time (days) of validity of model */
  ofloat endVMModelTime;
  /* Thread copy of Components list*/
  ObitFArray *VMComps;
  /* VMBeam class entries */
  /* Amp, phase interpolator for R/X pol Beam image array */
  ObitFInterpolate **BeamRXInterp, **BeamRXPhInterp;
  /* Amp, phase interpolator for L/Y pol Beam image array */
  ObitFInterpolate **BeamLYInterp, **BeamLYPhInterp;
  /* Amp, phase interpolator for RL/XY pol Beam image array */
  ObitFInterpolate **BeamRLInterp, **BeamRLPhInterp;
  /* Amp, phase interpolator for LR/YX pol Beam image array */
  ObitFInterpolate **BeamLRInterp, **BeamLRPhInterp;
  /** Number of antenna types */
  olong numAntType;
  /** Current uv channel number being processed.  */
  olong channel;
  /** Frequency of desired beam (Hz) corresponds to channel */
  odouble  BeamFreq;
  /** Dimension of Rgain...  */
  olong dimGain;
  /** Arrays of time/spatially variable R/X component gain, real, imag */
  ofloat **Rgain, **Rgaini;
  /** Arrays of time/spatially variable L/Y component gain, real, imag */
  ofloat **Lgain, **Lgaini;
  /** Arrays of time/spatially variable RL/XY component gain, real, imag */
  ofloat **RLgain, **RLgaini;
  /** Arrays of time/spatially variable LR/YX component gain, real, imag */
  ofloat **LRgain, **LRgaini;
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
  olong i, j;
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
  /* Component arrays per antenna type */
  out->dimGain    = in->dimGain;
  out->numAntType = in->numAntType; 
  if (out->Rgain==NULL) out->Rgain  = g_malloc0(in->numAntType*sizeof(ofloat*));
  if (out->Lgain==NULL) out->Lgain  = g_malloc0(in->numAntType*sizeof(ofloat*));
  if (out->Rgaini==NULL) out->Rgaini = g_malloc0(in->numAntType*sizeof(ofloat*));
  if (out->Lgaini==NULL) out->Lgaini = g_malloc0(in->numAntType*sizeof(ofloat*));
  for (j=0; j<in->numAntType; j++) {
    if (in->dimGain>0) {
      if (out->Rgain[j]) g_free(out->Rgain[j]); out->Rgain[j] = NULL;
      out->Rgain[j] = g_malloc0(in->dimGain*sizeof(ofloat));
      for (i=0; i<in->dimGain; i++) out->Rgain[j][i] = in->Rgain[j][i];
      if (out->Rgaini[j]) g_free(out->Rgaini[j]); out->Rgaini[j] = NULL;
      if (in->Rgaini[j]) {
	out->Rgaini[j] = g_malloc0(in->dimGain*sizeof(ofloat));
	for (i=0; i<in->dimGain; i++) out->Rgaini[j][i] = in->Rgaini[j][i];
      }
      
      if (out->Lgain[j]) g_free(out->Lgain[j]); out->Lgain[j] = NULL;
      out->Lgain[j] = g_malloc0(in->dimGain*sizeof(ofloat));
      for (i=0; i<in->dimGain; i++) out->Lgain[j][i] = in->Lgain[j][i];
      if (out->Lgaini[j]) g_free(out->Lgaini[j]); out->Lgaini[j] = NULL;
      if (in->Lgaini[j]) {
	out->Lgaini[j] = g_malloc0(in->dimGain*sizeof(ofloat));
	for (i=0; i<in->dimGain; i++) out->Lgaini[j][i] = in->Lgaini[j][i];
      }
      if (in->doCrossPol) {
	if (out->RLgain[j]) g_free(out->RLgain[j]); out->RLgain[j] = NULL;
	out->RLgain[j] = g_malloc0(in->dimGain*sizeof(ofloat));
	for (i=0; i<in->dimGain; i++) out->RLgain[j][i] = in->RLgain[j][i];
	if (in->RLgaini[j]) {
	  if (out->RLgaini[j]) g_free(out->RLgaini[j]); out->RLgaini[j] = NULL;
	  out->RLgaini[j] = g_malloc0(in->dimGain*sizeof(ofloat));
	  for (i=0; i<in->dimGain; i++) out->RLgaini[j][i] = in->RLgaini[j][i];
	}
	
	if (out->LRgain) g_free(out->LRgain); out->LRgain = NULL;
	out->LRgain = g_malloc0(in->dimGain*sizeof(ofloat));
	for (i=0; i<in->dimGain; i++) out->LRgain[i] = in->LRgain[i];
	if (in->LRgaini) {
	  if (out->LRgaini) g_free(out->LRgaini); out->LRgaini = NULL;
	  out->LRgaini = g_malloc0(in->dimGain*sizeof(ofloat));
	  for (i=0; i<in->dimGain; i++) out->LRgaini[i] = in->LRgaini[i];
	}
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
 * \param numAntType number of antenna types
 * \param RXBeam   R/X Beam normalized image array per type
 * \param LYBeam   L/Y Beam normalized image array per type
 * \param RLBeam   RL/XY Beam fractional image array per type if nonNULL
 * \param LRBeam   LR/YX Beam fractional image array per type if nonNULL
 * \param RXBeamIm R/X Beam phase image array per type if nonNULL
 * \param LYBeamIm L/Y Beam phase image array per type if nonNULL
 * \param RLBeamIm RL/XY Beam phase image array per type if nonNULL
 * \param LRBeamIm L/Y  Beam phase image array per type if nonNULL
 * \param Diams    Antenna diameters (m) per type
 * \return the new object.
 */
ObitSkyModelVMBeam* 
ObitSkyModelVMBeamCreate (gchar* name, ObitImageMosaic* mosaic,
			  ObitUV *uvData, olong numAntType,
			  ObitImage **RXBeam,   ObitImage **LYBeam, 
			  ObitImage **RLBeam,   ObitImage **LRBeam, 
			  ObitImage **RXBeamIm, ObitImage **LYBeamIm, 
			  ObitImage **RLBeamIm, ObitImage **LRBeamIm, 
			  ofloat *Diams, ObitErr *err)
{
  ObitSkyModelVMBeam* out=NULL;
  olong number, i, nchan, nif, refType;
  ofloat refDiam;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean doTab = TRUE, doCmplx;
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

  doCmplx = RXBeamIm && RXBeamIm[0];  /* Have phases? */
  out->doCmplx = doCmplx;
  /* Swallow Beam images */
  out->numPlane   = g_malloc0(numAntType*sizeof(olong));
  out->Diams      = g_malloc0(numAntType*sizeof(ofloat)); 
  out->AntType    = g_malloc0(uvData->myDesc->maxAnt*sizeof(olong)); 
  out->RXBeam     = g_malloc0(numAntType*sizeof(ObitImageInterp*)); 
  out->LYBeam     = g_malloc0(numAntType*sizeof(ObitImageInterp*)); 
  out->RLBeam     = g_malloc0(numAntType*sizeof(ObitImageInterp*)); 
  out->LRBeam     = g_malloc0(numAntType*sizeof(ObitImageInterp*)); 
  if (doCmplx) {
    out->RXBeamIm = g_malloc0(numAntType*sizeof(ObitImageInterp*));
    out->LYBeamIm = g_malloc0(numAntType*sizeof(ObitImageInterp*));
    out->RLBeamIm = g_malloc0(numAntType*sizeof(ObitImageInterp*));
    out->LRBeamIm = g_malloc0(numAntType*sizeof(ObitImageInterp*));
  }
  out->numAntType = numAntType;
  out->doCrossPol = TRUE;
  for (i=0; i<out->numAntType; i++) {
    out->RXBeam[i] = ObitImageInterpCreate("RXBeam", RXBeam[i], 2, err);
    out->LYBeam[i] = ObitImageInterpCreate("LYBeam", LYBeam[i], 2, err);
    if (RLBeam[i])
      out->RLBeam[i] = ObitImageInterpCreate("RLBeam", RLBeam[i], 2, err);
    else
      out->doCrossPol = FALSE;
    if (LRBeam[i])
      out->LRBeam[i] = ObitImageInterpCreate("LRBeam", LRBeam[i], 2, err);
    else
      out->doCrossPol = FALSE;
    if (err->error) Obit_traceback_val (err, routine, name, out);
    out->numPlane[i] = out->RXBeam[i]->nplanes;

    /* Phase beams */
    if (doCmplx && RXBeamIm && RXBeamIm[i])
      out->RXBeamIm[i] = ObitImageInterpCreate("RXBeamIm", RXBeamIm[i], 2, err);
    if (doCmplx && LYBeamIm && LYBeamIm[i])
      out->LYBeamIm[i] = ObitImageInterpCreate("LYBeamIm", LYBeamIm[i], 2, err);
    if (doCmplx && RLBeamIm && RLBeamIm[i])
      out->RLBeamIm[i] = ObitImageInterpCreate("RLBeamIm", RLBeamIm[i], 2, err);
    if (doCmplx && LRBeamIm && LRBeamIm[i])
      out->LRBeamIm[i] = ObitImageInterpCreate("LRBeamIm", LRBeamIm[i], 2, err);
    if (err->error) Obit_traceback_val (err, routine, name, out);
    
    /* Make sure they are all consistent */
    Obit_retval_if_fail ((ObitFArrayIsCompatable(out->RXBeam[i]->ImgPixels, 
						 out->LYBeam[i]->ImgPixels)), err, out,
			 "%s: Incompatable pp, qq beam arrays", routine);
    if (out->doCrossPol) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->LRBeam[i]->ImgPixels, 
						   out->RLBeam[i]->ImgPixels)), err, out,
			   "%s: Incompatable pq,qp, beam arrays", routine);
    }
    if (doCmplx && out->RXBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->RXBeam[i]->ImgPixels, 
						   out->RXBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable pp amp, phase beam arrays", routine);
    }
    if (doCmplx && out->LYBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->LYBeam[i]->ImgPixels, 
						   out->LYBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable qq amp, phase beam arrays", routine);
    }
    if (out->doCrossPol && out->RLBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->RLBeam[i]->ImgPixels, 
						   out->RLBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable pq amp, phase beam arrays", routine);
    }
    if (out->doCrossPol && out->LRBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->LRBeam[i]->ImgPixels, 
						   out->LRBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable qp amp, phase beam arrays", routine);
    }

    ObitInfoListAlwaysPut (RXBeam[i]->info, "doTab", OBIT_bool, dim, &doTab);
  
  } /* end loop over ant type */

  refType = 0; refDiam=0.0;  /* Reference type is type with largest diameter */
  /* Find largest Antenna type diameter */
  for (i=0; i<out->numAntType; i++) {
    /* Reference type */
    if (out->Diams[i] > refDiam) {
      refDiam = out->Diams[i]; refType = i;}
    }
  /* Reference Beam shape - Tabulated if possible */
  ObitInfoListAlwaysPut (RXBeam[refType]->info, "doTab", OBIT_bool, dim, &doTab);
  out->BeamShape = ObitBeamShapeCreate("Shape", RXBeam[refType], 0.01, 25.0, TRUE);

 /* Get list of planes per channel */
  nchan = uvData->myDesc->inaxes[uvData->myDesc->jlocf];
  if (uvData->myDesc->jlocif>=0) 
    nif = uvData->myDesc->inaxes[uvData->myDesc->jlocif];
  else nif = 1;
  out->numUVChann = nchan*nif;
  out->FreqPlane  = g_malloc0(out->numUVChann*sizeof(olong));
  for (i=0; i<out->numUVChann; i++) 
    out->FreqPlane[i] = MAX(0, MIN (out->numPlane[0]-1, 
				    ObitImageInterpFindPlane(out->RXBeam[0], uvData->myDesc->freqArr[i])));
  /* Release beam buffers */
  for (i=0; i<out->numAntType; i++) {
    if (RXBeam[i]  && (RXBeam[i]->image)) RXBeam[i]->image = ObitImageUnref(RXBeam[i]->image);
    if (RLBeam[i]  && (RLBeam[i]->image)) RLBeam[i]->image = ObitImageUnref(RLBeam[i]->image);
    if (LRBeam[i]  && (LRBeam[i]->image)) LRBeam[i]->image = ObitImageUnref(LRBeam[i]->image);
    if (LYBeam [i] && (LYBeam[i]->image)) LYBeam[i]->image = ObitImageUnref(LYBeam[i]->image);
    if (doCmplx) {
      if (RXBeamIm[i]  && (RXBeamIm[i]->image)) RXBeamIm[i]->image = ObitImageUnref(RXBeamIm[i]->image);
      if (RLBeamIm[i]  && (RLBeamIm[i]->image)) RLBeamIm[i]->image = ObitImageUnref(RLBeamIm[i]->image);
      if (LRBeamIm[i]  && (LRBeamIm[i]->image)) LRBeamIm[i]->image = ObitImageUnref(LRBeamIm[i]->image);
      if (LYBeamIm[i]  && (LYBeamIm[i]->image)) LYBeamIm[i]->image = ObitImageUnref(LYBeamIm[i]->image);
    }
  }
  /* Set antenna Types */
  ObitSkyModelVMBeamSetAnt (out, uvData, err);

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
  olong numAntList;
  ObitTableList *list=NULL;
  ObitTableAN *TableAN=NULL;
  ObitUVDesc *uvDesc;
  ofloat phase=0.5, cp, sp;
  /*gchar *blank="    ";*/
  olong i, j;
  VMBeamFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMBeamInitMod";

  if (err->error) return;

  /* Save/reset calibration state */
  /* It's not clear why it was doing this - not always wanted*/
  in->saveDoCalSelect = FALSE;
  ObitInfoListGetTest(uvdata->info, "doCalSelect", &type, dim, &in->saveDoCalSelect);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  /*btemp = FALSE;
    ObitInfoListAlwaysPut (uvdata->info, "doCalSelect", OBIT_bool, dim, &btemp);*/
  in->saveDoCalib = 0;
  ObitInfoListGetTest(uvdata->info, "doCalib", &type, dim, &in->saveDoCalib);
  dim[0] = dim[1] = dim[2] = 1;
  /*itemp = -1;
    ObitInfoListAlwaysPut (uvdata->info, "doCalib", OBIT_long, dim, &itemp);*/
  strncpy (in->saveStokes, "    ", 4);
  ObitInfoListGetTest(uvdata->info, "Stokes", &type, dim, in->saveStokes);
  dim[0] = 4; dim[1] = dim[2] = 1;
  /* ObitInfoListAlwaysPut (uvdata->info, "Stokes", OBIT_string, dim, blank);*/

  /* Open/close to reset */
  ObitUVOpen (uvdata, OBIT_IO_ReadOnly, err);
  ObitUVClose (uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(VMBeamFTFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(VMBeamFTFuncArg)); 
  
    for (i=0; i<in->nThreads; i++) {
      args = (VMBeamFTFuncArg*)in->threadArgs[i];
      strcpy (args->type, "vmbeam");  /* Enter type as first entry */
      args->in     = inn;
      args->uvdata = uvdata;
      args->ithread = i;
      args->err    = err;
      args->numAntType     = in->numAntType;
      args->BeamRXInterp   = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      args->BeamLYInterp   = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      args->BeamRLInterp   = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      args->BeamLRInterp   = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      args->BeamRXPhInterp = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      args->BeamLYPhInterp = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      args->BeamRLPhInterp = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      args->BeamLRPhInterp = g_malloc0(in->numAntType*sizeof(ObitImageInterp*));
      for (j=0; j<in->numAntType; j++) {
	if (in->RXBeam[j]) args->BeamRXInterp[j] = ObitImageInterpCloneInterp(in->RXBeam[j],err);
	if (in->LYBeam[j]) args->BeamLYInterp[j] = ObitImageInterpCloneInterp(in->LYBeam[j],err);
	if (in->RLBeam[j]) args->BeamRLInterp[j] = ObitImageInterpCloneInterp(in->RLBeam[j],err);
	if (in->LRBeam[j]) args->BeamLRInterp[j] = ObitImageInterpCloneInterp(in->LRBeam[j],err);
	if (in->doCmplx) {
	  if (in->RXBeamIm[j]) args->BeamRXPhInterp[j] = ObitImageInterpCloneInterp(in->RXBeamIm[j],err);
	  if (in->LYBeamIm[j]) args->BeamLYPhInterp[j] = ObitImageInterpCloneInterp(in->LYBeamIm[j],err);
	  if (in->RLBeamIm[j]) args->BeamRLPhInterp[j] = ObitImageInterpCloneInterp(in->RLBeamIm[j],err);
	  if (in->LRBeamIm[j]) args->BeamLRPhInterp[j] = ObitImageInterpCloneInterp(in->LRBeamIm[j],err);
	}
      } /* end antenna type loop */
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      args->begVMModelTime = -1.0e20;
      args->endVMModelTime = -1.0e20;
      args->VMComps = NULL;
      args->dimGain = 0;
      args->Rgain   = g_malloc0(in->numAntType*sizeof(ofloat*));
      args->Rgaini  = g_malloc0(in->numAntType*sizeof(ofloat*));
      args->Lgain   = g_malloc0(in->numAntType*sizeof(ofloat*));
      args->Lgaini  = g_malloc0(in->numAntType*sizeof(ofloat*));
      args->RLgain  = g_malloc0(in->numAntType*sizeof(ofloat*));
      args->RLgaini = g_malloc0(in->numAntType*sizeof(ofloat*));
      args->LRgain  = g_malloc0(in->numAntType*sizeof(ofloat*));
      args->LRgaini = g_malloc0(in->numAntType*sizeof(ofloat*));
    }
  } /* end initialize */

  /* Call parent initializer */
  ObitSkyModelVMInitMod(inn, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fourier transform routines - DFT only */
  /* Are phases given? */
  if (in->RXBeamIm) 
    in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamFTDFTPh;
  else /* No phase */
    in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamFTDFT;

  /* Check requested Stokes
  Obit_return_if_fail((!strncmp(in->stokes,"    ",4)), err,
		      "%s: Unsupported Stokes %s", routine, in->stokes); */

  /* Check that data contains RR, LL or (XX,YY)*/
  uvDesc = uvdata->myDesc;
  Obit_return_if_fail((((uvDesc->crval[uvDesc->jlocs]==-1.0) || 
			(uvDesc->crval[uvDesc->jlocs]==-5.0)) && 
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
				   OBIT_IO_ReadOnly, 0, 0, 0, err);
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
  ObitExpInit();

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
  olong i, j;
  VMBeamFTFuncArg *args;

  /* Call parent shutdown */
  ObitSkyModelVMShutDownMod(inn, uvdata, err);

  if (in->threadArgs) {
    /* Check type - only handle "vmbeam" */
    args = (VMBeamFTFuncArg*)in->threadArgs[0];
    if ((strlen(args->type)>6) || (!strncmp(args->type, "vmbeam", 6))) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMBeamFTFuncArg*)in->threadArgs[i];
	for (j=0; j<args->numAntType; j++) {
	  args->BeamRXInterp[j] = ObitFInterpolateUnref(args->BeamRXInterp[j]);
	  args->BeamLYInterp[j] = ObitFInterpolateUnref(args->BeamLYInterp[j]);
	  args->BeamRLInterp[j] = ObitFInterpolateUnref(args->BeamRLInterp[j]);
	  args->BeamLRInterp[j] = ObitFInterpolateUnref(args->BeamLRInterp[j]);
	  args->BeamRXPhInterp[j] = ObitFInterpolateUnref(args->BeamRXPhInterp[j]);
	  args->BeamLYPhInterp[j] = ObitFInterpolateUnref(args->BeamLYPhInterp[j]);
	  args->BeamRLPhInterp[j] = ObitFInterpolateUnref(args->BeamRLPhInterp[j]);
	  args->BeamLRPhInterp[j] = ObitFInterpolateUnref(args->BeamLRPhInterp[j]);
	  if (args->Rgain[j])   g_free(args->Rgain[j]);
	  if (args->Lgain[j])   g_free(args->Lgain[j]);
	  if (args->RLgain[j])  g_free(args->RLgain[j]);
	  if (args->LRgain[j])  g_free(args->LRgain[j]);
	  if (args->Rgaini[j])  g_free(args->Rgaini[j]);
	  if (args->Lgaini[j])  g_free(args->Lgaini[j]);
	  if (args->RLgaini[j]) g_free(args->RLgaini[j]);
	  if (args->LRgaini[j]) g_free(args->LRgaini[j]);
	} /* end ant type loop */
	if (args->BeamRXInterp)   {g_free(args->BeamRXInterp);}   args->BeamRXInterp = NULL;
	if (args->BeamLYInterp)   {g_free(args->BeamLYInterp);}   args->BeamLYInterp = NULL;
	if (args->BeamRLInterp)   {g_free(args->BeamRLInterp);}   args->BeamRLInterp = NULL;
	if (args->BeamLRInterp)   {g_free(args->BeamLRInterp);}   args->BeamLRInterp = NULL;
	if (args->BeamRXPhInterp) {g_free(args->BeamRXPhInterp);} args->BeamRXPhInterp = NULL;
	if (args->BeamLYPhInterp) {g_free(args->BeamLYPhInterp);} args->BeamLYPhInterp = NULL;
	if (args->BeamRLPhInterp) {g_free(args->BeamRLPhInterp);} args->BeamRLPhInterp = NULL;
	if (args->BeamLRPhInterp) {g_free(args->BeamLRPhInterp);} args->BeamLRPhInterp = NULL;
	if (args->Rgain)   g_free(args->Rgain);
	if (args->Lgain)   g_free(args->Lgain);
	if (args->RLgain)  g_free(args->RLgain);
	if (args->LRgain)  g_free(args->LRgain);
	if (args->Rgaini)  g_free(args->Rgaini);
	if (args->Lgaini)  g_free(args->Lgaini);
	if (args->RLgaini) g_free(args->RLgaini);
	if (args->LRgaini) g_free(args->LRgaini);
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
   /* in->modelMode = OBIT_SkyModel_DFT; Only can do DFT */

  /* Threading */
  if (in->threadArgs) {
    /* Check type - only handle "vmbeam" */
    args = (VMBeamFTFuncArg*)in->threadArgs[0];
    if ((strlen(args->type)>6) || (!strncmp(args->type, "vmbeam", 6))) {
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
 * Update model in Rgain, Lgain, RLgain, LRgain etc, members for comps member 
 * Sets begVMModelTime to current time
 * Sets endVMModelTime for when parallactic angle differs by 1 degree.
 * The parallel hand corrections should convert to a nominal value for radius.
 * The cross hand correction when multiplied by Stokes I will give the response.
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
  olong npos[2], lcomp, ncomp, i, iaty, ifield, lithread, iant=1, plane=0;
  ofloat **Rgain=NULL,  **Lgain=NULL,  **RLgain=NULL,  **LRgain=NULL, *ccData=NULL;
  ofloat **Rgaini=NULL, **Lgaini=NULL, **RLgaini=NULL, **LRgaini=NULL;
  ofloat curPA, tPA, tTime, bTime, fscale, PBCor, xx, yy, iPBCor, jPBCor;
  ofloat xr, xi, tr, ti, v, v1, v2;
  ofloat RXpol, LYpol, RXpolIm, LYpolIm;
  ofloat minPBCor=0.0, bmNorm, fblank = ObitMagicF();
  gboolean isCirc=TRUE, badBm=FALSE;
  odouble x, y;
  VMBeamFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMBeamUpdateModel";

  if (err->error) return;

  /* Thread to update */
  lithread = MAX (0, ithread);
  args = (VMBeamFTFuncArg*)in->threadArgs[lithread];
  /* g_assert (!strncmp(args->type,"vmbeam",6));  Test arg */

  /* Check subarray */
  Obit_return_if_fail(((suba>=0) && (suba<in->numAntList)), err,
		      "%s: bad subarray number %d", routine, suba+1);

  /* Check gain array dimension */
  Obit_return_if_fail((args->dimGain>=in->numComp), err,
    "%s:More components %d than allowed in gains  %d", 
    routine, in->numComp, args->dimGain);

  /* Find valid Antenna */
  iant = 1;
  for (i=0; i<in->AntList[suba]->number; i++) {
    if (in->AntList[suba]->ANlist[i]->AntID>0) {
      iant = in->AntList[suba]->ANlist[i]->AntID; break;
    }
  }

  /* Need Parallactic angle */
  tTime = time - 5.0*suba;  /* Correct time for subarray offset  from DBCON */
  bTime = tTime;
  curPA = ObitAntennaListParAng (in->AntList[suba], iant, tTime, in->curSource);

  /* Beginning of model validity */
  args->begVMModelTime = tTime + 5.0*suba;

  /* Find end time of validity - need Parallactic angle */
  tTime += 1.0/1440.0;
  tPA = ObitAntennaListParAng (in->AntList[suba], iant, tTime, in->curSource);
  /* Step by a min until the parallactic angle changes by 1 deg */
  while (fabs(tPA-curPA) < 1.0*DG2RAD) {
    tTime += 1.0/1440.0;
    tPA = ObitAntennaListParAng (in->AntList[suba], iant, tTime, in->curSource);
    /* But not forever */
    if (tTime-time>0.25) break;
  }

  /* Time for next update */
  args->endVMModelTime = tTime + 5.0*suba;

  /* Get parallactic angle half way between now and end */
  bTime = time + (tTime-bTime) / 2.0;
  curPA = ObitAntennaListParAng (in->AntList[suba], iant, bTime, in->curSource);
  /* Rotate cross pol by parallactic angle */
  xr  = cos (curPA);
  xi  = sin (curPA);
  curPA *= RAD2DG;  /* To deg */
  /* curPA += 180.0;   Not sure why */

  /* Local copies of pointers */
  Rgain  = args->Rgain; Lgain = args->Lgain; Rgaini = args->Rgaini; Lgaini = args->Lgaini;
  if (in->doCrossPol) {
    RLgain  = args->RLgain; LRgain = args->LRgain; RLgaini = args->RLgaini; LRgaini = args->LRgaini;
  } else {
    RLgain  = NULL; LRgain  = NULL; RLgaini = NULL; LRgaini = NULL;
  }
     /* Loop over antenna type */
  for (iaty=0; iaty<in->numAntType; iaty++) {
    npos[0] = 0; npos[1] = 0; 
    ccData = ObitFArrayIndex(in->comps, npos);
    lcomp = in->comps->naxis[0];  /* Length of row in comp table */
    ncomp = in->numComp;          /* number of components */
    
    /* Does this have circular polarization feeds? */
    isCirc =  uvdata->myDesc->crval[uvdata->myDesc->jlocs]==-1.0;
    
    /* Scale by ratio of frequency to beam image ref. frequency */
    plane = in->FreqPlane[MIN(args->channel, in->numUVChann-1)];
    fscale = args->BeamFreq / in->RXBeam[iaty]->freqs[plane];
    
    /* Beam normalization (center (power) defined to be 0.5 per parallel hand)*/
    v1 = ObitImageInterpValueInt (in->RXBeam[iaty], args->BeamRXInterp[iaty], 0.0, 0.0, curPA, plane, err);
    v2 = ObitImageInterpValueInt (in->LYBeam[iaty], args->BeamLYInterp[iaty], 0.0, 0.0, curPA, plane, err);
    badBm = (v1==fblank) || (v2==fblank);
    if (badBm) bmNorm = 0.5;
    else if (isCirc) bmNorm = 0.5 / (0.5*(v1 + v2)); /* Circular feedds */
    else             bmNorm = 1.0 / (v1 + v2);       /* Linear feeds */
    bmNorm=1.0; /* DEBUG */
    
    /* Compute antenna gains and put into Rgain, Lgain, RLgain, LRgain */
    for (i=0; i<ncomp; i++) {
      Lgain[iaty][i]  = 1.0;
      Rgain[iaty][i]  = 1.0;
      Lgaini[iaty][i] = 0.0;
      Rgaini[iaty][i] = 0.0;
      if (in->doCrossPol) {
	RLgain[iaty][i]  = 1.0;
	LRgain[iaty][i]  = 1.0;
	RLgaini[iaty][i] = 0.0;
	LRgaini[iaty][i] = 0.0;
      }
      /* Where in the beam? Offsets are from pointing position */
      xx = -ccData[i*lcomp+1];  /* AZ opposite of RA; offsets in beam images are of source */
      /* ?? xx = ccData[i*lcomp+1];  AZ opposite of RA; offsets in beam images are of source */
      yy =  ccData[i*lcomp+2];
      /* Scale by ratio of frequency to beam image ref. frequency */
      x = (odouble)xx * fscale;
      y = (odouble)yy * fscale;
      ifield = ccData[i*lcomp+0]+0.5;
      if (ifield<0) continue;
      
      /* Get symmetric primary (Power) beam correction for component */
      PBCor = getPBBeam(in->BeamShape, in->mosaic->images[ifield]->myDesc, xx, yy, in->antSize,  
			    args->BeamFreq, minPBCor); 
      iPBCor  = 1.0 / PBCor; 
      jPBCor = sqrtf(iPBCor);
       
      /* Interpolate gains - RR and LL (XX, YY) as voltage gains */
      v = ObitImageInterpValueInt (in->RXBeam[iaty], args->BeamRXInterp[iaty], x, y, curPA, plane, err);
      RXpol = bmNorm*v; RXpolIm = 0.0;
      if (badBm || (v==fblank) || (fabs(v)<0.001) || (fabs(v)>1.1)) RXpol = 1.0;
      if (in->doCmplx && in->RXBeamIm[iaty]) {
	v = ObitImageInterpValueInt (in->RXBeamIm[iaty], args->BeamRXPhInterp[iaty], x, y, curPA, plane, err);
	if (v!=fblank) RXpolIm = v;
      }
      v = ObitImageInterpValueInt (in->LYBeam[iaty], args->BeamLYInterp[iaty], x, y, curPA, plane, err);
      LYpol = bmNorm*v; LYpolIm = 0.0;
      if (badBm || (v==fblank) || (fabs(v)<0.001) || (fabs(v)>01.1)) LYpol = 1.0;
      if (in->doCmplx && in->LYBeamIm[iaty]) {
	v = ObitImageInterpValueInt (in->LYBeamIm[iaty], args->BeamLYPhInterp[iaty], x, y, curPA, plane, err);
	if (v!=fblank) LYpolIm = v;
      }
      /* Put voltage gain in arrays */
      if (in->doCmplx && in->RXBeamIm[iaty] && in->LYBeamIm[iaty]) {
	/* Circular feeds */
	if (isCirc) {
	  Rgain[iaty][i]  = jPBCor*RXpol;
	  Rgaini[iaty][i] = jPBCor*RXpolIm;
	  Lgain[iaty][i]  = jPBCor*LYpol;
	  Lgaini[iaty][i] = jPBCor*LYpolIm;
	} else {  /* linear feeds */
	  Rgain[iaty][i]  = jPBCor*RXpol;
	  Rgaini[iaty][i] = jPBCor*RXpolIm;
	  Lgain[iaty][i]  = jPBCor*LYpol;
	  Lgaini[iaty][i] = jPBCor*LYpolIm;
	} /* end linear feeds */
      } else { /* no phase, save power gains */
	/* Circular feeds */
	if (isCirc) {
	  /* Not sure about circ correction */
	  ti = 0.5*(iPBCor*(RXpol * RXpol) + iPBCor*(LYpol * LYpol));     /* Stokes I correction */
	  tr = 0.5*(iPBCor*(RXpol * RXpol) - iPBCor*(LYpol * LYpol));     /* Stokes V correction */
	  Rgain[iaty][i]  = ti + tr;
	  Lgain[iaty][i]  = ti - tr; 
	} else {  /* linear feeds */
	  Rgain[iaty][i]  = iPBCor*(RXpol * RXpol);
	  Lgain[iaty][i]  = iPBCor*(LYpol * LYpol);
	  Rgaini[iaty][i] = 0.0;
	  Lgaini[iaty][i] = 0.0;
	} /* end linear feeds */
      }
      if (fabs(PBCor)<0.01) {
	Rgain[iaty][i] = 1.0; Rgaini[iaty][i] = 0.0;
	Lgain[iaty][i] = 1.0; Lgaini[iaty][i] = 0.0;
      }
      /* If this is CLEANing use appropriate gains */
      if (in->doBeamCorClean && isCirc && (fabs(PBCor)>0.01)) {
	if (in->doCmplx && in->RXBeamIm[iaty] && in->LYBeamIm[iaty]) {
	  /* Using phase beam- really doesn't matter */
	  Rgain[iaty][i]  = iPBCor*(RXpol * RXpol);
	  Lgain[iaty][i]  = iPBCor*(LYpol * LYpol);
	  ti = 0.5*(Rgain[iaty][i] + Lgain[iaty][i]);  /* Stokes I correction */
	  tr = 0.5*(Rgain[iaty][i] - Lgain[iaty][i]);  /* Stokes V correction */
	  Rgain[iaty][i]  = ti + tr;
	  Lgain[iaty][i]  = ti - tr; 
	} else { /* no phase */
	  Rgain[iaty][i]  = iPBCor*(RXpol * RXpol);
	  Lgain[iaty][i]  = iPBCor*(LYpol * LYpol);
	  ti = 0.5*(Rgain[iaty][i] + Lgain[iaty][i]);     /* Stokes I correction */
	  tr = 0.5*(Rgain[iaty][i] - Lgain[iaty][i]);     /* Stokes V correction */
	  Rgain[iaty][i]  = ti + tr;
	  Lgain[iaty][i]  = ti - tr; 
	}
      } /* End force Stokes V for circular feeds */
      
      /* Cross pol corrections wanted? */
      if (in->doCrossPol) {
	/*  RL or XY */
	v =  ObitImageInterpValueInt (in->RLBeam[iaty], args->BeamRLInterp[iaty], x, y, curPA, plane, err);
	if (v==fblank) RLgain[iaty][i] = 0.0;
	else           RLgain[iaty][i] = jPBCor*v;
	if (in->doCmplx && in->RLBeamIm[iaty]) {
	  v = ObitImageInterpValueInt (in->RLBeamIm[iaty], args->BeamRLPhInterp[iaty], x, y, curPA, plane, err);
	  if (v==fblank) RLgaini[iaty][i] = 0.0;
	  else           RLgaini[iaty][i] = jPBCor*v;
	}
	if (isCirc) {
	  /* counter rotate paralactic angle */
	  tr = RLgain[iaty][i]; ti =  RLgaini[iaty][i];
	  RLgain[iaty][i]  =  tr*xr + xi*ti;
	  RLgaini[iaty][i] = -tr*xi + xr*ti;
	}
	/*  LR or YX */
	v = ObitImageInterpValueInt (in->LRBeam[iaty], args->BeamLRInterp[iaty], x, y, curPA, plane, err);
	if (LRgain[iaty][i]==fblank) LRgain[iaty][i] = 0.0;
	else                         LRgain[iaty][i] = jPBCor*v;
	if (in->doCmplx && in->LRBeamIm[iaty] ) {
	  v  = ObitImageInterpValueInt (in->LRBeamIm[iaty], args->BeamLRPhInterp[iaty], x, y, curPA, plane, err);
	  if (v==fblank) LRgaini[iaty][i] = 0.0;
	  else           LRgaini[iaty][i] = jPBCor*v;
	} 
	if (isCirc) {
	  /* counter rotate paralactic angle */
	  tr = LRgain[iaty][i]; ti = LRgaini[iaty][i];
	  LRgain[iaty][i]  =  tr*xr + xi*ti;
	  LRgaini[iaty][i] = -tr*xi + xr*ti;
	}
	/* Too low gain? */
	if (fabs(PBCor)<0.01) {
	  RLgain[iaty][i] = 1.0; RLgaini[iaty][i] = 0.0;
	  LRgain[iaty][i] = 1.0; LRgaini[iaty][i] = 0.0;
      }
     } /* end if crosspol */
    } /* end loop over components */
  } /* endloop over antenna type */

} /* end ObitSkyModelVMBeamUpdateModel */

/**
 * Determine antenna numbers
 * Only works for the 1st subarray (AN table 1)
 * \param in      SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamSetAnt (ObitSkyModelVMBeam* in, ObitUV *uvdata, 
			       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong iver;
  ObitTableAN *TableAN=NULL;
  ObitTableANRow *row;
  olong j;
  gchar *routine = "ObitSkyModelVMBeamSetAnt";

  if (err->error) return;

  /* Antenna Table */
  iver = 1;
  TableAN = newObitTableANValue ("AN", (ObitData*)uvdata, &iver, 
				   OBIT_IO_ReadOnly, 0, 0, 0, err);
  /* Create table row */
  row = newObitTableANRow (TableAN);

  /* Open table */
  retCode = ObitTableANOpen (TableAN, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,in->name);

  /* loop over table */
  while (retCode==OBIT_IO_OK) {
    retCode = ObitTableANReadRow (TableAN, -1, row, err);
    if (retCode == OBIT_IO_EOF) break;
    
    /* Check diameter */
    for (j=0; j<in->numAntType; j++) {
      if (row->diameter && (in->Diams[j]>0.0)) {
	if (abs(row->diameter-in->Diams[j])<0.001) in->AntType[row->noSta-1] = j;
      }
      else in->AntType[row->noSta-1] = 0; /* If no diameter given */
    }
    
    /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_msg (err, routine,in->name);
  } /* end of reading table */

    /* Close table */
retCode = ObitTableANClose (TableAN, err);
if ((retCode != OBIT_IO_OK) || (err->error))
  Obit_traceback_msg (err, routine, in->name);
/* Done - cleanup*/
row = ObitTableANRowUnref (row);
TableAN = ObitTableANUnref(TableAN);

} /* end ObitSkyModelVMBeamSetAnt */
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
 Initialize global ClassInfo Function pointers.
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
  in->Rgain        = NULL;
  in->Lgain        = NULL;
  in->RLgain       = NULL;
  in->LRgain       = NULL;
  in->Rgaini       = NULL;
  in->Lgaini       = NULL;
  in->RLgaini      = NULL;
  in->LRgaini      = NULL;
  in->BeamShape    = NULL;
  in->numAntType   = 0;
  in->AntType      = NULL;
  in->numPlane     = NULL;
  in->RXBeam       = NULL;
  in->RLBeam       = NULL;
  in->LRBeam       = NULL;
  in->LYBeam       = NULL;
  in->RXBeamIm     = NULL;
  in->RLBeamIm     = NULL;
  in->LRBeamIm     = NULL;
  in->LYBeamIm     = NULL;
  in->AntList      = NULL;
  in->curSource    = NULL;
  in->numAntList   = 0;
  in->Threshold    = 0.0;
  in->maxResid     = 0.0;
  in->doBeamCorClean = FALSE;
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
  for (i=0; i<in->numAntType; i++) {
    if ((in->Rgain!=NULL)  && (in->Rgain[i]!=NULL))   g_free(in->Rgain[i]);
    if ((in->Lgain!=NULL)  && (in->Lgain[i]!=NULL))   g_free(in->Lgain[i]);
    if ((in->LRgain!=NULL) && (in->RLgain[i]!=NULL))  g_free(in->RLgain[i]);
    if ((in->LRgain!=NULL) && (in->LRgain[i]!=NULL))  g_free(in->LRgain[i]);
    if ((in->Rgaini!=NULL) && (in->Rgaini[i]!=NULL))  g_free(in->Rgaini[i]);
    if ((in->Lgaini!=NULL) && (in->Lgaini[i]!=NULL))  g_free(in->Lgaini[i]);
    if ((in->RLgaini!=NULL)&& (in->RLgaini[i]!=NULL)) g_free(in->RLgaini[i]);
    if ((in->LRgaini!=NULL)&& (in->LRgaini[i]!=NULL)) g_free(in->LRgaini[i]);
  }
  if (in->Rgain)   g_free(in->Rgain);  in->Rgain   = NULL;
  if (in->Lgain)   g_free(in->Lgain);  in->Lgain   = NULL;
  if (in->RLgain)  g_free(in->RLgain); in->RLgain  = NULL;
  if (in->LRgain)  g_free(in->LRgain); in->LRgain  = NULL;
  if (in->Rgaini)  g_free(in->Rgaini); in->Rgaini  = NULL;
  if (in->Lgaini)  g_free(in->Lgaini); in->Lgaini  = NULL;
  if (in->RLgaini) g_free(in->RLgaini);in->RLgaini = NULL;
  if (in->LRgaini) g_free(in->LRgaini);in->LRgaini = NULL;
  if (in->numPlane)g_free(in->numPlane);in->numPlane = NULL;
  if (in->AntType) g_free(in->AntType);    in->AntType     = NULL;
  in->BeamShape = ObitBeamShapeUnref(in->BeamShape);
  for (i=0; i<in->numAntType; i++) {
    in->RXBeam[i]    = ObitImageInterpUnref(in->RXBeam[i]);
    in->RLBeam[i]    = ObitImageInterpUnref(in->RLBeam[i]);
    in->LRBeam[i]    = ObitImageInterpUnref(in->LRBeam[i]);
    in->LYBeam[i]    = ObitImageInterpUnref(in->LYBeam[i]);
    in->LYBeamIm[i]  = ObitImageInterpUnref(in->RXBeamIm[i]);
    in->RLBeamIm[i]  = ObitImageInterpUnref(in->RLBeamIm[i]);
    in->LRBeamIm[i]  = ObitImageInterpUnref(in->LRBeamIm[i]);
    in->LYBeamIm[i]  = ObitImageInterpUnref(in->LYBeamIm[i]);
  }
  g_free(in->RXBeam[i]);    in->RXBeam[i]    = NULL;
  g_free(in->RLBeam[i]);    in->RLBeam[i]    = NULL;
  g_free(in->LRBeam[i]);    in->LRBeam[i]    = NULL;
  g_free(in->LYBeam[i]);    in->LYBeam[i]    = NULL;
  g_free(in->RXBeamIm[i]);  in->LYBeamIm[i]  = NULL;
  g_free(in->RLBeamIm[i]);  in->RLBeamIm[i]  = NULL;
  g_free(in->LRBeamIm[i]);  in->LRBeamIm[i]  = NULL;
  g_free(in->LYBeamIm[i]);  in->LYBeamIm[i]  = NULL;
  
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
    args = (VMBeamFTFuncArg*)in->threadArgs[0];
    if ((strlen(args->type)>6) || (!strncmp(args->type, "vmbeam", 6))) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMBeamFTFuncArg*)in->threadArgs[i];
	for (i=0; i<in->numAntType; i++) {
	  if ((args->Rgain!=NULL)   && (args->Rgain[i]!=NULL))   g_free(args->Rgain[i]);
	  if ((args->Lgain!=NULL)   && (args->Lgain[i]!=NULL))   g_free(args->Lgain[i]);
	  if ((args->LRgain!=NULL)  && (args->RLgain[i]!=NULL))  g_free(args->RLgain[i]);
	  if ((args->LRgain!=NULL)  && (args->LRgain[i]!=NULL))  g_free(args->LRgain[i]);
	  if ((args->Rgaini!=NULL)  && (args->Rgaini[i]!=NULL))  g_free(args->Rgaini[i]);
	  if ((args->Lgaini!=NULL)  && (args->Lgaini[i]!=NULL))  g_free(args->Lgaini[i]);
	  if ((args->RLgaini!=NULL) && (args->RLgaini[i]!=NULL)) g_free(args->RLgaini[i]);
	  if ((args->LRgaini!=NULL) && (args->LRgaini[i]!=NULL)) g_free(args->LRgaini[i]);
	}
	if (args->Rgain)   g_free(args->Rgain);
	if (args->Lgain)   g_free(args->Lgain);
	if (args->RLgain)  g_free(args->RLgain);
	if (args->LRgain)  g_free(args->LRgain);
	if (args->Rgaini)  g_free(args->Rgaini);
	if (args->Lgaini)  g_free(args->Lgaini);
	if (args->RLgaini) g_free(args->RLgaini);
	if (args->LRgaini) g_free(args->LRgaini);
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
  ObitInfoListGetTest(in->info, "Threshold", &type, dim, &InfoReal);
  in->Threshold = InfoReal.flt;

  /* Is this part of a CLEAN? */
  in->doBeamCorClean = FALSE;
  ObitInfoListGetTest(in->info, "BeamCorClean", &type, dim,  &in->doBeamCorClean);

  /* Current maximum abs residual flux density */
  InfoReal.flt = -1.0; type = OBIT_float;
  lookup = !ObitInfoListGetTest(in->info, "maxResid", &type, (gint32*)dim, &InfoReal);
  if (lookup || (InfoReal.flt<0.0)) {
    /* Need to lookup from attached Mosaic */
    maxv = 0.0;
    for (i=0; i<in->mosaic->numberImages; i++) {
      /* Ignore tapered images */
      if (in->mosaic->BeamTaper[i]>0.0) continue;
      gotit = ObitInfoListGetTest (in->mosaic->images[i]->info, "maxAbsResid", 
				   &type, dim, &maxAbsResid); 
      if (gotit) maxv = MAX (maxv, maxAbsResid);
    } /* end loop over mosaic */

    /* If still nothing use max/min in image headers */
    if (maxv==0.0) {
      for (i=0; i<in->mosaic->numberImages; i++) {
	/* Ignore tapered images */
	if (in->mosaic->BeamTaper[i]>0.0) continue;
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
  olong i, j, mcomp, iComp, pos[2], nvis, lovis, hivis, nvisPerThread, nThreads;
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  VMBeamFTFuncArg *args;
  ofloat *ddata;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelVMBeamFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* Check */
  args = (VMBeamFTFuncArg*)in->threadArgs[0];
  if ((strlen(args->type)>6) || (strncmp(args->type, "vmbeam", 6))) {
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
      args->dimGain = in->numComp;
      for (j=0; j<in->numAntType; j++) {
	if ((args->Rgain!=NULL)   && (args->Rgain[j]!=NULL))   g_free(args->Rgain[j]);
	if ((args->Lgain!=NULL)   && (args->Lgain[j]!=NULL))   g_free(args->Lgain[j]);
	if ((args->LRgain!=NULL)  && (args->RLgain[j]!=NULL))  g_free(args->RLgain[j]);
	if ((args->LRgain!=NULL)  && (args->LRgain[j]!=NULL))  g_free(args->LRgain[j]);
	if ((args->Rgaini!=NULL)  && (args->Rgaini[j]!=NULL))  g_free(args->Rgaini[j]);
	if ((args->Lgaini!=NULL)  && (args->Lgaini[j]!=NULL))  g_free(args->Lgaini[j]);
	if ((args->RLgaini!=NULL) && (args->RLgaini[j]!=NULL)) g_free(args->RLgaini[j]);
	if ((args->LRgaini!=NULL) && (args->LRgaini[j]!=NULL)) g_free(args->LRgaini[j]);
	args->Rgain[j]  = g_malloc0(args->dimGain*sizeof(ofloat));
	args->Rgaini[j] = g_malloc0(args->dimGain*sizeof(ofloat));
	args->Lgain[j]  = g_malloc0(args->dimGain*sizeof(ofloat));
	args->Lgaini[j] = g_malloc0(args->dimGain*sizeof(ofloat));
	if (in->doCrossPol) {
	  args->RLgain[j]  = g_malloc0(args->dimGain*sizeof(ofloat));
	  args->RLgaini[j] = g_malloc0(args->dimGain*sizeof(ofloat));
	  args->LRgain[j]  = g_malloc0(args->dimGain*sizeof(ofloat));
	  args->LRgaini[j] = g_malloc0(args->dimGain*sizeof(ofloat));
	}
      } /* end antenna type */
    } /* end rebuild arrays */
    /* Update which vis */
    lovis += nvisPerThread;
    hivis += nvisPerThread;
    hivis = MIN (hivis, nvis);
  } /* end loop over threads */

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
 * \li Rgain   Float array of time/spatially variable R/X component gain
 * \li Lgain   Float array of time/spatially variable L/Y component gain
 * \li RLgain  Float array of time/spatially variable RL/XY component gain
 * \li LRgain  Float array of time/spatially variable LR/YX component gain
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
  ofloat **Rgain    = largs->Rgain;
  ofloat **Lgain    = largs->Lgain;
  ofloat **RLgain   = largs->RLgain;
  ofloat **LRgain   = largs->LRgain;
  ofloat **RLgaini  = largs->RLgaini;
  ofloat **LRgaini  = largs->LRgaini;

  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, iaty, jaty;
  olong lrec, nrparm, naxis[2], channel, plane;
  olong jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF, lim;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, suba, it1, it2, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale, exparg;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealRL, sumImagRL, modRealRL=0.0, modImagRL=0.0;
  ofloat sumRealLR, sumImagLR, modRealLR=0.0, modImagLR=0.0;
  ofloat **rgain1, **lgain1, tx, ty, tz, ll, lll;
  ofloat **rlgain1, **lrgain1, re, im, specFact;
  ofloat **rlgain1i, **lrgain1i;
  ofloat amp, ampr, ampl, arg, freq2=0.0, freqFact, wtRR=0.0, wtLL=0.0, temp;
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArrR[FazArrSize], AmpArrL[FazArrSize];
  ofloat AmpArr[FazArrSize],  FazArr[FazArrSize];
  ofloat CosArr[FazArrSize],  SinArr[FazArrSize];
  olong it, jt, kt, itcnt, nextCh, ch1; /*dbgcnt=0; DEBUG */
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

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF       = in->startIFPB-1;
  numberIF      = MAX (1, in->numberIFPB);
  jincif        = uvdata->myDesc->incif;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  nUVchan       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocf];
  jincf         = uvdata->myDesc->incf;
  nUVpoln       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocs];
  jincs         = uvdata->myDesc->incs;  /* increment in real array */
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocif>=0) {
    if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
      kincf = 1;
      kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
    } else { /* IF beforefreq  */
      kincif = 1;
      kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
    }
  } else {  /* NO IF axis */
      kincif = 1;
      kincf  = 1;
  }

  /* Cross or only parallel pol? */
  doCrossPol = (nUVpoln > 2) && in->doCrossPol;
  /* Only parallel for divide */
  if (in->doDivide) doCrossPol = FALSE;

  /* Get pointer for components */
  naxis[0] = 0; naxis[1] = 0; 
  Data = ObitFArrayIndex(in->comps, naxis);
  lcomp = in->comps->naxis[0];   /* Length of row in comp table */
  mcomp = in->numComp;           /* Actual number */

  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

  /* Current channel (0-rel) */
  channel = nextCh = 0;
  plane   = in->FreqPlane[largs->channel];  /* Which plane in correction cube */

  /* Outer loop over blocks of channels */
  /* Starting parameters this pass */
  lstartIF       = startIF;
  lstartChannel  = startChannel;
  while (nextCh<(numberIF*numberChannel)) {
    
    /* Loop over vis in buffer */
    lrec    = uvdata->myDesc->lrec;         /* Length of record */
    visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
    nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */
    
    for (iVis=loVis; iVis<hiVis; iVis++) {
      
      /* Is current model still valid? */
      if ((visData[iloct] > largs->endVMModelTime) || (visData[iloct] < largs->begVMModelTime)) {
	/* Current channel - which plane in correction to apply? */
	largs->channel  = channel;
	largs->BeamFreq = freqArr[channel];
	/* Subarray 0-rel */
	ObitUVDescGetAnts(uvdata->myDesc, visData, &it1, &it2, &suba);
	/* Update */
 	myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba-1, uvdata, ithread, err);
	if (err->error) {
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
			 routine);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}
      }
      
      /* Need antennas numbers */
      ObitUVDescGetAnts(uvdata->myDesc, visData, &ant1, &ant2, &it1);
      ant1--;    /* 0 rel */
      ant2--;    /* 0 rel */
      iaty = in->AntType[ant1];  /* Antenna type */
      jaty = in->AntType[ant2];
     
      /* Loop over IFs */
      ch1 = lstartChannel;
      for (iIF=lstartIF; iIF<startIF+numberIF; iIF++) {
	offsetIF = nrparm + iIF*jincif; 
	if (iIF>lstartIF) ch1 = 0;  /* beginning channel after first IF */
	for (iChannel=ch1; iChannel<startChannel+numberChannel; iChannel++) {
	  channel = iIF* nUVchan + iChannel; /* UV Channel */
	  offsetChannel = offsetIF + iChannel*jincf; 
	  freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */
	  freq2    = freqFact*freqFact;    /* Frequency factor squared */
	  
	  /* Sum over components */
	  /* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	  sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	  sumRealRL = sumImagRL = sumRealLR = sumImagLR = 0.0;
	  /* Set component gain lists by antenna and type */
	  ddata   = Data;
	  rgain1  = Rgain;
	  lgain1  = Lgain;
	  rlgain1 = RLgain;
	  lrgain1 = LRgain;
	  /* Imaginary parts - needed for X pol */
	  rlgain1i = RLgaini;
	  lrgain1i = LRgaini;
	  
	  /* Sum by model type - assume phase same for RR, LL */
	  kt = 0;
	  switch (in->modType) {
	  case OBIT_SkyModel_PointMod:     /* Point */
	    /* From the AIPSish QXXPTS.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArrR[itcnt] = ddata[3] * rgain1[iaty][iComp];
		AmpArrL[itcnt] = ddata[3] * lgain1[jaty][iComp];
		AmpArr[itcnt]  = ddata[3];
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_PointModSpec:     /* Point + spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[7+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  FazArr[itcnt]  = freqFact * (tx + ty + tz);
		  AmpArrR[itcnt] = specFact * ddata[3] * rgain1[iaty][iComp];
		  AmpArrL[itcnt] = specFact * ddata[3] * lgain1[jaty][iComp];
		  AmpArr[itcnt]  = specFact * ddata[3];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
		ddata += lcomp;  /* update pointer */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	    /* From the AIPSish QGASUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
				ddata[8]*visData[ilocv]*visData[ilocv] +
				ddata[9]*visData[ilocu]*visData[ilocv]);
		if (arg<-1.0e-5) exparg = exp (arg);
		else exparg = 1.0;
		ampr = ddata[3] * rgain1[iaty][iComp]*exparg;
		ampl = ddata[3] * lgain1[jaty][iComp]*exparg;
		AmpArrR[itcnt] = ampr;
		AmpArrL[itcnt] = ampl;
		AmpArr[itcnt]  = ddata[3]*exparg;
		itcnt++;          /* Count in amp/phase buffers */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[10+iterm] * lll;
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
		  AmpArrR[itcnt] = amp * rgain1[iaty][iComp];
		  AmpArrL[itcnt] = amp * lgain1[jaty][iComp];
		  AmpArr[itcnt]  = amp;
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	    /* From the AIPSish QSPSUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ddata[7];
		arg = MAX (arg, 0.1);
		arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		AmpArrR[itcnt] = ddata[3] * rgain1[iaty][iComp] * arg;
		AmpArrL[itcnt] = ddata[3] * lgain1[jaty][iComp] * arg;
		AmpArr[itcnt]  = ddata[3] * arg;
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereModSpec:    /* Uniform sphere + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
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
		  AmpArrR[itcnt] = amp * rgain1[iaty][iComp];
		  AmpArrL[itcnt] = amp * lgain1[jaty][iComp];
		  AmpArr[itcnt]  = amp;
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  default:
	    ObitThreadLock(in->thread);  /* Lock against other threads */
	    Obit_log_error(err, OBIT_Error,"%s Unknown Comp model type %d in %s",
			   routine, in->modType, in->name);
	    ObitThreadUnlock(in->thread); 
	    goto finish;
	  }; /* end switch by model type */
	  
	  /* Apply model */
	  modRealRR = sumRealRR;
	  modImagRR = sumImagRR;
	  modRealLL = sumRealLL;
	  modImagLL = sumImagLL;
	  if (doCrossPol) {
	    modRealRL = sumRealRL;
	    modImagRL = sumImagRL;
	    modRealLR = sumRealLR;
	    modImagLR = sumImagLR;
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
	      if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
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
	      if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealLL;
	      visData[offset+1] -= modImagLL;
	    }
	  } /* end LL not blanked */
	  
	  if (doCrossPol) {
	    /* RL */
	    iStoke = 2;
	    offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	    
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
                visData[offset]   = modRealRL;
                visData[offset+1] = modImagRL;
		if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	      } else {
		/* Subtract model */
                visData[offset]   -= modRealRL;
                visData[offset+1] -= modImagRL;
	      }
	    } /* end RL not blanked */
	    
	    /* LR  */
	    offset += jincs;
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
		visData[offset]   = modRealLR;
                visData[offset+1] = modImagLR;
		if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	      } else {
		/* Subtract model */
		visData[offset]   -= modRealRL;
                visData[offset+1] -= modImagRL;
	      }
	    } /* end LR not blanked */
	  } /* end crosspol */
	  
	  offsetChannel += jincf;
	  nextCh = channel+1;  /* Next channel */
	  /* Have we finished this plane in the correction cubes? 
	   If so break frequency looping but continue until visibility loop done.
	   Then, restart channel/IF loop ing in new plane */
	  if (plane!=in->FreqPlane[MIN(nextCh, (in->numUVChann-1))]) {
	    /* Reset gains & channels if this the last vis in loop */
	    if (iVis>=(hiVis-1)) {
	      plane   = in->FreqPlane[MIN(nextCh, (in->numUVChann-1))];  /* Which plane in correction cube */
	      largs->endVMModelTime = -1.0e20;  
	      lstartChannel = iChannel;
	      lstartIF      = iIF;
	      if (iChannel==(startChannel+numberChannel-1)) { /* end of channel loop */
	    	lstartChannel = 0;  /* start next IF */
		lstartIF++;
	      }
	    } /* end if last vis */
	    /* if ((iVis==(hiVis-1)) && (dbgcnt<20)) {
	      fprintf (stderr,"End  plane=%d ch=%d %d %d\n",plane,nextCh,lstartChannel,lstartIF);  DEBUG 
	      dbgcnt++;}*/
	    if (channel>1) goto newPlane;  /* Break looping except for first channel */
	  } /* end if new plane */
	  /* if ((dbgcnt<50) && (iVis==(hiVis-1)))
	    fprintf (stderr,"   ch=%d IF %d channel %d\n",iChannel,iIF,nextCh); DEBUG */
	} /* end loop over Channel */
	offsetIF += jincif;
      } /* end loop over IF */
      
      /* Come to here if finished with correction plane */
    newPlane:
      visData += lrec; /* Update vis pointer */
    } /* end loop over visibilities */
    /* if (dbgcnt<50)
      fprintf (stderr," New chann=%d ch=%d IF=%d\n",channel,lstartChannel,lstartIF); DEBUG
      dbgcnt++; */
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
 * \li Rgain   Float array of time/spatially variable R/X component gain
 * \li Rgaini  Float array of time/spatially variable R/X imaginary component gain
 * \li Lgain   Float array of time/spatially variable L/Y component gain
 * \li Lgaini  Float array of time/spatially variable L/Y imaginary component gain
 * \li RLgain  Float array of time/spatially variable RL/XY component gain
 * \li RLgaini Float array of time/spatially variable RL/XY imaginary component gain
 * \li LRgain  Float array of time/spatially variable LR/YX component gain
 * \li LRgaini Float array of time/spatially variable LR/YX imaginary component gain
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
  ofloat **Rgain    = largs->Rgain;
  ofloat **Lgain    = largs->Lgain;
  ofloat **RLgain    = largs->RLgain;
  ofloat **LRgain    = largs->LRgain;
  ofloat **Rgaini   = largs->Rgaini;
  ofloat **Lgaini   = largs->Lgaini;
  ofloat **RLgaini   = largs->RLgaini;
  ofloat **LRgaini   = largs->LRgaini;

  olong iVis, iIF, iChannel, iStoke, iComp, lcomp;
  olong lrec, nrparm, naxis[2], channel, plane;
  olong jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF, lim, iaty, jaty;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, suba, it1, it2, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealRL,  sumImagRL,  modRealRL=0.0,  modImagRL=0.0;
  ofloat sumRealLR,  sumImagLR,  modRealLR=0.0,  modImagLR=0.0;
  ofloat tx, ty, tz, ll, lll;
  ofloat re, im, specFact, ModR, ModI;
  ofloat **rgain1, **lgain1, **rgain1i, **lgain1i;
  ofloat **rlgain1, **lrgain1, **rlgain1i, **lrgain1i;
  ofloat amp, ampr, ampl, arg, freq2=0.0, freqFact, wtRR=0.0, wtLL=0.0, temp;
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArrRr[FazArrSize], AmpArrLr[FazArrSize], AmpArrRi[FazArrSize], AmpArrLi[FazArrSize];
  ofloat Arr1Rr[FazArrSize], Arr1Ri[FazArrSize], Arr1Lr[FazArrSize], Arr1Li[FazArrSize];
  ofloat Arr2Rr[FazArrSize], Arr2Ri[FazArrSize], Arr2Lr[FazArrSize], Arr2Li[FazArrSize];
  ofloat FazArr[FazArrSize], AmpArr[FazArrSize];
  ofloat CosArr[FazArrSize], SinArr[FazArrSize];
  olong it, jt, kt, itcnt, nextCh, ch1;
  gboolean doCrossPol;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMBeamFTDFTPh";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF       = in->startIFPB-1;
  numberIF      = MAX (1, in->numberIFPB);
  jincif        = uvdata->myDesc->incif;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  nUVchan       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocf];
  jincf         = uvdata->myDesc->incf;
  nUVpoln       = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocs];
  jincs         = uvdata->myDesc->incs;  /* increment in real array */
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocif>=0) {
    if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
      kincf = 1;
      kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
    } else { /* IF beforefreq  */
      kincif = 1;
      kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
    }
  } else {  /* NO IF axis */
    kincif = 1;
    kincf  = 1;
  }

  /* Cross or only parallel pol? */
  doCrossPol = (nUVpoln > 2) && in->doCrossPol;
  /* Only parallel for divide */
  if (in->doDivide) doCrossPol = FALSE;

  /* Get pointer for components */
  naxis[0] = 0; naxis[1] = 0; 
  Data = ObitFArrayIndex(in->comps, naxis);
  lcomp = in->comps->naxis[0];   /* Length of row in comp table */
  mcomp = in->numComp;           /* Actual number */

  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;

  /* Current channel (0-rel) */
  channel = nextCh = 0;
  plane   = in->FreqPlane[largs->channel];  /* Which plane in correction cube */

  /* Outer loop over blocks of channels */
  /* Starting parameters this pass */
  lstartIF       = startIF;
  lstartChannel  = startChannel;
  while (nextCh<(numberIF*numberChannel)) {
    
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
	ObitUVDescGetAnts(uvdata->myDesc, visData, &it1, &it2, &suba);
	/* Update */
	myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba-1, uvdata, ithread, err);
	if (err->error) {
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
			 routine);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}
      }
      
      /* Need antennas numbers */
      ObitUVDescGetAnts(uvdata->myDesc, visData, &ant1, &ant2, &it1);
      ant1--;    /* 0 rel */
      ant2--;    /* 0 rel */
      iaty = in->AntType[ant1];  /* Antenna type */
      jaty = in->AntType[ant2];
     
      /* Loop over IFs */
      channel = lstartIF* nUVchan + lstartChannel; /* UV Channel */
      ch1 = lstartChannel;
      for (iIF=lstartIF; iIF<startIF+numberIF; iIF++) {
	if (iIF>lstartIF) ch1 = 0;  /* beginning channel after first IF */
	offsetIF = nrparm + iIF*jincif; 
	for (iChannel=ch1; iChannel<startChannel+numberChannel; iChannel++) {
	  offsetChannel = offsetIF + iChannel*jincf; 
	  freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */
	  freq2    = freqFact * freqFact;   	           /* Frequency factor squared */

	  /* Sum over components */
	  /* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	  sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	  sumRealRL  = sumImagRL  = sumRealLR  = sumImagLR  = 0.0;
	  /* Set component gain lists by antenna and type */
	  ddata   = Data;
	  rgain1  = Rgain;
	  lgain1  = Lgain;
	  rlgain1 = RLgain;
	  lrgain1 = LRgain;
	  /* Imaginary parts */
	  rgain1i = Rgaini;
	  lgain1i = Lgaini;
	  rlgain1i = RLgaini;
	  lrgain1i = LRgaini;

	  /* Sum by model type - assume phase same for RR, LL */
	  kt = 0;
	  switch (in->modType) {
	  case OBIT_SkyModel_PointMod:     /* Point */
	    /* From the AIPSish QXXPTS.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArr[itcnt]   = ddata[3];
		AmpArrRr[itcnt] = ddata[3] * rgain1[iaty][iComp];
		AmpArrLr[itcnt] = ddata[3] * lgain1[jaty][iComp];
		AmpArrRi[itcnt] = ddata[3] * rgain1i[iaty][iComp];
		AmpArrLi[itcnt] = ddata[3] * lgain1i[jaty][iComp];
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_PointModSpec:     /* Point + spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
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
		  amp = ddata[3] * specFact;
		  FazArr[itcnt]   = freqFact * (tx + ty + tz);
		  AmpArr[itcnt]   = amp;
		  Arr1Rr[itcnt] = rgain1[iaty][iComp];
		  Arr1Lr[itcnt] = lgain1[iaty][iComp];
		  Arr1Ri[itcnt] = rgain1i[iaty][iComp];
		  Arr1Li[itcnt] = lgain1i[iaty][iComp];
		  Arr2Rr[itcnt] = rgain1[jaty][iComp];
		  Arr2Lr[itcnt] = lgain1[jaty][iComp];
		  Arr2Ri[itcnt] = rgain1i[jaty][iComp];
		  Arr2Li[itcnt] = lgain1i[jaty][iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
		ddata += lcomp;  /* update pointer */
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);

	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		/* Parallel  pol */
		re = AmpArr[jt]*CosArr[jt];
		im = AmpArr[jt]*SinArr[jt];
		ModR = Arr1Rr[jt]*Arr2Rr[jt] + Arr1Ri[jt]*Arr2Ri[jt];
		ModI = Arr1Ri[jt]*Arr2Rr[jt] - Arr1Rr[jt]*Arr2Ri[jt];
		sumRealRR += ModR*re - ModI*im;
		sumImagRR += ModR*im + ModI*re;

		ModR = Arr1Lr[jt]*Arr2Lr[jt] + Arr1Li[jt]*Arr2Li[jt];
		ModI = Arr1Li[jt]*Arr2Lr[jt] - Arr1Lr[jt]*Arr2Li[jt];
		sumRealLL += ModR*re - ModI*im;
		sumImagLL += ModR*im + ModI*re;
		if (doCrossPol) {
		  /* Cross pol */
		  ModR = Arr1Rr[jt]*Arr2Lr[jt] + Arr1Ri[jt]*Arr2Li[jt];
		  ModI = Arr1Ri[jt]*Arr2Lr[jt] - Arr1Rr[jt]*Arr2Li[jt];
		  sumRealRL += ModR*re - ModI*im;
		  sumImagRL += ModR*im + ModI*re;

		  ModR = Arr1Lr[jt]*Arr2Rr[jt] + Arr1Li[jt]*Arr2Ri[jt];
		  ModI = Arr1Li[jt]*Arr2Rr[jt] - Arr1Lr[jt]*Arr2Ri[jt];
		  sumRealRL += ModR*re - ModI*im;
		  sumImagRL += ModR*im + ModI*re;
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	    /* From the AIPSish QGASUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			       ddata[8]*visData[ilocv]*visData[ilocv] +
			       ddata[9]*visData[ilocu]*visData[ilocv]);
		if (arg<-1.0e-5) amp = ddata[3] * exp (arg);
		else             amp = ddata[3];
		ampl = amp * lgain1i[iaty][iComp];
		ampr = amp * rgain1i[jaty][iComp];
		AmpArrRi[itcnt] = ampr;
		AmpArrLi[itcnt] = ampl;

		ampl = amp * lgain1i[iaty][iComp];
		ampr = amp * rgain1i[jaty][iComp];
		AmpArrRi[itcnt] = ampr;
		AmpArrLi[itcnt] = ampl;
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
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
		  FazArr[itcnt]   = freqFact * (tx + ty + tz);
		  AmpArr[itcnt]   = ddata[3];
		  AmpArrRr[itcnt] = ddata[3] * rgain1[iaty][iComp];
		  AmpArrLr[itcnt] = ddata[3] * lgain1[jaty][iComp];
		  AmpArrRi[itcnt] = ddata[3] * rgain1i[iaty][iComp];
		  AmpArrLi[itcnt] = ddata[3] * lgain1i[jaty][iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	    /* From the AIPSish QSPSUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ddata[7];
		arg = MAX (arg, 0.1);
		arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		amp = ddata[3] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		AmpArr[itcnt]   = amp;
		AmpArrRr[itcnt] = amp * rgain1[iaty][iComp];
		AmpArrLr[itcnt] = amp * lgain1[jaty][iComp];
		AmpArrRi[itcnt] = amp * rgain1i[iaty][iComp];
		AmpArrLi[itcnt] = amp * lgain1i[jaty][iComp];
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereModSpec:    /* Uniform sphere + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
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
		  FazArr[itcnt]   = freqFact * (tx + ty + tz);
		  AmpArr[itcnt]   = amp;
		  AmpArrRr[itcnt] = amp * rgain1[iaty][iComp];
		  AmpArrLr[itcnt] = amp * lgain1[jaty][iComp];
		  AmpArrRi[itcnt] = amp * rgain1i[iaty][iComp];
		  AmpArrLi[itcnt] = amp * lgain1i[jaty][iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
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
		  re = AmpArr[jt]*CosArr[jt];
		  im = AmpArr[jt]*SinArr[jt];
		  sumRealRL += re * rlgain1[iaty][kt+jt]  - im * rlgain1i[jaty][kt+jt];
		  sumImagRL += re * rlgain1i[iaty][kt+jt] + im * rlgain1[jaty][kt+jt];
		  sumRealLR += re * lrgain1[jaty][kt+jt]  - im * lrgain1i[iaty][kt+jt];
		  sumImagLR += re * lrgain1i[jaty][kt+jt] + im * lrgain1[iaty][kt+jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	      kt = it+1;  /* offset in rlgain/lrgain */
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
	    modRealRL = sumRealRL;
	    modImagRL = sumImagRL;
	    modRealLR = sumRealLR;
	    modImagLR = sumImagLR;
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
	      if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
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
	      if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealLL;
	      visData[offset+1] -= modImagLL;
	    }
	  } /* end LL not blanked */
	  
	  if (doCrossPol) {
	    /* RL */
	    iStoke = 2;
	    offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	    
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
                visData[offset]   = modRealRL;
                visData[offset+1] = modImagRL;
		if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	      } else {
		/* Subtract model */
                visData[offset]   -= modRealLR;
		visData[offset+1] -= modImagLR;
	      }
	    } /* end RL not blanked */
	    
	    /* LR  */
	    offset += jincs;
	    /* Ignore blanked data unless replacing the data */
	    if ((visData[offset+2]>0.0) || in->doReplace) {
	      /* Apply model to data */
	      if (in->doReplace) {  /* replace data with model */
		visData[offset]   = modRealLR;
                visData[offset+1] = modImagLR;
		if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	      } else {
		/* Subtract model */
		visData[offset]   -= modRealLR;
		visData[offset+1] -= modImagLR;
	      }
	    } /* end LR not blanked */
	  } /* end crosspol */
	  
	  offsetChannel += jincf;
	  nextCh = channel+1;  /* Next channel */
	  /* Have we finished this plane in the correction cubes? 
	   If so break frequency looping but continue until visibility loop done.
	   Then, restart channel/IF loop ing in new plane */
	  if (plane!=in->FreqPlane[MIN(nextCh, (in->numUVChann-1))]) {
	    /* Reset gains & channels if this the last vis in loop */
	    if (iVis>=(hiVis-1)) {
	      plane   = in->FreqPlane[MIN(nextCh, (in->numUVChann-1))];  /* Which plane in correction cube */
	      largs->endVMModelTime = -1.0e20;  
	      lstartChannel = iChannel;
	      lstartIF      = iIF;
	      if (iChannel==(startChannel+numberChannel-1)) { /* end of channel loop */
	    	lstartChannel = 0;  /* start next IF */
		lstartIF++;
	      }
	    } /* end if last vis */
	    if (channel>1) goto newPlane;  /* Break looping except for first channel */
	  } /* end if new channel */
	} /* end loop over Channel */
	offsetIF += jincif;
      } /* end loop over IF */
      
      /* Come to here is finished with correction plane */
    newPlane:
      visData += lrec; /* Update vis pointer */
    } /* end loop over visibilities */
  } /* end outer frequency loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadSkyModelVMBeamFTDFTPh */

/** 
 *  Get model frequency primary beam 
 * \param beamSize Image BeamSize object
 * \param dec      Image Descriptor
 * \param x        x offset from pointing center (deg)
 * \param y        y offset from pointing center (deg)
 * \param antSize  Antenna diameter in meters. (defaults to 25.0)
 * \param freq     Frequency in Hz
 * \param pbmin    Minimum antenna gain Jinc 0=>0.05, poly 0=> 0.01
 * \return Relative gain at freq refFreq wrt average of Freq.
*/
static ofloat getPBBeam(ObitBeamShape *beamShape, ObitImageDesc *desc, 
			ofloat x, ofloat y, ofloat antSize, odouble freq, ofloat pbmin)
{
  ofloat PBBeam = 1.0;
  odouble Angle;
  /*odouble RAPnt, DecPnt, ra, dec, xx, yy, zz;*/

  if (antSize <=0.0) antSize = 25.0;  /* default antenna size */
  /* If antSize != 13.5, it's not MeerKAT 
  if (fabs(antSize-13.5)>0.5) beamShape->doMeerKAT = FALSE;*/
  /* DEBUG antSize=13.5; */
  ObitBeamShapeSetAntSize(beamShape, antSize);
  if (freq!=beamShape->refFreq) ObitBeamShapeSetFreq(beamShape, freq);

  /* Get pointing position
  ObitImageDescGetPoint(desc, &RAPnt, &DecPnt); */

  /* Convert offset to position 
  ObitSkyGeomXYShift (RAPnt, DecPnt, x, y, ObitImageDescRotate(desc), &ra, &dec);*/

  /* Angle from center
  RAPnt  *= DG2RAD;
  DecPnt *= DG2RAD;
  xx = DG2RAD * ra;
  yy = DG2RAD * dec;
  zz = sin (yy) * sin (DecPnt) + cos (yy) * cos (DecPnt) * cos (xx-RAPnt);
  zz = MIN (zz, 1.000);
  Angle = acos (zz) * RAD2DG; */
  Angle = sqrt (x*x + y*y);

  /* Get gain */
  PBBeam = ObitBeamShapeGainSym(beamShape, Angle);
 
  return PBBeam;
} /* end  getPBBeam */

