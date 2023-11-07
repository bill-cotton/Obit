/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2023                                          */
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
#include "ObitMatx.h"
#include "ObitExp.h"
#include "ObitComplex.h"

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

/** Private: Threaded full model FTDFT */
static gpointer ThreadSkyModelVMBeamFTDFT (gpointer arg);

/** Private: get model frequency primary beam */
static ofloat getPBBeam(ObitBeamShape *beamShape, ObitImageDesc *desc, 
			ofloat x, ofloat y, 
			ofloat antSize, odouble freq, ofloat pbmin);

/** Private: Fill Antenna Jones Matrix */
static void FillJones(ObitMatx *Jones, odouble x, odouble y, ofloat curPA, olong plane,
		      ofloat bmNorm, ofloat iPBCor, ofloat jPBCor, gboolean isCirc,
		      ObitMatx *JPA,
		      ObitImageInterp* P1BeamRe,  ObitFInterpolate* BeamP1ReInterp,
		      ObitImageInterp* P1BeamIm,  ObitFInterpolate* BeamP1ImInterp,
		      ObitImageInterp* P2BeamRe,  ObitFInterpolate* BeamP2ReInterp,
		      ObitImageInterp* P2BeamIm,  ObitFInterpolate* BeamP2InInterp,
		      ObitImageInterp* X1BeamRe,  ObitFInterpolate* BeamX1ReInterp,
		      ObitImageInterp* X1BeamIm,  ObitFInterpolate* BeamX1ImInterp,
		      ObitImageInterp* X2BeamRe,  ObitFInterpolate* BeamX2ReInterp,
		      ObitImageInterp* X2BeamIm,  ObitFInterpolate* BeamX2ImInterp,
		      ObitMatx *work, ObitErr *err);

/*---------------Private structures----------------*/
/* FT threaded function argument 
 Note: Derived classes MUST have the following entries at the beginning 
 of the corresponding structure */
typedef struct {
  /* type "vmbeam" in this class */
  gchar type[12];
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitSkyModelVMBeam *in;
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
  /** Dimension of JMatrix  */
  olong dimGain;
  /** Arrays of time/spatially variable Jones matrices per type, per component */
  ObitMatx ***JMatrix;
  /** Dimension of workVis  (256) */
  olong dimWork;
  /** Work matrix and array */
  ObitMatx **workVis, *work1, *work2, *sumVis;
  /** cos and sin of twice parallactic angle */
  ofloat cos2PA, sin2PA;
} VMBeamFTFuncArg;

/** Private: Make model args structure */
void ObitSkyModelVMBeamMakeArg (ObitSkyModelVMBeam *in, ObitUV *uvdata, ObitErr *err);

/** Private: Update model args structure */
void ObitSkyModelVMBeamUpdateArg (ObitSkyModelVMBeam *in, VMBeamFTFuncArg *args);

/** Private: Kill model arg structure */
gpointer ObitSkyModelVMBeamKillArg (ObitSkyModelVMBeam *in);

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
  out->RXBeamIm   = g_malloc0(numAntType*sizeof(ObitImageInterp*));
  out->LYBeamIm   = g_malloc0(numAntType*sizeof(ObitImageInterp*));
  out->RLBeamIm   = g_malloc0(numAntType*sizeof(ObitImageInterp*));
  out->LRBeamIm   = g_malloc0(numAntType*sizeof(ObitImageInterp*));
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

    /* Imaginary beams */
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
    if (doCmplx && out->RLBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->RLBeam[i]->ImgPixels, 
						   out->RLBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable pq amp, phase beam arrays", routine);
    }
    if (doCmplx && out->LRBeamIm[i]) {
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
    if (RXBeamIm[i]  && (RXBeamIm[i]->image)) RXBeamIm[i]->image = ObitImageUnref(RXBeamIm[i]->image);
    if (RLBeamIm[i]  && (RLBeamIm[i]->image)) RLBeamIm[i]->image = ObitImageUnref(RLBeamIm[i]->image);
    if (LRBeamIm[i]  && (LRBeamIm[i]->image)) LRBeamIm[i]->image = ObitImageUnref(LRBeamIm[i]->image);
    if (LYBeamIm[i]  && (LYBeamIm[i]->image)) LYBeamIm[i]->image = ObitImageUnref(LYBeamIm[i]->image);
  }
  /* Set antenna Types */
  ObitSkyModelVMBeamSetAnt (out, uvData, err);

  return out;
} /* end ObitSkyModelVMBeamCreate */

/**
 * Initializes Sky Model
 * Checks that data contain parallel hands, save calibration/selection request
 * and set uv data for no selection/calibration
 * \param in      SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err     Obit error stack object.
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
  strncpy (in->saveStokes, "   ", 4);
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
  if (in->threadArgs==NULL) ObitSkyModelVMBeamMakeArg(in, uvdata, err);

  /* Call parent initializer */
  ObitSkyModelVMInitMod(inn, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fourier transform routines - DFT only */
  /* Now only full Jones corrections */
  in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamFTDFT;

  /* Check requested Stokes
  Obit_return_if_fail((!strncmp(in->stokes,"    ",4)), err,
		      "%s: Unsupported Stokes %s", routine, in->stokes); */

  /* Check that data contains RR, LL or (XX,YY)*/
  uvDesc = uvdata->myDesc;
  Obit_return_if_fail((((uvDesc->crval[uvDesc->jlocs]==-1.0) || 
			(uvDesc->crval[uvDesc->jlocs]==-5.0)) && 
		       (uvDesc->inaxes[uvDesc->jlocs]>=2)), err,
		      "%s: RR,LL or XX,YY not in UV data", routine);

  /* Antenna Lists */
  /* Get TableList */
  list = uvdata->tableList;
  numAntList = ObitTableListGetHigh (list, "AIPS AN");  /* How many subarrays? */
  if (numAntList!=in->numAntList) { /* Rebuild Antenna Lists if needed */
    for (iver=1; iver<=in->numAntList; iver++) { 
      in->AntList[iver-1] = ObitAntennaListUnref(in->AntList[iver-1]);
    }
    if (in->AntList) {g_free(in->AntList);} in->AntList = NULL;
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

  /* Call parent shutdown */
  ObitSkyModelVMShutDownMod(inn, uvdata, err);
 
  /* Delete args if the correct type */
  if (in->threadArgs) in->threadArgs = ObitSkyModelVMBeamKillArg (in);
  ObitThreadPoolFree (in->thread);  /* Shut down any threading */
 
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
 * Update model in JMatrix, members for comps member 
 * Sets begVMModelTime to current time
 * Sets endVMModelTime for when parallactic angle differs by 0.3 degree.
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
  ofloat *ccData=NULL;
  ObitMatx*** JMatrix=NULL, *JPA=NULL;
  ofloat curPA, tPA, tTime, bTime, fscale, PBCor, xx, yy, iPBCor, jPBCor;
  ofloat v1, v2, minPBCor=0.0, bmNorm, fblank = ObitMagicF();
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
  /* Step by a min until the parallactic angle changes by 0.25 deg */
  while (fabs(tPA-curPA) < 0.25*DG2RAD) {
    tTime += 1.0/1440.0;
    tPA = ObitAntennaListParAng (in->AntList[suba], iant, tTime, in->curSource);
    /* But not forever */
    if (tTime-time>0.25) break;
  }

  /* Time for next update */
  args->endVMModelTime = tTime + 5.0*suba;

  /* Get parallactic angle half way between now and end */
  bTime = time + (tTime-bTime) / 2.0;
  bTime = time; /* at current time */
  curPA = ObitAntennaListParAng (in->AntList[suba], iant, bTime, in->curSource);
  /* save cos, sin of twice parallactic angle */
  args->cos2PA = cos(2.0*curPA);
  args->sin2PA = sin(2.0*curPA);
  curPA *= RAD2DG;  /* To deg */
  /* rotation matrix for linear feeds */
  JPA = args->work1;
  ObitMatxSet2C (JPA, args->cos2PA,0.0, args->sin2PA,0.0, -args->sin2PA,0.0, args->cos2PA,0.0);
 
  /* Local copies of pointers */
  JMatrix = args->JMatrix;
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
    
    /* Beam normalization (center (power) defined to be 0.5 per parallel hand - NO really 1.0 */
    v1 = ObitImageInterpValueInt (in->RXBeam[iaty], args->BeamRXInterp[iaty], 0.0, 0.0, curPA, plane, err);
    v2 = ObitImageInterpValueInt (in->LYBeam[iaty], args->BeamLYInterp[iaty], 0.0, 0.0, curPA, plane, err);
    badBm = (v1==fblank) || (v2==fblank);
    if (badBm) bmNorm = 0.5;
    else if (isCirc) bmNorm = 0.5 / (0.5*(v1 + v2)); /* Circular feedds */
    /* NO else             bmNorm = 1.0 / (v1 + v2);        Linear feeds */
    else             bmNorm = 2.0 / (v1 + v2);       /* Linear feeds */
    bmNorm=1.0; /* DEBUG */
    
    /* Compute antenna gains and put into JMatrix */
    for (i=0; i<ncomp; i++) {
      ObitMatxUnit(JMatrix[iaty][i]); /* Initialize with unit matrix */
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

      /* interpolate to Jones Matrix */
      if ((!badBm) && (fabs(PBCor)>0.01))
	/* DEBUG FillJones(JMatrix[iaty][i], x, y, curPA, plane, bmNorm, iPBCor, jPBCor, isCirc, args->s2X, args->c2X,*/
	FillJones(JMatrix[iaty][i], x, -y, -curPA, plane, bmNorm, iPBCor, jPBCor, isCirc, JPA,
		  in->RXBeam[iaty], args->BeamRXInterp[iaty], in->RXBeamIm[iaty], args->BeamRXPhInterp[iaty],
		  in->LYBeam[iaty], args->BeamLYInterp[iaty], in->RXBeamIm[iaty], args->BeamRXPhInterp[iaty],
		  in->RLBeam[iaty], args->BeamRLInterp[iaty], in->RLBeamIm[iaty], args->BeamRLPhInterp[iaty],
		  in->LRBeam[iaty], args->BeamLRInterp[iaty], in->LRBeamIm[iaty], args->BeamLRPhInterp[iaty], 
		  args->work2, err); 
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

/** Public: Jones correct and sum SkyModel components 
 * converts SkyModel comps to pseudovisibilities, multiple by Jones matrices and sum.
 * \param numComp  Number of SkyModel components
 * \param Stokes   Stokes type, 1,2,3,4 => I, Q, U, V
 * \param isCirc   True if visibility basis is circular, else linear
 * \param compFlux Component amplitude with any factors multiplied
 * \param sinArr   Array of sine of phase to apply, per numComp
 * \param cosArr   Array of cosine of phase to apply, per numComp
 * \param Jones1   Jones matrix array of first antenna of baseline, per numComp
 * \param Jones2   Jones matrix array of second antenna of baseline, per numComp
 * \param workVis  work array of numComp 2x2 complex matrices
 * \param work1    work 2x2 complex matrix
 * \param work2    work 2x2 complex matrix
 * \param sumVis   [out] 2x2 complex matrix for resultant sum
 */
void ObitSkyModelVMBeamJonesCorSum (olong numComp, olong Stokes,  gboolean isCirc, 
				    ofloat *compFlux, ofloat *sinArr, ofloat *cosArr,
				    ObitMatx **Jones1, ObitMatx **Jones2,
				    ObitMatx **workVis, ObitMatx *work1, ObitMatx *work2,
				    ObitMatx *sumVis)
{
  olong i;

  /* Create pseudovis -
     By feed type */
  if (isCirc) {  /* For circular feeds */
    /* By Stokes Type */
    switch (Stokes) {
    case 1: /* Circular Stokes I */
      /* vis: I+jV, Q+jU, Q-jU, I-V */
      for (i=0; i<numComp; i++) {
	ObitMatxSet2C (workVis[i],
		       compFlux[i]*cosArr[i],compFlux[i]*sinArr[i],
		       0.0,0.0, 0.0,0.0,
		       compFlux[i]*cosArr[i],compFlux[i]*sinArr[i]);
      } /* end loop over components */
      break;
    case 2: /* Circular Stokes Q */
      for (i=0; i<numComp; i++) {
 	ObitMatxSet2C (workVis[i], 0.0,0.0,
		       compFlux[i]*cosArr[i],compFlux[i]*sinArr[i],
		       compFlux[i]*cosArr[i],compFlux[i]*sinArr[i], 0.0,0.0);
      } /* end loop over components */
      break;
    case 3: /* Circular Stokes U */
      for (i=0; i<numComp; i++) {
 	ObitMatxSet2C (workVis[i], 0.0,0.0,
		       -compFlux[i]*sinArr[i],+compFlux[i]*cosArr[i],
		       +compFlux[i]*sinArr[i],-compFlux[i]*cosArr[i], 0.0,0.0);
      } /* end loop over components */
      break;
    case 4: /* Circular Stokes V */
      for (i=0; i<numComp; i++) {
  	ObitMatxSet2C (workVis[i],
		       +compFlux[i]*cosArr[i],+compFlux[i]*sinArr[i],
		       0.0,0.0, 0.0,0.0,
		       -compFlux[i]*cosArr[i],-compFlux[i]*sinArr[i]);
     } /* end loop over components */
      break;
    default:  /* Should NEVER get here */
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch element type */
  } else {  /* For linear feeds */
    /* vis : I+Q, U+jV, U-jV, I-Q */
    /* By Stokes Type */
    switch (Stokes) {
    case 1: /* Linear Stokes I */
      for (i=0; i<numComp; i++) {
	ObitMatxSet2C (workVis[i], compFlux[i]*cosArr[i],compFlux[i]*sinArr[i], 0.0,0.0, 0.0,0.0, 
		                   compFlux[i]*cosArr[i],compFlux[i]*sinArr[i]);
      } /* end loop over components */
      break;
    case 2: /* Linear Stokes Q */
      for (i=0; i<numComp; i++) {
 	ObitMatxSet2C (workVis[i],
		       +compFlux[i]*cosArr[i],+compFlux[i]*sinArr[i], 0.0,0.0, 0.0,0.0, 
		       -compFlux[i]*cosArr[i],-compFlux[i]*sinArr[i]);
      } /* end loop over components */
      break;
    case 3: /* Linear Stokes U */
      for (i=0; i<numComp; i++) {
 	ObitMatxSet2C (workVis[i], 0.0,0.0, 
		       +compFlux[i]*cosArr[i],+compFlux[i]*sinArr[i],
		       +compFlux[i]*cosArr[i],+compFlux[i]*sinArr[i], 0.0,0.0);
      } /* end loop over components */
      break;
    case 4: /* Linear Stokes V */
      for (i=0; i<numComp; i++) {
 	ObitMatxSet2C (workVis[i], 0.0,0.0,
		       -compFlux[i]*sinArr[i],+compFlux[i]*cosArr[i],
		       +compFlux[i]*sinArr[i],-compFlux[i]*cosArr[i], 0.0,0.0);
     } /* end loop over components */
      break;
    default:  /* Should NEVER get here */
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch element type */
  } /* end of circular/linear */

  /* Loop multiplying each Jones1*workVis*conjugate transpose(Jones2) and sum */
  ObitMatxSet2C(sumVis, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0); /* Init */
  for (i=0; i<numComp; i++) {
    ObitMatxMult(Jones1[i], workVis[i], work1);
    ObitMatxMultCT(work1, Jones2[i], work2);
    ObitMatxAdd (sumVis, work2, sumVis);
  } /* end loop over components */

  
} /* end ObitSkyModelVMBeamJonesCorSum */

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
  ObitSkyModelVMBeam *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->numPlane){g_free(in->numPlane);} in->numPlane = NULL;
  if (in->AntType) {g_free(in->AntType);}  in->AntType  = NULL;
  in->BeamShape = ObitBeamShapeUnref(in->BeamShape);
  for (i=0; i<in->numAntType; i++) {
    if (in->RXBeam && in->RXBeam[i]) in->RXBeam[i]  = ObitImageInterpUnref(in->RXBeam[i]);
    if (in->RLBeam && in->RLBeam[i]) in->RLBeam[i]  = ObitImageInterpUnref(in->RLBeam[i]);
    if (in->LRBeam && in->LRBeam[i]) in->LRBeam[i]  = ObitImageInterpUnref(in->LRBeam[i]);
    if (in->LYBeam && in->LYBeam[i]) in->LYBeam[i]  = ObitImageInterpUnref(in->LYBeam[i]);
    if (in->RXBeamIm && in->RXBeamIm[i]) in->RXBeamIm[i] = ObitImageInterpUnref(in->RXBeamIm[i]);
    if (in->RLBeamIm && in->RLBeamIm[i]) in->RLBeamIm[i] = ObitImageInterpUnref(in->RLBeamIm[i]);
    if (in->LRBeamIm && in->LRBeamIm[i]) in->LRBeamIm[i] = ObitImageInterpUnref(in->LRBeamIm[i]);
    if (in->LYBeamIm && in->LYBeamIm[i]) in->LYBeamIm[i] = ObitImageInterpUnref(in->LYBeamIm[i]);
  }
  if (in->RXBeam)   {g_free(in->RXBeam);    in->RXBeam    = NULL;}
  if (in->RLBeam)   {g_free(in->RLBeam);    in->RLBeam    = NULL;}
  if (in->LRBeam)   {g_free(in->LRBeam);    in->LRBeam    = NULL;}
  if (in->LYBeam)   {g_free(in->LYBeam);    in->LYBeam    = NULL;}
  if (in->RXBeamIm) {g_free(in->RXBeamIm);  in->RXBeamIm  = NULL;}
  if (in->RLBeamIm) {g_free(in->RLBeamIm);  in->RLBeamIm  = NULL;}
  if (in->LRBeamIm) {g_free(in->LRBeamIm);  in->LRBeamIm  = NULL;}
  if (in->LYBeamIm) {g_free(in->LYBeamIm);  in->LYBeamIm  = NULL;}
  
  in->curSource = ObitSourceUnref(in->curSource);
  if (in->AntList)  {
    for (i=0; i<in->numAntList; i++) { 
      in->AntList[i] = ObitAntennaListUnref(in->AntList[i]);
    }
    g_free(in->AntList); in->AntList = NULL;
  }
    
  /* Thread stuff */
  if (in->threadArgs) ObitSkyModelVMBeamKillArg(in);

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
 * this routine will divide the buffer up amoung the number of processors
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
    args->in     = in;
    args->field  = field;
    args->uvdata = uvdata;
    args->first  = lovis;
    args->last   = hivis;
    if (nThreads>1) args->ithread= i;
    else args->ithread = -1;
    args->err    = err;
    /* Update which vis */
    lovis += nvisPerThread;
    hivis += nvisPerThread;
    hivis = MIN (hivis, nvis);
    ObitSkyModelVMBeamUpdateArg (in, args);  /* Redo arrays if needed */
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
 * \li dimGain Dimension of JMatrix
 * \li JMatrix Complex 2x2 Jones matrices per type/component
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
  ObitMatx ***JMatrix = largs->JMatrix;
  ObitMatx **workVis  = largs->workVis;
  ObitMatx *work1     = largs->work1;
  ObitMatx *work2     = largs->work2;
  ObitMatx *sumVis    = largs->sumVis;
  
  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, iaty, jaty;
  olong lrec, nrparm, naxis[2], channel, plane;
  olong jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF, lim;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, suba, it1, it2, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale, exparg, ll, lll, specFact;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealRL, sumImagRL, modRealRL=0.0, modImagRL=0.0;
  ofloat sumRealLR, sumImagLR, modRealLR=0.0, modImagLR=0.0;
  ofloat amp, arg, freq2=0.0, freqFact, wtRR=0.0, wtLL=0.0, temp;
#define FazArrSize 256  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArr[FazArrSize],  FazArr[FazArrSize];
  ofloat CosArr[FazArrSize],  SinArr[FazArrSize];
  olong it, kt, itcnt, nextCh, ch1, Stokes; /*dbgcnt=0; DEBUG */
  gboolean doCrossPol, isCirc;
  odouble *freqArr, tx, ty, tz;
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
  nUVpoln       = uvdata->myDesc->inaxes[uvdata->myDesc->jlocs];
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
  
  isCirc = uvdata->myDesc->crval[uvdata->myDesc->jlocs]==-1.0;  /* Circular or linear feeds? */
  Stokes = (olong)(in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocs]+0.5);
  if ((Stokes<1) || (Stokes>4)) Stokes = 1;  /* Just to be sure, IQUV =>1,2,3,4 */
  
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
	  ddata = Data;
	  
	  /* Sum over components */
	  /* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	  sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	  sumRealRL = sumImagRL = sumRealLR = sumImagLR = 0.0;

	  /* Sum by model type - assume phase same for RR, LL */
	  kt = 0;
	  switch (in->modType) {
	  case OBIT_SkyModel_PointMod:     /* Point */
	    /* From the AIPSish QXXPTS.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		if (ddata[3]!=0.0) {  /* valid? */
		  FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					      ddata[5]*visData[ilocv] + 
					      ddata[6]*visData[ilocw]);
		  /* Parallel  pol */
		  AmpArr[itcnt]  = ddata[3];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
		ddata += lcomp;   /* update pointer */
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);

	      /* Accumulate real and imaginary parts to sumVis */
	      ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					    AmpArr, SinArr, CosArr, &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					    workVis, work1, work2, sumVis);
	      
	      sumRealRR += sumVis->flt[0]; sumImagRR += sumVis->flt[1];
	      sumRealLL += sumVis->flt[6]; sumImagLL += sumVis->flt[7];
	      if (doCrossPol) {
		/* Cross pol */
		sumRealRL += sumVis->flt[2]; sumImagRL += sumVis->flt[3];
		sumRealLR += sumVis->flt[4]; sumImagLR += sumVis->flt[5];
	      } /* end xpol */
	      kt = it+1;  /* offset in JMatrix */
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
		  AmpArr[itcnt]  = specFact * ddata[3];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
		ddata += lcomp;  /* update pointer */
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);

	      /* Accumulate real and imaginary parts to sumVis */
	      ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					    AmpArr, SinArr, CosArr, &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					    workVis, work1, work2, sumVis);
	      
	      sumRealRR += sumVis->flt[0]; sumImagRR += sumVis->flt[1];
	      sumRealLL += sumVis->flt[6]; sumImagLL += sumVis->flt[7];
	      if (doCrossPol) {
		/* Cross pol */
		sumRealRL += sumVis->flt[2]; sumImagRL += sumVis->flt[3];
		sumRealLR += sumVis->flt[4]; sumImagLR += sumVis->flt[5];
	      } /* end xpol */
	      kt = it+1;  /* offset in JMatrix */
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
		AmpArr[itcnt]  = ddata[3]*exparg;
		itcnt++;          /* Count in amp/phase buffers */
	      }  /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts to sumVis */
	      ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					    AmpArr, SinArr, CosArr, &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					    workVis, work1, work2, sumVis);
	      
	      sumRealRR += sumVis->flt[0]; sumImagRR += sumVis->flt[1];
	      sumRealLL += sumVis->flt[6]; sumImagLL += sumVis->flt[7];
	      if (doCrossPol) {
		/* Cross pol */
		sumRealRL += sumVis->flt[2]; sumImagRL += sumVis->flt[3];
		sumRealLR += sumVis->flt[4]; sumImagLR += sumVis->flt[5];
	      } /* end xpol */
	      kt = it+1;  /* offset in JMatrix */
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
		  AmpArr[itcnt]  = amp;
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
	      }  /* end inner loop over components */

	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts to sumVis */
	      ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					    AmpArr, SinArr, CosArr, &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					    workVis, work1, work2, sumVis);
	      
	      sumRealRR += sumVis->flt[0]; sumImagRR += sumVis->flt[1];
	      sumRealLL += sumVis->flt[6]; sumImagLL += sumVis->flt[7];
	      if (doCrossPol) {
		/* Cross pol */
		sumRealRL += sumVis->flt[2]; sumImagRL += sumVis->flt[3];
		sumRealLR += sumVis->flt[4]; sumImagLR += sumVis->flt[5];
	      } /* end xpol */
	      kt = it+1;  /* offset in JMatrix */
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
		AmpArr[itcnt]  = ddata[3] * arg;
		ddata += lcomp;   /* update pointer */
		itcnt++;          /* Count in amp/phase buffers */
	      }
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts to sumVis */
	      ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					    AmpArr, SinArr, CosArr, &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					    workVis, work1, work2, sumVis);
	      
	      sumRealRR += sumVis->flt[0]; sumImagRR += sumVis->flt[1];
	      sumRealLL += sumVis->flt[6]; sumImagLL += sumVis->flt[7];
	      if (doCrossPol) {
		/* Cross pol */
		sumRealRL += sumVis->flt[2]; sumImagRL += sumVis->flt[3];
		sumRealLR += sumVis->flt[4]; sumImagLR += sumVis->flt[5];
	      } /* end xpol */
	      kt = it+1;  /* offset in JMatrix */
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
		  AmpArr[itcnt]  = amp;
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		ddata += lcomp;  /* update pointer */
	      }  /* end inner loop over components */

	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      /* Accumulate real and imaginary parts to sumVis */
	      ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					    AmpArr, SinArr, CosArr, &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					    workVis, work1, work2, sumVis);
	      
	      sumRealRR += sumVis->flt[0]; sumImagRR += sumVis->flt[1];
	      sumRealLL += sumVis->flt[6]; sumImagLL += sumVis->flt[7];
	      if (doCrossPol) {
		/* Cross pol */
		sumRealRL += sumVis->flt[2]; sumImagRL += sumVis->flt[3];
		sumRealLR += sumVis->flt[4]; sumImagLR += sumVis->flt[5];
	      } /* end xpol */
	      kt = it+1;  /* offset in JMatrix */
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

/** Private: Fill Antenna Jones Matrix 
 * Interpolate requested beam value given the interpolator
 * Divides by the nominal beam gain
 * Beams with NULL pointers are given 0 value
 * \param Jones   Jones matrix to load
 * \param bmNorm eam normalization factor 
 * \param x       X or Right Ascension desired (deg)
 * \param y       Y or Declination deg) desired 
 * \param curPA   Current Parallactic Angle to rotate  (deg)
 * \param plane   Image plane 0-rel (frequency) to use
 * \param iPBCor  1/PB Correction
 * \param jPBCor  sqrt(1/PB Correction)
 * \param isCirc  Are the feeds circular?
 * \param JPA     Rotation matrix for Parallactic angle for linear feeds
 * \param P1BeamRe       P1 (RR or XX) Real beam to interpolate
 * \param BeamP1ReInterp P1 (RR or XX) Real Interpolator
 * \param P1BeamIm       P1 (RR or XX) Imag beam to interpolate
 * \param BeamP1ImInterp P1 (RR or XX) Imag Interpolator
 * \param P2BeamRe       P2 (LL or YY) Real beam to interpolate
 * \param BeamP2ReInterp P2 (LL or YY) Real Interpolator
 * \param P2BeamIm       P2 (LL or YY) Imag beam to interpolate
 * \param BeamP2ImInterp P2 (LL or YY) Imag Interpolator
 * \param X1BeamRe       X1 (RL or XY) Real beam to interpolate
 * \param BeamX1ImInterp X1 (RL or XY) Real Interpolator
 * \param X1BeamIm       X1 (RL or XY) Imag beam to interpolate
 * \param BeamX1ReInterp X1 (RL or XY) Imag Interpolator
 * \param X2BeamRe       X2 (LR or YX) Real beam to interpolate
 * \param BeamX2ReInterp X2 (LR or YX) Real Interpolator
 * \param X2BeamRe       X2 (LR or YX) Imag beam to interpolate
 * \param BeamX2ReInterp X2 (LR or YX) Imag Interpolator
 * \param err     Obit error stack object.
 * \return interpolated beam value -  may be fblank
 */
static void FillJones(ObitMatx *Jones, odouble x, odouble y, ofloat curPA, olong plane,
		      ofloat bmNorm, ofloat iPBCor, ofloat jPBCor, gboolean isCirc, ObitMatx *JPA,
		      ObitImageInterp* P1BeamRe,  ObitFInterpolate* BeamP1ReInterp,
		      ObitImageInterp* P1BeamIm,  ObitFInterpolate* BeamP1ImInterp,
		      ObitImageInterp* P2BeamRe,  ObitFInterpolate* BeamP2ReInterp,
		      ObitImageInterp* P2BeamIm,  ObitFInterpolate* BeamP2ImInterp,
		      ObitImageInterp* X1BeamRe,  ObitFInterpolate* BeamX1ReInterp,
		      ObitImageInterp* X1BeamIm,  ObitFInterpolate* BeamX1ImInterp,
		      ObitImageInterp* X2BeamRe,  ObitFInterpolate* BeamX2ReInterp,
		      ObitImageInterp* X2BeamIm,  ObitFInterpolate* BeamX2ImInterp,
		      ObitMatx *work, ObitErr *err)
{
  /* parallel and cross (on and off) diagonal terms. */
  ofloat v, P1r=0., P1i=0., P2r=0., P2i=0., X1r=0., X1i=0., X2r=0., X2i=0.; 
  ofloat fblank = ObitMagicF();
  olong i;

  if (err->error) return;  /* going wrong? */
  
  /* Interpolate voltage gains */
  /* P1 (RR or XX) */
  v = bmNorm*ObitImageInterpValueInt (P1BeamRe, BeamP1ReInterp, x, y, curPA, plane, err);
  if (v!=fblank) P1r = jPBCor*bmNorm*v;
  else           P1r = 0.0;
  if (P1BeamIm) {
    v = bmNorm*ObitImageInterpValueInt (P1BeamIm, BeamP1ImInterp, x, y, curPA, plane, err);
    if (v!=fblank) P1i = jPBCor*bmNorm*v;
  } else           P1i = 0.0;

  /* P2 (LL or YY) */
  v = bmNorm*ObitImageInterpValueInt (P2BeamRe, BeamP2ReInterp, x, y, curPA, plane, err);
  if (v!=fblank) P2r = jPBCor*bmNorm*v;
  else           P2r = 0.0;
  if (P2BeamIm) {
    v = bmNorm*ObitImageInterpValueInt (P2BeamIm, BeamP2ImInterp, x, y, curPA, plane, err);
    if (v!=fblank) P2i = jPBCor*bmNorm*v;
  } else           P2i = 0.0;

  /* X1 (RL or XY) */
  if (X1BeamRe) {
    v = bmNorm*ObitImageInterpValueInt (X1BeamRe, BeamX1ReInterp, x, y, curPA, plane, err);
    if (v!=fblank) X1r = jPBCor*bmNorm*v;
  }  else          X1r = 0.0;
  if (X1BeamIm) {
    v = bmNorm*ObitImageInterpValueInt (X1BeamIm, BeamX1ImInterp, x, y, curPA, plane, err);
    if (v!=fblank) P1i = jPBCor*bmNorm*v;
  } else           P1i = 0.0;

  /* X2 (LR or YX) */
  if (X1BeamRe) {
    v = bmNorm*ObitImageInterpValueInt (X2BeamRe, BeamX2ReInterp, x, y, curPA, plane, err);
    if (v!=fblank) X2r = jPBCor*bmNorm*v;
  } else           X2r = 0.0;
  if (P1BeamIm) {
    v = bmNorm*ObitImageInterpValueInt (X2BeamIm, BeamX2ImInterp, x, y, curPA, plane, err);
    if (v!=fblank) X2i = jPBCor*bmNorm*v;
  } else           X2i = 0.0;

  /* Load 'em up */
  /*DEBUG ObitMatxSet2C (Jones, P1r, P1i, X1r, X1i, X2r, X2i, P2r, P2i);*/
  /* Need to correct by parallactic angle *****************************/
  /* Largely using formalism of Smirnov 2011 */
  ObitMatxSet2C (Jones, P1r, -P1i, X1r, -X1i, X2r, -X2i, P2r, -P2i);  /* Conjugate =transmit */
  /*ObitMatxSet2C (Jones, P2r, -P2i, X2r, -X2i, X1r, -X1i, P1r, -P1i);  swap & Conjugate */
  /*Nope ObitMatxSet2C (work, P1r, P1i, X1r, X1i, X2r, X2i, P2r, P2i);
    ObitMatxInv2x2 (work, Jones);*/
  /* Multiply by parallactic angle term */
  if (isCirc) {
    /* At least for VLA data (AIPS/Obit), parallactic angle corrections were made to the data */
  } else {  /* Linear */
    ObitMatxMult(Jones, JPA, work);
    /* Hack, copy to Jones */
    for (i=0; i<8; i++) Jones->array[i] = work->array[i];
  }
} /* end FillJones */

/** Private: Make model args structure */
void ObitSkyModelVMBeamMakeArg (ObitSkyModelVMBeam *in, ObitUV *uvdata, ObitErr *err)
{
  olong i, j, k, dim[2]={2,2};
  VMBeamFTFuncArg* args;
  gchar *routine="VMBeamMakeArg";
  
  if (!in->threadArgs) in->threadArgs = g_malloc0(in->nThreads*sizeof(VMBeamFTFuncArg*));
  for (i=0; i<in->nThreads; i++) {
    if (!in->threadArgs[i]) in->threadArgs[i] = g_malloc0(sizeof(VMBeamFTFuncArg)); 
    args = (VMBeamFTFuncArg*)in->threadArgs[i];
    strcpy (args->type, "vmbeam");  /* Enter type as first entry */
    args->in     = in;
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
      if (in->RXBeamIm[j]) args->BeamRXPhInterp[j] = ObitImageInterpCloneInterp(in->RXBeamIm[j],err);
      if (in->LYBeamIm[j]) args->BeamLYPhInterp[j] = ObitImageInterpCloneInterp(in->LYBeamIm[j],err);
      if (in->RLBeamIm[j]) args->BeamRLPhInterp[j] = ObitImageInterpCloneInterp(in->RLBeamIm[j],err);
      if (in->LRBeamIm[j]) args->BeamLRPhInterp[j] = ObitImageInterpCloneInterp(in->LRBeamIm[j],err);
    } /* end antenna type loop */
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    args->begVMModelTime = -1.0e20;
    args->endVMModelTime = -1.0e20;
    args->VMComps = NULL;
    args->cos2PA  = 1.0;
    args->sin2PA  = 0.0;
    args->dimGain = 0;
    args->JMatrix = g_malloc0(in->numAntType*sizeof(ObitMatx***));
    args->dimGain = -1;   /* Nothing yet */
    args->dimWork = 256;  /* Same size as arrays in FTDFT functions (FazArrSize) */
    args->workVis = g_malloc0(args->dimWork*sizeof(ObitMatx**));
    for (k=0; k<args->dimWork; k++) args->workVis[k] = ObitMatxCreate(OBIT_Complex, 2, dim);
    args->work1   = ObitMatxCreate(OBIT_Complex, 2, dim);
    args->work2   = ObitMatxCreate(OBIT_Complex, 2, dim);
    args->sumVis  = ObitMatxCreate(OBIT_Complex, 2, dim);
  } /* end loop over threads */
} /* end  ObitSkyModelVMBeamMakeArg */

/** Private: Update model args structure if the number of components has changed
 * \param in     Pointer to the ObitSkyModel .
 * \param args   argument list to update
 */
void ObitSkyModelVMBeamUpdateArg (ObitSkyModelVMBeam *in, VMBeamFTFuncArg *args)
{
  olong j, k, dim[2]={2,2};
  
  /* Build JMatrix as needed */
  if (args->dimGain!=in->numComp) {
    /* Delete old JMatrix if it exists */
    if (args->dimGain>0) {
      for (j=0; j<in->numAntType; j++) {
	if ((args->JMatrix!=NULL)   && (args->JMatrix[j]!=NULL)) {
	  for (k=0; k<args->dimGain; k++) args->JMatrix[j][k] = ObitMatxUnref(args->JMatrix[j][k]);
	  g_free(args->JMatrix[j]);
	}
      } /* end antenna type */
    } /* end if delete old */
    /* (Re)build  */
    for (j=0; j<in->numAntType; j++) {
      args->JMatrix[j] = g_malloc0(in->numComp*sizeof(ObitMatx**));
      for (k=0; k<in->numComp; k++) args->JMatrix[j][k] = ObitMatxCreate(OBIT_Complex, 2, dim);
    } /* end antenna type */
  } /* end rebuild arrays */
  args->dimGain = in->numComp;
} /* end ObitSkyModelVMBeamUpdateArg */

/** Private: Kill model arg structure */
gpointer ObitSkyModelVMBeamKillArg (ObitSkyModelVMBeam *in)
{
  olong i, j, k;
  VMBeamFTFuncArg* args;
  
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
	if (args->JMatrix && args->JMatrix[j]) {
	  for (k=0; k<args->dimGain; k++) args->JMatrix[j][k] = ObitMatxUnref(args->JMatrix[j][k]);
	  g_free(args->JMatrix[j]); args->JMatrix[j] = NULL;
	}
      } /* end ant type loop */
      g_free(args->JMatrix); 
      if (args->workVis) {
	for (k=0; k<args->dimWork; k++) args->workVis[k] = ObitMatxUnref(args->workVis[k]);
	g_free(args->workVis); args->workVis = NULL;
      }
      g_free(args->work1);ObitMatxUnref(args->work1);
      g_free(args->work2);ObitMatxUnref(args->work2);
      if (args->BeamRXInterp)   {g_free(args->BeamRXInterp);}   args->BeamRXInterp = NULL;
      if (args->BeamLYInterp)   {g_free(args->BeamLYInterp);}   args->BeamLYInterp = NULL;
      if (args->BeamRLInterp)   {g_free(args->BeamRLInterp);}   args->BeamRLInterp = NULL;
      if (args->BeamLRInterp)   {g_free(args->BeamLRInterp);}   args->BeamLRInterp = NULL;
      if (args->BeamRXPhInterp) {g_free(args->BeamRXPhInterp);} args->BeamRXPhInterp = NULL;
      if (args->BeamLYPhInterp) {g_free(args->BeamLYPhInterp);} args->BeamLYPhInterp = NULL;
      if (args->BeamRLPhInterp) {g_free(args->BeamRLPhInterp);} args->BeamRLPhInterp = NULL;
      if (args->BeamLRPhInterp) {g_free(args->BeamLRPhInterp);} args->BeamLRPhInterp = NULL;
      args->work1  = ObitMatxUnref(args->work1); 
      args->work2  = ObitMatxUnref(args->work2); 
      args->sumVis = ObitMatxUnref(args->sumVis); 
      g_free(in->threadArgs[i]);
    } /* Loop over threads */
    g_free(in->threadArgs);
    in->threadArgs = NULL;
  } /* end if this a "vmbeam" threadArg */
  return in->threadArgs;
} /* end ObitSkyModelVMBeamKillArg */
