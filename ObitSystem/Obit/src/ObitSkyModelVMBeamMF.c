/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011-2026                                          */
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

#include "ObitAntenna.h"
#include "ObitTableANUtil.h"
#include "ObitUVDesc.h"
#include "ObitThread.h"
#if HAVE_GPU==1  /* GPU? */
#include "ObitGPUSkyModel.h"
#include "ObitGPUBeamModel.h"
#endif /* HAVE_GPU */
#include "ObitBeamShape.h"
#include "ObitSkyModelVMBeamMF.h"
#include "ObitSkyModelVMSquint.h"
#include "ObitSkyModelVMIon.h"
#include "ObitSkyModelMF.h"
#include "ObitImageMF.h"
#include "ObitBeamShape.h"
#include "ObitSpectrumFit.h"
#include "ObitTableCCUtil.h"
#include "ObitFFT.h"
#include "ObitUVUtil.h"
#include "ObitImageUtil.h"
#include "ObitPBUtil.h"
#include "ObitMem.h"
#include "ObitPrecess.h"
#include "ObitSkyGeom.h"
#include "ObitSinCos.h"
#include "ObitExp.h"
#include "ObitMatx.h"
#include "ObitSpectrumMF.h"
#include "ObitUtil.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVMBeamMF.c
 * ObitSkyModelVMBeamMF class function definitions.
 *
 * This class is derived from the #ObitSkyModelVM class
 *
 * This class represents sky models incorporating beam corrections and 
 * their Fourier transforms.
 * Wideband imaging version.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSkyModelVMBeamMF";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitSkyModelVMBeamGetClass;

/**
 * ClassInfo structure ObitSkyModelVMBeamMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSkyModelVMBeamMFClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*----------------- Macroes ---------------------------*/
/** Half width of gridded subtraction interpolation kernal */
#define HWIDTH 12

/*---------------Private function prototypes----------------*/
/** Private: FT by DFT, may be overridden in derived class */
void ObitSkyModelVMBeamMFFTDFT (ObitSkyModelVM *in, olong field, 
				ObitUV *uvdata, ObitErr *err);

/** Private: Initialize newly instantiated object. */
void  ObitSkyModelVMBeamMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSkyModelVMBeamMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSkyModelVMBeamMFClassInfoDefFn (gpointer inClass);

/** Private: Get Inputs. */
void  ObitSkyModelVMBeamMFGetInput (ObitSkyModel* inn, ObitErr *err);
/** Private: Load Components */
gboolean ObitSkyModelVMBeamMFLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
					ObitErr *err);


/** Private: Threaded FTDFT */
static gpointer ThreadSkyModelVMBeamMFFTDFT (gpointer arg);

/** Private: Primary beam and/or smoothing corrections */
static ObitTableCC* getPBCCTab (ObitSkyModelVMBeamMF* in, ObitUV* uvdata, 
				olong field, olong *inCCVer, olong *outCCver,
				olong *startCC, olong *endCC, ofloat range[2],
				ObitErr *err);

/** Private: Need new beam gains? */
static gboolean needNewPB(ObitSkyModelVMBeamMF* in, ObitUV* uvdata, 
			  olong iIF, olong iChannel,
			  ofloat oldPB, ofloat *newPB, ObitErr *err);

/*---------------Private structures----------------*/
/* FT threaded function argument 
 Note: Derived classes MUST have the following entries at the beginning 
 of the corresponding structure 
 Note: this is identical to and labeled as the parent class 
 NB": First entries MUST match those in ObitSkyModelVMBeam */
typedef struct {
  /* type "mfvmbeam" in this class */
  gchar type[12];
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitSkyModelVMBeamMF *in;
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
  /* Thread copy of Components list */
  ObitFArray *VMComps;
  /* VMBeam class entries */
  /* Amp, phase interpolator for R/X pol Beam image array*/
  ObitFInterpolate **BeamRXInterp, **BeamRXPhInterp;
  /* Amp, phase interpolator for L/Y pol Beam image array*/
  ObitFInterpolate **BeamLYInterp, **BeamLYPhInterp;
  /* Amp, phase interpolator for RL/XY pol Beam image array*/
  ObitFInterpolate **BeamRLInterp, **BeamRLPhInterp;
  /* Amp, phase interpolator for LR/YX pol Beam image array*/
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
   /* Beam interpolation parameters of first CC - for debugging */
  ofloat beamX, beamY, beamPA, beamGain, beamWork[10];
  olong beamPlane;
  /** Are the Jones matrices (JMatrix) unit matrices? */
  gboolean isUnitJones;
  /** MF Spectrum Evaluator */
  ObitSpectrumMF *spEval;
  /* ObitSkyModelVMBeamMF components */
  /** Number of spectral bins */
  olong        nSpec;
  /** Apply prior alpha correction? */
  gboolean doAlphaCorr;
  /** Prior spectral index correction */
  ofloat priorAlpha;
  /** Reference frequency (Hz) for Prior spectral index */
  odouble priorAlphaRefF;
  /** dimension of chFlux */
  olong nChFlux;
  /** Channel flux densities */
  ofloat *chFlux;
} VMBeamMFFTFuncArg;

/** Private: Make model args structure */
void ObitSkyModelVMBeamMFMakeArg (ObitSkyModelVMBeamMF *in, ObitUV *uvdata, ObitErr *err);

/** Private: Update model args structure */
void  ObitSkyModelVMBeamMFUpdateArg (ObitSkyModelVMBeamMF *in, VMBeamMFFTFuncArg *args);

/** Private: Kill model arg structure */
gpointer ObitSkyModelVMBeamMFKillArg (ObitSkyModelVMBeamMF *in);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSkyModelVMBeamMF* newObitSkyModelVMBeamMF (gchar* name)
{
  ObitSkyModelVMBeamMF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMBeamMFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSkyModelVMBeamMF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSkyModelVMBeamMFInit((gpointer)out);

 return out;
} /* end newObitSkyModelVMBeamMF */

/**
 * Initializes from ObitInfoList.
 * Initializes class if needed on first call.
 * \param out     the new object.to be initialized
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string SkyModel type, "Squint" for this class
 *      \li "xxxThreshold" ofloat Threshold flux density for doing high accuracy DFT model
 *      \li "xxxClassType" string SkyModelVMBeamMF type, "MF" for base class
 *      \li "xxxmosaic"    string prefix of ObitImageMosaic mosaic
 *      \li "xxxmodelType" olong Model type (ObitSkyModelType)
 *      \li "xxxmodType"   olong Component model type (ObitSkyModelCompType)
 *      \li "xxxmodelMode" olong Model calculation mode for components (ObitSkyModelCompType)
 *      \li "xxxCCver"     olong* List of AIPSCC table versions per image in mosaic 
 *                                there are mosaic->numberImages of these
 *      \li "xxxstartComp" olong* List of beginning component per image in mosaic (1-rel)
 *      \li "xxxendComp"   olong* List of highest component per image in mosaic (1-rel)
 *      \li "xxxfactor"    ofloat Factor to multiply times model
 *      \li "xxxminFlux"   ofloat Minimum flux density model or pixel
 *      \li "xxxstokFactor"ofloat Factor to multiply times second Stokes of model
 *      \li "xxxpointFlux" ofloat Point model flux density (Jy)
 *      \li "xxxpointXOff" ofloat Point, x (ra)offset in deg.
 *      \li "xxxpointYOff" ofloat Point, y (dec) offset in deg.
 *      \li "xxxpointParms"ofloat[10] Other (non-point)model components:
 *                                major_axis (deg),  minor_axis (deg),  position_angle (deg),
 *                                type (ObitSkyModelCompType as gint), spectral terms;
 *      \li "xxxantSize"   ofloat Antennna diameter (m) for rel. PB corrections
 *      \li "xxxdo3D"            boolean Apply 3D imaging corrections?
 *      \li "xxxdoDivide"        boolean Divide model into data?
 *      \li "xxxdoReplace"       boolean Replace data with model?
 *      \li "xxxdoPBCor"         boolean Make relative Primary Beam corrections?
 *      \li "xxxstartChannel"    olong   Selected start channel[1-rel]
 *      \li "xxxnumberChannel"   olong   Selected channel and number 
 *      \li "xxxstartIF"         olong   Selected start IF [1-rel]
 *      \li "xxxnumberIF"        olong   Selected IF number
 *      \li "xxxstartChannelPB"  olong   Selected start rel. PB correction channel[1-rel]
 *      \li "xxxnumberChannelPB" olong   Selected PB correction channel number
 *      \li "xxxstartIFPB"       olong   Selected start rel. PB correction IF[1-rel]
 *      \li "xxxnumberIFPB"      olong   Selected PB correction IF number
 *      \li "xxxnfreqPB"         olong   number of frequency channels for PB correction
 *      \li "xxxPBFreq"          odouble Reference frequency (Hz) for this block of channels 
 *                                       for PB corrections
 *      \li "xxxstokes"          gchar[5] Selected Stokes
 *      \li "xxxstartPoln"       olong   Selected start Poln [1-rel]
 *      \li "xxxnumberPoln"      olong   Selected Poln number
 *      \li "xxxdoFlip"          boolean True if need to multiply the FT by sqrt(-1) before applying
 *      \li "xxxnoNeg"           boolean True if only positive flux components are to be used
 *      \li "xxxminDFT"          ofloat  Minimum absolute component flux to use in DFT
 *      \li "xxxdoDFT"           boolean Something to do for DFT model?
 *      \li "xxxprtLv"           olong   message level for progress messages
 *      \li "xxxnSpecTerm"       olong   Number of spectral terms
 *      \li "xxxnThreads"        olong   Number of threads
 *      \li "xxxdoAlphaCorr"     boolean TRUE if prior spectral index corrections to be made
 *      \li "xxxpriorAlpha"      ofloat  prior spectral index applied to be corrected.
 *      \li "xxxpriorAlphaRefF"  odouble prior spectral index ref freq (Hz).
 *      \li "xxxdoSmoo"          boolean TRUE if tabulated spectra to be smooothed
 * \param err     ObitErr for reporting errors.
 */
void ObitSkyModelVMBeamMFFromInfo (ObitSkyModel *out, gchar *prefix, ObitInfoList *inList, 
				   ObitErr *err)
{ 
  ObitSkyModelVMBeamMF *myOut = (ObitSkyModelVMBeamMF*)out;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *value=NULL, *classType="        ";
  gboolean missing;
  ObitImageMosaic *mosaic=NULL;
  olong classCnt, otemp;
  gchar *Type = "BeamMF";
  gchar ctemp[50];
  gchar *routine = "ObitSkyModelVMBeamMFFromInfo";
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMBeamMFClassInit();

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(out, &myClassInfo));

  /* check class type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  if ((missing) || (type!=OBIT_string) || (!strncmp(Type,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s Wrong class type %s!=%s", routine, value, Type);
    return;
  }
  classCnt = dim[0]; /* How many characters in name */
  g_free(keyword);

  /* "xxxThreshold" ofloat Threshold flux density for doing high accuracy DFT model */
  if (prefix) keyword = g_strconcat (prefix, "Threshold", NULL);
  else        keyword = g_strdup("Threshold");
  ObitInfoListGetTest(inList, keyword, &type, dim, &myOut->Threshold);
  g_free(keyword);

  /* "xxxmosaic" string prefix of ObitImageMosaic mosaic */
  if (prefix) keyword = g_strconcat (prefix, "mosaic", NULL);
  else        keyword = g_strdup("mosaic");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  /* Does it exist? */
  if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s ImageMosaic not defined in %s", 
		   routine, keyword);
    return;
  } else { /* exists*/
    mosaic = (ObitImageMosaic*)ObitImageMosaicFromInfo(keyword, inList, err);
    if (err->error) Obit_traceback_msg (err, routine, keyword);
  }
  g_free(keyword);

  /* Create output - by type */
  if (!strncmp("MF", classType, classCnt)) {
    out = (ObitSkyModel*)ObitSkyModelMFCreate(prefix, mosaic);
  } else if (!strncmp("mfvmbeam", classType, classCnt)) {
    /* This need more work - more arguments for create */
    g_error("FIX ME");
    /*out = (ObitSkyModelVMBeamMF*)ObitSkyModelVMBeamMFCreate(prefix, mosaic);
      ObitSkyModelVMBeamMFFromInfo((ObitSkyModel*)out, prefix, inList, err);*/
  } else if (!strncmp("Squint", classType, classCnt)) {
    out = (ObitSkyModel*)ObitSkyModelVMSquintCreate(prefix, mosaic);
    ObitSkyModelVMSquintFromInfo((ObitSkyModel*)out, prefix, inList, err);
  } else if (!strncmp("Ion", classType, classCnt)) {
    out = (ObitSkyModel*)ObitSkyModelVMIonCreate(prefix, mosaic);
    ObitSkyModelVMIonFromInfo((ObitSkyModel*)out, prefix, inList, err);
 } else {  /* Assume base and hope for the best */
    out = (ObitSkyModel*)ObitSkyModelCreate(prefix, mosaic);
    /* Note problem in log */
    strncpy (ctemp, classType, MIN (48,classCnt)); ctemp[MIN (49,classCnt+1)] = 0;
    Obit_log_error(err, OBIT_InfoWarn, "%s: Unknown type %s using base class",
		   routine, ctemp);
  }

  /* Copy any InfoList Parameters */
  if (prefix) keyword = g_strconcat (prefix, "Info", NULL);
  else        keyword = g_strdup("Info");
  ObitInfoListCopyWithPrefix (inList, out->info, keyword, TRUE);
  
  /* "xxxmodelType" olong Model type (ObitSkyModelType) */
  if (prefix) keyword = g_strconcat (prefix, "modelType", NULL);
  else        keyword = g_strdup("modelType");
  otemp = 0;
  ObitInfoListGetTest(inList, keyword, &type, dim, &otemp);
  out->modelType = (ObitSkyModelType)otemp;
  g_free(keyword);

  /* "xxxmodType"   olong Component model type (ObitSkyModelCompType) */
  if (prefix) keyword = g_strconcat (prefix, "modType", NULL);
  else        keyword = g_strdup("modType");
  otemp = 0;
  ObitInfoListGetTest(inList, keyword, &type, dim, &otemp);
  out->modType = (ObitSkyModelCompType)otemp;
  g_free(keyword);

  /* "xxxmodelMode" olong Model calculation mode for components (ObitSkyModelCompType) */
  if (prefix) keyword = g_strconcat (prefix, "modelMode", NULL);
  else        keyword = g_strdup("modelMode");
  otemp = 0;
  ObitInfoListGetTest(inList, keyword, &type, dim, &otemp);
  out->modelMode = (ObitSkyModelCompType)otemp;
  g_free(keyword);

  /* "xxxCCver"     olong* List of AIPSCC table versions per image in mosaic 
     there are mosaic->numberImages of these */
  if (prefix) keyword = g_strconcat (prefix, "CCver", NULL);
  else        keyword = g_strdup("CCver");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->CCver);
  g_free(keyword);

  /* "xxxstartComp" olong* List of beginning component per image in mosaic (1-rel) */
  if (prefix) keyword = g_strconcat (prefix, "startComp", NULL);
  else        keyword = g_strdup("startComp");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->startComp);
  g_free(keyword);

  /* "xxxendComp"   olong* List of highest component per image in mosaic (1-rel) */
  if (prefix) keyword = g_strconcat (prefix, "endComp", NULL);
  else        keyword = g_strdup("endComp");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->endComp);
  g_free(keyword);

  /* "xxxfactor"    ofloat Factor to multiply times model */
  if (prefix) keyword = g_strconcat (prefix, "factor", NULL);
  else        keyword = g_strdup("factor");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->factor);
  g_free(keyword);

  /* "xxxminFlux"   ofloat Minimum flux density model or pixel */
  if (prefix) keyword = g_strconcat (prefix, "minFlux", NULL);
  else        keyword = g_strdup("minFlux");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->minFlux);
  g_free(keyword);

  /* "xxxstokFactor"ofloat Factor to multiply times second Stokes of model */
  if (prefix) keyword = g_strconcat (prefix, "stokFactor", NULL);
  else        keyword = g_strdup("stokFactor");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->stokFactor);
  g_free(keyword);

  /* "xxxpointFlux" ofloat Point model flux density (Jy) */
  if (prefix) keyword = g_strconcat (prefix, "pointFlux", NULL);
  else        keyword = g_strdup("pointFlux");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->pointFlux);
  g_free(keyword);

  /* "xxxpointXOff" ofloat Point, x (ra)offset in deg. */
  if (prefix) keyword = g_strconcat (prefix, "pointXOff", NULL);
  else        keyword = g_strdup("pointXOff");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->pointXOff);
  g_free(keyword);

  /* "xxxpointYOff" ofloat Point, y (dec) offset in deg. */
  if (prefix) keyword = g_strconcat (prefix, "pointYOff", NULL);
  else        keyword = g_strdup("pointYOff");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->pointYOff);
  g_free(keyword);

  /* "xxxpointParms"ofloat[10] Other (non-point)model components: */
  if (prefix) keyword = g_strconcat (prefix, "pointParms", NULL);
  else        keyword = g_strdup("pointParms");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->pointParms);
  g_free(keyword);

  /* "xxxantSize"   ofloat Antennna diameter (m) for rel. PB corrections */
  if (prefix) keyword = g_strconcat (prefix, "antSize", NULL);
  else        keyword = g_strdup("antSize");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->antSize);
  g_free(keyword);

  /* "xxxdo3D"            boolean Apply 3D imaging corrections? */
  if (prefix) keyword = g_strconcat (prefix, "do3D", NULL);
  else        keyword = g_strdup("do3D");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->do3D);
  g_free(keyword);

  /* "xxxdoDivide"        boolean Divide model into data? */
  if (prefix) keyword = g_strconcat (prefix, "doDivide", NULL);
  else        keyword = g_strdup("doDivide");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doDivide);
  g_free(keyword);

  /* "xxxdoReplace"       boolean Replace data with model? */
  if (prefix) keyword = g_strconcat (prefix, "doReplace", NULL);
  else        keyword = g_strdup("doReplace");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doReplace);
  g_free(keyword);

  /* "xxxdoPBCor"         boolean Make relative Primary Beam corrections? */
  if (prefix) keyword = g_strconcat (prefix, "doPBCor", NULL);
  else        keyword = g_strdup("doPBCor");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doPBCor);
  g_free(keyword);

  /* "xxxstartChannel"    olong   Selected start channel[1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startChannel", NULL);
  else        keyword = g_strdup("startChannel");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->startChannel);
  g_free(keyword);

  /* "xxxnumberChannel"   olong   Selected channel and number  */
  if (prefix) keyword = g_strconcat (prefix, "numberChannel", NULL);
  else        keyword = g_strdup("numberChannel");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->numberChannel);
  g_free(keyword);

  /* "xxxstartIF"         olong   Selected start IF [1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startIF", NULL);
  else        keyword = g_strdup("startIF");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->startIF);
  g_free(keyword);

  /* "xxxnumberIF"        olong   Selected IF number */
  if (prefix) keyword = g_strconcat (prefix, "numberIF", NULL);
  else        keyword = g_strdup("numberIF");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->numberIF);
  g_free(keyword);

  /* "xxxstartChannelPB"  olong   Selected start rel. PB correction channel[1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startChannelPB", NULL);
  else        keyword = g_strdup("startChannelPB");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->startChannelPB);
  g_free(keyword);

  /* "xxxnumberChannelPB" olong   Selected PB correction channel number */
  if (prefix) keyword = g_strconcat (prefix, "numberChannelPB", NULL);
  else        keyword = g_strdup("numberChannelPB");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->numberChannelPB);
  g_free(keyword);

  /* "xxxstartIFPB"       olong   Selected start rel. PB correction IF[1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startIFPB", NULL);
  else        keyword = g_strdup("startIFPB");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->startIFPB);
  g_free(keyword);

  /* "xxxnumberIFPB"      olong   Selected PB correction IF number */
  if (prefix) keyword = g_strconcat (prefix, "numberIFPB", NULL);
  else        keyword = g_strdup("numberIFPB");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->numberIFPB);
  g_free(keyword);

  /* "xxxnfreqPB"         olong   number of frequency channels for PB correction */
  if (prefix) keyword = g_strconcat (prefix, "nfreqPB", NULL);
  else        keyword = g_strdup("nfreqPB");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->nfreqPB);
  g_free(keyword);

  /* "xxxPBFreq"          odouble Reference frequency (Hz) for this block of 
     channels for PB corrections  */
  if (prefix) keyword = g_strconcat (prefix, "PBFreq", NULL);
  else        keyword = g_strdup("PBFreq");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->PBFreq);
  g_free(keyword);

  /* "xxxstokes"          gchar[5] Selected Stokes */
  if (prefix) keyword = g_strconcat (prefix, "stokes", NULL);
  else        keyword = g_strdup("stokes");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->stokes);
  g_free(keyword);

  /* "xxxstartPoln"       olong   Selected start Poln [1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startPoln", NULL);
  else        keyword = g_strdup("startPoln");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->startPoln);
  g_free(keyword);

  /* "xxxnumberPoln"      olong   Selected Poln number */
  if (prefix) keyword = g_strconcat (prefix, "numberPoln", NULL);
  else        keyword = g_strdup("numberPoln");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->numberPoln);
  g_free(keyword);

  /* "xxxdoFlip"          boolean True if need to multiply the FT by sqrt(-1) 
     before applying */
  if (prefix) keyword = g_strconcat (prefix, "doFlip", NULL);
  else        keyword = g_strdup("doFlip");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doFlip);
  g_free(keyword);

  /* "xxxnoNeg"           boolean True if only positive flux components are to be used */
  if (prefix) keyword = g_strconcat (prefix, "noNeg", NULL);
  else        keyword = g_strdup("noNeg");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->noNeg);
  g_free(keyword);

  /* "xxxminDFT"          ofloat  Minimum absolute component flux to use 
     in DFT */
  if (prefix) keyword = g_strconcat (prefix, "minDFT", NULL);
  else        keyword = g_strdup("minDFT");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->minDFT);
  g_free(keyword);

  /* "xxxdoDFT"           boolean Something to do for DFT model? */
  if (prefix) keyword = g_strconcat (prefix, "doDFT", NULL);
  else        keyword = g_strdup("doDFT");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doDFT);
  g_free(keyword);

  /* "xxxprtLv"           olong   message level for progress messages */
  if (prefix) keyword = g_strconcat (prefix, "prtLv", NULL);
  else        keyword = g_strdup("prtLv");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->prtLv);
  g_free(keyword);

  /* "xxxnSpecTerm"       olong   Number of spectral terms */
  if (prefix) keyword = g_strconcat (prefix, "nSpecTerm", NULL);
  else        keyword = g_strdup("nSpecTerm");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->nSpecTerm);
  g_free(keyword);

  /* "xxxnThreads"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "nThreads", NULL);
  else        keyword = g_strdup("nThreads");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->nThreads);
  g_free(keyword);

  /* ""xxxdoAlphaCorr"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "doAlphaCorr", NULL);
  else        keyword = g_strdup("doAlphaCorr");
  ObitInfoListGetTest(inList, keyword, &type, dim, &myOut->doAlphaCorr);
  g_free(keyword);

  /* ""xxxpriorAlpha"    ofloat  prior spectral index applied to be corrected. */
  if (prefix) keyword = g_strconcat (prefix, "priorAlpha", NULL);
  else        keyword = g_strdup("priorAlpha");
  ObitInfoListGetTest(inList, keyword, &type, dim, &myOut->priorAlpha);
  g_free(keyword);

  /* ""xxxpriorAlphaRefF"    odouble prior spectral index ref freq (Hz). */
  if (prefix) keyword = g_strconcat (prefix, "priorAlpha", NULL);
  else        keyword = g_strdup("priorAlphaRefF");
  ObitInfoListGetTest(inList, keyword, &type, dim, &myOut->priorAlphaRefF);
  g_free(keyword);

  /* ""xxxdoSmoo"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "doSmoo", NULL);
  else        keyword = g_strdup("doSmoo");
  ObitInfoListGetTest(inList, keyword, &type, dim, &myOut->doSmoo);
  g_free(keyword);

 /* Cleanup */
  mosaic = ObitImageMosaicUnref(mosaic);

  return;
} /* end ObitSkyModelVMBeamMFFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSkyModelVMBeamMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMBeamMFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSkyModelVMBeamMFGetClass */

/**
 * Make a deep copy of an ObitSkyModelVMBeamMF.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSkyModelVMBeamMF* 
ObitSkyModelVMBeamMFCopy  (ObitSkyModelVMBeamMF *in, 
			   ObitSkyModelVMBeamMF *out, ObitErr *err)
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
  out->doAlphaCorr    = in->doAlphaCorr;
  out->priorAlpha     = in->priorAlpha;
  out->priorAlphaRefF = in->priorAlphaRefF;
  out->doSmoo         = in->doSmoo;

  return out;
} /* end ObitSkyModelVMBeamMFCopy */

/**
 * Creates an ObitSkyModelVMBeamMF 
 * \param name     An optional name for the object.
 * \param mosaic   ObitImageMosaic giving one or more images/CC tables
 * \param uvData   UV data to be operated on
 * \param numAntType number of antenna types
 * \param RXBeam   R/X Beam normalized image array per type
 * \param LYBeam   L/Y Beam normalized image array per type
 * \param RLBeam   RL/XY Beam fractional image array per type if nonNULL
 * \param LRBeam   LR/YX Beam fractional image array per type if nonNULL
 * \param RXBeamIm R/X Beam phase image array per type  if nonNULL
 * \param LYBeamIm L/Y Beam phase image array per type if nonNULL
 * \param RLBeamIm RL/XY Beam phase image array per type if nonNULL
 * \param LRBeamIm LR/YX Beam phase image array per type if nonNULL
 * \param Diams    Antenna diameters (m) per type
 * \return the new object.
 */
ObitSkyModelVMBeamMF* 
ObitSkyModelVMBeamMFCreate (gchar* name, ObitImageMosaic* mosaic,
			    ObitUV *uvData, olong numAntType,
			    ObitImage **RXBeam,    ObitImage **LYBeam, 
			    ObitImage **RLBeam,    ObitImage **LRBeam, 
			    ObitImage **RXBeamIm,  ObitImage **LYBeamIm, 
			    ObitImage **RLBeamIm,  ObitImage **LRBeamIm, 
			    ofloat *Diams, ObitErr *err)
{
  ObitSkyModelVMBeamMF* out=NULL;
  ObitImageDesc* beamDesc=NULL;
  ObitTableAN *ANTable=NULL;
  olong number, i, j, nchan, nif, refType, numOrb, numPCal, numIF, iANver, nbchan;
  olong iAnt, iDiam, found;
  ofloat refDiam, delt_f, crp_f;
  odouble f_0;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean doTab=TRUE, doCmplx;
  gchar *routine = "ObitSkyModelVMBeamMFCreate";

  /* Error tests */
  if (err->error) return out;  /* Previous error */

  /* Tell SkyMdel type */
  Obit_log_error(err, OBIT_InfoErr, "Using Sky Model with beam corrections");

  /* Create basic structure */
  out = newObitSkyModelVMBeamMF (name);

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
    out->Diams[i]  = Diams[i];
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
    if (doCmplx && out->RXBeamIm && out->RXBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->RXBeam[i]->ImgPixels, 
						   out->RXBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable pp amp, phase beam arrays", routine);
    }
    if (doCmplx && out->LYBeamIm && out->LYBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->LYBeam[i]->ImgPixels, 
						   out->LYBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable qq amp, phase beam arrays", routine);
    }
    if (out->doCrossPol && out->RLBeamIm && out->RLBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->RLBeam[i]->ImgPixels, 
						   out->RLBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable pq amp, phase beam arrays", routine);
    }
    if (out->doCrossPol && out->LRBeamIm && out->LRBeamIm[i]) {
      Obit_retval_if_fail ((ObitFArrayIsCompatable(out->LRBeam[i]->ImgPixels, 
						   out->LRBeamIm[i]->ImgPixels)), err, out,
			   "%s: Incompatable qp amp, phase beam arrays", routine);
    }

    ObitInfoListAlwaysPut (RXBeam[i]->info, "doTab", OBIT_bool, dim, &doTab);
  
  } /* end loop over ant type */

  refType = 0; refDiam=0.0;  /* Reference type is type with largest diameter */
  /* Find largest Antenna type diameter */
  /* Find smallest Antenna type diameter */
  refType = 0; refDiam=1.0e9;  /* Reference type is type with largest diameter */
  for (i=0; i<out->numAntType; i++) {
    /* Reference type */
    //if (out->Diams[i] > refDiam) {
    if (out->Diams[i] < refDiam) {
      refDiam = out->Diams[i]; refType = i;}
  }
  // HACK try geometric mean
  //refDiam = 14.2;
  Obit_log_error(err, OBIT_InfoErr, "Using Reference antenna diam %6.2f",refDiam);
  out->antSize = refDiam;
  /* Save antSiae on Infoist */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (out->info, "antSize", OBIT_float, dim, &out->antSize);

  /* AntList and AntType */
  out->numAntList = uvData->myDesc->numSubA;
  out->AntList = g_malloc0(out->numAntList*sizeof(ObitAntennaList*));
  out->AntType = g_malloc0(uvData->myDesc->maxAnt*sizeof(olong));
  for (i=0; i<out->numAntList; i++) {
    /* Loop over AN tables (subarrays) */
    numOrb = 0; numPCal = 0; numIF = 0; iANver = i+1;
    ANTable = newObitTableANValue ("AN table", (ObitData*)uvData, 
				   &iANver, OBIT_IO_ReadOnly, numIF, numOrb, numPCal, err);
    if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    out->AntList[i] = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
    /* Cleanup */
    ANTable = ObitTableANUnref(ANTable);
    /* Search for antennas and get type */
    found = -1;
    for (iAnt=1; iAnt<=uvData->myDesc->maxAnt; iAnt++) {
      /* Is iAnt (1-rel) in this one? */
      for (j=0; j<out->AntList[i]->number-1; j++) {
	if (out->AntList[i]->ANlist[j]->AntID==iAnt) {
	  found = j; break;  /* Yep */
	}
	/* Which diameter? */
	for (iDiam=0; iDiam<out->numAntType; iDiam++) {
	  if (fabs(out->AntList[i]->ANlist[found]->Diam-out->Diams[iDiam])<0.001) {
	    out->AntType[iAnt-1] = iDiam; break;  /* Yep */
	  }
	} /* end check diameter */
      } /* End loop over this list */
    } /* end loop checking antennas */
  }

  /* Reference Beam shape*/
  out->BeamShape = ObitBeamShapeCreate("Shape", RXBeam[refType], 0.01, out->antSize, TRUE);

  /* Get list of frequencies per beam channel */
  beamDesc = out->RXBeam[0]->ImgDesc;
  nbchan = beamDesc->inaxes[beamDesc->jlocf];  /* Number of beam channels */
  out->BeamFreq = g_malloc0(nbchan*sizeof(ofloat*));
  f_0 = beamDesc->crval[beamDesc->jlocf]; delt_f = beamDesc->cdelt[beamDesc->jlocf];
  crp_f = beamDesc->crpix[beamDesc->jlocf];
  for (i=0; i<nbchan; i++) out->BeamFreq[i] = (float)(f_0+(i+crp_f-1)*delt_f);

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
  doCmplx = RXBeamIm[0]!=NULL;  /* Have phases? */
  for (i=0; i<out->numAntType; i++) {
    if ((RXBeam!=NULL) && (RXBeam[i]->image!=NULL)) RXBeam[i]->image = ObitImageUnref(RXBeam[i]->image);
    if ((RLBeam!=NULL) && (RLBeam[i]->image!=NULL)) RLBeam[i]->image = ObitImageUnref(RLBeam[i]->image);
    if ((LRBeam!=NULL) && (LRBeam[i]->image!=NULL)) LRBeam[i]->image = ObitImageUnref(LRBeam[i]->image);
    if ((LYBeam!=NULL) && (LYBeam[i]->image!=NULL)) LYBeam[i]->image = ObitImageUnref(LYBeam[i]->image);
    if ((RXBeamIm!=NULL) && (RXBeamIm[i]->image!=NULL)) RXBeamIm[i]->image = ObitImageUnref(RXBeamIm[i]->image);
    if ((RLBeamIm!=NULL) && (RLBeamIm[i]->image!=NULL)) RLBeamIm[i]->image = ObitImageUnref(RLBeamIm[i]->image);
    if ((LRBeamIm!=NULL) && (LRBeamIm[i]->image!=NULL)) LRBeamIm[i]->image = ObitImageUnref(LRBeamIm[i]->image);
    if ((LYBeamIm!=NULL) && (LYBeamIm[i]->image!=NULL)) LYBeamIm[i]->image = ObitImageUnref(LYBeamIm[i]->image);
  }
  /* Set antenna Types */
  ObitSkyModelVMBeamSetAnt ((ObitSkyModelVMBeam*)out, uvData, err);

  /* Possibly using GPU? */
#if HAVE_GPU==1  /*  GPU? */
    out->GPUBeamModel = ObitGPUBeamModelCreate ("GPU", "DFT");
#endif /* HAVE_GPU */

    return out;
} /* end ObitSkyModelVMBeamMFCreate */

/**
 * Initializes Sky Model
 * Checks that data contain RR, LL , save calibration/selection request
 * and set uv data for no selection/calibration
 * \param in      SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamMFInitMod (ObitSkyModel* inn, ObitUV *uvdata, 
				  ObitErr *err)
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
  /*gchar *blank="    ";*/
  ObitImageMF *image0;
  olong i, j, nSpec, nif, nfreq, n;
  odouble test;
  ofloat phase=0.5, cp, sp;
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar keyword[24];
  gchar *routine = "ObitSkyModelVMBeamMFInitMod";

  if (err->error) return;

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

   /* Initialize threadArg array */
  if (in->threadArgs==NULL) ObitSkyModelVMBeamMFMakeArg(in, uvdata, err);

  /* Call parent initializer */
  ObitSkyModelVMBeamInitMod(inn, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fourier transform threading routines */
  in->DFTFunc  = (ObitThreadFunc)ThreadSkyModelVMBeamMFFTDFT;
  in->GridFunc = (ObitThreadFunc)NULL;
  
   /* Init Sine/Cosine, exp calculator - just to be sure about threading */
  ObitSinCosCalc(phase, &sp, &cp);
  ObitExpInit();

  /* Create spectrum info arrays */
  nSpec = 1; image0 = NULL;
  if ((ObitImageMF*)in->mosaic) {
    image0 = (ObitImageMF*)in->mosaic->images[0];	  
    ObitInfoListGetTest(image0->myDesc->info, "NSPEC", &type, dim, &nSpec);
    in->nSpec   = nSpec;
    in->refFreq = image0->myDesc->crval[image0->myDesc->jlocf];
    /* get number of and channel frequencies for CC spectra from 
       CC table on first image in mosaic */
    if (nSpec>1) {
      in->specFreq = g_malloc0(nSpec*sizeof(odouble));
      for (i=0; i<nSpec; i++) {
	in->specFreq[i] = 1.0;
	sprintf (keyword, "FREQ%4.4d",i+1);
	ObitInfoListGetTest(image0->myDesc->info, keyword, &type, dim, &in->specFreq[i]);
      }
    } else { /* Bummer */
      Obit_log_error(err, OBIT_Error,"%s No Frequency info in Image header for %s", 
		     routine, in->mosaic->images[0]->name);
      return;
    }
  } else { /* end if mosaic */
    in->specFreq = g_malloc0(nSpec*sizeof(odouble));
    in->specFreq[0] = uvdata->myDesc->crval[uvdata->myDesc->jlocf];
  } /* End no mosaic */

  /* Prior spectral index */
  InfoReal.flt = 0.0;   type = OBIT_float;
  if (image0) ObitInfoListGetTest(image0->myDesc->info, "ALPHA", &type, dim, &InfoReal);
  if (type==OBIT_double) in->priorAlpha = (ofloat)InfoReal.dbl;
  if (type==OBIT_float)  in->priorAlpha = (ofloat)InfoReal.flt;

  in->priorAlphaRefF = in->refFreq;
   if (image0) ObitInfoListGetTest(image0->myDesc->info, "RFALPHA", &type, dim, &in->priorAlphaRefF);
  
  /* Make array of which coarse spectrum value is closest to each uv channel */
  nfreq = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  if (uvdata->myDesc->jlocif>=0) 
    nif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  else nif = 1;
  n = nfreq*nif;
  in->specIndex = g_malloc0(n*sizeof(olong)); 
  for (i=0; i<n; i++) {
    test = 1.0e20;
    in->specIndex[i] = -1;
    for (j=0; j<nSpec; j++) {
      if (fabs(uvdata->myDesc->freqArr[i]- in->specFreq[j])<test) {
	test = fabs(uvdata->myDesc->freqArr[i]- in->specFreq[j]);
	in->specIndex[i] = j;
      }
    }
  } /* End of loop making lookup table */

  /* Tell selected model info if prtLv>1 */
  /* Only DFT */
  in->currentMode = OBIT_SkyModel_DFT;
  if (err->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr, "SkyModelVMBeamMF using DFT calculation type");
  }
#if HAVE_GPU==1  /*  GPU?*/
  /* Init GPU, create object if necessary */
  ofloat delta_PA = 0.25;  /* Allowed range in parallactic angle in degrees */
  if (in->doGPU) {
    if (in->GPUBeamModel==NULL)
      in->GPUBeamModel = ObitGPUBeamModelCreate ("GPU", "DFT");
    // can't derive a sub class, another argument is added (delta_PA)
    //const ObitGPUBeamModelClassInfo *GPUClass;
    //GPUClass = (ObitGPUBeamModelClassInfo*)in->GPUBeamModel->ClassInfo;
    if (in->modelMode==OBIT_SkyModel_DFT)
      ObitGPUBeamModelDFTInit(in->GPUBeamModel, (Obit*)in, uvdata, delta_PA, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
#endif /* HAVE_GPU */
} /* end ObitSkyModelVMBeamMFInitMod */

/**
 * Any shutdown operations needed for a model
 * Restore calibration/selection state
 * \param in  SkyModel to shutdown
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamMFShutDownMod (ObitSkyModel* inn, ObitUV *uvdata,
				      ObitErr *err)
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;

  in->myInterp = ObitCInterpolateUnref(in->myInterp);
  /* Delete args if the correct type */
  if (in->threadArgs) in->threadArgs = ObitSkyModelVMBeamMFKillArg (in);
  in->plane    = ObitFArrayUnref(in->plane);
  ObitThreadPoolFree (in->thread);  /* Shut down any threading */

  /* Cleanup arrays */
  if (in->specFreq)  {g_free(in->specFreq);}  in->specFreq = NULL;
  if (in->specIndex) {g_free(in->specIndex);} in->specIndex = NULL;
  
  /* Call parent shutdown */
  ObitSkyModelVMBeamShutDownMod(inn, uvdata, err);
  
#if HAVE_GPU==1  /*  GPU? */
  gchar *routine = "ObitSkyModelVMBeamMFShutDownMod";
  ObitGPUBeamModelDFTShutdown(in->GPUBeamModel, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  in->GPUBeamModel = ObitGPUSkyModelUnref (in->GPUBeamModel);
#endif /* HAVE_GPU */
} /* end ObitSkyModelVMBeamMFShutDownMod */

/**
 * Initializes an ObitSkyModel for a pass through data in time order.
 * Copies CLEAN model and the beam images to GPU
 * \param in  SkyModel to initialize
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamMFInitModel (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;

  /* Call parent InitModel */
  ObitSkyModelVMBeamInitModel (inn, err);

  /* Fourier transform routines - DFT only */
  /* Now only does full Jones */
  in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamMFFTDFT;
  
#if HAVE_GPU==1  /*  GPU? */
  gchar *routine = "ObitSkyModelvMBeamMFInitModel";
  const ObitGPUBeamModelClassInfo *GPUClass;
  if (in->doGPU) {
    GPUClass = (ObitGPUBeamModelClassInfo*)in->GPUBeamModel->ClassInfo;
    if ((in->modelMode==OBIT_SkyModel_DFT) && 
	((in->modType==OBIT_SkyModel_PointMod)      || 
	 (in->modType==OBIT_SkyModel_PointModSpec)  || 
	 (in->modType==OBIT_SkyModel_PointModTSpec) ||
	 (in->modType==OBIT_SkyModel_GaussMod)      || 
	 (in->modType==OBIT_SkyModel_GaussModSpec)  || 
	 (in->modType==OBIT_SkyModel_GaussModTSpec))) {
      /* CLEAN Model */
      GPUClass->ObitGPUBeamModelDFTSetMod (in->GPUBeamModel, (Obit*)in, in->comps, err);
      /* Beam Models */
      GPUClass->ObitGPUBeamModelDFTSetBeam (in->GPUBeamModel, (Obit*)in, err);
    } /* end if this enables in GPU */
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end if using GPU */
#endif /* HAVE_GPU */
} /* end ObitSkyModelVMBeamMFInitModel */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string SkyModel type, "Squint" for this class
 *      \li "xxxThreshold" ofloat Threshold flux density for doing high accuracy DFT model
 *      \li "xxxmosaic"    string prefix of ObitImageMosaic mosaic
 *      \li "xxxmodelType" olong Model type (ObitSkyModelType)
 *      \li "xxxmodType"   olong Component model type (ObitSkyModelCompType)
 *      \li "xxxmodelMode" olong Model calculation mode for components (ObitSkyModelCompType)
 *      \li "xxxCCver"     olong* List of AIPSCC table versions per image in mosaic 
 *                                there are mosaic->numberImages of these
 *      \li "xxxstartComp" olong* List of beginning component per image in mosaic (1-rel)
 *      \li "xxxendComp"   olong* List of highest component per image in mosaic (1-rel)
 *      \li "xxxfactor"    ofloat Factor to multiply times model
 *      \li "xxxminFlux"   ofloat Minimum flux density model or pixel
 *      \li "xxxstokFactor"ofloat Factor to multiply times second Stokes of model
 *      \li "xxxpointFlux" ofloat Point model flux density (Jy)
 *      \li "xxxpointXOff" ofloat Point, x (ra)offset in deg.
 *      \li "xxxpointYOff" ofloat Point, y (dec) offset in deg.
 *      \li "xxxpointParms"ofloat[10] Other (non-point)model components:
 *                                major_axis (deg),  minor_axis (deg),  position_angle (deg),
 *                                type (ObitSkyModelCompType as gint), spectral terms;
 *      \li "xxxantSize"   ofloat Antennna diameter (m) for rel. PB corrections
 *      \li "xxxdo3D"            boolean Apply 3D imaging corrections?
 *      \li "xxxdoDivide"        boolean Divide model into data?
 *      \li "xxxdoReplace"       boolean Replace data with model?
 *      \li "xxxdoPBCor"         boolean Make relative Primary Beam corrections?
 *      \li "xxxstartChannel"    olong   Selected start channel[1-rel]
 *      \li "xxxnumberChannel"   olong   Selected channel and number 
 *      \li "xxxstartIF"         olong   Selected start IF [1-rel]
 *      \li "xxxnumberIF"        olong   Selected IF number
 *      \li "xxxstartChannelPB"  olong   Selected start rel. PB correction channel[1-rel]
 *      \li "xxxnumberChannelPB" olong   Selected PB correction channel number
 *      \li "xxxstartIFPB"       olong   Selected start rel. PB correction IF[1-rel]
 *      \li "xxxnumberIFPB"      olong   Selected PB correction IF number
 *      \li "xxxnfreqPB"         olong   number of frequency channels for PB correction
 *      \li "xxxPBFreq"          odouble Reference frequency (Hz) for this block of channels 
 *                                       for PB corrections
 *      \li "xxxstokes"          gchar[5] Selected Stokes
 *      \li "xxxstartPoln"       olong   Selected start Poln [1-rel]
 *      \li "xxxnumberPoln"      olong   Selected Poln number
 *      \li "xxxdoFlip"          boolean True if need to multiply the FT by sqrt(-1) before applying
 *      \li "xxxnoNeg"           boolean True if only positive flux components are to be used
 *      \li "xxxminDFT"          ofloat  Minimum absolute component flux to use in DFT
 *      \li "xxxdoDFT"           boolean Something to do for DFT model?
 *      \li "xxxprtLv"           olong   message level for progress messages
 *      \li "xxxnSpecTerm"       olong   Number of spectral terms
 *      \li "xxxnThreads"        olong   Number of threads
 *      \li "xxxdoAlphaCorr"     boolean TRUE if prior spectral index corrections to be made
 *      \li "xxxpriorAlpha"      ofloat  prior spectral index applied to be corrected.
 *      \li "xxxpriorAlphaRefF"  odouble prior spectral index ref freq (Hz).
 *      \li "xxxdoSmoo"          boolean TRUE if tabulated spectra to be smooothed
 * \param err     ObitErr for reporting errors.
 */
void ObitSkyModelVMBeamMFGetInfo (ObitSkyModel *inn, gchar *prefix, 
				  ObitInfoList *outList, ObitErr *err)
{ 
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong numberImages=1, otemp;
  gchar *keyword=NULL, *Type="BeamMF", *OK="OK", *None = "None";
  gchar *routine = "ObitSkyModelVMBeamMFGetInfo";

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

  /* "xxxmosaic" string prefix of ObitImageMosaic mosaic */
  if (prefix) keyword = g_strconcat (prefix, "mosaic", NULL);
  else        keyword = g_strdup("mosaic");
  if (in->mosaic) {
    ObitImageMosaicGetInfo(in->mosaic, keyword, outList, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    numberImages = in->mosaic->numberImages;
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, OK);
  } else {
    dim[0] = strlen(None);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, None);
  }
  g_free(keyword);

  /* "xxxmodelType" olong Model type (ObitSkyModelType) */
  if (prefix) keyword = g_strconcat (prefix, "modelType", NULL);
  else        keyword = g_strdup("modelType");
  dim[0] = 1;
  otemp = (olong)in->modelType;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &otemp);
  g_free(keyword);

  /* "xxxmodType"   olong Component model type (ObitSkyModelCompType) */
  if (prefix) keyword = g_strconcat (prefix, "modType", NULL);
  else        keyword = g_strdup("modType");
  dim[0] = 1;
  otemp = (olong)in->modType;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &otemp);
  g_free(keyword);

  /* "xxxmodelMode" olong Model calculation mode for components (ObitSkyModelCompType) */
  if (prefix) keyword = g_strconcat (prefix, "modelMode", NULL);
  else        keyword = g_strdup("modelMode");
  dim[0] = 1;
  otemp = (olong)in->modelMode;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &otemp);
  g_free(keyword);

  /* "xxxCCver"     olong* List of AIPSCC table versions per image in mosaic 
     there are mosaic->numberImages of these */
  if (prefix) keyword = g_strconcat (prefix, "CCver", NULL);
  else        keyword = g_strdup("CCver");
  dim[0] = numberImages;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, in->CCver);
  g_free(keyword);

  /* "xxxstartComp" olong* List of beginning component per image in mosaic (1-rel) */
  if (prefix) keyword = g_strconcat (prefix, "startComp", NULL);
  else        keyword = g_strdup("startComp");
  dim[0] = numberImages;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, in->startComp);
  g_free(keyword);

  /* "xxxendComp"   olong* List of highest component per image in mosaic (1-rel) */
  if (prefix) keyword = g_strconcat (prefix, "endComp", NULL);
  else        keyword = g_strdup("endComp");
  dim[0] = numberImages;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, in->endComp);
  g_free(keyword);

  /* "xxxfactor"    ofloat Factor to multiply times model */
  if (prefix) keyword = g_strconcat (prefix, "factor", NULL);
  else        keyword = g_strdup("factor");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->factor);
  g_free(keyword);

  /* "xxxminFlux"   ofloat Minimum flux density model or pixel */
  if (prefix) keyword = g_strconcat (prefix, "minFlux", NULL);
  else        keyword = g_strdup("minFlux");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->minFlux);
  g_free(keyword);

  /* "xxxstokFactor"ofloat Factor to multiply times second Stokes of model */
  if (prefix) keyword = g_strconcat (prefix, "stokFactor", NULL);
  else        keyword = g_strdup("stokFactor");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->stokFactor);
  g_free(keyword);

  /* "xxxpointFlux" ofloat Point model flux density (Jy) */
  if (prefix) keyword = g_strconcat (prefix, "pointFlux", NULL);
  else        keyword = g_strdup("pointFlux");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->pointFlux);
  g_free(keyword);

  /* "xxxpointXOff" ofloat Point, x (ra)offset in deg. */
  if (prefix) keyword = g_strconcat (prefix, "pointXOff", NULL);
  else        keyword = g_strdup("pointXOff");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->pointXOff);
  g_free(keyword);

  /* "xxxpointYOff" ofloat Point, y (dec) offset in deg. */
  if (prefix) keyword = g_strconcat (prefix, "pointYOff", NULL);
  else        keyword = g_strdup("pointYOff");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->pointYOff);
  g_free(keyword);

  /* "xxxpointParms"ofloat[10] Other (non-point)model components: */
  if (prefix) keyword = g_strconcat (prefix, "pointParms", NULL);
  else        keyword = g_strdup("pointParms");
  dim[0] = 10;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, in->pointParms);
  g_free(keyword);

  /* "xxxantSize"   ofloat Antennna diameter (m) for rel. PB corrections */
  if (prefix) keyword = g_strconcat (prefix, "antSize", NULL);
  else        keyword = g_strdup("antSize");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->antSize);
  g_free(keyword);

  /* "xxxdo3D"            boolean Apply 3D imaging corrections? */
  if (prefix) keyword = g_strconcat (prefix, "do3D", NULL);
  else        keyword = g_strdup("do3D");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->do3D);
  g_free(keyword);

  /* "xxxdoDivide"        boolean Divide model into data? */
  if (prefix) keyword = g_strconcat (prefix, "doDivide", NULL);
  else        keyword = g_strdup("doDivide");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doDivide);
  g_free(keyword);

  /* "xxxdoReplace"       boolean Replace data with model? */
  if (prefix) keyword = g_strconcat (prefix, "doReplace", NULL);
  else        keyword = g_strdup("doReplace");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doReplace);
  g_free(keyword);

  /* "xxxdoPBCor"         boolean Make relative Primary Beam corrections? */
  if (prefix) keyword = g_strconcat (prefix, "doPBCor", NULL);
  else        keyword = g_strdup("doPBCor");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doPBCor);
  g_free(keyword);

  /* "xxxstartChannel"    olong   Selected start channel[1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startChannel", NULL);
  else        keyword = g_strdup("startChannel");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->startChannel);
  g_free(keyword);

  /* "xxxnumberChannel"   olong   Selected channel and number  */
  if (prefix) keyword = g_strconcat (prefix, "numberChannel", NULL);
  else        keyword = g_strdup("numberChannel");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->numberChannel);
  g_free(keyword);

  /* "xxxstartIF"         olong   Selected start IF [1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startIF", NULL);
  else        keyword = g_strdup("startIF");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->startIF);
  g_free(keyword);

  /* "xxxnumberIF"        olong   Selected IF number */
  if (prefix) keyword = g_strconcat (prefix, "numberIF", NULL);
  else        keyword = g_strdup("numberIF");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->numberIF);
  g_free(keyword);

  /* "xxxstartChannelPB"  olong   Selected start rel. PB correction channel[1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startChannelPB", NULL);
  else        keyword = g_strdup("startChannelPB");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->startChannelPB);
  g_free(keyword);

  /* "xxxnumberChannelPB" olong   Selected PB correction channel number */
  if (prefix) keyword = g_strconcat (prefix, "numberChannelPB", NULL);
  else        keyword = g_strdup("numberChannelPB");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->numberChannelPB);
  g_free(keyword);

  /* "xxxstartIFPB"       olong   Selected start rel. PB correction IF[1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startIFPB", NULL);
  else        keyword = g_strdup("startIFPB");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->startIFPB);
  g_free(keyword);

  /* "xxxnumberIFPB"      olong   Selected PB correction IF number */
  if (prefix) keyword = g_strconcat (prefix, "numberIFPB", NULL);
  else        keyword = g_strdup("numberIFPB");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->numberIFPB);
  g_free(keyword);

  /* "xxxnfreqPB"         olong   number of frequency channels for PB correction */
  if (prefix) keyword = g_strconcat (prefix, "nfreqPB", NULL);
  else        keyword = g_strdup("nfreqPB");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->nfreqPB);
  g_free(keyword);

  /* "xxxPBFreq"          odouble Reference frequency (Hz) for this block of 
     channels for PB corrections  */
  if (prefix) keyword = g_strconcat (prefix, "PBFreq", NULL);
  else        keyword = g_strdup("PBFreq");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_double, dim, &in->PBFreq);
  g_free(keyword);

  /* "xxxstokes"          gchar[5] Selected Stokes */
  if (prefix) keyword = g_strconcat (prefix, "stokes", NULL);
  else        keyword = g_strdup("stokes");
  dim[0] = 5;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, &in->stokes);
  g_free(keyword);

  /* "xxxstartPoln"       olong   Selected start Poln [1-rel] */
  if (prefix) keyword = g_strconcat (prefix, "startPoln", NULL);
  else        keyword = g_strdup("startPoln");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->startPoln);
  g_free(keyword);

  /* "xxxnumberPoln"      olong   Selected Poln number */
  if (prefix) keyword = g_strconcat (prefix, "numberPoln", NULL);
  else        keyword = g_strdup("numberPoln");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->numberPoln);
  g_free(keyword);

  /* "xxxdoFlip"          boolean True if need to multiply the FT by sqrt(-1) 
     before applying */
  if (prefix) keyword = g_strconcat (prefix, "doFlip", NULL);
  else        keyword = g_strdup("doFlip");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doFlip);
  g_free(keyword);

  /* "xxxnoNeg"           boolean True if only positive flux components are to be used */
  if (prefix) keyword = g_strconcat (prefix, "noNeg", NULL);
  else        keyword = g_strdup("noNeg");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->noNeg);
  g_free(keyword);

  /* "xxxminDFT"          ofloat  Minimum absolute component flux to use 
     in DFT */
  if (prefix) keyword = g_strconcat (prefix, "minDFT", NULL);
  else        keyword = g_strdup("minDFT");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->minDFT);
  g_free(keyword);

  /* "xxxdoDFT"           boolean Something to do for DFT model? */
  if (prefix) keyword = g_strconcat (prefix, "doDFT", NULL);
  else        keyword = g_strdup("doDFT");
  dim[0] = 1;
  in->doDFT = TRUE;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doDFT);
  g_free(keyword);

  /* "xxxprtLv"           olong   message level for progress messages */
  if (prefix) keyword = g_strconcat (prefix, "prtLv", NULL);
  else        keyword = g_strdup("prtLv");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->prtLv);
  g_free(keyword);

  /* "xxxnSpecTerm"       olong   Number of spectral terms */
  if (prefix) keyword = g_strconcat (prefix, "nSpecTerm", NULL);
  else        keyword = g_strdup("nSpecTerm");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->nSpecTerm);
  g_free(keyword);

  /* "xxxnThreads"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "nThreads", NULL);
  else        keyword = g_strdup("nThreads");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->nThreads);
  g_free(keyword);

  /* "xxxdoAlphaCorr"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "doAlphaCorr", NULL);
  else        keyword = g_strdup("doAlphaCorr");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doAlphaCorr);
  g_free(keyword);

  /* "xxxpriorAlpha"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "priorAlpha", NULL);
  else        keyword = g_strdup("priorAlpha");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->priorAlpha);
  g_free(keyword);

  /* "xxxpriorAlphaRefF"        odouble prior spectral index ref freq (Hz) */
  if (prefix) keyword = g_strconcat (prefix, "priorAlphaRefF", NULL);
  else        keyword = g_strdup("priorAlphaRefF");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_double, dim, &in->priorAlphaRefF);
  g_free(keyword);

  /* "xxxdoSmoo"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "doSmoo", NULL);
  else        keyword = g_strdup("doSmoo");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doSmoo);
  g_free(keyword);

} /* end ObitSkyModelVMBeamMFGetInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSkyModelVMBeamMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSkyModelVMBeamMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSkyModelVMBeamMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSkyModelVMBeamMFClassInfoDefFn (gpointer inClass)
{
  ObitSkyModelVMBeamMFClassInfo *theClass = (ObitSkyModelVMBeamMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSkyModelVMBeamMFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSkyModelVMBeamMFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSkyModelVMBeamMFGetClass;
  theClass->newObit       = (newObitFP)newObitSkyModelVMBeamMF;
  theClass->ObitCopy      = (ObitCopyFP)ObitSkyModelVMBeamMFCopy;
  theClass->ObitClear     = (ObitClearFP)ObitSkyModelVMBeamMFClear;
  theClass->ObitInit      = (ObitInitFP)ObitSkyModelVMBeamMFInit;
  theClass->ObitSkyModelCreate       = (ObitSkyModelCreateFP)ObitSkyModelVMBeamMFCreate;
  theClass->ObitSkyModelInitMod      = (ObitSkyModelInitModFP)ObitSkyModelVMBeamMFInitMod;
  theClass->ObitSkyModelShutDownMod  = (ObitSkyModelShutDownModFP)ObitSkyModelVMBeamMFShutDownMod;
  theClass->ObitSkyModelInitModel    = (ObitSkyModelInitModelFP)ObitSkyModelVMBeamMFInitModel;
  theClass->ObitSkyModelLoadComps    = (ObitSkyModelLoadCompsFP)ObitSkyModelVMBeamMFLoadComps;
  theClass->ObitSkyModelFTDFT        = (ObitSkyModelFTDFTFP)ObitSkyModelVMBeamMFFTDFT;
  theClass->ObitSkyModelGetInput     = (ObitSkyModelGetInputFP)ObitSkyModelVMBeamMFGetInput;
  theClass->ObitSkyModelChose        = (ObitSkyModelChoseFP)ObitSkyModelVMBeamMFChose;
  theClass->ObitSkyModelGetInfo= (ObitSkyModelGetInfoFP)ObitSkyModelVMBeamMFGetInfo;

} /* end ObitSkyModelVMBeamMFClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSkyModelVMBeamMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModelVMBeamMF *in = inn;

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
} /* end ObitSkyModelVMBeamMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSkyModelVMBeamMF* cast to an Obit*.
 */
void ObitSkyModelVMBeamMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitSkyModelVMBeamMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->numPlane) {g_free(in->numPlane);} in->numPlane = NULL;
  if (in->AntType)  {g_free(in->AntType);}  in->AntType  = NULL;
  if (in->BeamFreq)  {g_free(in->BeamFreq);}  in->BeamFreq  = NULL;
  in->BeamShape = ObitBeamShapeUnref(in->BeamShape);
  for (i=0; i<in->numAntType; i++) {
    if (in->RXBeam && in->RXBeam[i]) in->RXBeam[i]    = ObitImageInterpUnref(in->RXBeam[i]);
    if (in->RLBeam && in->RLBeam[i]) in->RLBeam[i]    = ObitImageInterpUnref(in->RLBeam[i]);
    if (in->LRBeam && in->LRBeam[i]) in->LRBeam[i]    = ObitImageInterpUnref(in->LRBeam[i]);
    if (in->LYBeam && in->LYBeam[i]) in->LYBeam[i]    = ObitImageInterpUnref(in->LYBeam[i]);
    if (in->RXBeamIm && in->RXBeamIm[i]) in->RXBeamIm[i]  = ObitImageInterpUnref(in->RXBeamIm[i]);
    if (in->RLBeamIm && in->RLBeamIm[i]) in->RLBeamIm[i]  = ObitImageInterpUnref(in->RLBeamIm[i]);
    if (in->LRBeamIm && in->LRBeamIm[i]) in->LRBeamIm[i]  = ObitImageInterpUnref(in->LRBeamIm[i]);
    if (in->LYBeamIm && in->LYBeamIm[i]) in->LYBeamIm[i]  = ObitImageInterpUnref(in->LYBeamIm[i]);
  }
  if (in->RXBeam)   {g_free(in->RXBeam);    in->RXBeam    = NULL;}
  if (in->RLBeam)   {g_free(in->RLBeam);    in->RLBeam    = NULL;}
  if (in->LRBeam)   {g_free(in->LRBeam);    in->LRBeam    = NULL;}
  if (in->LYBeam)   {g_free(in->LYBeam);    in->LYBeam    = NULL;}
  if (in->RXBeamIm) {g_free(in->RXBeamIm);  in->RXBeamIm  = NULL;}
  if (in->RLBeamIm) {g_free(in->RLBeamIm);  in->RLBeamIm  = NULL;}
  if (in->LRBeamIm) {g_free(in->LRBeamIm);  in->LRBeamIm  = NULL;}
  if (in->LYBeamIm) {g_free(in->LYBeamIm);  in->LYBeamIm  = NULL;}
  in->curSource  = ObitSourceUnref(in->curSource);
  if (in->AntList)  {
    for (i=0; i<in->numAntList; i++) { 
      in->AntList[i] = ObitAntennaListUnref(in->AntList[i]);
    }
    g_free(in->AntList); in->AntList = NULL;
  }
    
  /* Thread stuff */
  if (in->threadArgs) ObitSkyModelVMBeamMFKillArg(in);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSkyModelVMBeamMFClear */

/**
 * Get input parameters from info member
 * If maxResid value not given or <0 then for it is determined
 * by examining each image in the Image mosaic.  
 * If any image has an info item "maxAbsResid" the the maximum of any
 * of these is used, else the MAX (fabs(minval), fabs(maxval)) is used.
 * \param in Pointer to the ObitSkyModelVMBeamMF .
 * \param err Obit error stack object.
 */
void  ObitSkyModelVMBeamMFGetInput (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitSkyModelVMBeamMFGetInput";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelVMBeamMFIsA(in));
  if (!ObitInfoListIsA(in->info)) return;
  InfoReal.itg = 0;type = OBIT_oint;

  /* Call base class version */
  ObitSkyModelVMBeamGetInput (inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Prior Alpha correction wanted? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doAlphaCorr", &type, (gint32*)dim, &InfoReal);
  in->doAlphaCorr = InfoReal.itg;

  /* Smoothing of flux densities wanted? */
  InfoReal.itg = (olong)FALSE; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "doSmoo", &type, (gint32*)dim, &InfoReal);
  in->doSmoo = InfoReal.itg;

  /* Turn off request for frequency dependent primary beam correction,
     as currently implemented, it is not correct as it should only make a 
     correction wrt the channels and IFs being imaged together and NOT the
     total set of frequencies. */
  in->doPBCor = FALSE;

} /* end ObitSkyModelVMBeamMFGetInput */

/**
 * Decide which method is the most appropriate to calculate the FT of a model
 * Sets currentMode member function
 * Only DFT supported
 * Threshold related decisioins made when loading components.
 * \param in     Pointer to the ObitSkyModel .
 * \param uvdata UV data set
 */
void  ObitSkyModelVMBeamMFChose (ObitSkyModel* inn, ObitUV* uvdata) 
{
  ObitSkyModelVMBeam *in = (ObitSkyModelVMBeam*)inn;
  /* Only OBIT_SkyModel_DFT supported */
  in->currentMode = OBIT_SkyModel_DFT;
} /* end ObitSkyModelVMBeamMFChose */

/**
 * Average subband flux densities
 * \param n   Number of subband channels
 * \param SBFlux Subband flux densities, 0=> no value
 * \return average, 1 if no valid data
 */
ofloat  ObitSkyModelVMBeamMFAvgSpec (olong n, ofloat* SBFlux) 
{
  ofloat sum;
  olong i, cnt;
  sum = 0.0; cnt = 0;
  for (i=0; i<n; i++) {
    if (SBFlux[i]!=0.0) {sum += SBFlux[i]; cnt++;}
  }
  if (cnt>0) sum /= cnt;
  else       sum = 1.0;
  return sum;
} /* end ObitSkyModelVMBeamMFAvgSpec */

/**
 * Load components model into in comps member.
 * If the frequency axis has ctype "SPECLNMF" and if "NSPEC" exists in the 
 * first image descriptor InfoList and is > 0 then there should be a spectrum 
 * in each CLEAN component.
 * Multiplies by factor member and any prior spectral index correction.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Allows mixed points and Gaussians, one per CC Table
 * The summ of the component flux densities is made and used to decide whether
 * Beam or Unit Jones matrices is made
 * \param inn  SkyModel 
 * \param n    Image number on mosaic, if -1 load all images
 * \param      uvdata UV data set to model
 * \param      err Obit error stack object.
 * Output is in member comps, the entries are
 * \li 0    Reserved to evaluated flux density flux density if no spectrum
 * \li 1    Field (0-rel)
 * \li 2    CC DeltaX for beam correction
 * \li 3    CC DeltaY
 * \li 4    -2*pi*x (radians),
 * \li 5    -2*pi*y (radians),
 * \li 6    -2*pi*z (radians),
 * \li 7->6+nTerm  nTerm Spectral terms: flux, si, curv1, curv2...
 * \li 7+nTerm->6+nTerm+nSpec Subband spectra, nSpec entries
 * \li    Gauss: 7+nadd->6+nadd+3 Gaussian terms  nadd = 7+nTerm+nSpec             
 * \li    Unif. Sp: 7+nadd->6+nadd+2 : unif. sp. factor, 0.1
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelVMBeamMFLoadComps (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
					ObitErr *err)
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *imDesc=NULL, *imIODesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  ObitFArray *CompArr=NULL;
  ObitSkyModelCompType modType, maxModType=OBIT_SkyModel_PointMod;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong warray, larray, iterm, nterm, nplane=1,  maxTerm=1, toff;
  ofloat *array, parms[20], range[2], gp1=0., gp2=0., gp3=0.;
  olong ver, i, j, k, kk, indxsb, indxfs, hi, lo, count, ncomp, startComp, endComp, irow, lrec;
  olong outCCVer, ndim, naxis[2], lenEntry, nadd, nspec, nsterm;
  ofloat *table, xxoff, yyoff, zzoff, maxCCFlux=0.0;
  ofloat konst, konst2, xyz[3], xp[3], umat[3][3], pmat[3][3];
  ofloat ccrot, ssrot, xpoff, ypoff, maprot, uvrot, wfact, allRMS=1.0e-6;
  ofloat dxyzc[3], cpa, spa, xmaj, xmin, *specCorr=NULL;
  gboolean doCheck=FALSE, want, do3Dmul, noNeg, noSIfit, badFit=FALSE;
  odouble *fitFreq = NULL, sumCCFlux=0.0;
  ofloat  *fitFlux = NULL, *PlnRMS = NULL;
  gpointer fitArg=NULL;
  ofloat *fitSigma=NULL, *fitParms=NULL;
  ofloat  fblank = ObitMagicF();
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitSkyModelVMBeamMFLoadComps";
  
  /* error checks */
  if (err->error) return gotSome;

  /* Don't bother if no components requested */
  if ((n>=0) && (in->startComp[n]>in->endComp[n])) return gotSome;

  /* Uv descriptor */
  uvDesc = uvdata->myDesc;

    /* Don't fit spectral index for polarization */
  noSIfit = uvDesc->crval[uvdata->myDesc->jlocs]>1.;
  noSIfit = noSIfit || (uvDesc->crval[uvdata->myDesc->jlocs]==-3.);
  noSIfit = noSIfit || (uvDesc->crval[uvdata->myDesc->jlocs]==-4.);
  noSIfit = noSIfit || (uvDesc->crval[uvdata->myDesc->jlocs]<-6.);
  imDesc  = in->mosaic->images[0]->myDesc;
  noSIfit = noSIfit || (imDesc->crval[imDesc->jlocs]>=2.);

  konst = DG2RAD * 2.0 * G_PI;
  /* konst2 converts FWHM(deg) to coefficients for u*u, v*v, u*v */
  /*konst2 = DG2RAD * 2.15169;*/
  /* konst2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5);*/
  konst2 = DG2RAD * sqrt(2.0) * G_PI / 2.35482044;

  /* Loop over images counting CCs */
  count = 0;
  in->modType = OBIT_SkyModel_Unknown; /* Model type not known */
  if (in->mosaic) {lo = 0; hi = in->mosaic->numberImages-1;}
  else {lo = 0; hi = 0;}
  if (n>=0) {lo = n; hi = n;}
  for (i=lo; i<=hi; i++) {

    /* Expect anything in this table? */
    if ((in->startComp[i]>in->endComp[i]) || (in->endComp[i]<=0)) continue;

    /* Get RMSes from first image with components */
    if ((PlnRMS==NULL) && !noSIfit) {
      fitSigma = g_malloc0(in->nSpec*sizeof(ofloat));  /* Allocate */
      ObitImageUtilGetRMS(in->mosaic->images[i], err);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      /* Where do subband spectra start? */
      if (!ObitInfoListGetTest(in->mosaic->images[i]->myDesc->info, "NTERM", &type, dim, &nsterm))
	nsterm = 2;/* Default Number of spectral terms in image */
      //ObitInfoListPrint (in->mosaic->images[i]->info,stderr); // SERIOUS DEBUG
      if (ObitInfoListGetP(in->mosaic->images[i]->info, "PlnRMS", &type, dim, (gpointer*)&PlnRMS)) {
	nplane = dim[0];
	for (kk=0; kk<in->nSpec; kk++) fitSigma[kk] = PlnRMS[nsterm+kk];
	allRMS =  PlnRMS[0];  /* Broadband  RMS */
	// DEBUG
	//fprintf (stderr,"Fit RMS[1]=%10.5g RMS[9]=%10.5g RMS[14]=%10.5g, nsterm=%d\n",
	//	 fitSigma[0],fitSigma[8],fitSigma[13],nsterm);
      } else {  /* Something wrong, use default */
	for (kk=0; kk<in->nSpec; kk++) fitSigma[kk] = 0.0001;
	// DEBUG
	//fprintf (stderr,"Default RMS[0]=%10.5g RMS[14]=%10.5g\n",fitSigma[0],fitSigma[13]);
      }
    }  /* end getting plane RMSes */

    /* Save PlnRMS on skymodel if present */
    if (PlnRMS) {
      dim[0] = nplane; dim[1] = dim[2] = dim[3] = dim[4] = 1;
      ObitInfoListAlwaysPut(in->info, "PlnRMS", OBIT_float, dim, PlnRMS);
    } /* end save PlnRMS */

    /* Get CC table */
    ver = in->CCver[i];
    tempTable = newObitImageTable (in->mosaic->images[i],OBIT_IO_ReadOnly, 
				   tabType, &ver, err);
    if ((tempTable==NULL) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    /* Open */
    retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);

    /* How many? */
    endComp = in->endComp[i];
    if (endComp<=0) endComp = CCTable->myDesc->nrow;
    count += MIN(CCTable->myDesc->nrow, endComp) - MAX(1, in->startComp[i]) + 1;

    /* Get model type in first with components */
    /* If only 3 col, or parmsCol 0 size then this is a point model */
    if ((CCTable->myDesc->nfield==3) || 
	(CCTable->parmsCol<0) ||
	(CCTable->myDesc->dim[CCTable->parmsCol]<=0)) {
      if (in->startComp[i]<=endComp) {
	in->modType = OBIT_SkyModel_PointMod;
	maxModType =  MAX (maxModType, in->modType);
      }
      /* Check type of all tables with components to subtract */
    } else if (in->startComp[i]<=endComp) {
      /* Create table row */
      CCRow = newObitTableCCRow (CCTable);
      /* Read first */
      irow = in->startComp[i];
      retCode = ObitTableCCReadRow (CCTable, irow, CCRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error)) 
	Obit_traceback_val (err, routine, in->name, retCode);

      /* Get model type */
      in->modType = CCRow->parms[3] + 0.5;
      maxModType  =  MAX (maxModType, in->modType);
      /* Release table row */
      CCRow = ObitTableCCRowUnref (CCRow);
    }

    /* Do we need to check model type */
    doCheck = doCheck || ((CCTable->myDesc->nfield>4) && (CCTable->parmsCol>=0) && 
      (CCTable->myDesc->dim[CCTable->parmsCol][0]>=3));
    
    /* Close */
    retCode = ObitTableCCClose (CCTable, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);

    /* release table  */
    CCTable = ObitTableCCUnref (CCTable);

    /* Is paramaterized spectral information included as opposed to a tabulated spectrum? */
    imDesc = in->mosaic->images[i]->myDesc; /* Image descriptor */
    if (!strncmp (imDesc->ctype[imDesc->jlocf],  "SPECLOGF", 8)) {
      nterm = imDesc->inaxes[imDesc->jlocf];
      ObitInfoListGetTest (imDesc->info, "NTERM", &type, dim, &nterm);
      maxTerm = MAX (maxTerm, nterm);
      in->nSpecTerm = nterm;     /* Number of terms in paramaterized spectrum */
      in->nTerm = nterm;         /* Number of terms in paramaterized spectrum */
      in->nSpec = 0;             /* No subband spectrum */
   }
    /* Tabulated spectrum? */
    if (!strncmp (imDesc->ctype[imDesc->jlocf],  "SPECLNMF", 8)) {
      /* IO descriptor gives true size */
      imIODesc = (ObitImageDesc*)in->mosaic->images[0]->myIO->myDesc;
      nspec = imIODesc->inaxes[imIODesc->jlocf];
      ObitInfoListGetTest (imDesc->info, "NSPEC", &type, dim, &nspec);
      in->nSpec = nspec;     /* Number of subbands in the spectrum */
      maxTerm = MAX (maxTerm, nspec);
      in->nTerm = 4;         /* Number of terms to use in spectral fit */
   }
  } /* end loop counting CCs */

  /* Use mode type, no subbands, spectral terms of the highest encountered */
  in->modType   = maxModType;
  if ((in->modType==OBIT_SkyModel_PointModSpec) || (in->modType==OBIT_SkyModel_GaussModSpec) || 
      (in->modType==OBIT_SkyModel_USphereModSpec)) {
    in->nSpecTerm = MAX (0, maxTerm-1);
  } else  in->nSpecTerm = 0;
  if ((in->modType==OBIT_SkyModel_PointModTSpec) || (in->modType==OBIT_SkyModel_GaussModTSpec) || 
      (in->modType==OBIT_SkyModel_USphereModTSpec)) {
    in->nSpec     = MAX (in->nSpec, maxTerm);
  } else in->nSpec = 0;

  /* (re)allocate structure */
  ndim = 2;
  naxis[0] = 7+in->nTerm+in->nSpec;
  naxis[1] = count;
  if (in->modType==OBIT_SkyModel_GaussMod)        naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_GaussModSpec)    naxis[0] += 3; /* Gaussian + spectrum */
  if (in->modType==OBIT_SkyModel_GaussModTSpec)   naxis[0] += 3; /* Gaussian  + tabuated spectrum */
  if (in->modType==OBIT_SkyModel_USphereMod)      naxis[0] += 2; /* Uniform sphere */
  if (in->modType==OBIT_SkyModel_USphereModSpec)  naxis[0] += 2; /* Uniform sphere + spectrum*/
  if (in->modType==OBIT_SkyModel_USphereModTSpec) naxis[0] += 2; /* Uniform sphere + tabuated spectrum */
  lenEntry = naxis[0];  /* Length of table entry */
  nadd = 7+in->nTerm+in->nSpec;   /* Beginning of additional parameters */

  /* Fitting component spectra? Create work arrays */
  if ((!noSIfit) && (in->nSpec>1)) {
    fitFreq  = g_malloc0(in->nSpec*sizeof(odouble));
    fitFlux  = g_malloc0(in->nSpec*sizeof(ofloat));
  }  /* End create work arrays */
  
  /* Spectral correction for prior alpha array (with any factor) */
  specCorr = g_malloc0(in->nSpec*sizeof(ofloat));
  if (in->doAlphaCorr && (in->priorAlpha!=0.0)) {
    for (i=0; i<in->nSpec; i++) {
      specCorr[i] = pow((in->specFreq[i]/in->priorAlphaRefF), in->priorAlpha) * in->factor;
    }
  } else { /* No correction */
    for (i=0; i<in->nSpec; i++) specCorr[i] = in->factor;
  }

  if (in->comps!=NULL) in->comps = ObitFArrayRealloc(in->comps, ndim, naxis);
  else in->comps = ObitFArrayCreate("Components", ndim, naxis);
  lrec = naxis[0]; /* Save size of entry */
  /* Get pointer */
  naxis[0] = 0; naxis[1]=0; 
  table = ObitFArrayIndex(in->comps, naxis);

  /* Loop over images loading CCs */
  ncomp = 0; sumCCFlux = 0.0;
  for (i=lo; i<=hi; i++) {

    /* Anything to do? */
    if ((in->endComp[i]<=0) || (in->endComp[i]<in->startComp[i])) continue;

    /* Get CC table */
    outCCVer = 0;
    ver = in->CCver[i];
    startComp = in->startComp[i];
    endComp = in->endComp[i];
    range[0] = in->minDFT;  /* Range of merged fluxes for DFT */
    range[1] = 1.0e20;
    CCTable = getPBCCTab (in, uvdata, (olong)i, &ver, &outCCVer, 
			  &startComp, &endComp, range, err); 
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    /* Anybody home? */
    if (CCTable==NULL) continue;
    
    /* Anything to do? */
    if (endComp<startComp) {
      /* No - free up tables */
      /* if outCCver>0 then the CCtable is temporary - Zap */
      if (outCCVer>0) {
	CCTable = ObitTableCCUnref (CCTable);
	ObitImageZapTable(in->mosaic->images[i], tabType, outCCVer, err);
      /* else simply release table  */
      } else CCTable = ObitTableCCUnref (CCTable);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      continue;
    }

    /* Field specific stuff */
    imDesc = in->mosaic->images[i]->myDesc; /* Image descriptor */
    /*  Set field center offsets. */
    maprot = ObitImageDescRotate(imDesc);
    uvrot  = ObitUVDescRotate(uvDesc);
    ssrot = sin (DG2RAD * (uvrot - maprot));
    ccrot = cos (DG2RAD * (uvrot - maprot));
    /* Get do3D from the image */
    in->do3D = imDesc->do3D;

    /* noNeg FALSE for Stokes != I */
    if ((fabs(imDesc->crval[imDesc->jlocs])-1.0)>0.01) noNeg = FALSE;
    else noNeg = in->noNeg;
    
    /* Get position phase shift parameters */
    ObitUVDescShiftPhase(uvDesc, imDesc, dxyzc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
    /*    Get reference pixel offsets from tangent point */
    if (in->do3D) {
      /* These should always be zero for 3D imaging? */
      xpoff = 0.0;
      ypoff = 0.0;
    } else { /** 2D - use offsets */
      xpoff = imDesc->xPxOff * imDesc->cdelt[imDesc->jlocr];
      ypoff = imDesc->yPxOff * imDesc->cdelt[imDesc->jlocd];
       /* ypoff = (imDesc->yPxOff+1.0) * imDesc->cdelt[imDesc->jlocd];DEBUG */
    }

    /* Set field center offsets */
    xxoff = dxyzc[0] * ccrot + dxyzc[1] * ssrot;
    yyoff = dxyzc[1] * ccrot - dxyzc[0] * ssrot;
    zzoff = dxyzc[2];

    /* rotation matrix if needed */
    do3Dmul = ObitUVDescShift3DMatrix (uvDesc, imDesc, umat, pmat);
    
    /* Convert table to merged array */
    CompArr = ObitTableCCUtilMergeSel (CCTable, startComp, endComp, parms, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    /* entries 0=flux, 1= deltaX 2=deltaY per merged CC, other parameters in parms */
    naxis[0] = 0; naxis[1]=0; 
    array = ObitFArrayIndex(CompArr, naxis);
    warray = CompArr->naxis[0];
    larray = CompArr->naxis[1];
    modType = (ObitSkyModelCompType)(parms[3]+0.5);  /* model type */
    in->modType = MAX (in->modType, modType);  /* Need highest number */
 
    /* Gaussian parameters - assumed all the same */
    if ((modType==OBIT_SkyModel_GaussMod) || (modType==OBIT_SkyModel_GaussModSpec) || 
	(modType==OBIT_SkyModel_GaussModTSpec)) {
      cpa = cos (DG2RAD * parms[2]);
      spa = sin (DG2RAD * parms[2]);
      xmaj = parms[0] * konst2;
      xmin = parms[1] * konst2;
      gp1 = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
      gp2 = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
      gp3 = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
    }

    /* Does the CC table have a DeltaZ column? */
    if (CCTable->DeltaZCol>=0) toff = 4;
    else                       toff = 3;

    /* loop over CCs */
    for (j=0; j<larray; j++) {
 
     /* Only down to first negative? */
      if (noNeg && (array[0]<=0.0)) break;

      /* Do we want this one? */
      want = (fabs(array[0])>0.0);
      want = want && (fabs(array[0])>in->minFlux);
      want = want && (ncomp<count);  /* don't overflow */
      if (want) {

	/* Field number */
	table[0] = 1.0; /* Save for later use ; 1.0 to indicate valid data */
	table[1] = i;
	/* Beam offset for in->mosaic->images[i] */
	/* Converted to offsets from pointing in ObitSkyModelVMBeamInitModel */
        table[2] = array[1];
	table[3] = array[2];

	/* Geometric terms */
	xp[0] = (array[1] + xpoff) * konst;
	xp[1] = (array[2] + ypoff) * konst;
	if (CCTable->DeltaZCol>=0) xp[2] = array[3] * konst;
	else                       xp[2] = 0.0;
	if (do3Dmul) {
	  xyz[0] = xp[0]*umat[0][0] + xp[1]*umat[1][0];
	  xyz[1] = xp[0]*umat[0][1] + xp[1]*umat[1][1];
	  xyz[2] = xp[0]*umat[0][2] + xp[1]*umat[1][2]; 
	  /* PRJMUL (2, XP, UMAT, XYZ); */
	} else {  /* no rotation  */
	  xyz[0] = ccrot * xp[0] + ssrot * xp[1];
	  xyz[1] = ccrot * xp[1] - ssrot * xp[0];
	  xyz[2] = xp[2];
 	}
	table[4] = xyz[0] + xxoff;
	table[5] = xyz[1] + yyoff;
	table[6] = xyz[2] + zzoff;
	wfact = (1.0 + 2*(xyz[2] + zzoff)); /* Huh? */
	    
	/* Zero rest in case */
	for (iterm=7; iterm<lenEntry; iterm++) table[iterm] = 0.0;

	/* Copy subband spectrum with correction to table */
	indxfs = 7;             /* Index in table of beginning of fitted spectrum */
	indxsb = 7+in->nTerm;   /* Index in table of beginning of subband spectrum */
	if (in->nSpec>0) {
	  for (iterm=0; iterm<in->nSpec; iterm++) table[indxsb+iterm] = array[iterm+toff]*specCorr[iterm]*wfact;
	  
	  /* Fitting component spectra? */
	  if ((!noSIfit) && (in->nSpec>1)) {
	    /* Remove any zeroes, blanks from spectrum */
	    k = 0;
	    for (kk=0; kk<in->nSpec; kk++) {
	      if ((table[indxsb+kk]!=0.0) && (table[indxsb+kk]!=fblank)) {
		fitFreq[k] = in->specFreq[kk]; /* Frequency for fit */
		fitFlux[k] = table[indxsb+kk]; /* Subband flux density */
		k++;
	      }
	    } /* end of clean arrays */
	    /* Number of entries may vary, have to create fitArg each component */
	    fitArg = ObitSpectrumFitMakeArg (k, in->nTerm, in->refFreq, fitFreq, 
					     FALSE, &fitParms, err);
	    if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
	    ObitSpectrumFitSingleArg (fitArg, fitFlux, fitSigma, fitParms);
	    /* Crazy solution filter - doesn't seem to do that much good */
	    badFit = fitParms[0]<=0.0;  /* Bad flux density */
	    badFit = badFit || fitParms[0]<(10.*allRMS); /* must be at least 10 sigma */
	    for (kk=1; kk<in->nTerm; kk++) badFit |= fabsf(fitParms[0])>10.0;  /* |Terms| > 10 */
	    if (badFit) {  /* Flat spectrum at average */
	      /*fprintf (stderr,"Crazy soln %8.5f, %8.5f,%8.5f,%8.5f\n",fitParms[0], fitParms[1],fitParms[3],fitParms[3]);*/
	      fitParms[0] = ObitSkyModelVMBeamMFAvgSpec(in->nSpec, &table[indxsb]);
	      for (kk=1; kk<in->nTerm; kk++) fitParms[kk] = 0.0;
	    }
	    for (kk=0; kk<in->nTerm; kk++) table[indxfs+kk] = fitParms[kk]; /* Fitted spectrum */
	    if (fitArg)   {g_free(fitArg);   fitArg = NULL;}
	    if (fitParms) {g_free(fitParms); fitParms = NULL;}
	    /* end fit component spectrum */
	  } else if ((in->nTerm>0) && noSIfit) {
	  /* Dummy spectral Fit */
	    table[indxfs]   =  ObitSkyModelVMBeamMFAvgSpec(in->nSpec, &table[indxsb]);  /* Average so component doesn't get dropped */
	    table[indxfs+1] = 0.0;  /* If it gets here, this is likely polarization data */
	    for (kk=2; kk<in->nTerm; kk++) table[indxfs+kk] = 0.0;
	  } /* End of spectral fit */
	} /* End of have subband spectrum */

	/* Sum CC flux densities  */
	if (fitFlux) sumCCFlux += table[indxfs];
	else {  /* Use max abs value */
	  maxCCFlux = -1.0;
	  for (kk=0; kk<in->nSpec; kk++) {
	    if (table[indxsb+kk]!=fblank) maxCCFlux = MAX( maxCCFlux, fabsf(table[indxsb+kk]));
	  }
	  sumCCFlux += maxCCFlux;
	} /* end summing flux densities */

	/* Point - no spectrum  
	table[0] = array[0] * in->factor;*/

	/* Only same type as highest */
	if  (in->modType==modType) {
	  /* Only Point */
	  if (modType==OBIT_SkyModel_PointMod){
	    /* Nothing special this case */
	    
	    /* Only Point with tabulated spectrum */
	  } else if (modType==OBIT_SkyModel_PointModTSpec) {
	    /*for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+4] = array[iterm+toff]*specCorr[iterm]*wfact;*/
	    
	    /* Only Point with spectrum */
	  } else if (modType==OBIT_SkyModel_PointModSpec) {
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+indxfs] = array[iterm+toff];
	    table[indxfs]*= wfact;
	    
	    /* Only Gaussian */
	  } else if (in->modType==OBIT_SkyModel_GaussMod) {
	    table[nadd+0] = gp1; table[nadd+1] = gp2; table[nadd+2] = gp3;
	    
	    /* Only Gaussian + tabulated spectrum */
	  } else if (in->modType==OBIT_SkyModel_GaussModTSpec) {
	    table[7+nadd] = gp1; table[8+nadd] = gp2; table[9+nadd] = gp3;
	    /*  subband spectrum 
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+indxfs] = array[iterm+toff]*specCorr[iterm]*wfact; */
	    
	    /* Only Gaussian + spectrum */
	  } else if (in->modType==OBIT_SkyModel_GaussModSpec) {
	    table[7+nadd] = gp1; table[8+nadd] = gp2; table[9+nadd] = gp2;
	    /*  parameterized spectrum */
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+indxfs] = array[iterm+toff];
	    table[indxfs]*= wfact; 
	     
	    /* Only Uniform sphere */
	  } else if (in->modType==OBIT_SkyModel_USphereMod) {
	    table[nadd+0] = 3.0 * array[0] * in->factor;
	    table[nadd+1] = parms[1]  * 0.109662271 * 2.7777778e-4;
	    table[nadd+2] = 0.1;
	    
	    /* Only Uniform sphere + tabulated spectrum */
	  } else if (in->modType==OBIT_SkyModel_USphereModTSpec) {
	    table[nadd+0] = 3.0 * array[0] * in->factor;
	    table[nadd+1] = parms[1]  * 0.109662271 * 2.7777778e-4;
	    table[nadd+2] = 0.1;
	    /*   subband spectrum
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+indxfs] = array[iterm+toff]*specCorr[iterm]*wfact; */

	    /* Only Uniform sphere+ spectrum */
	  } else if (in->modType==OBIT_SkyModel_USphereModSpec) {
	    table[nadd+0] = 3.0 * array[0] * in->factor;
	    table[nadd+1] = parms[1]  * 0.109662271 * 2.7777778e-4;
	    table[nadd+2] = 0.1;
	    /*  parameterized spectrum */
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+indxfs] = array[iterm+toff];
	    table[indxfs]*= wfact;
	  }
	} else { /* Mixed type - zero unused model components */

	  /* Only Point here but some Gaussian - zero Gauss comps */
	  if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_GaussMod)) {
	    table[nadd+0] = 0.0; table[nadd+1] = 0.0; table[nadd+2] = 0.0;
	  
	    /* Gauss here but also some points */
	  } else if ((modType==OBIT_SkyModel_GaussMod) && (in->modType==OBIT_SkyModel_PointMod)) {
	    table[nadd+0] = gp1; table[nadd+1] = gp2; table[nadd+2] = gp3;
	  
	  /* Only PointSpectrum here but some GaussianSpectrum - zero Gauss comps */
	  } else if ((modType==OBIT_SkyModel_PointModTSpec) && (in->modType==OBIT_SkyModel_GaussModTSpec)) {
	    table[nadd+0] = 0.0; table[nadd+1] = 0.0; table[nadd+2] = 0.0;
	    /*for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+indxfs] = array[iterm+toff]*specCorr[iterm];*/
	  
	    /* GaussianTSpectrum here but also some PointTSpectrum */
	  } else if ((in->modType==OBIT_SkyModel_PointModTSpec) && (modType==OBIT_SkyModel_GaussModTSpec)) {
	    table[nadd+0] = gp1; table[nadd+1] = gp2; table[nadd+2] = gp3;
	    /*  subband spectrum 
		for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+indxfs] = array[iterm+toff]*specCorr[iterm]*wfact;*/
	  
	    /* Only Point here but some with spectrum - zero spectra (Unlikely) */
	  } else if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_PointModTSpec)) {
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+indxfs] = 0.0;
	    
	    /* Only PointSpectrum here but some GaussianSpectrum - zero Gauss comps */
	  } else if ((modType==OBIT_SkyModel_PointModSpec) && (in->modType==OBIT_SkyModel_GaussModSpec)) {
	    table[nadd+0] = 0.0; table[nadd+1] = 0.0; table[nadd+2] = 0.0;
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+indxfs] = array[iterm+toff];
	    table[indxfs] *= wfact; 
	    	    
	    /* GaussianSpectrum here but also some PointSpectrum */
	  } else if ((in->modType==OBIT_SkyModel_PointModSpec) && (modType==OBIT_SkyModel_GaussModSpec)) {
	    table[nadd+0] = gp1; table[nadd+1] = gp2; table[nadd+2] = gp3;
	    /*  spectrum */
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+indxfs] = array[iterm+toff];
	    table[indxfs] *= wfact;
	    
	    /* Only Point here but some with spectrum - zero spectra (Unlikely) */
	  } else if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_PointModSpec)) {
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+indxfs] = 0.0;
	    
	  } else { /* Unsupported combination */
	    Obit_log_error(err, OBIT_Error,"%s Unsupported combination of model types %d %d  %s",
			   routine, modType, in->modType, CCTable->name);
	    Obit_traceback_val (err, routine, in->name, retCode);
	  }
	} /* end mixed type */
       
#ifdef DEBUG
	/* Debug, list components */
	//fprintf(stderr,"\n Comp %d\n",ncomp);
	//for (kk=0; kk< lrec; kk+=5) {
	//  fprintf(stderr,"  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f \n",table[kk+0],table[kk+1],table[kk+2],table[kk+3],table[kk+4]);
	//}
#endif /* DEBUG */ 
	    
	/* Update */
	table += lrec;
	ncomp++;
      } /* End only desired */
      array += warray;
    } /* end loop over components */

    /* Save sumCCFlux */
    in->sumCCFlux = (ofloat)sumCCFlux;

    /* Delete merged CC array */
    CompArr = ObitFArrayUnref(CompArr);

    /* if outCCver>0 then the CCtable is temporary - Zap */
    if (outCCVer>0) {
      CCTable = ObitTableCCUnref (CCTable);
      ObitImageZapTable(in->mosaic->images[i], tabType, outCCVer, err);
    /* else simply release table  */
    } else CCTable = ObitTableCCUnref (CCTable);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    /* Release table row */
    CCRow = ObitTableCCRowUnref (CCRow);
    
  } /* end loop loading CCs */

  /* Zero any extra entries in table. */
  for (i=ncomp; i<count; i++) {
    /* Zero entry */
    table[0] = 0.0; table[1] = -10.0;  /* Mark no field */
    for (kk=2; kk<lrec; kk++) table[kk] = 0.0;
    table += lrec;  /* Update pointer */
  } /* end loop zeroing extra components */

  in->numComp = ncomp;
  /* Cleanup */
  if (specCorr) g_free(specCorr); 
  if (fitSigma) g_free(fitSigma); 
  if (fitFreq)  g_free(fitFreq); 
  if (fitFlux)  g_free(fitFlux ); 

  /* Find anything */
  gotSome = ncomp>0;

  return gotSome;
} /* end ObitSkyModelVMBeamMFLoadComps */

/**
 * NYI
 * Sets the in->plane member to either the pixels from the image in the 
 * specified field in in->mosaic or this array with relative 
 * primary beam corrections if in->doPBCor.
 * \param in       SkyModelVMBeamMF
 * \param uvdata   UV data
 * \param field    Field number in in->mosaic
 * \param err      Obit error stack object.
 * \return ObitCCTable to use, this should be Unref when done and 
 *   Zapped if outCCver != 0
 */
void ObitSkyModelVMBeamMFgetPBImage (ObitSkyModel* in, ObitUV* uvdata, olong field, 
				     ObitErr *err)
{
  g_error("ObitSkyModelVMBeamMFgetPBImage NOT implemented");
} /* end ObitSkyModelVMBeamMFgetPBImage */
  
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
void ObitSkyModelVMBeamMFFTDFT (ObitSkyModelVM *inn, olong field, ObitUV *uvdata, ObitErr *err)
{
  olong i, mcomp, iComp=0, pos[2], nvis, lovis, hivis, nvisPerThread, nThreads;
  olong kamp = 7;
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
  VMBeamMFFTFuncArg *args;
  ofloat *ddata,suma=0.0;
  //ofloat  fblank = ObitMagicF();
  gboolean OK = TRUE, oldUnitJones = FALSE;
  gchar *routine = "ObitSkyModelVMBeamMFFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* Check */
  args = (VMBeamMFFTFuncArg*)in->threadArgs[0];
  if ((strlen(args->type)>8) || (strncmp(args->type, "mfvmbeam", 8))) {
    Obit_log_error(err, OBIT_Error,"%s: Wrong type FuncArg %s", routine,args->type);
    return;
  }

  /* does sumCCFlux exceed Threshold? if so use Unit Jones matrices rather than Beams */
  oldUnitJones = in->doUnitJones;
  in->doUnitJones = (in->sumCCFlux<in->Threshold);
  /* Tell about it first time*/
  if (in->doUnitJones && !oldUnitJones) {
    Obit_log_error(err, OBIT_InfoErr, "Switching off beam corrections;");
    Obit_log_error(err, OBIT_InfoErr, "Sum component flux %8.5f < Threshold %8.5f",
		   in->sumCCFlux, in->Threshold);
    ObitErrLog(err);
  }

  /* Count number of actual components */
  mcomp = 0; suma = 0.0;
  pos[0] = 0; pos[1] = 0;
  ddata = ObitFArrayIndex(in->comps, pos);
  for (iComp=0; iComp<in->comps->naxis[1]; iComp++) {
    if (ddata[kamp]!=0.0) {mcomp = iComp+1; suma += ddata[kamp];}
    if (suma>100.) {fprintf (stderr,"DAMN, sum=%f, add %f index %d\n",suma,ddata[kamp], iComp); break;}
    ddata += in->comps->naxis[0];  /* update pointer */
  } /* end loop over components */
  in->numComp = MIN(in->numComp, mcomp);  /* Number of actual components */

  /* Divide up work */
  nvis = uvdata->myDesc->numVisBuff;
  if (nvis<1000) nThreads = 1;
  else nThreads = in->nThreads;

  // Progress report if in->prtLv>=5
  if (err->prtLv>=5){
    gchar timeStr[24];
    day2dhms(uvdata->buffer[uvdata->myDesc->iloct], timeStr); /* Human readable time */
    Obit_log_error(err, OBIT_InfoErr,"Time %s nvis %6d, ncomp %6d, sum Mod %f", timeStr, nvis,in->numComp,suma);
    ObitErrLog(err);
  }

  nvisPerThread = nvis/nThreads;
  lovis = 1;
  hivis = nvisPerThread;
  hivis = MIN (hivis, nvis);

  /* Set up thread arguments */
  for (i=0; i<nThreads; i++) {
    if (i==(nThreads-1)) hivis = nvis;  /* Make sure do all */
    args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
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
    ObitSkyModelVMBeamMFUpdateArg (in, args);  /* Redo arrays if needed */
  } /* end loop over threads */

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nThreads, in->DFTFunc, in->threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
}  /* end ObitSkyModelVMBeamMFFTDFT */

/* gcc or icc */
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((aligned(32)))
# define ALIGN64_BEG
# define ALIGN64_END __attribute__((aligned(64)))
/* AVX512 =16 float vectors */
#if HAVE_AVX512==1
#include <immintrin.h>
typedef __m512   v16sf;
typedef __m512d  v8df;
typedef __m256   v8sf;
typedef ALIGN64_BEG union {
  double    f[8];
  long long i[8];
  v8df      v;
} ALIGN64_END V8DF;
#endif
#if HAVE_AVX==1
#include <immintrin.h>
typedef __m256  v8sf;
typedef __m256d v4df;
typedef __m128  v4sf;
typedef ALIGN32_BEG union {
  double    f[4];
  long long i[4];
  v4df      v;
} ALIGN32_END V4DF;
#endif
/**
 * Routine to get dot product of two vectors 
 * possibly with AVX512 implementation
 * \param n   Number of elements
 * \param v1  First vector
 * \param v2  Second vector
 * \return dot product
 */
static inline ofloat SkyModelDot (olong n, ofloat *v1, ofloat *v2) {
  olong i, ilast;
  ofloat outval=0.0;
  /* AVX512=16 floats */
#if HAVE_AVX512==1
  v16sf vv1, vv2;
#endif

  if (n<=0) return outval;
  ilast = 0;
  /* AVX512 in blocks of 16 */
#if HAVE_AVX512==1
  for (i=ilast; i<n; i+=16) {
  if (i+16>n) break;
  ilast = i+16;
  vv1 = _mm512_load_ps((float*)&v1[i]);
  vv2 = _mm512_load_ps((float*)&v2[i]);
  vv1 = _mm512_mul_ps(vv1, vv2); 
  outval += (ofloat)_mm512_reduce_add_ps(vv1);
  }
#endif
  for (i=ilast; i<n; i++) {
    outval += v1[i]*v2[i];
  }
  return outval;
} /* end SkyModelDot */

/**
 * Routine to get 2D dot product of two vectors 
 * each element of ov is v1*s1 + v2*s2 + v3*s3, in double
 * possibly with AVX & AVX512 implementation
 * \param n   Number of elements
 * \param v1  First vector
 * \param v2  Second vector
 * \param v3  Third vector
 * \param s1  First scalar
 * \param s2  Second scalar
 * \param s3  Third scalar
 * \param ov  output vector
 */
static inline void SkyModel2DDot (olong n, ofloat *v1, ofloat *v2, ofloat *v3, 
			   ofloat s1, ofloat s2, ofloat s3,
			   ofloat *ov) {
  olong i, ilast;
  double ds1, ds2, ds3;
  /* AVX512=16 floats */
#if HAVE_AVX512==1
  v8df vv1, vv2, vv3, vs1, vs2, vs3;
  v8sf vsp;
#elif HAVE_AVX==1
  v4df vv1, vv2, vv3, vs1, vs2, vs3;
  v4sf vsp;
#endif

  if (n<=0) return;
  ds1 = s1; ds2 = s2; ds3 = s3;
  ilast = 0;
  /* AVX512 double in blocks of 8 */
#if HAVE_AVX512==1
  for (i=0; i<n; i+=8) { 
    if (ilast+8>n) break;
    ilast = i+8;
    vsp = _mm256_load_ps((float*)&v1[i]);
    vv1 = _mm512_cvtps_pd(vsp);  /* to double */
    vs1 = _mm512_set1_pd((double)s1);
    vv1 = _mm512_mul_pd(vv1, vs1); 
    vsp = _mm256_load_ps((float*)&v2[i]);
    vv2 = _mm512_cvtps_pd(vsp);  /* to double */
    vs2 = _mm512_set1_pd((double)s2);
    vv1 = _mm512_fmadd_pd(vv2, vs2, vv1); 
    vsp = _mm256_load_ps((float*)&v3[i]);
    vv3 = _mm512_cvtps_pd(vsp);  /* to double */
    vs3 = _mm512_set1_pd((double)s3);
    vv1 = _mm512_fmadd_pd(vv3, vs3, vv1); 
    vsp = _mm512_cvtpd_ps(vv1); /* to single precision */		       
    _mm256_store_ps((float*)&ov[i], vsp);
  }
  /* rest as scalar */
  /* AVX  double in blocks of 4 */
#elif HAVE_AVX==1
  for (i=0; i<n; i+=4) {
    if (ilast+4>n) break;
    ilast = i+4;
    vsp = _mm_load_ps((float*)&v1[i]);
    vv1 = _mm256_cvtps_pd(vsp);  /* to double */
    vs1 = _mm256_set1_pd((double)s1);
    vv1 = _mm256_mul_pd(vv1, vs1); 
    vsp = _mm_load_ps((float*)&v2[i]);
    vv2 = _mm256_cvtps_pd(vsp);  /* to double */
    vs2 = _mm256_set1_pd((double)s2);
    vv2 = _mm256_mul_pd(vv2, vs2); 
    vv1 = _mm256_add_pd(vv1, vv2); 
    vsp = _mm_load_ps((float*)&v3[i]);
    vv3 = _mm256_cvtps_pd(vsp);  /* to double */
    vs3 = _mm256_set1_pd((double)s3);
    vv3 = _mm256_mul_pd(vv3, vs3); 
    vv1 = _mm256_add_pd(vv1, vv3); 
    vsp = _mm256_cvtpd_ps(vv1); /* to single precision */		       
    _mm_store_ps((float*)&ov[i], vsp);
  }
  /* rest as scalar */
#endif
  for (i=ilast; i<n; i++) {
    ov[i] = (ofloat)(v1[i]*ds1 + v2[i]*ds2 + v3[i]*ds3);
  }
  return;
}/* end SkyModel2DDot */

/**
 * Routine to multiply vectors
 * possibly with AVX implementation
 * \param n   Number of elements
 * \param v1  First vector
 * \param v2  Second vector
 * \param ov  output vector
 */
 static inline void SkyModelVMul (olong n, ofloat *v1, ofloat *v2, ofloat *ov) {
  olong i, ilast;
#if HAVE_AVX512XX==1
  v16sf vv1, vv2;
#elif HAVE_AVX==1
  v8sf vv1, vv2;
#endif

  if (n<=0) return;
  ilast = 0;
  /* AVX512 in blocks of 16 - gives segv - disable */
#if HAVE_AVX512XX==1
  for (i=ilast; i<n; i+=16) {
  if (ilast+16>n) break;
  ilast = i+16;
  vv1 = _mm512_load_ps((float*)&v1[i]);
  vv2 = _mm512_load_ps((float*)&v2[i]);
  vv1 = _mm512_mul_ps(vv1, vv2); 
  _mm512_store_ps((float*)&ov[i], vv1);
  }
  /* rest as scalar */
  /* AVX in blocks of 8 */
#elif HAVE_AVX==1
  for (i=ilast; i<n; i+=8) {
  if (ilast+8>n) break;
  ilast = i+8;
  vv1 = _mm256_load_ps((float*)&v1[i]);
  vv2 = _mm256_load_ps((float*)&v2[i]);
  vv1 = _mm256_mul_ps(vv1, vv2); 
  _mm256_store_ps((float*)&ov[i], vv1);
  }
  /* rest as scalar */
#endif
  for (i=ilast; i<n; i++) {
    ov[i] = v1[i]*v2[i];
  }
  return;
 } /* end SkyModelVMul */

/**
 * Calculate visibility model and sum over a set of components
 * \param sumMod   in/out full Stokes visibility accumulator
 * \param ncomp    number of components in AmpArr, SinArr, CosArr
 * \param isCirc   TRUE if circular feeds
 * \param polType  Model Stoke type, 1,2,3,4 => I,Q.U.V
 * \param AmpArr   Array of component amplitudes
 * \param SinArr   Array of component sine of phase
 * \param CosArr   Array of component cosine of phase
 * \param cos2PA   cosine twice parallactic angle
 * \param sin2PA   sine twice parallactic angle
 * \param beamArr1 Antenna 1 beam
 * \param beamArr2 Antenna 2 beam
 * \param beamCT   work matrix
 * \param modl     work matrix
 * \param matxt1   work matrix
 * \param matxt2   work matrix
 */
void calcMod(ObitMatx *sumMod, olong ncomp, gboolean isCirc, olong polType, 
	     ofloat *AmpArr, ofloat *SinArr, ofloat *CosArr, 
	     ofloat cos2PA, ofloat sin2PA, ObitMatx *beamArr1[], ObitMatx *beamArr2[], 
	     ObitMatx *beamCT, ObitMatx *modl, ObitMatx *matxt1, ObitMatx *matxt2) {
  olong jt;
  /* Accumulate model in sumMod */
  /* Use matrix formalism modl = beam1*model*beam2.conjugate_transpose */
  for (jt=0; jt<ncomp; jt++) {
    if (isCirc) {  /* Circular Feeds */
      /* from TMS2 p 109, eq 4.46 for XX,XY,YX,YY */
      switch (polType) {
      case 1:     /* Stokes I */
	ObitMatxSet2C(modl, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt], 0.0, 0.0,
		      0.0, 0.0, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt]);
	break;
      case 2:     /* Stokes Q */
	ObitMatxSet2C(modl, 0.0,0.0, 
		      -AmpArr[jt]*SinArr[jt], -AmpArr[jt]*CosArr[jt],
		      +AmpArr[jt]*SinArr[jt], -AmpArr[jt]*CosArr[jt],
		      0.0,0.0);
	break;
      case 3:     /* Stokes U */
	ObitMatxSet2C(modl, 0.0,0.0, 
		      +AmpArr[jt]*CosArr[jt], +AmpArr[jt]*SinArr[jt],
		      -AmpArr[jt]*CosArr[jt], -AmpArr[jt]*SinArr[jt],
		      0.0,0.0);
	break;
      case 4:     /* Stokes V */
	ObitMatxSet2C(modl, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt], 0.0, 0.0,
		      0.0, 0.0, -AmpArr[jt]*CosArr[jt], -AmpArr[jt]*SinArr[jt]);
	break;
      default:   /* Stokes I */
	ObitMatxSet2C(modl, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt], 0.0, 0.0,
		      0.0, 0.0, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt]);
      }; /* end polType switch */
    } else {       /* Linear feeds */
      /* from TMS2 p 107, eq 4.41 for XX,XY,YX,YY */
      switch (polType) {
      case 1:     /* Stokes I */
	ObitMatxSet2C(modl, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt], 0.0, 0.0,
		      0.0, 0.0, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt]);
	break;
      case 2:     /* Stokes Q */
	ObitMatxSet2C(modl, 
		      +cos2PA*AmpArr[jt]*CosArr[jt], +cos2PA*AmpArr[jt]*SinArr[jt],
		      -sin2PA*AmpArr[jt]*CosArr[jt], -sin2PA*AmpArr[jt]*SinArr[jt],
		      -sin2PA*AmpArr[jt]*CosArr[jt], -sin2PA*AmpArr[jt]*SinArr[jt],
		      -cos2PA*AmpArr[jt]*CosArr[jt], -cos2PA*AmpArr[jt]*SinArr[jt]);
	break;
      case 3:     /* Stokes U */
	ObitMatxSet2C(modl, 
		      +sin2PA*AmpArr[jt]*CosArr[jt], +sin2PA*AmpArr[jt]*SinArr[jt], 
		      +cos2PA*AmpArr[jt]*CosArr[jt], +cos2PA*AmpArr[jt]*SinArr[jt],
		      +cos2PA*AmpArr[jt]*CosArr[jt], +cos2PA*AmpArr[jt]*SinArr[jt],
		      -sin2PA*AmpArr[jt]*CosArr[jt], -sin2PA*AmpArr[jt]*SinArr[jt]);
	break;
      case 4:     /* Stokes V */
	ObitMatxSet2C(modl, 0.0, 0.0, -AmpArr[jt]*SinArr[jt], +AmpArr[jt]*CosArr[jt], 
		      +AmpArr[jt]*SinArr[jt], -AmpArr[jt]*CosArr[jt], 0.0, 0.0);
	break;
      default:   /* Stokes I */
	ObitMatxSet2C(modl, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt], 0.0, 0.0,
		      0.0, 0.0, AmpArr[jt]*CosArr[jt], AmpArr[jt]*SinArr[jt]);
      }; /* end polType switch */
    } /* end circ/lin feeds */
    ObitMatxCTrans(beamArr2[jt], beamCT);      /* Conjugate transpose of Beam2 */
    ObitMatxMult(modl, beamCT, matxt1);        /* model*beam2.conjugate_transpose */
    ObitMatxMult(beamArr1[jt], matxt1, matxt2);/* beam1*model*beam2.conjugate_transpose */
    ObitMatxAdd(matxt2, sumMod, sumMod);
  } /* end component loop */
} /* end calcMod */

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
 * Arguments are given in the VMBeamMFFTFuncArg structure passed as arg starting 
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
 * \li channel Current UV channel being processed (used in model update ).
 * \return NULL
 */
static gpointer ThreadSkyModelVMBeamMFFTDFT (gpointer args)
{
  /* Get arguments from structure */
  VMBeamMFFTFuncArg *largs = (VMBeamMFFTFuncArg*)args;
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)largs->in;
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
  ofloat *chFlux      = largs->chFlux;
 
  olong iVis=0, iIF, ifq, iChannel, iStoke, iComp=0, lcomp;
  olong lrec, nrparm, naxis[2], iaty, jaty;
  olong jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF, lim;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVIF, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, suba, it1, it2, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealRL,  sumImagRL,  modRealRL=0.0,  modImagRL=0.0;
  ofloat sumRealLR,  sumImagLR,  modRealLR=0.0,  modImagLR=0.0;
  ofloat ll, lll, logNuONu0;
  ofloat arg, freq2=0.0,freqFact, wtRR=0.0, wtLL=0.0, temp;
  ofloat  fblank = ObitMagicF();
#define FazArrSize 256  /* Size of the amp/phase/sine/cosine arrays */
#if HAVE_AVX512==1
  __attribute__((aligned(32))) ofloat AmpArr[FazArrSize], FazArr[FazArrSize], CosArr[FazArrSize], SinArr[FazArrSize];
  __attribute__((aligned(32))) ofloat ExpArg[FazArrSize], ExpVal[FazArrSize], ExpArg2[FazArrSize], ExpVal2[FazArrSize];
  __attribute__((aligned(32))) ofloat VecL[FazArrSize], VecM[FazArrSize], VecN[FazArrSize];
 #elif HAVE_AVX==1
  __attribute__((aligned(32))) ofloat AmpArr[FazArrSize], FazArr[FazArrSize], CosArr[FazArrSize], SinArr[FazArrSize];
  __attribute__((aligned(32))) ofloat ExpArg[FazArrSize], ExpVal[FazArrSize], ExpArg2[FazArrSize], ExpVal2[FazArrSize];
  __attribute__((aligned(32))) ofloat VecL[FazArrSize], VecM[FazArrSize], VecN[FazArrSize];
#else 
  ofloat AmpArr[FazArrSize], FazArr[FazArrSize], CosArr[FazArrSize], SinArr[FazArrSize];
  ofloat ExpArg[FazArrSize], ExpVal[FazArrSize], ExpArg2[FazArrSize], ExpVal2[FazArrSize];
  ofloat VecL[FazArrSize], VecM[FazArrSize], VecN[FazArrSize];
#endif
  olong it, kt, jt, itcnt, iSpec, nSpec, itab, Stokes;
  /* Offsets in ddata table */
  /*
    assignments
   0                      1.0
   1                      Field number, 0-rel
   2                      CC DeltaX for beam correction
   3                      CC DeltaY
   4->6                   l,m,n terms
   7->6+nTerm             nTerm Spectral terms: flux, si, curv1, curv2...
   7+nTerm->6+nTerm+nSpec Subband spectra, nSpec entries
   7+nTerm+nSpec->        other terms
       gaussian, uniform sphere
       ??? parameterized spectra use evaluator
   */
  olong kamp=0, ku=4, kv=5, kw=6, kfit=7, ksi=8,  kadd, ka1=0,ka2=1,ka3=2;
  gboolean doCrossPol, isCirc, doSpEval=FALSE, allBad=FALSE;
  odouble *freqArr, SMRefFreq, specFreqFact, Freq;
  odouble u, v, w;
  const ObitSkyModelVMClassInfo *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
#ifdef DEBUG
  ofloat gain, tre, tim;  /* Debugging values */
  olong imod;
#endif /* DEBUG */ 

  gchar *routine = "ThreadObitSkyModelVMBeamMFFTDFT";

         //DEBUG
  //fprintf(stderr,"in Thread %d t=%f\n",ithread,  24*uvdata->buffer[loVis*lrec+uvdata->myDesc->iloct]);
        // END DEBUG


  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  if ((largs->ithread<=0) && (uvdata->myDesc->firstVis<10)) {
    if (in->doUnitJones)
      Obit_log_error(err, OBIT_InfoErr, "SkyModel with %d components without Beam Cor", in->numComp);
    else
      Obit_log_error(err, OBIT_InfoErr, "SkyModel with %d components with Beam Cor", in->numComp);
    ObitErrLog(err);
  }

  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  nSpec         = in->nSpec - 1;  /* Offset in comp array using TSpec */
  /* PB corr should be mostly turned off higher up */
  startIF       = in->startIF-1;
  numberIF      = MAX (1, in->numberIF);
  jincif        = uvdata->myDesc->incif;
  if (uvdata->myDesc->jlocif>=0) nUVIF = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  else                           nUVIF = 1;
  startChannel  = in->startChannel-1;
  numberChannel = MAX (1, in->numberChannel);
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

  isCirc = uvdata->myDesc->crval[uvdata->myDesc->jlocs]==-1.0;  /* Circular or linear feeds? */
  /* Get stokes from image, else I */
  if (in->mosaic) Stokes = (olong)(in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocs]+0.5);
  else Stokes = 1; 
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
  kadd  = 7+in->nSpec+in->nTerm;
  
  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;
  /* Inverse of Sky model reference frequency */
  if (in->mosaic!=NULL) 
    SMRefFreq = 1.0 / in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocf];
  else 
    SMRefFreq = 1.0 / uvdata->myDesc->freq;

  lrec    = uvdata->myDesc->lrec;         /* Length of record */
  visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
  nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */

  /* Starting parameters this pass */
  lstartIF       = startIF;
  lstartChannel  = startChannel;
  iVis     = 0;
  
  /**************************** DEBUG **********************************/
#ifdef DEBUG0 
  /* Replace model with single component and set to replace the data with the model */
  in->doReplace = TRUE;
  mcomp         = 1;      /* One component */
  in->modType   = OBIT_SkyModel_GaussModTSpec;  /* Model type */
  ddata = Data;  /* replace model */
  ddata[0] = ddata[0];  /* Field number */
  ddata[1] = -0.124923;  /* X offset from center (deg) */
  ddata[2] = 0.0;  /* Y offset from center (deg) */
  amp = 1.0;    /* Amplitude all channels */
  for (iSpec=0; iSpec<=nSpec; iSpec++) ddata[kamp + iSpec] = amp;
  /* Position - no rotation - no field offset */
  ddata[4+nSpec] = ddata[1] * DG2RAD * 2.0 * G_PI;
  ddata[5+nSpec] = ddata[2] * DG2RAD * 2.0 * G_PI;
  ddata[6+nSpec] = 0.0;
  /* Gaussian - only round 0.2" */
  ddata[7+nSpec] = (0.20/3600.0) * DG2RAD * (G_PI / 1.17741022) * sqrt (0.5);
  ddata[7+nSpec] = -ddata[7+nSpec] * ddata[7+nSpec];
  ddata[8+nSpec] = ddata[7+nSpec];
  ddata[9+nSpec] = 0.0;
#endif /* DEBUG0 */ 
  /************************ end DEBUG **********************************/
  
  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

  /* Update spectral evaluator - only Stokes I */
  doSpEval = (largs->spEval!=NULL) && (in->nTerm>0) && (Stokes==1);  /* Have an interpolator? spectrum? */
  if (doSpEval) ObitSpectrumMFUpdate (largs->spEval, in->nSpec, NULL, in->refFreq, in->nTerm);

  /* Loop over IFs */
  for (iIF=lstartIF; iIF<lstartIF+numberIF; iIF++) {
    offsetIF = nrparm + iIF*jincif; 
    /* Loop over channels */
    for (iChannel=lstartChannel; iChannel<lstartChannel+numberChannel; iChannel++) {
      offsetChannel = offsetIF + iChannel*jincf; 
      ifq      = MIN (nUVchan*nUVIF, MAX (0, iIF*kincif + iChannel*kincf));
      iSpec    = MIN (nSpec, MAX (0,in->specIndex[ifq]));
      freqFact = fscale[ifq];  /* Frequency scaling factor */
      freq2    = freqFact*freqFact;    /* Frequency factor squared */
      specFreqFact = freqArr[ifq] * SMRefFreq;

      /* Setup to evaluate MF spectrum if in->spEval */
      itab      = in->specIndex[ifq];
      /* Log ratio of channel freq to Tabulated (subband) freq */
      logNuONu0 = (ofloat)log(freqArr[ifq]/in->specFreq[itab]);
      Freq = freqArr[ifq];  /* Channel freq */
      if (doSpEval)  {
	for (it=0; it<mcomp; it++) { /* interpolate component values for this channel */
	  ddata = Data + it*lcomp;
	  chFlux[it] = ObitSpectrumMFEval(largs->spEval, Freq, iSpec, &ddata[kfit+in->nTerm], &ddata[kfit]);
	}/* End loop over components */
      } else {
	/* Hack to trap single component passed as an argument */
	if ((in->modType==OBIT_SkyModel_PointMod) && (in->numComp=1) && (in->mosaic==NULL)) {
	  kamp = 0; ku = 1; kv = 2; kw = 3;
	}
      } /* end if (not) using evaluator */
	  
      if ((ifq==0) || (in->FreqPlane[ifq]!=in->FreqPlane[ifq-1]))
	largs->endVMModelTime = -1.0e20;  /* Force Reset gains for new beam channel */
      largs->channel  = ifq;
      largs->BeamFreq = freqArr[ifq];  /* need beam frequency at Beam channel */
      
      /* Loop over vis  */
      for (iVis=loVis; iVis<hiVis; iVis++) {
	visData = uvdata->buffer+iVis*lrec;    /* Buffer pointer with appropriate offset */
	
	//DEBUG
	//if ((ithread==1) && (iIF==6) && (iChannel==20)) {
	//if ((ithread==1) && (iIF==7) && (iVis==loVis+2) && (iChannel>230)) {
	//if ((ithread==0) && (iIF==0) && (iVis==loVis+2) && (iChannel==0)) {
	//  fprintf(stderr,"thread %d vis %d t=%f,iIF %d, iChannel %d doSpEval %d mcomp %d\n",ithread, iVis, 24*visData[iloct],iIF, iChannel,doSpEval,mcomp);
	//} // END DEBUG
	//if ((ithread==23) && (iIF==7) && (iVis==hiVis-2) && (iChannel==236)) {
	//  fprintf(stderr,"thread %d vis %d t=%f,iIF %d, iChannel %d doSpEval %d mcomp %d\n",ithread, iVis, 24*visData[iloct],iIF, iChannel,doSpEval,mcomp);
	//} // END DEBUG
	
	/* Exceed current time validity for gains? */
	if (visData[iloct] > largs->endVMModelTime) {
	  /* Subarray 0-rel */
	  ObitUVDescGetAnts(uvdata->myDesc, visData, &it1, &it2, &suba);
	  /* Update */
	  myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba-1, uvdata, ithread, err);
	  if (err->error) {
	    ObitThreadLock(in->thread);  /* Lock against other threads */
	    Obit_log_error(err, OBIT_Error,"%s Error updating Model",
			   routine);
	    ObitThreadUnlock(in->thread); 
	    goto finish;
	  }
	} /* end update */
	
	//*If all data flagged and not replacing, skip */
	offset = offsetChannel; /* Visibility offset */
	allBad = (visData[offset+2]<=0.0) && (visData[offset+5]<=0.0);
	if (doCrossPol) { // Also XY,YX
	  allBad = allBad && (visData[offset+8]<=0.0) && (visData[offset+11]<=0.0);
	}
	if (allBad && (!in->doReplace)) continue;
	
	/* Need antennas numbers */
	ObitUVDescGetAnts(uvdata->myDesc, visData, &ant1, &ant2, &it1);
	ant1--;    /* 0 rel */
	ant2--;    /* 0 rel */
	iaty = in->AntType[ant1];  /* Antenna type */
	jaty = in->AntType[ant2];
	
	ifq = iIF*kincif + iChannel*kincf;
	freqFact = fscale[ifq];  /* Frequency scaling factor */
	itab      = in->specIndex[ifq];

	/* u,v,w at frequency */
	u = (odouble)visData[ilocu]*freqFact;
	v = (odouble)visData[ilocv]*freqFact;
	w = (odouble)visData[ilocw]*freqFact;

	/* Sum over components */
	sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	sumRealRL = sumImagRL = sumRealLR = sumImagLR = 0.0;
	ddata = Data;
	
	/* Sum by model type - assume phase same for RR, LL */
	kt = 0;   /* Jones matrix index */
	switch (in->modType) {
	case OBIT_SkyModel_PointMod:     /* Point */
	  /***************************************************************
	                Point                                       
	  */
	  /* From the AIPSish QXXPTS.FOR  */
	  jt = 0;  /* index in chFlux */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    /* Separate loops for with spEval and without */
	    if (doSpEval)  {   /* Evaluated spectrum? */
	      for (iComp=it; iComp<mcomp; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if ((chFlux[jt]!=0.0) && (chFlux[jt]!=fblank)) {  /* valid? */
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  AmpArr[itcnt] = chFlux[jt];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		jt++;
		if (itcnt>=FazArrSize) break;
	      } /* end inner loop over components */
	    } else { /* without spEval */
	      for (iComp=it; iComp<mcomp; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[kamp]!=0.0) {
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  AmpArr[itcnt] = ddata[kamp];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
		if (itcnt>=FazArrSize) break;
	      } /* end inner loop over components */
	    } /* end w/o SpEval */
	    
	    /* Compute phases */
	    SkyModel2DDot(itcnt, VecL, VecM, VecN, u, v, w, FazArr);
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Accumulate real and imaginary parts to sumVis */
	    ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					  AmpArr, SinArr, CosArr,
					  &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					  largs->cos2PA, largs->sin2PA, 
					  workVis, work1, work2, sumVis);
	    
	    sumRealRR += sumVis->data.flt[0]; sumImagRR += sumVis->data.flt[1];
	    sumRealLL += sumVis->data.flt[6]; sumImagLL += sumVis->data.flt[7];
	    if (doCrossPol) {
	      /* Cross pol */
	      sumRealRL += sumVis->data.flt[2]; sumImagRL += sumVis->data.flt[3];
	      sumRealLR += sumVis->data.flt[4]; sumImagLR += sumVis->data.flt[5];
	    } /* end xpol */
	    /************************ DEBUG **********************************/
#ifdef DEBUG
	    /* DEBUG t=0.00/08:19:20.5 = 0.346765. within 5 sec), IF=7, Chan=40, bl 40-50 
	     2nd scan on J0408-6455  */
	    /* DEBUG J0402-6509  t=00/08:28:20.3 =0.3530127 in 2nd scan, IF=7, Chan=40, bl 40-50 QPol */
	    if (0 &&  (iIF==6) && (iChannel==39) && (ant1==39) && (ant2==49) && (kt==0) && (Stokes==2)) {
	      /* Desired time? */
	      //if (fabsf (visData[iloct]-0.352550)<5.787037037037037e-05) {
	      //if (fabsf (visData[iloct]-0.353012)<5.787037037037037e-05) {
	      //if (fabsf (visData[iloct]-0.346765)<5.787037037037037e-05) {
	      if (fabsf (visData[iloct]-0.3781875)<5.78703e-05) { // J0400-6544 0/9:4:35.4
		gain = largs->beamGain;
		fprintf (stderr,"t=%f, IF=%d, chan=%d bl=%d-%d iVis=%d\n",
			 visData[iloct],iIF+1, iChannel+1, ant1+1, ant2+1, iVis);
		fprintf (stderr,"mod 1, %8.5f, %8.5f, %8.5f, %8.5f,\n",AmpArr[0], VecL[0], VecM[0], VecN[0]);
		fprintf (stderr,"beamX %8.5f beamY %8.5f beamPA %8.5f beamPlane %d gain %6.4f\n",
			 largs->beamX,largs->beamY,largs->beamPA,largs->beamPlane, largs->beamGain); 
		fprintf (stderr,"Jones %8.5f, %8.5f,  %8.5f, %8.5f,\n", 
			 JMatrix[iaty][0]->data.cpx[0].real, JMatrix[iaty][0]->data.cpx[0].imag, 
			 JMatrix[iaty][0]->data.cpx[1].real, JMatrix[iaty][0]->data.cpx[1].imag);
		fprintf (stderr,"      %8.5f, %8.5f,  %8.5f, %8.5f,\n", 
			 JMatrix[iaty][0]->data.cpx[2].real, JMatrix[iaty][0]->data.cpx[2].imag, 
			 JMatrix[iaty][0]->data.cpx[3].real, JMatrix[iaty][0]->data.cpx[3].imag);
		offset = offsetChannel; /* Visibility offset */
		fprintf (stderr,"XX: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 visData[offset+0],visData[offset+1],visData[offset+2],
			 sqrtf(visData[offset+0]*visData[offset+0]+visData[offset+1]*visData[offset+1]),
			 57.296*atan2f(visData[offset+1],visData[offset+0]));
		fprintf (stderr,"YY: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 visData[offset+3],visData[offset+4],visData[offset+5],
			 sqrtf(visData[offset+3]*visData[offset+3]+visData[offset+4]*visData[offset+4]),
			 57.296*atan2f(visData[offset+4],visData[offset+3]));
		fprintf (stderr,"XY: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 visData[offset+6],visData[offset+7],visData[offset+8],
			 sqrtf(visData[offset+6]*visData[offset+6]+visData[offset+7]*visData[offset+7]),
			 57.296*atan2f(visData[offset+7],visData[offset+6]));
		fprintf (stderr,"YX: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 visData[offset+9],visData[offset+10],visData[offset+11],
			 sqrtf(visData[offset+9]*visData[offset+9]+visData[offset+10]*visData[offset+10]),
			 57.296*atan2f(visData[offset+10],visData[offset+9]));
		fprintf (stderr,"ModXX: %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 sumRealRR, sumImagRR,
			 sqrtf(sumRealRR*sumRealRR+sumImagRR*sumImagRR),
			 57.296*atan2f(sumImagRR,sumRealRR));
		fprintf (stderr,"ModYY: %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 sumRealLL, sumImagLL,
			 sqrtf(sumRealLL*sumRealLL+sumImagLL*sumImagLL),
			 57.296*atan2f(sumImagLL,sumRealLL));
		fprintf (stderr,"ModXY: %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 sumRealRL, sumImagRL,
			 sqrtf(sumRealRL*sumRealRL+sumImagRL*sumImagRL),
			 57.296*atan2f(sumImagRL,sumRealRL));
		fprintf (stderr,"ModYX: %8.5f, %8.5f, polar %8.5f %8.1f\n",
			 sumRealLR, sumImagLR,
			 sqrtf(sumRealLR*sumRealLR+sumImagLR*sumImagLR),
			 57.296*atan2f(sumImagLR,sumRealLR));
	      }}
#endif /* DEBUG */ 
	    /************************ end DEBUG **********************************/
	    
	    kt = it;  /* offset in JMatrix */
	    } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_PointModSpec:     /* Point + spectrum */
	  /********************************************************************
	                Point Spec                                      
	  */
	  itab = 7;  /* Start of parameterized spectrum */
	  jt = 0;  /* index in chFlux */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    /* Separate loops for with spEval and without */
	    if (doSpEval)  {   /* Evaluated spectrum? */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if ((chFlux[jt]!=0.0) && (chFlux[jt]!=fblank)){  /* valid? */
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  AmpArr[itcnt] = chFlux[jt];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
		jt++;
	      } /* end inner loop over components */
	    } else { /*  w/o SpEval */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[kfit]!=0.0) {  /* valid? */
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  /* Frequency dependent terms */
		  lll = ll = log(specFreqFact);
		  arg = 0.0;
		  for (iterm=1; iterm<nterm; iterm++) {
		    arg += ddata[kfit+iterm] * lll;
		    lll *= ll;
		  }
		  ExpArg2[itcnt] = arg;
		  AmpArr[itcnt]  = ddata[kfit];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
	      } /* end inner loop over components */
	    } 	/* end w/o  SpEval */
	    
	    /* Compute phases */
	    SkyModel2DDot(itcnt, VecL, VecM, VecN, u, v, w, FazArr);
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    if (!doSpEval) {
	      /* Evaluate spectrum */
	      ObitExpVec(itcnt, ExpArg2, ExpVal2);
	      SkyModelVMul(itcnt, AmpArr, ExpVal2, AmpArr);
	    }

	    /* Accumulate real and imaginary parts to sumVis */
	    ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					  AmpArr, SinArr, CosArr,
					  &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					  largs->cos2PA, largs->sin2PA, 
					  workVis, work1, work2, sumVis);
	    
	    sumRealRR += sumVis->data.flt[0]; sumImagRR += sumVis->data.flt[1];
	    sumRealLL += sumVis->data.flt[6]; sumImagLL += sumVis->data.flt[7];
	    if (doCrossPol) {
	      /* Cross pol */
	      sumRealRL += sumVis->data.flt[2]; sumImagRL += sumVis->data.flt[3];
	      sumRealLR += sumVis->data.flt[4]; sumImagLR += sumVis->data.flt[5];
	    } /* end xpol */
	    kt = it;  /* offset in JMatrix */
	  } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_PointModTSpec:     /* Point + tabulated spectrum */
	  /************************************************************
	                Point  TSpec                                     
	  */
	  itab = 7 + in->nTerm + iSpec;       /* Closest subband index */
	  jt = 0;  /* index in chFlux */
	  for (it=0; it<mcomp; it+=FazArrSize) {  /* Outer component loop */
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    /* Separate loops for with spEval and without */
	    if (doSpEval)  {   /* Evaluated spectrum? */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if ((chFlux[jt]!=0.0) && (chFlux[jt]!=fblank)) {  /* valid? */
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  AmpArr[itcnt] = chFlux[jt];
		  itcnt++;    /* Count valid in amp/phase buffers */
		}  /* end if valid */
	      jt++;
	      } /* end inner loop over components */
	    } else { /*  without spEval */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[itab]!=0.0) {  /* valid? */
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  /* Amplitude  - use spectral index */
		  AmpArr[itcnt] = ddata[itab];
		  ExpArg2[itcnt] = logNuONu0 * ddata[ksi];
		  itcnt++;    /* Count valid in amp/phase buffers */
		}  /* end if valid */
	      } /* end inner loop over components */
	  } /* end w/o SpEval */
	    
	    /* Compute phases */
	    SkyModel2DDot(itcnt, VecL, VecM, VecN, u, v, w, FazArr);
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Evaluate spectral index? */
	    if (!doSpEval)  {
	      ObitExpVec(itcnt, ExpArg2, ExpVal2);
	      /* Spectral index correction */
	      SkyModelVMul(itcnt, AmpArr,  ExpVal2, AmpArr);
	    }
	    /************************ DEBUG **********************************/
#ifdef DEBUG
	    if (1 &&  (iIF==0) && (iChannel==0) && (ant1==38) && (ant2==47) && (kt==0)) { // 
	      /* Desired time?*/
	      if (fabsf (visData[iloct]-0.27495259046)<2*5.78703e-05) { // 0/06:15:38.416
	    	fprintf (stderr,"DEBUG stop\n");
	      }}  // END DEBUG
#endif /* DEBUG */ 
	    /************************ end DEBUG **********************************/

	    /* Accumulate real and imaginary parts to sumVis */
	    ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					  AmpArr, SinArr, CosArr,
					  &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					  largs->cos2PA, largs->sin2PA, 
					  workVis, work1, work2, sumVis);
	    
	    sumRealRR += sumVis->data.array[0]; sumImagRR += sumVis->data.array[1];
	    sumRealLL += sumVis->data.array[6]; sumImagLL += sumVis->data.array[7];
	    if (doCrossPol) {
	      /* Cross pol */
	      sumRealRL += sumVis->data.array[2]; sumImagRL += sumVis->data.array[3];
	      sumRealLR += sumVis->data.array[4]; sumImagLR += sumVis->data.array[5];
	    } /* end xpol */
	    
	    /************************ DEBUG **********************************/
#ifdef DEBUG
	    /* 1st channel, bl 39-48 (1-rel)00/06:35:55.9  */
	    if (1 &&  (iIF==0) && (iChannel==0) && (ant1==38) && (ant2==47) && (kt==0)) { // 
	      /* Desired time?*/
	      if (fabsf (visData[iloct]-0.2749525904)<2*5.78703e-05) { // 0/06:15:38.416
		gain = largs->beamGain;
		fprintf (stderr,"t=%f, IF=%d, chan=%d bl=%d-%d iVis=%d Freq=%8.5f\n",
			 visData[iloct],iIF+1, iChannel+1, ant1+1, ant2+1, iVis, Freq);
		fprintf (stderr,"beamX %10.7f beamY %10.7f beamPA %10.7f beamPlane %d gain %8.5f\n",
			 largs->beamX,largs->beamY,largs->beamPA,largs->beamPlane, gain);
		gain = sqrtf(gain);
		/* Loop over model components */
		for (imod=0; imod<itcnt; imod++) {
		  if (itcnt>=1) fprintf (stderr,"\n mod %d, %10.7f, %10.7f, %10.7f, %10.7f,rot %8.5f %8.5f\n",
					 imod, AmpArr[imod], VecL[imod], VecM[imod], VecN[imod], CosArr[imod], SinArr[imod]);
		  if (imod<10) fprintf (stderr,"gain = %8.5f\n",largs->beamWork[imod]);
		  fprintf (stderr,"Jones1 %10.7f, %10.7f,  %10.7f, %10.7f,\n ",
			   JMatrix[iaty][imod]->data.cpx[0].real, JMatrix[iaty][imod]->data.cpx[0].imag, 
			   JMatrix[iaty][imod]->data.cpx[1].real, JMatrix[iaty][imod]->data.cpx[1].imag);
		  fprintf (stderr,"      %10.7f, %10.7f,  %10.7f, %10.7f,\n", 
			   JMatrix[iaty][imod]->data.cpx[2].real, JMatrix[iaty][imod]->data.cpx[2].imag, 
			   JMatrix[iaty][imod]->data.cpx[3].real, JMatrix[iaty][imod]->data.cpx[3].imag);
		  fprintf (stderr,"Jones1 raw %10.7f, %10.7f,  %10.7f, %10.7f,\n ", 
			   JMatrix[iaty][imod]->data.cpx[0].real*gain, JMatrix[iaty][0]->data.cpx[0].imag*gain, 
			   JMatrix[iaty][imod]->data.cpx[1].real*gain, JMatrix[iaty][0]->data.cpx[1].imag*gain);
		  fprintf (stderr,"          %10.7f, %10.7f,  %10.7f, %10.7f,\n", 
			   JMatrix[iaty][imod]->data.cpx[2].real*gain, JMatrix[iaty][imod]->data.cpx[2].imag*gain, 
			   JMatrix[iaty][imod]->data.cpx[3].real*gain, JMatrix[iaty][imod]->data.cpx[3].imag*gain);
		  fprintf (stderr,"Jones2 %10.7f, %10.7f,  %10.7f, %10.7f,\n ",
			   JMatrix[jaty][imod]->data.cpx[0].real, JMatrix[jaty][imod]->data.cpx[0].imag, 
			   JMatrix[jaty][imod]->data.cpx[1].real, JMatrix[jaty][imod]->data.cpx[1].imag);
		  fprintf (stderr,"      %10.7f, %10.7f,  %10.7f, %10.7f,\n", 
			   JMatrix[jaty][imod]->data.cpx[2].real, JMatrix[jaty][imod]->data.cpx[2].imag, 
			   JMatrix[jaty][imod]->data.cpx[3].real, JMatrix[jaty][imod]->data.cpx[3].imag);
		  fprintf (stderr,"Jones2 raw %10.7f, %10.7f,  %10.7f, %10.7f,\n ", 
			   JMatrix[jaty][imod]->data.cpx[0].real*gain, JMatrix[jaty][0]->data.cpx[0].imag*gain, 
			   JMatrix[jaty][imod]->data.cpx[1].real*gain, JMatrix[jaty][0]->data.cpx[1].imag*gain);
		  fprintf (stderr,"          %10.7f, %10.7f,  %10.7f, %10.7f,\n", 
			   JMatrix[jaty][imod]->data.cpx[2].real*gain, JMatrix[jaty][imod]->data.cpx[2].imag*gain, 
			   JMatrix[jaty][imod]->data.cpx[3].real*gain, JMatrix[jaty][imod]->data.cpx[3].imag*gain);
		} /* End loop giving model with Jones matrix1 */
		offset = offsetChannel; /* Visibility offset */
		tre = visData[offset+0]*CosArr[0]+visData[offset+1]*SinArr[0];
		tim = -visData[offset+0]*SinArr[0]+visData[offset+1]*CosArr[0];
		fprintf (stderr,"XX: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f \n",
			 visData[offset+0],visData[offset+1],visData[offset+2],
			 sqrtf(visData[offset+0]*visData[offset+0]+visData[offset+1]*visData[offset+1]),
			 57.296*atan2f(visData[offset+1],visData[offset+0]), tre, tim);
		tre = visData[offset+3]*CosArr[0]+visData[offset+4]*SinArr[0];
		tim = -visData[offset+3]*SinArr[0]+visData[offset+4]*CosArr[0];
		fprintf (stderr,"YY: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f \n",
			 visData[offset+3],visData[offset+4],visData[offset+5],
			 sqrtf(visData[offset+3]*visData[offset+3]+visData[offset+4]*visData[offset+4]),
			 57.296*atan2f(visData[offset+4],visData[offset+3]), tre, tim);
		tre = visData[offset+6]*CosArr[0]+visData[offset+7]*SinArr[0];
		tim = -visData[offset+6]*SinArr[0]+visData[offset+7]*CosArr[0];
		if (doCrossPol) fprintf (stderr,"XY: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f \n",
					 visData[offset+6],visData[offset+7],visData[offset+8],
					 sqrtf(visData[offset+6]*visData[offset+6]+visData[offset+7]*visData[offset+7]),
					 57.296*atan2f(visData[offset+7],visData[offset+6]), tre, tim);
		tre = visData[offset+9]*CosArr[0]+visData[offset+10]*SinArr[0];
		tim = -visData[offset+9]*SinArr[0]+visData[offset+10]*CosArr[0];
		if (doCrossPol) fprintf (stderr,"YX: %8.5f, %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f \n",
					 visData[offset+9],visData[offset+10],visData[offset+11],
					 sqrtf(visData[offset+9]*visData[offset+9]+visData[offset+10]*visData[offset+10]),
					 57.296*atan2f(visData[offset+10],visData[offset+9]), tre, tim);
		tre = sumRealRR*CosArr[0]+sumImagRR*SinArr[0];
		tim = -sumRealRR*SinArr[0]+sumImagRR*CosArr[0];
		fprintf (stderr,"ModXX: %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f\n",
			 sumRealRR, sumImagRR,
			 sqrtf(sumRealRR*sumRealRR+sumImagRR*sumImagRR),
			 57.296*atan2f(sumImagRR,sumRealRR), tre, tim);
		tre = sumRealLL*CosArr[0]+sumImagLL*SinArr[0];
		tim = -sumRealLL*SinArr[0]+sumImagLL*CosArr[0];
		fprintf (stderr,"ModYY: %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f\n",
			 sumRealLL, sumImagLL,
			 sqrtf(sumRealLL*sumRealLL+sumImagLL*sumImagLL),
			 57.296*atan2f(sumImagLL,sumRealLL), tre, tim);
		tre = sumRealRL*CosArr[0]+sumImagRL*SinArr[0]; // DEBUG
		tim = -sumRealRL*SinArr[0]+sumImagRL*CosArr[0];
		if (doCrossPol) fprintf (stderr,"ModXY: %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f\n",
					 sumRealRL, sumImagRL,
					 sqrtf(sumRealRL*sumRealRL+sumImagRL*sumImagRL),
					 57.296*atan2f(sumImagRL,sumRealRL), tre, tim);
		tre = sumRealLR*CosArr[0]+sumImagLR*SinArr[0];
		tim = -sumRealLR*SinArr[0]+sumImagLR*CosArr[0];
		if (doCrossPol) fprintf (stderr,"ModYX: %8.5f, %8.5f, polar %8.5f %8.1f unwind %8.5f  %8.5f\n",
					 sumRealLR, sumImagLR,
					 sqrtf(sumRealLR*sumRealLR+sumImagLR*sumImagLR),
					 57.296*atan2f(sumImagLR,sumRealLR), tre, tim);
		// rerun ObitSkyModelVMBeamJonesCorSum to see in the debugger
		//ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
		//			      AmpArr, SinArr, CosArr,
		//			      &JMatrix[iaty][kt], &JMatrix[jaty][kt],
		//                            largs->cos2PA, largs->sin2PA, 
		//			      workVis, work1, work2, sumVis);
	      }}
#endif /* DEBUG */ 
	    /************************ end DEBUG **********************************/


	    kt = it;  /* offset in JMatrix */
	  } /* end outer loop over components */

	  break;
	case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	  /************************************************************
	                Gauss                                     
	  */
	  itab = 7;  /* Start of parameterized spectrum */
	  jt = 0;  /* index in chFlux */
	  /* From the AIPSish QGASUB.FOR  */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    /* Separate loops for with spEval and without */
	    if (doSpEval)  {   /* Evaluated spectrum? */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if ((chFlux[jt]!=0.0) && (chFlux[jt]!=fblank)){  /* valid? */
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  /* Gaussian factor */
		  arg = freq2 * (ddata[kadd+ka1]*visData[ilocu]*visData[ilocu] +
				 ddata[kadd+ka2]*visData[ilocv]*visData[ilocv] +
				 ddata[kadd+ka3]*visData[ilocu]*visData[ilocv]);
		  ExpArg[itcnt] = arg;
		  AmpArr[itcnt] = chFlux[jt];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
	      jt++;
	      }  /* end inner loop over components */
	    } else { /*  w/o SpEval */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[kamp]!=0.0) {
		  VecL[itcnt] = ddata[ku];
		  VecM[itcnt] = ddata[kv];
		  VecN[itcnt] = ddata[kw];
		  /* Frequency dependent terms */
		  lll = ll = log(specFreqFact);
		  arg = 0.0;
		  for (iterm=1; iterm<nterm; iterm++) {
		    arg += ddata[kfit+iterm] * lll;
		    lll *= ll;
		  }
		  ExpArg2[itcnt] = arg;
		  /* Gaussian factor */
		  arg = freq2 * (ddata[kadd+ka1]*visData[ilocu]*visData[ilocu] +
				 ddata[kadd+ka2]*visData[ilocv]*visData[ilocv] +
				 ddata[kadd+ka3]*visData[ilocu]*visData[ilocv]);
		  ExpArg[itcnt]  = arg;
		  AmpArr[itcnt]  = ddata[kfit];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
	      }  /* end inner loop over components */
	    } /* end w/o  SpEval */
	    
	    /* Compute phases */
	    SkyModel2DDot(itcnt, VecL, VecM, VecN, u, v, w, FazArr);
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Convert Gaussian exp arguments */
	    ObitExpVec(itcnt, ExpArg, ExpVal);
	    if (!doSpEval) {
	      /* Evaluate spectrum */
	      ObitExpVec(itcnt, ExpArg2, ExpVal2);
	      SkyModelVMul(itcnt, AmpArr, ExpVal2, AmpArr);
	    }
	    /* Gaussian correction */
	    SkyModelVMul(itcnt, AmpArr,  ExpVal, AmpArr);

	    /* Accumulate real and imaginary parts to sumVis */
	    ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					  AmpArr, SinArr, CosArr,
					  &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					  largs->cos2PA, largs->sin2PA, 
					  workVis, work1, work2, sumVis);
	    
	    sumRealRR += sumVis->data.array[0]; sumImagRR += sumVis->data.array[1];
	    sumRealLL += sumVis->data.array[6]; sumImagLL += sumVis->data.array[7];
	    if (doCrossPol) {
	      /* Cross pol */
	      sumRealRL += sumVis->data.array[2]; sumImagRL += sumVis->data.array[3];
	      sumRealLR += sumVis->data.array[4]; sumImagLR += sumVis->data.array[5];
	    } /* end xpol */
	    kt = it;  /* offset in JMatrix */
	  } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	  /************************************************************
	                Gauss Spec                                    
	  */
	  itab = 7;  /* Start of parameterized spectrum */
	  jt = 0;  /* index in chFlux */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    /* Separate loops for with spEval and without */
	    if (doSpEval)  {   /* Evaluated spectrum? */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if ((chFlux[jt]!=0.0) && (chFlux[jt]!=fblank)) {  /* valid? */
		  VecL[itcnt] = ddata[ku+nterm];
		  VecM[itcnt] = ddata[kv+nterm];
		  VecN[itcnt] = ddata[kw+nterm];
		  /* Gaussian factor */
		  arg = freq2 * (ddata[kadd+ka1]*visData[ilocu]*visData[ilocu] +
				 ddata[kadd+ka2]*visData[ilocv]*visData[ilocv] +
				 ddata[kadd+ka3]*visData[ilocu]*visData[ilocv]);
		  ExpArg[itcnt] = arg;
		  AmpArr[itcnt] = chFlux[jt];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
		jt++;
	      }  /* end inner loop over components */
	    } else { /*  without spEval */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[itab]!=0.0) {  /* valid? */
		  VecL[itcnt] = ddata[ku+nterm];
		  VecM[itcnt] = ddata[kv+nterm];
		  VecN[itcnt] = ddata[kw+nterm];
		  /* Amplitude from component flux and two gains */
		  AmpArr[itcnt]  = ddata[kamp];
		  /* Frequency dependent term */
		  lll = ll = log(specFreqFact);
		  arg = 0.0;
		  for (iterm=1; iterm<nterm; iterm++) {
		    arg += ddata[kfit+iterm] * lll;
		    lll *= ll;
		  }
		  ExpArg2[itcnt] = arg;
		  /* Gaussian term */
		  arg = freq2 * (ddata[kadd+ka1]*visData[ilocu]*visData[ilocu] +
				 ddata[kadd+ka2]*visData[ilocv]*visData[ilocv] +
				 ddata[kadd+ka3]*visData[ilocu]*visData[ilocv]);
		  ExpArg[itcnt] = arg;
		  AmpArr[itcnt]  = ddata[kfit];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
	      }  /* end inner loop over components */
	    } /* end w/o SpEval */
	    
	    /* Compute phases */
	    SkyModel2DDot(itcnt, VecL, VecM, VecN, u, v, w, FazArr);
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    if (!doSpEval) {
	      /* Evaluate spectrum */
	      ObitExpVec(itcnt, ExpArg2, ExpVal2);
	      SkyModelVMul(itcnt, AmpArr, ExpVal2, AmpArr);
	    }
	    /* Convert Gaussian exp arguments */
	    ObitExpVec(itcnt, ExpArg, ExpVal);
	    /* Gaussian correction */
	    SkyModelVMul(itcnt, AmpArr,  ExpVal, AmpArr);
	    
	    /* Accumulate real and imaginary parts to sumVis */
	    ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					  AmpArr, SinArr, CosArr,
					  &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					  largs->cos2PA, largs->sin2PA, 
					  workVis, work1, work2, sumVis);
	    
	    sumRealRR += sumVis->data.array[0]; sumImagRR += sumVis->data.array[1];
	    sumRealLL += sumVis->data.array[6]; sumImagLL += sumVis->data.array[7];
	    if (doCrossPol) {
	      /* Cross pol */
	      sumRealRL += sumVis->data.array[2]; sumImagRL += sumVis->data.array[3];
	      sumRealLR += sumVis->data.array[4]; sumImagLR += sumVis->data.array[5];
	    } /* end xpol */
	    kt = it;  /* offset in JMatrix */
	  } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_GaussModTSpec:     /* Gaussian on sky + tabulated spectrum*/
	  /************************************************************
	                Gauss TSpec                                    
	  */
	  itab = kamp + iSpec;  /* Closest subband index */
	  jt = 0;  /* index in chFlux */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    /* Separate loops for with spEval and without */
	    if (doSpEval)  {   /* Evaluated spectrum? */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if ((chFlux[jt]!=0.0) && (chFlux[jt]!=fblank)) {  /* valid? */
		  VecL[itcnt] = ddata[ku+nSpec];
		  VecM[itcnt] = ddata[kv+nSpec];
		  VecN[itcnt] = ddata[kw+nSpec];
		  /* Gaussian gain term */
		  arg = freq2 * (ddata[kadd+ka1]*visData[ilocu]*visData[ilocu] +
				 ddata[kadd+ka2]*visData[ilocv]*visData[ilocv] +
				 ddata[kadd+ka3]*visData[ilocu]*visData[ilocv]);
		  ExpArg[itcnt] = arg;
		  AmpArr[itcnt] = chFlux[jt];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
	      jt++;
	      }  /* end inner loop over components */
	    } else { /*  without spEval */
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[itab]!=0.0) {  /* valid? */
		  VecL[itcnt] = ddata[ku+nSpec];
		  VecM[itcnt] = ddata[kv+nSpec];
		  VecN[itcnt] = ddata[kw+nSpec];
		  AmpArr[itcnt] = ddata[kamp]; 
		  /* Gaussian gain term */
		  arg = freq2 * (ddata[kadd+ka1]*visData[ilocu]*visData[ilocu] +
				 ddata[kadd+ka2]*visData[ilocv]*visData[ilocv] +
				 ddata[kadd+ka3]*visData[ilocu]*visData[ilocv]);
		  ExpArg[itcnt] = arg;
		  /* Amplitude  - use spectral index */
		  AmpArr[itcnt] = ddata[itab];
		  ExpArg2[itcnt] = logNuONu0 * ddata[ksi];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
	      }  /* end inner loop over components */
	    } /* end w/o SpEval */

	    /* Compute phases */
	    SkyModel2DDot(itcnt, VecL, VecM, VecN, u, v, w, FazArr);
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Convert Gaussian exp arguments */
	    ObitExpVec(itcnt, ExpArg, ExpVal);
	    /* Gaussian correction */
	    SkyModelVMul(itcnt, AmpArr,  ExpVal, AmpArr);
	    
	    /* Evaluate spectral index? */
	    if (!doSpEval)  {
	      ObitExpVec(itcnt, ExpArg2, ExpVal2);
	      /* Spectral index correction */
	      SkyModelVMul(itcnt, AmpArr,  ExpVal2, AmpArr);}
	    
	    /* Accumulate real and imaginary parts to sumVis */
	    ObitSkyModelVMBeamJonesCorSum(itcnt, Stokes, isCirc, 
					  AmpArr, SinArr, CosArr,
					  &JMatrix[iaty][kt], &JMatrix[jaty][kt],
					  largs->cos2PA, largs->sin2PA, 
					  workVis, work1, work2, sumVis);
	    
	    sumRealRR += sumVis->data.array[0]; sumImagRR += sumVis->data.array[1];
	    sumRealLL += sumVis->data.array[6]; sumImagLL += sumVis->data.array[7];
	    if (doCrossPol) {
	      /* Cross pol */
	      sumRealRL += sumVis->data.array[2]; sumImagRL += sumVis->data.array[3];
	      sumRealLR += sumVis->data.array[4]; sumImagLR += sumVis->data.array[5];
	    } /* end xpol */
	    kt = it;  /* offset in JMatrix */
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
	
	/* RR/XX */
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
	    // MAYBE NOT if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	  } else {
	    /* Subtract model */
	    visData[offset]   -= modRealRR;
	    visData[offset+1] -= modImagRR;
	  }
	} /* end RR/YY not blanked */
	
	/* LL/YY */
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
	    //Maybe NOT if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	  } else {
	    /* Subtract model */
	    visData[offset]   -= modRealLL;
	    visData[offset+1] -= modImagLL;
	  }
	} /* end LL/YY not blanked */
	
	if (doCrossPol) {
	  /* RL/XY */
	  iStoke = 2;
	  offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	  
	  /* Ignore blanked data unless replacing the data */
	  if ((visData[offset+2]>0.0) || in->doReplace) {
	    /* Apply model to data */
	    if (in->doReplace) {  /* replace data with model */
	      visData[offset]   = modRealRL;
	      visData[offset+1] = modImagRL;
	      //Maybe NOT if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealRL;
	      visData[offset+1] -= modImagRL;
	    }
	  } /* end RL/XY not blanked */
	  
	  /* LR/YX */
	  offset += jincs;
	  /* Ignore blanked data unless replacing the data */
	  if ((visData[offset+2]>0.0) || in->doReplace) {
	    /* Apply model to data */
	    if (in->doReplace) {  /* replace data with model */
	      visData[offset]   = modRealLR;
	      visData[offset+1] = modImagLR;
	      //Maybe NOT if (visData[offset+2]<=0.0) visData[offset+2] = 1.0;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modRealLR;
	      visData[offset+1] -= modImagLR;
	    }
	  } /* end LR/YX not blanked */
	} /* end crosspol */	      
      } /* end vis loop */
    } /* end channel loop */
  } /* end IF loop */

  
  /* Indicate completion */
 finish: 
  /* DEBUG
  ObitThreadLock(in->thread); 
  Obit_log_error(err, OBIT_InfoErr,"Thread %d finished update %d vcnt %d iVis %d",
  largs->ithread, nup, vcnt, iVis);
  ObitThreadUnlock(in->thread);  */

  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadSkyModelVMBeamMFFTDFT */

/**
 * Returns the CC table to use for the current set of channels/IF
 * Not making relative PB corrections, this is the input CC table
 * In the latter case, the table should be Zapped when use is finished.
 * If not making relative Primary Beam correctsions then all selected,
 * else the next block for which the primary beam correction 
 * varies by less than 1% at the edge of the FOV.
 * If in->currentMode=OBIT_SkyModel_Mixed then the output table will be merged
 * and only contain entries with abs. flux densities in the range range.
 * If there are no components selected to process, the input table is 
 * always returned.
 * \param in       SkyModel
 *                 If info member doSmoo is present and TRUE then tabulated flux 
 *                 densities are fitted by a spectrum and replaced byr the spectrum 
 *                 evaluated at that frequency.
 * \param uvdata   UV data
 * \param field    Field number in in->mosaic
 * \param inCCVer  input CC table version
 * \param outCCver output CC table version number, 
 *                 0=> create new in which case the actual value is returned
 * \param startCC  [in] the desired first CC number (1-rel)
 *                 [out] the actual first CC number in returned table
 * \param endCC    [in] the desired highest CC number, 0=> to end of table
 *                 [out] the actual highest CC number in returned table
 * \param range    Range of allowed, merged CC fluxes.
 * \param err      Obit error stack object.
 * \return ObitCCTable to use, this should be Unref and Zapped when done 
 * NULL if no data in input table.
 */
static ObitTableCC* getPBCCTab (ObitSkyModelVMBeamMF* in, ObitUV* uvdata, 
				olong field, olong *inCCVer, olong *outCCVer,
				olong *startCC, olong *endCC, ofloat range[2],
				ObitErr *err)
{
  ObitTableCC *CCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitImageMF *image=NULL;
  ObitBeamShape *BeamShape=NULL;
  ObitBeamShapeClassInfo *BSClass;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  ofloat *flux=NULL, *sigma=NULL, *fitResult=NULL, pbmin=0.01, *PBCorr=NULL, alpha=0.0;
  ofloat *FreqFact=NULL, *sigmaField=NULL, ll, lll, arg, specFact;
  odouble *Freq=NULL, refFreq;
  odouble Angle=0.0;
  gpointer fitArg=NULL;
  olong irow, row, i, iterm, nterm, offset, nSpec, tiver;
  gchar keyword[16];
  gchar *routine = "ObitSkyModelVMBeamMF:getPBCCTab";

  /* error checks */
  if (err->error) return CCTable;
  g_assert (ObitSkyModelIsA(in));

  /* Compress/select CC table to new table */
  *outCCVer =  0;  /* Create new one */
  tiver = *inCCVer;
  CCTable = ObitTableCCUtilMergeSel2Tab (in->mosaic->images[field], tiver, outCCVer, 
					    *startCC, *endCC, range, err);
  /* Anybody home? */
  if (CCTable==NULL) {
    *startCC = 1;
    *endCC   = 0;
    return CCTable;
    }
  *startCC = 1;   /* Want all of these */
  *endCC   = CCTable->myDesc->nrow ;
  
  /* Smooth table - loop through table, correcting the tabulated spectral points
     for the primary gain, set weights to 1/PB^2, fit spectrum, replace values by
     fit evaluated at frequency uncorrected by PBcorr.
   */
  if (in->doSmoo) {
    image = (ObitImageMF*)in->mosaic->images[field];
   /* Make sure this is an ObitImageMF */
    Obit_retval_if_fail((ObitImageMFIsA(image)), err, CCTable,
			"%s: Image %s NOT an ObitImageMF", 
			routine, image->name);

    /* Create spectrum info arrays */
    nSpec = 1;
    ObitInfoListGetTest(image->myDesc->info, "NSPEC", &type, dim, &nSpec);
    nterm = 1;
    ObitInfoListGetTest(image->myDesc->info, "NTERM", &type, dim, &nterm);
    refFreq = image->myDesc->crval[image->myDesc->jlocf];
    Freq    = g_malloc0(nSpec*sizeof(odouble));
    FreqFact= g_malloc0(nSpec*sizeof(ofloat));
    /* get number of and channel frequencies for CC spectra from 
       CC table on first image in mosaic */
    if (nSpec>1) {
      for (i=0; i<nSpec; i++) {
	Freq[i] = 1.0;
	sprintf (keyword, "FREQ%4.4d",i+1);
	ObitInfoListGetTest(image->myDesc->info, keyword, &type, dim, &Freq[i]);
      }
    }
    
    if (in->prtLv>1) Obit_log_error(err, OBIT_InfoErr, 
				    "Smooth CCs to %d term spectrum",
				    nterm);
    
    /* Prior spectral index */
    InfoReal.flt = 0.0;   type = OBIT_float;
    ObitInfoListGetTest(image->myDesc->info, "ALPHA", &type, dim, &InfoReal);
    if (type==OBIT_double) alpha = (ofloat)InfoReal.dbl;
    if (type==OBIT_float)  alpha = (ofloat)InfoReal.flt;
  
    /* Log Freq ratio */
    for (i=0; i<nSpec; i++)  FreqFact[i] = log(Freq[i]/refFreq);

    /* Open CC table */
    retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadWrite, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, image->name, CCTable);
    /* Make sure table has nSpec channels */
    Obit_retval_if_fail((CCTable->noParms >= (4+nSpec)), err, CCTable,
			"%s: CC table %d appears not to have tabulated spectra", 
			routine, *outCCVer);

    /* Setup */
    offset     = 4;
    flux       = g_malloc0(nSpec*sizeof(ofloat));
    sigma      = g_malloc0(nSpec*sizeof(ofloat));
    sigmaField = g_malloc0(nSpec*sizeof(ofloat));
    PBCorr     = g_malloc0(nSpec*sizeof(ofloat));
    BeamShape  = ObitBeamShapeCreate ("BS", (ObitImage*)image, pbmin, in->antSize, TRUE);
    BSClass    = (ObitBeamShapeClassInfo*)(BeamShape->ClassInfo);
    fitArg     = ObitSpectrumFitMakeArg (nSpec, nterm, refFreq, Freq, FALSE, 
					 &fitResult, err);
    for (i=0; i<nterm; i++) fitResult[i]  = 0.0;
    for (i=0; i<nSpec; i++) sigmaField[i] = -1.0;
    if  (err->error) goto cleanup;
    
    /* Create table row */
    CCRow = newObitTableCCRow (CCTable);
    
    /* loop over table */
    for (irow=(*startCC); irow<=(*endCC); irow++) {
      
      /* Read */
      row = irow;
      retCode = ObitTableCCReadRow (CCTable, row, CCRow, err);
      if  (err->error) goto cleanup;

      /* Set field sigma to 0.01 of first */
      if (sigmaField[0]<0.0) {
	for (i=0; i<nSpec; i++) {
	  sigmaField[i] = 0.01 * CCRow->parms[offset+i];
	}
      }
      
      /* Primary beam stuff - Distance from Center  */
      Angle = ObitImageDescAngle(image->myDesc, CCRow->DeltaX, CCRow->DeltaY);

      /* Loop over spectral channels get corrected flux, sigma */
      for (i=0; i<nSpec; i++) {
	BeamShape->refFreq = Freq[i];  /* Set frequency */
	PBCorr[i] = BSClass->ObitBeamShapeGainSym(BeamShape, Angle);
	flux[i]   = CCRow->parms[offset+i] / PBCorr[i];
	sigma[i]  = sigmaField[i] / (PBCorr[i]*PBCorr[i]);
      }
 
      /* Fit spectrum */
      ObitSpectrumFitSingleArg (fitArg, flux, sigma, fitResult);
      
      /* Prior spectral index correction */
      if (nterm>=2) fitResult[1] += alpha;

      /* Replace channel fluxes with fitted spectrum */
      for (i=0; i<nSpec; i++) {
	/* Frequency dependent term */
	lll = ll = FreqFact[i];
	arg = 0.0;
	for (iterm=1; iterm<nterm; iterm++) {
	  arg += fitResult[iterm] * lll;
	  lll *= ll;
	}
	specFact = exp(arg);
	CCRow->parms[offset+i] = specFact*fitResult[0]*PBCorr[i];
      }
 
     /* ReWrite output */
      row = irow;
      retCode = ObitTableCCWriteRow (CCTable, row, CCRow, err);
      if  (err->error) goto cleanup;
    } /* end loop over table */
    
     /* Close */
    retCode = ObitTableCCClose (CCTable, err);
    if  (err->error) goto cleanup;
 
    /* Cleanup */
  cleanup:
    if (flux)       g_free(flux);
    if (sigma)      g_free(sigma);
    if (sigmaField) g_free(sigmaField);
    if (PBCorr)     g_free(PBCorr);
    if (Freq)       g_free(Freq);
    if (FreqFact)   g_free(FreqFact);
    if (fitResult)  g_free(fitResult);
    BeamShape = ObitBeamShapeUnref(BeamShape);
    ObitSpectrumFitKillArg(fitArg);
    CCRow = ObitTableRowUnref(CCRow);
   if  (err->error) Obit_traceback_val (err, routine, image->name, CCTable);
 } /* End smoothing table */
  
  return CCTable;
} /* end getPBCCTab */

/**
 * Determines if a new beam correction is needed based on channel/IF and the 
 * general model of the beam (Jinc).
 * A change of 0.2 percent is deemed necessary for an update.
 * \param in       SkyModel
 * \param uvdata   UV data
 * \param iIF      IF number (0-rel)
 * \param iChannel Channel, 0-rel
 * \param oldPB    Old primary beam correction
 * \param newPB    New primary beam correction
 * \param err      Obit error stack object.
 * \return TRUE if beam update needed, else FALSE
 */
static gboolean needNewPB(ObitSkyModelVMBeamMF* in, ObitUV* uvdata, 
			  olong iIF, olong iChannel,
			  ofloat oldPB, ofloat *newPB, ObitErr *err)
{
  gboolean update = FALSE;
  gboolean doJinc;
  ofloat FOV, PBFact;
  olong nfreq, ifreq, incf, incif;
  odouble Angle;
  ObitUVDesc *uvDesc=NULL;
  gchar *routine = "needNewPB";

  /* error checks */
  if (err->error) return update;

  *newPB = oldPB;  /* Until proven otherwise */

    /* Info from uv descriptor */
  uvDesc = uvdata->myDesc;
  nfreq = uvDesc->inaxes[uvDesc->jlocf];
  incf = uvDesc->incf/uvDesc->inaxes[0];   /* frequency increment in freq table */
  incif = uvDesc->incif/uvDesc->inaxes[0]; /* IF increment in freq table */
  if (uvDesc->jlocs<uvDesc->jlocif) incif /= uvDesc->inaxes[uvDesc->jlocs];
  if (uvDesc->jlocs<uvDesc->jlocf)  incf  /= uvDesc->inaxes[uvDesc->jlocs];

  /* Field of view? */
  if (in->pointFlux > 0.0) {
    FOV = sqrt (in->pointXOff*in->pointXOff + in->pointYOff*in->pointYOff);
    FOV = MAX (FOV, 0.1/3600.0);
  } else if (in->mosaic==NULL) {
    /* Use offset of point or 0.1 asec which ever is larger */
    FOV = sqrt (in->pointXOff*in->pointXOff + in->pointYOff*in->pointYOff);
    FOV = MAX (FOV, 0.1/3600.0);
  } else {
    /* Get Actual current from mosaic */
    FOV = ObitImageMosaicMaxFOV(in->mosaic, in->startComp, in->endComp, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, update);
  }
  Angle = FOV;  /* Field of view = max angle from point */
  
  /* which frequency channel is this? */
  ifreq = iIF * incif + iChannel * incf;
  ifreq = MIN ((nfreq-1), ifreq);
  ifreq = MAX (0, ifreq);

  /* Which beam shape function to use? */
  doJinc = (uvDesc->freqArr[ifreq] >= 1.0e9);
  
  /* Gain */
  if (doJinc) PBFact = ObitPBUtilJinc(Angle, uvDesc->freqArr[ifreq], in->antSize, 0.01);
  else        PBFact = ObitPBUtilPoly(Angle, uvDesc->freqArr[ifreq], 0.01);
  
  /* Need update? */
  if (fabs(PBFact-oldPB)>0.002) {
    /* Yes */
    update = TRUE;
    *newPB = PBFact;
  }

 return update;
} /* end needNewPB */

/** Private: Make model args structure */
void ObitSkyModelVMBeamMFMakeArg (ObitSkyModelVMBeamMF *in, ObitUV *uvdata, ObitErr *err)
{
  olong i, j, k, dim[2]={2,2}, nFreq=0;
  gboolean doSpEval=FALSE;   /* Interpolate spectrum? */
  odouble *Freqs=NULL, *specFreqEff=NULL;
  ObitImageMF *imageMF=NULL;
  VMBeamMFFTFuncArg* args;
  ObitInfoType type;
  gint32 ddim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar keyword[24];
  gchar *routine = "VMBeamMFMakeArg";

  /* First attempt to get Array of Effective freqnencies from in->info */
  if (ObitInfoListGetP(in->info, "specFreqEff", &type, ddim, (gpointer*)&specFreqEff)) {
    nFreq = ddim[0];
    Freqs = g_malloc0(nFreq*sizeof(odouble));
    for (i=0; i<nFreq; i++) Freqs[i] = specFreqEff[i]; 
    doSpEval = TRUE;
 } /* end got "specFreqEff */
  else {
    /* Get Frequency info from ImageMF image[0], if any */
    if ((in->mosaic)&&(in->mosaic->images)&&(in->mosaic->images[0])) {
      imageMF = (ObitImageMF*)in->mosaic->images[0];
      ObitImageOpen(in->mosaic->images[0],OBIT_IO_ReadOnly, err);
      if (ObitInfoListGetTest(imageMF->myDesc->info, "NSPEC", &type, ddim, &nFreq)) {
	if (nFreq>1) {
	  Freqs = g_malloc0(nFreq*sizeof(odouble));
	  for (i=0; i<nFreq; i++) {
	    Freqs[i] = 1.0;
	    sprintf (keyword, "FEFF%4.4d",i+1);  /* Try FEFFnnnn */
	    if (ObitInfoListGetTest(imageMF->myDesc->info, keyword, &type, ddim, &Freqs[i])) {
	    } else { /* USE FREQnnnn */
	      sprintf (keyword, "FREQ%4.4d",i+1);  /* Try FEFFnnnn */
	      ObitInfoListGetTest(imageMF->myDesc->info, keyword, &type, ddim, &Freqs[i]);
	    } /* end get Freq */
	    doSpEval = TRUE;
	  } /* end loop over freq */
	  //fprintf (stderr,"DEBUG Freq 12 %lg\n",Freqs[11]); // DEBUG HACK
	} /* end >1 freq */
	ObitImageClose(in->mosaic->images[0],err);
      }
    } /* End get info from image */
  }/* End get frequency info */
  
  
  if (!in->threadArgs) in->threadArgs = g_malloc0(in->nThreads*sizeof(VMBeamMFFTFuncArg*));
  for (i=0; i<in->nThreads; i++) {
     if (!in->threadArgs[i]) in->threadArgs[i] = g_malloc0(sizeof(VMBeamMFFTFuncArg)); 
    args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
    strcpy (args->type, "mfvmbeam");  /* Enter type as first entry */
    args->in     = (ObitSkyModelVMBeamMF*)in;
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
    args->dimGain = -1; /* Nothing yet */
    args->dimWork = 256;  /* Same size as arrays in FTDFT functions (FazArrSize) */
    args->workVis = g_malloc0(args->dimWork*sizeof(ObitMatx**));
    for (k=0; k<args->dimWork; k++) args->workVis[k] = ObitMatxCreate(OBIT_Complex, 2, dim);
    args->work1   = ObitMatxCreate(OBIT_Complex, 2, dim);
    args->work2   = ObitMatxCreate(OBIT_Complex, 2, dim);
    args->sumVis  = ObitMatxCreate(OBIT_Complex, 2, dim);
    if (doSpEval) {
      args->spEval = newObitSpectrumMFCreate("SpEval", nFreq, Freqs, in->refFreq, in->nTerm);
      args->chFlux = NULL; args->nChFlux = 0;
    } else {  args->spEval = NULL; args->chFlux = NULL;  args->nChFlux = 0;}
    args->isUnitJones = FALSE;
  } /* end loop over threads */
  if (Freqs) g_free(Freqs); /* Cleanup */
} /* end  ObitSkyModelVMBeamMFMakeArg */

/** Private: Update model args structure if the number of components has changed
 * \param in     Pointer to the ObitSkyModel .
 * \param args   argument list to update
 */
void ObitSkyModelVMBeamMFUpdateArg (ObitSkyModelVMBeamMF *in, VMBeamMFFTFuncArg *args)
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
    /* (Re)Build  in->JMatrix */
    for (j=0; j<in->numAntType; j++) {
      args->JMatrix[j] = g_malloc0(in->numComp*sizeof(ObitMatx**));
      for (k=0; k<in->numComp; k++) args->JMatrix[j][k] = ObitMatxCreate(OBIT_Complex, 2, dim);
    } /* end antenna type */
  } /* end rebuild arrays */
  args->dimGain = in->numComp;

  /* chFlux if spEval given */
  if ((args->spEval) && (args->nChFlux != in->numComp)) {
    if (args->chFlux) g_free(args->chFlux);  /* Delete old */
    args->chFlux = g_malloc0(in->numComp*sizeof(ofloat));
    args->nChFlux = in->numComp;
  }
} /* end ObitSkyModelVMBeamMFUpdateArg */

/** Private: Kill model arg structure */
gpointer ObitSkyModelVMBeamMFKillArg (ObitSkyModelVMBeamMF *in)
{
  olong i, j, k;
  VMBeamMFFTFuncArg *args;
  
  /* Check type - only handle "mfvmbeam" */
  args = (VMBeamMFFTFuncArg*)in->threadArgs[0];
  if ((strlen(args->type)>6) || (!strncmp(args->type, "mfvmbeam", 8))) {
    for (i=0; i<in->nThreads; i++) {
      args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
      if (args) {
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
	if (args->spEval) args->spEval = ObitSpectrumMFUnref(args->spEval);
	if (args->chFlux) g_free(args->chFlux);
	g_free(in->threadArgs[i]);
      } /* end if args valid */
    } /* Loop over threads */
    g_free(in->threadArgs);
    in->threadArgs = NULL;
  } /* end if this a "mfvmbeam" threadArg */
  return in->threadArgs;
} /* end ObitSkyModelVMBeamMFKillArg */

