/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011-2014                                          */
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
#include "ObitTableANUtil.h"
#include "ObitSkyGeom.h"
#include "ObitSinCos.h"
#include "ObitExp.h"
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
 * Widebane imaging version.
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

/** Over sampling factor in uv plane  */
olong OverSampleVMBeamMF=4;

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

/** Private: Load Components */
gboolean ObitSkyModelVMBeamMFLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
					ObitErr *err);

/** Private: Get Inputs. */
void  ObitSkyModelVMBeamMFGetInput (ObitSkyModel* inn, ObitErr *err);

/** Private: Threaded FTDFT */
static gpointer ThreadSkyModelVMBeamMFFTDFT (gpointer arg);

/** Private: Threaded FTDFT with phase correction */
static gpointer ThreadSkyModelVMBeamMFFTDFTPh (gpointer arg);

/** Private: Threaded FTGrid */
static gpointer ThreadSkyModelVMBeamMFFTGrid (gpointer arg);

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
 Note: this is identical to and labeled as the parent class */
typedef struct {
  /* type "mfvmbeam" in this class */
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
  /* Number of spectral channels in Interp */
  /* thread number, <0 -> no threading  */
  olong        ithread;
  /* Obit error stack object */
  ObitErr      *err;
  /* UV Interpolator for FTGrid */
  ObitCInterpolate *Interp;
  /* Start time (days) of validity of model */
  ofloat begVMModelTime;
  /* End time (days) of validity of model */
  ofloat endVMModelTime;
  /* Thread copy of Components list */
  ObitFArray *VMComps;
  /* VMBeam class entries */
  /* Amp, phase interpolator for I pol Beam image */
  ObitFInterpolate *BeamIInterp, *BeamIPhInterp;
  /* Amp, phase interpolator for V pol Beam image */
  ObitFInterpolate *BeamVInterp, *BeamVPhInterp;
  /* Amp, phase interpolator for Q pol Beam image */
  ObitFInterpolate *BeamQInterp, *BeamQPhInterp;
  /* Amp, phase interpolator for U pol Beam image */
  ObitFInterpolate *BeamUInterp, *BeamUPhInterp;
  /** Current uv channel number being processed.  */
  olong channel;
  /** Frequency of desired beam (Hz) corresponds to channel */
  odouble  BeamFreq;
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
  /** Number of spectral bins */
  olong        nSpec;
  /** Apply prior alpha correction? */
  gboolean doAlphaCorr;
  /** Prior spectral index correction */
  ofloat priorAlpha;
  /* Reference frequency (Hz) for Prior spectral index */
  odouble priorAlphaRefF;
  /** UV Interpolator array for FTGrid */
  ObitCInterpolate **Interps;
} VMBeamMFFTFuncArg;
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
 *      \li "xxxmaxGrid"         ofloat  Maximum absolute component flux to use in Gridded model 
 *      \li "xxxdoDFT"           boolean Something to do for DFT model?
 *      \li "xxxdoGrid"          boolean Something to do for Grid model?
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
  gchar *keyword=NULL, *None = "None", *value=NULL, *classType=NULL;
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
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->maxGrid);
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

  /* "xxxmaxGrid"         ofloat  Maximum absolute component flux to use 
     in Gridded model  */
  if (prefix) keyword = g_strconcat (prefix, "maxGrid", NULL);
  else        keyword = g_strdup("maxGrid");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->maxGrid);
  g_free(keyword);

  /* "xxxdoDFT"           boolean Something to do for DFT model? */
  if (prefix) keyword = g_strconcat (prefix, "doDFT", NULL);
  else        keyword = g_strdup("doDFT");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doDFT);
  g_free(keyword);

  /* "xxxdoGrid"          boolean Something to do for Grid model? */
  if (prefix) keyword = g_strconcat (prefix, "doGrid", NULL);
  else        keyword = g_strdup("doGrid");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doGrid);
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
ObitSkyModelVMBeamMF* 
ObitSkyModelVMBeamMFCreate (gchar* name, ObitImageMosaic* mosaic,
			  ObitUV *uvData,
			  ObitImage *IBeam,  ObitImage *VBeam, 
			  ObitImage *QBeam,  ObitImage *UBeam, 
			  ObitImage *IBeamPh,  ObitImage *VBeamPh, 
			  ObitImage *QBeamPh,  ObitImage *UBeamPh, 
			  ObitErr *err)
{
  ObitSkyModelVMBeamMF* out=NULL;
  olong number, i, nchan, nif;
  gchar *routine = "ObitSkyModelVMBeamMFCreate";

  /* Error tests */
  if (err->error) return out;  /* Previous error */

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
  VMBeamMFFTFuncArg *args;
  gchar keyword[12];
  gchar *routine = "ObitSkyModelVMBeamMFInitMod";

  if (err->error) return;

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(VMBeamMFFTFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(VMBeamMFFTFuncArg)); 
  
    for (i=0; i<in->nThreads; i++) {
      args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
      strcpy (args->type, "mfvmbeam");  /* Enter type as first entry */
      args->in      = inn;
      args->uvdata  = uvdata;
      args->ithread = i;
      args->err     = err;
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
      args->Interps= NULL;
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
  } /* end initialize */

  /* Call parent initializer */
  ObitSkyModelVMBeamInitMod(inn, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fourier transform threading routines */
  in->DFTFunc  = (ObitThreadFunc)ThreadSkyModelVMBeamMFFTDFT;
  in->GridFunc = (ObitThreadFunc)ThreadSkyModelVMBeamMFFTGrid;
  
   /* Init Sine/Cosine, exp calculator - just to be sure about threading */
  ObitSinCosCalc(phase, &sp, &cp);

  /* Create spectrum info arrays */
  nSpec = 1;
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

  /* Prior spectral index */
  InfoReal.flt = 0.0;   type = OBIT_float;
  ObitInfoListGetTest(image0->myDesc->info, "ALPHA", &type, dim, &InfoReal);
  if (type==OBIT_double) in->priorAlpha = (ofloat)InfoReal.dbl;
  if (type==OBIT_float)  in->priorAlpha = (ofloat)InfoReal.flt;

  in->priorAlphaRefF = in->refFreq;
  ObitInfoListGetTest(image0->myDesc->info, "RFALPHA", &type, dim, &in->priorAlphaRefF);
  
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
  if (err->prtLv>1) {
    if (in->currentMode==OBIT_SkyModel_DFT)
      Obit_log_error(err, OBIT_InfoErr, "SkyModelVMBeamMF using DFT calculation type");
    else if (in->currentMode==OBIT_SkyModel_Grid)
      Obit_log_error(err, OBIT_InfoErr, "SkyModelVMBeamMF using Grid calculation type");
    else if (in->currentMode==OBIT_SkyModel_Fastest)
      Obit_log_error(err, OBIT_InfoErr, "SkyModelVMBeamMF using Fastest calculation type");
  }
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

 olong i, k;
  VMBeamMFFTFuncArg *args;

  in->myInterp = ObitCInterpolateUnref(in->myInterp);
  in->plane    = ObitFArrayUnref(in->plane);
  ObitThreadPoolFree (in->thread);  /* Shut down any threading */
  if (in->threadArgs) {
    /* Check type - only handle "mfvmbeam" */
    if (!strncmp((gchar*)in->threadArgs[0], "mfvmbeam", 8)) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
	args->BeamIInterp = ObitFInterpolateUnref(args->BeamIInterp);
	args->BeamVInterp = ObitFInterpolateUnref(args->BeamVInterp);
	args->BeamQInterp = ObitFInterpolateUnref(args->BeamQInterp);
	args->BeamUInterp = ObitFInterpolateUnref(args->BeamUInterp);
	args->BeamIPhInterp = ObitFInterpolateUnref(args->BeamIPhInterp);
	args->BeamVPhInterp = ObitFInterpolateUnref(args->BeamVPhInterp);
	args->BeamQPhInterp = ObitFInterpolateUnref(args->BeamQPhInterp);
	args->BeamUPhInterp = ObitFInterpolateUnref(args->BeamUPhInterp);
	args->VMComps       = ObitFArrayUnref(args->VMComps);
	if (args->Rgain)  g_free(args->Rgain);
	if (args->Lgain)  g_free(args->Lgain);
	if (args->Qgain)  g_free(args->Qgain);
	if (args->Ugain)  g_free(args->Ugain);
	if (args->Rgaini) g_free(args->Rgaini);
	if (args->Lgaini) g_free(args->Lgaini);
	if (args->Qgaini) g_free(args->Qgaini);
	if (args->Ugaini) g_free(args->Ugaini);
	if (args->Interps) {
	  for (k=0; k<args->nSpec; k++) 
	    args->Interps[k] = ObitCInterpolateUnref(args->Interps[k]);
	  g_free(args->Interps);
	}
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
      in->nThreads   = 0;
    } /* end if this a "mfvmbeam" threadArg */
  }


  /* Cleanup arrays */
  if (in->specFreq)  g_free(in->specFreq);  in->specFreq = NULL;
  if (in->specIndex) g_free(in->specIndex); in->specIndex = NULL;
  
  /* Call parent shutdown */
  ObitSkyModelVMBeamShutDownMod(inn, uvdata, err);
  
} /* end ObitSkyModelVMBeamMFShutDownMod */

/**
 * Initializes an ObitSkyModel for a pass through data in time order.
 * Resets current times, converts field offsets of components to pointing offsets
 * \param in  SkyModel to initialize
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamMFInitModel (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;

  /* Call parent InitModel */
  ObitSkyModelVMBeamInitModel (inn, err);

  /* Fourier transform routines - DFT only */
  /* Are phases given? */
  if (in->IBeamPh) 
    in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamMFFTDFTPh;
  else /* No phase */
    in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMBeamMFFTDFT;

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
 *      \li "xxxmaxGrid"         ofloat  Maximum absolute component flux to use in Gridded model 
 *      \li "xxxdoDFT"           boolean Something to do for DFT model?
 *      \li "xxxdoGrid"          boolean Something to do for Grid model?
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

  /* "xxxmaxGrid"         ofloat  Maximum absolute component flux to use 
     in Gridded model  */
  if (prefix) keyword = g_strconcat (prefix, "maxGrid", NULL);
  else        keyword = g_strdup("maxGrid");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->maxGrid);
  g_free(keyword);

  /* "xxxdoDFT"           boolean Something to do for DFT model? */
  if (prefix) keyword = g_strconcat (prefix, "doDFT", NULL);
  else        keyword = g_strdup("doDFT");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doDFT);
  g_free(keyword);

  /* "xxxdoGrid"          boolean Something to do for Grid model? */
  if (prefix) keyword = g_strconcat (prefix, "doGrid", NULL);
  else        keyword = g_strdup("doGrid");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doGrid);
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
  theClass->ObitSkyModelGridComps    = (ObitSkyModelGridCompsFP)ObitSkyModelVMBeamMFGridComps;
  theClass->ObitSkyModelLoadImage    = (ObitSkyModelLoadImageFP)ObitSkyModelVMBeamMFLoadImage;
  theClass->ObitSkyModelFTDFT        = (ObitSkyModelFTDFTFP)ObitSkyModelVMBeamMFFTDFT;
  theClass->ObitSkyModelFTGrid       = (ObitSkyModelFTGridFP)ObitSkyModelVMBeamMFFTGrid;
  theClass->ObitSkyModelLoadComps    = (ObitSkyModelLoadCompsFP)ObitSkyModelVMBeamMFLoadComps;
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
  VMBeamMFFTFuncArg *args;
  ObitSkyModelVMBeamMF *in = inn;

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
    /* Check type - only handle "mfvmbeam" */
    args = (VMBeamMFFTFuncArg*)in->threadArgs[0];
    if ((strlen(args->type)>8) || (!strncmp(args->type, "mfvmbeam", 8))) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
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
    } /* end if this a "vmbeammf" threadArg */
  }

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
 * \param in     Pointer to the ObitSkyModel .
 * \param uvdata UV data set
 */
void  ObitSkyModelVMBeamMFChose (ObitSkyModel* inn, ObitUV* uvdata) 
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
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
    ObitSkyModelVMBeamChose (inn, uvdata);
    in->maxGrid = 1.0e20;
    in->minDFT  = 0.0;
    return;
  }

} /* end ObitSkyModelVMBeamMFChose */

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
 * Grid components onto in->planes (zeroed arrays the twice the size 
 * of the image) and Fourier transformed to in->FTplanes.
 * Scaling of components and any tapering is applied.
 * Grid is double size for increased accuracy.
 * For convenience in interpolation, HWIDTH columns are added by 
 * copying from the positive half plane.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * Multiplies by factor member and any prior spectral index correction.
 * \param inn    Pointer to theObitSkyModelVMBeamMF .
 * \param field  field number (0-rel) in in->mosaic->images
 * \param uvdata UV data set to model
 * \param err    Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelVMBeamMFGridFTComps (ObitSkyModel* inn, olong field, ObitUV* uvdata, 
					  ObitErr *err)
{
  ObitSkyModelVMBeamMF *in  = (ObitSkyModelVMBeamMF*)inn;
  gboolean gotSome = FALSE;
  ObitImageDesc *imDesc = NULL;
  olong i, j, k, nx, ny;
  olong ncomp, ndim, naxis[2];
  ofloat gparm[3], dU, dV, UU, VV, texp;
  ofloat konst, xmaj, xmin, cpa, spa, b1, b2, b3, bb2, bb3;
  ofloat taper, *grid, factor[2];
  gboolean doGaus;
  ObitCArray *FFTImage = NULL;
  gchar *routine = "ObitSkyModelVMBeamMFGridFTComps";
  /* DEBUG 
  ObitFArray *tempFArray = NULL; */
  /* END DEBUG */

  /* error check */
  if (err->error) return gotSome ;

  /* Create grid, sum components into in->planes */
  ObitSkyModelVMBeamMFLoadGridComps (inn, field, uvdata, gparm, &ncomp, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  /* Don't bother if no components requested */
  gotSome = ncomp>=1;
  if (!gotSome) return gotSome;

  /* DEBUG 
     ObitImageUtilArray2Image ("DbugGriddedComps.fits", 0, in->planes[0], err);
     if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
     fprintf(stderr,"After ObitSkyModelLoadGridComps\n"); */
  /* END DEBUG */

  /* Create output arrays if needed */
  if (in->FTplanes==NULL)
    in->FTplanes = g_malloc0(in->nSpec*sizeof(ObitCArray*));
  if (in->myInterps==NULL)
    in->myInterps = g_malloc0(in->nSpec*sizeof(ObitCInterpolate*));

  /* Output of FFT */
  ndim = 2;
  naxis[0] = 1+in->planes[0]->naxis[0]/2; naxis[1] = in->planes[0]->naxis[1]; 
  FFTImage = ObitCArrayCreate ("FFT output", ndim, naxis);
  
  /* Loop over spectral planes */
  for (k=0; k<in->nSpec; k++) {
    
    /* Fourier Transform image */
    ObitSkyModelVMBeamMFFTImage (inn, in->planes[k], FFTImage);

    /* Release image plane */
    in->planes[k] = ObitFArrayUnref(in->planes[k]);
    
    /* DEBUG
       tempFArray = ObitCArrayMakeF(FFTImage);
       ObitCArrayReal (FFTImage, tempFArray); 
       ObitImageUtilArray2Image ("DbugFFTReal.fits", 0, tempFArray, err);
       tempFArray = ObitFArrayUnref(tempFArray);
       if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
       tempFArray = ObitCArrayMakeF(FFTImage);
       ObitCArrayImag (FFTImage, tempFArray); 
       ObitImageUtilArray2Image ("DbugFFTImag.fits", 0, tempFArray, err);
       tempFArray = ObitFArrayUnref(tempFArray);
       if (err->error) Obit_traceback_val (err, routine, in->name, gotSome); */
    /* END DEBUG */
    
    imDesc = in->mosaic->images[field]->myDesc; 
    
    /* Add taper if necessary */
    /* Are these Gaussians? */
    doGaus = (gparm[0]>0.0) || (gparm[1]>0.0);
    /* If tapering, create array, set constants */
    if (doGaus) {
      /* Image info - descriptor should still be valid */
      nx = OverSampleVMBeamMF*imDesc->inaxes[imDesc->jlocr];
      ny = OverSampleVMBeamMF*imDesc->inaxes[imDesc->jlocd];
      
      /* UV cell spacing */
      dU = RAD2DG /  (nx * fabs(imDesc->cdelt[imDesc->jlocr]));
      dV = RAD2DG /  (ny * fabs(imDesc->cdelt[imDesc->jlocd]));
      
      konst = DG2RAD * G_PI * sqrt (0.5) / 1.17741022;
      xmaj = gparm[0] * konst;
      xmin = gparm[1] * konst;
      cpa = cos (DG2RAD * (90.0+gparm[2])); /* FFTW grid different from AIPS */
      spa = sin (DG2RAD * (90.0+gparm[2]));
      b1 = -(((cpa*xmaj)*(cpa*xmaj)) + ((spa*xmin)*(spa*xmin)));
      b2 = -(((spa*xmaj)*(spa*xmaj)) + ((cpa*xmin)*(cpa*xmin)));
      b3 = - 2.0 * spa * cpa * (xmaj*xmaj - xmin*xmin);
      
      /* pointer to complex grid */
      ndim = 2; naxis[0] = 0; naxis[1] = 0; 
      grid = ObitCArrayIndex(FFTImage, naxis);
  
      /* loop over uv array */  
      for (i=0; i<ny; i++) {
	VV = dV * (i-nx/2);
	UU = 0.0;
	bb2 = b2 * VV * VV;
	bb3 = b3 * VV;
	/* Loop down row computing, applying taper */
	for (j=0; j<1+nx/2; j++) {
	  texp = b1 * UU * UU + bb2 + bb3 * UU;
	  if (texp>-14.0) taper = exp (texp);
	  else  taper = 0.0;
	  UU = UU + dU;
	  grid[2*j]   *= taper;
	  grid[2*j+1] *= taper;
	}
	grid += 2*FFTImage->naxis[0];
      }
    } /* end tapering */
    
    /* Add conjugate columns for interpolator */
    in->numConjCol = HWIDTH;  /* Number of columns on conjugate side of plane */
    in->FTplanes[k] = ObitCArrayUnref(in->FTplanes[k]);
    in->FTplanes[k] = ObitCArrayAddConjg(FFTImage, in->numConjCol);
    
    /* DEBUG */
    /*tempFArray = ObitCArrayMakeF(in->FTplane);*/  /* Temp FArray */
    /*ObitCArrayReal (in->FTplane, tempFArray);*/   /* Get real part */
    /*ObitImageUtilArray2Image ("DbugConjgReal.fits", 0, tempFArray, err);*/
    /*tempFArray = ObitFArrayUnref(tempFArray); */  /* delete temporary */
    /*if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);*/
    /*fprintf(stderr,"After ObitCArrayAddConjg\n");*/
    /* END DEBUG */
    
    /* (re)Create interpolator */
    factor[0] = OverSampleVMBeamMF; factor[1] = OverSampleVMBeamMF;
    in->myInterps[k] = ObitCInterpolateUnref(in->myInterps[k]);
    in->myInterps[k] = 
      newObitCInterpolateCreate("UV data interpolator", in->FTplanes[k], imDesc,
				factor[0], factor[1], in->numConjCol, HWIDTH, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  } /* end loop over planes */

  /* Cleanup */
  FFTImage  = ObitCArrayUnref(FFTImage);


  return gotSome;
} /* end ObitSkyModelVMBeamMFGridFTComps */

/**
 * Create arrays OverSampleVMBeamMF times the size of the input image (in->planes) 
 * and sum components onto them.
 * Grid is oversize for increased accuracy.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * \param inn    Pointer to the ObitSkyModelVMBeamMF .
 * \param field  field number (0-rel) in in->mosaic->images
 * \param uvdata UV data set to model
 * \param gparm  [out] the parameters of the Gaussians in the table
 *               [-1,-1,-1] => not Gaussian.
 * \param ncomp  Actual number of components in in->comps
 * \param err    Obit error stack object.
 */
void  ObitSkyModelVMBeamMFLoadGridComps (ObitSkyModel* inn, olong field, ObitUV* uvdata, 
					 ofloat gparm[3], olong *ncomp, ObitErr *err)
{
  ObitSkyModelVMBeamMF *in  = (ObitSkyModelVMBeamMF*)inn;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableCC *CCTable = NULL;
  ObitImageDesc *imDesc = NULL;
  ofloat range[2], specCorr;
  olong k;
  gchar *tabType = "AIPS CC";
  olong outCCVer, ver, first, last, startComp, endComp;
  gchar *routine = "ObitSkyModelVMBeamMFLoadGridComps";

  /* error check */
  if (err->error) return;

  /* Any components? */
  if ((in->endComp[field]<in->startComp[field]) || (in->endComp[field]<=0)) {
    *ncomp = 0;
    return;
  }

  /* Open Image */
  /* Use external buffer (Not actually reading image here) */
  in->mosaic->images[field]->extBuffer = TRUE;
  retCode = ObitImageOpen (in->mosaic->images[field], OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);

  /* Get CC table */
  outCCVer = 0;
  ver = in->CCver[field];
  startComp = in->startComp[field];
  endComp = in->endComp[field];
  range[0] = 0.0;  /* Range of merged fluxes for Grid */
  range[1] = in->maxGrid;
  CCTable = getPBCCTab (in, uvdata, field, &ver, &outCCVer, 
			&startComp, &endComp, range, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  in->CCver[field] = ver;  /* save if defaulted (0) */
  
  /* Grid planes */
  first = startComp;
  last  = endComp;
  /* If noNeg last = last before first negative */
  imDesc = in->mosaic->images[field]->myDesc;
  if (in->planes==NULL)
    in->planes = g_malloc0(in->nSpec*sizeof(ObitFArray*));
  for (k=0; k<in->nSpec; k++) {

    /* Spectral correction for prior alpha array */
     if (in->doAlphaCorr && (in->priorAlpha!=0.0)) {
       specCorr = pow((in->specFreq[k]/in->priorAlphaRefF), in->priorAlpha);
    } else { /* No correction */
      specCorr = 1.0;
    }
    
    retCode = ObitTableCCUtilGridSpect (CCTable, OverSampleVMBeamMF, k+1,
					&first, &last, in->noNeg,
					in->factor*specCorr, 
					in->minFlux, in->maxGrid,
					imDesc, &in->planes[k], gparm, 
					ncomp, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, in->name);
  } /* end loop oer planes */

  /* Save values of highest comp - probably bad*/
  if (outCCVer==0) {
    /* no translation of table */
    /* Bad in->startComp[field] = first;
       in->endComp[field] = last; */
  } else {
    /* Translated table with only selected values */
    /* Bad in->endComp[field] = in->startComp[field] + last-first; */
  }
  
  /* if outCCver>0 then the CCtable is temporary - Zap */
  if (outCCVer>0) {
    CCTable = ObitTableCCUnref (CCTable);
    ObitImageZapTable(in->mosaic->images[field], tabType, outCCVer, err);
  /* else simply release table  */
  } else CCTable = ObitTableCCUnref (CCTable);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Close Image */
  retCode = ObitImageClose (in->mosaic->images[field], err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);
  
  /* Unset use external buffer switch */
  in->mosaic->images[field]->extBuffer = FALSE;
  
} /* end ObitSkyModelVMBeamMFLoadGridComps */

/**
 * Fourier Transform image array in in->plane, 
 * Half plane complex returned in center-at-the-center order.
 * \param inn      the ObitSkyModelVMBeamMF .
 * \param inArray  Array to be Transformed.
 * \param outArray Output of FFT, half plane complex
 */
void  ObitSkyModelVMBeamMFFTImage (ObitSkyModel* inn, ObitFArray *inArray, 
				   ObitCArray *outArray)
{
  /*ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;*/
  olong naxis[2];
  ObitFFT *myFFT;

  /* Swaparoonie to FFT order */
  ObitFArray2DCenter (inArray);

  /* Create FFT */
  naxis[0] = inArray->naxis[0]; naxis[1] = inArray->naxis[1];
  myFFT = newObitFFT("FFT:FTImage", OBIT_FFT_Forward, 
		     OBIT_FFT_HalfComplex, 2, naxis);

  /* FFT */
  ObitFFTR2C (myFFT, inArray, outArray);

  /* Put the center at the center */
  ObitCArray2DCenter (outArray);

  /* Cleanup */
  myFFT     = ObitFFTUnref(myFFT);

} /* end ObitSkyModelVMBeamMFFTImage  */


/**
 * Load components model into in comps member.
 * If the frequency axis has ctype "SPECLNMF" and if "NSPEC" exists in the 
 * first image descriptor InfoList and is > 0 then there should be a spectrum 
 * in each CLEAN component.
 * Multiplies by factor member and any prior spectral index correction.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Allows mixed points and Gaussians, one per CC Table
 * Adapted from the AIPSish QNOT:IONDFT
 * \param inn  SkyModel 
 * \param n   Image number on mosaic, if -1 load all images
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * Output is in member comps, the entries are
 * \li 0 Field (0-rel)
 * \li 1 CC DeltaX
 * \li 2 CC DeltaY
 * \li 3... nSpec Amplitude (Jy) values
 * \li 4+nadd -2*pi*x (radians), nadd = nSpec - 1
 * \li 5+nadd -2*pi*y (radians)
 * \li 6+nadd -2*pi*z (radians)
 * \li Other model parameters depending on model type
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
  olong warray, larray, iterm, nterm, maxTerm=1, toff;
  ofloat *array, parms[20], range[2], gp1=0., gp2=0., gp3=0.;
  olong ver, i, j, hi, lo, count, ncomp, startComp, endComp, irow, lrec;
  olong outCCVer, ndim, naxis[2], lenEntry, nadd, nspec;
  ofloat *table, xxoff, yyoff, zzoff;
  ofloat konst, konst2, xyz[3], xp[3], umat[3][3], pmat[3][3];
  ofloat ccrot, ssrot, xpoff, ypoff, maprot, uvrot;
  ofloat dxyzc[3], cpa, spa, xmaj, xmin, *specCorr=NULL;
  gboolean doCheck=FALSE, want, do3Dmul;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitSkyModelVMBeamMFLoadComps";
  
  /* error checks */
  if (err->error) return gotSome;

  /* Don't bother if no components requested */
  if ((n>=0) && (in->startComp[n]>in->endComp[n])) return gotSome;

  /* Uv descriptor */
  uvDesc = uvdata->myDesc;

  konst = DG2RAD * 2.0 * G_PI;
  /* konst2 converts FWHM(deg) to coefficients for u*u, v*v, u*v */
  konst2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5);
  /*konst2 = DG2RAD * 2.15169;*/

  /* Loop over images counting CCs */
  count = 0;
  in->modType = OBIT_SkyModel_Unknown; /* Model type not known */
  if (in->mosaic) {lo = 0; hi = in->mosaic->numberImages-1;}
  else {lo = 0; hi = 0;}
  if (n>=0) {lo = n; hi = n;}
  for (i=lo; i<=hi; i++) {

    /* Expect anything in this table? */
    if ((in->startComp[i]>in->endComp[i]) || (in->endComp[i]<=0)) continue;

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

    /* Is spectral information included? */
    if (!strncmp (in->mosaic->images[i]->myDesc->ctype[in->mosaic->images[i]->myDesc->jlocf], 
		  "SPECLOGF", 8)) {
      nterm = in->mosaic->images[i]->myDesc->inaxes[in->mosaic->images[i]->myDesc->jlocf];
      ObitInfoListGetTest (in->mosaic->images[i]->myDesc->info, "NTERM", &type, dim, &nterm);
      maxTerm = MAX (maxTerm, nterm);
      in->nSpecTerm = nterm -1;  /* Only higher order terms */

    }
    if (!strncmp (in->mosaic->images[i]->myDesc->ctype[in->mosaic->images[i]->myDesc->jlocf], 
		  "SPECLNMF", 8)) {
      /* IO descriptor give true size */
      imIODesc = (ObitImageDesc*)in->mosaic->images[0]->myIO->myDesc;
      nspec = imIODesc->inaxes[imIODesc->jlocf];
      ObitInfoListGetTest (in->mosaic->images[i]->myDesc->info, "NSPEC", &type, dim, &nspec);
      in->nSpec = nspec;  /* Number of terms in the spectrum */
      maxTerm = MAX (maxTerm, nspec);
    }

  } /* end loop counting CCs */

  /* Use mode type, nterm of the highest encountered */
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
  naxis[0] = 7; naxis[1] = count;
  if (in->modType==OBIT_SkyModel_GaussMod)        naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_GaussModSpec)    naxis[0] += 3; /* Gaussian + spectrum */
  if (in->modType==OBIT_SkyModel_GaussModTSpec)   naxis[0] += 3; /* Gaussian  + tabuated spectrum */
  if (in->modType==OBIT_SkyModel_USphereMod)      naxis[0] += 2; /* Uniform sphere */
  if (in->modType==OBIT_SkyModel_USphereModSpec)  naxis[0] += 2; /* Uniform sphere + spectrum*/
  if (in->modType==OBIT_SkyModel_USphereModTSpec) naxis[0] += 2; /* Uniform sphere + tabuated spectrum */

  /* Any spectral terms */
  nadd = in->nSpec - 1;  /* Number of spectral points minus one */
  naxis[0] += nadd;
  lenEntry = naxis[0];  /* Length of table entry */

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
  ncomp = 0;
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
    
    /* Save values of highest comp - probably bad */
    if (outCCVer==0) {
      /* no translation of table */
      /*??? in->endComp[i] = endComp; */
    } else {
      /* Translated table with only selected values */
      /*??? in->endComp[i] = in->startComp[i] + endComp-startComp; */
    }
    
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
 
    /* DEBUG - replace model with fitted beam */
    if(modType==OBIT_SkyModel_GaussModTSpec) {
      parms[0] = imDesc->beamMaj;
      parms[1] = imDesc->beamMin;
      parms[2] = imDesc->beamPA;
    }

    /* Gaussian parameters */
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
      if (in->noNeg && (array[0]<=0.0)) break;

      /* Do we want this one? */
      want = (fabs(array[0])>0.0);
      want = want && (fabs(array[0])>in->minFlux);
      want = want && (ncomp<count);  /* don't overflow */
      if (want) {

	/* Field number */
	table[0] = i;
	table[1] = array[1];
	table[2] = array[2];

	/* Point */
	table[3] = array[0] * in->factor;
	xp[0] = (array[1] + xpoff) * konst;
	xp[1] = (array[2] + ypoff) * konst;
	if (CCTable->DeltaZCol>=0) xp[2] = array[3] * konst;
	else                       xp[2] = 0.0;
	if (do3Dmul) {
	  xyz[0] = xp[0]*umat[0][0] + xp[1]*umat[1][0];
	  xyz[1] = xp[0]*umat[0][1] + xp[1]*umat[1][1];
	  xyz[2] = xp[0]*umat[0][2] + xp[1]*umat[1][2]; 
	  /* DEBUG
	  xyz[0] = xp[0]*umat[0][0] + xp[1]*umat[0][1];
	  xyz[1] = xp[0]*umat[1][0] + xp[1]*umat[1][1];
	  xyz[2] = xp[0]*umat[2][0] + xp[1]*umat[2][1]; */
	  /* PRJMUL (2, XP, UMAT, XYZ); */
	} else {  /* no rotation  */
	  xyz[0] = ccrot * xp[0] + ssrot * xp[1];
	  xyz[1] = ccrot * xp[1] - ssrot * xp[0];
	  xyz[2] = xp[2];
 	}
	table[4+nadd] = xyz[0] + xxoff;
	table[5+nadd] = xyz[1] + yyoff;
	table[6+nadd] = xyz[2] + zzoff;

	/* Zero rest in case */
	for (iterm=7+nadd; iterm<lenEntry; iterm++) table[iterm] = 0.0;

	/* Only same type as highest */
	if  (in->modType==modType) {
	  /* Only Point */
	  if (modType==OBIT_SkyModel_PointMod){
	    /* Nothing special this case */
	    
	    /* Only Point with tabulated spectrum */
	  } else if (modType==OBIT_SkyModel_PointModTSpec) {
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+3] = array[iterm+toff]*specCorr[iterm];
	    
	    /* Only Point with spectrum */
	  } else if (modType==OBIT_SkyModel_PointModSpec) {
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+7+nadd] = array[iterm+toff];
	    
	    /* Only Gaussian */
	  } else if (in->modType==OBIT_SkyModel_GaussMod) {
	    table[7+nadd] = gp1;
	    table[8+nadd] = gp2;
	    table[9+nadd] = gp3;
	    
	    /* Only Gaussian + tabulated spectrum */
	  } else if (in->modType==OBIT_SkyModel_GaussModTSpec) {
	    table[7+nadd] = gp1;
	    table[8+nadd] = gp2;
	    table[9+nadd] = gp3;
	    /*  spectrum */
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+3] = array[iterm+toff]*specCorr[iterm];
	    
	    /* Only Gaussian + spectrum */
	  } else if (in->modType==OBIT_SkyModel_GaussModSpec) {
	    table[7+nadd] = gp1;
	    table[8+nadd] = gp2;
	    table[9+nadd] = gp2;
	    /*  spectrum */
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+3] = array[iterm+toff];
	    
	    /* Only Uniform sphere */
	  } else if (in->modType==OBIT_SkyModel_USphereMod) {
	    table[3] = 3.0 * array[0] * in->factor;
	    table[7+nadd] = parms[1]  * 0.109662271 * 2.7777778e-4;
	    table[8+nadd] = 0.1;
	    
	    /* Only Uniform sphere + tabulated spectrum */
	  } else if (in->modType==OBIT_SkyModel_USphereModTSpec) {
	    table[3] = 3.0 * array[0] * in->factor;
	    table[4+nadd] = parms[1]  * 0.109662271 * 2.7777778e-4;
	    table[5+nadd] = 0.1;
	    /*  spectrum */
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+3] = array[iterm+toff]*specCorr[iterm];

	    /* Only Uniform sphere+ spectrum */
	  } else if (in->modType==OBIT_SkyModel_USphereModSpec) {
	    table[3] = 3.0 * array[0] * in->factor;
	    table[7+nadd] = parms[1]  * 0.109662271 * 2.7777778e-4;
	    table[8+nadd] = 0.1;
	    /*  spectrum */
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+3] = array[iterm+toff];

	  }
	} else { /* Mixed type - zero unused model components */

	  /* Only Point here but some Gaussian - zero Gauss comps */
	  if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_GaussMod)) {
	    table[7+nadd] = 0.0;
	    table[8+nadd] = 0.0;
	    table[9+nadd] = 0.0;
	  
	    /* Gauss here but also some points */
	  } else if ((modType==OBIT_SkyModel_GaussMod) && (in->modType==OBIT_SkyModel_PointMod)) {
	    table[7+nadd] = gp1;
	    table[8+nadd] = gp2;
	    table[9+nadd] = gp3;
	  
	  /* Only PointSpectrum here but some GaussianSpectrum - zero Gauss comps */
	  } else if ((modType==OBIT_SkyModel_PointModTSpec) && (in->modType==OBIT_SkyModel_GaussModTSpec)) {
	    table[7+nadd] = 0.0;
	    table[8+nadd] = 0.0;
	    table[9+nadd] = 0.0;
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+3] = array[iterm+toff]*specCorr[iterm];
	  
	    /* GaussianTSpectrum here but also some PointTSpectrum */
	  } else if ((in->modType==OBIT_SkyModel_PointModTSpec) && (modType==OBIT_SkyModel_GaussModTSpec)) {
	    table[7+nadd] = gp1;
	    table[8+nadd] = gp2;
	    table[9+nadd] = gp3;
	    /*  spectrum */
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+3] = array[iterm+toff]*specCorr[iterm];
	  
	    /* Only Point here but some with spectrum - zero spectra (Unlikely) */
	  } else if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_PointModTSpec)) {
	    for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+3] = 0.0;
	    
	    /* Only PointSpectrum here but some GaussianSpectrum - zero Gauss comps */
	  } else if ((modType==OBIT_SkyModel_PointModSpec) && (in->modType==OBIT_SkyModel_GaussModSpec)) {
	    table[7+nadd] = 0.0;
	    table[8+nadd] = 0.0;
	    table[9+nadd] = 0.0;
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+3] = array[iterm+toff];
	    
	    /* GaussianSpectrum here but also some PointSpectrum */
	  } else if ((in->modType==OBIT_SkyModel_PointModSpec) && (modType==OBIT_SkyModel_GaussModSpec)) {
	    table[7+nadd] = gp1;
	    table[8+nadd] = gp2;
	    table[9+nadd] = gp3;
	    /*  spectrum */
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+3] = array[iterm+toff];
	    
	    /* Only Point here but some with spectrum - zero spectra (Unlikely) */
	  } else if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_PointModSpec)) {
	    for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+3] = 0.0;
	    
	  } else { /* Unsupported combination */
	    Obit_log_error(err, OBIT_Error,"%s Unsupported combination of model types %d %d  %s",
			   routine, modType, in->modType, CCTable->name);
	    Obit_traceback_val (err, routine, in->name, retCode);
	  }
	} /* end mixed type */
	    
	/* Update */
	table += lrec;
	ncomp++;
      } /* End only desired */
      array += warray;
    } /* end loop over components */

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
    table[0] = -10.0;
    table[1] = 0.0;
    table[2] = 0.0;
    table[3] = 0.0;
    table[4] = 0.0;
    table[5] = 0.0;
    table[6] = 0.0;
    table += lrec;  /* Update pointer */
  } /* end loop zeroing extra components */

  in->numComp = ncomp;

  if (specCorr) g_free(specCorr); /* Cleanup */

  /* Find anything */
  gotSome = ncomp>0;

  return gotSome;
} /* end ObitSkyModelVMBeamMFLoadComps */

/**
 * Grid components model into in plane member and Fourier transform to
 * FTplanes and apply Gaussian taper if needed.
 * Multiplies by factor member and any prior spectral index correction.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * \param inn  SkyModelVMBeamMF 
 * \param n   Image number on mosaic, 0-rel
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelVMBeamMFGridComps (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
					ObitErr *err)
{
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
  gboolean gotSome = FALSE;
  gchar *routine = "ObitSkyModelVMBeamMFGridComps";
  
  /* error checks */
  if (err->error) return gotSome;
  if ((n<0) || (n>in->mosaic->numberImages-1)) {
    Obit_log_error(err, OBIT_Error,"%s requested field %d out of range [0,%d]",
		   routine, n, in->mosaic->numberImages-1);
      return gotSome;
  }

  /* Load/FT Grid CC table */
  gotSome = ObitSkyModelVMBeamMFGridFTComps (inn, n, uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  return gotSome;
} /* end ObitSkyModeVMBeamlMFGridComps */

/**
 * NYI
 * Load image model into in plane member and Fourier transform.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * \param inn SkyModelVMBeamMF 
 * \param n   Image number on mosaic
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model
 */
gboolean ObitSkyModelVMBeamMFLoadImage (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
					ObitErr *err)
{
  g_error ("ObitSkyModelVMBeamMFLoadImage not implemented");
  return FALSE;
} /* end ObitSkyModelVMBeamMFLoadImage */

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
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)inn;
  VMBeamMFFTFuncArg *args;
  ofloat *ddata;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelVMBeamMFFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* Check */
  args = (VMBeamMFFTFuncArg*)in->threadArgs[0];
  if ((strlen(args->type)>8) || (strncmp(args->type, "mfvmbeam", 8))) {
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
    args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
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
}  /* end ObitSkyModelVMBeamMFFTDFT */

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
 * \li dimGain Dimension of Rgain
 * \li Rgain   Float array of time/spatially variable R component gain
 * \li Lgain   Float array of time/spatially variable L component gain
 * \li Qgain   Float array of time/spatially variable Q component gain
 * \li Ugain   Float array of time/spatially variable U component gain
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
  ofloat *Rgain    = largs->Rgain;
  ofloat *Lgain    = largs->Lgain;
  ofloat *Qgain    = largs->Qgain;
  ofloat *Ugain    = largs->Ugain;

  olong iVis=0, iIF, ifq, iChannel, iStoke, iComp=0, lcomp, ncomp;
  olong lrec, nrparm, naxis[2], channel, plane, sVis, eVis;
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF, lim;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVIF, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, ilocb, suba, itemp, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale, oldPB, newPB;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealQ,  sumImagQ,  modRealQ=0.0,  modImagQ=0.0;
  ofloat sumRealU,  sumImagU,  modRealU=0.0,  modImagU=0.0;
  ofloat cbase, *rgain1, *lgain1, *rgain2, *lgain2, tx, ty, tz, ll, lll;
  ofloat *qgain1, *ugain1, *qgain2, *ugain2, re, im, specFact;
  ofloat amp, ampr, ampl, arg, freq2=0.0,freqFact, wtRR=0.0, wtLL=0.0, temp;
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArrR[FazArrSize+5], AmpArrL[FazArrSize+5], FazArr[FazArrSize+5];
  ofloat CosArr[FazArrSize+5], SinArr[FazArrSize+5];
  ofloat ExpArg[FazArrSize],  ExpVal[FazArrSize];
  olong it, jt, itcnt, iSpec, nSpec, itab;
  olong sIF, sChannel;
  gboolean doCrossPol, updatePB, reGain;
  odouble *freqArr, SMRefFreq, specFreqFact;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMBeamMFFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;
  ilocb =  uvdata->myDesc->ilocb;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  nSpec         = in->nSpec - 1;  /* Offset in comp array using TSpec */
  /* PB corr should be mostly turned off higher up */
  startIF       = in->startIF-1;
  numberIF      = MAX (1, in->numberIF);
  jincif        = uvdata->myDesc->incif;
  nUVIF         = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocif];
  startChannel  = in->startChannel-1;
  numberChannel = MAX (1, in->numberChannel);
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
  /* Inverse of Sky model reference frequency */
  if (in->mosaic!=NULL) 
    SMRefFreq = 1.0 / in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocf];
  else 
    SMRefFreq = 1.0 / uvdata->myDesc->freq;

  /* Current channel (0-rel) */
  channel = 0;
  plane   = in->FreqPlane[largs->channel];  /* Which plane in correction cube */
  oldPB   = -1.0;   /* Initial primary beam correction */

  /* Set component gain lists by antenna and type 
     assume all antennas the same */
  ddata = Data;
  rgain1 = Rgain;
  rgain2 = Rgain;
  lgain1 = Lgain;
  lgain2 = Lgain;
  qgain1 = Qgain;
  qgain2 = Qgain;
  ugain1 = Ugain;
  ugain2 = Ugain;
  
  lrec    = uvdata->myDesc->lrec;         /* Length of record */
  visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
  nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */

  /* Outer loop over blocks of channels */
  /* Starting parameters this pass */
  lstartIF       = startIF;
  lstartChannel  = startChannel;

  sVis     = loVis;
  eVis     = hiVis;
  iVis     = 0;
  sIF      = lstartIF;
  sChannel = lstartChannel;

  /**************************** DEBUG **********************************/
#ifdef DEBUG 
  /* Replace model with single component and set to replace the data with the model */
  in->doReplace = TRUE;
  mcomp         = 1;      /* One component */
  in->modType   = OBIT_SkyModel_GaussModTSpec;  /* Model type */
  ddata = Data;  /* replace model */
  ddata[0] = ddata[0];  /* Field number */
  ddata[1] = 0.0;  /* X offset from center (deg) */
  ddata[2] = 0.0;  /* Y offset from center (deg) */
  amp = 1.0;    /* Amplitude all channels */
  for (iSpec=0; iSpec<=nSpec; iSpec++) ddata[3 + iSpec] = amp;
  /* Position - no rotation - no field offset */
  ddata[4+nSpec] = ddata[1] * DG2RAD * 2.0 * G_PI;
  ddata[5+nSpec] = ddata[2] * DG2RAD * 2.0 * G_PI;
  ddata[6+nSpec] = 0.0;
  /* Gaussian - only round 60" */
  ddata[7+nSpec] = (20.0/3600.0) * DG2RAD * (G_PI / 1.17741022) * sqrt (0.5);
  ddata[7+nSpec] = -ddata[7+nSpec] * ddata[7+nSpec];
  ddata[8+nSpec] = ddata[7+nSpec];
  ddata[9+nSpec] = 0.0;
#endif /* DEBUG */ 
  /************************ end DEBUG **********************************/
  
  /* Loop over IFs */
  channel = lstartIF * nUVchan + lstartChannel; /* UV Channel */
  for (iIF=sIF; iIF<lstartIF+numberIF; iIF++) {
    offsetIF = nrparm + iIF*jincif; 
    /* Loop over channels */
    for (iChannel=sChannel; iChannel<lstartChannel+numberChannel; iChannel++) {
      offsetChannel = offsetIF + iChannel*jincf; 
      ifq      = MIN (nUVchan*nUVIF, MAX (0, iIF*kincif + iChannel*kincf));
      channel  = ifq;  
      iSpec    = MIN (nSpec, MAX (0,in->specIndex[ifq]));
      freqFact = fscale[ifq];  /* Frequency scaling factor */
      freq2    = freqFact*freqFact;    /* Frequency factor squared */
      specFreqFact = freqArr[ifq] * SMRefFreq;
      
      /* New PB correction? */
      updatePB = needNewPB (in, uvdata, iIF, iChannel, oldPB, &newPB, err);
      /* New plane in beam image? or PB update? */
      if (updatePB || (plane!=in->FreqPlane[MIN(channel, (in->numUVChann-1))])) {
	oldPB = newPB;
	plane   = in->FreqPlane[MIN(channel, (in->numUVChann-1))];  /* Which plane in correction cube */
	largs->channel  = ifq;
	largs->BeamFreq = freqArr[ifq];
	/* Subarray 0-rel */
	itemp = (olong)visData[ilocb];
	suba = 100.0 * (visData[ilocb]-itemp) + 0.5; 
	/* Update antenna gains */
	myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba, uvdata, ithread, err);
	/* DEBUG
	ObitThreadLock(in->thread);
	fprintf (stderr,"ch %d IF %d plane %d ifq %d rgain %f lgain %f\n",
		 iChannel, iIF, plane, ifq, Rgain[0], Lgain[0]);
	ObitThreadUnlock(in->thread);  */
	/* DEBUG */
      } /* end new plane */
      if (err->error) {  /* Error? */
	ObitThreadLock(in->thread);  /* Lock against other threads */
	Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
		       routine);
	ObitThreadUnlock(in->thread); 
	goto finish;
      }
      
      /* Loop over vis  */
      reGain = FALSE;
      for (iVis=sVis; iVis<eVis; iVis++) {
	visData = uvdata->buffer+iVis*lrec;    /* Buffer pointer with appropriate offset */
	/* Exceed current time validity for gains */
	if (visData[iloct] > largs->endVMModelTime) {
	  /* Subarray 0-rel */
	  itemp = (olong)visData[ilocb];
	  suba = 100.0 * (visData[ilocb]-itemp) + 0.5; 
	  /* Update */
	  reGain = TRUE;
	  myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba, uvdata, ithread, err);
	  if (err->error) {
	    ObitThreadLock(in->thread);  /* Lock against other threads */
	    Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
			   routine);
	    ObitThreadUnlock(in->thread); 
	    goto finish;
	  }
	} /* end update */
	
	  /* Need antennas numbers */
	cbase = visData[ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	ant1--;    /* 0 rel */
	ant2--;    /* 0 rel */
	
	/* Sum over components */
	/* Table values 3...=Amp, 4+nSpec=-2*pi*x, 5+nSpec =-2*pi*y, 6+nSpec=-2*pi*z */
	sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	sumRealQ  = sumImagQ  = sumRealU  = sumImagU  = 0.0;
	ddata = Data;
	
	/* Sum by model type - assume phase same for RR, LL */
	switch (in->modType) {
	case OBIT_SkyModel_PointMod:     /* Point */
	  /* From the AIPSish QXXPTS.FOR  */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    for (iComp=it; iComp<mcomp; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					  ddata[5]*visData[ilocv] + 
					  ddata[6]*visData[ilocw]);
	      /* Parallel  pol */
	      /* Amplitude from component flux and two gains */
	      AmpArrR[itcnt] = ddata[3] * rgain1[iComp];
	      AmpArrL[itcnt] = ddata[3] * lgain1[iComp];
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
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      if (ddata[3]!=0.0) {  /* valid? */
		tx = ddata[4]*(odouble)visData[ilocu];
		ty = ddata[5]*(odouble)visData[ilocv];
		tz = ddata[6]*(odouble)visData[ilocw];
		/* Frequency dependent term */
		lll = ll = log(specFreqFact);
		arg = 0.0;
		for (iterm=0; iterm<nterm; iterm++) {
		  arg += ddata[6+iterm] * lll;
		  lll *= ll;
		}
		specFact = exp(arg);
		FazArr[itcnt]  = freqFact * (tx + ty + tz);
		AmpArrR[itcnt] = specFact * ddata[3] * rgain1[iComp];
		AmpArrL[itcnt] = specFact * ddata[3] * lgain1[iComp];
		itcnt++;          /* Count in amp/phase buffers */
	      }  /* end if valid */
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
	case OBIT_SkyModel_PointModTSpec:     /* Point + tabulated spectrum */
	  itab = 3 + iSpec;
	  for (it=0; it<mcomp; it+=FazArrSize) {  /* Outer component loop */
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      if (ddata[itab]!=0.0) {  /* valid? */
		tx = ddata[4+nSpec]*(odouble)visData[ilocu];
		ty = ddata[5+nSpec]*(odouble)visData[ilocv];
		tz = ddata[6+nSpec]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArrR[itcnt] = ddata[itab] * rgain1[iComp];
		AmpArrL[itcnt] = ddata[itab] * lgain1[iComp];
		itcnt++;          /* Count in amp/phase buffers */
	      }  /* end if valid */
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
	  } /* End outer component loop */
	  break;
	case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	  /* From the AIPSish QGASUB.FOR  */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					  ddata[5]*visData[ilocv] + 
					  ddata[6]*visData[ilocw]);
	      /* Parallel  pol */
	      arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			     ddata[8]*visData[ilocv]*visData[ilocv] +
			     ddata[9]*visData[ilocu]*visData[ilocv]);
	      ampr = ddata[3] * rgain1[iComp];
	      ampl = ddata[3] * lgain1[iComp];
	      ExpArg[itcnt]  = arg;
	      AmpArrR[itcnt] = ampr;
	      AmpArrL[itcnt] = ampl;
	      itcnt++;          /* Count in amp/phase buffers */
	    }  /* end inner loop over components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Convert Gaussian exp arguments */
	    ObitExpVec(itcnt, ExpArg, ExpVal);
	    
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      amp = AmpArrR[jt]*ExpVal[jt];
	      sumRealRR += amp*CosArr[jt];
	      sumImagRR += amp*SinArr[jt];
	      amp = AmpArrL[jt]*ExpVal[jt];
	      sumRealLL += amp*CosArr[jt];
	      sumImagLL += amp*SinArr[jt];
	      if (doCrossPol) {
		/* Cross pol */
		re = 0.5*(AmpArrR[it]+AmpArrL[jt])*CosArr[jt]*ExpVal[jt];
		im = 0.5*(AmpArrR[it]+AmpArrL[jt])*SinArr[jt]*ExpVal[jt];
		sumRealQ += re * qgain1[jt];
		sumImagQ += im * qgain1[jt];
		sumRealU += re * ugain1[jt];
		sumImagU += im * ugain1[jt];
	      } /* end xpol */
	    } /* End loop over amp/phase buffer */
	  } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      if (ddata[3]!=0.0) {  /* valid? */
		/* Frequency dependent term */
		lll = ll = log(specFreqFact);
		arg = 0.0;
		for (iterm=0; iterm<nterm; iterm++) {
		  arg += ddata[9+iterm] * lll;
		  lll *= ll;
		}
		specFact = exp(arg);
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			       ddata[8]*visData[ilocv]*visData[ilocv] +
			       ddata[9]*visData[ilocu]*visData[ilocv]);
		amp = specFact * ddata[3];
		tx = ddata[4]*(odouble)visData[ilocu];
		ty = ddata[5]*(odouble)visData[ilocv];
		tz = ddata[6]*(odouble)visData[ilocw];
		ExpArg[itcnt] = arg;
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		AmpArrR[itcnt] = amp * rgain1[iComp];
		AmpArrL[itcnt] = amp * lgain1[iComp];
		itcnt++;          /* Count in amp/phase buffers */
	      }  /* end if valid */
	    }  /* end inner loop over components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Convert Gaussian exp arguments */
	    ObitExpVec(itcnt, ExpArg, ExpVal);
	    
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      amp = AmpArrR[jt]*ExpVal[jt];
	      sumRealRR += amp*CosArr[jt];
	      sumImagRR += amp*SinArr[jt];
	      amp = AmpArrL[jt]*ExpVal[jt];
	      sumRealLL += amp*CosArr[jt];
	      sumImagLL += amp*SinArr[jt];
	      if (doCrossPol) {
		/* Cross pol */
		re = 0.5*(AmpArrR[jt]+AmpArrL[jt])*CosArr[jt]*ExpVal[jt];
		im = 0.5*(AmpArrR[jt]+AmpArrL[jt])*SinArr[jt]*ExpVal[jt];
		sumRealQ += re * qgain1[jt];
		sumImagQ += im * qgain1[jt];
		sumRealU += re * ugain1[jt];
		sumImagU += im * ugain1[jt];
	      } /* end xpol */
	    } /* End loop over amp/phase buffer */
	  } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_GaussModTSpec:     /* Gaussian on sky + tabulated spectrum*/
	  itab = 3 + iSpec;
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      if (ddata[itab]!=0.0) {  /* valid? */
		arg = freq2 * (ddata[7+nSpec]*visData[ilocu]*visData[ilocu] +
			       ddata[8+nSpec]*visData[ilocv]*visData[ilocv] +
			       ddata[9+nSpec]*visData[ilocu]*visData[ilocv]);
		amp = ddata[itab];
		tx = ddata[4+nSpec]*(odouble)visData[ilocu];
		ty = ddata[5+nSpec]*(odouble)visData[ilocv];
		tz = ddata[6+nSpec]*(odouble)visData[ilocw];
		ExpArg[itcnt] = arg;
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArrR[itcnt] = amp * rgain1[iComp];
		AmpArrL[itcnt] = amp * lgain1[iComp];
		itcnt++;          /* Count in amp/phase buffers */
	      }  /* end if valid */
	    }  /* end inner loop over components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Convert Gaussian exp arguments */
	    ObitExpVec(itcnt, ExpArg, ExpVal);
	    
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      /* Parallel  pol */
	      amp = AmpArrR[jt]*ExpVal[jt];
	      sumRealRR += amp*CosArr[jt];
	      sumImagRR += amp*SinArr[jt];
	      amp = AmpArrL[jt]*ExpVal[jt];
	      sumRealLL += amp*CosArr[jt];
	      sumImagLL += amp*SinArr[jt];
	      if (doCrossPol) {
		/* Cross pol */
		re = 0.5*(AmpArrR[jt]+AmpArrL[jt])*CosArr[jt]*ExpVal[jt];
		im = 0.5*(AmpArrR[jt]+AmpArrL[jt])*SinArr[jt]*ExpVal[jt];
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
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					  ddata[5]*visData[ilocv] + 
					  ddata[6]*visData[ilocw]);
	      arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				    visData[ilocv]*visData[ilocv]) * ddata[7];
	      arg = MAX (arg, 0.1);
	      arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
	      AmpArrR[itcnt] = ddata[3] * rgain1[iComp] * arg;
	      AmpArrL[itcnt] = ddata[3] * lgain1[iComp] * arg;

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
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      if (ddata[3]!=0.0) {  /* valid? */
		/* Frequency dependent term */
		lll = ll = log(specFreqFact);
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
		itcnt++;          /* Count in amp/phase buffers */
	      } /* end if valid */
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
	case OBIT_SkyModel_USphereModTSpec:    /* Uniform sphere + tabulated spectrum*/
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    lim = MIN (mcomp, it+FazArrSize);
	    for (iComp=it; iComp<lim; iComp++) {
	      itab = 3 + iSpec;
	      ddata = Data + iComp*lcomp;  /* Component pointer */
	      if (ddata[itab]!=0.0) {  /* valid? */
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ddata[7+nSpec];
		arg = MAX (arg, 0.1);
		amp = ddata[itab] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		tx = ddata[4+nSpec]*(odouble)visData[ilocu];
		ty = ddata[5+nSpec]*(odouble)visData[ilocv];
		tz = ddata[6+nSpec]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArrR[itcnt] = amp * rgain1[iComp];
		AmpArrL[itcnt] = amp * lgain1[iComp];
		itcnt++;          /* Count in amp/phase buffers */
	      } /* end if valid */
	    }  /* end inner loop over components */
	    
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
	  } /* End outer componentloop */
	  break;
	default:
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Unknown Comp model type %d in %s",
			 routine, in->modType, in->name);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}; /* end switch by model type */

	/* DEBUG 
	if ((fabs(modRealRR)>1.0)||(fabs(modRealLL)>1.0)) {
	  ObitThreadLock(in->thread);
	  fprintf (stderr,"iVis %d ch %d IF %d Model %f %f %f %f\n",
		   iVis, iChannel, iIF, sumRealRR, sumImagRR, sumRealLL, sumImagLL);
	  ObitThreadUnlock(in->thread); 
	}*/
	/* DEBUG */
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
	
      } /* end vis loop */
      if (reGain) largs->endVMModelTime = -1.0e20;  /* Reset gains? */
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
static gpointer ThreadSkyModelVMBeamMFFTDFTPh (gpointer args)
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
  ofloat *Rgain    = largs->Rgain;
  ofloat *Lgain    = largs->Lgain;
  ofloat *Qgain    = largs->Qgain;
  ofloat *Ugain    = largs->Ugain;
  ofloat *Rgaini   = largs->Rgaini;
  ofloat *Lgaini   = largs->Lgaini;
  ofloat *Qgaini   = largs->Qgaini;
  ofloat *Ugaini   = largs->Ugaini;

  olong iVis=0, iIF, ifq, iChannel, iStoke, iComp, lcomp, ncomp;
  olong lrec, nrparm, naxis[2], channel, plane, sVis, eVis;
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF, lim;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVIF, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, ilocb, suba, itemp, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale, oldPB, newPB;
  ofloat sumRealRR, sumImagRR, modRealRR=0.0, modImagRR=0.0;
  ofloat sumRealLL, sumImagLL, modRealLL=0.0, modImagLL=0.0;
  ofloat sumRealQ,  sumImagQ,  modRealQ=0.0,  modImagQ=0.0;
  ofloat sumRealU,  sumImagU,  modRealU=0.0,  modImagU=0.0;
  ofloat *rgain1, *lgain1, *rgain2, *lgain2, *rgain1i, *lgain1i, *rgain2i, *lgain2i;
  ofloat *qgain1, *ugain1, *qgain2, *ugain2, *qgain1i, *ugain1i, *qgain2i, *ugain2i;
  ofloat cbase, tx, ty, tz, ll, lll;
  ofloat re, im, specFact;
  ofloat amp, ampr, ampl, arg, freq2=0.0,freqFact, wtRR=0.0, wtLL=0.0, temp;
#ifndef FazArrSize
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
#endif /* FazArrSize */ 
  ofloat AmpArrRr[FazArrSize], AmpArrLr[FazArrSize], AmpArrRi[FazArrSize], AmpArrLi[FazArrSize];
  ofloat FazArr[FazArrSize];
  ofloat CosArr[FazArrSize], SinArr[FazArrSize];
  ofloat ExpArg[FazArrSize],  ExpVal[FazArrSize];
  olong it, jt, itcnt, iSpec, nSpec, itab;
  gboolean doCrossPol, updatePB, reGain;
  odouble *freqArr, SMRefFreq, specFreqFact;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMBeamMFFTDFTPh";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;
  ilocb =  uvdata->myDesc->ilocb;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  nSpec         = in->nSpec - 1;  /* Offset in comp array using TSpec */
  /* PB corr should be mostly turned off higher up */
  startIF       = in->startIF-1;
  numberIF      = MAX (1, in->numberIF);
  jincif        = uvdata->myDesc->incif;
  nUVIF         = uvdata->myDesc->inaxes[ uvdata->myDesc->jlocif];
  startChannel  = in->startChannel-1;
  numberChannel = MAX (1, in->numberChannel);
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
  /* Inverse of Sky model reference frequency */
  if (in->mosaic!=NULL) 
    SMRefFreq = 1.0 / in->mosaic->images[0]->myDesc->crval[in->mosaic->images[0]->myDesc->jlocf];
  else 
    SMRefFreq = 1.0 / uvdata->myDesc->freq;

  /* Current channel (0-rel) */
  channel = 0;
  plane   = in->FreqPlane[largs->channel];  /* Which plane in correction cube */
  oldPB   = -1.0;   /* Initial primary beam correction */

  /* Set component gain lists by antenna and type 
     assume all antennas the same */
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
  
  lrec    = uvdata->myDesc->lrec;         /* Length of record */
  visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
  nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */

  /* Outer loop over blocks of channels */
  /* Starting parameters this pass */
  lstartIF       = startIF;
  lstartChannel  = startChannel;

  /* Outer loop over vis */
  sVis = loVis;
  eVis = hiVis;
  iVis = 0;
  while (sVis<eVis) {
    
    /* Loop over IFs */
    channel = lstartIF* nUVchan + lstartChannel; /* UV Channel */
    for (iIF=lstartIF; iIF<lstartIF+numberIF; iIF++) {
      offsetIF = nrparm + iIF*jincif; 

      /* Loop over channels */
      for (iChannel=lstartChannel; iChannel<lstartChannel+numberChannel; iChannel++) {
	offsetChannel = offsetIF + iChannel*jincf; 
	ifq      = MIN (nUVchan*nUVIF, MAX (0, iIF*kincif + iChannel*kincf));
	iSpec    = MIN (nSpec, MAX (0,in->specIndex[ifq]));
	freqFact = fscale[ifq];  /* Frequency scaling factor */
	freq2    = freqFact*freqFact;    /* Frequency factor squared */
	specFreqFact = freqArr[iIF*kincif + iChannel*kincf] * SMRefFreq;

	/* New PB correction? */
	updatePB = needNewPB (in, uvdata, iIF, iChannel, oldPB, &newPB, err);
	/* New plane in beam image? or PB update? */
	if (updatePB || (plane!=in->FreqPlane[MIN(channel, (in->numUVChann-1))])) {
	  oldPB = newPB;
	  plane   = in->FreqPlane[MIN(channel, (in->numUVChann-1))];  /* Which plane in correction cube */
	  largs->channel  = ifq;
	  largs->BeamFreq = freqArr[ifq];
	  /* Subarray 0-rel */
	  itemp = (olong)visData[ilocb];
	  suba = 100.0 * (visData[ilocb]-itemp) + 0.5; 
	  /* Update antenna gains */
	  /* Update */
	  myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba, uvdata, ithread, err);
	} /* end new plane */
	if (err->error) {  /* Error? */
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
			 routine);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}
	
	/* Loop over vis  */
	reGain = FALSE;
	for (iVis=sVis; iVis<eVis; iVis++) {
	  visData = uvdata->buffer+iVis*lrec;    /* Buffer pointer with appropriate offset */
	  /* Exceed current time validity for gains */
	  if (visData[iloct] > largs->endVMModelTime) {
	    /* Subarray 0-rel */
	    itemp = (olong)visData[ilocb];
	    suba = 100.0 * (visData[ilocb]-itemp) + 0.5; 
	    /* Update */
	    reGain = TRUE;
	    myClass->ObitSkyModelVMUpdateModel ((ObitSkyModelVM*)in, visData[iloct], suba, uvdata, ithread, err);
	    if (err->error) {
	      ObitThreadLock(in->thread);  /* Lock against other threads */
	      Obit_log_error(err, OBIT_Error,"%s Error updating VMComps",
			     routine);
	      ObitThreadUnlock(in->thread); 
	      goto finish;
	    }
	  } /* end update */
      
	  /* Need antennas numbers */
	  cbase = visData[ilocb]; /* Baseline */
	  ant1 = (cbase / 256.0) + 0.001;
	  ant2 = (cbase - ant1 * 256) + 0.001;
	  ant1--;    /* 0 rel */
	  ant2--;    /* 0 rel */
	  
	  /* Sum over components */
	  /* Table values 3...=Amp, 4+nSpec=-2*pi*x, 5+nSpec =-2*pi*y, 6+nSpec=-2*pi*z */
	  sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	  sumRealQ  = sumImagQ  = sumRealU  = sumImagU  = 0.0;
	  ddata = Data;
	  
	  /* Sum by model type - assume phase same for RR, LL */
	  switch (in->modType) {
	  case OBIT_SkyModel_PointMod:     /* Point */
	    /* From the AIPSish QXXPTS.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		/* Amplitude from component flux and two gains */
		AmpArrRr[itcnt] = ddata[3] * rgain1[iComp];
		AmpArrLr[itcnt] = ddata[3] * lgain1[iComp];
		AmpArrRi[itcnt] = ddata[3] * rgain1i[iComp];
		AmpArrLi[itcnt] = ddata[3] * lgain1i[iComp];
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
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[3]!=0.0) {  /* valid? */
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  /* Frequency dependent term */
		  lll = ll = log(specFreqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[6+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  amp = ddata[3] * specFact;
		  FazArr[itcnt]  = freqFact * (tx + ty + tz);
		  AmpArrRr[itcnt] = amp * rgain1[iComp];
		  AmpArrLr[itcnt] = amp * lgain1[iComp];
		  AmpArrRi[itcnt] = amp * rgain1i[iComp];
		  AmpArrLi[itcnt] = amp * lgain1i[iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
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
	  case OBIT_SkyModel_PointModTSpec:     /* Point + tabulated spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {  /* Outer component loop */
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		itab = 3 + iSpec;
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[itab]!=0.0) {  /* valid? */
		  tx = ddata[4+nSpec]*(odouble)visData[ilocu];
		  ty = ddata[5+nSpec]*(odouble)visData[ilocv];
		  tz = ddata[6+nSpec]*(odouble)visData[ilocw];
		  FazArr[itcnt] = freqFact * (tx + ty + tz);
		  /* Parallel  pol */
		  /* Amplitude from component flux and two gains */
		  amp = ddata[itab];
		  AmpArrRr[itcnt] = amp * rgain1[iComp];
		  AmpArrLr[itcnt] = amp * lgain1[iComp];
		  AmpArrRi[itcnt] = amp * rgain1i[iComp];
		  AmpArrLi[itcnt] = amp * lgain1i[iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		}  /* end if valid */
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
	    } /* End outer componentloop */
	    break;
	  case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	    /* From the AIPSish QGASUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		/* Parallel  pol */
		arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
				ddata[8]*visData[ilocv]*visData[ilocv] +
				ddata[9]*visData[ilocu]*visData[ilocv]);
		ExpArg[itcnt]  = arg;
		amp  = ddata[3];
		ampr = amp * rgain1[iComp];
		ampl = amp * lgain1[iComp];
		AmpArrRr[itcnt] = ampr;
		AmpArrLr[itcnt] = ampl;

		ampl = amp * lgain1i[iComp];
		ampr = amp * rgain1i[iComp];
		AmpArrRi[itcnt] = ampr;
		AmpArrLi[itcnt] = ampl;
		itcnt++;          /* Count in amp/phase buffers */
	      }  /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Convert Gaussian exp arguments */
	      ObitExpVec(itcnt, ExpArg, ExpVal);

	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt]*ExpVal[jt] - AmpArrRi[jt]*SinArr[jt]*ExpVal[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt]*ExpVal[jt] + AmpArrRi[jt]*CosArr[jt]*ExpVal[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt]*ExpVal[jt] - AmpArrLi[jt]*SinArr[jt]*ExpVal[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt]*ExpVal[jt] + AmpArrLi[jt]*CosArr[jt]*ExpVal[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  re *= ExpVal[jt];
		  im *= ExpVal[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(specFreqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += ddata[9+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
				 ddata[8]*visData[ilocv]*visData[ilocv] +
				 ddata[9]*visData[ilocu]*visData[ilocv]);
		  amp = specFact * ddata[3];
		  tx = ddata[4]*(odouble)visData[ilocu];
		  ty = ddata[5]*(odouble)visData[ilocv];
		  tz = ddata[6]*(odouble)visData[ilocw];
		  ExpArg[itcnt]   = arg;
		  FazArr[itcnt]   = freqFact * (tx + ty + tz);
		  AmpArrRr[itcnt] = amp * rgain1[iComp];
		  AmpArrLr[itcnt] = amp * lgain1[iComp];
		  AmpArrRi[itcnt] = amp * rgain1i[iComp];
		  AmpArrLi[itcnt] = amp * lgain1i[iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
	      }  /* end inner loop over components */

	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Convert Gaussian exp arguments */
	      ObitExpVec(itcnt, ExpArg, ExpVal);

	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt]*ExpVal[jt] - AmpArrRi[jt]*SinArr[jt]*ExpVal[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt]*ExpVal[jt] + AmpArrRi[jt]*CosArr[jt]*ExpVal[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt]*ExpVal[jt] - AmpArrLi[jt]*SinArr[jt]*ExpVal[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt]*ExpVal[jt] + AmpArrLi[jt]*CosArr[jt]*ExpVal[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  re *= ExpVal[jt];
		  im *= ExpVal[jt];
		  sumRealQ += re * qgain1[jt]  - im * qgain1i[jt];
		  sumImagQ += re * qgain1i[jt] + im * qgain1[jt];
		  sumRealU += re * ugain1[jt]  - im * ugain1i[jt];
		  sumImagU += re * ugain1i[jt] + im * ugain1[jt];
		} /* end xpol */
	      } /* End loop over amp/phase buffer */
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModTSpec:     /* Gaussian on sky + tabulated spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		itab = 3 + iSpec;
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[itab]!=0.0) {  /* valid? */
		  arg = freq2 * (ddata[7+nSpec]*visData[ilocu]*visData[ilocu] +
				 ddata[8+nSpec]*visData[ilocv]*visData[ilocv] +
				 ddata[9+nSpec]*visData[ilocu]*visData[ilocv]);
		  amp = ddata[itab];
		  tx = ddata[4+nSpec]*(odouble)visData[ilocu];
		  ty = ddata[5+nSpec]*(odouble)visData[ilocv];
		  tz = ddata[6+nSpec]*(odouble)visData[ilocw];
		  ExpArg[itcnt]  = arg;
		  FazArr[itcnt] = freqFact * (tx + ty + tz);
		  /* Parallel  pol */
		  /* Amplitude from component flux and two gains */
		  AmpArrRr[itcnt] = amp * rgain1[iComp];
		  AmpArrLr[itcnt] = amp * lgain1[iComp];
		  AmpArrRi[itcnt] = amp * rgain1i[iComp];
		  AmpArrLi[itcnt] = amp * lgain1i[iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
	      }  /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	      /* Convert Gaussian exp arguments */
	      ObitExpVec(itcnt, ExpArg, ExpVal);

	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealRR += AmpArrRr[jt]*CosArr[jt]*ExpVal[jt] - AmpArrRi[jt]*SinArr[jt]*ExpVal[jt];
		sumImagRR += AmpArrRr[jt]*SinArr[jt]*ExpVal[jt] + AmpArrRi[jt]*CosArr[jt]*ExpVal[jt];
		sumRealLL += AmpArrLr[jt]*CosArr[jt]*ExpVal[jt] - AmpArrLi[jt]*SinArr[jt]*ExpVal[jt];
		sumImagLL += AmpArrLr[jt]*SinArr[jt]*ExpVal[jt] + AmpArrLi[jt]*CosArr[jt]*ExpVal[jt];
		if (doCrossPol) {
		  /* Cross pol */
		  re = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*CosArr[jt] - 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*SinArr[jt];  
		  im = 0.5*(AmpArrRr[jt]+AmpArrLr[jt])*SinArr[jt] + 
		    0.5*(AmpArrRi[jt]+AmpArrLi[jt])*CosArr[jt];
		  re *= ExpVal[jt];
		  im *= ExpVal[jt];
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
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		FazArr[itcnt] = freqFact * (ddata[4]*visData[ilocu] + 
					    ddata[5]*visData[ilocv] + 
					    ddata[6]*visData[ilocw]);
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ddata[7];
		arg = MAX (arg, 0.1);
		arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		amp = ddata[3] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		AmpArrRr[itcnt] = amp * rgain1[iComp];
		AmpArrLr[itcnt] = amp * lgain1[iComp];
		AmpArrRi[itcnt] = amp * rgain1i[iComp];
		AmpArrLi[itcnt] = amp * lgain1i[iComp];
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
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(specFreqFact);
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
		  AmpArrRr[itcnt] = amp * rgain1[iComp];
		  AmpArrLr[itcnt] = amp * lgain1[iComp];
		  AmpArrRi[itcnt] = amp * rgain1i[iComp];
		  AmpArrLi[itcnt] = amp * lgain1i[iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
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
	  case OBIT_SkyModel_USphereModTSpec:    /* Uniform sphere + tabulated spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		itab = 3 + iSpec;
		ddata = Data + iComp*lcomp;  /* Component pointer */
		if (ddata[itab]!=0.0) {  /* valid? */
		  arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
					visData[ilocv]*visData[ilocv]) * ddata[7+nSpec];
		  arg = MAX (arg, 0.1);
		  amp = ddata[itab] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  tx = ddata[4+nSpec]*(odouble)visData[ilocu];
		  ty = ddata[5+nSpec]*(odouble)visData[ilocv];
		  tz = ddata[6+nSpec]*(odouble)visData[ilocw];
		  FazArr[itcnt] = freqFact * (tx + ty + tz);
		  /* Parallel  pol */
		  /* Amplitude from component flux and two gains */
		  AmpArrRr[itcnt] = amp * rgain1[iComp];
		  AmpArrLr[itcnt] = amp * lgain1[iComp];
		  AmpArrRi[itcnt] = amp * rgain1i[iComp];
		  AmpArrLi[itcnt] = amp * lgain1i[iComp];
		  itcnt++;          /* Count in amp/phase buffers */
		} /* end if valid */
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
	    } /* End outer componentloop */
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

	} /* end inner vis loop */
	if (reGain) largs->endVMModelTime = -1.0e20;  /* Reset gains? */
      } /* end channel loop */
    } /* end IF loop */
    goto finish;  /* Do not need outer loop */
  } /* End outer vis loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadSkyModelVMBeamMFFTDFTPh */

/**
 * Do Fourier transform using the a gridded image or set of components 
 * for a buffer of data.
 * If threading has been enabled by a call to ObitThreadAllowThreads 
 * this routine will divide the buffer up amount the number of processors
 * returned by ObitThreadNumProc.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * Adapted from the AIPSish ALGSTB, QUVINT, QINTP
 * Note: Unlike AIPS, FFTw produces nontransposed images with half
 * the first (U) axis.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * \param inn    SkyModelVMBeamMF with model components loaded (ObitSkyModelVMBeamMFLoad)
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to model and subtract from current buffer
 * \param err Obit error stack object.
 */
void ObitSkyModelVMBeamMFFTGrid (ObitSkyModel *inn, olong field, ObitUV *uvdata, 
				 ObitErr *err)
{
  ObitSkyModelVMBeamMF *in  = (ObitSkyModelVMBeamMF*)inn;
  olong i, k, nvis, lovis, hivis, nvisPerThread, nThreads;
  VMBeamMFFTFuncArg *args;
  gboolean OK = TRUE, resetInterp=FALSE;
  gchar *routine = "ObitSkyModelVMBeamMFFTGrid";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array on first call */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(VMBeamMFFTFuncArg*));
    for (i=0; i<in->nThreads; i++) {
      in->threadArgs[i] = g_malloc0(sizeof(VMBeamMFFTFuncArg)); 
      args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
      strcpy (args->type, "mfvmbeam");  /* Enter type as first entry */
     }
  } /* end initialize */
  
  /* Divide up work - single threaded if too little data per call */
  nvis = uvdata->myDesc->numVisBuff;
  if (nvis<1000) nThreads = 1;
  else nThreads = in->nThreads;
  nvisPerThread = MAX (1, nvis/nThreads);
  lovis = 1;
  hivis = nvisPerThread;
  hivis = MIN (hivis, nvis);

  /* Set up thread arguments */
  for (i=0; i<nThreads; i++) {
    if (i==(nThreads-1)) hivis = nvis;  /* Make sure do all */
    args = (VMBeamMFFTFuncArg*)in->threadArgs[i];
    args->in     = inn;
    resetInterp = FALSE;
    if (args->Interps) {
      if (args->field!=field) {
	resetInterp = TRUE;
	for (k=0; k<in->nSpec; k++)
	  args->Interps[k] = ObitCInterpolateUnref(args->Interps[k]);
      }
    }
    args->field  = field;
    args->uvdata = uvdata;
    args->first  = lovis;
    args->last   = hivis;
    if (nThreads>1) args->ithread= i;
    else args->ithread = -1;
    args->err    = err;
    /* local copies of interpolators if needed */
    if ((!args->Interps) || resetInterp) {
      args->nSpec   = in->nSpec;
      args->Interps = g_malloc0(args->nSpec*sizeof(ObitCInterpolate*));
      resetInterp   = TRUE;
    }
    if (resetInterp) {
      for (k=0; k<in->nSpec; k++) {
	if (i>0) {
	  args->Interps[k] = ObitCInterpolateClone(in->myInterps[k], NULL);
	} else {
	  args->Interps[k] = ObitCInterpolateRef(in->myInterps[k]);
	}
	if (err->error) Obit_traceback_msg (err, routine, in->name);
      }
    } /* end local copy of interpolator */
    /* Update which vis */
    lovis += nvisPerThread;
    hivis += nvisPerThread;
    hivis = MIN (hivis, nvis);
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nThreads, in->GridFunc, in->threadArgs);
  
  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);

}  /* end ObitSkyModelVMBeamMFFTGrid */

/**
 * Do Fourier transform using the a gridded image or set of components 
 * for a buffer of data.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Adapted from the AIPSish ALGSTB, QUVINT, QINTP
 * Note: Unlike AIPS, FFTw produces nontransposed images with half
 * the first (U) axis.
 * Arguments are given in the structure passed as arg
 * \param arg  Pointer to FTFuncArg argument with elements
 * \li type   String identifying structure
 * \li in     SkyModelVMBeamMF with model components loaded (ObitSkyModelVMBeamMFLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number, <0-> no threads
 * \li err Obit error stack object.
 * \li Interp UV Interpolator
 * \return NULL
 */
gpointer ThreadSkyModelVMBeamMFFTGrid (gpointer args)
{
  /* Get arguments from structure */
  VMBeamMFFTFuncArg *largs = (VMBeamMFFTFuncArg*)args;
  ObitSkyModelVMBeamMF *in = (ObitSkyModelVMBeamMF*)largs->in;
  olong field        = largs->field;
  ObitUV *uvdata     = largs->uvdata;
  olong loVis        = largs->first-1;
  olong hiVis        = largs->last;
  ObitErr *err       = largs->err;
  ObitCInterpolate **Interp = largs->Interps;

  ObitImageDesc *imDesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  olong iVis, iIF, iChannel, iStoke;
  olong i, j, k, lrec, nrparm;
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF;
  olong ilocu, ilocv, ilocw, ifq, itab;
  ofloat *visData, *fscale, vis[2], flip;
  ofloat sumReal, sumImag, modReal, modImag;
  ofloat freqFact, wt=0.0, temp;
  ofloat dxyzc[3],  uvw[3], ut, vt, rt, it, fblank = ObitMagicF();
  ofloat umat[3][3], pmat[3][3], rmat[3][3], dmat[3][3];
  ofloat PC, cosPC, sinPC, konst, maprot, uvrot, ssrot, ccrot;
  odouble *freqArr;
  gboolean doRot, doConjg, isBad, do3Dmul, doPC;
  gchar *routine = "ThreadSkyModelVMBeamMFFTGrid";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Any "type" arg list allowed */

  /* Visibility pointers */
  uvDesc = uvdata->myDesc;
  ilocu =  uvDesc->ilocu;
  ilocv =  uvDesc->ilocv;
  ilocw =  uvDesc->ilocw;

  /* Set channel, IF and Stokes ranges */
  startIF  = in->startIFPB-1;
  numberIF = MAX (1, in->numberIFPB);
  jincif   = uvDesc->incif;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  jincf         = uvDesc->incf;
  startPoln  = in->startPoln-1;
  numberPoln = in->numberPoln;
  jincs      = uvDesc->incs;  /* increment in real array */
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF before freq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  /* Get pointer for frequency correction tables */
  fscale  = uvDesc->fscale;
  freqArr = uvDesc->freqArr;

  /* Field specific stuff */
  imDesc = in->mosaic->images[field]->myDesc; /* Image descriptor */
  /*  Set field center offsets. */
  maprot = ObitImageDescRotate(imDesc);
  uvrot  = ObitUVDescRotate(uvDesc);
  ssrot = sin (DG2RAD * (uvrot - maprot));
  ccrot = cos (DG2RAD * (uvrot - maprot));
  konst = DG2RAD * 2.0 * G_PI;

  /* Which way does RA go with pixel? */
  if (imDesc->cdelt[imDesc->jlocr]>0.0) flip = -1;
  else flip = 1.0;

  /* Get position phase shift parameters */
  ObitUVDescShiftPhase(uvDesc, imDesc, dxyzc, err);
  if (err->error) {
    ObitThreadLock(in->thread);  /* Lock against other threads */
    Obit_log_error(err, OBIT_Error,"%s: Error phase shifting %s",
		   routine, uvdata->name);
    ObitThreadUnlock(in->thread); 
    goto finish;
  }
  /* Phase shift for field offset? */
  doPC = (fabs(dxyzc[0])>1.0e-12) || (fabs(dxyzc[1])>1.0e-12) || 
    (fabs(dxyzc[2])>1.0e-12);
    
  /* 3D rotation matrix if needed */
  if (in->do3D) {
    do3Dmul = ObitUVDescShift3DMatrix (uvDesc, imDesc, umat, pmat);

    /* Correct field shift */
    if (doPC) {
      /* Rotation matrix for relative rotation */
      rmat[0][0] = ccrot; rmat[1][0] = ssrot; rmat[2][0] = 0.0;
      rmat[0][1] =-ssrot; rmat[1][1] = ccrot; rmat[2][1] = 0.0;
      rmat[0][2] =   0.0; rmat[1][2] =   0.0; rmat[2][2] = 1.0;
      for (i=0; i<3; i++) {
	for (j=0; j<3; j++) {
	  dmat[j][i] = 0.0;
	  for (k=0; k<3; k++) dmat[j][i] += pmat[k][i]*rmat[j][k];
	}
      }
      /* Rotate field offset XXXX*/
      ut = dxyzc[0]*dmat[0][0] + dxyzc[1]*dmat[1][0] + dxyzc[2]*dmat[2][0];
      vt = dxyzc[0]*dmat[0][1] + dxyzc[1]*dmat[1][1] + dxyzc[2]*dmat[2][1];
      wt = dxyzc[0]*dmat[0][2] + dxyzc[1]*dmat[1][2] + dxyzc[2]*dmat[2][2];
      dxyzc[0] = ut;
      dxyzc[1] = vt;
      dxyzc[2] = wt;
      /* PRJMUL (2, DDX, DMAT, DDX) */
    } /* end field shift */
  } else {do3Dmul = FALSE;}
  
  /* Rotation needed? */
  doRot = (fabs (ssrot)>1.0e-10) || (fabs (ccrot-1.0)>1.0e-4);

  /* Loop over vis in buffer */
  lrec    = uvdata->myDesc->lrec;         /* Length of record */
  visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
  nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */
  for (iVis=loVis; iVis<hiVis; iVis++) {
    /* Loop over IFs */
    for (iIF=startIF; iIF<startIF+numberIF; iIF++) {
      offsetIF = nrparm + iIF*jincif; 
      for (iChannel=startChannel; iChannel<startChannel+numberChannel; iChannel++) {
	offsetChannel = offsetIF + iChannel*jincf; 
	ifq = iIF*kincif + iChannel*kincf;
	freqFact = fscale[ifq];  /* Frequency scaling factor */
	
	/* Get u, v, w at wavelength */
	uvw[0] = freqFact * visData[ilocu];
	uvw[1] = freqFact * visData[ilocv];
	uvw[2] = freqFact * visData[ilocw];
	
	if (do3Dmul) {       /* 3D reprojection */
	  ut = (uvw[0])*umat[0][0] + (uvw[1])*umat[0][1] + (uvw[2])*umat[0][2];
	  vt = (uvw[0])*umat[1][0] + (uvw[1])*umat[1][1] + (uvw[2])*umat[1][2];
	  wt = (uvw[0])*umat[2][0] + (uvw[1])*umat[2][1] + (uvw[2])*umat[2][2];
	  uvw[0] = ut;
	  uvw[1] = vt;
	  uvw[2] = wt;
	  /* PRJMUL (1, UVW, UMAT, UVW); */
	} else if (doRot) {  /* Only rotate in u,v */
	  ut = ccrot * uvw[0] - ssrot * uvw[1];
	  vt = ccrot * uvw[1] + ssrot * uvw[0];
	  uvw[0] = ut;
	  uvw[1] = vt;
	}
	
	/* need to conjugate? (only one half U plane available) */
	doConjg = flip*uvw[0] < 0.0;
	if (doConjg) {
	  uvw[0] = -uvw[0];
	  uvw[1] = -uvw[1];
	  uvw[2] = -uvw[2];
	}
	
	/* Interpolate from UV grid */
	itab = in->specIndex[ifq];
	ObitCInterpolateOffset (Interp[itab], uvw, vis, err);
	if (err->error) {
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s: Error interpolatingFT of model",
			 routine);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}
      
	/* Blanked if outside grid  - zero data and weight */
	isBad = (vis[0]==fblank);
	
	/* Phase correction for field offset? */
	if (doPC && !isBad) {
	  PC = uvw[0]*dxyzc[0] + uvw[1]*dxyzc[1] + uvw[2]*dxyzc[2];
	  cosPC = cos(PC);
	  sinPC = sin(PC);
	  rt = cosPC * vis[0] - sinPC * vis[1];
	  it = cosPC * vis[1] + sinPC * vis[0];
	  vis[0] = rt;
	  vis[1] = it;
	}
	
	/* Conjugate? */
	if (doConjg) {
	  sumReal =  vis[0];
	  sumImag = -vis[1];
	} else {
	  sumReal =  vis[0];
	  sumImag =  vis[1];
	}
	
	/* Need to multiply model by sqrt(-1)? */
	if (in->doFlip) {
	  modReal = -sumImag;
	  modImag =  sumReal;
	} else {
	  modReal =  sumReal;
	  modImag =  sumImag;
	}
	
	/* Dividing? */
	if (in->doDivide) {
	  /* Divide model - also correct weight */
	  wt = modReal * modReal + modImag * modImag;
	  modReal /= wt;
	  modImag /= wt;
	  wt = sqrt (wt);
	}
	
	/* Stokes Loop */
	for (iStoke=startPoln; iStoke<startPoln+numberPoln; iStoke++) {
	  offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	  
	  /* Ignore blanked data */
	  if ((visData[offset+2]<=0.0) && !in->doReplace) continue;
	  
	  /* Apply model to data */
	  if (isBad) { /* Bad model (outside grid) Blank */
	    visData[offset+1] = 0.0;
	    visData[offset]   = 0.0;
	    visData[offset+2] = 0.0; /* flag weight */
	    
	  } else {   /* Model OK */
	    
	    if (in->doDivide) {
	      temp = modReal * visData[offset] + modImag * visData[offset+1];
	      visData[offset+1] = modReal * visData[offset+1] - modImag * visData[offset];
	      visData[offset]   = temp;
	      visData[offset+2] *= wt;  /* correct weight */
	    } else if (in->doReplace) {  /* replace data with model */
	      visData[offset]   = modReal;
	      visData[offset+1] = modImag;
	    } else {
	      /* Subtract model */
	      visData[offset]   -= modReal;
	      visData[offset+1] -= modImag;
	    }
	  }
	  /* Factor for next Stokes */
	  modReal *= in->stokFactor;
	  modImag *= in->stokFactor;
	  
	  offset += jincs;
	} /* end loop over Stokes */
	offsetChannel += jincf;
      } /* end loop over Channel */
      offsetIF += jincif;
    } /* end loop over IF */
    
    visData += lrec; /* Update vis pointer */
  } /* end loop over visibilities */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadSkyModelVMBeamMFFTGrid */

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
  gchar keyword[12];
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
 * A change of 1 percent is deemed necessary for an update.
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
  olong nfreq, nif, niffreq, ifreq, incf, incif;
  odouble Angle;
  ObitUVDesc *uvDesc=NULL;
  gchar *routine = "needNewPB";

  /* error checks */
  if (err->error) return update;

  *newPB = oldPB;  /* Until proven otherwise */

    /* Info from uv descriptor */
  uvDesc = uvdata->myDesc;
  nfreq = uvDesc->inaxes[uvDesc->jlocf];
  if (uvDesc->jlocif>=0) nif = uvDesc->inaxes[uvDesc->jlocif];
  else nif = 1;
  niffreq = nfreq * nif;
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
  if (fabs(PBFact-oldPB)>0.01) {
    /* Yes */
    update = TRUE;
    *newPB = PBFact;
  }

 return update;
} /* end needNewPB */
