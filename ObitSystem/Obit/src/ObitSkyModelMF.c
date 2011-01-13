/* $Id: ObitSkyModelMF.c 155 2010-02-04 13:17:17Z bill.cotton $      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010,2011                                          */
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

#include <math.h>
#include "ObitThread.h"
#include "ObitSkyModelMF.h"
#include "ObitSkyModelVMSquint.h"
#include "ObitSkyModelVMIon.h"
#include "ObitImageMosaic.h"
#include "ObitTableCCUtil.h"
#include "ObitFFT.h"
#include "ObitUVUtil.h"
#include "ObitImageUtil.h"
#include "ObitPBUtil.h"
#include "ObitMem.h"
#include "ObitSinCos.h"
#include "ObitImageMF.h"
#include "ObitBeamShape.h"
#include "ObitSpectrumFit.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelMF.c
 * ObitSkyModelMF class function definitions.
 * This class represents sky models and their Fourier transforms and is
 * derived from the ObitSkyModel base class.
 * SkyModelMF components have coarse, tabulated spectra.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSkyModelMF";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitSkyModelGetClass;

/**
 * ClassInfo structure ObitSkyModelMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSkyModelMFClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/** Over sampling factor in uv plane */
olong OverSampleMF=4; 

/*----------------- Macroes ---------------------------*/
/** Half width of gridded subtraction interpolation kernal */
#define HWIDTH 12

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSkyModelMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSkyModelMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSkyModelMFClassInfoDefFn (gpointer inClass);

/** Private: Load point model, may be overridden in derived class */
gboolean ObitSkyModelMFLoadPoint (ObitSkyModel *in, ObitUV *uvdata, ObitErr *err);

/** Private: Load  Components model, may be overridden in derived class */
gboolean ObitSkyModelMFLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);

/** Private: Threaded FTDFT */
static gpointer ThreadSkyModelMFFTDFT (gpointer arg);

/** Private: Threaded FTGrid */
static gpointer ThreadSkyModelMFFTGrid (gpointer arg);

/** Private: Primary beam and/or smoothing corrections */
ObitTableCC* ObitSkyModelMFgetPBCCTab (ObitSkyModelMF* in, ObitUV* uvdata, 
				       olong field, olong *inCCVer, olong *outCCver,
				       olong *startCC, olong *endCC, ofloat range[2],
				       ObitErr *err);
/*---------------Private structures----------------*/
/* FT threaded function argument */
typedef struct {
  /* type "MF" in this class */
  gchar type[12];
  /* SkyModelMF with model components loaded (ObitSkyModelMFLoad) */
  ObitSkyModelMF *in;
  /* Field number being processed (-1 => all) */
  olong        field;
  /* UV data set to model and subtract from current buffer */
  ObitUV       *uvdata;
  /* First (1-rel) vis in uvdata buffer to process this thread */
  olong        first;
  /* Highest (1-rel) vis in uvdata buffer to process this thread  */
  olong        last;
  /* Number of spectral channels in Interp */
  olong        nSpec;
  /* Apply prior alpha correction? */
  gboolean doAlphaCorr;
  /* Prior spectral index correction */
  ofloat priorAlpha;
  /* Reference frequency (Hz) for Prior spectral index */
  odouble priorAlphaRefF;
  /* thread number, >0 -> no threading  */
  olong        ithread;
  /* Obit error stack object */
  ObitErr      *err;
  /* UV Interpolator for each MF channel for FTGrid */
  ObitCInterpolate **Interp;
} FTFuncArg;
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSkyModelMF* newObitSkyModelMF (gchar* name)
{
  ObitSkyModelMF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelMFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSkyModelMF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSkyModelMFInit((gpointer)out);

 return out;
} /* end newObitSkyModelMF */

/**
 * Constructor from ObitInfoList.
 * Initializes class if needed on first call.
 * Also works for derived classes.
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string SkyModelMF type, "MF" for base class
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
 *                                type (ObitSkyModelMFCompType as gint), spectral terms;
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
 * \return the new object.
 */
ObitSkyModelMF* ObitSkyModelMFFromInfo (gchar *prefix, ObitInfoList *inList, 
				    ObitErr *err)
{ 
  ObitSkyModelMF *out = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *value=NULL, *classType=NULL;
  ObitImageMosaic *mosaic=NULL;
  olong classCnt, otemp;
  gboolean missing;
  gchar ctemp[50];
  gchar *routine = "ObitSkyModelMFFromInfo";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelMFClassInit();

  /* error checks */
  if (err->error) return out;

  /* check class type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&classType);
  if ((missing) || (type!=OBIT_string)) {
    Obit_log_error(err, OBIT_Error,"%s No class type", routine);
    return out;
  }
  classCnt = dim[0]; /* How many characters in name */
  g_free(keyword);

  /* "xxxmosaic" string prefix of ObitImageMosaic mosaic */
  if (prefix) keyword = g_strconcat (prefix, "mosaic", NULL);
  else        keyword = g_strdup("mosaic");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  /* Does it exist? */
  if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s ImageMosaic not defined in %s", 
		   routine, keyword);
    return out;
  } else { /* exists*/
    mosaic = (ObitImageMosaic*)ObitImageMosaicFromInfo(keyword, inList, err);
    if (err->error) Obit_traceback_val (err, routine, keyword, out);
  }
  g_free(keyword);

  /* Create output - by type */
  if (!strncmp("MF", classType, classCnt)) {
    out = ObitSkyModelMFCreate(prefix, mosaic);
  } else if (!strncmp("Squint", classType, classCnt)) {
    out = (ObitSkyModelMF*)ObitSkyModelVMSquintCreate(prefix, mosaic);
    ObitSkyModelVMSquintFromInfo((ObitSkyModel*)out, prefix, inList, err);
  } else if (!strncmp("Ion", classType, classCnt)) {
    out = (ObitSkyModelMF*)ObitSkyModelVMIonCreate(prefix, mosaic);
    ObitSkyModelVMIonFromInfo((ObitSkyModel*)out, prefix, inList, err);
 } else {  /* Assume MF and hope for the best */
    out = (ObitSkyModelMF*)ObitSkyModelCreate(prefix, mosaic);
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
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doAlphaCorr);
  g_free(keyword);

  /* ""xxxpriorAlpha"    ofloat  prior spectral index applied to be corrected. */
  if (prefix) keyword = g_strconcat (prefix, "priorAlpha", NULL);
  else        keyword = g_strdup("priorAlpha");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->priorAlpha);
  g_free(keyword);

  /* ""xxxpriorAlphaRefF"    odouble prior spectral index ref freq (Hz). */
  if (prefix) keyword = g_strconcat (prefix, "priorAlpha", NULL);
  else        keyword = g_strdup("priorAlphaRefF");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->priorAlphaRefF);
  g_free(keyword);

  /* ""xxxdoSmoo"        olong   Number of threads */
  if (prefix) keyword = g_strconcat (prefix, "doSmoo", NULL);
  else        keyword = g_strdup("doSmoo");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doSmoo);
  g_free(keyword);

 /* Cleanup */
  mosaic = ObitImageMosaicUnref(mosaic);

  return out;
} /* end ObitSkyModelMFFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSkyModelMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelMFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSkyModelMFGetClass */

/**
 * Make a deep copy of an ObitSkyModelMF.
 * \param inn  The object to copy
 * \param outt An existing object pointer for output or NULL if none exists.
 * \param err  Obit error stack object.
 * \return pointer to the new object.
 */
ObitSkyModel* ObitSkyModelMFCopy  (ObitSkyModel *inn, ObitSkyModel *outt, ObitErr *err)
{
  ObitSkyModelMF *in  = (ObitSkyModelMF*)inn;
  ObitSkyModelMF *out = (ObitSkyModelMF*)outt;
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, number;
  gchar *routine = "ObitSkyModelMFCopy";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return (ObitSkyModel*)out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitSkyModelMF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (inn, outt, err);

  /*  copy this class */
  out->info = ObitInfoListUnref(out->info);  /* if it exists */
  if (in->info) out->info = ObitInfoListCopy (in->info);
  out->mosaic = ObitImageMosaicUnref(out->mosaic);  /* if it exists */
  if (in->mosaic) out->mosaic = ObitImageMosaicCopy (in->mosaic, out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, (ObitSkyModel*)out);
  out->plane = ObitFArrayUnref(out->plane);  /* if it exists */
  if (in->plane) out->plane = ObitFArrayCopy (in->plane, out->plane, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, (ObitSkyModel*)out);
  out->FTplane = ObitCArrayUnref(out->FTplane);  /* if it exists */
  if (in->FTplane) out->FTplane = ObitCArrayCopy (in->FTplane, out->FTplane, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, (ObitSkyModel*)out);
  out->myInterp = ObitCInterpolateUnref(out->myInterp);  /* if it exists */
  if (in->myInterp) out->myInterp = ObitCInterpolateCopy (in->myInterp, out->myInterp, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, (ObitSkyModel*)out);
  out->comps = ObitFArrayUnref(out->comps);  /* if it exists */
  if (in->comps) out->comps = ObitFArrayCopy (in->comps, out->comps, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, (ObitSkyModel*)out);
  out->modelType = in->modelType;
  out->modelMode = in->modelMode;
  out->minDFT    = in->minDFT;
  out->maxGrid   = in->maxGrid;
  out->doDFT     = in->doDFT;
  out->doGrid   = in->doGrid;
  if ((out->mosaic) && (out->mosaic->numberImages>0)) {
    number = out->mosaic->numberImages;
    out->CCver = ObitMemRealloc (out->CCver, sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = in->CCver[i];
    out->startComp = ObitMemRealloc (out->startComp, sizeof(olong)*number);
    for (i=0; i<number; i++) out->startComp[i] = in->startComp[i];
    out->endComp = ObitMemRealloc (out->endComp, sizeof(olong)*number);
    for (i=0; i<number; i++) out->endComp[i] = in->endComp[i];
  }

  return (ObitSkyModel*)out;
} /* end ObitSkyModelMFCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an SkyModelMF similar to the input one.
 * \param inn  The object to copy
 * \param outt  An existing object pointer for output, must be defined.
 * \param err  Obit error stack object.
 */
void ObitSkyModelMFClone  (ObitSkyModel *inn, ObitSkyModel *outt, ObitErr *err)
{
  ObitSkyModelMF *in  = (ObitSkyModelMF*)inn;
  ObitSkyModelMF *out = (ObitSkyModelMF*)outt;
  olong number, i;
  const ObitClassInfo *ParentClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  shallow copy this class */
  if (in->info)     out->info     = ObitInfoListRef(in->info);
  if (in->mosaic)   out->mosaic   = ObitImageMosaicRef(in->mosaic);
  if (in->plane)    out->plane    = ObitFArrayRef(in->plane);
  if (in->FTplane)  out->FTplane  = ObitCArrayRef(in->FTplane);
  if (in->myInterp) out->myInterp = ObitCInterpolateRef(in->myInterp);
  if (in->comps)    out->comps    = ObitFArrayRef(in->comps);

  /* Actual copy of some */
  out->modelType = in->modelType;
  out->modelMode = in->modelMode;
  if ((out->mosaic) && (out->mosaic->numberImages>0)) {
    number = out->mosaic->numberImages;
    out->CCver = ObitMemRealloc (out->CCver, sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = in->CCver[i];
    out->startComp = ObitMemRealloc (out->startComp, sizeof(olong)*number);
    for (i=0; i<number; i++) out->startComp[i] = in->startComp[i];
    out->endComp = ObitMemRealloc (out->endComp, sizeof(olong)*number);
    for (i=0; i<number; i++) out->endComp[i] = in->endComp[i];
  }

} /* end ObitSkyModelMFClone */

/**
 * Creates an ObitSkyModelMF 
 * \param name  An optional name for the object.
 * \param mosaic ObitImageMosaic giving one or more images/CC tables
 * \return the new object.
 */
ObitSkyModelMF* ObitSkyModelMFCreate (gchar* name, ObitImageMosaic* mosaic)
{
  ObitSkyModelMF* out;
  olong number, i;

  /* Create basic structure */
  out = newObitSkyModelMF (name);

  /* Modify for input mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);
  if ((out->mosaic) && (out->mosaic->numberImages>0)) {
    number = out->mosaic->numberImages;
    out->CCver = ObitMemRealloc (out->CCver, sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = 0;
    out->startComp = ObitMemRealloc (out->startComp, sizeof(olong)*number);
    out->endComp   = ObitMemRealloc (out->endComp, sizeof(olong)*number);
    for (i=0; i<number; i++) out->startComp[i] = 1;
    for (i=0; i<number; i++) out->endComp[i] = 0;
  }

  return out;
} /* end ObitSkyModelMFCreate */

/**
 * Initializes an ObitSkyModelMF 
 * \param inn  SkyModelMF to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelMFInitMod (ObitSkyModel* inn, ObitUV *uvdata, ObitErr *err)
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  ObitImageMF *image0;
  olong i, j, k, nSpec, nif, nfreq, n;
  ofloat phase=0.5, cp, sp;
  odouble test;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  FTFuncArg *args;
  gchar keyword[12];
  gchar *routine="SkyModelMFInitMod";
  
  /* Fourier transform threading routines */
  in->DFTFunc  = (ObitThreadFunc)ThreadSkyModelMFFTDFT;
  in->GridFunc = (ObitThreadFunc)ThreadSkyModelMFFTGrid;
  
  /* Setup for threading */
  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));
  
  /* Initialize threadArg array */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(FTFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(FTFuncArg)); 
  
    for (i=0; i<in->nThreads; i++) {
      args = (FTFuncArg*)in->threadArgs[i];
      strcpy (args->type, "MF");  /* Enter type as first entry */
      args->in             = in;
      args->uvdata         = uvdata;
      args->ithread        = i;
      args->doAlphaCorr    = in->doAlphaCorr;
      args->priorAlpha     = in->priorAlpha;
      args->priorAlphaRefF = in->priorAlphaRefF;
      args->err            = err;
      if (args->Interp) 
	for (k=0; k<args->nSpec; k++) 
	  args->Interp[k] = ObitCInterpolateUnref(args->Interp[k]);
    }
  } /* end initialize */

  /* Init Sine/Cosine calculator - just to be sure about threading */
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
  ObitInfoListGetTest(image0->myDesc->info, "ALPHA", &type, dim, &in->priorAlpha);
  in->priorAlphaRefF = in->refFreq;
  ObitInfoListGetTest(image0->myDesc->info, "ALPHARF", &type, dim, &in->priorAlphaRefF);
  
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
  if (in->prtLv>1) {
    if (in->currentMode==OBIT_SkyModel_DFT)
      Obit_log_error(err, OBIT_InfoErr, "SkyModelMF using DFT calculation type");
    else if (in->currentMode==OBIT_SkyModel_Grid)
      Obit_log_error(err, OBIT_InfoErr, "SkyModelMF using Grid calculation type");
    else if (in->currentMode==OBIT_SkyModel_Fastest)
      Obit_log_error(err, OBIT_InfoErr, "SkyModelMF using Fastest calculation type");
  }

} /* end ObitSkyModelMFInitMod */

/**
 * Any shutdown operations needed for a model
 * Cleanup structures no longer needed
 * \param inn     SkyModelMF to initialize
 * \param uvdata  uv data being modeled.
 * \param err     Obit error stack object.
 */
void ObitSkyModelMFShutDownMod (ObitSkyModel* inn, ObitUV *uvdata, ObitErr *err)
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  olong i, k;
  FTFuncArg *args;

  in->myInterp = ObitCInterpolateUnref(in->myInterp);
  in->plane    = ObitFArrayUnref(in->plane);
  ObitThreadPoolFree (in->thread);  /* Shut down any threading */
  if (in->threadArgs) {
    /* Check type - only handle "base" */
    if (!strncmp((gchar*)in->threadArgs[0], "MF", 2)) {
      for (i=0; i<in->nThreads; i++) {
	args = (FTFuncArg*)in->threadArgs[i];
	if (args->Interp) {
	  for (k=0; k<args->nSpec; k++) 
	    args->Interp[k] = ObitCInterpolateUnref(args->Interp[k]);
	  g_free(args->Interp);
	}
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
      in->nThreads   = 0;
    } /* end if this a "MF" threadArg */
  }

  /* Cleanup arrays */
  if (in->specFreq)  g_free(in->specFreq);  in->specFreq = NULL;
  if (in->specIndex) g_free(in->specIndex); in->specIndex = NULL;
} /* end ObitSkyModelMFShutDownMod */

/**
 * Initializes an ObitSkyModelMF for a pass through data in time order
 * No op in this class
 * \param inn  SkyModel to initialize
 * \param err  Obit error stack object.
 */
void ObitSkyModelMFInitModel (ObitSkyModel* inn, ObitErr *err)
{
  /* nada */
} /* end ObitSkyModelMFInitModel */

/**
 * Load point model into in comps member.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Adapted from the AIPSish QNOT:VISDFT
 * Output is in member comps with a single row, the entries are
 * \li Amplitude (Jy)
 * \li -2*pi*x (radians)
 * \li -2*pi*y (radians)
 * \li -2*pi*z (radians)
 * \param inn    SkyModelMF 
 * \param uvdata UV data set to model
 * \param err    Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelMFLoadPoint (ObitSkyModel *inn, ObitUV *uvdata, ObitErr *err)
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, cnt, ndim, naxis[2];
  ofloat *table, const2, ccrot, ssrot, cpa, spa, uvrot, xmaj, xmin;
  ofloat dxyzc[3], xxoff, yyoff, zzoff;
  gchar *routine = "ObitSkyModelMFLoadPoint";
  
  /* error checks */
  if (err->error) return gotSome;
  retCode = OBIT_IO_OK;

  /* NEEDS MORE WORK TO BE USEFUL */

  gotSome = (in->pointFlux * in->factor!=0.0);  /* Non zero model? */
  if (!gotSome) return gotSome;

  /* Get position phase shift parameters */
  ObitUVDescShiftPosn(uvdata->myDesc, -in->pointXOff, -in->pointYOff, 
		      dxyzc, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
 
  /* Set field center offsets */
  uvrot  = ObitUVDescRotate(uvdata->myDesc);
  ssrot = sin (DG2RAD * uvrot);
  ccrot = cos (DG2RAD * uvrot);
  xxoff = dxyzc[0] * ccrot + dxyzc[1] * ssrot;
  yyoff = dxyzc[1] * ccrot - dxyzc[0] * ssrot;
  zzoff = dxyzc[2];
   
  in->modType = in->pointParms[3] + 0.5;  /* Model component type */

  /* Include spectrum? */
  if ((in->modType>=10) && (in->modType<=19)) {
    /* Get number */
    cnt = 0;
    for (i=4; i<10; i++) if (in->pointParms[i]!=0.0) cnt++;
    in->nSpecTerm = cnt;
  }

  /* (re)allocate structure */
  ndim = 2;
  naxis[0] = 4; naxis[1] = 1;
  if (in->modType==OBIT_SkyModel_GaussMod)        naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_GaussModTSpec)   naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_USphereMod)      naxis[0] += 2; /* Uniform sphere */
  if (in->modType==OBIT_SkyModel_USphereModTSpec) naxis[0] += 2; /* Uniform sphere */
  /* Any spectral terms */
  naxis[0] += in->nSpecTerm;
  if (in->comps!=NULL) in->comps = ObitFArrayRealloc(in->comps, ndim, naxis);
  else in->comps = ObitFArrayCreate("Components", ndim, naxis);

  /* Fill values */
  naxis[0] = 0; naxis[1]=0;
  table = ObitFArrayIndex(in->comps, naxis);
  table[0] = in->pointFlux * in->factor;
  table[1] = xxoff;
  table[2] = yyoff;
  table[3] = zzoff;

  /* Point + spectrum */
  if (in->modType==OBIT_SkyModel_PointModSpec) {
    for (i=0; i<in->nSpecTerm; i++) table[i+4] = in->pointParms[i+4];
  }

  /* Gaussian */
  if (in->modType==OBIT_SkyModel_GaussMod) {
    /* const2 converts FWHM(asec) to coefficients for u*u, v*v, u*v */
    const2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5) * 2.77777778e-4;
    cpa = cos (DG2RAD * in->pointParms[2]);
    spa = sin (DG2RAD * in->pointParms[2]);
    xmaj = in->pointParms[0] * const2;
    xmin = in->pointParms[1] * const2;
    table[4] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
    table[5] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
    table[6] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
  }

  /* Gaussian + spectrum */
  if (in->modType==OBIT_SkyModel_GaussModSpec) {
    /* const2 converts FWHM(asec) to coefficients for u*u, v*v, u*v */
    const2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5) * 2.77777778e-4;
    cpa = cos (DG2RAD * in->pointParms[2]);
    spa = sin (DG2RAD * in->pointParms[2]);
    xmaj = in->pointParms[0] * const2;
    xmin = in->pointParms[1] * const2;
    table[4] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
    table[5] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
    table[6] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
    for (i=0; i<in->nSpecTerm; i++) table[i+7] = in->pointParms[i+4];
  }

  /* Uniform sphere */
  if (in->modType==OBIT_SkyModel_USphereMod) {
    table[0] = 3.0 * in->pointFlux * in->factor;
    table[4] = in->pointParms[1]  * 0.109662271 * 2.7777778e-4;
    table[5] = 0.1;
 }
  /* Uniform sphere + spectrum*/
  if (in->modType==OBIT_SkyModel_USphereModSpec) {
    table[0] = 3.0 * in->pointFlux * in->factor;
    table[4] = in->pointParms[1]  * 0.109662271 * 2.7777778e-4;
    table[5] = 0.1;
    for (i=0; i<in->nSpecTerm; i++) table[i+6] = in->pointParms[i+4];
 }

  return gotSome;
} /* end ObitSkyModelMFLoadPoint */

/**
 * Load components model into in comps member.
 * If the frequency axis has ctype "SPECLNMF" and if "NSPEC" exists in the 
 * first image descriptor InfoList and is > 0 then there should be a spectrum 
 * in each CLEAN component.
 * Multiplies by factor member and any prior spectral index correction.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Adapted from the AIPSish QNOT:VISDFT
 * \param inn    SkyModelMF 
 * \param n      Image number on mosaic, if -1 load all images
 * \param uvdata UV data set to model
 * \param err    Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelMFLoadComps (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
				  ObitErr *err)
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *imDesc=NULL, *imIODesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  ObitFArray *CompArr=NULL;
  ObitSkyModelCompType modType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong warray, larray;
  ofloat *array, parms[20], *specCorr=NULL;
  olong ver, i, j, hi, lo, count, ncomp, startComp, endComp, irow, lrec;
  olong nspec, iterm, outCCVer, ndim, naxis[2];
  ofloat *table, xxoff, yyoff, zzoff;
  ofloat konst, konst2, xyz[3], xp[3], umat[3][3], pmat[3][3];
  ofloat ccrot, ssrot, xpoff, ypoff, maprot, uvrot;
  ofloat dxyzc[3], cpa, spa, xmaj, xmin, range[2];
  gboolean doCheck=FALSE, want, do3Dmul;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitSkyModelMFLoadComps";
  
  /* error checks */
  if (err->error) return gotSome;

  /* Don't bother if no components requested */
  if ((n>=0) && (in->startComp[n]>in->endComp[n])) return gotSome;

  /* UV descriptor */
  uvDesc = uvdata->myDesc;

  konst = DG2RAD * 2.0 * G_PI;
  /* konst2 converts FWHM(deg) to coefficients for u*u, v*v, u*v */
  konst2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5);

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
    /* DEBUG */
    if (ver<=0) {
      Obit_log_error(err, OBIT_Error,"%s  CC table %d in %s",
		     routine, ver, in->name);
    } /* END DEBUG */

    tempTable = newObitImageTable (in->mosaic->images[i],OBIT_IO_ReadOnly, 
				   tabType, &ver, err);
    if ((tempTable==NULL) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    /* Open */
    retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) {
      Obit_log_error(err, OBIT_Error,"%s problem %d with CC table %d in %s",
		     routine, retCode, ver, in->name);
      Obit_traceback_val (err, routine, in->name, retCode);
    }

    /* How many? */
    endComp = in->endComp[i];
    if (endComp<=0) endComp = CCTable->myDesc->nrow;
    count += MIN(CCTable->myDesc->nrow, endComp) - MAX(1, in->startComp[i]) + 1;

    /* Get model type in first with components */
    /* If only 3 col, or parmsCol 0 size then this is a point model */
    if ((CCTable->myDesc->nfield==3) || 
	(CCTable->parmsCol<0) ||
	(CCTable->myDesc->dim[CCTable->parmsCol]<=0))
      in->modType = OBIT_SkyModel_PointMod;
    if ((in->modType == OBIT_SkyModel_Unknown) && (in->startComp[i]<=endComp)) {
      /* Create table row */
      CCRow = newObitTableCCRow (CCTable);
      /* Read first */
      irow = in->startComp[i];
      retCode = ObitTableCCReadRow (CCTable, irow, CCRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error)) {
	Obit_log_error(err, OBIT_Error,"%s problem %d with CC table %d  %s",
		       routine, retCode, ver, CCTable->name);
	Obit_traceback_val (err, routine, in->name, retCode);
      }

      /* Get model type */
      in->modType = CCRow->parms[3] + 0.5;
      /* Release table row */
      CCRow = ObitTableCCRowUnref (CCRow);
    }

    /* Do we need to check model type */
    doCheck = (CCTable->myDesc->nfield>4) && (CCTable->parmsCol>=0) && 
      (CCTable->myDesc->dim[CCTable->parmsCol][0]>=3);
    
    /* Close */
    retCode = ObitTableCCClose (CCTable, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);

    /* release table  */
    CCTable = ObitTableCCUnref (CCTable);

    /* Is spectral information included? */
    if (!strncmp (in->mosaic->images[i]->myDesc->ctype[in->mosaic->images[i]->myDesc->jlocf], 
		  "SPECLNMF", 8)) {
      /* IO descriptor give true size */
      imIODesc = (ObitImageDesc*)in->mosaic->images[0]->myIO->myDesc;
      nspec = imIODesc->inaxes[imIODesc->jlocf];
      ObitInfoListGetTest (in->mosaic->images[i]->myDesc->info, "NSPEC", &type, dim, &nspec);
      in->nSpec = nspec;  /* Number of terms in the spectrum */
    }

  } /* end loop counting CCs */

  /* (re)allocate structure */
  ndim = 2;
  naxis[0] = 4; naxis[1] = count;
  if (in->modType==OBIT_SkyModel_GaussMod)        naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_GaussModTSpec)   naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_USphereMod)      naxis[0] += 2; /* Uniform sphere */
  if (in->modType==OBIT_SkyModel_USphereModTSpec) naxis[0] += 2; /* Uniform sphere */

  /* Any spectral terms */
  naxis[0] += in->nSpec;

  /* Spectral correction for prior alpha array */
  specCorr = g_malloc0(in->nSpec*sizeof(ofloat));
  if (in->doAlphaCorr && (in->priorAlpha!=0.0)) {
    for (i=0; i<in->nSpec; i++) {
      specCorr[i] = pow((in->specFreq[i]/in->priorAlphaRefF), in->priorAlpha);
    }
  } else { /* No correction */
    for (i=0; i<in->nSpec; i++) specCorr[i] = 1.0;
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
    if ((in->endComp[i]>0) && (in->endComp[i]<in->startComp[i])) continue;

    /* Get CC table */
    outCCVer = 0;
    ver = in->CCver[i];
    startComp = in->startComp[i];
    endComp = in->endComp[i];
    range[0] = in->minDFT;  /* Range of merged fluxes for DFT */
    range[1] = 1.0e20;
    if (endComp>=startComp) {
      CCTable = ObitSkyModelMFgetPBCCTab (in, uvdata, (olong)i, &ver, 
					  &outCCVer, &startComp, &endComp, range, err); 
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    }      

    /* Save values of highest comp - probably bad */
    if (outCCVer==0) {
      /* no translation of table */
      /*??? in->endComp[i] = endComp;*/
    } else {
      /* Translated table with only selected values */
      /*??? in->endComp[i] = in->startComp[i] + endComp-startComp;*/
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
    }
    
    /* Set field center offsets */
    xxoff = dxyzc[0] * ccrot + dxyzc[1] * ssrot;
    yyoff = dxyzc[1] * ccrot - dxyzc[0] * ssrot;
    zzoff = dxyzc[2];

    /* projection rotation matrix if needed */
    do3Dmul = ObitUVDescShift3DMatrix (uvDesc, imDesc, umat, pmat);
    
    /* DEBUG
       fprintf (stderr,"%s: subtracting components %d to %d \n",
       routine, startComp, endComp); */

    /* Convert table to merged array */
    CompArr = ObitTableCCUtilMergeSel (CCTable, startComp, endComp, parms, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    /* entries 0=flux, 1= deltaX 2=deltaY + spectra per merged CC, other parameters in parms */
    naxis[0] = 0; naxis[1]=0; 
    array = ObitFArrayIndex(CompArr, naxis);
    warray = CompArr->naxis[0];
    larray = CompArr->naxis[1];
    modType = (ObitSkyModelCompType)(parms[3]+0.5);  /* model type +20 */

    /* Test for consistency for spectral terms */
    Obit_retval_if_fail ((CompArr->naxis[0]-3>=in->nSpec), err, retCode,
			 "%s Inconsistent request for spectral corrections, %d < %d",  
		       routine,CompArr->naxis[0]-3, in->nSpec);  
   
    /* loop over CCs */
    for (j=0; j<larray; j++) {

      /* Only down to first negative? */
      if (in->noNeg && (array[0]<=0.0)) break;

      /* Do we want this one? Only accept components of in->modType */
      if (doCheck) want = in->modType==modType;
	else want = TRUE;
      want = want && (fabs(array[0])>in->minDFT);
      want = want && (array[0]>in->minFlux);
      want = want && (ncomp<count);  /* don't overflow */
      if (want) {

	/* Point */
	table[0] = array[0] * in->factor;
	xp[0] = (array[1] + xpoff) * konst;
	xp[1] = (array[2] + ypoff) * konst;
	xp[2] = 0.0;
	if (do3Dmul) {
	  xyz[0] = xp[0]*umat[0][0] + xp[1]*umat[1][0];
	  xyz[1] = xp[0]*umat[0][1] + xp[1]*umat[1][1];
	  xyz[2] = xp[0]*umat[0][2] + xp[1]*umat[1][2];
	  /* PRJMUL (2, XP, UMAT, XYZ); */
	} else {  /* no rotn matrix */
	  xyz[0] = ccrot * xp[0] + ssrot * xp[1];
	  xyz[1] = ccrot * xp[1] - ssrot * xp[0];
	  xyz[2] = 0.0;
 	}
	table[1] = xyz[0] + xxoff;
	table[2] = xyz[1] + yyoff;
	table[3] = xyz[2] + zzoff;

	/* DEBUG
        fprintf (stderr,"%s: comp %d   pos %g %g %g xyz %g %g %g poff %g %g\n",
        routine, irow,table[1],table[2],table[3], xxoff, yyoff, zzoff, xpoff, ypoff);  */

	/* Point with tabulated spectrum */
	if (in->modType==OBIT_SkyModel_PointModTSpec) {
	  for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+4] = array[iterm+3]*specCorr[iterm];

	/* Gaussian */
	} else if (in->modType==OBIT_SkyModel_GaussMod) {
	  cpa = cos (DG2RAD * parms[2]);
	  spa = sin (DG2RAD * parms[2]);
	  xmaj = parms[0] * konst2;
	  xmin = parms[1] * konst2;
	  table[4] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
	  table[5] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
	  table[6] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
	  
	/* Gaussian + tabulated spectrum */
	} else if (in->modType==OBIT_SkyModel_GaussModTSpec) {
	  cpa = cos (DG2RAD * parms[2]);
	  spa = sin (DG2RAD * parms[2]);
	  xmaj = parms[0] * konst2;
	  xmin = parms[1] * konst2;
	  table[4] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
	  table[5] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
	  table[6] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
	  /*  spectrum */
	  for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+7] = array[iterm+3]*specCorr[iterm];
	  
	/* Uniform sphere */
	} else if (in->modType==OBIT_SkyModel_USphereMod) {
	  table[0] = 3.0 * array[0] * in->factor;
	  table[4] = parms[1]  * 0.109662271 * 2.7777778e-4;
	  table[5] = 0.1;

	/* Uniform sphere + tabulated spectrum */
	} else if (in->modType==OBIT_SkyModel_USphereModTSpec) {
	  table[0] = 3.0 * array[0] * in->factor;
	  table[4] = parms[1]  * 0.109662271 * 2.7777778e-4;
	  table[5] = 0.1;
	  /*  spectrum */
	  for (iterm=0; iterm<in->nSpec; iterm++) table[iterm+6] = array[iterm+3]*specCorr[iterm];
	}
	    
	/* Update */
	table += lrec;
	array += warray;
	ncomp++;
      } /* End only desired */
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
    
  } /* end loop over tables loading CCs */

  /* Zero any extra entries in table. */
  for (i=ncomp; i<count; i++) {
    /* Zero entry */
    table[0] = 0.0;
    table[1] = 0.0;
    table[2] = 0.0;
    table[3] = 0.0;
    table += lrec;  /* Update pointer */
  } /* end loop zeroing extra components */

  if (specCorr) g_free(specCorr); /* Cleanup */

  /* Find anything */
  gotSome = ncomp>0;

  return gotSome;
} /* end ObitSkyModelMFLoadComps */

/**
 * Grid components model into in plane member and Fourier transform to
 * FTplanes and apply Gaussian taper if needed.
 * Multiplies by factor member and any prior spectral index correction.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * \param inn  SkyModelMF 
 * \param n   Image number on mosaic, 0-rel
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelMFGridComps (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
				  ObitErr *err)
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  gboolean gotSome = FALSE;
  gchar *routine = "ObitSkyModelMFGridComps";
  
  /* error checks */
  if (err->error) return gotSome;
  if ((n<0) || (n>in->mosaic->numberImages-1)) {
    Obit_log_error(err, OBIT_Error,"%s requested field %d out of range [0,%d]",
		   routine, n, in->mosaic->numberImages-1);
      return gotSome;
  }

  /* Load/FT Grid CC table */
  gotSome = ObitSkyModelMFGridFTComps (inn, n, uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  return gotSome;
} /* end ObitSkyModelMFGridComps */

/**
 * NYI
 * Load image model into in plane member and Fourier transform.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * \param inn SkyModelMF 
 * \param n   Image number on mosaic
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (generally true here).
 */
gboolean ObitSkyModelMFLoadImage (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
				  ObitErr *err)
{
  g_error ("ObitSkyModelMFLoadImage not implemented");
  return FALSE;
} /* end ObitSkyModelMFLoadImage */

/**
 * Do Fourier transform using a DFT for a buffer of data.
 * If threading has been enabled by a call to ObitThreadAllowThreads 
 * this routine will divide the buffer up amount the number of processors
 * returned by ObitThreadNumProc.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * If doFlip member is true the Fourier transform is multiplied by sqrt(-1)
 * (for Stokes RL and LR)
 * After the AIPSish QXXPTS, QPTDIV and friends
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * \param inn    SkyModelMF with model components loaded (ObitSkyModelMFLoad)
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to model and subtract from current buffer
 * \param err Obit error stack object.
 */
void ObitSkyModelMFFTDFT (ObitSkyModel *inn, olong field, ObitUV *uvdata, ObitErr *err)
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  olong i, nvis, lovis, hivis, nvisPerThread, nThreads;
  FTFuncArg *args;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelMFFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

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
    args = (FTFuncArg*)in->threadArgs[i];
    strcpy (args->type, "MF");  /* Enter type as first entry */
    args->in     = in;
    args->field  = field;
    args->uvdata = uvdata;
    args->first  = lovis;
    args->last   = hivis;
    if (nThreads>1) args->ithread = i;
    else args->ithread = -1;
    args->err    = err;
    args->Interp = NULL;
    /* Update which vis */
    lovis += nvisPerThread;
    hivis += nvisPerThread;
    hivis = MIN (hivis, nvis);
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nThreads, in->DFTFunc, in->threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
}  /* end ObitSkyModelMFFTDFT */

/**
 * Do Fourier transform using a DFT for a buffer of data.
 * Callable as thread
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * If doFlip member is true the Fourier transform is multiplied by sqrt(-1)
 * (for Stokes RL and LR)
 * After the AIPSish QXXPTS, QPTDIV and friends
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to FTFuncArg argument with elements:
 * \li type   String identifying structure
 * \li in     SkyModelMF with model components loaded (ObitSkyModelMFLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number, <0-> no threads
 * \li err Obit error stack object.
 * \return NULL
 */
static gpointer ThreadSkyModelMFFTDFT (gpointer args)
{
  /* Get arguments from structure */
  FTFuncArg *largs = (FTFuncArg*)args;
  ObitSkyModelMF *in = largs->in;
  /*olong field      = largs->field;*/
  ObitUV *uvdata   = largs->uvdata;
  olong loVis      = largs->first-1;
  olong hiVis      = largs->last;
  ObitErr *err     = largs->err;

  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, ncomp, mcomp, nspec;
  olong lrec, nrparm, naxis[2];
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF;
  olong ilocu, ilocv, ilocw;
  ofloat *visData, *ccData, *data, *fscale;
  ofloat modReal, modImag;
  ofloat amp, arg, freq2, freqFact, wt=0.0, temp;
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArr[FazArrSize], FazArr[FazArrSize], CosArr[FazArrSize], SinArr[FazArrSize];
  olong it, jt, itcnt, ifq, itab;
  odouble tx, ty, tz, sumReal, sumImag, *freqArr;
  gchar *routine = "ThreadSkyModelMFFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Check */
  if (strncmp (largs->type, "MF", 2)) {
    ObitThreadLock(in->thread);  /* Lock against other threads */
    Obit_log_error(err, OBIT_Error,"%s: Wrong type FuncArg %s", routine, largs->type);
    ObitThreadUnlock(in->thread);
    goto finish;
  }

  /* Get pointer for components */
  naxis[0] = 0; naxis[1] = 0; 
  data  = ObitFArrayIndex(in->comps, naxis);
  lcomp = in->comps->naxis[0];  /* Length of row in comp table */
  ncomp = in->comps->naxis[1];  /* number of components */
  if (ncomp<=0) goto finish; /* Anything? */
  nspec = in->nSpecTerm;

  /* Count number of actual components */
  mcomp = 0;
  ccData = data;
  for (iComp=0; iComp<ncomp; iComp++) {
    if (ccData[0]!=0.0) mcomp = iComp+1;
    ccData += lcomp;  /* update pointer */
  } /* end loop over components */
  
  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF  = in->startIFPB-1;
  numberIF = MAX (1, in->numberIFPB);
  jincif   = uvdata->myDesc->incif;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  jincf         = uvdata->myDesc->incf;
  startPoln  = in->startPoln-1;
  numberPoln = in->numberPoln;
  jincs      = uvdata->myDesc->incs;  /* increment in real array */
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocif>=0) {
    if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
      kincf = 1;
      kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
    } else { /* IF before freq  */
      kincif = 1;
      kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
    } 
  } else {  /* NO IF axis */
      kincif = 1;
      kincf  = 1;
  }

  /* DEBUG 
  if (startChannel+numberChannel>7) {
    fprintf (stderr, "DEBUG BAD channel %d st %d no %d stin %d ncin %d\n", 
	     iChannel, startChannel, numberChannel, in->startChannel, in->numberChannel);
  }*/

  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

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

	/* Sum over components */
	sumReal = sumImag = 0.0;
	ccData  = data;
	
	/* Sum by model type */
	switch (in->modType) {
	case OBIT_SkyModel_PointMod:     /* Point */
	  /* From the AIPSish QXXPTS.FOR  */
	  /* outer loop */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    for (iComp=it; iComp<mcomp; iComp++) {
	      if (ccData[0]!=0.0) {  /* valid? */
		tx = ccData[1]*(odouble)visData[ilocu];
		ty = ccData[2]*(odouble)visData[ilocv];
		tz = ccData[3]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		AmpArr[itcnt] = ccData[0];
	      }  /* end if valid */
	      ccData += lcomp;  /* update pointer */
	      itcnt++;          /* Count in amp/phase buffers */
	      if (itcnt>=FazArrSize) break;
	    } /* end inner loop over components */

	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      sumReal += AmpArr[jt]*CosArr[jt];
	      sumImag += AmpArr[jt]*SinArr[jt];
	    }
	  } /*end outer loop over components */
	  break;
	case OBIT_SkyModel_PointModTSpec:     /* Point + tabulated spectrum */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    for (iComp=it; iComp<mcomp; iComp++) {
	      if (ccData[0]!=0.0) {  /* valid? */
		tx = ccData[1]*(odouble)visData[ilocu];
		ty = ccData[2]*(odouble)visData[ilocv];
		tz = ccData[3]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		itab = 4 + in->specIndex[ifq];
		AmpArr[itcnt] = ccData[itab];
	      }  /* end if valid */
	      ccData += lcomp;  /* update pointer */
	      itcnt++;          /* Count in amp/phase buffers */
	      if (itcnt>=FazArrSize) break;
	    } /* end inner loop over components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      sumReal += AmpArr[jt]*CosArr[jt];
	      sumImag += AmpArr[jt]*SinArr[jt];
	    }
	  } /*end outer loop over components */
	  break;
	case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	  /* From the AIPSish QGASUB.FOR  */
	  freq2 = freqFact*freqFact;    /* Frequency factor squared */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    for (iComp=it; iComp<mcomp; iComp++) {
	      if (ccData[0]!=0.0) {  /* valid? */
		arg = freq2 * (ccData[4]*visData[ilocu]*visData[ilocu] +
			       ccData[5]*visData[ilocv]*visData[ilocv] +
			       ccData[6]*visData[ilocu]*visData[ilocv]);
		if (arg<-1.0e-5) amp = ccData[0] * exp (arg);
		else amp = ccData[0];
		tx = ccData[1]*(odouble)visData[ilocu];
		ty = ccData[2]*(odouble)visData[ilocv];
		tz = ccData[3]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		AmpArr[itcnt] = amp;
	      } /* end if valid */
	      ccData += lcomp;  /* update pointer */
	      itcnt++;          /* Count in amp/phase buffers */
	      if (itcnt>=FazArrSize) break;
	    }  /* end inner loop over components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      sumReal += AmpArr[jt]*CosArr[jt];
	      sumImag += AmpArr[jt]*SinArr[jt];
	    }
	  } /*end outer loop over components */
	  break;
	case OBIT_SkyModel_GaussModTSpec:     /* Gaussian on sky + tabulated spectrum*/
	  freq2 = freqFact*freqFact;    /* Frequency factor squared */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    for (iComp=it; iComp<mcomp; iComp++) {
	      if (ccData[0]!=0.0) {  /* valid? */
		itab = 7 + in->specIndex[ifq];
		arg = freq2 * (ccData[4]*visData[ilocu]*visData[ilocu] +
			       ccData[5]*visData[ilocv]*visData[ilocv] +
			       ccData[6]*visData[ilocu]*visData[ilocv]);
		if (arg<-1.0e-5) amp = ccData[itab] * exp (arg);
		else amp = ccData[0];
		tx = ccData[1]*(odouble)visData[ilocu];
		ty = ccData[2]*(odouble)visData[ilocv];
		tz = ccData[3]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		AmpArr[itcnt] = amp;
	      } /* end if valid */
	      ccData += lcomp;  /* update pointer */
	      itcnt++;          /* Count in amp/phase buffers */
	      if (itcnt>=FazArrSize) break;
	    }  /* end inner loop over components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      sumReal += AmpArr[jt]*CosArr[jt];
	      sumImag += AmpArr[jt]*SinArr[jt];
	    }
	  } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	  /* From the AIPSish QSPSUB.FOR  */
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    for (iComp=it; iComp<mcomp; iComp++) {
	      if (ccData[0]!=0.0) {  /* valid? */
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ccData[4];
		arg = MAX (arg, 0.1);
		amp = ccData[0] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		tx = ccData[1]*(odouble)visData[ilocu];
		ty = ccData[2]*(odouble)visData[ilocv];
		tz = ccData[3]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		AmpArr[itcnt] = amp;
	      } /* end if valid */
	      ccData += lcomp;  /* update pointer */
	      itcnt++;          /* Count in amp/phase buffers */
	      if (itcnt>=FazArrSize) break;
	    }  /* end inner ver components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      sumReal += AmpArr[jt]*CosArr[jt];
	      sumImag += AmpArr[jt]*SinArr[jt];
	    }
	  } /* end outer loop over components */
	  break;
	case OBIT_SkyModel_USphereModTSpec:    /* Uniform sphere + tabulated spectrum*/
	  for (it=0; it<mcomp; it+=FazArrSize) {
	    itcnt = 0;
	    for (iComp=it; iComp<mcomp; iComp++) {
	      if (ccData[0]!=0.0) {  /* valid? */
		itab = 6 + in->specIndex[ifq];
		arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				      visData[ilocv]*visData[ilocv]) * ccData[4];
		arg = MAX (arg, 0.1);
		amp = ccData[itab] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		tx = ccData[1]*(odouble)visData[ilocu];
		ty = ccData[2]*(odouble)visData[ilocv];
		tz = ccData[3]*(odouble)visData[ilocw];
		FazArr[itcnt] = freqFact * (tx + ty + tz);
		AmpArr[itcnt] = amp;
	      } /* end if valid */
	      ccData += lcomp;  /* update pointer */
	      itcnt++;          /* Count in amp/phase buffers */
	      if (itcnt>=FazArrSize) break;
	    }  /* end inner loop over components */
	    
	    /* Convert phases to sin/cos */
	    ObitSinCosVec(itcnt, FazArr, SinArr, CosArr);
	    /* Accumulate real and imaginary parts */
	    for (jt=0; jt<itcnt; jt++) {
	      sumReal += AmpArr[jt]*CosArr[jt];
	      sumImag += AmpArr[jt]*SinArr[jt];
	    }
	  } /*end outer loop over components */
	  break;
	default:
	  ObitThreadLock(in->thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s Unknown Comp model type %d in %s",
			 routine, in->modType, in->name);
	  ObitThreadUnlock(in->thread); 
	  goto finish;
	}; /* end switch by model type */

	/* Need to multiply model by sqrt(-1)? */
	if (in->doFlip) {
	  modReal = -(ofloat)sumImag;
	  modImag =  (ofloat)sumReal;
	} else {
	  modReal =  (ofloat)sumReal;
	  modImag =  (ofloat)sumImag;
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

	  /* DEBUG Look for crazy value of data
	  if ((fabs(visData[offset])>10000.0) || (fabs(modReal)>10000.0)) {
	    fprintf (stderr, "DEBUG DATA %g %g chann %d IF %d Stok %d vis %d\n", 
		     visData[offset], modReal, iChannel, iIF, iStoke, iVis);
	  } */

	  /* Ignore blanked data */
	  if ((visData[offset+2]<=0.0) && !in->doReplace) continue;
	  
 	  /* Apply model to data */
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
} /* ThreadSkyModelMFFTDFT */

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
 * \param inn    SkyModelMF with model components loaded (ObitSkyModelMFLoad)
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to model and subtract from current buffer
 * \param err Obit error stack object.
 */
void ObitSkyModelMFFTGrid (ObitSkyModel *inn, olong field, ObitUV *uvdata, ObitErr *err)
{
  ObitSkyModelMF *in  = (ObitSkyModelMF*)inn;
  olong i, k, nvis, lovis, hivis, nvisPerThread, nThreads;
  FTFuncArg *args;
  gboolean OK = TRUE, resetInterp=FALSE;
  gchar *routine = "ObitSkyModelMFFTGrid";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array on first call */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(FTFuncArg*));
    for (i=0; i<in->nThreads; i++) {
      in->threadArgs[i] = g_malloc0(sizeof(FTFuncArg)); 
      args = (FTFuncArg*)in->threadArgs[i];
      strcpy (args->type, "base");  /* Enter type as first entry */
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
    args = (FTFuncArg*)in->threadArgs[i];
    args->in     = in;
    resetInterp = FALSE;
    if (args->Interp) {
      if (args->field!=field) {
	resetInterp = TRUE;
	for (k=0; k<in->nSpec; k++)
	  args->Interp[k] = ObitCInterpolateUnref(args->Interp[k]);
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
    if ((!args->Interp) || resetInterp) {
      args->nSpec = in->nSpec;
      args->Interp = g_malloc0(args->nSpec*sizeof(ObitCInterpolate*));
      resetInterp = TRUE;
    }
    if (resetInterp) {
      for (k=0; k<in->nSpec; k++) {
	if (i>0) {
	  args->Interp[k] = ObitCInterpolateClone(in->myInterps[k], NULL);
	} else {
	  args->Interp[k] = ObitCInterpolateRef(in->myInterps[k]);
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

}  /* end ObitSkyModelMFFTGrid */

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
 * \li in     SkyModelMF with model components loaded (ObitSkyModelMFLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number, <0-> no threads
 * \li err Obit error stack object.
 * \li Interp UV Interpolator
 * \return NULL
 */
gpointer ThreadSkyModelMFFTGrid (gpointer args)
{
  /* Get arguments from structure */
  FTFuncArg *largs = (FTFuncArg*)args;
  ObitSkyModelMF *in = largs->in;
  olong field        = largs->field;
  ObitUV *uvdata     = largs->uvdata;
  olong loVis        = largs->first-1;
  olong hiVis        = largs->last;
  ObitErr *err       = largs->err;
  ObitCInterpolate **Interp = largs->Interp;

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
  gchar *routine = "ThreadSkyModelMFFTGrid";

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
} /* ThreadSkyModelMFFTGrid */

/**
 * Sum the fluxes of components defined by CCVer and endComp
 * \param in  SkyModelMF Checks InfoList member noNeg
 * \param err Obit error stack object.
 * \return Sum of Clean components
 */
ofloat ObitSkyModelMFSum (ObitSkyModel *in, ObitErr *err)
{
  ofloat sum = 0.0;
  olong field;
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitIOCode retCode;
  gchar *tabType = "AIPS CC";
  olong ver, irow;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitSkyModelMFSum";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return sum;
  g_assert (ObitSkyModelMFIsA(in));

  /* Want only positive flux components?? */
  InfoReal.itg = (olong)in->noNeg; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "noNeg", &type, (gint32*)dim, &InfoReal);
  in->noNeg = InfoReal.itg;

  /* Loop over fields */
  for (field=0; field<in-> mosaic->numberImages; field++) {

    /* Check input table to see if there are any selected components */
    /* Get CC table */
    ver = in->CCver[field];
    tempTable = newObitImageTable (in->mosaic->images[field],OBIT_IO_ReadOnly, 
				   tabType, &ver, err);
    Obit_retval_if_fail (((tempTable!=NULL) && (!err->error)), err, sum,
		       "%s: CANNOT Find CC table %d on %s",  
		       routine, ver, in->name);  

    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);
    if (err->error) Obit_traceback_val (err, routine, in->name, sum);
    
    /* Open CC table */
    retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, sum);
    
    /* Create table row */
    CCRow = newObitTableCCRow (CCTable);

    /* If no end CC given - use all */
    if (in->endComp[field]<=0) in->endComp[field] = CCTable->myDesc->nrow;
    /* No more than all */
    in->endComp[field] = MIN (CCTable->myDesc->nrow, in->endComp[field]);
    /* Loop over table summing */
    for (irow=1; irow<=in->endComp[field]; irow++) {
      retCode = ObitTableCCReadRow (CCTable, irow, CCRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error)) 
	Obit_traceback_val (err, routine, in->name, sum);
      if (in->noNeg && (CCRow->Flux<=0.0)) break;
      sum += CCRow->Flux;  /* Sum components */
    }
    
    /* Close Table */
    retCode = ObitTableCCClose (CCTable, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, sum);

    /* Release table row */
    CCRow = ObitTableCCRowUnref (CCRow);
    /* release table  */
    CCTable = ObitTableCCUnref (CCTable);
  } /* end loop over fields */
  
 return sum;
} /* end ObitSkyModelMFSum */

/**
 * Get input parameters from info member
 * \param inn Pointer to the ObitSkyModelMF .
 * \param err Obit error stack object.
 */
void ObitSkyModelMFGetInput (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitSkyModelMFGetInput";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelMFIsA(in));
  if (!ObitInfoListIsA(in->info)) return;
  InfoReal.itg = 0;type = OBIT_oint;

 /* Call base class version */
  ObitSkyModelGetInput (inn, err);
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

} /* end ObitSkyModelMFGetInput */

/**
 * Decide which method is the most appropriate to calculate the FT of a model
 * Sets currentMode member function
 * Sets currentMode member function
 * Adopted from the AIPSish QMTYP
 * \param inn    Pointer to theObitSkyModelMF .
 * \param uvdata UV data set
 */
void  ObitSkyModelMFChose (ObitSkyModel *inn, ObitUV* uvdata) 
{
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  olong nfield, ncc, nx, ny, nvis, nchan, sumcc, i, timff1, timff2, timff3, x;
  olong startComp, endComp;
  ofloat timdft, timfft;

  /* Constants to pick relative times Last determined for FPS 120Bs (Oh yeah!) */
  ofloat tpvgrd = 1.0e-5; /* Time/vis to interpolate (ALGSUB) */
  ofloat tfft=0.8e-6;     /* Time/NX/NY for GRID (CCSGRD) Dependency on grid size.*/
  ofloat tpcgrd=1.0e-4;   /* Time/comp to grid (CCSGRD) dependency on no. comp. */
  ofloat tpvpc=6.6e-7;    /* Time/vis/comp DFT (VISDFT) */

  in->currentMode = in->modelMode;  /* default */
  if (in->currentMode!=OBIT_SkyModel_Fastest) return;

  /* If using point model use DFT */
  if (in->modelType == OBIT_SkyModel_Point) {
    in->currentMode = OBIT_SkyModel_DFT;
    return;
  }

  /* Particulars */
  nfield = in->mosaic->numberImages;
  nvis = uvdata->myDesc->nvis;
  nchan = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  if (uvdata->myDesc->jlocif>=0) 
    nchan *= uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];

  /* Loop over fields */
  sumcc = 0;
  timff1 = timff2 = timff3 = 0;
  for (i=0; i<nfield; i++) {
    startComp = MAX (1, in->startComp[i]);
    endComp = MAX (1, in->endComp[i]);
    ncc = MAX (0, (endComp - startComp + 1));
    if (ncc>0) {
      nx = in->mosaic->nx[i];
      ny = in->mosaic->ny[i];
      x = 4;    /* Small images oversampled */
      if ((nx>2048) || (ny>2048)) x /= 2;
      sumcc = sumcc + ncc;
      timff1 = timff1 + nvis;
      timff2 = timff2 + x * nx * ny;
      timff3 = timff3 + ncc;
     }
  }

  /* How long for gridded method? */
  timfft = (timff1 * tpvgrd + timff2 * tfft*in->nSpec + timff3 * tpcgrd) * nchan;
  /* Ad Hoc hack*/
  timfft *= 2.0;

  /* How long for a DFT? */
  timdft = tpvpc * nvis * sumcc * nchan;

  if (timdft<=timfft) in->currentMode = OBIT_SkyModel_DFT;
  else in->currentMode = OBIT_SkyModel_Grid;

  /* Must do Grid for Image input model - not really supported here */
  if (in->modelType==OBIT_SkyModel_Image) in->currentMode = OBIT_SkyModel_Grid;
} /* end ObitSkyModelMFChose */


/**
 * NYI
 * Sets the in->plane member to either the pixels from the image in the 
 * specified field in in->mosaic or this array with relative 
 * primary beam corrections if in->doPBCor.
 * \param in       SkyModelMF
 * \param uvdata   UV data
 * \param field    Field number in in->mosaic
 * \param err      Obit error stack object.
 * \return ObitCCTable to use, this should be Unref when done and 
 *   Zapped if outCCver != 0
 */
void ObitSkyModelMFgetPBImage (ObitSkyModel* in, ObitUV* uvdata, olong field, 
			       ObitErr *err)
{
  g_error("ObitSkyModelMFgetPBImage NOT implemented");
} /* end ObitSkyModelMFgetPBImage */
  
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
 * \param inn    Pointer to theObitSkyModelMF .
 * \param field  field number (0-rel) in in->mosaic->images
 * \param uvdata UV data set to model
 * \param err    Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelMFGridFTComps (ObitSkyModel* inn, olong field, ObitUV* uvdata, 
				    ObitErr *err)
{
  ObitSkyModelMF *in  = (ObitSkyModelMF*)inn;
  gboolean gotSome = FALSE;
  ObitImageDesc *imDesc = NULL;
  olong i, j, k, nx, ny;
  olong ncomp, ndim, naxis[2];
  ofloat gparm[3], dU, dV, UU, VV, texp;
  ofloat konst, xmaj, xmin, cpa, spa, b1, b2, b3, bb2, bb3;
  ofloat taper, *grid, factor[2];
  gboolean doGaus;
  ObitCArray *FFTImage = NULL;
  gchar *routine = "ObitSkyModelMFGridFTComps";
  /* DEBUG 
  ObitFArray *tempFArray = NULL; */
  /* END DEBUG */

  /* error check */
  if (err->error) return gotSome ;

  /* Create grid, sum components into in->planes */
  ObitSkyModelMFLoadGridComps (inn, field, uvdata, gparm, &ncomp, err);
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
    ObitSkyModelMFFTImage (inn, in->planes[k], FFTImage);

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
      nx = OverSampleMF*imDesc->inaxes[imDesc->jlocr];
      ny = OverSampleMF*imDesc->inaxes[imDesc->jlocd];
      
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
	grid += FFTImage->naxis[0];
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
    factor[0] = OverSampleMF; factor[1] = OverSampleMF;
    in->myInterps[k] = ObitCInterpolateUnref(in->myInterps[k]);
    in->myInterps[k] = 
      newObitCInterpolateCreate("UV data interpolator", in->FTplanes[k], imDesc,
				factor[0], factor[1], in->numConjCol, HWIDTH, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  } /* end loop over planes */

  /* Cleanup */
  FFTImage  = ObitCArrayUnref(FFTImage);


  return gotSome;
} /* end ObitSkyModelMFGridFTComps */

/**
 * Create arrays OverSampleMF times the size of the input image (in->planes) 
 * and sum components onto them.
 * Grid is oversize for increased accuracy.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * \param inn    Pointer to the ObitSkyModelMF .
 * \param field  field number (0-rel) in in->mosaic->images
 * \param uvdata UV data set to model
 * \param gparm  [out] the parameters of the Gaussians in the table
 *               [-1,-1,-1] => not Gaussian.
 * \param ncomp  Actual number of components in in->comps
 * \param err    Obit error stack object.
 */
void  ObitSkyModelMFLoadGridComps (ObitSkyModel* inn, olong field, ObitUV* uvdata, 
				   ofloat gparm[3], olong *ncomp, ObitErr *err)
{
  ObitSkyModelMF *in  = (ObitSkyModelMF*)inn;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableCC *CCTable = NULL;
  ObitImageDesc *imDesc = NULL;
  ofloat range[2], specCorr;
  olong k;
  gchar *tabType = "AIPS CC";
  olong outCCVer, ver, first, last, startComp, endComp;
  gchar *routine = "ObitSkyModelMFLoadGridComps";

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
  CCTable = ObitSkyModelMFgetPBCCTab (in, uvdata, field, &ver, &outCCVer, 
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
    
    retCode = ObitTableCCUtilGridSpect (CCTable, OverSampleMF, k+1,
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
  
} /* end ObitSkyModelMFLoadGridComps */

/**
 * Fourier Transform image array in in->plane, 
 * Half plane complex returned in center-at-the-center order.
 * \param inn      the ObitSkyModelMF .
 * \param inArray  Array to be Transformed.
 * \param outArray Output of FFT, half plane complex
 */
void  ObitSkyModelMFFTImage (ObitSkyModel* inn, ObitFArray *inArray, 
			     ObitCArray *outArray)
{
  /*ObitSkyModelMF *in = (ObitSkyModelMF*)inn;*/
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

} /* end ObitSkyModelMFFTImage  */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param inn     Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string SkyModelMF type, "MF" for this class
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
void ObitSkyModelMFGetInfo (ObitSkyModel *inn, gchar *prefix, ObitInfoList *outList, 
			    ObitErr *err)
{ 
  ObitSkyModelMF *in = (ObitSkyModelMF*)inn;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *OK="OK", *Type="MF";
  olong otemp, numberImages=0;
  gchar *routine = "ObitSkyModelMFGetInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Copy any InfoList keywords */
  if (prefix) keyword = g_strconcat (prefix, "Info", NULL);
  else       keyword = g_strdup("Info");
  ObitInfoListCopyAddPrefix (in->info, outList, keyword);

  /* Class Type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  dim[0] = strlen(Type);
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, Type);
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

} /* end ObitSkyModelMFGetInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSkyModelMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSkyModelMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSkyModelMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSkyModelMFClassInfoDefFn (gpointer inClass)
{
  ObitSkyModelMFClassInfo *theClass = (ObitSkyModelMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSkyModelMFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSkyModelMFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSkyModelMFGetClass;
  theClass->newObit       = (newObitFP)newObitSkyModelMF;
  theClass->ObitSkyModelFromInfo = (ObitSkyModelFromInfoFP)ObitSkyModelMFFromInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitSkyModelMFCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSkyModelMFClear;
  theClass->ObitInit      = (ObitInitFP)ObitSkyModelMFInit;
  theClass->ObitSkyModelInitMod = (ObitSkyModelInitModFP)ObitSkyModelMFInitMod;
  theClass->ObitSkyModelShutDownMod= (ObitSkyModelShutDownModFP)ObitSkyModelMFShutDownMod;
  theClass->ObitSkyModelInitModel= (ObitSkyModelInitModelFP)ObitSkyModelMFInitModel;
  theClass->ObitSkyModelLoadPoint = (ObitSkyModelLoadPointFP)ObitSkyModelMFLoadPoint;
  theClass->ObitSkyModelLoadComps = (ObitSkyModelLoadCompsFP)ObitSkyModelMFLoadComps;
  theClass->ObitSkyModelGridComps = (ObitSkyModelGridCompsFP)ObitSkyModelMFGridComps;
  theClass->ObitSkyModelLoadImage = (ObitSkyModelLoadImageFP)ObitSkyModelMFLoadImage;
  theClass->ObitSkyModelFTDFT     = (ObitSkyModelFTDFTFP)ObitSkyModelMFFTDFT;
  theClass->ObitSkyModelFTGrid    = (ObitSkyModelFTGridFP)ObitSkyModelMFFTGrid;
  theClass->ObitSkyModelGetInput      = (ObitSkyModelGetInputFP)ObitSkyModelMFGetInput;
  theClass->ObitSkyModelChose         = (ObitSkyModelChoseFP)ObitSkyModelMFChose;
  theClass->ObitSkyModelgetPBImage    = (ObitSkyModelgetPBImageFP)ObitSkyModelMFgetPBImage;
  theClass->ObitSkyModelGridFTComps   = (ObitSkyModelGridFTCompsFP)ObitSkyModelMFGridFTComps;
  theClass->ObitSkyModelLoadGridComps = (ObitSkyModelLoadGridCompsFP)ObitSkyModelMFLoadGridComps;
  theClass->ObitSkyModelFTImage       = (ObitSkyModelFTImageFP)ObitSkyModelMFFTImage;
  theClass->ObitSkyModelGetInfo       = (ObitSkyModelGetInfoFP)ObitSkyModelMFGetInfo;

} /* end ObitSkyModelMFClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSkyModelMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModelMF *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->planes      =  NULL;
  in->FTplanes    =  NULL;
  in->myInterps   =  NULL;
  in->specFreq    =  NULL;
  in->specIndex   =  NULL;
  in->refFreq     =  1.0;
  in->nSpec       =  0;
  in->doAlphaCorr =  FALSE;
  in->priorAlphaRefF = 1.0;

} /* end ObitSkyModelMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSkyModelMF* cast to an Obit*.
 */
void ObitSkyModelMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i, k;
  FTFuncArg *args=NULL;
  ObitSkyModelMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->threadArgs) {
    /* Check type - only handle "MF" */
    if (!strncmp((gchar*)in->threadArgs[0], "MF", 2)) {
      for (i=0; i<in->nThreads; i++) {
	args = (FTFuncArg*)in->threadArgs[i];
	if (args->Interp) {
	  for (k=0; k<args->nSpec; k++) 
	    args->Interp[k] = ObitCInterpolateUnref(args->Interp[k]);
	  g_free(args->Interp);
	}
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
    } /* end if this a "MF" threadArg */
  }

  if (in->planes) {
    for (k=0; k<in->nSpec; k++) 
      in->planes[k] = ObitFArrayUnref(in->planes[k]);
    g_free(in->planes);
  }
  if (in->FTplanes) {
    for (k=0; k<in->nSpec; k++) 
      in->FTplanes[k] = ObitCArrayUnref(in->FTplanes[k]);
    g_free(in->FTplanes);
  }
  if (in->myInterps) {
    for (k=0; k<in->nSpec; k++) 
      in->myInterps[k] = ObitCInterpolateUnref(in->myInterps[k]);
    g_free(in->myInterps);
  }
  if (in->specFreq)  g_free(in->specFreq);  in->specFreq = NULL;
  if (in->specIndex) g_free(in->specIndex); in->specIndex = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSkyModelMFClear */

/**
 * Returns the CC table to use for the current set of channels/IF
 * If not making relative PB corrections, this is the input CC table
 * else it is one generated making relative PB corrections.
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
 */
ObitTableCC* ObitSkyModelMFgetPBCCTab (ObitSkyModelMF* in, ObitUV* uvdata, 
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
  ofloat *flux=NULL, *sigma=NULL, *fitResult=NULL, pbmin=0.01, *PBCorr=NULL, alpha;
  ofloat *FreqFact=NULL, *sigmaField=NULL, ll, lll, arg, specFact;
  odouble *Freq=NULL, refFreq;
  odouble Angle=0.0;
  gpointer fitArg=NULL;
  olong irow, row, i, iterm, nterm, offset, nSpec, tiver;
  gchar keyword[12];
  gchar *routine = "ObitSkyModelMFgetPBCCTab";

  /* error checks */
  if (err->error) return CCTable;
  g_assert (ObitSkyModelIsA(in));

  /* Compress/select CC table to new table */
  *outCCVer =  0;  /* Create new one */
  tiver = *inCCVer;
  CCTable = ObitTableCCUtilMergeSel2Tab (in->mosaic->images[field], tiver, outCCVer, 
					    *startCC, *endCC, range, err);
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
    ObitInfoListGetTest(image->myDesc->info, "ALPHA", &type, dim, &alpha);
  
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
   	BeamShape->refFreq = Freq[i];  /* Set frequency */
	PBCorr  = BSClass->ObitBeamShapeGainSym(BeamShape, Angle);
	/* Frequency dependent term */
	lll = ll = FreqFact[i];
	arg = 0.0;
	for (iterm=1; iterm<nterm; iterm++) {
	  arg += fitResult[iterm] * lll;
	  lll *= ll;
	}
	specFact = exp(arg);
	CCRow->parms[offset+i] = specFact*fitResult[0]*PBCorr;
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
} /* end ObitSkyModelMFgetPBCCTab */

