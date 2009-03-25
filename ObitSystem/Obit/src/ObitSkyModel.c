/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-200                                          */
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
#include "ObitSkyModel.h"
#include "ObitSkyModelVMSquint.h"
#include "ObitSkyModelVMIon.h"
#include "ObitImageMosaic.h"
#include "ObitTableCCUtil.h"
#include "ObitFFT.h"
#include "ObitUVUtil.h"
#include "ObitImageUtil.h"
#include "ObitPBUtil.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModel.c
 * ObitSkyModel class function definitions.
 * This class represents sky models and their Fourier transforms and is
 * derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSkyModel";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitSkyModelClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSkyModelClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/** Over sampling factor in uv plane */
olong OverSample=2; 

/*----------------- Macroes ---------------------------*/
/** Half width of gridded subtraction interpolation kernal */
#define HWIDTH 12

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSkyModelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSkyModelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSkyModelClassInfoDefFn (gpointer inClass);

/** Private: Load point model, may be overridden in derived class */
gboolean ObitSkyModelLoadPoint (ObitSkyModel *in, ObitUV *uvdata, ObitErr *err);

/** Private: Load  Components model, may be overridden in derived class */
gboolean ObitSkyModelLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);

/** Private: Threaded FTDFT */
static gpointer ThreadSkyModelFTDFT (gpointer arg);

/** Private: Threaded FTGrid */
static gpointer ThreadSkyModelFTGrid (gpointer arg);

/*---------------Private structures----------------*/
/* FT threaded function argument */
typedef struct {
  /* type "base" in this class */
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
  /* thread number, >0 -> no threading  */
  olong        ithread;
  /* Obit error stack object */
  ObitErr      *err;
  /* UV Interpolator for FTGrid */
  ObitCInterpolate *Interp;
} FTFuncArg;
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSkyModel* newObitSkyModel (gchar* name)
{
  ObitSkyModel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSkyModel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSkyModelInit((gpointer)out);

 return out;
} /* end newObitSkyModel */

/**
 * Constructor from ObitInfoList.
 * Initializes class if needed on first call.
 * Also works for derived classes.
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string SkyModel type, "Base" for base class
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
 * \param err     ObitErr for reporting errors.
 * \return the new object.
 */
ObitSkyModel* ObitSkyModelFromInfo (gchar *prefix, ObitInfoList *inList, 
				    ObitErr *err)
{ 
  ObitSkyModel *out = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *value=NULL, *classType=NULL;
  ObitImageMosaic *mosaic=NULL;
  olong classCnt, otemp;
  gboolean missing;
  gchar ctemp[50];
  gchar *routine = "ObitSkyModelFromInfo";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelClassInit();

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
  if (!strncmp("Base", classType, classCnt)) {
    out = ObitSkyModelCreate(prefix, mosaic);
  } else if (!strncmp("Squint", classType, classCnt)) {
    out = (ObitSkyModel*)ObitSkyModelVMSquintCreate(prefix, mosaic);
    ObitSkyModelVMSquintFromInfo(out, prefix, inList, err);
  } else if (!strncmp("Ion", classType, classCnt)) {
    out = (ObitSkyModel*)ObitSkyModelVMIonCreate(prefix, mosaic);
    ObitSkyModelVMIonFromInfo(out, prefix, inList, err);
 } else {  /* Assume base and hope for the best */
    out = ObitSkyModelCreate(prefix, mosaic);
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

  /* Cleanup */
  mosaic = ObitImageMosaicUnref(mosaic);

  return out;
} /* end ObitSkyModelFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSkyModelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSkyModelGetClass */

/**
 * Make a deep copy of an ObitSkyModel.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSkyModel* ObitSkyModelCopy  (ObitSkyModel *in, ObitSkyModel *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, number;
  gchar *routine = "ObitSkyModelCopy";

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
    out = newObitSkyModel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->info = ObitInfoListUnref(out->info);  /* if it exists */
  if (in->info) out->info = ObitInfoListCopy (in->info);
  out->mosaic = ObitImageMosaicUnref(out->mosaic);  /* if it exists */
  if (in->mosaic) out->mosaic = ObitImageMosaicCopy (in->mosaic, out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  out->plane = ObitFArrayUnref(out->plane);  /* if it exists */
  if (in->plane) out->plane = ObitFArrayCopy (in->plane, out->plane, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  out->FTplane = ObitCArrayUnref(out->FTplane);  /* if it exists */
  if (in->FTplane) out->FTplane = ObitCArrayCopy (in->FTplane, out->FTplane, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  out->myInterp = ObitCInterpolateUnref(out->myInterp);  /* if it exists */
  if (in->myInterp) out->myInterp = ObitCInterpolateCopy (in->myInterp, out->myInterp, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  out->comps = ObitFArrayUnref(out->comps);  /* if it exists */
  if (in->comps) out->comps = ObitFArrayCopy (in->comps, out->comps, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
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

  return out;
} /* end ObitSkyModelCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an SkyModel similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitSkyModelClone  (ObitSkyModel *in, ObitSkyModel *out, ObitErr *err)
{
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

} /* end ObitSkyModelClone */

/**
 * Creates an ObitSkyModel 
 * \param name  An optional name for the object.
 * \param mosaic ObitImageMosaic giving one or more images/CC tables
 * \return the new object.
 */
ObitSkyModel* ObitSkyModelCreate (gchar* name, ObitImageMosaic* mosaic)
{
  ObitSkyModel* out;
  olong number, i;

  /* Create basic structure */
  out = newObitSkyModel (name);

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
} /* end ObitSkyModelCreate */

/**
 * Initializes an ObitSkyModel 
 * \param in  SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelInitMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err)
{
  olong i;
  FTFuncArg *args;

  /* Fourier transform threading routines */
  in->DFTFunc  = (ObitThreadFunc)ThreadSkyModelFTDFT;
  in->GridFunc = (ObitThreadFunc)ThreadSkyModelFTGrid;

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
      strcpy (args->type, "base");  /* Enter type as first entry */
      args->in     = in;
      args->uvdata = uvdata;
      args->ithread= i;
      args->err    = err;
      if (args->Interp) args->Interp = ObitCInterpolateUnref(args->Interp);
    }
  } /* end initialize */

} /* end ObitSkyModelInitMod */

/**
 * Any shutdown operations needed for a model
 * Cleanup structures no longer needed
 * \param in  SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelShutDownMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err)
{
  olong i;
  FTFuncArg *args;

  in->myInterp = ObitCInterpolateUnref(in->myInterp);
  in->plane    = ObitFArrayUnref(in->plane);
  ObitThreadPoolFree (in->thread);  /* Shut down any threading */
  if (in->threadArgs) {
    /* Check type - only handle "base" */
    if (!strncmp((gchar*)in->threadArgs[0], "base", 4)) {
      for (i=0; i<in->nThreads; i++) {
	args = (FTFuncArg*)in->threadArgs[i];
	if (args->Interp) args->Interp = ObitCInterpolateUnref(args->Interp);
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
    } /* end if this a "base" threadArg */
  }
  in->nThreads   = 0;
} /* end ObitSkyModelShutDownMod */

/**
 * Initializes an ObitSkyModel for a pass through data in time order
 * No op in this class
 * \param in  SkyModel to initialize
 * \param err Obit error stack object.
 */
void ObitSkyModelInitModel (ObitSkyModel* in, ObitErr *err)
{
  /* nada */
} /* end ObitSkyModelInitModel */

/**
 * Calculates the Fourier transform of the model and subtracts from UV data
 * \param in      SkyModel to Fourier transform
 * \param indata  UV data set to subtract model from
 * \param outdata UV data set to write to
 * \param err     Obit error stack object.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitSkyModelSubUV (ObitSkyModel *in, ObitUV *indata, ObitUV *outdata, 
			      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitSkyModelClassInfo 
    *myClass=(const ObitSkyModelClassInfo*)in->ClassInfo;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUV *inUV;
  ObitIOAccess access;
  gboolean same, doCalSelect, gotSome, done, isFirst;
  olong i, image, nimage, nload;
  olong firstVis, bufSize;
  ofloat *Buffer;
  gchar *routine = "ObitSkyModelSubUV";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitSkyModelIsA(in));
  g_assert (ObitUVIsA(indata));
  g_assert (ObitUVIsA(outdata));

  /* Get inputs */
  myClass->ObitSkyModelGetInput (in, err);
  if (err->error) goto cleanup;

   /* initialize model */
  myClass->ObitSkyModelInitMod(in, indata, err);
  if (err->error) goto cleanup;

 /* Unless using "point" model, check mosaic */
  if (in->pointFlux==0.0) {
    if (!ObitImageMosaicIsA(in->mosaic)) {
      Obit_log_error(err, OBIT_Error,"%s mosaic member not defined in %s",
		     routine, in->name);
      goto cleanup;
    }
    
    /* valid images? */
    for (i=0; i<in->mosaic->numberImages; i++) {
      if (!ObitImageIsA(in->mosaic->images[i])) {
	Obit_log_error(err, OBIT_Error,"%s mosaic image %d not defined in %s",
		       routine, i+1, in->name);
	goto cleanup;
      }
    }
  } /* end check images */

  /* DEBUG
  err->error = 1;*/
  if (err->error) goto cleanup;

  /* Do we need to calibrate/select input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(indata->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;  /* Output may be same as input */

  /* If input and output are the same, doCalSelect MUST be FALSE*/
  /* Are input and output the same file? */
  same = ObitUVSame(indata, outdata, err);
  if (err->error) goto cleanup;
  Obit_retval_if_fail ((!(same && doCalSelect)), err, retCode,
		       "%s CANNOT rewrite files with doCalSelect=TRUE",  
		       routine);  

  /* Open input uv data */
  retCode = ObitUVOpen (indata, access, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Update frequency tables on indata */
  ObitUVGetFreq (indata, err);
  if (err->error) goto cleanup;

  /* Fill in data selection */
  myClass->ObitSkyModelSetSelect (in, indata, err);
  if (err->error) goto cleanup;
 
  /* Choose mode if requested */
  myClass->ObitSkyModelChose (in, indata);

  /* Prepare output to be compatable with input */
  if (!ObitUVSame(indata, outdata, err)) {
    ObitUVClone (indata, outdata, err);
  } else { /* Copy descriptor to be sure */
    outdata->myDesc = (gpointer)ObitUVDescCopy(indata->myDesc, outdata->myDesc, err);
  }
  if (err->error) goto cleanup;

  /* use same data buffer on input and output.
     If multiple passes are made the input files will be closed
     which deallocates the buffer, use output buffer.
     so free input buffer */
  if (!same) {
    if (indata->buffer) ObitIOFreeBuffer(indata->buffer); /* free existing */

    /* Open output */
    retCode = ObitUVOpen (outdata, OBIT_IO_ReadWrite, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
    bufSize = outdata->bufferSize;  /* Save buffer size */
    Buffer  = outdata->buffer; 
    indata->bufferSize =  bufSize; 
    indata->buffer = Buffer; 
  }


  inUV = indata;  /* First pass read input */
  /* Loop over images */
  if (in->mosaic!=NULL) nimage = in->mosaic->numberImages;  /* How many in mosaic? */
  else nimage = 1;
  isFirst = !same;  /* keep track of first of several passes */
  if (in->currentMode==OBIT_SkyModel_DFT) nimage = 1;  /* All at once with DFT */
  for (image = 0; image<nimage; image++) {

    /* Loop over blocks of channels for PB corrections */
    in->startIFPB = in->startChannelPB = -1;  /* to initialize setPBChans */
    done = myClass->ObitSkyModelsetPBChans(in, indata, err);
    if (err->error) goto cleanup;

    while (!done) {

      done = TRUE;  /* until proven otherwise */
      /* Load image model */
      nload = image;   /* Which image/images to load */
      if (in->currentMode==OBIT_SkyModel_DFT) nload = -1;
      gotSome = myClass->ObitSkyModelLoad (in, nload, indata, err);
      if (err->error) goto cleanup;
      if (!gotSome) continue;   /* Anything to do? */

      /* Any model initialization at beginning of pass */
      myClass->ObitSkyModelInitModel(in, err);

      /* Are input and output the same file? */
      same = ObitUVSame(inUV, outdata, err);
      if (err->error) goto cleanup;
      
      /* Loop over data file */
      retCode = OBIT_IO_OK;
      while (retCode==OBIT_IO_OK) {
	if (doCalSelect) retCode = ObitUVReadSelect (inUV, inUV->buffer, err);
	else retCode = ObitUVRead (inUV, inUV->buffer, err);
	if (retCode!=OBIT_IO_OK) break;
	firstVis = inUV->myDesc->firstVis;
	
	/* Calculate and subtract model */
	myClass->ObitSkyModelFT (in, image, inUV, err);
	
	/* How many */
	outdata->myDesc->numVisBuff = inUV->myDesc->numVisBuff;
	
	/* rewrite */
	retCode = ObitUVWrite (outdata, inUV->buffer, err);
	/* suppress vis number update if rewriting the same file */
	if (same) {
	  outdata->myDesc->firstVis = firstVis;
	  ((ObitUVDesc*)(inUV->myIO->myDesc))->firstVis = firstVis;
	}
	/* check for errors */
	if ((retCode > OBIT_IO_EOF) || (err->error)) {
	  /* First unset input buffer (may be multiply deallocated ;'{ ) */
	  inUV->buffer     = NULL;
	  inUV->bufferSize = 0;
	  goto cleanup;
	}
      } /* end loop over uv data file */
      
      /* Further passes reread output */
      if (isFirst) {
	Buffer  = outdata->buffer;      /* Save buffer pointer */
	bufSize = outdata->bufferSize;  /* Save buffer size */
	/* Keep from deallocating buffer */
	inUV->buffer     = NULL;
	inUV->bufferSize = 0;
	retCode = ObitUVClose (inUV, err);
	if (err->error) goto cleanup;

	/* To be sure, restore output buffer info */
	outdata->buffer     = Buffer;
	outdata->bufferSize = bufSize;
	inUV = outdata;
	doCalSelect = FALSE;
	/* Update frequency tables on outdata */
	ObitUVGetFreq (outdata, err);
	if (err->error) goto cleanup;
	isFirst = FALSE;
      }
      
      retCode = ObitUVIOSet (inUV, err); /* reset to beginning of uv data */
      if (err->error) goto cleanup;
      
      done = myClass->ObitSkyModelsetPBChans(in, indata, err);  /* Update channels selected */
      if (err->error) goto cleanup;
      if (done) break;

    } /* end loop over PB channel blocks */

    /* if doReplace, have already thrown data away - just accumulate model */
    if (in->doReplace) {
      in->doReplace = FALSE; 
      in->factor = -fabs(in->factor);
    }
  } /* end loop over images */
  
  /* close files */
  Buffer  = outdata->buffer;      /* Save buffer pointer */
  bufSize = outdata->bufferSize;  /* Save buffer size */

  /* keep from multiple deallocations of buffer */
  indata->buffer     = NULL;
  indata->bufferSize = 0;
  retCode = ObitUVClose (indata, err);
  if (err->error) goto cleanup;

  /* To be sure, restore output buffer info */
  outdata->buffer     = Buffer;
  outdata->bufferSize = bufSize;
  retCode = ObitUVClose (outdata, err);
  if (err->error) goto cleanup;

  /* Cleanup */
 cleanup:
  myClass->ObitSkyModelShutDownMod(in, indata, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
  
  return retCode;
}  /* end ObitSkyModelSubUV */

/**
 * Calculates the Fourier transform of the model and divides into UV data
 * \param in      SkyModel to Fourier transform
 * \param indata  UV data set to subtract model from
 * \param outdata UV data set to write to
 * \param err Obit error stack object.
 * \return return code, OBIT_IO_OK=> OK
 */
ObitIOCode ObitSkyModelDivUV (ObitSkyModel *in, ObitUV *indata, ObitUV *outdata,
			      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  const ObitSkyModelClassInfo 
    *myClass=(const ObitSkyModelClassInfo*)in->ClassInfo;
  ObitUV *scratchUV=NULL, *workUV, *inUV;
  gboolean saveDoDivide, gotSome;
  ofloat saveFactor, *Buffer;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOAccess access;
  gboolean same, doCalSelect, doScratch, done, isFirst;
  olong i, image, nimage, nload;
  olong firstVis, bufSize;
  gchar *routine = "ObitSkyModelDivUV";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitSkyModelIsA(in));
  g_assert (ObitUVIsA(indata));
	    
  /* Get inputs */
  myClass->ObitSkyModelGetInput (in, err);
  if (err->error) goto cleanup;

  /* initialize model */
  myClass->ObitSkyModelInitMod(in, indata, err);
  if (err->error) goto cleanup;

  /* Unless using "point" model, check mosaic */
  if (in->pointFlux==0.0) {
    if (!ObitImageMosaicIsA(in->mosaic)) {
      Obit_log_error(err, OBIT_Error,"%s mosaic member not defined in %s",
		     routine, in->name);
      goto cleanup;
    }
    
    /* valid images? */
    for (i=0; i<in->mosaic->numberImages; i++) {
      if (!ObitImageIsA(in->mosaic->images[i])) {
	Obit_log_error(err, OBIT_Error,"%s mosaic image %d not defined in %s",
		       routine, i+1, in->name);
	goto cleanup;
      }
    }
  } /* end check images */

  /* Save previous Divide start */
  saveDoDivide = in->doDivide;
  saveFactor   = in->factor;

  /* Do we need to calibrate/select input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(indata->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;  /* Output may be same as input */

  /* If input and output are the same, doCalSelect MUST be FALSE*/
  /* Are input and output the same file? */
  same = ObitUVSame(indata, outdata, err);
  if (err->error) goto cleanup;
  Obit_retval_if_fail ((!(same && doCalSelect)), err, retCode,
		       "%s CANNOT rewrite files with doCalSelect=TRUE",  
		       routine);  

  /* Open input uv data */
  retCode = ObitUVOpen (indata, access, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Update frequency tables on indata */
  ObitUVGetFreq (indata, err);
  if (err->error) goto cleanup;

  /* Fill in data selection */
  myClass->ObitSkyModelSetSelect (in, indata, err);
  if (err->error) goto cleanup;

  /* Choose mode if requested */
  myClass->ObitSkyModelChose (in, indata);

  /* Prepare output to be compatable with input */
  if (!ObitUVSame(indata, outdata, err)) {
    ObitUVClone (indata, outdata, err);
  } else { /* Copy descriptor to be sure */
    outdata->myDesc = (gpointer)ObitUVDescCopy(indata->myDesc, outdata->myDesc, err);
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  /* If there are more than one fields and not doing DFT then create 
     zeroed scratch file version of input, accumulate the model and then divide at
     the end */
  doScratch = (in->mosaic!=NULL) && (in->mosaic->numberImages>1) && 
    (in->currentMode!=OBIT_SkyModel_DFT);
  /* For Mixed model always have to accumulate then divide */
  doScratch = doScratch || in->currentMode==OBIT_SkyModel_Mixed;
  if (doScratch) {
    
    retCode = ObitUVClose (indata, err);  /* Close indata */
    if (err->error) Obit_traceback_val (err, routine,indata->name, retCode);
    /* Get zeroed scatch version of indata */
    scratchUV = ObitUVUtilCopyZero (indata, TRUE, NULL, err);
    if (err->error) goto cleanup;
    workUV = scratchUV; /* Accumulate model to output */
    inUV   = scratchUV; /* All passes read scratchfile */
    in->doDivide = FALSE;               /* Accumulate model by summing */
    in->factor = -1.0*fabs(in->factor); /* want sum */
    
    /* Open new input uv data */
    retCode = ObitUVOpen (scratchUV, OBIT_IO_ReadWrite, err);
  } else {  /* One pass division */
    inUV   = indata;
    workUV = outdata;  
    in->doDivide = TRUE;
  } /* end setup for scratch file */
  
  /* Are input and output the same file? */
  same = ObitUVSame(inUV, workUV, err);
  if (err->error) goto cleanup;

  /* use same data buffer on input and output 
     If multiple passes are made the input files will be closed
     which deallocates the buffer, use output buffer.
     so free input buffer */
  if (!same) {
    if (inUV->buffer) ObitIOFreeBuffer(inUV->buffer); /* free existing */
    /* Open output */
    retCode = ObitUVOpen (workUV, OBIT_IO_ReadWrite, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  }

  bufSize = workUV->bufferSize;  /* Save buffer size */
  Buffer  = workUV->buffer; 
  inUV->bufferSize =  bufSize; 
  inUV->buffer = Buffer; 

  /* Loop over images */
  if (in->mosaic!=NULL) nimage = in->mosaic->numberImages;  /* How many in mosaic? */
  else nimage = 1;
  isFirst = !same;  /* keep track of first of several passes */
  if (in->currentMode==OBIT_SkyModel_DFT) nimage = 1;  /* All at once with DFT */
  for (image = 0; image<nimage; image++) {

    /* Loop over blocks of channels for PB corrections */
    in->startIFPB = in->startChannelPB = -1;  /* to initialize setPBChans */
    done = myClass->ObitSkyModelsetPBChans(in, indata, err);
    if (err->error) goto cleanup;

    while (!done) {

      done = TRUE;  /* until proven otherwise */
      /* Load image model */
      nload = image;   /* Which image/images to load */
      if (in->currentMode==OBIT_SkyModel_DFT) nload = -1;
      gotSome = myClass->ObitSkyModelLoad (in, nload, indata, err);
      if (err->error) goto cleanup;
      if (!gotSome) continue;   /* Anything to do? */

      /* Any model initialization at beginning of pass */
      myClass->ObitSkyModelInitModel(in, err);

      /* Are input and output the same file? */
      same = ObitUVSame(inUV, workUV, err);
      if (err->error) goto cleanup;
      
      /* Loop over data file */
      retCode = OBIT_IO_OK;
      while (retCode==OBIT_IO_OK) {
	if (doCalSelect) retCode = ObitUVReadSelect (inUV, inUV->buffer, err);
	else retCode = ObitUVRead (inUV, inUV->buffer, err);
	if (err->error) goto cleanup;
	if (retCode!=OBIT_IO_OK) break;
	firstVis = inUV->myDesc->firstVis;
	
	/* Calculate and subtract model */
	myClass->ObitSkyModelFT (in, image, inUV, err);
	if (err->error) goto cleanup;
	
	/* How many */
	workUV->myDesc->numVisBuff = inUV->myDesc->numVisBuff;
	
	/* rewrite */
	retCode = ObitUVWrite (workUV, Buffer, err);
	if (err->error) goto cleanup;
	/* suppress vis number update if rewriting the same file */
	if (same) {
	  workUV->myDesc->firstVis = firstVis;
	  ((ObitUVDesc*)(inUV->myIO->myDesc))->firstVis = firstVis;
	}
	/* check for errors */
	if ((retCode > OBIT_IO_EOF) || (err->error)) {
	  /* First unset input buffer (may be multiply deallocated ;'{ ) */
	  inUV->buffer     = NULL;
	  inUV->bufferSize = 0;
	  goto cleanup;
	}
      } /* end loop over uv data file */
      
      /* Further passes reread output */
      if (isFirst) {
	Buffer  = workUV->buffer;      /* Save buffer pointer */
	bufSize = workUV->bufferSize;  /* Save buffer size */
	/* Keep from deallocating buffer */
	inUV->buffer     = NULL;
	inUV->bufferSize = 0;
	retCode = ObitUVClose (inUV, err);
	if (err->error) goto cleanup;

	/* To be sure, restore output buffer info */
	doCalSelect = FALSE;
	workUV->buffer     = Buffer;
	workUV->bufferSize = bufSize;
	inUV = workUV;
	/* Update frequency tables on workUV */
	ObitUVGetFreq (workUV, err);
	if (err->error) goto cleanup;
	isFirst = FALSE;
      }
      
      retCode = ObitUVIOSet (inUV, err); /* reset to beginning of uv data */
      if (err->error) goto cleanup;
      
      done = myClass->ObitSkyModelsetPBChans(in, indata, err);  /* Update channels selected */
      if (err->error) goto cleanup;
      if (done) break;
      
    } /* end loop over PB channel blocks */ 
  } /* end loop over images */
  
  /* close files */
  Buffer  = workUV->buffer;      /* Save buffer pointer */
  bufSize = workUV->bufferSize;  /* Save buffer size */

  /* keep from multiple deallocations of buffer */
  inUV->buffer     = NULL;
  inUV->bufferSize = 0;
  retCode = ObitUVClose (inUV, err);
  if (err->error) goto cleanup;

  /* First unset output buffer (may be multiply deallocated ;'{ ) */
  workUV->buffer     =  Buffer;
  workUV->bufferSize = bufSize;
  retCode = ObitUVClose (workUV, err);
  if (err->error) goto cleanup;
  
  /* If model was accumulated in scratch file, do division */
  if (doScratch) {
    ObitUVUtilVisDivide (indata, scratchUV, outdata, err);
    if (err->error) goto cleanup;
 }

  /* Restore Divide request */
  in->doDivide = saveDoDivide;
  in->factor   = saveFactor;

  /* Cleanup */
 cleanup:
  scratchUV    = ObitUVUnref(scratchUV);
  myClass->ObitSkyModelShutDownMod(in, indata, err);
  if (err->error) Obit_traceback_val (err, routine,in->name, retCode);
  
  return OBIT_IO_OK;
}  /* end ObitSkyModelDivUV */

/**
 * Loads the model for a specified image on mosaic member
 * \param in      SkyModel to Fourier transform
 * \param image   Image number in mosaic member to load model, 0-rel
 *                If  <0 , load all
 * \param uvdata UV data set to model
 * \param err     Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelLoad (ObitSkyModel *in, olong image, ObitUV *uvdata, 
			   ObitErr *err)
{
  gboolean gotSome = FALSE;
  const ObitSkyModelClassInfo *myClass;
  olong n;
  gchar *routine = "ObitSkyModelLoad";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return gotSome;
  g_assert (ObitSkyModelIsA(in));

  myClass = in->ClassInfo; /* Need class structure */

 /* Do by model type */
  switch (in->modelType) {
  case OBIT_SkyModel_Point:
    gotSome = myClass->ObitSkyModelLoadPoint(in, uvdata, err);
    break;
  case OBIT_SkyModel_Comps:
    n = image;  /* image or images */
    /* Load components or grid/FT */
    if (in->currentMode==OBIT_SkyModel_Grid)
       gotSome = myClass->ObitSkyModelGridComps(in, n, uvdata, err);
    else if (in->currentMode==OBIT_SkyModel_DFT)
       gotSome = myClass->ObitSkyModelLoadComps(in, n, uvdata, err);
    else if (in->currentMode==OBIT_SkyModel_Mixed) {
       in->doDFT  = myClass->ObitSkyModelLoadComps(in, n, uvdata, err);
       in->doGrid = myClass->ObitSkyModelGridComps(in, n, uvdata, err);
       gotSome = in->doDFT || in->doGrid;
    }

    break;
  case OBIT_SkyModel_Image:
     gotSome = myClass->ObitSkyModelLoadImage(in, image, uvdata, err);
    break;
  default:
    Obit_log_error(err, OBIT_Error,"%s Unknown model type %d in %s",
		   routine, in->modelType, in->name);
    return  FALSE;
  }; /* end switch by model type */
  if (err->error) Obit_traceback_val (err, routine, in->name, FALSE);

  return gotSome;
}  /* end ObitSkyModelLoad */

/**
 * Calculates the Fourier transform for current UV buffer.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * If doFlip member is true the Fourier transform is multiplied by sqrt(-1)
 * (for Stokes RL and LR)
 * Processes all data selected.
 * \param in     SkyModel to Fourier transform
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to subtract model from
 * \param err    Obit error stack object.
 */
void ObitSkyModelFT (ObitSkyModel *in, olong field, ObitUV *uvdata, ObitErr *err)
{
  const ObitSkyModelClassInfo *myClass;
  gchar *routine = "ObitSkyModelFT";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelIsA(in));
  g_assert (ObitUVIsA(uvdata));
	    
  myClass = in->ClassInfo; /* Need class structure */
  
 /* Do by model type */
  switch (in->currentMode) {
  case OBIT_SkyModel_DFT:
    myClass->ObitSkyModelFTDFT(in, field, uvdata, err);
    break;
  case OBIT_SkyModel_Grid:
    myClass->ObitSkyModelFTGrid(in, field, uvdata, err);
    break;
  case OBIT_SkyModel_Mixed:
    if (in->doDFT)  myClass->ObitSkyModelFTDFT(in, field, uvdata, err);
    if (in->doGrid) myClass->ObitSkyModelFTGrid(in, field, uvdata, err);
    break;
  default:
    Obit_log_error(err, OBIT_Error,"%s Unknown mode type %d in %s",
		   routine, in->currentMode, in->name);
    return;
  }; /* end switch by model type */
	    
  if (err->error) Obit_traceback_msg (err, routine, in->name);

}  /* end ObitSkyModelFT */

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
 * \param in  SkyModel 
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelLoadPoint (ObitSkyModel *in, ObitUV *uvdata, ObitErr *err)
{
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, cnt, ndim, naxis[2];
  ofloat *table, const2, ccrot, ssrot, cpa, spa, uvrot, xmaj, xmin;
  ofloat dxyzc[3], xxoff, yyoff, zzoff;
  gchar *routine = "ObitSkyModelLoadPoint";
  
  /* error checks */
  if (err->error) return gotSome;
  retCode = OBIT_IO_OK;

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
  if (in->modType>=10) {
    /* Get number */
    cnt = 0;
    for (i=4; i<10; i++) if (in->pointParms[3+i]!=0.0) cnt++;
    in->nSpecTerm = cnt;
  }

 /* (re)allocate structure */
  ndim = 2;
  naxis[0] = 4; naxis[1] = 1;
  if (in->pointParms[3]==1.) naxis[0] += 3; /* Gaussian */
  if (in->pointParms[3]==2.) naxis[0] += 2; /* Uniform sphere */
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
  for (i=0; i<in->nSpecTerm; i++) table[i+4] = in->pointParms[i+4];

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

  /* Uniform sphere */
  if (in->modType==OBIT_SkyModel_USphereMod) {
    table[0] = 3.0 * in->pointFlux * in->factor;
    table[4] = in->pointParms[1]  * 0.109662271 * 2.7777778e-4;
    table[5] = 0.1;
 }

  return gotSome;
} /* end ObitSkyModelLoadPoint */

/**
 * Load components model into in comps member.
 * If the frequency axis has ctype "SPECLOGF" and if "NTERM" exists in the first image 
 * descriptor InfoList and is > 0 and there are at least nterm planes in the image, 
 * then there should be spectral information in the CC tables which is used.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Adapted from the AIPSish QNOT:VISDFT
 * \param in  SkyModel 
 * \param n   Image number on mosaic, if -1 load all images
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				ObitErr *err)
{
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *imDesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  ObitFArray *CompArr=NULL;
  ObitSkyModelCompType modType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong warray, larray;
  ofloat *array, parms[20];
  olong ver, i, j, hi, lo, count, ncomp, startComp, endComp, irow, lrec;
  olong nterm, iterm, outCCVer, ndim, naxis[2];
  ofloat *table, xxoff, yyoff, zzoff;
  ofloat konst, konst2, xyz[3], xp[3], umat[3][3], pmat[3][3];
  ofloat ccrot, ssrot, xpoff, ypoff, maprot, uvrot;
  ofloat dxyzc[3], cpa, spa, xmaj, xmin, range[2];
  gboolean doCheck=FALSE, want, do3Dmul;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitSkyModelLoadComps";
  
  /* error checks */
  if (err->error) return gotSome;

  /* Don't bother if no components requested */
  if ((n>=0) && (in->startComp[n]>in->endComp[n])) return gotSome;

  /* Uv descriptor */
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
		  "SPECLOGF", 8)) {
      nterm = in->mosaic->images[i]->myDesc->inaxes[in->mosaic->images[i]->myDesc->jlocf];
      ObitInfoListGetTest (in->mosaic->images[i]->myDesc->info, "NTERM", &type, dim, &nterm);
      in->nSpecTerm = nterm -1;  /* Only higher order terms */
    }

  } /* end loop counting CCs */

  /* (re)allocate structure */
  ndim = 2;
  naxis[0] = 4; naxis[1] = count;
  if (in->modType==OBIT_SkyModel_GaussMod)   naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_USphereMod) naxis[0] += 2; /* Uniform sphere */

  /* Any spectral terms */
  naxis[0] += in->nSpecTerm;

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
    CCTable = ObitSkyModelgetPBCCTab (in, uvdata, (olong)i, &ver, &outCCVer, 
				      &startComp, &endComp, range, err); 
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
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
    xpoff = (imDesc->crpix[imDesc->jlocr] - 
	     (imDesc->inaxes[imDesc->jlocr]/2)) * 
      imDesc->cdelt[imDesc->jlocr];
    ypoff = (imDesc->crpix[imDesc->jlocd] - 
	     (imDesc->inaxes[imDesc->jlocd]/2) - 1) * 
      imDesc->cdelt[imDesc->jlocd];
    /* These should always be zero for 3D imaging? */
    xpoff = 0.0;
    ypoff = 0.0;
    
    /* Set field center offsets */
    xxoff = dxyzc[0] * ccrot + dxyzc[1] * ssrot;
    yyoff = dxyzc[1] * ccrot - dxyzc[0] * ssrot;
    zzoff = dxyzc[2];

    /* 3D rotation matrix if needed */
    if (in->do3D) {
      do3Dmul = ObitUVDescShift3DMatrix (uvDesc, imDesc, umat, pmat);
    } else {do3Dmul = FALSE;}
    
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
    modType = (ObitSkyModelCompType)(parms[3]+0.5);  /* model type */

    /* Test for consistency for spectral terms */
    Obit_retval_if_fail ((CompArr->naxis[0]-3>=in->nSpecTerm), err, retCode,
			 "%s Inconsistent request for spectral corrections, %d < %d",  
		       routine,CompArr->naxis[0]-3, in->nSpecTerm);  
   
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
	} else {  /* no 3D */
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

	/* Point with spectrum */
	if (in->modType==OBIT_SkyModel_PointModSpec) {
	  for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+4] = array[iterm+3];

	/* Gaussian */
	} else if (in->modType==OBIT_SkyModel_GaussMod) {
	  cpa = cos (DG2RAD * parms[2]);
	  spa = sin (DG2RAD * parms[2]);
	  xmaj = parms[0] * konst2;
	  xmin = parms[1] * konst2;
	  table[4] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
	  table[5] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
	  table[6] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
	  
	/* Gaussian + spectrum */
	} else if (in->modType==OBIT_SkyModel_GaussModSpec) {
	  cpa = cos (DG2RAD * parms[2]);
	  spa = sin (DG2RAD * parms[2]);
	  xmaj = parms[0] * konst2;
	  xmin = parms[1] * konst2;
	  table[4] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
	  table[5] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
	  table[6] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
	  /*  spectrum */
	  for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+7] = array[iterm+3];
	  
	/* Uniform sphere */
	} else if (in->modType==OBIT_SkyModel_USphereMod) {
	  table[0] = 3.0 * array[0] * in->factor;
	  table[4] = parms[1]  * 0.109662271 * 2.7777778e-4;
	  table[5] = 0.1;

	/* Uniform sphere+ spectrum */
	} else if (in->modType==OBIT_SkyModel_USphereModSpec) {
	  table[0] = 3.0 * array[0] * in->factor;
	  table[4] = parms[1]  * 0.109662271 * 2.7777778e-4;
	  table[5] = 0.1;
	  /*  spectrum */
	  for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+6] = array[iterm+3];
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

  /* Find anything */
  gotSome = ncomp>0;

  return gotSome;
} /* end ObitSkyModelLoadComps */

/**
 * Grid components model into in plane member and Fourier transform to
 * FTplane and apply Gaussian taper if needed.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * \param in  SkyModel 
 * \param n   Image number on mosaic, 0-rel
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelGridComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err)
{
  gboolean gotSome = FALSE;
  gchar *routine = "ObitSkyModelGridComps";
  
  /* error checks */
  if (err->error) return gotSome;
  if ((n<0) || (n>in->mosaic->numberImages-1)) {
    Obit_log_error(err, OBIT_Error,"%s requested field %d out of range [0,%d]",
		   routine, n, in->mosaic->numberImages-1);
      return gotSome;
  }

  /* Load/FT Grid CC table */
  gotSome = ObitSkyModelGridFTComps (in, n, uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  return gotSome;
} /* end ObitSkyModelGridComps */

/**
 * Load image model into in plane member and Fourier transform.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * \param in  SkyModel 
 * \param n   Image number on mosaic
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (generally true here).
 */
gboolean ObitSkyModelLoadImage (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err)
{
  gboolean gotSome = FALSE;
  ofloat factor[MAXFARRAYDIM];
  olong ndim, naxis[2];
  ObitFArray *PadImage = NULL;
  ObitCArray *FFTImage = NULL;
  gchar *routine = "ObitSkyModelLoadImage";
  /* DEBUG
  ObitFArray *tempFArray = NULL; */
  /* END DEBUG */
  
  /* error checks - assume most done at higher level */
  if (err->error) return gotSome;

  /* Check Mosaic */
  if (!ObitImageMosaicIsA(in->mosaic)) {
    Obit_log_error(err, OBIT_Error,"%s mosaic member not defined in %s",
		   routine, in->name);
    return gotSome;
  }

  /* Check Image number */
  if ((n<0) || (n>in->mosaic->numberImages)) {
    Obit_log_error(err, OBIT_Error,"%s: Image %d not defined in %s",
		   routine, n, in->mosaic->name);
    return gotSome;
  }

  /* Get Image into in->plane */
  ObitSkyModelgetPBImage (in, uvdata, n, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  /* Clip below in->minFlux */
  ObitFArrayClip (in->plane, in->minFlux, 1.0e20, 0.0);

  /* Scale, Zero pad, deblank image */
  ndim = 2;
  factor[0] = OverSample; factor[1] = OverSample; factor[3] = OverSample;
  naxis[0] = (olong)(factor[0]*in->plane->naxis[0]+0.5);
  naxis[1] = (olong)(factor[1]*in->plane->naxis[1]+0.5);
  PadImage = ObitFArrayCreate ("Padded Image", ndim, naxis);
  ObitFArrayPad (in->plane, PadImage, in->factor);

  /* DEBUG */
  /* ObitImageUtilArray2Image ("DbugPaddedImage.fits", 1, PadImage, err);*/
  /* if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);*/
  /* fprintf(stderr,"After ReadingPaddingImage\n");*/
  /* END DEBUG */
  
  /* Output of FFT */
  ndim = 2;
  naxis[0] = 1+PadImage->naxis[1]/2; naxis[1] = PadImage->naxis[0]; 
  FFTImage = ObitCArrayCreate ("FFT output", ndim, naxis);

  /* Fourier transform to in->FTplane */
  ObitSkyModelFTImage(in, PadImage, FFTImage);

  /* DEBUG */
  /*tempFArray = ObitCArrayMakeF(FFTImage);*/     /* Temp FArray */
  /*ObitCArrayReal (FFTImage, tempFArray);*/      /* Get real part */
  /*ObitImageUtilArray2Image ("DbugImageFFTReal.fits", 1, tempFArray, err);*/
  /*tempFArray = ObitFArrayUnref(tempFArray);*/   /* delete temporary */
  /*if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);*/
  /*fprintf(stderr,"After Image FFT\n");*/
  /* END DEBUG */

  /* Add conjugate columns for interpolator */
  in->numConjCol = HWIDTH;  /* Number of columns on conjugate side of plane */
  in->FTplane = ObitCArrayUnref(in->FTplane);
  in->FTplane = ObitCArrayAddConjg(FFTImage, in->numConjCol);
  
  /* DEBUG */
  /*tempFArray = ObitCArrayMakeF(in->FTplane);*/  /* Temp FArray */
  /*ObitCArrayReal (in->FTplane, tempFArray); */  /* Get real part */
  /*ObitImageUtilArray2Image ("DbugImageConjgReal.fits", 1, tempFArray, err);*/
  /*tempFArray = ObitFArrayUnref(tempFArray); */  /* delete temporary */
  /*if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);*/
  /*fprintf(stderr,"After Image ObitCArrayAddConjg\n");*/
  /* END DEBUG */

  /* Create interpolator */
  /* (re)Create interpolator */
  in->myInterp = ObitCInterpolateUnref(in->myInterp);
  in->myInterp = 
    newObitCInterpolateCreate("UV data interpolator", in->FTplane, 
			      in->mosaic->images[n]->myDesc,
			      factor[0], factor[1], in->numConjCol, HWIDTH, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
  
  /* Cleanup */
  PadImage  = ObitFArrayUnref(PadImage);
  FFTImage  = ObitCArrayUnref(FFTImage);

  gotSome = TRUE;  /* Hard to miss with an image */
  return gotSome;
} /* end ObitSkyModelLoadImage */

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
 * \param in     SkyModel with model components loaded (ObitSkyModelLoad)
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to model and subtract from current buffer
 * \param err Obit error stack object.
 */
void ObitSkyModelFTDFT (ObitSkyModel *in, olong field, ObitUV *uvdata, ObitErr *err)
{
  olong i, nvis, lovis, hivis, nvisPerThread, nThreads;
  FTFuncArg *args;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelFTDFT";

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
    strcpy (args->type, "base");  /* Enter type as first entry */
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
}  /* end ObitSkyModelFTDFT */

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
 * \li in     SkyModel with model components loaded (ObitSkyModelLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number, <0-> no threads
 * \li err Obit error stack object.
 * \return NULL
 */
static gpointer ThreadSkyModelFTDFT (gpointer args)
{
  /* Get arguments from structure */
  FTFuncArg *largs = (FTFuncArg*)args;
  ObitSkyModel *in = largs->in;
  /*olong field      = largs->field;*/
  ObitUV *uvdata   = largs->uvdata;
  olong loVis      = largs->first-1;
  olong hiVis      = largs->last;
  ObitErr *err     = largs->err;

  /* Need to deal with:
     1) nterm = number of spectral terms, iin->nSpecTerm
     2) coef = spectral coefficients, in CCtable as entries at end
     3) any difference in reference frequency for data and spectra
  */

  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, ncomp, mcomp, iterm, nterm;
  olong lrec, nrparm, naxis[2];
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF;
  olong ilocu, ilocv, ilocw;
  ofloat *visData, *ccData, *data, *fscale, specFact;
  ofloat  modReal, modImag;
  ofloat amp, arg, freq2, freqFact, wt=0.0, temp, ll, lll;
  odouble phase, tx, ty, tz, sumReal, sumImag, *freqArr;
  gchar *routine = "ThreadSkyModelFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

  /* Check */
  if (strncmp (largs->type, "base", 4)) {
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
  nterm = in->nSpecTerm;

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
	freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */

	/* Sum over components */
	sumReal = sumImag = 0.0;
	ccData  = data;
	
	/* Sum by model type */
	switch (in->modType) {
	case OBIT_SkyModel_PointMod:     /* Point */
	  /* From the AIPSish QXXPTS.FOR  */
	  for (iComp=0; iComp<mcomp; iComp++) {
	    if (ccData[0]!=0.0) {  /* valid? */
	      tx = ccData[1]*(odouble)visData[ilocu];
	      ty = ccData[2]*(odouble)visData[ilocv];
	      tz = ccData[3]*(odouble)visData[ilocw];
	      phase = freqFact * (tx + ty + tz);
	      sumReal += ccData[0]*cos(phase);
	      sumImag += ccData[0]*sin(phase);
	      /* DEBUG
		 if (iVis==0) {
		 fprintf (stderr,"%s: comp %d real %f imag %f phase %f sum %f %f \n",
		 routine, iComp, ccData[0]*cos(phase), ccData[0]*sin(phase),
		 57.296*phase, sumReal, sumImag); 
		 } */
	    }  /* end if valid */
	    ccData += lcomp;  /* update pointer */
	  } /* end loop over components */
	  break;
	case OBIT_SkyModel_PointModSpec:     /* Point + spectrum */
	  for (iComp=0; iComp<mcomp; iComp++) {
	    if (ccData[0]!=0.0) {  /* valid? */
	      tx = ccData[1]*(odouble)visData[ilocu];
	      ty = ccData[2]*(odouble)visData[ilocv];
	      tz = ccData[3]*(odouble)visData[ilocw];
	      /* Frequency dependent term */
	      lll = ll = log(freqFact);
	      arg = 0.0;
	      for (iterm=0; iterm<nterm; iterm++) {
		arg += ccData[4+iterm] * lll;
		lll *= ll;
	      }
	      specFact = exp(arg);
	      phase = freqFact * (tx + ty + tz);
	      amp = specFact * ccData[0];
	      sumReal += amp*cos(phase);
	      sumImag += amp*sin(phase);
	    }  /* end if valid */
	    ccData += lcomp;  /* update pointer */
	  } /* end loop over components */
	  break;
	case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	  /* From the AIPSish QGASUB.FOR  */
	  freq2 = freqFact*freqFact;    /* Frequency factor squared */
	  for (iComp=0; iComp<mcomp; iComp++) {
	    if (ccData[0]!=0.0) {  /* valid? */
	      arg = freq2 * (ccData[4]*visData[ilocu]*visData[ilocu] +
			     ccData[5]*visData[ilocv]*visData[ilocv] +
			     ccData[6]*visData[ilocu]*visData[ilocv]);
	      if (arg<-1.0e-5) amp = ccData[0] * exp (arg);
	      else amp = ccData[0];
	      tx = ccData[1]*(odouble)visData[ilocu];
	      ty = ccData[2]*(odouble)visData[ilocv];
	      tz = ccData[3]*(odouble)visData[ilocw];
	      phase = freqFact * (tx + ty + tz);
	      sumReal += amp*cos(phase);
	      sumImag += amp*sin(phase);
	    } /* end if valid */
	    ccData += lcomp;  /* update pointer */
	  }  /* end loop over components */
	  break;
	case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum*/
	  freq2 = freqFact*freqFact;    /* Frequency factor squared */
	  for (iComp=0; iComp<mcomp; iComp++) {
	    if (ccData[0]!=0.0) {  /* valid? */
	      /* Frequency dependent term */
	      lll = ll = log(freqFact);
	      arg = 0.0;
	      for (iterm=0; iterm<nterm; iterm++) {
		arg += ccData[7+iterm] * lll;
		lll *= ll;
	      }
	      specFact = exp(arg);
	      arg = freq2 * (ccData[4]*visData[ilocu]*visData[ilocu] +
			     ccData[5]*visData[ilocv]*visData[ilocv] +
			     ccData[6]*visData[ilocu]*visData[ilocv]);
	      if (arg<-1.0e-5) amp = specFact * ccData[0] * exp (arg);
	      else amp = specFact * ccData[0];
	      tx = ccData[1]*(odouble)visData[ilocu];
	      ty = ccData[2]*(odouble)visData[ilocv];
	      tz = ccData[3]*(odouble)visData[ilocw];
	      phase = freqFact * (tx + ty + tz);
	      sumReal += amp*cos(phase);
	      sumImag += amp*sin(phase);
	    } /* end if valid */
	    ccData += lcomp;  /* update pointer */
	  }  /* end loop over components */
	  break;
	case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	  /* From the AIPSish QSPSUB.FOR  */
	  for (iComp=0; iComp<mcomp; iComp++) {
	    if (ccData[0]!=0.0) {  /* valid? */
	      arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				    visData[ilocv]*visData[ilocv]) * ccData[4];
	      arg = MAX (arg, 0.1);
	      amp = ccData[0] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
	      tx = ccData[1]*(odouble)visData[ilocu];
	      ty = ccData[2]*(odouble)visData[ilocv];
	      tz = ccData[3]*(odouble)visData[ilocw];
	      phase = freqFact * (tx + ty + tz);
	      sumReal += amp*cos(phase);
	      sumImag += amp*sin(phase);
	    } /* end if valid */
	    ccData += lcomp;  /* update pointer */
	  }  /* end loop over components */
	  break;
	case OBIT_SkyModel_USphereModSpec:    /* Uniform sphere + spectrum*/
	  for (iComp=0; iComp<mcomp; iComp++) {
	    if (ccData[0]!=0.0) {  /* valid? */
	      /* Frequency dependent term */
	      lll = ll = log(freqFact);
	      arg = 0.0;
	      for (iterm=0; iterm<nterm; iterm++) {
		arg += ccData[4+iterm] * lll;
		lll *= ll;
	      }
	      specFact = exp(arg);
	      arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				    visData[ilocv]*visData[ilocv]) * ccData[4];
	      arg = MAX (arg, 0.1);
	      amp = specFact * ccData[0] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
	      tx = ccData[1]*(odouble)visData[ilocu];
	      ty = ccData[2]*(odouble)visData[ilocv];
	      tz = ccData[3]*(odouble)visData[ilocw];
	      phase = freqFact * (tx + ty + tz);
	      sumReal += amp*cos(phase);
	      sumImag += amp*sin(phase);
	    } /* end if valid */
	    ccData += lcomp;  /* update pointer */
	  }  /* end loop over components */
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
} /* ObitSkyModelFTDFT */

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
 * \param in     SkyModel with model components loaded (ObitSkyModelLoad)
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to model and subtract from current buffer
 * \param err Obit error stack object.
 */
void ObitSkyModelFTGrid (ObitSkyModel *in, olong field, ObitUV *uvdata, ObitErr *err)
{
  olong i, nvis, lovis, hivis, nvisPerThread, nThreads;
  FTFuncArg *args;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelFTGrid";

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
    if (args->field!=field) args->Interp = ObitCInterpolateUnref(args->Interp);
    args->field  = field;
    args->uvdata = uvdata;
    args->first  = lovis;
    args->last   = hivis;
    if (nThreads>1) args->ithread= i;
    else args->ithread = -1;
    args->err    = err;
    /* local copy of interpolator if needed */
    if (!args->Interp) {
      if (i>0) {
	args->Interp = ObitCInterpolateClone(in->myInterp, NULL);
      } else {
	args->Interp = ObitCInterpolateRef(in->myInterp);
      }
      if (err->error) Obit_traceback_msg (err, routine, in->name);
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
}  /* end ObitSkyModelFTGrid */

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
 * \li in     SkyModel with model components loaded (ObitSkyModelLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number, <0-> no threads
 * \li err Obit error stack object.
 * \li Interp UV Interpolator
 * \return NULL
 */
gpointer ThreadSkyModelFTGrid (gpointer args)
{
  /* Get arguments from structure */
  FTFuncArg *largs = (FTFuncArg*)args;
  ObitSkyModel *in = largs->in;
  olong field      = largs->field;
  ObitUV *uvdata   = largs->uvdata;
  olong loVis      = largs->first-1;
  olong hiVis      = largs->last;
  ObitErr *err     = largs->err;
  ObitCInterpolate *Interp = largs->Interp;

  ObitImageDesc *imDesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  olong iVis, iIF, iChannel, iStoke;
  olong i, j, k, lrec, nrparm;
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF;
  olong ilocu, ilocv, ilocw;
  ofloat *visData, *fscale, vis[2], flip;
  ofloat sumReal, sumImag, modReal, modImag;
  ofloat freqFact, wt=0.0, temp;
  ofloat dxyzc[3],  uvw[3], ut, vt, rt, it, fblank = ObitMagicF();
  ofloat umat[3][3], pmat[3][3], rmat[3][3], dmat[3][3];
  ofloat PC, cosPC, sinPC, konst, maprot, uvrot, ssrot, ccrot;
  odouble *freqArr;
  gboolean doRot, doConjg, isBad, do3Dmul, doPC;
  gchar *routine = "ThreadSkyModelFTGrid";

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
	freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */
	
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
	ObitCInterpolateOffset (Interp, uvw, vis, err);
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
} /* ThreadSkyModelFTGrid */

/**
 * Sum the fluxes of components defined by CCVer and endComp
 * \param in  SkyModel Checks InfoList member noNeg
 * \param err Obit error stack object.
 * \return Sum of Clean components
 */
ofloat ObitSkyModelSum (ObitSkyModel *in, ObitErr *err)
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
  gchar *routine = "ObitSkyModelSum";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return sum;
  g_assert (ObitSkyModelIsA(in));

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
    if ((tempTable==NULL) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, sum);
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
} /* end ObitSkyModelSum */

/**
 * Compress CC tables
 * \param in  SkyModel
 * \param err Obit error stack object.
 * \return Sum of Clean components
 */
void ObitSkyModelCompressCC (ObitSkyModel *in, ObitErr *err)
{
  olong field;
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL;
  ObitIOCode retCode;
  gchar *tabType = "AIPS CC";
  olong ver, num, *iptr;
  gchar *routine = "ObitSkyModelCompressCC";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelIsA(in));

  /* Update CCVer on in */
  if (ObitInfoListGetP(in->info, "CCVer",  &type, (gint32*)dim, (gpointer)&iptr)) {
    num = MIN ( in->mosaic->numberImages, dim[0]);
    for (i=0; i<num; i++) in->CCver[i] = iptr[i];
    if (dim[0]<in->mosaic->numberImages) /* If one, use for all */
      for (i=num; i<in->mosaic->numberImages; i++) in->CCver[i] = iptr[0];
  }
  
  /* Loop over fields */
  for (field=0; field<in->mosaic->numberImages; field++) {

    /* Check input table to see if there are any components */
    /* Get CC table */
    ver = in->CCver[field];
    tempTable = newObitImageTable (in->mosaic->images[field],OBIT_IO_ReadOnly, 
				   tabType, &ver, err);
    if ((tempTable==NULL) || (err->error)) Obit_traceback_msg (err, routine, in->name);
    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Open CC table */
    retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);

    /* Close Table */
    retCode = ObitTableCCClose (CCTable, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);

    /* Anything in it? */
    if (CCTable->myDesc->nrow<=0) {
      /* release table  */
      CCTable = ObitTableCCUnref (CCTable);
      continue;
    }

    /* Merge */
    retCode = ObitTableCCUtilMerge (CCTable, CCTable, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
    CCTable = ObitTableCCUnref (CCTable);

    /* Use all components */
    in->startComp[field] = 1;
    in->endComp[field]   = 0; 

  } /* end loop over fields */
  
  return;
} /* end ObitSkyModelCompressCC  */

/**
 * Get input parameters from info member
 * \param in Pointer to the ObitSkyModel .
 * \param err Obit error stack object.
 */
void ObitSkyModelGetInput (ObitSkyModel* in, ObitErr *err)
{
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat rtemp[10];
  olong itemp, *iptr, num;
  gchar tempStr[5];
  union ObitInfoListEquiv InfoReal; 

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelIsA(in));
  if (!ObitInfoListIsA(in->info)) return;
  InfoReal.itg = 0;type = OBIT_oint;

  /* Channel range */
  ObitInfoListGetTest(in->info, "BChan", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  in->startChannel  = itemp;

  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(in->info, "EChan", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (itemp>0) in->numberChannel = itemp - in->startChannel+1;
  else  in->numberChannel = 0;

  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(in->info, "BIF", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  in->startIF  = itemp;

  InfoReal.itg = 0; type = OBIT_oint;
  ObitInfoListGetTest(in->info, "EIF", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (itemp>0) in->numberIF = itemp - in->startIF+1;
  else in->numberIF = 0;

  /* Stokes */
  for (i=0; i<4; i++) tempStr[i] = ' '; tempStr[4] = 0;
  ObitInfoListGetTest(in->info, "Stokes", &type, (gint32*)dim, &tempStr);
  for (i=0; i<4; i++) in->stokes[i] = tempStr[i]; in->stokes[4] = 0;

  /* 3D wanted? */
  InfoReal.itg = (olong)in->do3D; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "do3D", &type, (gint32*)dim, &InfoReal);
  in->do3D = InfoReal.itg;

  /* Division wanted? */
  InfoReal.itg = (olong)in->doDivide; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "DIVIDE", &type, (gint32*)dim, &InfoReal);
  in->doDivide = InfoReal.itg;

  /* Data replacement wanted? */
  InfoReal.itg = (olong)in->doReplace; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "REPLACE", &type, (gint32*)dim, &InfoReal);
  in->doReplace = InfoReal.itg;

  /* Want only positive flux components?? */
  InfoReal.itg = (olong)in->noNeg; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "noNeg", &type, (gint32*)dim, &InfoReal);
  in->noNeg = InfoReal.itg;

  /* Relative Primary Beam correction wanted? */
  InfoReal.itg = (olong)in->doPBCor; type = OBIT_bool;
  ObitInfoListGetTest(in->info, "PBCor", &type, (gint32*)dim, &InfoReal);
  in->doPBCor = InfoReal.itg;

  /* Antenna size for rel. Primary beam correction */
  InfoReal.flt = in->antSize; type = OBIT_float;
  ObitInfoListGetTest(in->info, "antSize", &type, (gint32*)dim, &InfoReal);
  in->antSize = InfoReal.flt;

  /* Model scaling factor */
  InfoReal.flt = 1.0; type = OBIT_float;
  ObitInfoListGetTest(in->info, "Factor", &type, (gint32*)dim, &InfoReal);
  in->factor = InfoReal.flt;

  /* Minimum flux density */
  InfoReal.flt = -1.0e20; type = OBIT_float;
  ObitInfoListGetTest(in->info, "minFlux", &type, (gint32*)dim, &InfoReal);
  in->minFlux = InfoReal.flt;

  /* Model type */
  InfoReal.itg = (olong)OBIT_SkyModel_Comps; type = OBIT_long;
  ObitInfoListGetTest(in->info, "ModelType", &type, (gint32*)dim, &InfoReal);
  in->modelType = InfoReal.itg;

  /* Model mode */
  InfoReal.itg = (olong)OBIT_SkyModel_Fastest; type = OBIT_long;
  ObitInfoListGetTest(in->info, "Mode", &type, (gint32*)dim, &InfoReal);
  in->modelMode = InfoReal.itg;

  /* Point model flux density */
  InfoReal.flt = 0.0; type = OBIT_float;
  ObitInfoListGetTest(in->info, "MODPTFLX", &type, (gint32*)dim, &InfoReal);
  in->pointFlux = InfoReal.flt;
  if (in->pointFlux > 0.0) in->modelType = OBIT_SkyModel_Point;
  if (in->pointFlux > 0.0) in->modelMode = OBIT_SkyModel_DFT;

  /* Point model X offset (deg) */
  InfoReal.flt = 0.0; type = OBIT_float;
  ObitInfoListGetTest(in->info, "MODPTXOF", &type, (gint32*)dim, &InfoReal);
  in->pointXOff = InfoReal.flt;

  /* Point model Y offset (deg) */
  InfoReal.flt = 0.0; type = OBIT_float;
  ObitInfoListGetTest(in->info, "MODPTYOF", &type, (gint32*)dim, &InfoReal);
  in->pointYOff = InfoReal.flt;

  /* Point model other parameters */
  for (i=0; i<10; i++) rtemp[i] = 0;
  type = OBIT_float; dim[0] = 4;
  ObitInfoListGetTest(in->info, "MODPTYPM", &type, (gint32*)dim, &rtemp);
  for (i=0; i<dim[0]; i++) in->pointParms[i] = rtemp[i];

  /* CC tables version, components */
  if ((in->mosaic) && (in->mosaic->numberImages>0)) {
    /* CC version number */
    if (ObitInfoListGetP(in->info, "CCVer",  &type, (gint32*)dim, (gpointer)&iptr)) {
      num = MIN ( in->mosaic->numberImages, dim[0]);
      for (i=0; i<num; i++) in->CCver[i] = iptr[i];
      if (dim[0]<in->mosaic->numberImages) /* If one, use for all */
	for (i=num; i<in->mosaic->numberImages; i++) in->CCver[i] = iptr[0];
    }

    /* Start CC number */
    if (ObitInfoListGetP(in->info, "BComp",  &type, (gint32*)dim, (gpointer)&iptr)) {
      num = MIN ( in->mosaic->numberImages, dim[0]);
      for (i=0; i<num; i++) in->startComp[i] = MAX (1,iptr[i]);
    }

    /* End CC number */
    if (ObitInfoListGetP(in->info, "EComp",  &type, (gint32*)dim, (gpointer)&iptr)) {
      num = MIN ( in->mosaic->numberImages, dim[0]);
      for (i=0; i<num; i++) in->endComp[i] = iptr[i];
    }
  }
  /* Print level */
  in->prtLv = 0;  /* default = none */
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &in->prtLv);

} /* end ObitSkyModelGetInput */

/**
 * Decide which method is the most appropriate to calculate the FT of a model
 * Sets currentMode member function
 * Adopted from the AIPSish QMTYP
 * \param in     Pointer to theObitSkyModel .
 * \param uvdata UV data set
 */
void  ObitSkyModelChose (ObitSkyModel* in, ObitUV* uvdata) 
{
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
  timfft = (timff1 * tpvgrd + timff2 * tfft + timff3 * tpcgrd) * nchan;

  /* How long for a DFT? */
  timdft = tpvpc * nvis * sumcc * nchan;

  if (timdft<=timfft) in->currentMode = OBIT_SkyModel_DFT;
  else in->currentMode = OBIT_SkyModel_Grid;

  /* Must do Grid for Image input model */
  if (in->modelType==OBIT_SkyModel_Image) in->currentMode = OBIT_SkyModel_Grid;
} /* end ObitSkyModelChose */


/**
 * Make sure all data selection values are fully filled in.
 * Especially the polarization controls are difficult.
 * All channel, IF, Stokes axis values are 1-rel.
 * Does a number of consistency checks.
 * Looks up the numbers of Clean components in each image.
 * Recognize the following Stokes'
 * \li "    " use I, RR, LL, RR+LL as given in data
 * \li "RR  " use RR only 
 * \li "LL  " use LL only 
 * \li "RL  " use RL only 
 * \li "LR  " use LR only 
 * \li "RLLR" use RL and LR
 * \li "I???" use I, RR, LL, RR+LL as given in data
 * \li "Q???" use Q, or RL&LR
 * \li "U???" use U, or RL&LR
 * \li "V???" use V, or RR&LL
 * Default action is equivalent to Stokes="    ";
 * Also sets doFlip as necessary for requested combination of
 * Stokes and uv data.
 * \param in     Pointer to theObitSkyModel .
 * \param uvdata UV data set
 * \param err    Obit error stack object.
 */
void ObitSkyModelSetSelect (ObitSkyModel* in, ObitUV* uvdata, ObitErr *err)
{
  ObitUVDesc *uvDesc=NULL;
  ObitTable *tempTable=NULL;
  ObitImageDesc *imDesc=NULL;
  ObitTableCC *CCTable = NULL;
  gchar *tabType = "AIPS CC";
  olong ver;
  ObitIOCode retCode;
  olong field, maxNumber, startPoln, endPoln, iPoln, nPoln;
  gboolean bad, doStok;
  odouble imStoke;
  gchar *routine = "SetSelect";

   /* error checks */
  if (err->error) return;

  uvDesc = uvdata->myDesc;                /* Uv descriptor */
 
  /* IF selection */
  in->startIF = MAX (1, in->startIF);
  if (uvDesc->jlocif>=0) 
    maxNumber = uvDesc->inaxes[uvDesc->jlocif] - in->startIF + 1;
  else {maxNumber = 1;  in->startIF=1;}
  in->numberIF = MIN (maxNumber, in->numberIF);
  if (in->numberIF<=0) in->numberIF = maxNumber;

  /* Channel selection */
  in->startChannel = MAX (1, in->startChannel);
  maxNumber = uvDesc->inaxes[uvDesc->jlocf] - in->startChannel + 1;
  in->numberChannel = MIN (maxNumber, in->numberChannel);
  if (in->numberChannel<=0) in->numberChannel = maxNumber;

  /* If using "point" model assume it is the same as the UV data */
  if (in->mosaic!=NULL) {
    imDesc = in->mosaic->images[0]->myDesc; /* Image descriptor */

    /* Image Stokes */
    imStoke = imDesc->crval[imDesc->jlocs];
    /* Check that Stokes type consistent with Image, i.e. both Stokes  or correlator */
    /* Image Stokes I */
    bad = (imDesc->crval[imDesc->jlocs]==1.0) && (in->stokes[0]=='Q');
    bad = bad || ((imDesc->crval[imDesc->jlocs]==1.0) && (in->stokes[0]=='U'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==1.0) && (in->stokes[0]=='V'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==1.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='L')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==1.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='R')));
    /* Image Stokes Q */
    bad = bad || ((imDesc->crval[imDesc->jlocs]==2.0) && (in->stokes[0]==' '));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==2.0) && (in->stokes[0]=='I'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==2.0) && (in->stokes[0]=='U'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==2.0) && (in->stokes[0]=='V'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==2.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='R')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==2.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='L')));
    /* Image Stokes U */
    bad = bad || ((imDesc->crval[imDesc->jlocs]==3.0) && (in->stokes[0]==' '));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==3.0) && (in->stokes[0]=='I'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==3.0) && (in->stokes[0]=='Q'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==3.0) && (in->stokes[0]=='V'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==3.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='R')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==3.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='L')));
    /* Image Stokes V */
    bad = bad || ((imDesc->crval[imDesc->jlocs]==4.0) && (in->stokes[0]=='Q'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==4.0) && (in->stokes[0]=='U'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==4.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='L')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==4.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='R')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==4.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='R')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==4.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='L')));
    /* Image Stokes R */
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-1.0) && (in->stokes[0]=='Q'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-1.0) && (in->stokes[0]=='U'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-1.0) && (in->stokes[0]=='V'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-1.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='L')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-1.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='R')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-1.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='L')));
    /* Image Stokes L */
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-2.0) && (in->stokes[0]=='Q'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-2.0) && (in->stokes[0]=='U'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-2.0) && (in->stokes[0]=='V'));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-2.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='L')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-2.0) && ((in->stokes[0]=='L') && (in->stokes[1]=='R')));
    bad = bad || ((imDesc->crval[imDesc->jlocs]==-2.0) && ((in->stokes[0]=='R') && (in->stokes[1]=='R')));
    
    if (bad) {
      Obit_log_error(err, OBIT_Error,"%s: Stokes request %s incompatable with %s",
		     routine, in->stokes, in->name);
      return;
    }
  } /* end check image */
  else {
    /* Using point model, assume Stokes Q - doesn't really matter if it's otherwise */
    imStoke = 3.0;
  }

  /* Check that Stokes type consistent with uvdata, i.e. actually in data */
  /* Which first Stokes requested; */
  doStok = FALSE; startPoln = 1; endPoln = MIN (2, uvDesc->inaxes[uvDesc->jlocs]);
  if (!strncmp (in->stokes, "    ",4)) {doStok=FALSE;startPoln=1;endPoln=MIN (2, uvDesc->inaxes[uvDesc->jlocs]);}
  if (!strncmp (in->stokes, "RR  ",4)) {doStok=FALSE;startPoln=1;endPoln=1;}
  if (!strncmp (in->stokes, "LL  ",4)) {doStok=FALSE;startPoln=2;endPoln=2;}
  if (!strncmp (in->stokes, "RL  ",4)) {doStok=FALSE;startPoln=3;endPoln=3;}
  if (!strncmp (in->stokes, "LR  ",4)) {doStok=FALSE;startPoln=4;endPoln=4;}
  if (!strncmp (in->stokes, "RLLR",4)) {doStok=FALSE;startPoln=3;endPoln=4;}
  if (!strncmp (in->stokes, "I",1)) {doStok=TRUE;startPoln=1;endPoln=1;}
  if (!strncmp (in->stokes, "Q",1)) {doStok=TRUE;startPoln=2;endPoln=2;}
  if (!strncmp (in->stokes, "U",1)) {doStok=TRUE;startPoln=3;endPoln=3;}
  if (!strncmp (in->stokes, "V",1)) {doStok=TRUE;startPoln=4;endPoln=4;}
  /* Treat Stokes 'I' and data in correlation mode as Stokes ' ' */
  if (!strncmp (in->stokes, "I",1) || (uvDesc->crval[uvDesc->jlocs]<0.0)) 
    {doStok=FALSE;startPoln=1;endPoln=MIN (2, uvDesc->inaxes[uvDesc->jlocs]);}

  /* Actual poln range in uv data */
  nPoln = uvDesc->inaxes[uvDesc->jlocs];
  if (uvDesc->crval[uvDesc->jlocs]>0.0) {  /* Stokes */
    iPoln = uvDesc->crval[uvDesc->jlocs] + 0.5;
    nPoln = 1;
  } else {                                 /* Correlation */
    iPoln = uvDesc->crval[uvDesc->jlocs] - 0.5;
    nPoln = MIN (2,uvDesc->inaxes[uvDesc->jlocs]);
    /* Treat  single LL as RR */
    if ((iPoln==-2) && (nPoln==1)) iPoln = -1;
  }

  /* Check that Stokes type consistent with uvdata, i.e. actually in data */
  if (doStok) {  /* Request in Stokes */
    if (uvDesc->crval[uvDesc->jlocs]>0.0) {
      /* request in Stokes and data in Stokes */
      bad = (startPoln>iPoln) || (endPoln>iPoln+nPoln-1);
    } else {
      /* request in Stokes and data in correlation */
      bad = (startPoln==1) && (iPoln<-2);
      bad = bad || ((startPoln==2) && ((iPoln<-3) || ((iPoln-nPoln+1)<-4)));
      bad = bad || ((startPoln==3) && ((iPoln<-3) || ((iPoln-nPoln+1)<-4)));
      bad = bad || ((startPoln==4) && ((iPoln<-1) || ((iPoln-nPoln+1)<-2)));
    }
  } else { /* Request in correlation */
    if (uvDesc->crval[uvDesc->jlocs]<0.0) {
      /* request in Correlation and data in correlation */
      bad = (startPoln<abs(iPoln)) || (endPoln>abs(iPoln)+uvDesc->inaxes[uvDesc->jlocs]-1);
    } else {
      /* request in Correlation and data in Stokes */
      bad = (startPoln==1) && (iPoln>2);
      bad = bad || ((startPoln==2) && ((iPoln<-3) || ((iPoln-nPoln+1)<-4)));
      /* Most possibilities don't make much sense */
      bad = bad || (startPoln>2);
    }
  }
  
  if (bad) {
    Obit_log_error(err, OBIT_Error,"%s: Stokes request %s incompatable with %s",
		   routine, in->stokes, uvdata->name);
    return;
  }

  /* Stokes selection - based on Stokes */
  in->startPoln  = startPoln-abs(iPoln)+1;
  in->numberPoln = nPoln;
  
  maxNumber = uvDesc->inaxes[uvDesc->jlocs] - in->startPoln + 1;
  maxNumber = MAX (1, maxNumber);
  in->numberPoln = MIN (maxNumber, in->numberPoln);
  
  /* Need to multiply by sqrt(-1)? UPol model and RL, LR requested */
  in->doFlip = (imStoke==3.0) && (!strncmp (in->stokes, "RLLR",4));

  /* Factor for second Stokes, 1.0 unless model Stokes U and RL,LR or 
     Stokes V and RR,LL data */
  in->stokFactor = 1.0;
  if (in->doFlip) in->stokFactor = -1.0;
  if ((imStoke>=3.0) && /* Image U or V */
      (!strncmp (in->stokes, "RLLR",4)) &&
      (uvDesc->crval[uvDesc->jlocs]<0.0))    /* Data RR,LL,RL,LR */
    in->stokFactor = -1.0;

  /* Find actual number of Clean components */
  if ((in->pointFlux == 0.0) && (in->mosaic!=NULL)) {
    imDesc = in->mosaic->images[0]->myDesc; /* Image descriptor */
    for (field=0; field<in->mosaic->numberImages; field++) {
      /* Open Image */
      /* Use external buffer (Not actually reading image here) */
      in->mosaic->images[field]->extBuffer = TRUE;
      retCode = ObitImageOpen (in->mosaic->images[field], OBIT_IO_ReadOnly, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* Get CC table */
      ver = in->CCver[field];
      tempTable = newObitImageTable (in->mosaic->images[field],OBIT_IO_ReadOnly, 
				     tabType, &ver, err);
      if ((tempTable==NULL) || (err->error)) {
	Obit_log_error(err, OBIT_Error,"%s: Cannot find %s table %d on %s",
		       routine, tabType, ver, in->mosaic->images[field]->name);
	return;
      } /*Obit_traceback_msg (err, routine, in->name);*/
      CCTable = ObitTableCCConvert(tempTable);
      tempTable = ObitTableUnref(tempTable);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      in->CCver[field] = ver;  /* save if defaulted (0) */
      
      /* Open CC table */
      retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* How many components to deal with? */
      if (in->endComp[field]<=0) in->endComp[field] = CCTable->myDesc->nrow;
      
      /* Image size? */
      imDesc = in->mosaic->images[field]->myDesc;
      in->mosaic->nx[field] = imDesc->inaxes[imDesc->jlocr];
      in->mosaic->ny[field] = imDesc->inaxes[imDesc->jlocd];
      
      /* Close Table */
      retCode = ObitTableCCClose (CCTable, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* Release table object */
      CCTable = ObitTableCCUnref (CCTable);
      
      /* Close Image */
      retCode = ObitImageClose (in->mosaic->images[field], err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* Unset use external buffer switch */
      in->mosaic->images[field]->extBuffer = FALSE;
    } /* End checking CCs */
  } /* end if mosaic */

} /* ObitSkyModelSetSelect  */

/**
 * Determines block of channels to use in this pass
 * If not making Primary Beam correctsion then all selected,
 * else the next block for which the primary beam correction 
 * varies by less than 1% at the edge of the FOV.
 * Works by modifying the followint ObitSkyModel members:
 * \li startIFPB   First IF selected in current pass, 
 *     if this is initally -1 then initialize.
 * \li numberIFPB  Number of IFs selected in current pass
 * \li startChannelPB   First channel selected in current pass
 * \li numberChannelPB  Number of channels selected in current pass
 * \param in      SkyModel
 * \param uvdata  UV data
 * \param err      Obit error stack object.
 * \return TRUE if there are more channels to do, else FALSE
 */
gboolean ObitSkyModelsetPBChans(ObitSkyModel* in, ObitUV* uvdata, ObitErr *err)
{
  gboolean done = TRUE;
  ofloat FOV, PBStart, PBFact;
  odouble Angle, sumFreq;
  olong nfreq, nif, niffreq, ifreq, incf, incif, iIF, iChan, countFreq;
  ObitUVDesc *uvDesc = NULL;
  gboolean found;
  gchar *routine = "setPBChans";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitSkyModelIsA(in));
  g_assert (ObitUVIsA(uvdata));

  /* See if primary Beam rel. corrections requested */
  if (!in->doPBCor) {
    /* No - use full set of IFs/channels */
    done = in->startIFPB>0;  /* have we done this before? */
    in->startIFPB       = in->startIF;
    in->numberIFPB      = in->numberIF;
    in->startChannelPB  = in->startChannel;
    in->numberChannelPB = in->numberChannel;
    return done;
  }

  /* everything done? */
  done = (in->startIFPB >= (in->startIF+in->numberIF-1)) && 
    (in->startChannelPB >= (in->startChannel+in->numberChannel-1));
  if (done) return done;

  /* Determine which channels to do this pass */
  /* If in->startIFPB = -1 initialize */
  if (in->startIFPB==-1) {
    in->startIFPB       = in->startIF;
    in->numberIFPB      = 0;
    in->startChannelPB  = in->startChannel;
    in->numberChannelPB = 0;
  } else {

    /* New starting values */
    in->startIFPB      += in->numberIFPB;
    in->startChannelPB += in->numberChannelPB;
    if (in->startChannelPB>in->numberChannel) in->startChannelPB=1;
    in->numberIFPB      = 0;
    in->numberChannelPB = 0;
  }

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
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
  }
  Angle = FOV;  /* Field of view = max angle from point */
  
  /* which frequency channel is the start? */
  ifreq = (in->startIFPB-1) * incif + (in->startChannelPB-1) * incf;
  
  /* Primary beam correction factor at first IF/channel */
  uvDesc = uvdata->myDesc;
  PBStart = ObitPBUtilRelPB (Angle, niffreq, uvDesc->freqArr, in->antSize, 0.0,
			     uvDesc->freqArr[ifreq]);
  
  /* Loop through the remaining channels until finding one whose PB 
     correction differs from PBStart by more than 1% */
  found = FALSE;  /* haven't found the high channel number yet */
  countFreq = 0;
  sumFreq = 0.0;
  for (iIF = in->startIFPB; iIF<=in->startIF+in->numberIF-1; iIF++) {
    for (iChan = in->startChannelPB; iChan<=in->startChannel+in->numberChannel-1; 
	 iChan++) {
      /* which frequency channel is this? */
      ifreq = (iIF-1) * incif + (iChan-1) * incf;
      PBFact = ObitPBUtilRelPB (Angle, niffreq, uvDesc->freqArr, in->antSize, 0.0,
				uvDesc->freqArr[ifreq]);
      /* Does this one differ by more than 1% from PBStart? */
      found = fabs(PBFact-PBStart) > 0.01;
      if (found) break;
      countFreq++;
      sumFreq += uvDesc->freqArr[ifreq];
      in->numberChannelPB++;
    } /* end loop over channels */
    if (found) break;
    in->numberIFPB++;
  } /* end loop over IFs */

  /* Either whole IFs or a block of channels in a single IF */
  in->numberChannelPB = MIN (in->numberChannelPB, in->numberChannel);

  /* if the limit not found counts too high
     if (!found) {
     in->numberIFPB--;
     in->numberChannelPB--;
     } */

  /* Set useful values on object */
  in->nfreqPB = niffreq;             /* Number of channels */
  if (countFreq>0) in->PBFreq  = sumFreq / countFreq; /* average frequency */
  else in->PBFreq  = 1.0;  /* Something very wrong */
  
  /* everything done? */
  done = (in->startIFPB > (in->startIF+in->numberIF-1));
  return done;
} /* end setPBChans */


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
 * \return ObitCCTable to use, this should be Unref when done and 
 *                 Zapped if outCCver != 0.
 */
ObitTableCC* ObitSkyModelgetPBCCTab (ObitSkyModel* in, ObitUV* uvdata, 
				     olong field, olong *inCCVer, olong *outCCver,
				     olong *startCC, olong *endCC, ofloat range[2],
				     ObitErr *err)
{
  ObitTable *tempTable = NULL;
  ObitTableCC *CCTable = NULL, *newCCTable = NULL;
  ObitIOCode retCode;
  gchar *tabType = "AIPS CC";
  olong ver, tiver=0;
  gchar *routine = "ObitSkyModelgetPBCCTab";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return CCTable;
  g_assert (ObitSkyModelIsA(in));

  /* Check input table to see if there are any selected components */
  /* Get CC table */
  ver = *inCCVer;
  tempTable = newObitImageTable (in->mosaic->images[field],OBIT_IO_ReadOnly, 
				 tabType, &ver, err);
  if ((tempTable==NULL) || (err->error)) 
    Obit_traceback_val (err, routine, in->name, CCTable);
  CCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_val (err, routine, in->name, CCTable);
  *inCCVer = ver;  /* save if defaulted (0) */
  tiver = ver;
  
  /* Open CC table */
  retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, CCTable);
  
  /* How many components to deal with? */
  if (*endCC<=0) *endCC = CCTable->myDesc->nrow;
  
  /* Close Table */
  retCode = ObitTableCCClose (CCTable, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, CCTable);
  
  /* See if primary Beam rel. corrections requested or anything to do */
  if ((!in->doPBCor) || (*endCC < *startCC)) {
    /* No - just return input table */

  } else {  /* Want PB corrected table */
    
    /* Release table object - recreate */
    CCTable = ObitTableCCUnref (CCTable);
    
    /* Get PB corrected table */
    CCTable = ObitPBUtilCCCor (in->mosaic->images[field], *inCCVer, outCCver, 
			       in->nfreqPB, uvdata->myDesc->freqArr, in->antSize, 
			       0.0, in->PBFreq, startCC, endCC, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, CCTable);
    tiver = *outCCver; /* which version number is this? */
  } 

  /* Need to compress, select? */
  if (in->currentMode==OBIT_SkyModel_Mixed ) {
    *outCCver =  0;  /* Create new one */
    newCCTable = ObitTableCCUtilMergeSel2Tab (in->mosaic->images[field], tiver, outCCver, 
					      *startCC, *endCC, range, err);
    *startCC = 1;   /* Want all of these */
    *endCC   = newCCTable->myDesc->nrow ;
 
    /* Release old table object */
    CCTable = ObitTableCCUnref (CCTable);
    CCTable = newCCTable;
  }
  
  return CCTable;
} /* end ObitSkyModelgetPBCCTab */

/**
 * Sets the in->plane member to either the pixels from the image in the 
 * specified field in in->mosaic or this array with relative 
 * primary beam corrections if in->doPBCor.
 * \param in       SkyModel
 * \param uvdata   UV data
 * \param field    Field number in in->mosaic
 * \param err      Obit error stack object.
 * \return ObitCCTable to use, this should be Unref when done and 
 *   Zapped if outCCver != 0
 */
void ObitSkyModelgetPBImage (ObitSkyModel* in, ObitUV* uvdata, olong field, 
				    ObitErr *err)
{
  olong   blc[IM_MAXDIM]={1,1,1,1,1}, trc[IM_MAXDIM]={0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong inPlane[IM_MAXDIM]={1,1,1,1,1};
  ObitIOCode retCode;
  ObitIOSize IOSize = OBIT_IO_byPlane;
  olong ndim, naxis[2];
  gchar *routine = "ObitSkyModelgetPBImage";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelIsA(in));

  /* See if primary Beam rel. corrections requested */
  if (!in->doPBCor) {
    /* No - just read input image */
    /* Set BLC,TRC, read whole first plane */
    dim[0] = IM_MAXDIM;
    ObitInfoListPut (in->mosaic->images[field]->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (in->mosaic->images[field]->info, "TRC", OBIT_long, dim, trc, err); 
    dim[0] = 1;
    ObitInfoListPut (in->mosaic->images[field]->info, "IOBy", OBIT_long, dim, &IOSize, err);
    if (err->error) Obit_traceback_msg (err, routine, in->mosaic->name);
    
    /* Use external buffer */
    in->mosaic->images[field]->extBuffer = TRUE;
    
    /* Open Image */
    retCode = ObitImageOpen (in->mosaic->images[field], OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
    
    /* (re)allocate memory for plane */
    ndim = 2;
    naxis[0] = in->mosaic->images[field]->myDesc->inaxes[0];
    naxis[1] = in->mosaic->images[field]->myDesc->inaxes[1];
    if (in->plane!=NULL) in->plane = ObitFArrayRealloc(in->plane, ndim, naxis);
    else in->plane = ObitFArrayCreate("ModelImage", ndim, naxis);
    
    /* Read plane */
    retCode = ObitImageRead (in->mosaic->images[field], in->plane->array, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
    
    /* Close Image */
    retCode = ObitImageClose (in->mosaic->images[field], err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
    
    /* Unset use external buffer switch */
    in->mosaic->images[field]->extBuffer = FALSE;
  } else { /* Get PB corrected image array */
    in->plane = ObitFArrayUnref(in->plane);
    in->plane = 
      ObitPBUtilImageCor (in->mosaic->images[field], inPlane,
			  in->nfreqPB, uvdata->myDesc->freqArr, in->antSize, 
			  0.0, in->PBFreq, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end get PB corrected array */  
} /* end ObitSkyModelgetPBImage */
  
/**
 * Grid components onto in->plane (a zeroed array the twice the size 
 * of the image) and Fourier transformed to in->FTplane.
 * Scaling of components and any tapering is applied.
 * Grid is double size for increased accuracy.
 * For convenience in interpolation, HWIDTH columns are added by 
 * copying from the positive half plane.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * \param in     Pointer to theObitSkyModel .
 * \param field  field number (0-rel) in in->mosaic->images
 * \param uvdata UV data set to model
 * \param err    Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelGridFTComps (ObitSkyModel* in, olong field, ObitUV* uvdata, 
				  ObitErr *err)
{
  gboolean gotSome = FALSE;
  ObitImageDesc *imDesc = NULL;
  olong i, j, nx, ny;
  olong ncomp, ndim, naxis[2];
  ofloat gparm[3], dU, dV, UU, VV, texp;
  ofloat konst, xmaj, xmin, cpa, spa, b1, b2, b3, bb2, bb3;
  ofloat taper, *grid, factor[2];
  gboolean doGaus;
  ObitCArray *FFTImage = NULL;
  gchar *routine = "ObitSkyModelGridFTComps";
  /* DEBUG 
  ObitFArray *tempFArray = NULL;*/
  /* END DEBUG */

  /* error check */
  if (err->error) return gotSome ;

  /* Create grid, sum components into in->plane */
  ObitSkyModelLoadGridComps (in, field, uvdata, gparm, &ncomp, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);

  /* Don't bother if no components requested */
  gotSome = ncomp>=1;
  if (!gotSome) return gotSome;

  /* DEBUG
     ObitImageUtilArray2Image ("DbugGriddedComps.fits", 1, in->plane, err);
     if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
     fprintf(stderr,"After ObitSkyModelLoadGridComps\n"); */
  /* END DEBUG */

  /* Output of FFT */
  ndim = 2;
  naxis[0] = 1+in->plane->naxis[0]/2; naxis[1] = in->plane->naxis[1]; 
  FFTImage = ObitCArrayCreate ("FFT output", ndim, naxis);
  
  /* Fourier Transform image */
  ObitSkyModelFTImage (in, in->plane, FFTImage);

  /* DEBUG
     tempFArray = ObitCArrayMakeF(FFTImage);
     ObitCArrayReal (FFTImage, tempFArray); 
     ObitImageUtilArray2Image ("DbugFFTReal.fits", 1, tempFArray, err);
     tempFArray = ObitFArrayUnref(tempFArray);
     if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
     tempFArray = ObitCArrayMakeF(FFTImage);
     ObitCArrayImag (FFTImage, tempFArray); 
     ObitImageUtilArray2Image ("DbugFFTImag.fits", 1, tempFArray, err);
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
    nx = OverSample*imDesc->inaxes[imDesc->jlocr];
    ny = OverSample*imDesc->inaxes[imDesc->jlocd];

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
  in->FTplane = ObitCArrayUnref(in->FTplane);
  in->FTplane = ObitCArrayAddConjg(FFTImage, in->numConjCol);
  
  /* DEBUG */
  /*tempFArray = ObitCArrayMakeF(in->FTplane);*/  /* Temp FArray */
  /*ObitCArrayReal (in->FTplane, tempFArray);*/   /* Get real part */
  /*ObitImageUtilArray2Image ("DbugConjgReal.fits", 1, tempFArray, err);*/
  /*tempFArray = ObitFArrayUnref(tempFArray); */  /* delete temporary */
  /*if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);*/
  /*fprintf(stderr,"After ObitCArrayAddConjg\n");*/
  /* END DEBUG */

  /* (re)Create interpolator */
  factor[0] = OverSample; factor[1] = OverSample;
  in->myInterp = ObitCInterpolateUnref(in->myInterp);
  in->myInterp = 
    newObitCInterpolateCreate("UV data interpolator", in->FTplane, imDesc,
			      factor[0], factor[1], in->numConjCol, HWIDTH, err);					   
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
  
  /* Cleanup */
  FFTImage  = ObitCArrayUnref(FFTImage);

  return gotSome;
} /* end ObitSkyModelGridFTComps */

/**
 * Create array OverSample times the size of the input image (in->plane) 
 * and sum components onto it.
 * Grid is oversize for increased accuracy.
 * Due to the difference with the FFT ordering for half plane complex 
 * in AIPS and using FFTW, the method here is different.
 * Components are added to a grid which is then FFTed.
 * \param in     Pointer to theObitSkyModel .
 * \param field  field number (0-rel) in in->mosaic->images
 * \param uvdata UV data set to model
 * \param gparm  [out] the parameters of the Gaussians in the table
 *               [-1,-1,-1] => not Gaussian.
 * \param ncomp  Actual number of components in in->comps
 * \param err    Obit error stack object.
 */
void  ObitSkyModelLoadGridComps (ObitSkyModel* in, olong field, ObitUV* uvdata, 
				 ofloat gparm[3], olong *ncomp, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableCC *CCTable = NULL;
  ObitImageDesc *imDesc = NULL;
  ofloat range[2];
  gchar *tabType = "AIPS CC";
  olong outCCVer, ver, first, last, startComp, endComp;
  gchar *routine = "ObitSkyModelLoadGridComps";

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
  CCTable = ObitSkyModelgetPBCCTab (in, uvdata, field, &ver, &outCCVer, 
				    &startComp, &endComp, range, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  in->CCver[field] = ver;  /* save if defaulted (0) */
  
  /* Grid */
  first = startComp;
  last  = endComp;
  /* If noNeg last = last before first negative */
  imDesc = in->mosaic->images[field]->myDesc;
  retCode = ObitTableCCUtilGrid (CCTable, OverSample, 
				 &first, &last, in->noNeg,
				 in->factor, 
				 in->minFlux, in->maxGrid,
				 imDesc, &in->plane, gparm, 
				 ncomp, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) Obit_traceback_msg (err, routine, in->name);

  /* Save values of highest comp - probably bad*/
  if (outCCVer==0) {
    /* no translation of table */
    in->startComp[field] = first;
    in->endComp[field] = last;
  } else {
    /* Translated table with only selected values */
    in->endComp[field] = in->startComp[field] + last-first;
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
  
} /* end ObitSkyModelLoadGridComps */

/**
 * Fourier Transform image array in in->plane, 
 * Half plane complex returned in center-at-the-center order.
 * \param in       the ObitSkyModel .
 * \param inArray  Array to be Transformed.
 * \param outArray Output of FFT, half plane complex
 */
void  ObitSkyModelFTImage (ObitSkyModel* in, ObitFArray *inArray, 
			   ObitCArray *outArray)
{
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

} /* end ObitSkyModelFTImage */

/**
 * Add a new field to a SkyModel
 * Resizes arrays, assumes last in list is new one.
 * \param in       the ObitSkyModel .
 * \param err    Obit error stack object.
 */
void  ObitSkyModelAddField (ObitSkyModel* in, ObitErr *err)
{
  olong newField = 0;
  olong i, oldField, *itemp;

  /* error checks */
  if (err->error) return;
  g_assert (ObitSkyModelIsA(in));

  /* field to add - assume 1 being added and mosaic is already correct */
  oldField = in->mosaic->numberImages-1;
  newField = oldField+1;

  /* Resize/copy arrays */
  itemp = ObitMemAlloc0Name (sizeof(olong)*newField, "SkyModel CCver");
  for (i=0; i<oldField; i++) itemp[i] = in->CCver[i]; itemp[i] = itemp[0]; 
  in->CCver = ObitMemFree(in->CCver);
  in->CCver = itemp;

  itemp = ObitMemAlloc0Name (sizeof(olong)*newField, "SkyModel startComp");
  for (i=0; i<oldField; i++) itemp[i] = in->startComp[i]; itemp[i] = 1; 
  in->startComp = ObitMemFree(in->startComp);
  in->startComp = itemp;

  itemp = ObitMemAlloc0Name (sizeof(olong)*newField, "SkyModel endComp");
  for (i=0; i<oldField; i++) itemp[i] = in->endComp[i]; itemp[i] = 0; 
  in->endComp = ObitMemFree(in->endComp);
  in->endComp = itemp;

} /* end ObitSkyModelAddField */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string SkyModel type, "Base" for base class
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
 * \param err     ObitErr for reporting errors.
 */
void ObitSkyModelGetInfo (ObitSkyModel *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err)
{ 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *OK="OK", *Type="Base";
  olong otemp, numberImages=0;
  gchar *routine = "ObitSkyModelGetInfo";

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


} /* end ObitSkyModelGetInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSkyModelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSkyModelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSkyModelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSkyModelClassInfoDefFn (gpointer inClass)
{
  ObitSkyModelClassInfo *theClass = (ObitSkyModelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSkyModelClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSkyModelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSkyModelGetClass;
  theClass->newObit       = (newObitFP)newObitSkyModel;
  theClass->ObitSkyModelFromInfo = (ObitSkyModelFromInfoFP)ObitSkyModelFromInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitSkyModelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSkyModelClear;
  theClass->ObitInit      = (ObitInitFP)ObitSkyModelInit;
  theClass->ObitSkyModelCreate  = (ObitSkyModelCreateFP)ObitSkyModelCreate;
  theClass->ObitSkyModelInitMod = (ObitSkyModelInitModFP)ObitSkyModelInitMod;
  theClass->ObitSkyModelShutDownMod= (ObitSkyModelShutDownModFP)ObitSkyModelShutDownMod;
  theClass->ObitSkyModelInitModel= (ObitSkyModelInitModelFP)ObitSkyModelInitModel;
  theClass->ObitSkyModelLoad    = (ObitSkyModelLoadFP)ObitSkyModelLoad;
  theClass->ObitSkyModelSubUV   = (ObitSkyModelSubUVFP)ObitSkyModelSubUV;
  theClass->ObitSkyModelDivUV   = (ObitSkyModelDivUVFP)ObitSkyModelDivUV;
  theClass->ObitSkyModelFT      = (ObitSkyModelFTFP)ObitSkyModelFT;
  theClass->ObitSkyModelLoadPoint = (ObitSkyModelLoadPointFP)ObitSkyModelLoadPoint;
  theClass->ObitSkyModelLoadComps = (ObitSkyModelLoadCompsFP)ObitSkyModelLoadComps;
  theClass->ObitSkyModelGridComps = (ObitSkyModelGridCompsFP)ObitSkyModelGridComps;
  theClass->ObitSkyModelLoadImage = (ObitSkyModelLoadImageFP)ObitSkyModelLoadImage;
  theClass->ObitSkyModelFTDFT     = (ObitSkyModelFTDFTFP)ObitSkyModelFTDFT;
  theClass->ObitSkyModelFTGrid    = (ObitSkyModelFTGridFP)ObitSkyModelFTGrid;
  theClass->ObitSkyModelSum           = (ObitSkyModelSumFP)ObitSkyModelSum;
  theClass->ObitSkyModelCompressCC    = (ObitSkyModelCompressCCFP)ObitSkyModelCompressCC;
  theClass->ObitSkyModelGetInput      = (ObitSkyModelGetInputFP)ObitSkyModelGetInput;
  theClass->ObitSkyModelChose         = (ObitSkyModelChoseFP)ObitSkyModelChose;
  theClass->ObitSkyModelSetSelect     = (ObitSkyModelSetSelectFP)ObitSkyModelSetSelect;
  theClass->ObitSkyModelsetPBChans    = (ObitSkyModelsetPBChansFP)ObitSkyModelsetPBChans;
  theClass->ObitSkyModelgetPBCCTab    = (ObitSkyModelgetPBCCTabFP)ObitSkyModelgetPBCCTab;
  theClass->ObitSkyModelgetPBImage    = (ObitSkyModelgetPBImageFP)ObitSkyModelgetPBImage;
  theClass->ObitSkyModelGridFTComps   = (ObitSkyModelGridFTCompsFP)ObitSkyModelGridFTComps;
  theClass->ObitSkyModelLoadGridComps = (ObitSkyModelLoadGridCompsFP)ObitSkyModelLoadGridComps;
  theClass->ObitSkyModelFTImage       = (ObitSkyModelFTImageFP)ObitSkyModelFTImage;
  theClass->ObitSkyModelAddField      = (ObitSkyModelAddFieldFP)ObitSkyModelAddField;
  theClass->ObitSkyModelGetInfo       = (ObitSkyModelGetInfoFP)ObitSkyModelGetInfo;

} /* end ObitSkyModelClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSkyModelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModel *in = inn;
  olong number, i;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread = newObitThread();
  in->info = newObitInfoList();
  in->do3D      = TRUE;
  in->doDivide  = FALSE;
  in->doReplace = FALSE;
  in->doPBCor   = TRUE;
  in->noNeg     = FALSE;
  in->PBFreq    = 1.0e9;
  in->nfreqPB   = 1;
  in->antSize   = 25.0;
  in->startChannel = 1;
  in->numberChannel = 0;
  in->startIF = 1;
  in->numberIF = 0;
  in->startChannelPB = -1;
  in->numberChannelPB = 0;
  in->startIFPB = -1;
  in->numberIFPB = 0;
  in->stokes[0] = ' '; in->stokes[1] = ' '; in->stokes[2] = ' ';
  in->stokes[3] = ' '; in->stokes[4] = 0;
  in->mosaic  = NULL;
  in->plane   = NULL;
  in->FTplane = NULL;
  in->myInterp= NULL;
  in->comps   = NULL;
  in->modelType = OBIT_SkyModel_Comps;
  in->modelMode = OBIT_SkyModel_Fastest;
  in->currentMode = OBIT_SkyModel_DFT;
  in->factor      = 1.0;
  in->stokFactor  = 1.0;
  number = 1;
  in->CCver     = ObitMemAllocName (sizeof(olong)*number, "SkyModel CCver");
  in->startComp = ObitMemAllocName (sizeof(olong)*number, "SkyModel startComp");
  in->endComp   = ObitMemAllocName (sizeof(olong)*number, "SkyModel endComp");
  in->minDFT    = 0.0;
  in->maxGrid   = 1.0e20;
  in->doDFT     = FALSE;
  in->doGrid    = FALSE;
  in->nThreads    = 0;
  in->threadArgs= NULL;
  in->DFTFunc   = NULL;
  in->GridFunc  = NULL;
  in->nSpecTerm = 0;
  for (i=0; i<10; i++) in->pointParms[i] = 0.0;

} /* end ObitSkyModelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSkyModel* cast to an Obit*.
 */
void ObitSkyModelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  FTFuncArg *args;
  ObitSkyModel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = freeObitThread(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->mosaic    = ObitImageMosaicUnref(in->mosaic);
  in->plane     = ObitFArrayUnref(in->plane);
  in->FTplane   = ObitCArrayUnref(in->FTplane);
  in->myInterp  = ObitCInterpolateUnref(in->myInterp);
  in->comps     = ObitFArrayUnref(in->comps); 
  in->CCver     = ObitMemFree(in->CCver); 
  in->startComp = ObitMemFree(in->startComp); 
  in->endComp   = ObitMemFree(in->endComp); 
  if (in->threadArgs) {
    /* Check type - only handle "base" */
    if (!strncmp((gchar*)in->threadArgs[0], "base", 4)) {
      for (i=0; i<in->nThreads; i++) {
	args = (FTFuncArg*)in->threadArgs[i];
	if (args->Interp) ObitCInterpolateUnref(args->Interp);
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
    } /* end if this a "base" threadArg */
  }
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSkyModelClear */

