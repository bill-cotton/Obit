/* To DO
 */
/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
#include "ObitSkyModelVMPoln.h"
#include "ObitTableCCUtil.h"
#include "ObitTablePD.h"
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
#include "ObitSpectrumFit.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVMPoln.c
 * ObitSkyModelVMPoln class function definitions.
 *
 * This class is derived from the #ObitSkyModelVM class
 *
 * This class represents sky models incorporating beam corrections and 
 * their Fourier transforms.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSkyModelVMPoln";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitSkyModelVMGetClass;

/**
 * ClassInfo structure ObitSkyModelVMPolnClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSkyModelVMPolnClassInfo myClassInfo = {FALSE};

/*---------------Private structures----------------*/
/* FT threaded function argument 
 Note: Derived classes MUST have the following entries at the beginning 
 of the corresponding structure */
typedef struct {
  /* type "vmpoln" in this class */
  gchar type[12];
  /* SkyModel with model components loaded (ObitSkyModelLoad) */
  ObitSkyModelVMPoln *in;
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
  /* UV Interpolator for FTGrid NYI */
  ObitCInterpolate *Interp;
  /* VM class entries */
  /* Start time (days) of validity of model */
  ofloat begVMModelTime;
  /* End time (days) of validity of model */
  ofloat endVMModelTime;
  /* Thread copy of parallactic angle arrays */
  ofloat *antParang, *sinPA, *cosPA;
} VMPolnFTFuncArg;
/*---------------Private function prototypes----------------*/
/** Private: FT by DFT, may be overridden in derived class */
void ObitSkyModelVMPolnFTDFT (ObitSkyModelVM *in, olong field, 
				ObitUV *uvdata, ObitErr *err);

/** Private: Initialize newly instantiated object. */
void  ObitSkyModelVMPolnInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSkyModelVMPolnClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSkyModelVMPolnClassInfoDefFn (gpointer inClass);

/** Private: Get Inputs. */
void  ObitSkyModelVMPolnGetInput (ObitSkyModel* inn, ObitErr *err);

/** Private: Threaded FTDFT */
static gpointer ThreadSkyModelVMPolnFTDFT (gpointer arg);

/** Private: Load Components */
gboolean ObitSkyModelVMPolnLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				      ObitErr *err);

/** Private: Get MF frequency information */
static olong PolnMFFreqInfo (ObitSkyModelVMPoln *in, ObitImage *image, ObitUV *uvdata, 
			     ofloat *Alpha,  odouble *AlphaRefF,  
			     odouble **specFreq, olong **specIndex, 
			     gboolean *forceDFT);

/** Private: Calculate visibility model */
static void calcModel (VMPolnFTFuncArg* arg, olong ant1, olong ant2, 
		       olong iChan, ofloat stokes[8], ofloat model[8]);

/** Private: Digest PD Table */
static void digestPDTable (ObitSkyModelVMPoln *in, ObitUV *uvdata, olong PDVer, 
			   ObitErr *err);

/** Private: Fetch/compress CC Table */
static ObitTableCC* 
ObitSkyModelVMPolngetPBCCTab (ObitSkyModel* in, ObitImageMosaic *mosaic, ObitUV* uvdata, 
			      olong field, olong *inCCVer, olong *outCCver,
			      olong *startCC, olong *endCC, ofloat range[2],
			      ObitErr *err);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSkyModelVMPoln* newObitSkyModelVMPoln (gchar* name)
{
  ObitSkyModelVMPoln* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMPolnClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSkyModelVMPoln));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSkyModelVMPolnInit((gpointer)out);

 return out;
} /* end newObitSkyModelVMPoln */

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
void ObitSkyModelVMPolnFromInfo (ObitSkyModel *out, gchar *prefix, ObitInfoList *inList, 
				   ObitErr *err)
{ 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *value=NULL;
  gboolean missing;
  gchar *Type = "Poln";
  gchar *routine = "ObitSkyModelVMPolnFromInfo";
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMPolnClassInit();

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

} /* end ObitSkyModelVMPolnFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSkyModelVMPolnGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMPolnClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSkyModelVMPolnGetClass */

/**
 * Make a deep copy of an ObitSkyModelVMPoln.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSkyModelVMPoln* 
ObitSkyModelVMPolnCopy  (ObitSkyModelVMPoln *in, 
			   ObitSkyModelVMPoln *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gchar *routine = "ObitSkyModelVMCopy";

  /* Copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL) && 
	    /* Don't call yourself */
	    (ParentClass!=(const ObitClassInfo*)&myClassInfo));
  out = ParentClass->ObitCopy (in, out, err);

  /* This class */
  /* Not properly inmplemented */
  g_error("%s not implemented", routine);
  return out;
} /* end ObitSkyModelVMPolnCopy */

/**
 * Creates an ObitSkyModelVMPoln 
 * \param name     An optional name for the object.
 * \param uvData   UV data to be operated on
 *                 info member needs:
 * \li PDVer       PD table version
 * \param mosaicI  IPol mosaic giving one or more images/CC tables
 * \param mosaicQ  QPol mosaic giving one or more images/CC tables
 *                 NULL => not given
 * \param mosaicU  UPol mosaic giving one or more images/CC tables
 *                 NULL => not given
 * \param mosaicV  VPol mosaic giving one or more images/CC tables
 *                 NULL => not given
 * \param err      Obit error/message object
 * \return the new object.
 */
ObitSkyModelVMPoln* 
ObitSkyModelVMPolnCreate (gchar* name, ObitUV *uvdata,
			  ObitImageMosaic* mosaicI,
			  ObitImageMosaic* mosaicQ,
			  ObitImageMosaic* mosaicU,
			  ObitImageMosaic* mosaicV,
			  ObitErr *err)
{
  ObitSkyModelVMPoln* out=NULL;
  ObitImage *image0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, number;
  gboolean forceDFT, forceTemp;
  gchar *routine = "ObitSkyModelVMPolnCreate";

  /* Error tests */
  if (err->error) return out;  /* Previous error */

  /* Create basic structure */
  out = newObitSkyModelVMPoln (name);

  /* Modify for input mosaics */
  /* I */
  out->mosaic  = ObitImageMosaicRef(mosaicI);
  out->mosaicI = ObitImageMosaicRef(mosaicI);
  if ((out->mosaicI) && (out->mosaicI->numberImages>0)) {
    out->nModelStokes = 1;
    number = out->mosaicI->numberImages;
    out->CCver = ObitMemAlloc0 (sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = 0;
    out->startIComp = ObitMemAlloc0 (sizeof(olong)*number);
    out->endIComp   = ObitMemAlloc0 (sizeof(olong)*number);
  }
  /* Q */
  out->mosaicQ = ObitImageMosaicRef(mosaicQ);
  if ((out->mosaicQ) && (out->mosaicQ->numberImages>0)) {
    out->nModelStokes = 2;
    number = out->mosaicQ->numberImages;
    out->CCver = ObitMemAlloc0 (sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = 0;
    out->startQComp = ObitMemAlloc0 (sizeof(olong)*number);
    out->endQComp   = ObitMemAlloc0 (sizeof(olong)*number);
  }
  /* U */
  out->mosaicU = ObitImageMosaicRef(mosaicU);
  if ((out->mosaicU) && (out->mosaicU->numberImages>0)) {
    out->nModelStokes = 3;
    number = out->mosaicU->numberImages;
    out->CCver = ObitMemAlloc0 (sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = 0;
    out->startUComp = ObitMemAlloc0 (sizeof(olong)*number);
    out->endUComp   = ObitMemAlloc0 (sizeof(olong)*number);
  }
  /* V */
  out->mosaicV = ObitImageMosaicRef(mosaicV);
  if ((out->mosaicV) && (out->mosaicV->numberImages>0)) {
    out->nModelStokes = 4;
    number = out->mosaicV->numberImages;
    out->CCver = ObitMemAlloc0 (sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = 0;
    out->startVComp = ObitMemAlloc0 (sizeof(olong)*number);
    out->endVComp   = ObitMemAlloc0 (sizeof(olong)*number);
  }

  /* Ensure uvdata fully instantiated and OK */
  ObitUVFullInstantiate (uvdata, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);

  /* swallow PD table */
  out->PDVer = 1;
  ObitInfoListGetTest(uvdata->info, "PDVer", &type, dim, &out->PDVer);
  digestPDTable (out, uvdata, out->PDVer, err);

  /* Fourier transform routines - DFT only */
  out->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMPolnFTDFT;

  /* Create spectrum info arrays */
  image0 = (ObitImage*)out->mosaicI->images[0];	  
  /* Open/close image to get header */
  ObitImageOpen(image0, OBIT_IO_ReadOnly, err);
  ObitImageClose(image0, err);
  if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);
  image0->image =  ObitFArrayUnref(image0->image);  /* free buffer */
  out->nSpecI = PolnMFFreqInfo (out, image0, uvdata, &out->AlphaI, &out->AlphaIRefF, 
			       &out->specFreqI, &out->specIndexI, &forceDFT);
  if (out->mosaicQ) {
    image0 = (ObitImage*)out->mosaicQ->images[0];	  
  /* Open/close image to get header */
  ObitImageOpen(image0, OBIT_IO_ReadOnly, err);
  ObitImageClose(image0, err);
  if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);
  image0->image =  ObitFArrayUnref(image0->image);  /* free buffer */
  out->nSpecQ = PolnMFFreqInfo (out, image0, uvdata, &out->AlphaQ, &out->AlphaQRefF, 
				 &out->specFreqQ, &out->specIndexQ, &forceTemp);
  }
  if (out->mosaicU) {
    image0 = (ObitImage*)out->mosaicQ->images[0];	  
    /* Open/close image to get header */
    ObitImageOpen(image0, OBIT_IO_ReadOnly, err);
    ObitImageClose(image0, err);
    if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);
    image0->image =  ObitFArrayUnref(image0->image);  /* free buffer */
    out->nSpecU = PolnMFFreqInfo (out, image0, uvdata, &out->AlphaU, &out->AlphaURefF, 
				  &out->specFreqU, &out->specIndexU, &forceTemp);
  }
  if (out->mosaicV) {
    image0 = (ObitImage*)out->mosaicV->images[0];	  
    /* Open/close image to get header */
    ObitImageOpen(image0, OBIT_IO_ReadOnly, err);
    ObitImageClose(image0, err);
    if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);
    image0->image =  ObitFArrayUnref(image0->image);  /* free buffer */
    out->nSpecQ = PolnMFFreqInfo (out, image0, uvdata, &out->AlphaV, &out->AlphaVRefF, 
				  &out->specFreqV, &out->specIndexV, &forceTemp);
  }
  
  /* Force DFT? */
  if (forceDFT) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Input freq out of table range, force DFT method",
		   routine);
    out->modelMode = OBIT_SkyModel_DFT;
  }

  /* Grid not implemented, force DFT */
  if (out->modelMode != OBIT_SkyModel_DFT)
    Obit_log_error(err, OBIT_InfoWarn, "%s: Forceing DFT method, Grid not implemented",
		   routine);
  out->modelMode = OBIT_SkyModel_DFT;

  return out;
} /* end ObitSkyModelVMPolnCreate */

/**
 * Initializes Sky Model
 * Checks that data contain RR, LL , save calibration/selection request
 * and set uv data for no selection/calibration
 * \param in      SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMPolnInitMod (ObitSkyModel* inn, ObitUV *uvdata, 
				ObitErr *err)
{
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong iver;
  gboolean btemp;
  olong itemp, numAntList;
  ObitTableList *list=NULL;
  ObitTableAN *TableAN=NULL;
  ObitTableCC *CCTable = NULL;
  ObitTable *tempTable=NULL;
  ObitImageMosaic *mosaic;
  ObitIOCode retCode;
  gchar *tabType = "AIPS CC";
  olong *StartComp, *EndComp;
  ofloat phase=0.5, cp, sp;
  /*gchar *blank="    ";*/
  olong i, j;
  VMPolnFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMPolnInitMod";

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

  /* How many threads? */
  in->nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array */
  if (in->threadArgs==NULL) {
    in->threadArgs = g_malloc0(in->nThreads*sizeof(VMPolnFTFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(VMPolnFTFuncArg)); 
  
    for (i=0; i<in->nThreads; i++) {
      args = (VMPolnFTFuncArg*)in->threadArgs[i];
      strcpy (args->type, "vmpoln");  /* Enter type as first entry */
      args->in     = in;
      args->uvdata = uvdata;
      args->ithread = i;
      args->err    = err;
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      args->begVMModelTime = -1.0e20;
      args->endVMModelTime = -1.0e20;
      args->antParang      = g_malloc0(in->numAnt*sizeof(ofloat));
      args->sinPA          = g_malloc0(in->numAnt*sizeof(ofloat));
      args->cosPA          = g_malloc0(in->numAnt*sizeof(ofloat));
    }
  } /* end initialize */

  /* Call parent initializer */
  ObitSkyModelVMInitMod(inn, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Check requested Stokes
  Obit_return_if_fail((!strncmp(in->stokes,"    ",4)), err,
		      "%s: Unsupported Stokes %s", routine, in->stokes); */

  /* Loop over Stokes */
  for (j=0; j<in->nModelStokes; j++) {
    /* Set pointers for this Stokes' */
    if (j==0) {           /* Stokes I */
      mosaic    = in->mosaicI;
      StartComp = in->startIComp;
      EndComp   = in->endIComp;
    } else if (j==1) {    /* Stokes Q */
      mosaic    = in->mosaicQ;
      StartComp = in->startQComp;
      EndComp   = in->endQComp;
    } else if (j==2) {    /* Stokes U */
      mosaic    = in->mosaicU;
      StartComp = in->startUComp;
      EndComp   = in->endUComp;
    } else if (j==4) {    /* Stokes V */
      mosaic    = in->mosaicV;
      StartComp = in->startVComp;
      EndComp   = in->endVComp;
    } /* End select Stokes */
    /* Check start and end component numbers */
    for (i=0; i<mosaic->numberImages; i++) {
      
      /* Make sure fully instantiated */
      ObitImageFullInstantiate(mosaic->images[i], TRUE, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      /* Get CC table */
      iver = in->CCver[i];
      tempTable = newObitImageTable (mosaic->images[i], OBIT_IO_ReadOnly, 
				     tabType, &iver, err);
      if (tempTable==NULL) {
	Obit_log_error(err, OBIT_Error,"%s: Error finding %s Table %d for %s", 
		       routine, tabType, iver,  mosaic->images[i]->name);
	return;
      }
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      CCTable = ObitTableCCConvert(tempTable);
      tempTable = ObitTableUnref(tempTable);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      /* Open */
      retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      
      /* How many? */
      EndComp[i] = MIN (EndComp[i], CCTable->myDesc->nrow);
      if (EndComp[i]<=0) EndComp[i] = CCTable->myDesc->nrow;
      StartComp[i] = MIN (StartComp[i], CCTable->myDesc->nrow+1);
      StartComp[i] = MAX (StartComp[i], 1);
      
      /* Close */
      retCode = ObitTableCCClose (CCTable, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_msg (err, routine, in->name);
      CCTable = ObitTableUnref(CCTable);
    } /* End loop over images */
  } /* end loop over stokes */

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

  /* Fourier transform routines - DFT only */
  in->DFTFunc   = (ObitThreadFunc)ThreadSkyModelVMPolnFTDFT;

} /* end ObitSkyModelVMPolnInitMod */

/**
 * Any shutdown operations needed for a model
 * Restore calibration/selection state
 * \param in  SkyModel to shutdown
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMPolnShutDownMod (ObitSkyModel* inn, ObitUV *uvdata,
				    ObitErr *err)
{
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i;
  VMPolnFTFuncArg *args;

  /* Call parent shutdown */
  ObitSkyModelVMShutDownMod(inn, uvdata, err);

  if (in->threadArgs) {
    /* Check type - only handle "vmpoln" */
    args = (VMPolnFTFuncArg*)in->threadArgs[0];
    if ((strlen(args->type)>6) || (!strncmp(args->type, "vmpoln", 6))) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMPolnFTFuncArg*)in->threadArgs[i];
	if (args->antParang) g_free(args->antParang);
	if (args->sinPA)     g_free(args->sinPA);
	if (args->cosPA)     g_free(args->cosPA);
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
    } /* end if this a "vmpoln" threadArg */
  }

  /* Restore calibration state */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "doCalSelect", OBIT_bool, dim, &in->saveDoCalSelect);
  ObitInfoListAlwaysPut (uvdata->info, "doCalib", OBIT_long, dim, &in->saveDoCalib);
  dim[0] = 4; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "Stokes", OBIT_string, dim, in->saveStokes);

} /* end ObitSkyModelVMPolnShutDownMod */

/**
 * Initializes an ObitSkyModel for a pass through data in time order.
 * Resets current times, converts field offsets of components to pointing offsets
 * \param in  SkyModel to initialize
 * \param err Obit error stack object.
 */
void ObitSkyModelVMPolnInitModel (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  olong i;
  odouble RAPnt, DecPnt;
  VMPolnFTFuncArg *args;

  /*  Reset time of current model */
  in->begVMModelTime = -1.0e20;
  in->endVMModelTime = -1.0e20;
  in->curVMModelTime = -1.0e20;
   /* in->modelMode = OBIT_SkyModel_DFT; Only can do DFT */

  /* Threading */
  if (in->threadArgs) {
    /* Check type - only handle "vmpoln" */
    args = (VMPolnFTFuncArg*)in->threadArgs[0];
    if ((strlen(args->type)>6) || (!strncmp(args->type, "vmpoln", 6))) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMPolnFTFuncArg*)in->threadArgs[i];
	args->begVMModelTime = -1.0e20;
	args->endVMModelTime = -1.0e20;
      }
    } /* end if this a "vmpoln" threadArg */
  }

  /* Get pointing position */
  ObitImageDescGetPoint(in->mosaicI->images[0]->myDesc, &RAPnt, &DecPnt);

} /* end ObitSkyModelVMPolnInitModel */

/**
 * Update VM model with time and spatial modifications to model
 * Update parallactic angle 
 * Sets begVMModelTime to current time
 * Sets endVMModelTime for when parallactic angle differs by 1 degree.
 * \param in      SkyModelVM 
 * \param time    current time (d)
 * \param suba    0-rel subarray number
 * \param uvdata  uv data being modeled.
 * \param ithread which thread (0-rel) , <0-> no threads
 * \param err Obit error stack object.
 */
void ObitSkyModelVMPolnUpdateModel (ObitSkyModelVM *inn, 
				    ofloat time, olong suba,
				    ObitUV *uvdata, olong ithread, 
				    ObitErr *err)
{
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  olong lithread, ia;
  ofloat curPA, tPA, tTime, bTime;
  VMPolnFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMPolnUpdateModel";

  if (err->error) return;

  /* Thread to update */
  lithread = MAX (0, ithread);
  args = (VMPolnFTFuncArg*)in->threadArgs[lithread];
  /* g_assert (!strncmp(args->type,"vmpoln",6));  Test arg */

   /* Check subarray */
  Obit_return_if_fail(((suba>=0) && (suba<in->numAntList)), err,
		      "%s: bad subarray number %d", routine, suba+1);

  /* Need Parallactic angle */
  tTime = time - 5.0*suba;  /* Correct time for subarray offset*/
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
  for (ia=0; ia<in->numAnt; ia++) {
    args->antParang[ia] = 
      ObitAntennaListParAng (in->AntList[suba], ia+1, bTime, in->curSource);
    args->sinPA[ia] = sin(args->antParang[ia]);
    args->cosPA[ia] = cos(args->antParang[ia]);
  }
} /* end ObitSkyModelVMPolnUpdateModel */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string SkyModel type, "Squint" for this class
 * \param err     ObitErr for reporting errors.
 */
void ObitSkyModelVMPolnGetInfo (ObitSkyModel *inn, gchar *prefix, 
				  ObitInfoList *outList, ObitErr *err)
{ 
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *Type="Squint";
  gchar *routine = "ObitSkyModelVMPolnGetInfo";

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

} /* end ObitSkyModelVMPolnGetInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSkyModelVMPolnClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSkyModelVMPolnClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSkyModelVMPolnClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSkyModelVMPolnClassInfoDefFn (gpointer inClass)
{
  ObitSkyModelVMPolnClassInfo *theClass = (ObitSkyModelVMPolnClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSkyModelVMPolnClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSkyModelVMPolnClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSkyModelVMPolnGetClass;
  theClass->newObit       = (newObitFP)newObitSkyModelVMPoln;
  theClass->ObitCopy      = (ObitCopyFP)ObitSkyModelVMPolnCopy;
  theClass->ObitClear     = (ObitClearFP)ObitSkyModelVMPolnClear;
  theClass->ObitInit      = (ObitInitFP)ObitSkyModelVMPolnInit;
  theClass->ObitSkyModelCreate       = (ObitSkyModelCreateFP)ObitSkyModelVMPolnCreate;
  theClass->ObitSkyModelInitMod      = (ObitSkyModelInitModFP)ObitSkyModelVMPolnInitMod;
  theClass->ObitSkyModelShutDownMod  = (ObitSkyModelShutDownModFP)ObitSkyModelVMPolnShutDownMod;
  theClass->ObitSkyModelInitModel    = (ObitSkyModelInitModelFP)ObitSkyModelVMPolnInitModel;
  theClass->ObitSkyModelFTDFT        = (ObitSkyModelFTDFTFP)ObitSkyModelVMPolnFTDFT;
  theClass->ObitSkyModelGetInput     = (ObitSkyModelGetInputFP)ObitSkyModelVMPolnGetInput;
  theClass->ObitSkyModelChose        = (ObitSkyModelChoseFP)ObitSkyModelVMPolnChose;
  theClass->ObitSkyModelLoadComps    = (ObitSkyModelLoadCompsFP)ObitSkyModelVMPolnLoadComps;
  theClass->ObitSkyModelVMUpdateModel=
    (ObitSkyModelVMUpdateModelFP)ObitSkyModelVMPolnUpdateModel;
  theClass->ObitSkyModelGetInfo= (ObitSkyModelGetInfoFP)ObitSkyModelVMPolnGetInfo;

} /* end ObitSkyModelVMPolnClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Load components model for all selected I,Q,U,V into in 
 * VM?Comps members.
 * Does nModelStokes of (I, Q, U, V)
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Allows mixed points and Gaussians, one per CC Table
 * Adapted from the AIPSish QNOT:IONDFT
 * \param inn  SkyModel 
 * \param n   Image number on mosaic, if -1 load all images
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * Output is in members VM?Comps, the entries are
 * \li 0 Field (0-rel)
 * \li 1 CC DeltaX
 * \li 2 CC DeltaY
 * \li 3 Amplitude (Jy)
 * \li 4 -2*pi*x (radians)
 * \li 5 -2*pi*y (radians)
 * \li 6 -2*pi*z (radians)
 * \li Other model parameters depending on model type
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelVMPolnLoadComps (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
				      ObitErr *err)
{
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  ObitImageDesc* imIODesc;
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *imDesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  ObitFArray *CompArr=NULL;
  ObitSkyModelCompType modType, maxModType=OBIT_SkyModel_PointMod;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong warray, larray, iterm, nterm, nspec=0, maxTerm=1, toff, iModStok;
  ofloat *array, parms[20], range[2], gp1=0., gp2=0., gp3=0.;
  olong ver, i, j, hi, lo, count, ncomp, startComp, endComp, irow, lrec;
  olong outCCVer, ndim, naxis[2], lenEntry;
  ofloat *table, xxoff, yyoff, zzoff;
  ofloat konst, konst2, xyz[3], xp[3], umat[3][3], pmat[3][3];
  ofloat ccrot, ssrot, xpoff, ypoff, maprot, uvrot;
  ofloat dxyzc[3], cpa, spa, xmaj, xmin;
  gboolean doCheck=FALSE, want, do3Dmul;
  gchar *tabType = "AIPS CC";
  olong *StartComp, *EndComp, *numComp, *nSpec=NULL;
  ObitFArray **VMComps;
  ObitImageMosaic *mosaic;
  gpointer fitArg=NULL;
  ofloat *fitSigma=NULL, *fitParms=NULL;
  ofloat *specCorr=NULL, priorAlpha;
  odouble priorAlphaRefF, *specFreq;
  gchar *routine = "ObitSkyModelVMPolnLoadComps";
  
  /* error checks */
  if (err->error) return gotSome;

  /* UV descriptor */
  uvDesc = uvdata->myDesc;
  
  konst = DG2RAD * 2.0 * G_PI;
  /* konst2 converts FWHM(deg) to coefficients for u*u, v*v, u*v */
  /*konst2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5);*/
  konst2 = DG2RAD * sqrt(2)*G_PI / 2.35482044;
  
  /* Loop over model Stokes */
  for (iModStok=0; iModStok<in->nModelStokes; iModStok++) {

    /* Set pointers for this Stokes' */
    if (iModStok==0) {           /* Stokes I */
      mosaic    = in->mosaicI;
      StartComp = in->startIComp;
      EndComp   = in->endIComp;
      VMComps   = &in->VMIComps;
      numComp   = &in->numIComp;
      nSpec     = &in->nSpecI;
      specFreq  = in->specFreqI;
      priorAlpha = in->AlphaI; priorAlphaRefF = in->AlphaIRefF;
    } else if (iModStok==1) {    /* Stokes Q */
      mosaic    = in->mosaicQ;
      StartComp = in->startQComp;
      EndComp   = in->endQComp;
      VMComps   = &in->VMQComps;
      numComp   = &in->numQComp;
      nSpec     = &in->nSpecQ;
      specFreq  = in->specFreqQ;
      priorAlpha = in->AlphaQ; priorAlphaRefF = in->AlphaQRefF;
    } else if (iModStok==2) {    /* Stokes U */
      mosaic    = in->mosaicU;
      StartComp = in->startUComp;
      EndComp   = in->endUComp;
      VMComps   = &in->VMUComps;
      numComp   = &in->numUComp;
      nSpec     = &in->nSpecU;
      specFreq  = in->specFreqU;
      priorAlpha = in->AlphaU; priorAlphaRefF = in->AlphaURefF;
    } else if (iModStok==4) {    /* Stokes V */
      mosaic    = in->mosaicV;
      StartComp = in->startVComp;
      EndComp   = in->endVComp;
      VMComps   = &in->VMVComps;
      numComp   = &in->numVComp;
      nSpec     = &in->nSpecV;
      specFreq  = in->specFreqV;
      priorAlpha = in->AlphaV; priorAlphaRefF = in->AlphaVRefF;
    } /* End select Stokes */

    /* Don't bother if no components requested */
    if ((n>=0) && (StartComp[n]>EndComp[n])) continue;
    
    /* Loop over images counting CCs */
    count = 0;
    in->modType = OBIT_SkyModel_Unknown; /* Model type not known */
    if (mosaic) {lo = 0; hi = mosaic->numberImages-1;}
    else {lo = 0; hi = 0;}
    if (n>=0) {lo = n; hi = n;}
    for (i=lo; i<=hi; i++) {
      
      /* Expect anything in this table? */
      if ((StartComp[i]>EndComp[i]) || (EndComp[i]<=0)) continue;
      
      /* Get CC table */
      ver = in->CCver[i];
      tempTable = newObitImageTable (mosaic->images[i],OBIT_IO_ReadOnly, 
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
      endComp = EndComp[i];
      if (endComp<=0) endComp = CCTable->myDesc->nrow;
      count += MIN(CCTable->myDesc->nrow, endComp) - MAX(1, in->startComp[i]) + 1;
      
      /* Get model type in first with components */
      /* If only 3 col, or parmsCol 0 size then this is a point model */
      if ((CCTable->myDesc->nfield==3) || 
	  (CCTable->parmsCol<0) ||
	  (CCTable->myDesc->dim[CCTable->parmsCol]<=0)) {
	if (StartComp[i]<=endComp) {
	  in->modType = OBIT_SkyModel_PointMod;
	  maxModType =  MAX (maxModType, in->modType);
	}
	/* Check type of all tables with components to subtract */
      } else if (StartComp[i]<=endComp) {
	/* Create table row */
	CCRow = newObitTableCCRow (CCTable);
	/* Read first */
	irow = StartComp[i];
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
      
      /* Is spectral information included  */
      if (!strncmp (mosaic->images[i]->myDesc->ctype[in->mosaic->images[i]->myDesc->jlocf], 
		    "SPECLOGF", 8)) {
	nterm = in->mosaic->images[i]->myDesc->inaxes[in->mosaic->images[i]->myDesc->jlocf];
	ObitInfoListGetTest (in->mosaic->images[i]->myDesc->info, "NTERM", &type, dim, &nterm);
	maxTerm = MAX (maxTerm, nterm);
	in->nSpecTerm = nterm -1;  /* Only higher order terms */
      }
      
      /* Is tabulated spectral information included? */
      if (!strncmp (mosaic->images[i]->myDesc->ctype[mosaic->images[i]->myDesc->jlocf], 
		    "SPECLNMF", 8)) {
	/* IO descriptor give true size */
	imIODesc = (ObitImageDesc*)mosaic->images[0]->myIO->myDesc;
	nspec = imIODesc->inaxes[imIODesc->jlocf];
	ObitInfoListGetTest (in->mosaic->images[i]->myDesc->info, "NSPEC", &type, dim, &nspec);
	maxTerm = MAX (nspec, maxTerm);
	*nSpec = nspec;  /* Number of terms in the spectrum */
      }
    } /* end loop counting CCs */
    
    /* Use mode type, nterm of the highest encountered */
    in->modType   = maxModType;
    in->nSpecTerm = MAX (0, maxTerm-1);
    *nSpec   = MAX (0, nspec);
    
    /* (re)allocate structure */
    ndim = 2;
    naxis[0] = 7; naxis[1] = count;
    if (in->modType==OBIT_SkyModel_GaussMod)       naxis[0] += 3; /* Gaussian */
    if (in->modType==OBIT_SkyModel_GaussModSpec)   naxis[0] += 3; /* Gaussian + spectrum */
    if (in->modType==OBIT_SkyModel_USphereMod)     naxis[0] += 2; /* Uniform sphere */
    if (in->modType==OBIT_SkyModel_USphereModSpec) naxis[0] += 2; /* Uniform sphere + spectrum*/
    
    /* Any spectral terms */
    naxis[0] += in->nSpecTerm;
    naxis[0] += *nSpec;
    /* Plus component spectral index if needed */
    if (*nSpec>0) naxis[0]++;
    lenEntry = naxis[0];  /* Length of table entry */
    
    if ((*VMComps)!=NULL) *VMComps = ObitFArrayRealloc(*VMComps, ndim, naxis);
    else                  *VMComps = ObitFArrayCreate("Components", ndim, naxis);
    lrec = naxis[0]; /* Save size of entry */
    /* Get pointer */
    naxis[0] = 0; naxis[1]=0; 
    table = ObitFArrayIndex(*VMComps, naxis);
    
    /* Fitting component spectra? */
    if (*nSpec>1) {
      fitSigma = g_malloc0(*nSpec*sizeof(ofloat));
      for (i=0; i<*nSpec; i++) fitSigma[i] = 0.0001;  /* Comp spectrum fitting sigmas (Jy/bm) */
      fitArg = ObitSpectrumFitMakeArg (*nSpec, 2, specFreq[0], specFreq, 
				       FALSE, &fitParms, err);
      if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
    } else {   /* Default spectral index = 0 */
      fitParms = g_malloc(2*sizeof(ofloat));
      fitParms[0] = fitParms[1] = 0.0;
    }
    
    /* Spectral corrections for prior alpha arrays */
    specCorr = g_malloc0(*nSpec*sizeof(ofloat));
    if (priorAlpha!=0.0) {
      for (i=0; i<*nSpec; i++) {
	specCorr[i] = pow((in->specFreqI[i]/priorAlphaRefF), priorAlpha);
      }
    } else { /* No correction */
      for (i=0; i<*nSpec; i++) specCorr[i] = 1.0;
    }

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
      CCTable = 
	ObitSkyModelVMPolngetPBCCTab (inn, mosaic, uvdata, (olong)i, &ver, &outCCVer, 
				      &startComp, &endComp, range, err); 
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      
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
      imDesc = mosaic->images[i]->myDesc; /* Image descriptor */
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
      
      /* Gaussian parameters */
      if ((modType==OBIT_SkyModel_GaussMod) || (modType==OBIT_SkyModel_GaussModSpec)) {
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
	  
	  /* Fitting component spectra? */
	  if (*nSpec>1) {
	    for (iterm=0; iterm<*nSpec; iterm++) table[iterm] = array[iterm+toff]*specCorr[iterm];
	    ObitSpectrumFitSingleArg (fitArg, table, fitSigma, fitParms);
	    /* Sanity check */
	    if (fitParms[1]<-2.0) fitParms[1] = 0.0;
	    if (fitParms[1]> 2.0) fitParms[1] = 0.0;
	  } /* end fit component spectrum */

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
	  } else {  /* no rotation  */
	    xyz[0] = ccrot * xp[0] + ssrot * xp[1];
	    xyz[1] = ccrot * xp[1] - ssrot * xp[0];
	    xyz[2] = xp[2];
	  }
	  table[4] = xyz[0] + xxoff;
	  table[5] = xyz[1] + yyoff;
	  table[6] = xyz[2] + zzoff;
	  
	  /* Zero rest in case */
	  for (iterm=7; iterm<lenEntry; iterm++) table[iterm] = 0.0;
	  
	  /* Only same type as highest */
	  if  (in->modType==modType) {
	    /* Only Point */
	    if (modType==OBIT_SkyModel_PointMod){
	      /* Nothing special this case */
	      
	      /* Only Point with spectrum */
	    } else if (modType==OBIT_SkyModel_PointModSpec) {
	      for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+7] = array[iterm+toff];
	      
	      /* Only Point with tabulated spectrum */
	    } else if (modType==OBIT_SkyModel_PointModTSpec) {
	      table[7] = fitParms[1];   /* Component spectral index */
	      for (iterm=0; iterm<*nSpec; iterm++) 
		table[iterm+8] = array[iterm+toff]*specCorr[iterm] * in->factor;
	    
	      /* Only Gaussian */
	    } else if (in->modType==OBIT_SkyModel_GaussMod) {
	      table[7] = gp1;
	      table[8] = gp2;
	      table[9] = gp3;
	      
	      /* Only Gaussian + spectrum */
	    } else if (in->modType==OBIT_SkyModel_GaussModSpec) {
	      table[7] = gp1;
	      table[8] = gp2;
	      table[9] = gp2;
	      /*  spectrum */
	      for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+10] = array[iterm+toff];
	      
	      /* Only Gaussian with tabulated spectrum */
	    } else if (modType==OBIT_SkyModel_GaussModTSpec) {
	      table[10] = fitParms[1];   /* Component spectral index */
	      for (iterm=0; iterm<*nSpec; iterm++) 
		table[iterm+11] = array[iterm+toff]*specCorr[iterm] * in->factor;
	    
	      /* Only Uniform sphere */
	    } else if (in->modType==OBIT_SkyModel_USphereMod) {
	      table[0] = 3.0 * array[0] * in->factor;
	      table[7] = parms[1]  * 0.109662271 * 2.7777778e-4;
	      table[8] = 0.1;
	      
	      /* Only Uniform sphere+ spectrum */
	    } else if (in->modType==OBIT_SkyModel_USphereModSpec) {
	      table[0] = 3.0 * array[0] * in->factor;
	      table[7] = parms[1]  * 0.109662271 * 2.7777778e-4;
	      table[8] = 0.1;
	      /*  spectrum */
	      for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+9] = array[iterm+toff];
	    }
	  } else { /* Mixed type - zero unused model components */
	    
	    /* Only Point here but some Gaussian - zero Gauss comps */
	    if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_GaussMod)) {
	      table[7] = 0.0;
	      table[8] = 0.0;
	      table[9] = 0.0;
	      
	      /* Gauss here but also some points */
	    } else if ((modType==OBIT_SkyModel_GaussMod) && (in->modType==OBIT_SkyModel_PointMod)) {
	      table[7] = gp1;
	      table[8] = gp2;
	      table[9] = gp3;
	      
	      /* Only PointSpectrum here but some GaussianSpectrum - zero Gauss comps */
	    } else if ((modType==OBIT_SkyModel_PointModSpec) && (in->modType==OBIT_SkyModel_GaussModSpec)) {
	      table[7] = 0.0;
	      table[8] = 0.0;
	      table[9] = 0.0;
	      for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+10] = array[iterm+toff];
	      
	      /* GaussianSpectrum here but also some PointSpectrum */
	    } else if ((in->modType==OBIT_SkyModel_PointModSpec) && (modType==OBIT_SkyModel_GaussModSpec)) {
	      table[7] = gp1;
	      table[8] = gp2;
	      table[9] = gp3;
	      /*  spectrum */
	      for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+10] = array[iterm+toff];
	      
	      /* Only Point here but some with spectrum - zero spectra (Unlikely) */
	    } else if ((modType==OBIT_SkyModel_PointMod) && (in->modType==OBIT_SkyModel_PointModSpec)) {
	      for (iterm=0; iterm<in->nSpecTerm; iterm++) table[iterm+7] = 0.0;
	      
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
    
    *numComp = ncomp;
    
    /* Find anything */
    gotSome = gotSome || (ncomp>0);
    /* Cleanup */
    if (specCorr) g_free(specCorr); specCorr = NULL;
    if (fitSigma) g_free(fitSigma); fitSigma = NULL;
    if (fitParms) g_free(fitParms); fitParms = NULL;
  } /* End loop over model Stokes types */

  /* Point in->comps at Stokes I model */
  in->comps = ObitFArrayRef(in->VMIComps);
  return gotSome;
} /* end ObitSkyModelVMPolnLoadComps */

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSkyModelVMPolnInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModelVMPoln *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->isEVLA       = NULL;
  in->AntList      = NULL;
  in->curSource    = NULL;
  in->mosaicI    = in->mosaicQ    = in->mosaicU    = in->mosaicV    = NULL;
  in->VMIComps   = in->VMQComps   = in->VMUComps   = in->VMVComps   = NULL;
  in->specIndexI = in->specIndexQ = in->specIndexU = in->specIndexV = NULL;
  in->specFreqI  = in->specFreqQ  = in->specFreqU  = in->specFreqV  = NULL;
  in->startIComp = in->startQComp = in->startUComp = in->startVComp = NULL;
  in->endIComp   = in->endQComp   = in->endUComp   = in->endVComp   = NULL;
  in->RS   = in->RD   = in->LS   = in->LD  = NULL;
  in->RSc  = in->RDc  = in->LSc  = in->LDc = NULL;
  in->CX   = in->CY   = in->SX   = in->SY  = NULL;
  in->CXc  = in->CYc  = in->SYc  = in->SYc = NULL;
  in->numAntList   = 0;
  in->numAnt       = 0;
  in->PDrefAnt     = -999;
  in->numIComp = in->numQComp = in->numUComp = in->numVComp = 0;
  in->nSpecI   = in->nSpecQ   = in->nSpecU   = in->nSpecV   = 0;
} /* end ObitSkyModelVMPolnInit */


/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSkyModelVMPoln* cast to an Obit*.
 */
void ObitSkyModelVMPolnClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  VMPolnFTFuncArg *args;
  ObitSkyModelVMPoln *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->isEVLA)      g_free(in->isEVLA);      in->isEVLA      = NULL;
  if (in->specIndexI)  g_free(in->specIndexI);  in->specIndexI  = NULL;
  if (in->specIndexQ)  g_free(in->specIndexQ);  in->specIndexQ  = NULL;
  if (in->specIndexU)  g_free(in->specIndexU);  in->specIndexU  = NULL;
  if (in->specIndexV)  g_free(in->specIndexV);  in->specIndexV  = NULL;
  if (in->startIComp)  g_free(in->startIComp);  in->startIComp  = NULL;
  if (in->startQComp)  g_free(in->startQComp);  in->startQComp  = NULL;
  if (in->startUComp)  g_free(in->startUComp);  in->startUComp  = NULL;
  if (in->startVComp)  g_free(in->startVComp);  in->startVComp  = NULL;
  if (in->endIComp)    g_free(in->endIComp);    in->endIComp    = NULL;
  if (in->endQComp)    g_free(in->endQComp);    in->endQComp    = NULL;
  if (in->endUComp)    g_free(in->endUComp);    in->endUComp    = NULL;
  if (in->endVComp)    g_free(in->endVComp);    in->endVComp    = NULL;
  in->VMIComps  = ObitFArrayUnref(in->VMIComps);
  in->VMQComps  = ObitFArrayUnref(in->VMQComps);
  in->VMUComps  = ObitFArrayUnref(in->VMUComps);
  in->VMVComps  = ObitFArrayUnref(in->VMUComps);
  in->mosaicI   = ObitImageMosaicUnref(in->mosaicI);
  in->mosaicQ   = ObitImageMosaicUnref(in->mosaicQ);
  in->mosaicU   = ObitImageMosaicUnref(in->mosaicU);
  in->mosaicV   = ObitImageMosaicUnref(in->mosaicV);
  in->curSource = ObitSourceUnref(in->curSource);
  if (in->AntList)  {
    for (i=0; i<in->numAntList; i++) { 
      in->AntList[i] = ObitAntennaListUnref(in->AntList[i]);
    }
    g_free(in->AntList); in->AntList = NULL;
  }
  if (in->RS) {
    for (i=0; i<in->numUVChan; i++ ) {
      if (in->RS[i])  g_free(in->RS[i]);
      if (in->RD[i])  g_free(in->RD[i]);
      if (in->LS[i])  g_free(in->LS[i]);
      if (in->LD[i])  g_free(in->LD[i]);
      if (in->RSc[i]) g_free(in->RSc[i]);
      if (in->RDc[i]) g_free(in->RDc[i]);
      if (in->LSc[i]) g_free(in->LSc[i]);
      if (in->LDc[i]) g_free(in->LDc[i]);
    }
    g_free(in->RS);  in->RS  = NULL;
    g_free(in->RD);  in->RD  = NULL;
    g_free(in->LS);  in->LS  = NULL;
    g_free(in->LD);  in->LD  = NULL;
    g_free(in->RSc); in->RSc = NULL;
    g_free(in->RDc); in->RDc = NULL;
    g_free(in->LSc); in->LSc = NULL;
    g_free(in->LDc); in->LDc = NULL;
  }
  if (in->PPRL) g_free(in->PPRL); in->PPRL = NULL;
  if (in->PPLR) g_free(in->PPLR); in->PPLR = NULL;
  if (in->CX) {
    for (i=0; i<in->numUVChan; i++ ) {
      if (in->CX[i])  g_free(in->CX[i]);
      if (in->SX[i])  g_free(in->SX[i]);
      if (in->CY[i])  g_free(in->CY[i]);
      if (in->SY[i])  g_free(in->SY[i]);
      if (in->CXc[i]) g_free(in->CXc[i]);
      if (in->SXc[i]) g_free(in->SXc[i]);
      if (in->CYc[i]) g_free(in->CYc[i]);
      if (in->SYc[i]) g_free(in->SYc[i]);
    }
    g_free(in->CX);  in->CX  = NULL;
    g_free(in->SX);  in->SX  = NULL;
    g_free(in->CY);  in->CY  = NULL;
    g_free(in->SY);  in->SY  = NULL;
    g_free(in->CXc); in->CXc = NULL;
    g_free(in->SXc); in->SXc = NULL;
    g_free(in->CYc); in->CYc = NULL;
    g_free(in->SYc); in->SYc = NULL;
  }
  if (in->PPXY) g_free(in->PPXY); in->PPXY = NULL;
  if (in->PPYX) g_free(in->PPYX); in->PPYX = NULL;
    
  /* Thread stuff */
  if (in->threadArgs) {
    /* Check type - only handle "vmpoln" */
    args = (VMPolnFTFuncArg*)in->threadArgs[0];
    if ((strlen(args->type)>6) || (!strncmp(args->type, "vmpoln", 6))) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMPolnFTFuncArg*)in->threadArgs[i];
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
    } /* end if this a "vmpoln" threadArg */
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSkyModelVMPolnClear */

/**
 * Get input parameters from info member
 * If maxResid value not given or <0 then for it is determined
 * by examining each image in the Image mosaic.  
 * \param in Pointer to the ObitSkyModelVMPoln .
 * \param err Obit error stack object.
 */
void  ObitSkyModelVMPolnGetInput (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong *iptr, num;
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitSkyModelVMPolnGetInput";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelVMPolnIsA(in));
  if (!ObitInfoListIsA(in->info)) return;
  InfoReal.itg = 0;type = OBIT_oint;

  /* Call base class version */
  ObitSkyModelVMGetInput (inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);


  /* Make sure doPBCor turned off, will always use beam model instead */
  in->doPBCor = FALSE;
  
    /* Start CC number */
    if (ObitInfoListGetP(in->info, "BComp",  &type, (gint32*)dim, (gpointer)&iptr)) {
      num = MIN ( in->mosaic->numberImages, dim[0]);
      for (i=0; i<num; i++)   in->startIComp[i] = MAX (1,iptr[i]);
      if (in->startQComp)
	for (i=0; i<num; i++) in->startQComp[i] = MAX (1,iptr[i]);
      if (in->startUComp)
	for (i=0; i<num; i++) in->startUComp[i] = MAX (1,iptr[i]);
      if (in->startVComp)
	for (i=0; i<num; i++) in->startVComp[i] = MAX (1,iptr[i]);
    }

    /* End CC number */
    if (ObitInfoListGetP(in->info, "EComp",  &type, (gint32*)dim, (gpointer)&iptr)) {
      num = MIN ( in->mosaic->numberImages, dim[0]);
      for (i=0; i<num; i++)   in->endIComp[i] = iptr[i];
      if (in->endQComp)
	for (i=0; i<num; i++) in->endQComp[i] = iptr[i];
      if (in->endUComp)
	for (i=0; i<num; i++) in->endUComp[i] = iptr[i];
      if (in->endVComp)
	for (i=0; i<num; i++) in->endVComp[i] = iptr[i];
   }
 } /* end ObitSkyModelVMPolnGetInput */

/**
 * Decide which method is the most appropriate to calculate the FT of a model
 * Sets currentMode member function
 * Only DFT supported
 * \param in     Pointer to the ObitSkyModel .
 * \param uvdata UV data set
 */
void  ObitSkyModelVMPolnChose (ObitSkyModel* inn, ObitUV* uvdata) 
{
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  /* Only option implemented */
  in->currentMode = OBIT_SkyModel_DFT;

} /* end ObitSkyModelVMPolnChose */


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
void ObitSkyModelVMPolnFTDFT (ObitSkyModelVM *inn, olong field, ObitUV *uvdata, ObitErr *err)
{
  olong i, nvis, lovis, hivis, nvisPerThread, nThreads;
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)inn;
  VMPolnFTFuncArg *args;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelVMPolnFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* Check */
  args = (VMPolnFTFuncArg*)in->threadArgs[0];
  if ((strlen(args->type)>6) || (strncmp(args->type, "vmpoln", 6))) {
    Obit_log_error(err, OBIT_Error,"%s: Wrong type FuncArg %s", routine,args->type);
    return;
  }

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
    args = (VMPolnFTFuncArg*)in->threadArgs[i];
    args->in     = (ObitSkyModelVMPoln*)inn;
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
    if (i==(nThreads-1)) hivis = nvis;
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nThreads, in->DFTFunc, in->threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
}  /* end ObitSkyModelVMPolnFTDFT */

/**
 * Do Fourier transform using a DFT for a buffer of data.
 * Version for instrumental polarization.
 * All Stokes are calculated and modified.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * If doFlip member is true the Fourier transform is multiplied by sqrt(-1)
 * (for Stokes RL and LR)
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Method assumes same correction to all antennas.
 * Arguments are given in the VMPolnFTFuncArg structure passed as arg starting 
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
 * \li VM?Comps Stokes ? Component list
 * \li channel Current UV channel being processed (used in model update ).
 * \return NULL
 */
static gpointer ThreadSkyModelVMPolnFTDFT (gpointer args)
{
  /* Get arguments from structure */
  VMPolnFTFuncArg *largs = (VMPolnFTFuncArg*)args;
  ObitSkyModelVMPoln *in = (ObitSkyModelVMPoln*)largs->in;
  /*olong field      = largs->field;*/
  ObitUV *uvdata   = largs->uvdata;
  olong loVis      = largs->first-1;
  olong hiVis      = largs->last;
  olong ithread    = MAX (0, largs->ithread);
  ObitErr *err     = largs->err;

  olong iVis, iIF, iChannel, iStoke, iComp;
  olong lrec, nrparm, naxis[2], channel;
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong lstartChannel, lstartIF, lim, ifq;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF, iterm, nterm=0, nUVchan, nUVpoln;
  olong ilocu, ilocv, ilocw, iloct, ilocb, suba, itemp, ant1, ant2;
  olong lcompI, lcompQ, lcompU, lcompV, ncompI, ncompQ, ncompU, ncompV;
  olong mcomp, mcompI, mcompQ, mcompU, mcompV;
  ofloat *visData, *IData, *QData=NULL, *UData=NULL, *VData=NULL;
  ofloat *idata, *qdata, *udata, *vdata, *fscale, exparg;
  ofloat sumRealI,  sumImagI,  modRealRR=0.0,  modImagRR=0.0;
  ofloat sumRealQ,  sumImagQ,  modRealLL=0.0,  modImagLL=0.0;
  ofloat sumRealU,  sumImagU,  modRealRL=0.0,  modImagRL=0.0;
  ofloat sumRealV,  sumImagV,  modRealLR=0.0,  modImagLR=0.0;
  ofloat cbase, tx, ty, tz, ll, lll,  specFact;
  ofloat amp,  arg, freq2=0.0, freqFact, wtI=0.0, temp;
  ofloat logNuONu0I=0.0, logNuONu0Q=0.0, logNuONu0U=0.0, logNuONu0V=0.0; 
#define FazArrSize 100  /* Size of the amp/phase/sine/cosine arrays */
  ofloat AmpArrI[FazArrSize], AmpArrQ[FazArrSize];
  ofloat AmpArrU[FazArrSize], AmpArrV[FazArrSize];
  ofloat FazArrI[FazArrSize], FazArrQ[FazArrSize];
  ofloat FazArrU[FazArrSize], FazArrV[FazArrSize];
  ofloat CosArr[FazArrSize], SinArr[FazArrSize];
  ofloat cosSumPA, sinSumPA, cosDifPA, sinDifPA;
  ofloat model[8], stokes[8], u, v, w;
  olong it, jt, itcnt, qtcnt, utcnt, vtcnt, itab=0, qtab=0, utab=0, vtab=0;
  odouble *freqArr;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMPolnFTDFT";

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

  /* Get pointers for components */
  naxis[0] = 0; naxis[1] = 0; 
  IData = ObitFArrayIndex(in->VMIComps, naxis);
  lcompI = in->comps->naxis[0];   /* Length of row in comp table */
  ncompI = in->comps->naxis[1];   /* Number of components */
  mcompI = in->numIComp;           /* Actual number */
  mcomp  = mcompI;
  if (in->VMQComps) QData = ObitFArrayIndex(in->VMQComps, naxis);
  if (QData) {
    lcompQ = in->comps->naxis[0];   /* Length of row in comp table */
    ncompQ = in->comps->naxis[1];   /* Number of components */
    mcompQ = in->numQComp;           /* Actual number */
    mcomp  = MAX (mcomp, mcompQ);
  } else lcompQ = ncompQ = mcompQ = 0;
  if (in->VMUComps) UData = ObitFArrayIndex(in->VMUComps, naxis);
  if (UData) {
    lcompU = in->comps->naxis[0];   /* Length of row in comp table */
    ncompU = in->comps->naxis[1];   /* Number of components */
    mcompU = in->numUComp;           /* Actual number */
    mcomp  = MAX (mcomp, mcompU);
  } else lcompU = ncompU = mcompU = 0;
  if (in->VMVComps) VData = ObitFArrayIndex(in->VMVComps, naxis);
  if (VData) {
    lcompV = in->comps->naxis[0];   /* Length of row in comp table */
    ncompV = in->comps->naxis[1];   /* Number of components */
    mcompV = in->numVComp;           /* Actual number */
    mcomp  = MAX (mcomp, mcompV);
  } else lcompV = ncompV = mcompV = 0;

  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

  /* Current channel (0-rel) */
  channel = 0;

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
      
      /* Is current parallactic angles still valid? */
      if ((visData[iloct] > largs->endVMModelTime) || (visData[iloct] < largs->begVMModelTime)) {
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
      /* Baseline sum and difference of parallactic angles */
      cosSumPA = largs->cosPA[ant1]*largs->cosPA[ant2] - largs->sinPA[ant1]*largs->sinPA[ant2];
      sinSumPA = largs->sinPA[ant1]*largs->cosPA[ant2] + largs->cosPA[ant1]*largs->sinPA[ant2];
      cosDifPA = largs->cosPA[ant1]*largs->cosPA[ant2] + largs->sinPA[ant1]*largs->sinPA[ant2];
      sinDifPA = largs->sinPA[ant1]*largs->cosPA[ant2] - largs->cosPA[ant1]*largs->sinPA[ant2];
      
      /* Loop over IFs */
      channel = lstartIF* nUVchan + lstartChannel; /* UV Channel */
      for (iIF=lstartIF; iIF<startIF+numberIF; iIF++) {
	offsetIF = nrparm + iIF*jincif; 
	for (iChannel=lstartChannel; iChannel<startChannel+numberChannel; iChannel++) {
	  offsetChannel = offsetIF + iChannel*jincf; 
	  ifq = iIF*kincif + iChannel*kincf;
	  freqFact = fscale[ifq];  /* Frequency scaling factor */
	  freq2    = freqFact*freqFact;    /* Frequency factor squared */
	  
	  /* Log ratio of channel freq to Tabulated freq */
	  itab      = in->specIndexI[ifq];
	  logNuONu0I = (ofloat)log(freqArr[ifq]/in->specFreqI[itab]);
	  if (in->nModelStokes>=2) {
	    itab      = in->specIndexQ[ifq];
	    logNuONu0Q = (ofloat)log(freqArr[ifq]/in->specFreqQ[itab]);
	  }
	  if (in->nModelStokes>=3) {
	    itab      = in->specIndexU[ifq];
	    logNuONu0U = (ofloat)log(freqArr[ifq]/in->specFreqU[itab]);
	  }
	  if (in->nModelStokes>=4) {
	    itab      = in->specIndexV[ifq];
	    logNuONu0V = (ofloat)log(freqArr[ifq]/in->specFreqV[itab]);
	  } 

	  /* u,v,w at frequency */
	  u = (odouble)visData[ilocu]*freqFact;
	  v = (odouble)visData[ilocv]*freqFact;
	  w = (odouble)visData[ilocw]*freqFact;
	  
	  /* Sum over components */
	  /* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	  sumRealI = sumImagI = sumRealV = sumImagV = 0.0;
	  sumRealQ = sumImagQ = sumRealU = sumImagU = 0.0;
	  idata = IData; qdata = QData; udata = UData; vdata = VData;
	  
	  /* Sum by model type */
	  switch (in->modType) {
	  case OBIT_SkyModel_PointMod:     /* Point */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = qtcnt = utcnt = vtcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArrI[itcnt] = (idata[4]*u + idata[5]*v + idata[6]*w);
		/* Amplitude from component flux */
		AmpArrI[itcnt] = idata[3];  itcnt++;     
		if (iComp<mcompQ) {
		  FazArrQ[qtcnt] = (qdata[4]*u + qdata[5]*v + qdata[6]*w);
		  AmpArrQ[qtcnt] = qdata[3]; qtcnt++;}
		if (iComp<mcompU) {
		  FazArrU[utcnt] = (udata[4]*u + udata[5]*v + udata[6]*w);
		  AmpArrU[utcnt] = udata[3]; utcnt++;}
		if (iComp<mcompV) {
		  FazArrV[vtcnt] = (vdata[4]*u + vdata[5]*v + vdata[6]*w);
		  AmpArrV[vtcnt] = vdata[3]; vtcnt++;}
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_PointModTSpec:     /* Point + tabulated spectrum */
	    itab = 8 + in->specIndexI[ifq];
	    if (in->specIndexQ) qtab = 8 + in->specIndexQ[ifq];
	    if (in->specIndexU) utab = 8 + in->specIndexU[ifq];
	    if (in->specIndexV) vtab = 8 + in->specIndexV[ifq];
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = qtcnt = utcnt = vtcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		if (idata[itab]!=0.0) {  /* valid? */
		  tx = idata[4]*(odouble)u;
		  ty = idata[5]*(odouble)v;
		  tz = idata[6]*(odouble)w;
		  FazArrI[itcnt] = (tx + ty + tz);
		  specFact = exp(-logNuONu0I * idata[7]);
		  AmpArrI[itcnt] = specFact * idata[itab];  itcnt++;
		} /* end I OK */
		if ((iComp<mcompQ) && (qdata[qtab]!=0.0)) {  /* valid? */
		  tx = qdata[4]*(odouble)u;
		  ty = qdata[5]*(odouble)v;
		  tz = qdata[6]*(odouble)w;
		  FazArrQ[qtcnt] = (tx + ty + tz);
		  specFact = exp(-logNuONu0Q * qdata[7]);
		  AmpArrQ[qtcnt] = specFact * qdata[qtab];  qtcnt++;
		} /* end Q OK */
		if ((iComp<mcompU) && (udata[utab]!=0.0)) {  /* valid? */
		  tx = udata[4]*(odouble)u;
		  ty = udata[5]*(odouble)v;
		  tz = udata[6]*(odouble)w;
		  FazArrU[utcnt] = (tx + ty + tz);
		  specFact = exp(-logNuONu0U * udata[7]);
		  AmpArrU[utcnt] = specFact * udata[utab];  utcnt++;
		} /* end U OK */
		if ((iComp<mcompV) && (vdata[vtab]!=0.0)) {  /* valid? */
		  tx = vdata[4]*(odouble)u;
		  ty = vdata[5]*(odouble)v;
		  tz = vdata[6]*(odouble)w;
		  FazArrV[vtcnt] = (tx + ty + tz);
		  specFact = exp(-logNuONu0V * vdata[7]);
		  AmpArrV[vtcnt] = specFact * vdata[vtab];  vtcnt++;
		} /* end V OK */
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_PointModSpec:     /* Point + spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		if (idata[3]!=0.0) {  /* valid? */
		  tx = idata[4]*(odouble)u;
		  ty = idata[5]*(odouble)v;
		  tz = idata[6]*(odouble)w;
		  FazArrI[itcnt] = (tx + ty + tz);
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += idata[7+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  AmpArrI[itcnt] = specFact * idata[3];  itcnt++;
		} /* end I OK */
		if ((iComp<mcompQ) && (qdata[3]!=0.0)) {  /* valid? */
		  tx = qdata[4]*(odouble)u;
		  ty = qdata[5]*(odouble)v;
		  tz = qdata[6]*(odouble)w;
		  FazArrQ[qtcnt] = (tx + ty + tz);
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += qdata[7+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  AmpArrQ[qtcnt] = specFact * qdata[3];  qtcnt++;
		} /* end Q OK */
		if ((iComp<mcompU) && (udata[3]!=0.0)) {  /* valid? */
		  tx = udata[4]*(odouble)u;
		  ty = udata[5]*(odouble)v;
		  tz = udata[6]*(odouble)w;
		  FazArrU[utcnt] =(tx + ty + tz);
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += udata[7+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  AmpArrU[utcnt] = specFact * udata[3];  utcnt++;
		} /* end U OK */
		if ((iComp<mcompV) && (vdata[3]!=0.0)) {  /* valid? */
		  tx = vdata[4]*(odouble)u;
		  ty = vdata[5]*(odouble)v;
		  tz = vdata[6]*(odouble)w;
		  FazArrV[vtcnt] = (tx + ty + tz);
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += vdata[7+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  AmpArrV[vtcnt] = specFact * vdata[3];  vtcnt++;
		} /* end V OK */
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = qtcnt = utcnt = vtcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArrI[itcnt] = (idata[4]*u + idata[5]*v + idata[6]*w);
		/* Amplitude from component flux, Gaussian */
		arg = (idata[7]*u*u + idata[8]*v*v + idata[9]*u*v);
		if (arg<-1.0e-5) exparg = exp (arg);
		else exparg = 1.0;
		AmpArrI[itcnt] = exparg*idata[3];  itcnt++;   /* End I */  
		if (iComp<mcompQ) {
		  FazArrQ[qtcnt] = (qdata[4]*u + qdata[5]*v + qdata[6]*w);
		  arg = freq2 * (qdata[7]*u*u + qdata[8]*v*v + qdata[9]*u*v);
		  if (arg<-1.0e-5) exparg = exp (arg);
		  else exparg = 1.0;
		  AmpArrQ[qtcnt] = exparg*qdata[3]; qtcnt++;} /* End Q */
		if (iComp<mcompU) {
		  FazArrU[utcnt] = (udata[4]*u + udata[5]*v + udata[6]*w);
		  arg = (udata[7]*u*u + udata[8]*v*v + udata[9]*u*v);
		  if (arg<-1.0e-5) exparg = exp (arg);
		  else exparg = 1.0;
		  AmpArrU[utcnt] = exparg*udata[3]; utcnt++;} /* End U */
		if (iComp<mcompV) {
		  FazArrV[vtcnt] = (vdata[4]*u + vdata[5]*v + vdata[6]*w);
		  arg = (vdata[7]*u*u + vdata[8]*v*v + vdata[9]*u*v);
		  if (arg<-1.0e-5) exparg = exp (arg);
		  else exparg = 1.0;
		  AmpArrV[vtcnt] = exparg*vdata[3]; vtcnt++;} /* End V */
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModTSpec:     /* Gaussian on sky + tabulated spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = qtcnt = utcnt = vtcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      ll = log(freqFact); 
	      itab = 11 + in->specIndexI[ifq];
	      if (in->specIndexQ) qtab = 11 + in->specIndexQ[ifq];
	      if (in->specIndexU) utab = 11 + in->specIndexU[ifq];
	      if (in->specIndexV) vtab = 11 + in->specIndexV[ifq];
	      for (iComp=it; iComp<lim; iComp++) {
		if (idata[3]!=0.0) {  /* valid? */
		  FazArrI[itcnt] = (idata[4]*u + idata[5]*v + idata[6]*w);
		  specFact = exp(-logNuONu0I * idata[10]);
		  arg = (idata[7]*u*u + idata[8]*v*v + idata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrI[itcnt] = amp*idata[3];  itcnt++;
		} /* End I OK */
		if ((iComp<mcompQ) && (qdata[3]!=0.0)) {
		  FazArrQ[qtcnt] = (qdata[4]*u + qdata[5]*v + qdata[6]*w);
		  specFact = exp(-logNuONu0Q * qdata[10]);
		  arg = (qdata[7]*u*u + qdata[8]*v*v + qdata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrQ[qtcnt] = amp*qdata[3]; qtcnt++;}  /* End Q */
	      
		if ((iComp<mcompU) && (udata[3]!=0.0)) {
		  FazArrU[utcnt] = (udata[4]*u + udata[5]*v + udata[6]*w);
		  specFact = exp(-logNuONu0U * udata[10]);
		  arg = (udata[7]*u*u + udata[8]*v*v + udata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrU[utcnt] = amp*udata[3]; utcnt++;} /* End U */
		if ((iComp<mcompV) && (vdata[3]!=0.0)) {
		  FazArrV[vtcnt] = (vdata[4]*u + vdata[5]*v + vdata[6]*w);
		  specFact = exp(-logNuONu0V * vdata[10]);
		  arg = (vdata[7]*u*u + vdata[8]*v*v + vdata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrV[vtcnt] = amp*vdata[3]; vtcnt++;}  /* End V */
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_GaussModSpec:     /* Gaussian on sky + spectrum */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = qtcnt = utcnt = vtcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      ll = log(freqFact); 
	      for (iComp=it; iComp<lim; iComp++) {
		if (idata[3]!=0.0) {  /* valid? */
		  FazArrI[itcnt] = (idata[4]*u + idata[5]*v + idata[6]*w);
		  /* Frequency dependent term */
		  lll = ll;  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += idata[10+iterm] * lll;  lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = (idata[7]*u*u + idata[8]*v*v + idata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrI[itcnt] = amp*idata[3];  itcnt++;
		} /* End I OK */
		if ((iComp<mcompQ) && (qdata[3]!=0.0)) {
		  FazArrQ[qtcnt] = (qdata[4]*u + qdata[5]*v + qdata[6]*w);
		  lll = ll; arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += qdata[10+iterm] * lll;  lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = (qdata[7]*u*u + qdata[8]*v*v + qdata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrQ[qtcnt] = amp*qdata[3]; qtcnt++;}  /* End Q */
	      
		if ((iComp<mcompU) && (udata[3]!=0.0)) {
		  FazArrU[utcnt] = (udata[4]*u + udata[5]*v + udata[6]*w);
		  lll = ll;  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += udata[10+iterm] * lll;  lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = (udata[7]*u*u + udata[8]*v*v + udata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrU[utcnt] = amp*udata[3]; utcnt++;} /* End U */
		if ((iComp<mcompV) && (vdata[3]!=0.0)) {
		  FazArrV[vtcnt] = (vdata[4]*u + vdata[5]*v + vdata[6]*w);
		  lll = ll; arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += vdata[10+iterm] * lll; lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = (vdata[7]*u*u + vdata[8]*v*v + vdata[9]*u*v);
		  if (arg<-1.0e-5) amp = specFact * exp (arg);
		  else amp = specFact;
		  AmpArrV[vtcnt] = amp*vdata[3]; vtcnt++;}  /* End V */
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	    /* From the AIPSish QSPSUB.FOR  */
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = qtcnt = utcnt = vtcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      for (iComp=it; iComp<lim; iComp++) {
		FazArrI[itcnt] = (idata[4]*u + idata[5]*v + idata[6]*w);
		/* Amplitude from component flux, uniform sphere, spectrum */
		if (idata[3]!=0.0) {  /* valid? */
		  /* Frequency dependent term */
		  lll = ll = log(freqFact);
		  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += idata[8+iterm] * lll;
		    lll *= ll;
		  }
		  specFact = exp(arg);
		  arg =  sqrt(u*u + v*v) * idata[7];
		  arg = MAX (arg, 0.1);
		  arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  AmpArrI[itcnt] = arg*idata[3];  itcnt++;     
		} /* End I OK */
		  if (iComp<mcompQ) {
		    FazArrQ[qtcnt] = (qdata[4]*u + qdata[5]*v + qdata[6]*w);
		    arg = sqrt(u*u + v*v) * qdata[7];
		    arg = MAX (arg, 0.1);
		    arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		    AmpArrQ[qtcnt] = arg*qdata[3]; qtcnt++;}  /* End Q */
		if (iComp<mcompU) {
		  FazArrU[utcnt] = (udata[4]*u + udata[5]*v + udata[6]*w);
		  arg = sqrt(u*u + v*v) * udata[7];
		  arg = MAX (arg, 0.1);
		  arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  AmpArrU[utcnt] = arg*udata[3]; utcnt++;} /* End U */
		if (iComp<mcompV) {
		  FazArrV[vtcnt] = (vdata[4]*u + vdata[5]*v + vdata[6]*w);
		  arg = sqrt(u*u + v*v) * vdata[7];
		  arg = MAX (arg, 0.1);
		  arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  AmpArrV[vtcnt] = arg*vdata[3]; vtcnt++;} /* End V */
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  case OBIT_SkyModel_USphereModSpec:    /* Uniform sphere + spectrum*/
	    for (it=0; it<mcomp; it+=FazArrSize) {
	      itcnt = qtcnt = utcnt = vtcnt = 0;
	      lim = MIN (mcomp, it+FazArrSize);
	      ll = log(freqFact); 
	      for (iComp=it; iComp<lim; iComp++) {
		if (idata[3]!=0.0) {  /* valid? */
		  FazArrI[itcnt] = (idata[4]*u + idata[5]*v + idata[6]*w);
		  /* Frequency dependent term */
		  lll = ll;  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += idata[10+iterm] * lll;  lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = sqrt(u*u + v*v) * idata[9];
		  arg = MAX (arg, 0.1);
		  arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  AmpArrI[itcnt] = specFact * arg * idata[3];  itcnt++;     
		} /* End I OK */
		if ((iComp<mcompQ) && (qdata[3]!=0.0)) {
		  FazArrQ[qtcnt] = (qdata[4]*u + qdata[5]*v + qdata[6]*w);
		  lll = ll; arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += qdata[10+iterm] * lll;  lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = sqrt(u*u + v*v) * qdata[9];
		  arg = MAX (arg, 0.1);
		  arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  AmpArrQ[qtcnt] = specFact * arg * qdata[3];  qtcnt++;}  /* End Q */     
		if ((iComp<mcompU) && (udata[3]!=0.0)) {
		  FazArrU[utcnt] = (udata[4]*u + udata[5]*v + udata[6]*w);
		  lll = ll;  arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += udata[10+iterm] * lll;  lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = sqrt(u*u + v*v) * udata[9];
		  arg = MAX (arg, 0.1);
		  arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  AmpArrU[utcnt] = specFact * arg * udata[3];  utcnt++;}  /* End U */     
		if ((iComp<mcompV) && (vdata[3]!=0.0)) {
		  FazArrV[vtcnt] = (vdata[4]*u + vdata[5]*v + vdata[6]*w);
		  lll = ll; arg = 0.0;
		  for (iterm=0; iterm<nterm; iterm++) {
		    arg += udata[10+iterm] * lll; lll *= ll;
		  }
		  specFact = exp(arg);
		  arg = sqrt(u*u + v*v) * vdata[9];
		  arg = MAX (arg, 0.1);
		  arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
		  AmpArrV[vtcnt] = specFact * arg * vdata[3];  vtcnt++;}  /* End V */     
		idata += lcompI;  qdata += lcompQ; /* update indices */
		udata += lcompU;  vdata += lcompV;
	      } /* end inner loop over components */
	      
	      /* Convert phases to sin/cos */
	      ObitSinCosVec(itcnt, FazArrI, SinArr, CosArr);
	      /* Accumulate real and imaginary parts */
	      for (jt=0; jt<itcnt; jt++) {
		sumRealI += AmpArrI[jt]*CosArr[jt]; sumImagI += AmpArrI[jt]*SinArr[jt];
	      } 
	      /* Q */
	      ObitSinCosVec(qtcnt, FazArrQ, SinArr, CosArr);
	      for (jt=0; jt<qtcnt; jt++) {
		sumRealQ += AmpArrQ[jt]*CosArr[jt]; sumImagQ += AmpArrQ[jt]*SinArr[jt];
	      } 
	      /* U */
	      ObitSinCosVec(utcnt, FazArrU, SinArr, CosArr);
	      for (jt=0; jt<utcnt; jt++) {
		sumRealU += AmpArrU[jt]*CosArr[jt]; sumImagU += AmpArrU[jt]*SinArr[jt];
	      } 
	      /* V */
	      ObitSinCosVec(vtcnt, FazArrV, SinArr, CosArr);
	      for (jt=0; jt<vtcnt; jt++) {
		sumRealV += AmpArrV[jt]*CosArr[jt]; sumImagV += AmpArrV[jt]*SinArr[jt];
	      } 
	    } /* end outer loop over components */
	    break;
	  default:
	    ObitThreadLock(in->thread);  /* Lock against other threads */
	    Obit_log_error(err, OBIT_Error,"%s Unknown Comp model type %d in %s",
			   routine, in->modType, in->name);
	    ObitThreadUnlock(in->thread); 
	    goto finish;
	  }; /* end switch by model type */
	  
	  /* Get corrected values for correlation by type */
	  if (in->isCirc) { /* Circular feeds - prior parallactic angle correction assumed */
	    if (in->polType==OBIT_UVPoln_ELORI) {
	      stokes[0] = sumRealI - sumRealV; stokes[1] =  sumImagI - sumImagV;
	      stokes[2] = sumRealI + sumRealV; stokes[3] =  sumImagI + sumImagV;
	      stokes[4] = sumImagQ + sumRealU; stokes[5] = -sumRealQ + sumImagU;
	      stokes[6] = sumImagQ - sumRealU; stokes[7] = -sumRealQ - sumImagU;
	      calcModel(largs, ant1, ant2, ifq, stokes, model);
	      modRealRR = model[0];  modImagRR = model[1];
	      modRealLL = model[2];  modImagLL = model[3];
	      modRealRL = model[4];  modImagRL = model[5];
	      modRealLR = model[6];  modImagLR = model[7];
	    } else {  /* No correction */
	      modRealRR =  sumRealI - sumRealV;
	      modImagRR =  sumImagI - sumImagV;
	      modRealLL =  sumRealI + sumRealV;
	      modImagLL =  sumImagI + sumImagV;
	      modRealRL =  sumImagQ + sumRealU;
	      modImagRL = -sumRealQ + sumImagU;
	      modRealLR =  sumImagQ - sumRealU;
	      modImagLR = -sumRealQ - sumImagU;
	    }
	  } else { /* Linear Feeds (really XX, YY, XY, YX ) 
		    no prior parallactic angle correction */
	    if (in->polType==OBIT_UVPoln_ELORI) {
	      stokes[0] = sumRealI - sumRealV; stokes[1] =  sumImagI - sumImagV;
	      stokes[2] = sumRealI + sumRealV; stokes[3] =  sumImagI + sumImagV;
	      stokes[4] = sumImagQ + sumRealU; stokes[5] = -sumRealQ + sumImagU;
	      stokes[6] = sumImagQ - sumRealU; stokes[7] = -sumRealQ - sumImagU;
	      calcModel(largs, ant1, ant2, ifq, stokes, model);
	      modRealRR = model[0];  modImagRR = model[1];
	      modRealLL = model[2];  modImagLL = model[3];
	      modRealRL = model[4];  modImagRL = model[5];
	      modRealLR = model[6];  modImagLR = model[7];
	    } else {  /* No correction (really XX, YY, XY, YX ) */
	      modRealRR =  sumRealI + sumRealQ*cosSumPA + sumRealU*sinSumPA;
	      modImagRR =  sumImagI + sumImagQ*cosSumPA + sumImagU*sinSumPA;
	      modRealLL =  sumRealI - sumRealQ*cosSumPA - sumRealU*sinSumPA;
	      modImagLL =  sumImagI - sumImagQ*cosSumPA - sumImagU*sinSumPA;
	      modRealRL = -sumRealQ*sinSumPA + sumRealU*cosSumPA - sumImagV;
	      modImagRL = -sumImagQ*sinSumPA + sumImagU*cosSumPA + sumRealV;
	      modRealLR = -sumRealQ*sinSumPA + sumRealU*cosSumPA + sumImagV;
	      modImagLR = -sumImagQ*sinSumPA + sumImagU*cosSumPA - sumRealV;
	    }
	  }

	  /* Replacing data? */
	  if (in->doReplace) {
	    iStoke = 0;   /* RR (XX) */
	    offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	    visData[offset]   = modRealRR;
	    visData[offset+1] = modImagRR;
	    offset += jincs;   /* LL (YY) */
	    visData[offset]   = modRealLL;
	    visData[offset+1] = modImagLL;
	    offset += jincs;   /* RL (XY) */
	    visData[offset]   = modRealRL;
	    visData[offset+1] = modImagRL;
	    offset += jincs;   /* LR (YX) */
	    visData[offset]   = modRealLR;
	    visData[offset+1] = modImagLR;
	    goto doneChan;
	  } /* end relpace */
	  
	  /* Dividing? - also correct weight */
	  if (in->doDivide) {
	    iStoke = 0;   /* RR (XX) */
	    wtI = modRealRR * modRealRR + modImagRR * modImagRR;
	    offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	    temp = modRealRR * visData[offset] + modImagRR * visData[offset+1];
	    visData[offset+1] = (modRealRR * visData[offset+1] - modImagRR * visData[offset])/wtI;
	    visData[offset]   = temp/wtI;
	    visData[offset+2] *= sqrt(wtI);  /* correct weight */

	    offset += jincs;   /* LL (YY) */
	    wtI = modRealLL * modRealLL + modImagLL * modImagLL;
	    temp = modRealLL * visData[offset] + modImagLL * visData[offset+1];
	    visData[offset+1] = (modRealLL * visData[offset+1] - modImagLL * visData[offset])/wtI;
	    visData[offset]   = temp/wtI;
	    visData[offset+2] *= sqrt(wtI);  /* correct weight */

	    offset += jincs;   /* RL (XY) */
	    wtI = modRealRL * modRealRL + modImagRL * modImagRL;
	    temp = modRealRL * visData[offset] + modImagRL * visData[offset+1];
	    visData[offset+1] = (modRealRL * visData[offset+1] - modImagRL * visData[offset])/wtI;
	    visData[offset]   = temp/wtI;
	    visData[offset+2] *= sqrt(wtI);  /* correct weight */

	    offset += jincs;   /* LR (YX) */
	    wtI = modRealLR * modRealLR + modImagLR * modImagLR;
	    temp = modRealLR * visData[offset] + modImagLR * visData[offset+1];
	    visData[offset+1] = (modRealLR * visData[offset+1] - modImagLR * visData[offset])/wtI;
	    visData[offset]   = temp/wtI;
	    visData[offset+2] *= sqrt(wtI);  /* correct weight */

	    goto doneChan;
	  }
	  
	  /* Subtract */
	  iStoke = 0;   /* RR (XX) */
	  offset = offsetChannel + iStoke*jincs; /* Visibility offset */
	  visData[offset]   -= modRealRR;
	  visData[offset+1] -= modImagRR;
	  
	  offset += jincs;   /* LL (YY) */
	  visData[offset]   -= modRealLL;
	  visData[offset+1] -= modImagLL;
	  
	  offset += jincs;   /* RL (XY) */
	  visData[offset]   -= modRealRL;
	  visData[offset+1] -= modImagRL;
	  
	  offset += jincs;   /* LR (YX) */
	  visData[offset]   -= modRealLR;
	  visData[offset+1] -= modImagLR;
	  	  
	doneChan:
	  offsetChannel += jincf;
	  channel++; /* Finished another channel */
	} /* end loop over Channel */
	offsetIF += jincif;
      } /* end loop over IF */
      
      visData += lrec; /* Update vis pointer */
    } /* end loop over visibilities */
  } /* end outer frequency loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadSkyModelVMPolnFTDFT */

/** 
 *  Get MF frequency information
 * \param in        Sky model object
 * \param image     Image with MF info
 * \param uvdata    UVdata with freq chann/IFs
 * \param Alpha     [out] prior spectral index applied
 * \param AlphaRefF [out] prior spectral index ref freq (Hz)
 * \param specFreq  [out] Array of MF frequencies, g_free when done
 * \param specIndex [out] Array of MF freq per uvdata chan/IF, g_free
 * \param forceDFT  [out] TRUE if DFT needed due to frequency incompatabilities
 * \return number of MF frequency bins, 1=> no MF info.
*/
static olong PolnMFFreqInfo (ObitSkyModelVMPoln *in, ObitImage *image, ObitUV *uvdata, ofloat *Alpha,  
			     odouble *AlphaRefF, odouble **specFreq, olong **specIndex,
			     gboolean *forceDFT)
{
  olong nSpec=1;
  olong i, j, nfreq, nif, n;
  gboolean lsb;
  odouble refFreq, *specFreqLo=NULL, *specFreqHi=NULL;
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar keyword[12];

  *forceDFT = FALSE;   /* In case */
  ObitInfoListGetTest(image->myDesc->info, "NSPEC", &type, dim, &nSpec);
  if (nSpec<=1) return nSpec;

  refFreq = image->myDesc->crval[image->myDesc->jlocf];
  /* get number of and channel frequencies for CC spectra from 
     CC table on first image in mosaic */
  *specFreq  = g_malloc0(nSpec*sizeof(odouble));
  specFreqLo = g_malloc0(nSpec*sizeof(odouble));
  specFreqHi = g_malloc0(nSpec*sizeof(odouble));
  for (i=0; i<nSpec; i++) {
    (*specFreq)[i] = 1.0;
    sprintf (keyword, "FREQ%4.4d",i+1);
    ObitInfoListGetTest(image->myDesc->info, keyword, &type, dim, &((*specFreq)[i]));
    sprintf (keyword, "FREL%4.4d",i+1);
    ObitInfoListGetTest(image->myDesc->info, keyword, &type, dim, &specFreqLo[i]);
    sprintf (keyword, "FREH%4.4d",i+1);
    ObitInfoListGetTest(image->myDesc->info, keyword, &type, dim, &specFreqHi[i]);
  }

  /* Make sure overlap */
  for (i=1; i<nSpec-1; i++) {
    if (specFreqHi[i]>specFreqLo[i-1]) {
      specFreqHi[i]   = 0.5*(specFreqHi[i]+specFreqLo[i-1]);
      specFreqLo[i-1] = specFreqHi[i];
    }
  }

  /* Prior spectral index */
  InfoReal.flt = 0.0;   type = OBIT_float;
  ObitInfoListGetTest(image->myDesc->info, "ALPHA", &type, dim, &InfoReal);
  if (type==OBIT_double) *Alpha = (ofloat)InfoReal.dbl;
  if (type==OBIT_float)  *Alpha = (ofloat)InfoReal.flt;

  *AlphaRefF = refFreq;
  ObitInfoListGetTest(image->myDesc->info, "RFALPHA", &type, dim, AlphaRefF);
  
  /* Make array of which coarse spectrum value corresponds to each uv channel */
  nfreq = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  if (uvdata->myDesc->jlocif>=0) 
    nif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  else nif = 1;
  lsb = uvdata->myDesc->cdelt[uvdata->myDesc->jlocf]<0.0; /* Lower sideband? */
  n = nfreq*nif;
  *specIndex = g_malloc0(n*sizeof(olong)); 
  *forceDFT = FALSE;   /* Need to force DFT? */
  for (i=0; i<n; i++) {
    (*specIndex)[i] = -1;
    for (j=0; j<nSpec; j++) {
      if (lsb) { /* Lower sideband */
	if ((uvdata->myDesc->freqArr[i] <= specFreqLo[j]) &&
	    (uvdata->myDesc->freqArr[i] >= specFreqHi[j]))
	  (*specIndex)[i] = j;
      } else {  /* Upper sideband */
	if ((uvdata->myDesc->freqArr[i] >= specFreqLo[j]) &&
	    (uvdata->myDesc->freqArr[i] <= specFreqHi[j]))
	  (*specIndex)[i] = j;
      }
    }
    /* If frequency out of range, force DFT */
    if ((*specIndex)[i]<0) {
      *forceDFT = TRUE;
      /* Set to closest */
      if ( lsb && (uvdata->myDesc->freqArr[i]<specFreqLo[nSpec-1]))      (*specIndex)[i] = nSpec-1;
      else if ( lsb && (uvdata->myDesc->freqArr[i]>specFreqHi[0]))       (*specIndex)[i] = 0;
      else if (!lsb && (uvdata->myDesc->freqArr[i]<specFreqLo[0]))       (*specIndex)[i] = 0;
      else if (!lsb && (uvdata->myDesc->freqArr[i]>specFreqHi[nSpec-1])) (*specIndex)[i] = nSpec-1;
      else (*specIndex)[i] = 0;  /* Shouldn't happen */
    }
  } /* End of loop making lookup table */
  if (specFreqLo) g_free(specFreqLo);
  if (specFreqHi) g_free(specFreqHi);
  return nSpec;
} /* end PolnMFFreqInfo */

/**
 * Digest PD Table
 * \param in     SkyModelVMPoln
 * \param err    Obit error/message object
 */
static void digestPDTable (ObitSkyModelVMPoln *in, ObitUV *uvdata, olong PDVer, ObitErr *err)
{
  olong i, iRow, numPol=0, numIF=0, numChan=0;
  olong kndx, jndx, indx, iif, ich, ia;
  ObitIOCode retCode;
  ObitTablePD *PDTab;
  ObitTablePDRow *PDRow=NULL;
  ObitUVDesc *uvDesc=uvdata->myDesc;
  dcomplex ct1, ct2, Jp, Jm;
  /* dcomplex PRref, PLref;*/
  gboolean haveRLPhase;
  ofloat PD, root2;
  gchar *routine = "digestPDTable";

  /* Complex constants */
  COMPLEX_SET (Jp,  0.0, 1.0);
  COMPLEX_SET (Jm,  0.0,-1.0);

  /* Get table */
  PDTab = newObitTablePDValue ("Instrum", (ObitData*)uvdata, &PDVer, 
			       OBIT_IO_ReadOnly, numPol, numIF, numChan, 
			       err);
  /* Open input table */
  retCode = ObitTablePDOpen (PDTab, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, PDTab->name);
  /* Set row */
  PDRow  = newObitTablePDRow (PDTab);

  /* Sizes of things */
  in->numAnt = PDTab->numAnt;

  /* Check compatability with UV data */
  Obit_return_if_fail(((uvDesc->jlocif<0) || (PDTab->numIF==uvDesc->inaxes[uvDesc->jlocif])), err,
		      "%s: UV, PD IFs inconsistent %d %d", 
		      routine, PDTab->numIF, uvDesc->inaxes[uvDesc->jlocif]);
  Obit_return_if_fail((PDTab->numChan==uvDesc->inaxes[uvDesc->jlocf]), err,
		      "%s: UV, PD Channels inconsistent %d %d", 
		      routine, PDTab->numChan, uvDesc->inaxes[uvDesc->jlocf]);
  in->numUVChan = PDTab->numIF*PDTab->numChan;

  /* Circular or linear feeds? */
  in->isCirc = (uvdata->myDesc->crval[uvdata->myDesc->jlocs]<0.0) &&
                (uvdata->myDesc->crval[uvdata->myDesc->jlocs]>-4.0);

  /* Create arrays */
  if (in->isCirc) {   /* Circular */
    in->PPRL = g_malloc0(in->numUVChan*sizeof(dcomplex));
    in->PPLR = g_malloc0(in->numUVChan*sizeof(dcomplex));
    in->RS   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->RD   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->LS   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->LD   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->RSc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->RDc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->LSc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->LDc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    for (i=0; i<in->numUVChan; i++) {
      in->RS[i]   = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->RD[i]   = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->LS[i]   = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->LD[i]   = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->RSc[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->RDc[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->LSc[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->LDc[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
    }
  } else { /* Linear */
    in->PPXY = g_malloc0(in->numUVChan*sizeof(dcomplex));
    in->PPYX = g_malloc0(in->numUVChan*sizeof(dcomplex));
    in->SYc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->CX   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->SX   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->CY   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->SY   = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->CXc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->SXc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    in->CYc  = g_malloc0(in->numUVChan*sizeof(dcomplex*));
    for (i=0; i<in->numUVChan; i++) {
      in->CX[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->SX[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->CY[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->SY[i]  = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->CXc[i] = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->SXc[i] = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->CYc[i] = g_malloc0(in->numAnt*sizeof(dcomplex));
      in->SYc[i] = g_malloc0(in->numAnt*sizeof(dcomplex));
    }
  }
  /* Parameterization type */
  in->polType = ObitAntennaListGetPolType (PDTab->polType);

  /* Loop over table */
  haveRLPhase = PDTab->myDesc->repeat[PDTab->RLPhaseCol]>0;
  for (iRow=1; iRow<=PDTab->myDesc->nrow; iRow++) {
    retCode = ObitTablePDReadRow (PDTab, iRow, PDRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, PDTab->name);
    if (PDRow->status==-1) continue;

    if (in->PDrefAnt<0) in->PDrefAnt = PDRow->RefAnt;

    ia = PDRow->antNo - 1;
    kndx = jndx = indx = 0;
    root2 = 1.0 / sqrt(2.0);
    for (iif=0; iif<PDTab->numIF; iif++) {
      for (ich=0; ich<PDTab->numChan; ich++) {
	/* See if RL (XY) Phase difference given */
	if (haveRLPhase) 
	  PD = PDRow->RLPhase[kndx];
	else 
	  PD = 0.0;
	/* Digest by type */
	if (in->isCirc) {   /* Circular feeds */
	  COMPLEX_SET(in->RS[indx][ia], root2*(cos(PDRow->Real1[kndx]) + sin(PDRow->Real1[kndx])), 0);
	  COMPLEX_SET(ct1,               root2*(cos(PDRow->Real1[kndx]) - sin(PDRow->Real1[kndx])), 0);
	  COMPLEX_EXP(ct2, 2*PDRow->Imag1[kndx]);
	  COMPLEX_MUL2 (in->RD[indx][ia], ct1, ct2);
	  COMPLEX_SET (ct1, root2*(cos(PDRow->Real2[kndx]) + sin(PDRow->Real2[kndx])), 0);
	  COMPLEX_EXP (ct2, -2*PDRow->Imag2[kndx]);
	  COMPLEX_MUL2 (in->LS[indx][ia], ct1, ct2);
	  COMPLEX_SET (in->LD[indx][ia], root2*(cos(PDRow->Real2[kndx]) - sin(PDRow->Real2[kndx])), 0);
	  COMPLEX_CONJUGATE (in->RSc[indx][ia], in->RS[indx][ia]);
	  COMPLEX_CONJUGATE (in->RDc[indx][ia], in->RD[indx][ia]);
	  COMPLEX_CONJUGATE (in->LSc[indx][ia], in->LS[indx][ia]);
	  COMPLEX_CONJUGATE (in->LDc[indx][ia], in->LD[indx][ia]);
	  if (ia==in->PDrefAnt) {
	    /* This should already have been done to the data
	       COMPLEX_EXP (PRref,  PDRow->Imag1[kndx]);
	       COMPLEX_EXP (PLref, -PDRow->Imag2[kndx]+PD);
	       COMPLEX_CONJUGATE (ct1, PLref);
	       COMPLEX_MUL2 (in->PPRL[indx], PRref, ct1);
	       COMPLEX_CONJUGATE (ct1, PRref);
	       COMPLEX_MUL2 (in->PPLR[indx], PLref, ct1);
	    */
	    COMPLEX_SET(in->PPRL[indx], 1.0, 0.0);
	    COMPLEX_SET(in->PPLR[indx], 1.0, 0.0);
	  }
	} else {             /* Linear feeds */

	  COMPLEX_EXP (ct1, -PDRow->Imag1[kndx]);
	  COMPLEX_SET (ct2, cos(G_PI*0.25+PDRow->Real1[kndx]), 0.0);
	  COMPLEX_MUL3 (in->CX[indx][ia], Jp, ct1, ct2);
	  COMPLEX_EXP (ct1, PDRow->Imag1[kndx]);
	  COMPLEX_SET (ct2, sin(G_PI*0.25+PDRow->Real1[kndx]), 0.0);
	  COMPLEX_MUL2 (in->SX[indx][ia], ct1, ct2);
	  COMPLEX_EXP (ct1, PDRow->Imag2[kndx]);
	  COMPLEX_SET (ct2, cos(G_PI*0.25-PDRow->Real2[kndx]), 0.0);
	  COMPLEX_MUL2 (in->CY[indx][ia], ct1, ct2);
	  COMPLEX_EXP (ct1, -PDRow->Imag2[kndx]);
	  COMPLEX_SET (ct2, sin(G_PI*0.25-PDRow->Real2[kndx]), 0.0);
	  COMPLEX_MUL3 (in->SY[indx][ia], Jp, ct1, ct2);
	  COMPLEX_CONJUGATE (in->CXc[indx][ia], in->CX[indx][ia]);
	  COMPLEX_CONJUGATE (in->SXc[indx][ia], in->SX[indx][ia]);
	  COMPLEX_CONJUGATE (in->CYc[indx][ia], in->CY[indx][ia]);
	  COMPLEX_CONJUGATE (in->SYc[indx][ia], in->SY[indx][ia]);
	  COMPLEX_EXP (in->PPXY[indx],  PD);
	  COMPLEX_EXP (in->PPYX[indx], -PD);
	}
	indx++; kndx++;
      } /* end Channel loop */
    } /* end IF loop */
  } /* End loop over table */

  /* Close table */
  retCode = ObitTablePDClose (PDTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, PDTab->name);
  
  /* release objects */
  PDRow  = ObitTablePDRowUnref(PDRow);
  PDTab  = ObitTablePDUnref(PDTab);
} /* end digestPDTable*/

/**
 * Calculate visibility model using a thread argument
 * \param args   argument
 * \param ant1   first antenna (0-rel)
 * \param ant2   second antenna (0-rel)
 * \param iChan  Channel/IF number (0-rel]
 * \param stokes Stokes visibility (VRR, VLL, VRL, VLR)
 * \param model  output array (VRR, VLL, VRL, VLR)
 */
static void calcModel (VMPolnFTFuncArg* arg, olong ant1, olong ant2, 
		       olong iChan, ofloat stokes[8], ofloat model[8])
{
  ObitSkyModelVMPoln *in =  arg->in;
  dcomplex *RS, *RD, *LS, *LD, *RSc, *RDc, *LSc, *LDc;
  dcomplex *CX, *SX, *CY, *SY, *CXc, *CYc, *SXc, *SYc;
  dcomplex PPRL, PPLR, PPXY, PPYX;
  dcomplex PA1, PA2, PA1c, PA2c, ct1;
  dcomplex S[4], VRR, VRL, VLR, VLL, MC1, MC2, MC3, MC4;
  dcomplex VXX, VYY, VXY, VYX, Jp, Jm, SPA, DPA, SPAc, DPAc;
  dcomplex SM1, SM2, SM3, SM4;

  /* Stokes model as RR,RL,LR,LL */
  COMPLEX_SET (S[0], stokes[0], stokes[1]);
  COMPLEX_SET (S[3], stokes[2], stokes[3]);
  COMPLEX_SET (S[1], stokes[4], stokes[5]);
  COMPLEX_SET (S[2], stokes[6], stokes[7]);
 
  COMPLEX_SET (MC1, 0.0, 0.0);  /* Muller matrix elements */
  COMPLEX_SET (MC2, 0.0, 0.0);
  COMPLEX_SET (MC3, 0.0, 0.0);
  COMPLEX_SET (MC4, 0.0, 0.0);

  /* Circular or linear */
  if (in->isCirc) {
    RS  = in->RS[iChan];
    RD  = in->RD[iChan];
    LS  = in->LS[iChan];
    LD  = in->LD[iChan];
    RSc = in->RSc[iChan];
    RDc = in->RDc[iChan];
    LSc = in->LSc[iChan];
    LDc = in->LDc[iChan];
    COMPLEX_SET (PPRL, in->PPRL[iChan].real, in->PPRL[iChan].imag);
    COMPLEX_SET (PPLR, in->PPLR[iChan].real, in->PPLR[iChan].imag);
    
    COMPLEX_SET (PA1, arg->cosPA[ant1], arg->sinPA[ant1]); /* Parallactic angle */
    COMPLEX_SET (PA2, arg->cosPA[ant2], arg->sinPA[ant2]);
    COMPLEX_CONJUGATE (PA1c, PA1);
    COMPLEX_CONJUGATE (PA2c, PA2);
    
    /* VRR = S[0] * RS[ant1] * RSc[ant2] +        
             S[1] * RS[ant1] * RDc[ant2] * PA2c + 
             S[2] * RD[ant1] * RSc[ant2] * PA1  + 
             S[3] * RD[ant1] * RDc[ant2] * PA1  * PA2c; */
    COMPLEX_MUL2 (MC1, RS[ant1], RSc[ant2]);
    COMPLEX_MUL2 (VRR, S[0], MC1);
    COMPLEX_MUL3 (MC2, RS[ant1], RDc[ant2],  PA2c);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VRR, VRR,  ct1);
    COMPLEX_MUL3 (MC3, RD[ant1], RSc[ant2], PA1);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VRR, VRR,  ct1);
    COMPLEX_MUL4 (MC4, RD[ant1], RDc[ant2], PA1, PA2c);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VRR, VRR,  ct1);
    model[0] = VRR.real;
    model[1] = VRR.imag;
    
    /* VLL = S[0] * LS[ant1] * LSc[ant2] * PA1c * PA2 +	
             S[1] * LS[ant1] * LDc[ant2] * PA1c +
             S[2] * LD[ant1] * LSc[ant2] * PA2  +
             S[3] * LD[ant1] * LDc[ant2]; */
    COMPLEX_MUL4 (MC1, LS[ant1], LSc[ant2], PA1c, PA2);
    COMPLEX_MUL2 (VLL, S[0], MC1);
    COMPLEX_MUL3 (MC2, LS[ant1], LDc[ant2], PA1c);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VLL, VLL,  ct1);
    COMPLEX_MUL3 (MC3, LD[ant1], LSc[ant2], PA2);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VLL, VLL,  ct1);
    COMPLEX_MUL2 (MC4, LD[ant1], LDc[ant2]);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VLL, VLL,  ct1);
    model[2] = VLL.real;
    model[3] = VLL.imag;
    
    /* 	    RL */
    /* VRL = PPRL * S[0] * RS[ant1] * LSc[ant2] * PA2 +
             PPRL * S[1] * RS[ant1] * LDc[ant2] + 
             PPRL * S[2] * RD[ant1] * LSc[ant2] * PA1 * PA2 +
             PPRL * S[3] * RD[ant1] * LDc[ant2] * PA1; */
    COMPLEX_MUL4 (MC1, PPRL, RS[ant1], LSc[ant2], PA2);
    COMPLEX_MUL2 (VRL, S[0], MC1);
    COMPLEX_MUL3 (MC2, PPRL, RS[ant1], LDc[ant2]);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VRL, VRL,  ct1);
    COMPLEX_MUL5 (MC3, PPRL, RD[ant1], LSc[ant2],  PA1,  PA2);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VRL, VRL,  ct1);
    COMPLEX_MUL4 (MC4, PPRL, RD[ant1], LDc[ant2],  PA1);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VRL, VRL,  ct1);
    model[4] = VRL.real;
    model[5] = VRL.imag;
    
    /*        LR */
    /* VLR = PPLR * S[0] * LS[ant1] * RSc[ant2] * PA1c +
             PPLR * S[1] * LS[ant1] * RDc[ant2] * PA1c * PA2c +
             PPLR * S[2] * LD[ant1] * RSc[ant2] +
             PPLR * S[3] * LD[ant1] * RDc[ant2] * PA2c */
    COMPLEX_MUL4 (MC1, PPLR, LS[ant1], RSc[ant2], PA1c);
    COMPLEX_MUL2 (VLR, S[0], MC1);
    COMPLEX_MUL5 (MC2, PPLR, LS[ant1], RDc[ant2], PA1c,  PA2c);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VLR, VLR,  ct1);
    COMPLEX_MUL3 (MC3, PPLR, LD[ant1], RSc[ant2]);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VLR, VLR,  ct1);
    COMPLEX_MUL4 (MC4, PPLR, LD[ant1], RDc[ant2],  PA2c);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VLR, VLR, ct1);
    model[6] = VLR.real;
    model[7] = VLR.imag;
  }  else {  /* Linear feeds */

    /* Init working variables */
    CX  = in->CX[iChan];
    SX  = in->SX[iChan];
    CY  = in->CY[iChan];
    SY  = in->SY[iChan];
    CXc = in->CXc[iChan];
    SXc = in->SXc[iChan];
    CYc = in->CYc[iChan];
    SYc = in->SYc[iChan];
    COMPLEX_SET (PPXY, in->PPXY[iChan].real, in->PPXY[iChan].imag);
    COMPLEX_SET (PPYX, in->PPYX[iChan].real, in->PPYX[iChan].imag);

    COMPLEX_SET (Jp,  0.0, 1.0);
    COMPLEX_SET (Jm,  0.0,-1.0);
    /* Sum and difference or parallactic angle */
    COMPLEX_SET(SPA, arg->cosPA[ant1]*arg->cosPA[ant2] - arg->sinPA[ant1]*arg->sinPA[ant2],
	  	     arg->sinPA[ant1]*arg->cosPA[ant2] + arg->cosPA[ant1]*arg->sinPA[ant2]);
    COMPLEX_SET(DPA, arg->cosPA[ant1]*arg->cosPA[ant2] + arg->sinPA[ant1]*arg->sinPA[ant2],
	  	     arg->sinPA[ant1]*arg->cosPA[ant2] - arg->cosPA[ant1]*arg->sinPA[ant2]);
    COMPLEX_CONJUGATE (SPAc, SPA);
    COMPLEX_CONJUGATE (DPAc, DPA);
    
    /* Note: no correction for antenna gain error */
    /* VXX = {S[0] * CX[ant1] * CXc[ant2] * DPAc  +
              S[1] * CX[ant1] * SXc[ant2] * SPAc  +
              S[2] * SX[ant1] * CXc[ant2] * SPA   + 
              S[3] * SX[ant1] * SXc[ant2] * DPA};
    */
    COMPLEX_MUL3 (MC1, CX[ant1], CXc[ant2], DPAc);
    COMPLEX_MUL3 (MC2, CX[ant1], SXc[ant2], SPAc);
    COMPLEX_MUL3 (MC3, SX[ant1], CXc[ant2], SPA);
    COMPLEX_MUL3 (MC4, SX[ant1], SXc[ant2], DPA);
    COMPLEX_MUL2 (VXX, S[0], MC1);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VXX, VXX,  ct1);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VXX, VXX,  ct1);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VXX, VXX,  ct1);
    model[0] = VXX.real;
    model[1] = VXX.imag;
    
    /* VYY = {S[0] * SY[ant1] * SYc[ant2] * DPAc +       
              S[1] * SY[ant1] * CYc[ant2] * SPAc +
              S[2] * CY[ant1] * SYc[ant2] * SPA  + 
              S[3] * CY[ant1] * CYc[ant2] * DPA};
    */
    COMPLEX_MUL3 (MC1, SY[ant1], SYc[ant2], DPAc);
    COMPLEX_MUL3 (MC2, SY[ant1], CYc[ant2], SPAc);
    COMPLEX_MUL3 (MC3, CY[ant1], SYc[ant2], SPA);
    COMPLEX_MUL3 (MC4, CY[ant1], CYc[ant2], DPA);
    COMPLEX_MUL2 (VYY, S[0], MC1);
    COMPLEX_MUL2 (ct1, S[1], MC2);
    COMPLEX_ADD2 (VYY, VYY,  ct1);
    COMPLEX_MUL2 (ct1, S[2], MC3);
    COMPLEX_ADD2 (VYY, VYY,  ct1);
    COMPLEX_MUL2 (ct1, S[3], MC4);
    COMPLEX_ADD2 (VYY, VYY,  ct1);
    model[6] = VYY.real;
    model[7] = VYY.imag;
    
    /* VXY = {S[0] * CX[ant1] * SYc[ant2] * DPAc +       
              S[1] * CX[ant1] * CYc[ant2] * SPAc +
              S[2] * SX[ant1] * SYc[ant2] * SPA  + 
              S[3] * SX[ant1] * CYc[ant2] * DPA}} * PPXY;
    */
    COMPLEX_MUL3 (MC1, CX[ant1], SYc[ant2], DPAc);
    COMPLEX_MUL3 (MC2, CX[ant1], CYc[ant2], SPAc);
    COMPLEX_MUL3 (MC3, SX[ant1], SYc[ant2], SPA);
    COMPLEX_MUL3 (MC4, SX[ant1], CYc[ant2], DPA);
    COMPLEX_MUL2 (SM1, S[0], MC1);
    COMPLEX_MUL2 (SM2, S[1], MC2);
    COMPLEX_MUL2 (SM3, S[2], MC3);
    COMPLEX_MUL2 (SM4, S[3], MC4);
    COMPLEX_ADD4 (VXY, SM1, SM2, SM3, SM4);
    COMPLEX_MUL2 (VXY, VXY, PPXY);
    model[2] = VXY.real;
    model[3] = VXY.imag;
    
    /* VYX = {S[0] * SY[ant1] * CXc[ant2] * DPAc +       
              S[1] * SY[ant1] * SXc[ant2] * SPAc +
              S[2] * CY[ant1] * CXc[ant2] * SPA  + 
              S[3] * CY[ant1] * SXc[ant2] * DPA} * PPYX;
    */
    COMPLEX_MUL3 (MC1, SY[ant1], CXc[ant2], DPAc);
    COMPLEX_MUL3 (MC2, SY[ant1], SXc[ant2], SPAc);
    COMPLEX_MUL3 (MC3, CY[ant1], CXc[ant2], SPA);
    COMPLEX_MUL3 (MC4, CY[ant1], SXc[ant2], DPA);
    COMPLEX_MUL2 (SM1, S[0], MC1);
    COMPLEX_MUL2 (SM2, S[1], MC2);
    COMPLEX_MUL2 (SM3, S[2], MC3);
    COMPLEX_MUL2 (SM4, S[3], MC4);
    COMPLEX_ADD4 (VYX, SM1, SM2, SM3, SM4);
    COMPLEX_MUL2 (VYX, VYX, PPYX);
    model[4] = VYX.real;
    model[5] = VYX.imag;
  } /* end linear */
} /* end calCmodel */

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
 * \param mosaic   Mosaic for image
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
static ObitTableCC* 
ObitSkyModelVMPolngetPBCCTab (ObitSkyModel* in, ObitImageMosaic *mosaic, ObitUV* uvdata, 
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
  tempTable = newObitImageTable (mosaic->images[field],OBIT_IO_ReadOnly, 
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
} /* end ObitSkyModelVMPolngetPBCCTab */

