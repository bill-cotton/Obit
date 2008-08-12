/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitSkyModelVMSquint.h"
#include "ObitTableCCUtil.h"
#include "ObitFFT.h"
#include "ObitUVUtil.h"
#include "ObitImageUtil.h"
#include "ObitPBUtil.h"
#include "ObitMem.h"
#include "ObitPrecess.h"
#include "ObitTableANUtil.h"
#include "ObitSkyGeom.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVMSquint.c
 * ObitSkyModelVMSquint class function definitions.
 *
 * This class is derived from the #ObitSkyModelVM class
 *
 * This class represents sky models incorporating beam squint corrections and 
 * their Fourier transforms.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSkyModelVMSquint";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitSkyModelVMGetClass;

/**
 * ClassInfo structure ObitSkyModelVMSquintClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSkyModelVMSquintClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: FT by DFT, may be overridden in derived class */
void ObitSkyModelVMSquintFTDFT (ObitSkyModelVM *in, olong field, 
				ObitUV *uvdata, ObitErr *err);

/** Private: Initialize newly instantiated object. */
void  ObitSkyModelVMSquintInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSkyModelVMSquintClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSkyModelVMSquintClassInfoDefFn (gpointer inClass);

/** Private: Get Inputs. */
void  ObitSkyModelVMSquintGetInput (ObitSkyModel* inn, ObitErr *err);

/** Private: Threaded FTDFT */
static gpointer ThreadSkyModelVMSquintFTDFT (gpointer arg);

/** Private: Get Azimuth of feed for VLA */
static ofloat FeedAz (ObitSkyModelVMSquint* in, ObitUV *uvdata, ofloat *squint);

/** Private: Get Azimuth of feed for EVLA. */
static ofloat FeedAzE (ObitSkyModelVMSquint* in,  ObitUV *uvdata, ofloat *squint);

/** Private: Get Angle from center of the beam */
ofloat BeamAngle (ObitImageDesc *in, ofloat x, ofloat y, ofloat offx, ofloat offy);

/*---------------Private structures----------------*/
/* FT threaded function argument 
 Note: Derived classes MUST have the following entries at the beginning 
 of the corresponding structure */
typedef struct {
  /* type "squint" in this class */
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
  /* thread number  */
  olong        ithread;
  /* Obit error stack object */
  ObitErr      *err;
  /* End time (days) of validity of model */
  ofloat endVMModelTime;
  /* Thread copy of Components list - not used here */
  ObitFArray *VMComps;
  /** Dimension of Rgain...  */
  olong dimGain;
  /** Array of time/spatially variable VLA R component gain */
  ofloat *Rgain;
  /** Array of time/spatially variable VLA L component gain */
  ofloat *Lgain;
  /** Array of time/spatially variable EVLA R component gain */
  ofloat *REgain;
  /** Array of time/spatially variable EVLA L component gain */
  ofloat *LEgain;
} VMSquintFTFuncArg;
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSkyModelVMSquint* newObitSkyModelVMSquint (gchar* name)
{
  ObitSkyModelVMSquint* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMSquintClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSkyModelVMSquint));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSkyModelVMSquintInit((gpointer)out);

 return out;
} /* end newObitSkyModelVMSquint */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSkyModelVMSquintGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMSquintClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSkyModelVMSquintGetClass */

/**
 * Make a deep copy of an ObitSkyModelVMSquint.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSkyModelVMSquint* 
ObitSkyModelVMSquintCopy  (ObitSkyModelVMSquint *in, 
			   ObitSkyModelVMSquint *out, ObitErr *err)
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
  out->Threshold = in->Threshold;
  out->maxResid  = in->maxResid;
  /* Component arrays */
  out->dimGain = in->dimGain;
  if (in->dimGain>0) {
    if (out->Rgain) g_free(out->Rgain); out->Rgain = NULL;
    out->Rgain = g_malloc0(in->dimGain*sizeof(ofloat));
    for (i=0; i<in->dimGain; i++) out->Rgain[i] = in->Rgain[i];
    if (out->Lgain) g_free(out->Lgain); out->Lgain = NULL;
    out->Lgain = g_malloc0(in->dimGain*sizeof(ofloat));
    for (i=0; i<in->dimGain; i++) out->Lgain[i] = in->Lgain[i];
    if (out->REgain) g_free(out->REgain); out->REgain = NULL;
    out->REgain = g_malloc0(in->dimGain*sizeof(ofloat));
    for (i=0; i<in->dimGain; i++) out->Rgain[i] = in->REgain[i];
    if (out->LEgain) g_free(out->LEgain); out->LEgain = NULL;
    out->LEgain = g_malloc0(in->dimGain*sizeof(ofloat));
    for (i=0; i<in->dimGain; i++) out->LEgain[i] = in->LEgain[i];
  }

  return out;
} /* end ObitSkyModelVMSquintCopy */

/**
 * Creates an ObitSkyModelVMSquint 
 * \param name  An optional name for the object.
 * \param mosaic ObitImageMosaic giving one or more images/CC tables
 * \return the new object.
 */
ObitSkyModelVMSquint* 
ObitSkyModelVMSquintCreate (gchar* name, ObitImageMosaic* mosaic)
{
  ObitSkyModelVMSquint* out;
  olong number, i;

  /* Create basic structure */
  out = newObitSkyModelVMSquint (name);

  /* Modify for input mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);
  if ((out->mosaic) && (out->mosaic->numberImages>0)) {
    number = out->mosaic->numberImages;
    out->CCver = ObitMemAlloc0 (sizeof(olong)*number);
    for (i=0; i<number; i++) out->CCver[i] = 0;
    out->startComp = ObitMemAlloc0 (sizeof(olong)*number);
    out->endComp   = ObitMemAlloc0 (sizeof(olong)*number);
  }

  return out;
} /* end ObitSkyModelVMSquintCreate */

/**
 * Initializes Sky Model
 * Checks that data contain RR, LL , save calibration/selection request
 * and set uv data for no selection/calibration
 * \param in      SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMSquintInitMod (ObitSkyModel* inn, ObitUV *uvdata, 
				  ObitErr *err)
{
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong iver;
  gboolean btemp;
  olong itemp, numAntList;
  ObitTableList *list=NULL;
  ObitTableAN *TableAN=NULL;
  ObitUVDesc *uvDesc;
  gchar *blank="    ";
  olong i;
  VMSquintFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMSquintInitMod";

  if (err->error) return;

  /* Save/reset calibration state */
  in->saveDoCalSelect = FALSE;
  ObitInfoListGetTest(uvdata->info, "doCalSelect", &type, dim, &in->saveDoCalSelect);
  dim[0] = dim[1] = dim[2] = 1;
  btemp = FALSE;
  ObitInfoListAlwaysPut (uvdata->info, "doCalSelect", OBIT_bool, dim, &btemp);
  in->saveDoCalib = FALSE;
  ObitInfoListGetTest(uvdata->info, "doCalib", &type, dim, &in->saveDoCalib);
  dim[0] = dim[1] = dim[2] = 1;
  itemp = -1;
  ObitInfoListAlwaysPut (uvdata->info, "doCalib", OBIT_long, dim, &itemp);
  strncpy (in->saveStokes, "    ", 4);
  ObitInfoListGetTest(uvdata->info, "Stokes", &type, dim, in->saveStokes);
  dim[0] = 4; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "Stokes", OBIT_string, dim, blank);

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
    in->threadArgs = g_malloc0(in->nThreads*sizeof(VMSquintFTFuncArg*));
    for (i=0; i<in->nThreads; i++) 
      in->threadArgs[i] = g_malloc0(sizeof(VMSquintFTFuncArg)); 
  } /* end initialize */
  
  for (i=0; i<in->nThreads; i++) {
    args = (VMSquintFTFuncArg*)in->threadArgs[i];
    strcpy (args->type, "squint");  /* Enter type as first entry */
    args->in     = inn;
    args->uvdata = uvdata;
    args->ithread= i;
    args->err    = err;
    args->endVMModelTime = -1.0e20;
    args->VMComps= NULL;
    args->dimGain= 0;
    args->Rgain  = NULL;
    args->Lgain  = NULL;
    args->REgain = NULL;
    args->LEgain = NULL;
  }

  /* Call parent initializer */
  ObitSkyModelVMInitMod(inn, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Fourier transform routines - DFT only */
  in->DFTFunc  = (ObitThreadFunc)ThreadSkyModelVMSquintFTDFT;

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

} /* end ObitSkyModelVMSquintInitMod */

/**
 * Any shutdown operations needed for a model
 * Restore calibration/selection state
 * \param in  SkyModel to shutdown
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMSquintShutDownMod (ObitSkyModel* inn, ObitUV *uvdata,
				      ObitErr *err)
{
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)inn;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i;
  VMSquintFTFuncArg *args;

  /* Call parent shutdown */
  ObitSkyModelVMShutDownMod(inn, uvdata, err);

  if (in->threadArgs) {
    /* Check type - only handle "squint" */
    if (!strncmp((gchar*)in->threadArgs[0], "squint", 6)) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMSquintFTFuncArg*)in->threadArgs[i];
	if (args->Rgain)  g_free(args->Rgain);
	if (args->Lgain)  g_free(args->Lgain);
	if (args->REgain) g_free(args->REgain);
	if (args->LEgain) g_free(args->LEgain);
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
    } /* end if this a "squint" threadArg */
  }

  /* Restore calibration state */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "doCalSelect", OBIT_bool, dim, &in->saveDoCalSelect);
  ObitInfoListAlwaysPut (uvdata->info, "doCalib", OBIT_long, dim, &in->saveDoCalib);
  dim[0] = 4; dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "Stokes", OBIT_string, dim, in->saveStokes);

} /* end ObitSkyModelVMSquintShutDownMod */

/**
 * Initializes an ObitSkyModel for a pass through data in time order.
 * Resets current times
 * \param in  SkyModel to initialize
 * \param err Obit error stack object.
 */
void ObitSkyModelVMSquintInitModel (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)inn;
  olong i;
  VMSquintFTFuncArg *args;

  /*  Reset time of current model */
  in->endVMModelTime = -1.0e20;
  in->curVMModelTime = -1.0e20;

  /* Threading */
  if (in->threadArgs) {
    /* Check type - only handle "squint" */
    if (!strncmp((gchar*)in->threadArgs[0], "squint", 6)) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMSquintFTFuncArg*)in->threadArgs[i];
	args->endVMModelTime = -1.0e20;
      }
    } /* end if this a "squint" threadArg */
  }
} /* end ObitSkyModelVMSquintInitModel */

/**
 * Update VM model with time and spatial modifications to model
 * Update model in Rgain, Lgain, REgain, LEgain members from comps member 
 * Sets endVMModelTime for when parallactic angle differs by 1 degree.
 * \param in      SkyModelVM 
 *                Parameters on info:
 *                  PBmin = min. beam gain for Squint correction (Jinc 0.05, Poly 0.01)
 * \param time    current time (d)
 * \param suba    0-rel subarray number
 * \param uvdata  uv data being modeled.
 * \param ithread which thread (0-rel)
 * \param err Obit error stack object.
 */
void ObitSkyModelVMSquintUpdateModel (ObitSkyModelVM *inn, 
				      ofloat time, olong suba,
				      ObitUV *uvdata, olong ithread, 
				      ObitErr *err)
{
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)inn;
  odouble Freq, Angle, AngleRR, AngleLL;
  olong npos[2], lcomp, ncomp, i, ifield;
  ofloat *ccData, *Rgain, *Lgain, curPA, tPA, tTime;
  ofloat feedPA, squint, dx, dy, x, y, antsize = 24.5, pbmin = 0.0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  VMSquintFTFuncArg *args;
  gchar *routine = "ObitSkyModelVMSquintUpdateModel";

  /* DEBUG
  ofloat RFact, LFact; */
  if (err->error) return;

  /* Thread to update */
  args = (VMSquintFTFuncArg*)in->threadArgs[ithread];
  g_assert (!strncmp(args->type,"squint",6));  /* Test arg */

   /* Check subarray */
  Obit_return_if_fail(((suba>=0) && (suba<in->numAntList)), err,
		      "%s: bad subarray number %d", routine, suba+1);

  /* Check gain array dimension */
  Obit_return_if_fail((args->dimGain>=in->numComp), err,
    "%s:More components %d than allowed in gains  %d", 
    routine, in->numComp, args->dimGain);

  /* Min gain for squint correction - don't go below 1% */
  ObitInfoListGetTest(in->info, "PBmin", &type, dim, &pbmin);
  pbmin = MAX (pbmin, 0.01);

  /* Need Parallactic angle */
  tTime = time - 5.0*suba;  /* Correct time for subarray offset  from DBCON */
  curPA = ObitAntennaListParAng (in->AntList[suba], 1, tTime, in->curSource);
  Freq = uvdata->myDesc->crval[uvdata->myDesc->jlocf];

  /* Which antennas are EVLA ? */
  if (in->isEVLA==NULL) {
    ObitThreadLock(in->thread);  /* Lock against other threads */
    if (in->isEVLA==NULL)   /* Still NULL? */
      in->isEVLA = g_malloc((100+in->AntList[suba]->number)*sizeof(gboolean));
    ObitThreadUnlock(in->thread); 
  }
  for (i=0; i<in->AntList[suba]->number; i++) {
    in->isEVLA[i] = 
      (in->AntList[suba]->ANlist[i]->AntName[0]=='E') &&
      (in->AntList[suba]->ANlist[i]->AntName[1]=='V') &&
      (in->AntList[suba]->ANlist[i]->AntName[2]=='L') &&
      (in->AntList[suba]->ANlist[i]->AntName[3]=='A');
  }

  /* Compute VLA corrections */
  /* What azimuth is the feed at? */
  feedPA = DG2RAD*FeedAz (in, uvdata, &squint);
  feedPA -= curPA; /* Correct to on sky */

  /* DEBUG
  squint = 0.0;
  squint = -squint; *//* ?? */
  /* feedPA -= 2.0*curPA; Correct to on sky */
  /* The sign of this angle is not well tested as it only matters to x
     and the source in the test data is due south (y only).*/
  feedPA = -feedPA; /* Correct to on sky */

  /* Offsets due to beam squint - for LL + for RR */
  dx = squint * sin(feedPA);
  dy = squint * cos(feedPA);

  npos[0] = 0; npos[1] = 0; 
  ccData = ObitFArrayIndex(in->comps, npos);
  Rgain = args->Rgain;
  Lgain = args->Lgain;
  lcomp = in->comps->naxis[0];  /* Length of row in comp table */
  ncomp = in->numComp;  /* number of components */

  /* Compute antenna gains and put in to Rgain, Lgain */
  for (i=0; i<ncomp; i++) {
    Lgain[i] = 1.0;
    Rgain[i] = 1.0;
    /* Where in the beam? */
    x = ccData[i*lcomp+1];
    y = ccData[i*lcomp+2];
    ifield = ccData[i*lcomp+0]+0.5;
    if (ifield<0) continue;
    /* Determine angle wrt pointing position */
    /*Angle = ObitImageDescAngle(in->mosaic->images[ifield]->myDesc, y, x);*/
    Angle = BeamAngle(in->mosaic->images[ifield]->myDesc, x, y, 0.0, 0.0);

    /* Angle with squint LL */
    /*AngleLL = ObitImageDescAngle(in->mosaic->images[ifield]->myDesc, 
      y-dy, x-dx);*/
    AngleLL = BeamAngle(in->mosaic->images[ifield]->myDesc, x, y, -dx, -dy);
    Lgain[i] = sqrt(ObitPBUtilPntErr (Angle, AngleLL, antsize, pbmin, Freq));

    /* Angle with squint RR */
    /* AngleRR = ObitImageDescAngle(in->mosaic->images[ifield]->myDesc, 
       y+dy, x+dx);*/
    AngleRR = BeamAngle(in->mosaic->images[ifield]->myDesc, x, y, dx, dy);
    Rgain[i] = sqrt(ObitPBUtilPntErr (Angle, AngleRR, antsize, pbmin, Freq));

    /* DEBUG 
    if (i==0) {
      RFact  =  ObitPBUtilPntErr (Angle, AngleRR, antsize, pbmin, Freq);
      LFact  =  ObitPBUtilPntErr (Angle, AngleLL, antsize, pbmin, Freq);
      fprintf (stdout, "%f %f %f %f %f %f %f\n",time*24.0, feedPA*RAD2DG, 
	       Angle, AngleRR, AngleLL, RFact, LFact);
    } */

  } /* end loop over components for VLA */

  /* Compute EVLA corrections */
  /* What azimuth is the feed at? */
  feedPA = DG2RAD*FeedAzE (in, uvdata, &squint);
  feedPA -= curPA; /* Correct to on sky */

  /* DEBUG
  squint = 0.0;
  squint = -squint; *//* ?? */
  /* feedPA -= 2.0*curPA; Correct to on sky */
  /* The sign of this angle is not well tested as it only matters to x
     and the source in the test data is due south (y only).*/
  feedPA = -feedPA; /* Correct to on sky */

  /* Offsets due to beam squint - for LL + for RR */
  dx = squint * sin(feedPA);
  dy = squint * cos(feedPA);

  npos[0] = 0; npos[1] = 0; 
  ccData = ObitFArrayIndex(in->comps, npos);
  Rgain = args->REgain;
  Lgain = args->LEgain;
  lcomp = in->comps->naxis[0];  /* Length of row in comp table */
  ncomp = in->numComp;          /* number of components */

  /* Compute antenna gains and put in to Rgain, Lgain */
  for (i=0; i<ncomp; i++) {
    Lgain[i] = 1.0;
    Rgain[i] = 1.0;
    /* Where in the beam? */
    x = ccData[i*lcomp+1];
    y = ccData[i*lcomp+2];
    ifield = ccData[i*lcomp+0]+0.5;
    if (ifield<0) continue;
    /* Determine angle wrt pointing position */
    /*Angle = ObitImageDescAngle(in->mosaic->images[ifield]->myDesc, y, x);*/
    Angle = BeamAngle(in->mosaic->images[ifield]->myDesc, x, y, 0.0, 0.0);

    /* Angle with squint LL */
    /*AngleLL = ObitImageDescAngle(in->mosaic->images[ifield]->myDesc, 
      y-dy, x-dx);*/
    AngleLL = BeamAngle(in->mosaic->images[ifield]->myDesc, x, y, -dx, -dy);
    Lgain[i] = sqrt(ObitPBUtilPntErr (Angle, AngleLL, antsize, pbmin, Freq));

    /* Angle with squint RR */
    /* AngleRR = ObitImageDescAngle(in->mosaic->images[ifield]->myDesc, 
       y+dy, x+dx);*/
    AngleRR = BeamAngle(in->mosaic->images[ifield]->myDesc, x, y, dx, dy);
    Rgain[i] = sqrt(ObitPBUtilPntErr (Angle, AngleRR, antsize, pbmin, Freq));

    /* DEBUG 
    if (i==0) {
      RFact  =  ObitPBUtilPntErr (Angle, AngleRR, antsize, pbmin, Freq);
      LFact  =  ObitPBUtilPntErr (Angle, AngleLL, antsize, pbmin, Freq);
      fprintf (stdout, "%f %f %f %f %f %f %f\n",time*24.0, feedPA*RAD2DG, 
	       Angle, AngleRR, AngleLL, RFact, LFact);
    } */

  } /* end loop over components for EVLA */

  /* end time - need Parallactic angle */
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
  
} /* end ObitSkyModelVMSquintUpdateModel */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSkyModelVMSquintClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSkyModelVMSquintClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSkyModelVMSquintClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSkyModelVMSquintClassInfoDefFn (gpointer inClass)
{
  ObitSkyModelVMSquintClassInfo *theClass = (ObitSkyModelVMSquintClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSkyModelVMSquintClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSkyModelVMSquintClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSkyModelVMSquintGetClass;
  theClass->newObit       = (newObitFP)newObitSkyModelVMSquint;
  theClass->ObitCopy      = (ObitCopyFP)ObitSkyModelVMSquintCopy;
  theClass->ObitClear     = (ObitClearFP)ObitSkyModelVMSquintClear;
  theClass->ObitInit      = (ObitInitFP)ObitSkyModelVMSquintInit;
  theClass->ObitSkyModelCreate       = (ObitSkyModelCreateFP)ObitSkyModelVMSquintCreate;
  theClass->ObitSkyModelInitMod      = (ObitSkyModelInitModFP)ObitSkyModelVMSquintInitMod;
  theClass->ObitSkyModelShutDownMod  = (ObitSkyModelShutDownModFP)ObitSkyModelVMSquintShutDownMod;
  theClass->ObitSkyModelInitModel    = (ObitSkyModelInitModelFP)ObitSkyModelVMSquintInitModel;
  theClass->ObitSkyModelFTDFT        = (ObitSkyModelFTDFTFP)ObitSkyModelVMSquintFTDFT;
  theClass->ObitSkyModelGetInput     = (ObitSkyModelGetInputFP)ObitSkyModelVMSquintGetInput;
  theClass->ObitSkyModelChose        = (ObitSkyModelChoseFP)ObitSkyModelVMSquintChose;
  theClass->ObitSkyModelVMUpdateModel=
    (ObitSkyModelVMUpdateModelFP)ObitSkyModelVMSquintUpdateModel;

} /* end ObitSkyModelVMSquintClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSkyModelVMSquintInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModelVMSquint *in = inn;

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
  in->REgain       = NULL;
  in->LEgain       = NULL;
  in->AntList      = NULL;
  in->curSource    = NULL;
  in->numAntList   = 0;
  in->Threshold    = 0.0;
  in->maxResid     = 0.0;
} /* end ObitSkyModelVMSquintInit */


/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSkyModelVMSquint* cast to an Obit*.
 */
void ObitSkyModelVMSquintClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  VMSquintFTFuncArg *args;
  ObitSkyModelVMSquint *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->isEVLA) g_free(in->isEVLA); in->isEVLA = NULL;
  if (in->Rgain)  g_free(in->Rgain);  in->Rgain  = NULL;
  if (in->Lgain)  g_free(in->Lgain);  in->Lgain  = NULL;
  if (in->REgain) g_free(in->REgain); in->REgain = NULL;
  if (in->LEgain) g_free(in->LEgain); in->LEgain = NULL;
  in->curSource = ObitSourceUnref(in->curSource);
  if (in->AntList)  {
    for (i=0; i<in->numAntList; i++) { 
      in->AntList[i] = ObitAntennaListUnref(in->AntList[i]);
    }
    g_free(in->AntList); in->AntList = NULL;
  }
    
  /* Thread stuff */
  if (in->threadArgs) {
    /* Check type - only handle "squint" */
    if (!strncmp((gchar*)in->threadArgs[0], "squint", 6)) {
      for (i=0; i<in->nThreads; i++) {
	args = (VMSquintFTFuncArg*)in->threadArgs[i];
	if (args->Rgain)  g_free(args->Rgain);
	if (args->Lgain)  g_free(args->Lgain);
	if (args->REgain) g_free(args->REgain);
	if (args->LEgain) g_free(args->LEgain);
	g_free(in->threadArgs[i]);
      }
      g_free(in->threadArgs);
      in->threadArgs = NULL;
    } /* end if this a "ion" threadArg */
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSkyModelVMSquintClear */

/**
 * Get input parameters from info member
 * If maxResid value not given or <0 then for it is determined
 * by examining each image in the Image mosaic.  
 * If any image has an info item "maxAbsResid" the the maximum of any
 * of these is used, else the MAX (fabs(minval), fabs(maxval)) is used.
 * \param in Pointer to the ObitSkyModelVMSquint .
 * \param err Obit error stack object.
 */
void  ObitSkyModelVMSquintGetInput (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)inn;
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM];
  ObitImageDesc *desc;
  gboolean lookup, gotit;
  ofloat maxv, maxAbsResid;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitSkyModelVMSquintGetInput";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitSkyModelVMSquintIsA(in));
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
  
} /* end ObitSkyModelVMSquintGetInput */

/**
 * Decide which method is the most appropriate to calculate the FT of a model
 * Sets currentMode member function
 * Only DFT supported
 * \param in     Pointer to the ObitSkyModel .
 * \param uvdata UV data set
 */
void  ObitSkyModelVMSquintChose (ObitSkyModel* inn, ObitUV* uvdata) 
{
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)inn;
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

} /* end ObitSkyModelVMSquintChose */


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
void ObitSkyModelVMSquintFTDFT (ObitSkyModelVM *inn, olong field, ObitUV *uvdata, ObitErr *err)
{
  olong i, mcomp, iComp, pos[2], nvis, lovis, hivis, nvisPerThread;
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)inn;
  VMSquintFTFuncArg *args;
  ofloat *ddata;
  gboolean OK = TRUE;
  gchar *routine = "ObitSkyModelVMSquintFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

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
  nvisPerThread = nvis/in->nThreads;
  lovis = 1;
  hivis = nvisPerThread;
  hivis = MIN (hivis, nvis);

  /* Set up thread arguments */
  for (i=0; i<in->nThreads; i++) {
    if (i==(in->nThreads-1)) hivis = nvis;  /* Make sure do all */
    args = (VMSquintFTFuncArg*)in->threadArgs[i];
    args->in     = (ObitSkyModel*)inn;
    args->field  = field;
    args->uvdata = uvdata;
    args->first  = lovis;
    args->last   = hivis;
    args->ithread= i;
    args->err    = err;
    if (args->dimGain!=in->numComp) {
      if (args->Rgain)  g_free(args->Rgain);
      if (args->Lgain)  g_free(args->Lgain);
      if (args->REgain) g_free(args->REgain);
      if (args->LEgain) g_free(args->LEgain);
      args->dimGain = in->numComp;
      args->Rgain  = g_malloc0(args->dimGain*sizeof(ofloat));
      args->Lgain  = g_malloc0(args->dimGain*sizeof(ofloat));
      args->REgain = g_malloc0(args->dimGain*sizeof(ofloat));
      args->LEgain = g_malloc0(args->dimGain*sizeof(ofloat));
    }
    /* Update which vis */
    lovis += nvisPerThread;
    hivis += nvisPerThread;
    hivis = MIN (hivis, nvis);
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, in->nThreads, in->DFTFunc, in->threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
}  /* end ObitSkyModelVMSquintFTDFT */

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
 * Arguments are given in the VMSquintFTFuncArg structure passed as arg starting 
 * with the following:
 * \li type   String identifying structure
 * \li in     SkyModelVM with model components loaded (ObitSkyModelLoad)
 * \li field  Field number being processed (-1 => all)
 * \li uvdata UV data set to model and subtract from current buffer
 * \li first  First (1-rel) vis in uvdata buffer to process this thread
 * \li last   Highest (1-rel) vis in uvdata buffer to process this thread
 * \li ithread thread number
 * \li err    Obit error stack object.
 * \li endVMModelTime End time (days) of validity of model
 * \li VMComps Thread copy of Components list - not used here
 * \li dimGain Dimension of Rgain
 * \li Rgain   Float array of time/spatially variable VLA R component gain
 * \li Lgain   Float array of time/spatially variable VLA L component gain
 * \li REgain  Float array of time/spatially variable EVLA R component gain
 * \li LEgain  Float array of time/spatially variable EVLA L component gain
 * \return NULL
 */
static gpointer ThreadSkyModelVMSquintFTDFT (gpointer args)
{
  /* Get arguments from structure */
  VMSquintFTFuncArg *largs = (VMSquintFTFuncArg*)args;
  ObitSkyModelVMSquint *in = (ObitSkyModelVMSquint*)largs->in;
  /*olong field      = largs->field;*/
  ObitUV *uvdata   = largs->uvdata;
  olong loVis      = largs->first-1;
  olong hiVis      = largs->last;
  olong ithread    = largs->ithread;
  ObitErr *err     = largs->err;
  /*olong dimGain    = largs->dimGain;*/
  ofloat *Rgain    = largs->Rgain;
  ofloat *Lgain    = largs->Lgain;
  ofloat *REgain   = largs->REgain;
  ofloat *LEgain   = largs->LEgain;

  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, ncomp;
  olong lrec, nrparm, naxis[2];
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF;
  olong ilocu, ilocv, ilocw, iloct, ilocb, suba, itemp, ant1, ant2, mcomp;
  ofloat *visData, *Data, *ddata, *fscale;
  ofloat sumRealRR, sumImagRR, modRealRR, modImagRR;
  ofloat sumRealLL, sumImagLL, modRealLL, modImagLL;
  ofloat cbase, *rgain1, *lgain1, *rgain2, *lgain2;
  ofloat amp, phase, sp, cp, arg, freq2, freqFact, wtRR=0.0, wtLL=0.0, temp;
  odouble *freqArr;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMSquintFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) goto finish;

 /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;
  ilocb =  uvdata->myDesc->ilocb;

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
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF beforefreq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  /* Get pointer for components */
  naxis[0] = 0; naxis[1] = 0; 
  Data = ObitFArrayIndex(in->comps, naxis);
  lcomp = in->comps->naxis[0];   /* Length of row in comp table */
  ncomp = in->comps->naxis[1];   /* Number of components */
  mcomp = in->numComp;           /* Actual number */

  /* Get pointer for frequency correction tables */
  fscale  = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

  /* Loop over vis in buffer */
  lrec    = uvdata->myDesc->lrec;         /* Length of record */
  visData = uvdata->buffer+loVis*lrec;    /* Buffer pointer with appropriate offset */
  nrparm  = uvdata->myDesc->nrparm;       /* Words of "random parameters" */

  for (iVis=loVis; iVis<hiVis; iVis++) {

    /* Is current model still valid? */
    if (visData[iloct] > largs->endVMModelTime) {
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
    for (iIF=startIF; iIF<startIF+numberIF; iIF++) {
      offsetIF = nrparm + iIF*jincif; 
      for (iChannel=startChannel; iChannel<startChannel+numberChannel; iChannel++) {
	offsetChannel = offsetIF + iChannel*jincf; 
	freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */

	/* Sum over components */
	/* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	sumRealRR = sumImagRR = sumRealLL = sumImagLL = 0.0;
	/* Set component gain lists by antenna and type */
	ddata = Data;
	if (in->isEVLA[ant1]) {
	  rgain1 = REgain;
	  lgain1 = LEgain;
	} else {
	  rgain1 = Rgain;
	  lgain1 = Lgain;
	}
	if (in->isEVLA[ant2]) {
	  rgain2 = REgain;
	  lgain2 = LEgain;
	} else {
	  rgain2 = Rgain;
	  lgain2 = Lgain;
	}
	
	/* Sum by model type - assume phase same for RR, LL */
	switch (in->modType) {
	case OBIT_SkyModel_PointMod:     /* Point */
	  /* From the AIPSish QXXPTS.FOR  */
	  for (iComp=0; iComp<mcomp; iComp++) {
	    phase = freqFact * (ddata[4]*visData[ilocu] + 
				ddata[5]*visData[ilocv] + 
				ddata[6]*visData[ilocw]);
	    cp = cos(phase);
	    sp = sin(phase);
	    /* Amplitude from component flux and two gains */
	    amp = ddata[3] * rgain1[iComp] * rgain2[iComp];
	    sumRealRR += amp*cp;
	    sumImagRR += amp*sp;
	    amp = ddata[3] * lgain1[iComp] * lgain2[iComp];
	    sumRealLL += amp*cp;
	    sumImagLL += amp*sp;
	    /* DEBUG
	       if (iVis==0) {
	       fprintf (stderr,"%s: comp %d real %f imag %f phase %f sum %f %f \n",
	       routine, iComp, ccData[0]*cos(phase), ccData[0]*sin(phase),
	       57.296*phase, sumReal, sumImag); 
	       } */
	    ddata += lcomp;   /* update pointer */
	  } /* end loop over components */
	  break;
	case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	  /* From the AIPSish QGASUB.FOR  */
	  freq2 = freqArr[iIF*kincif + iChannel*kincf];    /* Frequency squared */
	  freq2 *= freq2;
	  for (iComp=0; iComp<mcomp; iComp++) {
	    phase = freqFact * (ddata[4]*visData[ilocu] + 
				ddata[5]*visData[ilocv] + 
				ddata[6]*visData[ilocw]);
	    cp = cos(phase);
	    sp = sin(phase);
	    arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			   ddata[8]*visData[ilocv]*visData[ilocv] +
			   ddata[9]*visData[ilocu]*visData[ilocv]);
	    amp = ddata[3] * rgain1[iComp] * rgain2[iComp];
	    if (arg>1.0e-5) amp *= exp (arg);
	    sumRealRR += amp*cp;
	    sumImagRR += amp*sp;
	    arg = freq2 * (ddata[7]*visData[ilocu]*visData[ilocu] +
			   ddata[8]*visData[ilocv]*visData[ilocv] +
			   ddata[9]*visData[ilocu]*visData[ilocv]);
	    amp = ddata[3] * lgain1[iComp] * lgain2[iComp];
	    if (arg>1.0e-5) amp *= exp (arg);
	    sumRealLL += amp*cp;
	    sumImagLL += amp*sp;
	    ddata += lcomp;   /* update pointer */
	  }  /* end loop over components */
	  break;
	case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	  /* From the AIPSish QSPSUB.FOR  */
	  for (iComp=0; iComp<mcomp; iComp++) {
	    phase = freqFact * (ddata[4]*visData[ilocu] + 
				ddata[5]*visData[ilocv] + 
				ddata[6]*visData[ilocw]);
	    cp = cos(phase);
	    sp = sin(phase);
	    arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				  visData[ilocv]*visData[ilocv]);
	    arg = MAX (arg, 0.1);
	    arg = ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
	    amp = ddata[3] * rgain1[iComp] * rgain2[iComp] * arg;
	    sumRealRR += amp*cp;
	    sumImagRR += amp*sp;
	    amp = ddata[3] * lgain1[iComp] * lgain2[iComp] * arg;
	    sumRealLL += amp*cp;
	    sumImagLL += amp*sp;
	    ddata += lcomp;   /* update pointer */
	  }
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

	  offsetChannel += jincf;
      } /* end loop over Channel */
 	  offsetIF += jincif;
   } /* end loop over IF */

    visData += lrec; /* Update vis pointer */
  } /* end loop over visibilities */

  /* Indicate completion */
 finish: ObitThreadPoolDone (in->thread, (gpointer)&ithread);
  
  return NULL;
} /* ObitSkyModelVMquintFTDFT */

/**
 * Return azimuth (antenna coords) and throw of VLA feed being used.
 * \param in     ObitSkyModel giving selected PB IFs and channels
 * \param uvdata ObitUV data with Descriptor including frequency
 * \param squint [out] beam squint (deg) for average PB frequency
 * \return antenna azimuth (deg) for reference freq
 */
static ofloat FeedAz (ObitSkyModelVMSquint* in, ObitUV *uvdata, 
		      ofloat *squint)
{
  ofloat out;
  olong iIF, iChannel;
  olong startIF, numberIF, startChannel, numberChannel, kincf, kincif;
  odouble Freq=1.0, PBFreq;

  /* Beam squint 0.065 times HPBW = (1.5 lambda(cm) in amin) */
  Freq = uvdata->myDesc->crval[uvdata->myDesc->jlocf];  /* Ref. freq */
  /* Get center channel for PB corr - 
     Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF  = in->startIFPB-1;
  numberIF = MAX (1, in->numberIFPB);
  iIF      = startIF + numberIF/2;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  iChannel      = startChannel + numberChannel/2;
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF before freq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  PBFreq = uvdata->myDesc->freqArr[iIF*kincif + iChannel*kincf]; /* center PB Freq */
  /* correct squint - empirically from VLA test data (0.82) */
  (*squint) =  0.82 * 0.5 * 0.065 * (VELIGHT/PBFreq) * 1.5e2  / 60.0;
  out = 0.0;

  /* Get angle by reference frequency */
  if ((Freq>1.0e9) && (Freq<2.0e9)) {
    out = -136.8;   /* VLA L band */
    /* DEBUG = best (out = -132.0; */  /* VLA L band */
    /* At 1.4174 GHz basic formula => 61.87 " Juan likes 51.3" @ 1.400 GHz */
  (*squint) =  0.819 * 0.5 * 0.065 * (VELIGHT/PBFreq) * 1.5e2  / 60.0;
  /* DEBUG = best (*squint) =  0.75 * 0.5 * 0.065 * (VELIGHT/PBFreq) * 1.5e2  / 60.0; */

  } else if ((Freq>2.0e9) && (Freq<3.0e9)) {
    out = 0.0;   /* VLA S band not on VLA */
  } else if ((Freq>3.0e9) && (Freq<7.0e9)) {
    out = 61.6;   /* VLA C band */
  } else if ((Freq>7.0e9) && (Freq<12.0e9)) {
    out = -48.6;   /* VLA X band, ant 6,14,22 = 9.7 */
  } else if ((Freq>12.0e9) && (Freq<20.0e9)) {
    /* out = 113.5; VLA Ku ant 1,2,3,4,6,7,9,10,12,14 */
    /* out = 127.8; VLA Ku band ant 5,8,11,16,17,18,19,20,21,22,23,24,25,26,27,28*/
    out = 120.0;   /* VLA Ku band compromise - poor for everything */
  } else if ((Freq>20.0e9) && (Freq<28.0e9)) {
    /* out =  91.1; VLA K band ant 6,14,22 */
    /* out =  -4.5; VLA K band ant 1,2,3,16,18,24,27,28 */
    out = -25.4; /* VLA K band ant 4,7,8,9,10,11,12,15,17,19,20,21,23,25,26 */
  } else if ((Freq>28.0e9) && (Freq<35.0e9)) {
    out = 0.0; /* VLA KA band - not on VLA */
  } else if ((Freq>35.0e9) && (Freq<50.0e9)) {
    /* out = -20.8; VLA Q band ant 6,14,22 */
    out =   6.0; /* VLA Q band all other antennas */
  } else {
    out = 0.0;
  }
  return out;
} /* end FeedAz */

/**
 * Return azimuth (antenna coords) and throw of EVLA feed being used.
 * \param in     ObitSkyModel giving selected PB IFs and channels
 * \param uvdata ObitUV data with Descriptor including frequency
 * \param squint [out] beam squint (deg) for average PB frequency
 * \return antenna azimuth (deg) for reference freq
 */
static ofloat FeedAzE (ObitSkyModelVMSquint* in, ObitUV *uvdata, 
		       ofloat *squint)
{
  ofloat out;
  olong iIF, iChannel;
  olong startIF, numberIF, startChannel, numberChannel, kincf, kincif;
  odouble Freq=1.0, PBFreq;

  /* Beam squint 0.065 times HPBW = (1.5 lambda(cm) in amin) */
  Freq = uvdata->myDesc->crval[uvdata->myDesc->jlocf];  /* Ref. freq */
  /* Get center channel for PB corr - 
     Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF  = MAX (0, in->startIFPB-1);
  numberIF = MAX (1, in->numberIFPB);
  iIF      = startIF + numberIF/2;
  startChannel  = MAX (0, in->startChannelPB-1);
  numberChannel = MAX (1, in->numberChannelPB);
  iChannel      = startChannel + numberChannel/2;
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF before freq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  PBFreq = uvdata->myDesc->freqArr[iIF*kincif + iChannel*kincf]; /* center PB Freq */
  /* correct squint - empirically from VLA test data (0.82) */
  (*squint) =  0.82 * 0.5 * 0.065 * (VELIGHT/PBFreq) * 1.5e2  / 60.0;
  out = 0.0;

  /* Get angle by reference frequency */
  if ((Freq>1.0e9) && (Freq<2.0e9)) {
    out = -84.1;   /* EVLA L band measured */
  } else if ((Freq>2.0e9) && (Freq<3.0e9)) {
    out = 101.6;   /* EVLA S band  nominal*/
  } else if ((Freq>3.0e9) && (Freq<7.0e9)) {
    out = 165.6;   /* EVLA C band measured*/
  } else if ((Freq>7.0e9) && (Freq<12.0e9)) {
    out = -156.3;   /* EVLA X band, ant 6,14,22 = 9.7 */
  } else if ((Freq>12.0e9) && (Freq<20.0e9)) {
    out = 47.6;   /* EVLA Ku nominal */
  } else if ((Freq>20.0e9) && (Freq<28.0e9)) {
    out =  25.9;  /*  EVLA K band nominal */
  } else if ((Freq>28.0e9) && (Freq<35.0e9)) {
    out = -16.9; /* EVLA KA band nominal */
  } else if ((Freq>35.0e9) && (Freq<50.0e9)) {
    out =   4.5; /* EVLA Q band  nominal */
  } else {
    out = 0.0;
  }
  return out;
} /* end FeedAzE */

/**
 * Determine angle on sky from a position to the antenna pointing position
 * \param in  Image descriptor
 * \param x   X offset in deg from ref pixel in in
 *            This should be in the same system and equinox as in in.
 * \param y   Y offset in deg from ref pixel in in
 * \param err Error stack, returns if not empty.
 * \return angular distance on the sky (deg)
 */
ofloat BeamAngle (ObitImageDesc *in, ofloat x, ofloat y, ofloat offx, ofloat offy)
{
  ofloat dist = 0.0;
  odouble RAPnt, DecPnt, RAPntc, DecPntc, ra, dec, xx, yy, zz;

  /* Get pointing position */
  ObitImageDescGetPoint(in, &RAPnt, &DecPnt);

  /* Actual pointing position with squint */
  ObitSkyGeomXYShift (RAPnt, DecPnt, offx, offy, 0.0, &RAPntc, &DecPntc);
  RAPntc  *= DG2RAD;
  DecPntc *= DG2RAD;

  /* Convert offset to position */
  ObitSkyGeomXYShift (in->crval[in->jlocr], in->crval[in->jlocd],
		      x, y, ObitImageDescRotate(in), &ra, &dec);

  /* Compute distance */
  xx = DG2RAD * (ra);
  yy = DG2RAD * (dec);
  zz = sin (yy) * sin (DecPntc) + cos (yy) * cos (DecPntc) * cos (xx-RAPntc);
  zz = MIN (zz, 1.000);
  dist = acos (zz) * RAD2DG;
  return dist;
} /* end BeamAngle */

