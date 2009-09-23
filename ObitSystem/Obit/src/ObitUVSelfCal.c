/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2009                                          */
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

#include "ObitUVSelfCal.h"
#include "ObitMem.h"
#include "ObitTableUtil.h"
#include "ObitUVUtil.h"
#include "ObitUVSoln.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSelfCal.c
 * ObitUVSelfCal class function definitions.
 * This class enables self calibration of ObitUV data sets
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVSelfCal";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVSelfCalClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVSelfCalClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVSelfCalInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVSelfCalClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVSelfCalClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVSelfCal* newObitUVSelfCal (gchar* name)
{
  ObitUVSelfCal* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVSelfCalClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVSelfCal));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVSelfCalInit((gpointer)out);

 return out;
} /* end newObitUVSelfCal */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVSelfCalGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVSelfCalClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVSelfCalGetClass */

/**
 * Make a deep copy of an ObitUVSelfCal.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVSelfCal* ObitUVSelfCalCopy  (ObitUVSelfCal *in, ObitUVSelfCal *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;

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
    out = newObitUVSelfCal(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->modelMode = in->modelMode;
  out->skyModel = ObitSkyModelUnref(out->skyModel);
  out->skyModel = ObitSkyModelRef(in->skyModel);
  out->SCData   = ObitUVUnref(out->SCData);
  out->SCData   = ObitSkyModelRef(in->SCData);
  out->display  = ObitDisplayUnref(out->display);
  out->display  = ObitDisplayRef(in->display);
  out->hist     = ObitMemFree (out->hist);
  out->hist     = ObitMemAlloc0Name(in->numHist*sizeof(ofloat),"Flux Histogram");
  for (i=0; i<in->numHist; i++) out->hist[i] = in->hist[i];
  out->histRMS  = ObitMemFree (out->histRMS);
  out->histRMS  = ObitMemAlloc0Name(in->numHist*sizeof(ofloat),"Flux Histogram RMS");
  for (i=0; i<in->numHist; i++) out->histRMS[i] = in->histRMS[i];
  out->numHist  = in->numHist;
  out->HistInc  = in->HistInc;

  return out;
} /* end ObitUVSelfCalCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVSelfCal similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVSelfCalClone  (ObitUVSelfCal *in, ObitUVSelfCal *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->modelMode = in->modelMode;
  out->skyModel = ObitSkyModelUnref(out->skyModel);
  out->skyModel = ObitSkyModelRef(in->skyModel);
  out->SCData   = ObitUVUnref(out->SCData);
  out->SCData   = ObitSkyModelRef(in->SCData);
  out->display  = ObitDisplayUnref(out->display);
  out->display  = ObitDisplayRef(in->display);
  out->hist     = ObitMemFree (out->hist);
  out->hist     = ObitMemAlloc0Name(in->numHist*sizeof(ofloat),"Flux Histogram");
  for (i=0; i<in->numHist; i++) out->hist[i] = in->hist[i];
  out->histRMS  = ObitMemFree (out->histRMS);
  out->histRMS  = ObitMemAlloc0Name(in->numHist*sizeof(ofloat),"Flux Histogram RMS");
  for (i=0; i<in->numHist; i++) out->histRMS[i] = in->histRMS[i];
  out->numHist  = in->numHist;
  out->HistInc  = in->HistInc;

} /* end ObitUVSelfCalClone */

/**
 * Creates an ObitUVSelfCal 
 * \param name      An optional name for the object.
 * \param skyModel  Sky model [optional] to normalize uv data
 * \return the new object.
 */
ObitUVSelfCal* ObitUVSelfCalCreate (gchar* name, ObitSkyModel *skyModel)
{
  ObitUVSelfCal* out;

  /* Create basic structure */
  out = newObitUVSelfCal (name);

  /* Save SkyModel */
  if (skyModel!=NULL) out->skyModel = ObitSkyModelRef(skyModel);

  /* Calibration solution object */
  out->mySolver = ObitUVGSolveCreate("Calibration solution");

  return out;
} /* end ObitUVSelfCalCreate */

/**
 * Determine Self calibration for a uv dataset 
 * inUV is divided by the skyModel and an antenna gain solution is made.
 * On output, the dataset has a new SN table with this calibration
 * and is set up to apply them.
 * Routine determines if self calibration is converged, if so TRUE is returned 
 * (else FALSE) and the best SN table is set to be applied.
 * \param in      Input self cal object. 
 * Control parameters are on the info member.
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "solInt"  OBIT_float (1,1,1) Solution interval (min). (default 1 sec)
 * \li "refAnt"  OBIT_int   (1,1,1) Ref ant to use. (default 1)
 * \li "avgPol"  OBIT_bool  (1,1,1) True if RR and LL to be averaged (false)
 * \li "avgIF"   OBIT_bool  (1,1,1) True if all IFs to be averaged (false)
 * \li "SNRMin"  OBIT_float (1,1,1) Minimum acceptable SNR (5)
 * \li "doMGM"   OBIT_bool  (1,1,1) True then find the mean gain modulus (true)
 * \li "solType" OBIT_string (4,1,1 Solution type '  ', 'L1',  (' ')
 * \li "solMode" OBIT_string (4,1,1 Solution mode: 'A&P', 'P', 'P!A', 'GCON' ('P')
 * \li "minNo"   OBIT_int   (1,1,1) Min. no. antennas. (default 4)
 * \li "antWt"   OBIT_float (*,1,1) Antenna weights. (default 1.0)
 * \li "UVR_Full"OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (kilolamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "prtLv"   OBIT_int   (1,1,1) Print level (default no print)
 * \li "minFluxPSC" OBIT_float (1,1,1) min peak flux for phase selfcal               
 * \li "minFluxASC" OBIT_float (1,1,1) min peak flux for A&P selfcal
 * \li "doSmoo"  OBIT_bool  (1,1,1) True then interpolate failed solutions [F]
 * \li "dispURL" OBIT_string scalar = URL of display server
 * \li "peakFlux" OBIT_float scalar If present and > 0.0, then this is the highest
 *                                  image pixel value to use to determine if SC needed
 *                                  else use sum of CC in SkyModel.
 * \param inUV     Input UV data. 
 * \param init     If True, this is the first SC in a series.
 * \param noSCNeed If True, no self calibration was needed
 * \param window   If nonNULL, the CLEAN window to be edited.
 * \param err      Error/message stack, returns if error.
 * \return True if SC converged, else False.  No gain solution if True
 */
gboolean ObitUVSelfCalSelfCal (ObitUVSelfCal *in, ObitUV *inUV, gboolean init, 
			       gboolean *noSCNeed, ObitDConCleanWindow *window,
			       ObitErr *err)
{
  gboolean converged = FALSE;
  ObitTableSN *TabSN=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, jj, nfield, *iatemp, itemp, isuba, refant, bestSN=0, oldDoCal;
  ofloat best, posFlux, peakFlux;
  ofloat minFluxPSC, minFluxASC, minFlux;
  gchar *dispURL=NULL, tname[129], solmod[5], tbuff[128];
  gboolean Tr=TRUE, Fl=FALSE, diverged, quit, doSmoo;
  gchar        *SCParms[] = {  /* Self parameters */
    "refAnt", "solInt", "solType", "solMode", "WtUV", "avgPol", "avgIF", 
    "doMGM", "minSNR", "minNo", "doSmoo", "prtLv", 
    NULL
  };
  gchar *routine = "ObitUVSelfCalSelfCal";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return converged;
  g_assert (ObitUVSelfCalIsA(in));
  g_assert (ObitUVIsA(inUV));

  *noSCNeed = FALSE;
  /* First this series? */
  if (init) {
    in->numLast = 0;
    for (i=0; i<5; i++) {in->lastQual[i] = 0.0; in->lastSNVer[i] = -1;}
  }

  /* Copy gain solution control info */
  ObitInfoListCopyList (in->info, in->mySolver->info, SCParms);

  /* Image display? */
  if (!in->display) {
    ObitInfoListGetP(in->info, "dispURL", &type, dim, (gpointer)&dispURL);
    /* dispURL not NULL terminated */
    if (dispURL) {for (i=0; i<dim[0]; i++) tname[i] = dispURL[i]; tname[i]=0;}
    if (dispURL && (strncmp(tname, "None", 4))) 
      in->display = ObitDisplayCreate("Display", tname, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, converged);
  }

  /* Quality for convergence test of SkyModel */
  /* Total sum of CCs */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut(in->skyModel->info, "noNeg", OBIT_bool, dim, &Fl);
  in->sumCC = ObitSkyModelSum (in->skyModel, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, converged);

  /* Compress CC files */
  ObitSkyModelCompressCC (in->skyModel, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, converged);

  /* Sum of only positive components */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut(in->skyModel->info, "noNeg", OBIT_bool, dim, &Tr);
  posFlux = ObitSkyModelSum (in->skyModel, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, converged);

  /* Self cal needed? */
  peakFlux = posFlux;
  ObitInfoListGetTest(in->info, "peakFlux", &type, dim, &peakFlux);
  /* Peak must be positive */
  if (peakFlux<=0.0) peakFlux = posFlux;
  solmod[0] = solmod[1] = solmod[2] = solmod[3] = '0'; solmod[4] = 0;
  ObitInfoListGetTest(in->info, "solMode", &type, dim, solmod);
  minFluxPSC = 1.0e20;
  ObitInfoListGetTest(in->info, "minFluxPSC", &type, dim, &minFluxPSC);
  if (minFluxPSC==0.0) minFluxPSC = 1.0e20;
  minFluxASC = 1.0e20;
  ObitInfoListGetTest(in->info, "minFluxASC", &type, dim, &minFluxASC);
  if (minFluxASC==0.0) minFluxASC = 1.0e20;
  if (solmod[0] == 'A') minFlux = minFluxASC;
  else minFlux = minFluxPSC;
  if (peakFlux<minFlux) {
    Obit_log_error(err, OBIT_InfoErr, "Peak Flux %g < %g - no Self cal needed", 
		   peakFlux, minFlux);
    *noSCNeed = TRUE;
    return TRUE;
  }
  /* Need model total flux density at least  minFlux) */
  if (posFlux<minFlux) {
    Obit_log_error(err, OBIT_InfoErr, "Model Flux %g < %g - no Self cal needed", 
		   posFlux, minFlux);
    *noSCNeed = TRUE;
    return TRUE;
  }

  /* Other info */
  refant = 0;
  ObitInfoListGetTest(in->info, "refAnt", &type, dim, &refant);
  isuba = 1;
  ObitInfoListGetTest(in->info, "subA", &type, dim, &isuba);
  doSmoo = FALSE;
  ObitInfoListGetTest(in->info, "doSmoo", &type, dim, &doSmoo);
  /* save prior doCal state */
  oldDoCal = -1;
  ObitInfoListGetTest(in->info, "doCalib", &type, dim, &oldDoCal);

  /* Tell user */
  Obit_log_error(err, OBIT_InfoErr, "Self cal with model flux = %g, peak=%g", 
		 posFlux, peakFlux);

  /* RMS field 1 */
  in->RMSFld1 = ObitImageMosaicGetImageRMS (in->skyModel->mosaic, 0, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, converged);

  /* Quality */
  if (in->sumCC>0.0) in->totalQual = in->RMSFld1 / in->sumCC;
  else in->totalQual = 1.0e20;

  /* Divergence test - call it diverged if solution worse that 1.3*last */
  diverged  = (in->numLast>21) && (in->totalQual > 1.3*in->lastQual[in->numLast-1]);
  converged = diverged;

  /* Convergence test - not improved in 2 cycles */
  converged = 
    converged || ((in->numLast>=2) &&
		  (in->totalQual >= in->lastQual[in->numLast-1]) &&
		  (in->lastQual[in->numLast] >= in->lastQual[in->numLast-1]));

  /* Tell user */
  Obit_log_error(err, OBIT_InfoErr, "RMS = %g Jy, RMS/SumCC = %g",
		 in->RMSFld1, in->totalQual);
  /* Tell user last few */
  if (in->numLast>0) {
    jj = 0;
    for (j=0; j<in->numLast; j++) {
      sprintf (&tbuff[jj], " %12.5g",in->lastQual[j]);
      jj = strlen(tbuff);
    }
    Obit_log_error(err, OBIT_InfoErr, "last %d  RMS/SumCC = %s",in->numLast, tbuff);
  }
  if (diverged) Obit_log_error(err, OBIT_InfoErr, "Solution diverging");
  else if (converged) Obit_log_error(err, OBIT_InfoErr, "Solution converged");
  ObitErrLog(err);
 
  /* Show user and allow judgement */
  if (in->display) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Display residual from CLEAN, continue with self-cal?");
    ObitErrLog(err);
    quit = ObitDisplayShow (in->display, (Obit*)in->skyModel->mosaic, window, 
			    1, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, converged);
    if (quit && init) *noSCNeed = TRUE; /* User turned off all self cal */
    if (quit) {converged=TRUE;}
    if (!quit) {  /* Give some sign of life */
      Obit_log_error(err, OBIT_InfoErr, "Continuing");
      ObitErrLog(err);
    }
  }

  if (converged) { /* Converged, find best SN table */
    best = 1.0e20;
    for (i=0; i<in->numLast; i++) {
      if (best<in->lastQual[i]) {
	best = in->lastQual[i];
	bestSN = in->lastSNVer[i];
      }
    }
    
  } else { /* Do self cal if not converged */

    oldDoCal = 2;  /* Apply this SN table later */

    /* Keep track of this result */
    if (in->numLast>=5) {  /* slide window if necessary */
      in->numLast = 4;
      for (i=0; i<4; i++) {
	in->lastQual[i]  = in->lastQual[i+1];
	in->lastSNVer[i] = in->lastSNVer[i+1];
      }
    }
    in->lastQual[in->numLast]  = in->totalQual;
    in->lastSNVer[in->numLast] = 0;
    in->numLast++;
  
    /* Reset SkyModel to use all components */
    nfield = in->skyModel->mosaic->numberImages;
    iatemp = ObitMemAlloc(nfield*sizeof(olong));  /* temp. array */
    dim[0] = nfield;
    for (i=0; i<nfield; i++) iatemp[i] = 1;
    ObitInfoListAlwaysPut(in->skyModel->info, "BComp", OBIT_long, dim, iatemp);
    for (i=0; i<nfield; i++) iatemp[i] = in->skyModel->endComp[i];
    ObitInfoListAlwaysPut(in->skyModel->info, "EComp", OBIT_long, dim, iatemp);
    iatemp = ObitMemFree(iatemp);  /* Deallocate memory */
    
    /* Scratch file to use for self calibration */
    if (in->SCData==NULL) {
      in->SCData = newObitUVScratch (inUV, err);
      if (err->error) Obit_traceback_val (err, routine, inUV->name, converged);
    }
    
    /* Don't apply any  calibration to inUV */
    dim[0] = 1; dim[1] = 1;
    ObitInfoListAlwaysPut(inUV->info, "doCalSelect", OBIT_bool, dim, &Fl);
    itemp = -1;
    ObitInfoListAlwaysPut(inUV->info, "doCalib", OBIT_long, dim, &itemp);
    
    /* Divide data by model */
    ObitSkyModelDivUV (in->skyModel, inUV, in->SCData, err);
    if (err->error) Obit_traceback_val (err, routine, in->skyModel->name, converged);
    
    /* Self Calibrate */
    in->sumCC = posFlux;  /* Model flux actually being used */
    in->mySolver->UVFullRange[0] = in->UVFullRange[0];
    in->mySolver->UVFullRange[1] = in->UVFullRange[1];
    /* Force new solution table */
    itemp = 0;
    ObitInfoListAlwaysPut(in->mySolver->info, "solnVer", OBIT_long, dim, &itemp);
    TabSN = ObitUVGSolveCal (in->mySolver, in->SCData, inUV, inUV->mySel, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, converged);

    /* Get SN version number */
    in->lastSNVer[in->numLast-1] = ObitTableGetVersion ((ObitTable*)TabSN, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, converged);

    /* Make sure all references to same antenna */
    ObitUVSolnRefAnt(TabSN, isuba, &refant, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, converged);

    /* interpolate failed solutions if requested */
    if (doSmoo) {
      dim[0] = 1; dim[1] = 1;
      ObitInfoListAlwaysPut(TabSN->info, "doBlank", OBIT_bool, dim, &Tr);
      ObitUVSolnSNSmo (TabSN, isuba, err);
    }
    if (err->error) Obit_traceback_val (err, routine, in->name, converged);
    
    TabSN =  ObitTableSNUnref (TabSN);
  } /* end of self calibrate if not converged */

  /* Apply Calibration */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut(inUV->info, "doCalSelect", OBIT_bool, dim, &Tr);
  itemp = oldDoCal;  /* Restore prior if no self cal here */
  ObitInfoListAlwaysPut(inUV->info, "doCalib", OBIT_long, dim, &itemp);
  itemp = bestSN;
  ObitInfoListAlwaysPut(inUV->info, "gainUse", OBIT_long, dim, &itemp);
  /* Use all components */
  ObitInfoListAlwaysPut(in->skyModel->info, "noNeg", OBIT_bool, dim, &Fl);

  return converged;
} /* end ObitUVSelfCalSelfCal */

/**
 * Determine Self calibration for a uv dataset from a specified model
 * inUV is divided by the point model and an antenna gain solution is made.
 * On output, the dataset has a new SN table with this calibration
 * and is set up to apply them.
 * \param in      Input self cal object. 
 * Control parameters are on the info member.
 * \li "modelFlux"  OBIT_float (1,1,1) Model flux density [def 1.0]
 * \li "modelPos"   OBIT_float (2,1,1) Model posn. offset (X, y) asec [def 0]
 * \li "modelParm"  OBIT_float (?,1,1) Other model parameters:
 *      major_axis (deg),  minor_axis (deg),  position_angle (deg),
 *      type (ObitSkyModelCompType as gint); [def 0]
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "solInt"  OBIT_float (1,1,1) Solution interval (min). (default 1 sec)
 * \li "refAnt"  OBIT_int   (1,1,1) Ref ant to use. (default 1)
 * \li "avgPol"  OBIT_bool  (1,1,1) True if RR and LL to be averaged (false)
 * \li "avgIF"   OBIT_bool  (1,1,1) True if all IFs to be averaged (false)
 * \li "SNRMin"  OBIT_float (1,1,1) Minimum acceptable SNR (5)
 * \li "doMGM"   OBIT_bool  (1,1,1) True then find the mean gain modulus (true)
 * \li "solType" OBIT_string (4,1,1 Solution type '  ', 'L1',  (' ')
 * \li "solMode" OBIT_string (4,1,1 Solution mode: 'A&P', 'P', 'P!A', 'GCON' ('P')
 * \li "minNo"   OBIT_int   (1,1,1) Min. no. antennas. (default 4)
 * \li "antWt"   OBIT_float (*,1,1) Antenna weights. (default 1.0)
 * \li "UVR_Full"OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (kilolamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "prtLv"   OBIT_int   (1,1,1) Print level (default no print)
 * \li "doSmoo"  OBIT_bool  (1,1,1) True then interpolate failed solutions [F]
 * \param inUV     Input UV data. 
 * \param err      Error/message stack, returns if error.
 */
void ObitUVSelfCalModel (ObitUVSelfCal *in, ObitUV *inUV, ObitErr *err)
{
  ObitTableSN *TabSN=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong itemp, isuba, refant, nparm, bestSN=0;
  ofloat pointFlux, pointPos[2], pointParms[20];
  gchar solmod[5];
  gboolean Tr=TRUE, Fl=FALSE, doSmoo;
  gchar        *SCParms[] = {  /* Self parameters */
    "refAnt", "solInt", "solType", "solMode", "WtUV", "avgPol", "avgIF", 
    "doMGM", "minSNR", "minNo", "doSmoo", "prtLv", 
    NULL
  };
  gchar *routine = "ObitUVSelfCalModel";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVSelfCalIsA(in));
  g_assert (ObitUVIsA(inUV));

  /* Control parameters? */
  pointFlux = 1.0;
  ObitInfoListGetTest(in->info, "modelFlux", &type, dim, &pointFlux);
  pointPos[0] = pointPos[1] = 0.0;
  ObitInfoListGetTest(in->info, "modelPos", &type, dim, &pointPos);
  pointParms[0] = pointParms[3] = 0.0; dim[0] = 1;
  ObitInfoListGetTest(in->info, "modelParm", &type, dim, &pointParms);
  nparm = dim[0];
  solmod[0] = solmod[1] = solmod[2] = solmod[3] = '0'; solmod[4] = 0;
  ObitInfoListGetTest(in->info, "solMode", &type, dim, solmod);
  refant = 0;
  ObitInfoListGetTest(in->info, "refAnt", &type, dim, &refant);
  isuba = 1;
  ObitInfoListGetTest(in->info, "subA", &type, dim, &isuba);
  doSmoo = FALSE;
  ObitInfoListGetTest(in->info, "doSmoo", &type, dim, &doSmoo);

  /* Tell user */
  Obit_log_error(err, OBIT_InfoErr, "Self cal with point model ");

  /* Scratch file to use for self calibration */
  if (in->SCData==NULL) {
    in->SCData = newObitUVScratch (inUV, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  }
  
  /* Copy gain solution control info */
  ObitInfoListCopyList (in->info, in->mySolver->info, SCParms);

  /* Don't apply any  calibration to inUV */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut(inUV->info, "doCalSelect", OBIT_bool, dim, &Fl);
  itemp = -1;
  ObitInfoListAlwaysPut(inUV->info, "doCalib", OBIT_long, dim, &itemp);
  
  /* Model to correct units */
  pointPos[0] /= 3600.0;
  pointPos[1] /= 3600.0;
  if (pointParms[3]==1.0) { /* Gaussian Size to degrees */
    pointParms[1] /= 3600.0;  
    pointParms[2] /= 3600.0;
  } else if (pointParms[3]==2.0) { /* Unif. sphere Size to degrees */
    pointParms[1] /= 3600.0;  
  }

  /* Set Model on SkyModel */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(in->skyModel->info, "MODPTFLX", OBIT_float, dim, &pointFlux);
  ObitInfoListAlwaysPut(in->skyModel->info, "MODPTXOF", OBIT_float, dim, &pointPos[0]);
  ObitInfoListAlwaysPut(in->skyModel->info, "MODPTYOF", OBIT_float, dim, &pointPos[1]);
  dim[0] = nparm;
  ObitInfoListAlwaysPut(in->skyModel->info, "MODPTYPM", OBIT_float, dim, &pointParms);

  /* Divide data by model */
  ObitSkyModelDivUV (in->skyModel, inUV, in->SCData, err);
  if (err->error) Obit_traceback_msg (err, routine, in->skyModel->name);
  
  /* Self Calibrate */
  in->sumCC = pointFlux;  /* Model flux actually being used */
  in->mySolver->UVFullRange[0] = in->UVFullRange[0];
  in->mySolver->UVFullRange[1] = in->UVFullRange[1];
  TabSN =  ObitUVGSolveCal (in->mySolver, in->SCData, inUV, inUV->mySel, err);
  if (err->error) {
    TabSN =  ObitTableSNUnref (TabSN);
    Obit_traceback_msg (err, routine, in->name);
  }
  
  /* Get SN version number */
  in->lastSNVer[in->numLast-1] = ObitTableGetVersion ((ObitTable*)TabSN, err);
  if (err->error) {
    TabSN =  ObitTableSNUnref (TabSN);
    Obit_traceback_msg (err, routine, in->name);
  }
  
  /* Make sure all references to same antenna */
  ObitUVSolnRefAnt(TabSN, isuba, &refant, err);
  if (err->error) {
    TabSN =  ObitTableSNUnref (TabSN);
    Obit_traceback_msg (err, routine, in->name);
  }
  
  /* interpolate failed solutions if requested */
  if (doSmoo) ObitUVSolnSNSmo (TabSN, isuba, err);
  if (err->error) {
    TabSN =  ObitTableSNUnref (TabSN);
    Obit_traceback_msg (err, routine, in->name);
  }
  
  TabSN =  ObitTableSNUnref (TabSN);

  /* Clear Model on SkyModel */
  dim[0] = dim[1] = dim[2] = 1; pointFlux = 0.0;
  ObitInfoListAlwaysPut(in->skyModel->info, "MODPTFLX", OBIT_float, dim, &pointFlux);
  in->skyModel->pointFlux = 0.0;

  /* Apply Calibration next time */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut(inUV->info, "doCalSelect", OBIT_bool, dim, &Tr);
  itemp = 2;
  ObitInfoListAlwaysPut(inUV->info, "doCalib", OBIT_long, dim, &itemp);
  itemp = bestSN;
  ObitInfoListAlwaysPut(inUV->info, "gainUse", OBIT_long, dim, &itemp);
  /* Use all components */
  ObitInfoListAlwaysPut(in->skyModel->info, "noNeg", OBIT_bool, dim, &Fl);

} /* end ObitUVSelfCalModel */

/**
 * Determine the flux density histogram of amp vs baseline length.
 * Histogram forced to be monitonically decreasing with baseline.
 * Also computes RMS of average values as histRMS
 * \param in    Input Self Cal
 * \param inUV  UV data to examine
 * \param err   Error stack, returns if not empty.
 */
void ObitUVSelfCalFluxHist (ObitUVSelfCal *in, ObitUV *inUV, ObitErr *err)
{
  ObitIOCode retCode=OBIT_IO_SpecErr;
  olong i, j, k, inc, icell;
  ofloat MaxBL, MaxW, *count, bl, *u, *v, *w, *vis, icells, last, amp;
  ofloat arg, sum=0.0, sum2=0.0, mavg, mrms;
  olong cnt;
  gchar *routine = "ObitUVSelfCalFluxHist";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVSelfCalIsA(in));
  g_assert (ObitUVIsA(inUV));

  /* Get baseline range */
  ObitUVUtilUVWExtrema (inUV, &MaxBL, &MaxW, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Create histogram, count arrays */
  /* Size of histogram Max (1000, 1% of vis) */
  in->numHist = MIN (1000, (olong)(0.01*inUV->myDesc->nvis)); 
  in->HistInc = MaxBL / in->numHist;  /* Increment between cells */
  if (in->hist) in->hist    = ObitMemFree(in->hist); /* Free old if exists */
  in->hist    = ObitMemAlloc0Name(in->numHist*sizeof(ofloat),"Flux Histogram");
  if (in->histRMS) in->histRMS = ObitMemFree(in->histRMS); /* Free old if exists */
  in->histRMS = ObitMemAlloc0Name(in->numHist*sizeof(ofloat),"Flux Histogram");
  count       = ObitMemAlloc0Name(in->numHist*sizeof(ofloat),"Histogram Counts");
  for (i=0; i<in->numHist; i++) {in->hist[i] = 0.0; in->histRMS[i]=0.0; count[i]=0.0;}
  icells = 1.0 / in->HistInc;  /* inverse of cell spacing */

  /* Open uv data if not already open */
  if (inUV->myStatus==OBIT_Inactive) {
    retCode = ObitUVOpen (inUV, OBIT_IO_ReadOnly, err);
    if (err->error) goto cleanup;
  }
  
  /* increment in frequency array = number of Stokes */
  inc = inUV->myDesc->inaxes[inUV->myDesc->jlocs];  

  /* Loop through data */
  while (retCode==OBIT_IO_OK) {
    /* read buffer full */
    retCode = ObitUVRead (inUV, NULL, err);
    if (err->error) goto cleanup;
    
    /* initialize data pointers */
    u   = inUV->buffer+inUV->myDesc->ilocu;
    v   = inUV->buffer+inUV->myDesc->ilocv;
    w   = inUV->buffer+inUV->myDesc->ilocw;
    vis = inUV->buffer+inUV->myDesc->nrparm;
    for (i=0; i<inUV->myDesc->numVisBuff; i++) { /* loop over buffer */
      
      /* Get statistics */
      bl = sqrt ((*u)*(*u) + (*v)*(*v));

      /* Loop over correlations */
      for (j=0; j<inUV->myDesc->ncorr; j++) {
	k = j / inc;  /* Index in frequency scaling array */
	if (vis[j*3+2]>0.0) {
	  icell = bl * icells * inUV->myDesc->fscale[k] + 0.5;
	  icell = MAX (0, MIN (icell, in->numHist-1));
	  amp = sqrt (vis[j*3]*vis[j*3] + vis[j*3+1]*vis[j*3+1]);
	  count[icell]++;
	  in->hist[icell]    += amp;
	  in->histRMS[icell] += amp*amp;
	} /* end if valid data */
      } /* end loop over correlations */
      
      /* update data pointers */
      u += inUV->myDesc->lrec;
      v += inUV->myDesc->lrec;
      w += inUV->myDesc->lrec;
      vis += inUV->myDesc->lrec;
    } /* end loop over buffer */
  } /* end loop over file */
  
    /* Close */
  retCode = ObitUVClose (inUV, err);
  if (err->error) goto cleanup;

  /* Normalize */
  last = 1.0e20;
  for (i=0; i<in->numHist; i++) {
    if (count[i] > 3.0) {
      sum  = in->hist[i];
      sum2 = in->histRMS[i];
      in->hist[i] /= count[i];
      last = in->hist[i];
      /* Get RMS */
      if (count[i]>9.0) {
	arg = (sum2 - ((sum*sum)/count[i])) / (count[i]-1.0);
	in->histRMS[i] = sqrt (fabs(arg));
      } else {
	/* If too little data use average */
	in->histRMS[i] = in->hist[i];	
      }
    }  else in->hist[i] = last;
  }

  /* Average histogram RMS of cells with at least 10 counts */
  sum = 0.0; sum2 = 0.0; cnt = 0;
  for (i=1; i<in->numHist; i++) {
    if (count[i]>=10.0) {
      cnt++;
      sum  += in->histRMS[i];
      sum2 += in->hist[i];
    }
  }
  mrms = sum/cnt;
  mavg = sum2/cnt;  /* In case not enough in inner 10 % */

  /* Average histogram values in first 10 % of histogram */
  sum = 0.0; cnt = 0; i=1;
  while (i<0.1*in->numHist) {
    if (count[i]>=10.0) {
      cnt++;
      sum += in->hist[i];
    }
    i++;
  }
  if (cnt>1)  mavg = sum/cnt;

  /*  Reset RMSes to no less than the average and amplitudes to no 
      more than the average of the first 10 % of the uv range. */

  for (i=0; i<in->numHist; i++) {
    in->hist[i]    = MIN (in->hist[i],    mavg);
    in->histRMS[i] = MAX (in->histRMS[i], mrms);
  }

  /* Cleanup */
  cleanup: count  = ObitMemFree(count);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
} /* end  ObitUVSelfCalFluxHist */

/**
 * Find UV range that encloses a given flux density..
 * Returns a uv range from the longest baseline whose average flux 
 * amplitude - 1 X RMS exceeds the sum of the clean components in the 
 * mosaic member of the SkyModel member to infinity.
 * Results saved in the UVFullRange member.
 * If there is no histogram all baselines are accepted.
 * \param in       Input Self Cal
 * \param inUV     UV data to examine
 * \param err      Error stack, returns if not empty.
 */
void ObitUVSelfCalBLRange (ObitUVSelfCal *in, ObitErr *err)
{
  olong i, last;
  /*gchar *routine = "ObitUVSelfCalBLRange";*/

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVSelfCalIsA(in));

  /* If no histogram take all */
  if (in->numHist<=0) {
    in->UVFullRange[0] = 0.0;
    in->UVFullRange[1] = 1.0e20;
    return;
  }

  /* Find first element in hist < in->sumCC */
  last = 0;
  for (i=0; i<in->numHist; i++) {
    if ((in->hist[i]-in->histRMS[i]) > in->sumCC) last = i;
  }
  last = MIN (last, in->numHist-1);

  /* Save result */
  in->UVFullRange[0] = last * in->HistInc;
  in->UVFullRange[1] = 1.0e20;

  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, "Full weight  from %g lambda in SelfCal", 
		 in->UVFullRange[0]*0.001);

} /* end ObitUVSelfCalBLRange */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVSelfCalClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVSelfCalClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVSelfCalClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVSelfCalClassInfoDefFn (gpointer inClass)
{
  ObitUVSelfCalClassInfo *theClass = (ObitUVSelfCalClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVSelfCalClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVSelfCalClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVSelfCalGetClass;
  theClass->newObit       = (newObitFP)newObitUVSelfCal;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVSelfCalCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVSelfCalClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVSelfCalInit;
  theClass->ObitUVSelfCalCreate = (ObitUVSelfCalCreateFP)ObitUVSelfCalCreate;

} /* end ObitUVSelfCalClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVSelfCalInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVSelfCal *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList(); 
  in->modelMode = OBIT_SkyModel_Fastest;
  in->mySolver  = NULL;
  in->skyModel  = NULL;
  in->display   = NULL;
  in->SCData    = NULL;
  in->numHist   = 0;
  in->HistInc   = 1.0;
  in->hist      = NULL;
  in->UVFullRange[0] = 0.0;
  in->UVFullRange[1] = 0.0;

} /* end ObitUVSelfCalInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVSelfCal* cast to an Obit*.
 */
void ObitUVSelfCalClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVSelfCal *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread   = ObitThreadUnref(in->thread);
  in->info     = ObitInfoListUnref(in->info);
  in->mySolver = ObitUVGSolveUnref(in->mySolver);
  in->display  = ObitDisplayUnref(in->display);
  in->skyModel = ObitSkyModelUnref(in->skyModel);
  in->SCData   = ObitUVUnref(in->SCData);
  in->hist     = ObitMemFree(in->hist);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVSelfCalClear */


