/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2010                                          */
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

#include "ObitUVSoln.h"
#include "ObitTableUtil.h"
#include "ObitTableSUUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSoln.c
 * ObitUVSoln class function definitions.
 * This clas allows manipulation (interpolation) of Solutions in an ObitTableSN.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVSoln";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVSolnClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVSolnClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVSolnInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVSolnClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVSolnClassInfoDefFn (gpointer inClass);

/** Private:  Read calibration for a new time into the internal arrays. */
static void ObitUVSolnNewTime (ObitUVSoln *in, ofloat time,
			       ObitErr *err);

/** Private:  Determine number of times and usage of each antenna as reference */
static olong 
*refCount (ObitTableSN *SNTab, olong isub, olong *numtime, ObitErr* err);

/** Private: Refer phases to a common reference antenna  */
static void 
refPhase (ObitTableSN *SNTab, olong isub, olong iif, olong refant, olong ant, 
	  olong mxtime, ofloat* wrktim, 
	  ofloat* work1, ofloat* work2, ofloat* work3, ofloat* work4, 
	  ofloat* work5, ObitErr* err);

/** Private: Refer MB delays to a common reference antenna  */
static void 
refMBDelay (ObitTableSN *SNTab, olong isub, olong refant, olong ant, 
	    olong mxtime, ofloat* wrktim, 
	    ofloat* work1, ofloat* work2, ofloat* work4, ObitErr* err);

/** Private: Refer delays to a common reference antenna  */
static void 
refDelay (ObitTableSN *SNTab, olong isub, olong iif, olong refant, olong ant, 
	  olong mxtime, ofloat* wrktim, 
	  ofloat* work1, ofloat* work2, ofloat* work4, ObitErr* err);

/** Private: Refer rates to a common reference antenna  */
static void 
refRate (ObitTableSN *SNTab, olong isub, olong iif, olong refant, olong ant, 
	 olong mxtime, ofloat* wrktim, 
	 ofloat* work1, ofloat* work2, ofloat* work4, ObitErr* err);

/** Private: Boxcar smoothing of an irregularly spaced array  */
static void 
boxsmo (ofloat width, ofloat* x, ofloat* y, olong n, ofloat* ys);

/** Private: Average rates in an open SN table */
static void avgRate (ObitTableSN *SNTab, ObitErr* err);

/** Private: Average prior rates for Antenna and poln  over IF */
static ofloat avgRatePrior (ObitUVSoln *in, olong ant, olong ipol);

/** Private: Average following rates for Antenna and poln  over IF */
static ofloat avgRateFollow (ObitUVSoln *in, olong ant, olong ipol);

/** Private: Make selector for calibrators */
static ObitUVSel* MakeCalSelect (ObitUV *in, ObitInfoList *info, ObitErr *err);

/** Coherent amplitude/phase smoothing */
void smoAmpPh (ObitTableSN *SNTab, ObitUVSel *sel, gchar* smoFunc, gchar* smoType, 
	       ofloat alpha, ofloat *smoParm,
	       olong sub, olong ifbeg, olong ifend, odouble* freqs, gboolean doBlank, 
	       ofloat* gncnt, ofloat* gnsum, ObitErr* err);

/** Private: Generic smoothing */
void smoIt (gchar* smmeth, ofloat width, ofloat alpha, 
	    ofloat* x, ofloat *t, ofloat *w, olong n, 
	    ofloat* xs, ofloat* ws, ofloat *wrk1, ofloat *wrk2, gboolean doBlank);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVSoln* newObitUVSoln (gchar* name)
{
  ObitUVSoln* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVSolnClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVSoln));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVSolnInit((gpointer)out);

 return out;
} /* end newObitUVSoln */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVSolnGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVSolnClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVSolnGetClass */

/**
 * Creates an ObitUVSoln 
 * \param name  An optional name for the object.
 * \param inUV  UV data with input SN table
 * \return the new object.
 */
ObitUVSoln* ObitUVSolnCreate (gchar* name, ObitUV *inUV)
{
  ObitUVSoln* out;

  /* Create basic structure */
  out = newObitUVSoln (name);

  /* save input uv data pointer */
  out->myUV = ObitUVRef(inUV);

  return out;
} /* end ObitUVSolnCreate */

/**
 * Initialize structures solution interpolation .
 * \param in   Solution Object.
 * \param err  ObitError stack.
 */
void ObitUVSolnStartUp (ObitUVSoln *in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitUVSel  *sel  = in->myUV->mySel;
  ObitUVDesc *desc = in->myUV->myDesc;
  ObitTableSN *newTableSN=NULL;
  ObitInfoType type;
  gboolean doBlank;
  gint32 dim[MAXINFOELEMDIM];
  olong size, i, iif, iant, SNver, highVer;
  olong itemp;
  gchar interMode[9];
  gchar *routine="ObitUVSolnStartUp";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVSolnIsA(in));

  /* Get calibrator selector */
  in->CalSel = MakeCalSelect (in->myUV, in->info, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Get parameters */
  strcpy (interMode, "    ");
  ObitInfoListGetTest(in->info, "interMode", &type, dim, interMode);
  in->interNPoly = 2;
  ObitInfoListGetTest(in->info, "interNPoly", &type, dim, &in->interNPoly);
  in->maxInter = 1440.0; /* default 1 day */
  ObitInfoListGetTest(in->info, "maxInter", &type, dim, &in->maxInter);
  if (in->maxInter<=0.0) in->maxInter = 1.0e20;
  in->maxInter /= 1440.0;  /* to days */
  for (i=0; i<10; i++) in->interParm[i] = 0.0;
  ObitInfoListGetTest(in->info, "interParm", &type, dim, &in->interParm);
  itemp = 0;
  ObitInfoListGetTest(in->info, "solnVer", &type, dim, &itemp);
  SNver = itemp;
  doBlank = FALSE;
  ObitInfoListGetTest(in->info, "doBlank", &type, dim, &doBlank);


  /* Which interpolation mode */
  in->interMode = OBIT_UVSolnInterUnknown;
  if (!strncmp(interMode, "    ",4) || !strncmp(interMode, "2PT",3)) {
    in->interMode = OBIT_UVSolnInter2PT;
  } else if (!strncmp(interMode, "SELF",4)) {
    in->interMode = OBIT_UVSolnInterSELF;
  } else if (!strncmp(interMode, "POLY",4)) {
    in->interMode = OBIT_UVSolnInterPOLY;
    /* NYI */
    Obit_log_error(err, OBIT_Error,"%s: Unsupported interpolation type %s",
		   routine, interMode);
    return;
  } else if (!strncmp(interMode, "SIMP",4)) {
    in->interMode = OBIT_UVSolnInterSIMP;
  } else if (!strncmp(interMode, "AMBG",4)) {
    in->interMode = OBIT_UVSolnInterAMBG;
  } else if (!strncmp(interMode, "CUBE",4)) {
    in->interMode = OBIT_UVSolnInterCUBE;
  } else if (!strncmp(interMode, "MWF",3)) {
    in->interMode = OBIT_UVSolnInterMWF;
  } else if (!strncmp(interMode, "BOX",3)) {
    in->interMode = OBIT_UVSolnInterBOX;
  } else if (!strncmp(interMode, "GAUS",4)) {
    in->interMode = OBIT_UVSolnInterGAUS;
  } else { /* Unknown */
    Obit_log_error(err, OBIT_Error,"%s: Unknown interpolation type %s",
		   routine, interMode);
    return;
  }

  /* Open and close UV data to set selection */
  ObitUVOpen(in->myUV, OBIT_IO_ReadCal, err);
  ObitUVClose(in->myUV, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Copy Selector information */
  in->SubA        = sel->SubA;
  in->FreqID      = sel->FreqID;

  /* Copy descriptor information */
  in->numAnt    = desc->maxAnt;
  in->numSubA   = desc->numSubA;
  in->DeltaTime = desc->DeltaTime;
  in->numIF     = desc->inaxes[desc->jlocif];

  /* Nothing read yet */
  in->LastRowRead = 0;

  /* In case restarting */
  if (in->PriorAntTime)  g_free (in->PriorAntTime);  in->PriorAntTime  = NULL;
  if (in->FollowAntTime) g_free (in->FollowAntTime); in->FollowAntTime = NULL;
  if (in->CalApply)      g_free (in->CalApply);      in->CalApply      = NULL;
  if (in->CalPrior)      g_free (in->CalPrior);      in->CalPrior      = NULL;
  if (in->CalFollow)     g_free (in->CalFollow);     in->CalFollow     = NULL;
  if (in->IFR)           g_free (in->IFR);           in->IFR           = NULL;
  if (in->PriorIFR)      g_free (in->PriorIFR);      in->PriorIFR      = NULL;
  if (in->FollowIFR)     g_free (in->FollowIFR);     in->FollowIFR     = NULL;
  if (in->MBDelay)       g_free (in->MBDelay);       in->MBDelay       = NULL;
  if (in->PriorMBDelay)  g_free (in->PriorMBDelay);  in->PriorMBDelay  = NULL;
  if (in->FollowMBDelay) g_free (in->FollowMBDelay); in->FollowMBDelay = NULL;
  if (in->RefAnt)        g_free (in->RefAnt);        in->RefAnt        = NULL;
  if (in->RateFact)      g_free (in->RateFact);      in->RateFact      = NULL;
  if (in->MissAnt)       g_free (in->MissAnt);       in->MissAnt       = NULL;

  /* Which SN table? */
  highVer = ObitTableListGetHigh (in->myUV->tableList, "AIPS SN");
  if (SNver==0) SNver = highVer;
  in->SNTable =
    newObitTableSNValue (in->name, (ObitData*)in->myUV, &SNver, 
			 OBIT_IO_ReadOnly, 0, 0, err);
  Obit_return_if_fail((in->SNTable!=NULL), err,
		      "%s: Cannot find AIPS SN table", routine);
  
  /* Open solution table, get numPol  */
  retCode = ObitTableSNOpen (in->SNTable, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  in->numPol = in->SNTable->numPol;
  in->numRow = in->SNTable->myDesc->nrow;

  /* Smoothing needed? */
  in->isSNSmoo = FALSE;
  if ((in->interMode == OBIT_UVSolnInterMWF) || 
      (in->interMode == OBIT_UVSolnInterBOX) || 
      (in->interMode == OBIT_UVSolnInterGAUS)) {

    /* Close SN table */
    ObitTableSNClose (in->SNTable, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Create new, temporary SN table */
    SNver = highVer+1;
    newTableSN = 
      newObitTableSNValue (in->name, (ObitData*)in->myUV, &SNver, 
			   OBIT_IO_WriteOnly, in->SNTable->numPol, 
			   in->SNTable->numIF, err);
    newTableSN = ObitTableSNCopy (in->SNTable, newTableSN, err);
    in->SNTable = ObitTableSNUnref(in->SNTable);
    in->SNTable = newTableSN;
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    in->isSNSmoo = TRUE;  /* Remember and zap in shutdown */

    /* Smooth table */
    dim[0] = 4; dim[1] = 1;
    ObitInfoListAlwaysPut(in->SNTable->info, "smoType", OBIT_string, dim, interMode);
    dim[0] = 1; dim[1] = 1;
    ObitInfoListAlwaysPut(in->SNTable->info, "smoAmp",   OBIT_float, dim, &in->interParm[0]);
    ObitInfoListAlwaysPut(in->SNTable->info, "smoPhase", OBIT_float, dim, &in->interParm[1]);
    ObitInfoListAlwaysPut(in->SNTable->info, "doBlank", OBIT_bool, dim, &doBlank);
    ObitUVSolnSNSmo (in->SNTable, in->SubA, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Open smoothed table */
    retCode = ObitTableSNOpen (in->SNTable, OBIT_IO_ReadOnly, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    in->interMode = OBIT_UVSolnInter2PT; /* Now straight interpolate */
  }  /* end smooth */

  /* row to read SN table */  
  in->SNTableRow = newObitTableSNRow((ObitTableSN*)(in->SNTable));

  /* Allocate calibration arrays */
  in->lenCalArrayEntry = 6; /* length of cal array entry */
  size = in->numAnt * in->numIF * in->numPol * in->lenCalArrayEntry;
  in->CalApply     = g_malloc0(size*sizeof(ofloat));
  in->CalPrior     = g_malloc0(size*sizeof(ofloat));
  in->CalFollow    = g_malloc0(size*sizeof(ofloat));
  in->IFR          = g_malloc0(in->numAnt*sizeof(ofloat));
  in->PriorIFR     = g_malloc0(in->numAnt*sizeof(ofloat));
  in->FollowIFR    = g_malloc0(in->numAnt*sizeof(ofloat));
  in->MBDelay      = g_malloc0(2*in->numAnt*sizeof(ofloat));
  in->PriorMBDelay = g_malloc0(2*in->numAnt*sizeof(ofloat));
  in->FollowMBDelay= g_malloc0(2*in->numAnt*sizeof(ofloat));
  in->PriorAntTime = g_malloc0(in->numAnt*sizeof(ofloat));
  in->FollowAntTime= g_malloc0(in->numAnt*sizeof(ofloat));
  in->RefAnt       = g_malloc0(2*in->numAnt*sizeof(olong));
  in->RateFact     = g_malloc0(in->numIF*sizeof(ofloat));
  in->MissAnt      = g_malloc0(in->numAnt*sizeof(gboolean));

  /* Initial times to trigger update of calibration arrays */
  in->CalTime       = -1.0e20;
  in->PriorCalTime  = -1.0e20;
  in->FollowCalTime = -1.0e20;

  /* Fill frequency related arrays */
  for (iif=0; iif<in->numIF; iif++) {
    in->RateFact[iif] = 2.0 * G_PI * 86400.0 * desc->freqIF[iif];
  }

  /* No antennas known to be missing yet */
  for (iant=0; iant<in->numAnt; iant++) in->MissAnt[iant] = FALSE;

} /*  end ObitUVSolnStartUp */

/**
 * Interpolate calibration for a given time, antenna, source...
 * Input values on an ObitTableSNRow into which the results are added.
 * Data selection is enabled.
 * \param in    Solution Object.
 * \param SNRow SN row giving desired time, antenna..., and to receive output
 * \param err   ObitError stack.
 * \return TRUE if source, time etc. selected, else FALSE.
 */
gboolean ObitUVSolnGetSN (ObitUVSoln *in, ObitTableSNRow *SNrow, ObitErr *err)
{
  gboolean out = FALSE;
  gboolean   wanted;
  olong asize, jndxa, indxa, ant, numAnt, iif, ioff, ipol;
  ofloat fblank = ObitMagicF();
  ObitUVDesc *desc;
  ObitUVSel *sel;
  gchar *routine="ObitUVSolnGetSN ";

  /* error checks */
  if (err->error) return out;

  /* local pointers for structures */
  desc = in->myUV->myDesc;
  sel  = in->myUV->mySel;

  /* Check if this data wanted */
  wanted = (SNrow->Time>=sel->timeRange[0]) && (SNrow->Time<=sel->timeRange[1]);
  wanted = wanted && ObitUVSelWantSour (sel, SNrow->SourID);
  wanted = wanted && ObitUVSelWantAnt (sel, SNrow->antNo);
  wanted = wanted && ((SNrow->FreqID==sel->FreqID) || (sel->FreqID<=0));
  wanted = wanted && ((SNrow->SubA==sel->SubA) || (sel->SubA<=0));
  if (!wanted) return wanted;
  out = wanted;

  /* see if new time - update tables. */
  if (SNrow->Time > in->CalTime) {
    ObitUVSolnUpdate (in, SNrow->Time, SNrow->SourID, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
  }

  /* set antenna index */
  asize = in->numPol * in->numIF * in->lenCalArrayEntry;
  jndxa = (SNrow->antNo - 1) * asize;
  ant = SNrow->antNo - 1;
  numAnt = in->numAnt;

  /* loop over IF */
  for (iif=0; iif<in->numIF; iif++) { /* loop 300 */
    ioff = (iif-1) * desc->incif;
    
    /* loop over polarization */
    for (ipol= 1; ipol<=MIN (2, in->numPol); ipol++) { /* loop 200 */
      indxa = jndxa + (ipol-1) * in->lenCalArrayEntry;

      /* if cal Bad - blank */
      if (in->CalApply[indxa]==fblank)  {
	if (ipol==1) {  /* First pol */
	  SNrow->IFR = 0.0;
	  SNrow->MBDelay1     = fblank;
	  SNrow->Real1[iif]   = fblank;
	  SNrow->Imag1[iif]   = fblank;
	  SNrow->Delay1[iif]  = fblank;
	  SNrow->Rate1[iif]   = fblank;
	  SNrow->Weight1[iif] = 0.0;
	  SNrow->RefAnt1[iif] = 0;
	} else {       /* second pol */
	  SNrow->MBDelay2     = fblank;
	  SNrow->Real2[iif]   = fblank;
	  SNrow->Imag2[iif]   = fblank;
	  SNrow->Delay2[iif]  = fblank;
	  SNrow->Rate2[iif]   = fblank;
	  SNrow->Weight2[iif] = 0.0;
	  SNrow->RefAnt2[iif] = 0;
	}
	continue;
      }

      /* Set SN row values */
      if (ipol==1) {  /* First pol */
	SNrow->IFR          = in->IFR[ant];
	SNrow->MBDelay1     = in->MBDelay[ant];
	SNrow->Real1[iif]   = in->CalApply[indxa];
	SNrow->Imag1[iif]   = in->CalApply[indxa+1];
	SNrow->Delay1[iif]  = in->CalApply[indxa+2];
	SNrow->Rate1[iif]   = in->CalApply[indxa+3];
	SNrow->Weight1[iif] = in->CalApply[indxa+4];
	SNrow->RefAnt1[iif] = in->RefAnt[ant];
      } else {       /* second pol */
	SNrow->MBDelay2     = in->MBDelay[ant+numAnt];
	SNrow->Real2[iif]   = in->CalApply[indxa];
	SNrow->Imag2[iif]   = in->CalApply[indxa+1];
	SNrow->Delay2[iif]  = in->CalApply[indxa+2];
	SNrow->Rate2[iif]   = in->CalApply[indxa+3];
	SNrow->Weight2[iif] = in->CalApply[indxa+4];
	SNrow->RefAnt2[iif] = in->RefAnt[ant+numAnt];
      }

    } /* end loop loop over Stokes  L200 */;

    /* setup for next IF */
    jndxa += in->lenCalArrayEntry * in->numPol;
  } /* end loop  L300 - loop over IF */;

  return out;
} /* end ObitUVSolnGetSN */


/**
 * Shutdown Solution interpolation.
 * Close any open file and destroy structures.
 * \param in   Calibrate Object.
 * \param err  ObitError stack.
 */
void ObitUVSolnShutDown (ObitUVSoln *in, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *tabtype=NULL;
  olong lver;
  gchar *routine="ObitUVSolnShutDown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVSolnIsA(in));

  /* Close calibration table, release row structure  */
  retCode = ObitTableSNClose ((ObitTableSN*)in->SNTable, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  in->SNTableRow = ObitTableSNRowUnref(in->SNTableRow);

  /* Smoothed table to be zapped? */
  if (in->isSNSmoo) {
    tabtype = g_strdup(in->SNTable->tabType);
    lver = in->SNTable->tabVer;
    ObitUVZapTable (in->myUV, tabtype, lver, err);
    if (tabtype) g_free (tabtype); tabtype = NULL;
    ObitUVUpdateTables (in->myUV, err);  /* Update disk header */
  }
  in->isSNSmoo = FALSE;

  if (in->PriorAntTime)  g_free (in->PriorAntTime);  in->PriorAntTime  = NULL;
  if (in->FollowAntTime) g_free (in->FollowAntTime); in->FollowAntTime = NULL;
  if (in->CalApply)      g_free (in->CalApply);      in->CalApply      = NULL;
  if (in->CalPrior)      g_free (in->CalPrior);      in->CalPrior      = NULL;
  if (in->CalFollow)     g_free (in->CalFollow);     in->CalFollow     = NULL;
  if (in->IFR)           g_free (in->IFR);           in->IFR           = NULL;
  if (in->PriorIFR)      g_free (in->PriorIFR);      in->PriorIFR      = NULL;
  if (in->FollowIFR)     g_free (in->FollowIFR);     in->FollowIFR     = NULL;
  if (in->MBDelay)       g_free (in->MBDelay);       in->MBDelay       = NULL;
  if (in->PriorMBDelay)  g_free (in->PriorMBDelay);  in->PriorMBDelay  = NULL;
  if (in->FollowMBDelay) g_free (in->FollowMBDelay); in->FollowMBDelay = NULL;
  if (in->RefAnt)        g_free (in->RefAnt);        in->RefAnt        = NULL;
  if (in->RateFact)      g_free (in->RateFact);      in->RateFact      = NULL;
  if (in->MissAnt)       g_free (in->MissAnt);       in->MissAnt       = NULL;

} /*  end ObitUVSolnShutDown */

/**
 * Update ObitUVSoln calibration tables for time time.
 * The current table is interpolated between the previous and following
 * sets of solutions.
 * If a new set of entries is needed from the SN/CL table they are read.
 * Adopted from AIPS CGASET.FOR
 * A calibration CalApply entry  consists of:
 * \li real part of gain
 * \li imaginary of gain
 * \li group delay (sec)
 * \li fringe rate (sec/sec)
 * \li weight
 * \param in   Calibrate Object.
 * \param time desired time in days
 * \param err  Error stack for messages and errors.
 */
void ObitUVSolnUpdate (ObitUVSoln *in, ofloat time, olong SourID, ObitErr *err)
{
  olong iant, iif, ipol, index, jndex, na, nturn;
  gboolean   good1, good2, bad;
  ofloat      wtt1, wtt2, phase1, phase2, amp, phase=0.0,
    g1a, g1p, g2a, g2p, delta, wt1, wt2, dt, dt1, dt2, rate1, rate2, temp, 
    phi0, phi1, phi2, phi3;
  ofloat fblank = ObitMagicF();
  gboolean newcal;
  gchar *routine="ObitUVSolnUpdate";
 
 
  /* see if time for new table entry */
  if ((in->LastRowRead <= in->numRow)  &&  (time > in->FollowCalTime)) {
    ObitUVSolnNewTime (in, time, err);
    if (err->error) Obit_traceback_msg (err, routine, "unspecified");
    newcal = TRUE;
  } else {
    newcal = FALSE;
  }

  /* see if calibration needs update; every 0.03 of solution interval. */
  delta = (time - in->CalTime);  
  if ((!newcal) &&  (delta <= 0.03*(in->FollowCalTime-in->PriorCalTime))) return;

  /* interpolate current calibration to time */
  in->CalTime = time;

  /* loop thru antennas */
  for (iant=0; iant<in->numAnt; iant++) { /* loop 500 */
    /* initialize indices for CalApply (index),  CalPrior and CalFollow (jndex)*/
    index = iant * in->numIF * in->numPol * in->lenCalArrayEntry;
    jndex = iant * in->numIF * in->numPol * in->lenCalArrayEntry;

    if (in->MissAnt[iant]) continue;  /* Nothing for this antenna */

    /* set interpolation weights proportional to time difference. */
    dt  = in->FollowAntTime[iant] - in->PriorAntTime[iant] + 1.0e-20;
    dt1 = time - in->PriorAntTime[iant];
    dt2 = time - in->FollowAntTime[iant];
    wtt1 = 0.0;
    if ((time < in->FollowAntTime[iant]) && (dt>1.0e-19)) wtt1 = -dt2 / dt;
    wtt2 = 1.0 - wtt1;

    /* Set Faraday rotation - interpolate if both good */
    if ((in->PriorIFR[iant]  !=  fblank) && (in->FollowIFR[iant] != fblank)) {
      in->IFR[iant] = wtt1 * in->PriorIFR[iant] + wtt2 * in->FollowIFR[iant];
    } else if (in->PriorIFR[iant] != fblank) {
      in->IFR[iant] = in->PriorIFR[iant];
    } else if (in->FollowIFR[iant] != fblank) {
      in->IFR[iant] = in->FollowIFR[iant];
    } else {
      in->IFR[iant] = fblank;
    } 

    /* Set multiband delay - interpolate if both good */
    /* Average rates */
    rate1 =  avgRatePrior (in, iant+1, 1);
    rate2 =  avgRateFollow (in, iant+1, 1);
    if ((in->PriorMBDelay[iant]  !=  fblank) && (in->FollowMBDelay[iant] != fblank)) {
      in->MBDelay[iant] = 
	wtt1 * (in->PriorMBDelay[iant]  + dt1 * rate1) + 
	wtt2 * (in->FollowMBDelay[iant] + dt2 * rate2);
    } else if (in->PriorMBDelay[iant] != fblank) {
      in->MBDelay[iant] = in->PriorMBDelay[iant] + dt1 * rate1;
    } else if (in->FollowMBDelay[iant] != fblank) {
      in->MBDelay[iant] = in->FollowMBDelay[iant] + dt2 * rate2;
    } else {
      in->MBDelay[iant] = fblank;
    } 
    /* Second poln MB delay */
    if (in->numPol>1) {
      /* Average rates */
      rate1 =  avgRatePrior (in, iant+1, 2);
      rate2 =  avgRateFollow (in, iant+1,2);
      na = in->numAnt;
      if ((in->PriorMBDelay[iant+na]  !=  fblank) && (in->FollowMBDelay[iant+na] != fblank)) {
	in->MBDelay[iant+na] = 
	  wtt1 * (in->PriorMBDelay[iant+na]  + dt1 * rate1) + 
	  wtt2 * (in->FollowMBDelay[iant+na] + dt2 * rate2);
      } else if (in->PriorMBDelay[iant+na] != fblank) {
	in->MBDelay[iant+na] = in->PriorMBDelay[iant+na] + dt1 * rate1;
      } else if (in->FollowMBDelay[iant+na] != fblank) {
	in->MBDelay[iant+na] = in->FollowMBDelay[iant+na] + dt2 * rate2;
      } else {
	in->MBDelay[iant+na] = fblank;
      } 
    }
    
    /* loop thru IF */
    for (iif=0; iif<in->numIF; iif++) { /* loop 400 */

      /* loop thru polarization */
      for (ipol= 0; ipol<in->numPol; ipol++) { /* loop 300 */

	/* initialize soln with blanks */
	in->CalApply[index]   = fblank;
	in->CalApply[index+1] = fblank;
	in->CalApply[index+2] = fblank;
	in->CalApply[index+3] = fblank;
	in->CalApply[index+4] = fblank;
	in->CalApply[index+5] = fblank;

	/* check for blanked soln. */
	good1 = (in->CalPrior[jndex] != fblank)  && 
	  (in->CalPrior[jndex+1] != fblank)  && 
	  (in->CalPrior[jndex+2] != fblank)  && 
	  (in->CalPrior[jndex+3] != fblank)  && 
	  (in->PriorAntTime[iant]>-100.0)    && 
	  (wtt1 > 0.0);
	good2 = (in->CalFollow[jndex] != fblank)  && 
	  (in->CalFollow[jndex+1] != fblank)  && 
	  (in->CalFollow[jndex+2] != fblank)  && 
	  (in->CalFollow[jndex+3] != fblank)  && 
	  (in->FollowAntTime[iant]>-100.0)    && 
	  (wtt2 > 0.0);

	/* If 'SELF' interpolation make sure it's the same source */
	if (in->interMode==OBIT_UVSolnInterSELF) {
	  if ((SourID!=in->PriorSourID)  && 
	      (SourID>0) && (in->PriorSourID>0)) good1 = FALSE;
	  if ((SourID!=in->FollowSourID)   && 
	      (SourID>0) && (in->FollowSourID>0)) good2 = FALSE;

	  /* If both OK pick closest */
	  if (good1 && good2) {
	    good1 = wtt1  >=  wtt2;
	    good2 = wtt1  <  wtt2;
	  }
	} /* end 'SELF' interpolation */

	/* Impose maximum interpolation time */
	if (fabs(dt1) > in->maxInter)  good1 = FALSE;
	if (fabs(dt2) > in->maxInter)  good2 = FALSE;

	/* solution all flagged? */
	bad = !(good1 || good2);

	/* nothing more if both prior and following are bad */
	if (!bad) {

	  /* different reference antennas  use closest */
	  if ((fabs (in->CalPrior[jndex+5] - in->CalFollow[jndex+5]) >= 0.5)
	      &&  good1  &&  good2) {
	    good1 = wtt1  >=  wtt2;
	    good2 = wtt1  <  wtt2;
	  } 
	  
	  /* initial weights */
	  wt1 = wtt1;
	  wt2 = wtt2;
	  
	  /* Only Following good */
	  if (!good1) {
	    wt1 = 0.0;
	    wt2 = 1.0;
	  } 
	  
	  /* Only Prior good */
	  if (!good2) {
	    wt1 = 1.0;
	    wt2 = 0.0;
	  } 
	  
	  /* Get gains from tables */
	  if (good1) {
	    g1a = in->CalPrior[jndex];
	    g1p = in->CalPrior[jndex+1];
	  } else {  /* if not good use following */
	    g1a = in->CalFollow[jndex];
	    g1p = in->CalFollow[jndex+1];
	  }
	  if (good2) {
	    g2a = in->CalFollow[jndex];
	    g2p = in->CalFollow[jndex+1];
	  } else {  /* if not good use preceeding */
	    g2a = in->CalPrior[jndex];
	    g2p = in->CalPrior[jndex+1];
	  }
	  
	  /* check if fringe rates given - if so add accumulated phase */
	  if ((in->CalPrior[jndex+3]  != 0.0)  || (in->CalFollow[jndex+3] != 0.0)) {
	    phase1 = 0.0;
	    phase2 = 0.0;
	    if (good1) phase1 = in->CalPrior[jndex+3]  * dt1 * in->RateFact[iif];
	    if (good2) phase2 = in->CalFollow[jndex+3] * dt2 * in->RateFact[iif];

	    /* rotate phase by accumulated fringe rate */
	    g1p += phase1;
	    g2p += phase2;

	    /* Put in same turn */
	    g1p = fmod(g1p, 2.0*G_PI);
	    if (g1p<0.0) g1p += 2.0*G_PI;
	    g2p = fmod(g2p, 2.0*G_PI);
	    if (g2p<0.0) g2p += 2.0*G_PI;
	  } 
	  
	  /* interpolate phase by mode */
	  switch (in->interMode) {
	  case OBIT_UVSolnInterSELF:
	  case OBIT_UVSolnInter2PT:
	  case OBIT_UVSolnInterPOLY:
	  case OBIT_UVSolnInterSIMP:
	    /* Both good? */
	    if (good1 && good2 ) {
	      /* resolve phase ambiguity by disallowing phase jumps of >180 deg */
	      if (fabs(g2p-g1p)> G_PI) {
		if ((g2p-g1p)> 0.0) g2p -=  2.0 * G_PI;
		else g2p +=  2.0 * G_PI;
	      }
	      phase = g1p + (g2p - g1p) * dt1 / dt;
	    } else if (good1) { /* Only prior */
	      phase = g1p;
	    } else if (good2) { /* Only following */
	      phase = g2p;
	    }
	    break;
	  case OBIT_UVSolnInterAMBG: /* Use rates to resolve ambiguities */
	    /* Both good? */
	    if (good1 && good2 ) {
	      /* rates in radians per day
	      rate1 = in->CalPrior[jndex+3]  * in->RateFact[iif];
	      rate2 = in->CalFollow[jndex+3] * in->RateFact[iif];
	      temp = (g1p + 0.5*(rate1*dt1+rate2*dt2) - g2p) / (2.0*G_PI); */
	      /* How many turns of phase 
	      if (temp>=0.0) nturn = (olong)(temp+0.5);
	      else nturn = (olong)(temp-0.5);
	      g2p += (ofloat)(nturn) * 2.0 * G_PI;*/
	      phase = g1p + (g2p - g1p) * dt1 / dt;
 	    } else if (good1) { /* Only prior */
	      phase = g1p;
	    } else if (good2) { /* Only following */
	      phase = g2p;
	    }
	    break;
	  case OBIT_UVSolnInterCUBE:  /* CUBIC method */
	    /* Both good? */
	    if (good1 && good2 ) {
	      /* rates in radians per second */
	      rate1 = in->CalPrior[jndex+3]  * in->RateFact[iif];
	      rate2 = in->CalFollow[jndex+3] * in->RateFact[iif];
	      temp = (g1p + 0.5*(rate1+rate2) - g2p) / (2.0*G_PI);
	      /* How many turns pf phase */
	      if (temp>=0.0) nturn = (olong)(temp+0.5);
	      else nturn = (olong)(temp-0.5);
	      g2p += (ofloat)(nturn) * 2.0 * G_PI;
	      phi0 = g1p;
	      phi1 = rate1;
	      phi2 = (3.0*(g2p-g1p) - 2.0*rate1*dt-rate2*dt)/(dt*dt);
	      phi3 = -(2.0*(g2p-g1p) - rate1*dt-rate2*dt)/(dt*dt*dt);
	      phase = phi0 + phi1*dt1 + phi2*dt1*dt1 + phi3*dt1*dt1*dt1;
	    } else if (good1) { /* Only prior */
	      phase = g1p;
	    } else if (good2) { /* Only following */
	      phase = g2p;
	    }
	    break;
	  default:
	    break;
	  }; /* end switch by interpolation mode */

	  /* interpolate amplitude */
	  amp   = wt1 * g1a + wt2 * g2a;

	  /* set amplitude and phase in output array */
	  in->CalApply[index]   = amp * cos(phase);
	  in->CalApply[index+1] = amp * sin(phase);

	  /* interpolate delay */
	  in->CalApply[index+2] = (wt1 * in->CalPrior[jndex+2] +
				   wt2 * in->CalFollow[jndex+2]);
	  /* interpolate rate */
	  in->CalApply[index+3] = (wt1 * in->CalPrior[jndex+3] +
				   wt2 * in->CalFollow[jndex+3]);
	  /* interpolate weight */
	  in->CalApply[index+4] = (wt1 * in->CalPrior[jndex+4] +
				   wt2 * in->CalFollow[jndex+4]);

	  /* reference antenna */
	  if (good1) in->RefAnt[iant+ipol*in->numAnt] = 
		       (olong)(in->CalPrior[jndex+5] + 0.5);
	  else if (good2) in->RefAnt[iant+ipol*in->numAnt] = 
			    (olong)(in->CalFollow[jndex+5] + 0.5);
	} /* end of only valid solutions section */

	/* update indices */
        index += in->lenCalArrayEntry;
	jndex += in->lenCalArrayEntry;

      } /* end poln loop  L300: */;
    } /* end IF loop  L400: */;
  } /* end antenna loop  L500: */

} /* end ObitUVSolnUpdate */
    
/**
 * References the phases to a common reference antenna in a  
 * polarization coherent fashion.  
 * Leaves the output table sorted in antenna-time order.  
 * Routine translated from the AIPSish UVUTIL.FOR/SLFREF  
 * \param SNTab  SN table object 
 * \param isuba  Desired subarray, 0=> 1 
 * \param refant Reference antenna, if 0 then the most commonly used 
 *               reference antenna is picked. 
 *               If <0 then only sort and not reference table
 * \param err    Error/message stack, returns if error.
 */
void ObitUVSolnRefAnt (ObitTableSN *SNTab, olong isuba, olong* refant, ObitErr* err) 
{
  ObitIOCode retCode;
  olong   ant, iif, isub, numif, numpol, numant, numtime, iref, crefa, *antuse=NULL;
  ofloat *wrkTime=NULL, *work1=NULL, *work2=NULL, *work3=NULL, *work4=NULL, *work5=NULL;
  gboolean noReRef = FALSE;
  gchar msgtxt[81];
  gchar *routine = "ObitUVSolnRefAnt";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
 
  /* Subarray */
  isub = MAX (1, isuba);

  /* Must be antenna-time order */
  ObitTableUtilSort2f ((ObitTable*)SNTab , "TIME  ", 1, FALSE, "ANTENNA",
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Do we really want to rereference? */
  noReRef = ((*refant)<0);
  if (noReRef) return;
  
  /* Open table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
 
  /* Get descriptive info */
  numif  = SNTab->numIF;
  numpol = SNTab->numPol;
  numant = SNTab->numAnt;

  /* Determine which antennas used as reference antennas, number of times */
  antuse = refCount (SNTab, isub, &numtime, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Create work arrays */
  numtime += 10;  /* Fudge a bit */
  wrkTime = g_malloc0(numtime*sizeof(ofloat));
  work1   = g_malloc0(numtime*sizeof(ofloat));
  work2   = g_malloc0(numtime*sizeof(ofloat));
  work3   = g_malloc0(numtime*sizeof(ofloat));
  work4   = g_malloc0(numtime*sizeof(ofloat));
  work5   = g_malloc0(numtime*sizeof(ofloat));

  /* Determine reference antenna if  necessary. */
  if (*refant <= 0) {
    (*refant) = 1;
    crefa = antuse[*refant-1];
    /* Pick most commonly used */
    for (ant=1; ant<numant; ant++) { /* loop 30 */
      if (antuse[ant] > crefa) {
	*refant = ant+1;
	crefa = antuse[(*refant)-1];
      } 
    } /* end loop  L30:  */;
  } /* end picking reference antenna */
  iref = *refant;

  /* Message about rereferencing. */
  g_snprintf (msgtxt,80, "Rereferencing phases to antenna %3d", *refant);
  Obit_log_error(err, OBIT_InfoErr, msgtxt);

  /* Loop through antennas used as  secondary reference antennas. */
  for (ant= 1; ant<=numant; ant++) { /* loop 500 */
    if ((antuse[ant-1] <= 0)  ||  (ant == *refant)) continue; /*goto L500;*/
    /* Multiband delays */
    refMBDelay (SNTab, isub, iref, ant, numtime, wrkTime, 
		work1, work2, work4, err);
    
    for (iif=0; iif<numif; iif++) { /* loop 300 */

      /* Delay */
      refDelay (SNTab, isub, iif, iref, ant, numtime, wrkTime, 
		work1, work2, work4, err);
      /* Rate */
      refRate (SNTab, isub, iif, iref, ant, numtime, wrkTime, 
		work1, work2, work4, err);
      /* Phases - this one will change the reference antenna - must be last */
      refPhase (SNTab, isub, iif, iref, ant, numtime, wrkTime, 
		work1, work2, work3, work4, work5, err);
      if (err->error) {
	Obit_log_error(err, OBIT_InfoErr,
		       "%s: ERROR ref. ant %3d phases to antenna %3d IF %d in %s", 
		       routine, ant, *refant, iif+1, SNTab->name);
	goto cleanup;
      }
    } /* end IF loop  L300: */;
  } /* end of antenna loop  L500: */;
  
  /* Close table */
  retCode = ObitTableSNClose (SNTab, err);
  if (err->error) goto cleanup;
  
  /* Cleanup */
  cleanup: if (antuse)  g_free (antuse);
  if (wrkTime) g_free (wrkTime);
  if (work1)   g_free (work1);
  if (work2)   g_free (work2);
  if (work3)   g_free (work3);
  if (work4)   g_free (work4);
  if (work5)   g_free (work5);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

} /* end of routine ObitUVSolnRefAnt */ 

/**
 * Averages the rates over poln, IF in an SN table
 * \param SNTab  SN table object 
 * \param err    Error/message stack, returns if error.
 */
void ObitUVSolnAvgRate (ObitTableSN *SNTab, ObitErr* err) 
{
  gchar *routine = "ObitUVSolnAvgRate";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));

  /* Open table */
  ObitTableSNOpen (SNTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Average rates */
  avgRate (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Close output table */
  ObitTableSNClose (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
} /* end of routine ObitUVSolnAvgRate */ 

/**
 * Smooths the SN table which should already be referenced to a single  
 * reference antenna.  Failed solutions are optionally interpolated.  
 * Leaves the output table sorted in antenna-time order.  
 * If SmoType='VLMB the host ObitUV of SNTab should have a 
 * valid selector and descriptor.
 * All rates are averaged.
 * Routine adapted from the AIPSish UVUTIL.FOR/SLFSMO  
 * Table MUST be in Time order when called and will be returned in antenna order
 * Controls on SNTab:
 * \li smoFunc   OBIT_string (4,1,1) smoothing function 
 *               'MWF' (median window filter), 
 *               "GAUS' (Gaussian) else "BOX", [def "BOX"]
 * \li smoParm   OBIT_float (5,1,1) Amplitude smoothing time in hr. [def 0.0]
 *               ampl, phase, rate, singleband delay, multiband delay
 * \li smoType   OBIT_string (4,1,1) Data to be smoothed
 *               'AMPL', 'PHAS', 'BOTH'[def], 'DELA', 'DERA', 
 *               'VLBI','VLMB', 'FULL', '    '='AMPL'
 * \li doBlank   OBIT_bool (1,1,1) Replace blanked values with interpolated? [def false]
 * \param SNTab  SN table object 
 * \param isuba  Desired subarray, 0=> 1 
 * \param err    Error/message stack, returns if error.
 */
void ObitUVSolnSNSmo (ObitTableSN *SNTab, olong isuba, ObitErr* err) 
{
  ObitIOCode retCode;
  gint32 dim[UV_MAXDIM];
  ObitInfoType type;
  ObitUV *host;
  gchar  smtype[5], smfunc[5], doType[5];
  olong   i, iif, isub, numant, numpol, numif;
  gboolean doBlank;
  ofloat smparm[5], gncnt, gnsum;
  olong   mxtime, ntmp, *antuse;
  ofloat *work1=NULL, *work2=NULL;
  gchar *routine = "ObitUVSolnSNSmo";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));

  /* Control Info */
  strcpy (smfunc, "BOX ");dim[0] = strlen(smfunc);
  ObitInfoListGetTest(SNTab->info, "smoFunc", &type, dim, smfunc);
  smfunc[dim[0]] = 0;
  strcpy (smtype, "BOTH");dim[0] = strlen(smtype);
  ObitInfoListGetTest(SNTab->info, "smoType", &type, dim, smtype);
  smtype[dim[0]] = 0;
  if (!strncmp (smtype, "    ",4)) strcpy (smtype, "AMPL");

  /* Smoothing times */
  for (i=0; i<5; i++) smparm[i] = 0.0;
  ObitInfoListGetTest(SNTab->info, "smoParm", &type, dim, smparm);
  for (i=0; i<5; i++) smparm[i] /= 24.0;  /* To days */
  doBlank = FALSE;
  ObitInfoListGetTest(SNTab->info, "doBlank", &type, dim, &doBlank);

  /* Subarray */
  isub = MAX (1, isuba);
  
  /* Determine which antennas used as reference antennas, number of times */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
  antuse = refCount (SNTab, isub, &mxtime, err);
  mxtime += 10;  /* Fudge a bit on the number of times */
  g_free(antuse);  /* don't need */
  retCode = ObitTableSNClose (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* MUST be antenna-time order */
  ObitTableUtilSort2f ((ObitTable*)SNTab, "ANTENNA", 1, FALSE, "TIME  ", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Open table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
 
  /* Get descriptive info */
  numif  = SNTab->numIF;
  numpol = SNTab->numPol;
  numant = SNTab->numAnt;

  /* Average rates */
  avgRate (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);

  /* Reset mean gain modulus info  */
  gncnt = 0.0;
  gnsum = 0.0;

  /* VLMB handled separately - delay smoothed first */
  if (!strncmp(smtype, "VLMB", 4)) strncpy (doType, "DERA", 4);
  else strncpy (doType, smtype, 4);
  doType[4] = 0;

  /* Primary work array */
  work1 = g_malloc(16*mxtime*sizeof(ofloat));
  
  /* How many secondary work arrays */
  if (!strncmp(smfunc, "GAUS", 4)) ntmp = 3;
  else if (!strncmp(smfunc, "MWF", 3)) ntmp = 4;
  else ntmp = 2;
  work2 = g_malloc(ntmp*mxtime*sizeof(ofloat));
  
  for (iif=0; iif<numif; iif++) { /* loop over IFs 300 */
    /* Smooth this IF */
    ObitUVSolnSNSmooth (SNTab, smfunc, doType, 0.5, smparm, iif, isub, 
			&gncnt, &gnsum, mxtime, work1, work2, doBlank, err);
    if (err->error) goto cleanup;
  } /* end loop  L300: */
  if (work1) g_free(work1); work1 = NULL;
  if (work2) g_free(work2); work2 = NULL;

  /* If VLMB smooth phase coherently */
      if (!strncmp(smtype, "VLMB", 4)) {
	host = (ObitUV*)SNTab->myHost;
	smoAmpPh (SNTab, host->mySel, smfunc, "PHAS", 0.5, smparm, isub, 1, numif, 
		  host->myDesc->freqIF, doBlank, &gncnt, &gnsum, err);
	if (err->error) goto cleanup;
      }

  /* Save mean gain modulus if requested */
  if ((fabs (SNTab->mGMod-1.0) > 1.0e-5)  &&  (gncnt > 0.1)) {
    SNTab->mGMod  = gnsum / gncnt;
  } 
  
  /* Close output table */
cleanup: 
  retCode = ObitTableSNClose (SNTab, err);
  
  /* Cleanup */
  if (work1) g_free(work1);
  if (work2) g_free(work2);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
} /* end of routine ObitUVSolnSNSmo */ 

/**
 * Routine to deselect records in an SN table if they match a given
 * subarray, have a selected FQ id, appear on a list of antennas and
 * are in a given timerange.
 * \param SNTab      SN table object 
 * \param isub       Subarray number, <=0 -> any
 * \param fqid       Selected FQ id, <=0 -> any
 * \param nantf      Number of antennas in ants
 * \param ants       List of antennas, NULL or 0 in first -> flag all
 * \param nsou       Number of source ids in sources
 * \param sources    List of sources, NULL or 0 in first -> flag all
 * \param timerange  Timerange to flag, 0s -> all
 * \param err        Error/message stack, returns if error.
 */
void ObitUVSolnDeselSN (ObitTableSN *SNTab, olong isuba, olong fqid, 
			olong nantf, olong *ants, olong nsou, olong *sources,
			ofloat timerange[2], ObitErr* err)
{
  ObitIOCode retCode;
  ObitTableSNRow *row=NULL;
  olong   i;
  gboolean allAnt, allSou, flagIt;
  ofloat tr[2];
  olong  loop;
  gchar  *routine = "ObitUVSolnDeselSN";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
 
 /* Open table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
 
  /* All antennas to be flagged? */
  allAnt = ((ants==NULL) || (ants[0]<=0));
  /* All sources to be flagged? */
  allSou = ((sources==NULL) || (sources[0]<=0));

  /* Timerange */
  tr[0] = timerange[0];
  tr[1] = timerange[1];
  if ((tr[0]==0.0) && (tr[1]==0.0)) {
    tr[0] = -1.0e20;
    tr[1] =  1.0e20;
  }

  /* Create Row */
  row = newObitTableSNRow (SNTab);
  /* Loop through table */
  for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop 20 */

    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    if (row->status<0) continue;  /* Skip deselected record */

    /* Flag this one? */
    flagIt = (row->SubA==isuba) || (isuba<=0);                     /* by subarray */
    flagIt = flagIt || (row->FreqID==fqid) || (fqid<=0);           /* by FQ id */
    flagIt = flagIt || ((row->Time>=tr[0]) && (row->Time<=tr[1])); /* by time */
    flagIt = flagIt && allAnt;                          /* Check Antenna */
    if (!flagIt) for (i=0; i<nantf; i++) {
      if (row->antNo==ants[i]) {flagIt = TRUE; break;}
    }
    flagIt = flagIt && allSou;                          /* Check Antenna */
    if (!flagIt) for (i=0; i<nsou; i++) {
      if (row->SourID==sources[i]) {flagIt = TRUE; break;}
    }

    if (flagIt) { /* flag */
      /* Rewrite record flagged */
      row->status = -1;
      retCode = ObitTableSNWriteRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    } /* end if flag */
  } /* End loop over table */

  /* Close table */
  retCode = ObitTableSNClose (SNTab, err);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
  row = ObitTableSNRowUnref (row);  /* Cleanup */
} /* end ObitUVSolnDeselSN */

/**
 * Routine to deselect records in an CL table if they match a given
 * subarray, have a selected FQ id, appear on a list of antennas and
 * are in a given timerange.
 * \param CLTab      CL table object 
 * \param isub       Subarray number, <=0 -> any
 * \param fqid       Selected FQ id, <=0 -> any
 * \param nantf      Number of antennas in ants
 * \param ants       List of antennas, NULL or 0 in first -> flag all
 * \param nsou       Number of source ids in sources
 * \param sources    List of sources, NULL or 0 in first -> flag all
 * \param timerange  Timerange to flag, 0s -> all
 * \param err        Error/message stack, returns if error.
 */
void ObitUVSolnDeselCL (ObitTableCL *CLTab, olong isuba, olong fqid, 
			   olong nantf, olong *ants, olong nsou, olong *sources,
			   ofloat timerange[2], ObitErr* err)
{
  ObitIOCode retCode;
  ObitTableCLRow *row=NULL;
  olong   i;
  gboolean allAnt, allSou, flagIt;
  ofloat tr[2];
  olong  loop;
  gchar  *routine = "ObitUVSolnDeselCL";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitTableCLIsA(CLTab));
 
 /* Open table */
  retCode = ObitTableCLOpen (CLTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, CLTab->name);
 
  /* All antennas to be flagged? */
  allAnt = ((ants==NULL) || (ants[0]<=0));
  /* All sources to be flagged? */
  allSou = ((sources==NULL) || (sources[0]<=0));

  /* Timerange */
  tr[0] = timerange[0];
  tr[1] = timerange[1];
  if ((tr[0]==0.0) && (tr[1]==0.0)) {
    tr[0] = -1.0e20;
    tr[1] =  1.0e20;
  }

  /* Create Row */
  row = newObitTableCLRow (CLTab);
  /* Loop through table */
  for (loop=1; loop<=CLTab->myDesc->nrow; loop++) { /* loop 20 */

    retCode = ObitTableCLReadRow (CLTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, CLTab->name);
    if (row->status<0) continue;  /* Skip deselected record */

    /* Flag this one? */
    flagIt = (row->SubA==isuba) || (isuba<=0);             /* by subarray */
    flagIt = flagIt || (row->FreqID==fqid) || (fqid<=0);      /* by FQ id */
    flagIt = flagIt || ((row->Time>=tr[0]) && (row->Time<=tr[1])); /* by time */
    flagIt = flagIt && allAnt;                          /* Check Antenna */
    if (!flagIt) for (i=0; i<nantf; i++) {
      if (row->antNo==ants[i]) {flagIt = TRUE; break;}
    }
    flagIt = flagIt && allSou;                          /* Check Antenna */
    if (!flagIt) for (i=0; i<nsou; i++) {
      if (row->SourID==sources[i]) {flagIt = TRUE; break;}
    }

    if (flagIt) { /* flag */
      /* Rewrite record flagged */
      row->status = -1;
      retCode = ObitTableCLWriteRow (CLTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, CLTab->name);
    } /* end if flag */
  } /* End loop over table */

  /* Close table */
  retCode = ObitTableCLClose (CLTab, err);
  if (err->error) Obit_traceback_msg (err, routine, CLTab->name);
  row = ObitTableCLRowUnref (row);  /* Cleanup */
} /* end ObitUVSolnDeselCL */

/**
 * Routine to smooth amplitudes and/or phases rates in an open SN  
 * table.  All poln present are smoothed but only one IF.  The KOL  
 * pointers are presumed to point at the desired IF.  An error is  
 * returned if there are any non-zero delays, rates, or multi-band  
 * delays.  If the reference antenna changes and phase is being  
 * smoothed, an error is returned.  
 * Multiband delays only smoothed if iif==0;
 * Input table must be in antenna-time order.  
 * Routine translated from the AIPSish SNSMOO.FOR/SNSMOO 
 * \param SNTab  SN table object; must be opened/closed externally
 * \param smoFunc  Smoothing function: 'MWF', 'GAUS', else BOX 
 * \param smoType  Type of data to smooth
 *                 'AMPL', 'PHAS', 'BOTH'[def], 'VLBI','FULL' 
 *                 'DERA' = Delay & Rate only
 * \param alpha    Alpha clip for MWF (0 -> box, 1 -> pure MWF) 
 * \param smoParm  Smoothing time in days for:
 *                 ampl, phase, rate, singleband delay, multiband delay
 *                 0=>fill in for blanked only. 
 * \param iif      Desired IF (0-rel)
 * \param sub      Desired subarray (1-rel)
 * \param gncnt    [In/Out] Count for gain normalization 
 * \param gnsum    [In/Out] Sum of gain modulii 
 * \param nxt      Number of times allowed in WRK 
 * \param work1    Work buffer (NXT*16) 
 * \param work2    Work area >= (NXT*m)  (m=2 BOX, 3 GAUS, 4 MWF) 
 * \param doBlank  replace blanked values with interpolated values?
 * \param err      Error/message stack, returns if error.
 */
void 
ObitUVSolnSNSmooth (ObitTableSN *SNTab, gchar* smoFunc, gchar* smoType, ofloat alpha, 
		    ofloat *smoParm, olong iif, olong sub, ofloat* gncnt, ofloat* gnsum, 
		    olong nxt, ofloat* work1, ofloat* work2, gboolean doBlank, ObitErr* err) 
{
  olong   loopa, numtim, ant, numrec, nleft, isnrno=0, itime, 
    refa1, refa2, n1good, n2good, i, numif, numpol, numant;
  olong loopr, fstrec, save;
  ofloat    amp, stamp, stph, stMB, stdela, strate, fblank =  ObitMagicF();
  gboolean  need2, doamp, doph, doMB, dodela, dorate, isVLBI, isVLBIMB;
  odouble   timoff=0.0;
  ObitIOCode retCode;
  ObitTableSNRow *row=NULL;
  gchar *routine = "ObitUVSolnSNSmooth";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
  
  /* Get number of records in table */
  numrec = SNTab->myDesc->nrow;
  if (numrec <= 0) return;   /* bail if nothing */
  
  /* Get descriptive info */
  numif  = SNTab->numIF;
  numpol = SNTab->numPol;
  numant = SNTab->numAnt;
  
  /* averaging times */
  stamp   = smoParm[0];
  stph    = smoParm[1];
  stdela  = smoParm[2];
  strate  = smoParm[3];
  stMB    = smoParm[4];

  /* Are there 2 polarizations? */
  need2 = numpol>1;
  refa1 = 0;
  refa2 = 0;

  /* Any delay/rate smoothing? */
  isVLBI = (!strncmp(smoType, "VLBI",4)) || (!strncmp(smoType, "VLMB",4)) ||
    (!strncmp(smoType, "FULL",4));
  isVLBIMB =  (!strncmp(smoType, "VLMB",4));

  doamp = FALSE;
  doph  = FALSE;
  doMB  = FALSE;
  dodela  = FALSE;
  dorate  = FALSE;
  /* Only amp/phase? */
  if (!strncmp(smoType, "BOTH",4)) {
    doamp   = TRUE;
    doph    = TRUE;
    stdela  = 0.0;
    strate  = 0.0;
    stMB    = 0.0;

    /* Only amplitude? */
  } else if (!strncmp(smoType, "AMPL",4)) {
    doamp   = TRUE;
    stph    = 0.0;
    stdela  = 0.0;
    strate  = 0.0;
    stMB    = 0.0;

  /* Only phase? */
  } else if (!strncmp(smoType, "PHAS",4)) {
    doph    = TRUE;
    stamp   = 0.0;
    stdela  = 0.0;
    strate  = 0.0;
    stMB    = 0.0;

  /* Only delay? */
  } else if (!strncmp(smoType, "DELA",4)) {
    dodela  = TRUE;
    stamp   = 0.0;
    stph    = 0.0;
    strate  = 0.0;

  /* Only delay & Rate? */
  } else if (!strncmp(smoType, "DERA",4)) {
    stamp   = 0.0;
    stph    = 0.0;
    dodela  = TRUE;
    dorate  = TRUE;
 
    /* VLBI or FULL */
 } else if ((!strncmp(smoType, "VLBI",4)) || (!strncmp(smoType, "FULL",4))) {
    stamp   = 0.0;
    doph    = TRUE;
    dodela  = TRUE;
    dorate  = TRUE;
 }

  /* Create Row */
  row = newObitTableSNRow (SNTab);
  /* Attach row to output buffer */
  ObitTableSNSetRow (SNTab, row, err);
  fstrec = 0;  /* Record number read in table */
  
  /* Loop over antenna */
  for (loopa= 1; loopa<=numant; loopa++) { /* loop 600 */
    ant = loopa;
    /* Set pointers, counters */
    numtim = 0;
    nleft = SNTab->myDesc->nrow - fstrec;  /* How many rows? */
    n1good = 0;
    n2good = 0;
    /* Loop in time, reading */
    for (loopr=1; loopr<=nleft; loopr++) { /* loop 100 */
      isnrno = fstrec + loopr;
      retCode = ObitTableSNReadRow (SNTab, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Finished antenna? */
      if (row->antNo < ant) continue; /* Shouldn't happen */
      if (row->antNo > ant) break;

      /* Want this record */
      if ((row->SubA == sub)  &&  (row->antNo == ant)) {

	/* take first reference antenna and check that it never changes */
	if (refa1 <= 0) {
	  refa1 = row->RefAnt1[iif];
	} else {
	  if (row->RefAnt1[iif] == 0) row->RefAnt1[iif] = refa1;
	} 
	if (refa1 != row->RefAnt1[iif]) { /* goto L980;*/
	  Obit_log_error(err, OBIT_Error, 
			 "%s: Reference antenna varies, cannot smooth %s", 
			 routine, SNTab->name);
	  row = ObitTableSNRowUnref(row); /* delete row object */
	  return;
	}
	
	if (need2) {  /* Get second polarization */
	  /* take first reference antenna and check that it never changes */
	  if (refa2 <= 0) {
	    refa2 = row->RefAnt2[iif];
	  } else {
	    if (row->RefAnt2[iif] == 0) row->RefAnt2[iif] = refa2;
	  } 
	  if (refa2 != row->RefAnt2[iif]) { /* goto L980;*/
	    Obit_log_error(err, OBIT_Error, 
			   "%s: Reference antenna varies, cannot smooth %s", 
			   routine, SNTab->name);
	    row = ObitTableSNRowUnref(row); /* delete row object */
	    return;
	  }
	}  /* end get second polarization */ 
	
	/* Put in buffer */
	if (numtim >= nxt) {
	  Obit_log_error(err, OBIT_Error, 
			 "%s: Exceed time limit of %d for %s", 
			 routine, nxt, SNTab->name);
	  row = ObitTableSNRowUnref(row); /* delete row object */
	  return;
	} 
	/* Work1 usage :
	   0 = Real pol 1
	   1 = Imag pol 1
	   2 = Amplitude pol 1
	   3 = Weight pol 1
	   4 = Real pol 2
	   5 = Imag pol 2
	   6 = Amplitude pol 2
	   7 = Weight pol 2
	   8 = Time(day) relative to first
	   9 = row number
	   10 = MB delay1
	   11 = SB delay1
	   12 = Rate1
	   13 = MB delay2
	   14 = SB delay2
	   15 = Rate2
	*/
	
	if (numtim == 0) timoff = row->Time;  /* First time */
	work1[8*nxt+numtim] = row->Time - timoff;
	work1[9*nxt+numtim] = (ofloat)isnrno;
	if ((row->Weight1[iif] > 0.0)  &&  (row->Real1[iif] != fblank)) {
	  work1[2*nxt+numtim] = sqrt (row->Real1[iif]*row->Real1[iif] + 
				      row->Imag1[iif]*row->Imag1[iif]);
	} else {
	  work1[2*nxt+numtim] = fblank;
	} 
	/* First polarization */
	if ((work1[2*nxt+numtim] > 0.0) && (work1[2*nxt+numtim]!=fblank)) {
	  work1[0*nxt+numtim]  = row->Real1[iif] / work1[2*nxt+numtim];
	  work1[1*nxt+numtim]  = row->Imag1[iif] / work1[2*nxt+numtim];
	  work1[3*nxt+numtim]  = row->Weight1[iif];
	  work1[10*nxt+numtim] = row->MBDelay1;
	  work1[11*nxt+numtim] = row->Delay1[iif];
	  work1[12*nxt+numtim] = row->Rate1[iif];
	  n1good = n1good + 1;
	} else {
	  work1[0*nxt+numtim]  = fblank;
	  work1[1*nxt+numtim]  = fblank;
	  work1[3*nxt+numtim]  = fblank;
	  work1[10*nxt+numtim] = fblank;
	  work1[11*nxt+numtim] = fblank;
	  work1[12*nxt+numtim] = fblank;
	}
	
	if (need2) {  /* Second polarization */
	  if ((row->Weight2[iif] > 0.)  &&  (row->Real2[iif] != fblank)) {
	    work1[6*nxt+numtim] = sqrt (row->Real2[iif]*row->Real2[iif] + 
					row->Imag2[iif]*row->Imag2[iif]);
	  } else {
	    work1[6*nxt+numtim] = fblank;
	  } 
	  if ((work1[6*nxt+numtim] > 0.0) && (work1[6*nxt+numtim]!=fblank)) {
	    work1[4*nxt+numtim]  = row->Real2[iif] / work1[6*nxt+numtim];
	    work1[5*nxt+numtim]  = row->Imag2[iif] / work1[6*nxt+numtim];
	    work1[7*nxt+numtim]  = row->Weight2[iif];
	    work1[13*nxt+numtim] = row->MBDelay2;
	    work1[14*nxt+numtim] = row->Delay2[iif];
	    work1[15*nxt+numtim] = row->Rate2[iif];
	    n2good = n2good + 1;
	  } else {
	    work1[4*nxt+numtim]  = fblank;
	    work1[5*nxt+numtim]  = fblank;
	    work1[7*nxt+numtim]  = fblank;
	    work1[13*nxt+numtim] = fblank;
	    work1[14*nxt+numtim] = fblank;
	    work1[15*nxt+numtim] = fblank;
	  } 
	} /* end second polarization */
      } /* end if want record */ 
      numtim++;   /* count times */
    } /* end loop  L100: */

    save = isnrno - 1; /* How far did we get? */
    if (numtim <= 0) goto endAnt;  /* Catch anything? */
    
    /* Smooth as requested */
    if (n1good > 0) { /* First polarization */
      if (doamp) {
	smoIt (smoFunc, stamp, alpha, &work1[8*nxt], &work1[2*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[2*nxt+i] = work2[i];
      } else {  /* Not smoothing amp, replace fblank with 1.0 if doBlank */
	if (doBlank) {
	  for (i=0; i<numtim; i++) {
	    if (work1[2*nxt+i]==fblank) work1[2*nxt+i] = 1.0;
	  }
	}
      }
      if (doph) {
	smoIt (smoFunc, stph, alpha, &work1[8*nxt], &work1[0*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[i] = work2[i];
	smoIt (smoFunc, stph, alpha, &work1[8*nxt], &work1[1*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[1*nxt+i] = work2[i];
      }
      if (doMB) {
	smoIt (smoFunc, stMB, alpha, &work1[8*nxt], &work1[10*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[10*nxt+i] = work2[i];
      }
      if (dodela) {
	smoIt (smoFunc, stdela, alpha, &work1[8*nxt], &work1[11*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[11*nxt+i] = work2[i];
      }
      if (dorate) {
	smoIt (smoFunc, strate, alpha, &work1[8*nxt], &work1[12*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[12*nxt+i] = work2[i];
      }

      /* Save deblanked weights if doBlank */
      if (doBlank) for (i=0; i<numtim; i++) work1[3*nxt+i] = work2[1*nxt+i];
    } /* end first polarization */
    
    
    if (n2good > 0) {  /* Second polarization */
      if (doamp) {
	smoIt (smoFunc, stamp, alpha, &work1[8*nxt], &work1[6*nxt], &work1[7*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[6*nxt+i] = work2[i];
      } else {  /* Not smoothing amp, replace fblank with 1.0 if doBlank */
	if (doBlank) {
	  for (i=0; i<numtim; i++) {
	    if (work1[6*nxt+i]==fblank) work1[6*nxt+i] = 1.0;
	  }
	}
      }
      if (doph) {
	smoIt (smoFunc, stph, alpha, &work1[8*nxt], &work1[4*nxt], &work1[7*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[4*nxt+i] = work2[i];
	smoIt (smoFunc, stph, alpha, &work1[8*nxt], &work1[5*nxt], &work1[7*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[5*nxt+i] = work2[i];
      }
      if (doMB) {
	smoIt (smoFunc, stMB, alpha, &work1[8*nxt], &work1[13*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[13*nxt+i] = work2[i];
      }
      if (dodela) {
	smoIt (smoFunc, stdela, alpha, &work1[8*nxt], &work1[14*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[14*nxt+i] = work2[i];
      }
      if (dorate) {
	smoIt (smoFunc, strate, alpha, &work1[8*nxt], &work1[15*nxt], &work1[3*nxt], numtim, 
	       &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	for (i=0; i<numtim; i++) work1[15*nxt+i] = work2[i];
      }
      /* Save deblanked weights if doBlank */
      if (doBlank) for (i=0; i<numtim; i++) work1[7*nxt+i] = work2[1*nxt+i];
    } /* end second polarization */
    
    /* Replace with smoothed values */
    for (itime=0; itime<numtim; itime++) { /* loop 200 */
      isnrno = (olong)(work1[9*nxt+itime]+0.5);
      retCode = ObitTableSNReadRow (SNTab, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Update */
      if (iif==0) {  /* Multiband delays only on first IF */
	row->MBDelay1 = work1[10*nxt+itime];
	if (need2) row->MBDelay2 = work1[13*nxt+itime];
      }

      /* weights zero rather than fblank */
      if (work1[3*nxt+itime]==fblank) work1[3*nxt+itime] = 0.0;
      if (work1[3*nxt+itime] > 0.0) {
	amp = sqrt (work1[0*nxt+itime]*work1[0*nxt+itime] + 
		    work1[1*nxt+itime]*work1[1*nxt+itime]);
	if (amp <= 0.0) amp = 1.0;
	row->Real1[iif]   = work1[0*nxt+itime] * work1[2*nxt+itime] / amp;
	row->Imag1[iif]   = work1[1*nxt+itime] * work1[2*nxt+itime] / amp;
	row->Weight1[iif] = work1[3*nxt+itime];
	if (work1[11*nxt+itime]!=fblank)  row->Delay1[iif]  = work1[11*nxt+itime];
	else row->Delay1[iif] = 0.0;
	if (work1[12*nxt+itime]!=fblank)  row->Rate1[iif]   = work1[12*nxt+itime];
	else row->Rate1[iif] = 0.0;
	(*gncnt) += + 1.0;
	(*gnsum) += work1[2*nxt+itime];
	if (row->RefAnt1[iif] == 0) row->RefAnt1[iif] = refa1;
      } else {  /* Datum bad */
	row->Real1[iif]   = fblank;
	row->Imag1[iif]   = fblank;
	row->Weight1[iif] = 0.0;
	row->Delay1[iif]  = fblank;
	row->Rate1[iif]   = fblank;
      }
      if (need2) {
	/* weights zero rather than fblank */
	if (work1[7*nxt+itime]==fblank) work1[7*nxt+itime] = 0.0;
	if (work1[7*nxt+itime] > 0.0) {
	  amp = sqrt (work1[4*nxt+itime]*work1[4*nxt+itime] + 
		      work1[5*nxt+itime]*work1[5*nxt+itime]);
	  if (amp <= 0.0) amp = 1.0;
	  row->Real2[iif]   = work1[4*nxt+itime] * work1[6*nxt+itime] / amp;
	  row->Imag2[iif]   = work1[5*nxt+itime] * work1[6*nxt+itime] / amp;
	  row->Weight2[iif] = work1[7*nxt+itime];
	  if (work1[14*nxt+itime]!=fblank)  row->Delay2[iif]  = work1[14*nxt+itime];
	  else row->Delay2[iif] = 0.0;
	  if (work1[14*nxt+itime]!=fblank)  row->Rate2[iif]   = work1[15*nxt+itime];
	  else row->Rate2[iif] = 0.0;
	  (*gncnt) += 1.0;
	  (*gnsum) += work1[6*nxt+itime];
	  if (row->RefAnt2[iif] == 0) row->RefAnt2[iif] = refa2;
	} else {  /* Datum bad */
	  row->Real2[iif]   = fblank;
	  row->Imag2[iif]   = fblank;
	  row->Weight2[iif] = 0.0;
	  row->Delay2[iif]  = fblank;
	  row->Rate2[iif]   = fblank;
	}
      }
      
      /* Rewrite record */
      retCode = ObitTableSNWriteRow (SNTab, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    } /* end loop rewriting smoothed solutions L200: */;
    /* First SN number of next antenna */
    
    /* End of antenna loop */
  endAnt: fstrec = save;
  } /* end loop over antennas  L600: */;

  row = ObitTableSNRowUnref(row); /* delete row object */
} /* end of routine ObitUVSolnSNSmooth */ 

/**
 * Does a box car (running mean) smoothing of weighted  
 * irregularly spaced points possibly with blanked values.  Only  
 * returns blanked values if no valid data found.  First good value  
 * used for all previous points, last good value used for all  
 * subsequent points in which all data are blanked in the boxcar.  A  
 * datum is blanked if its weight is <= 0 or fblank.  
 * Routine translated from the AIPSish SMBOX.FOR/SMBOX  
 * \param width   Width of boxcar in same units as X: 0 => replace 
 *                blanks with interpolated closest 2, < 0 => replace 
 *                only blanks with the boxcar smoothed values (all 
 *                others remain unchanged) 
 * \param x       Abscissas of points to be smoothed in increasing 
 *                order 
 * \param y       Values to be smoothed. 
 * \param w       Weights of data. 
 * \param n       Number of points to smooth. 
 * \param ys      Smoothed values. 
 * \param ws      Smoothed weights 
 * \param doBlank replace blanked values with interpolated values?
 */
void 
ObitUVSolnSmooBox (ofloat width, ofloat* x, ofloat* y, ofloat* w, olong n, 
		   ofloat* ys, ofloat* ws, gboolean doBlank) 
{
  olong   i , j, k, l, i1, i2;
  gboolean wasb, onlyb, blnkd;
  ofloat  hw, d, temp, fblank =  ObitMagicF();

  if (n <= 0) return;    /* any data? */
  d = fabs (width);      /* Smoothing width */
  onlyb = width <= 0.0;  /* Only interpolate blanked values? */
  hw = d / 2.0;          /* Half width of smoothing */
  wasb = FALSE;          /* any blanked data? */

  /* width = 0.0 => interp only - here copy good and interpolate later */
  if (d <= 0.0) {
    for (i=0; i<n; i++) { /* loop 10 */
      blnkd = (y[i] == fblank)  ||  (w[i] <= 0.0)  ||  (w[i] == fblank);
      if (blnkd) {
	ys[i] = fblank;
	ws[i] = 0.0;
	wasb = TRUE;
      } else {
	ys[i] = y[i];
	ws[i] = w[i];
      } 
    } /* end loop  L10:  */;

    
  } else {  /* Smooth */
    for (i=0; i<n; i++) { /* loop 100 */
      blnkd = (y[i] == fblank)  ||  (w[i] <= 0.0)  ||  (w[i] == fblank);
      if (blnkd) {
	ys[i] = 0.0;
	ws[i] = 0.0;
	wasb = TRUE;
      } else {
	ys[i] = y[i] * w[i];
	ws[i] = w[i];
      } 

      if ((blnkd)  ||  (!onlyb)) {  /* smoothing data */
	for (k=1; k<=i; k++) { /* weighted sum of previous data loop 20 */
	  j = i - k;
	  if (fabs(x[i]-x[j]) > hw) {
	    break;   /* out of window? goto L25;*/
	  } else {
	    if ((y[j] != fblank)  &&  (w[j] > 0.0)  &&  (w[j] != fblank)) {
	      ys[i] += y[j] * w[j];
	      ws[i] += w[j];
	    } 
	  } 
	} /* end loop  L20:  */;

	for (j= i+1; j<n; j++) { /* weighted sum of subsequent data loop 30 */
	  if (fabs(x[i]-x[j]) > hw) {
	    break;   /* out of window? goto L35;*/  
	  } else {
	    if ((y[j] != fblank)  &&  (w[j] > 0.0)  &&  (w[j] != fblank)) {
	      ys[i] += y[j] * w[j];
	      ws[i] += w[j];
	    } 
	  } 
	} /* end loop  L30:  */;
      } /* end of smoothing */
      
      if (ws[i] > 0.0) {  /*  normalize by sum of weights L35: */
	ys[i] /= ws[i];
      } else {
	wasb = TRUE;
	ys[i] = fblank;
      } 
    } /* end loop  L100: */;
  } /* End of smoothing */


  /* fill in any remaining blanks if needed */
  if (doBlank && wasb) {
    /* extrapolate to ends */
    i1 = n+1;
    i2 = 0;
    for (i=0; i<n; i++) { /* loop 110 */
      if (ws[i] > 0.0) {
	i1 = MIN (i, i1);
	i2 = MAX (i, i2);
      } 
    } /* end loop  L110: */;

    if (i1 > 1) {  /* Fill to beginning */
      j = i1 - 1;
      for (l=0; l<=j; l++) 
	{ys[l] = ys[i1]; ws[l] = ws[i1];}
    } 
    
    if (i2 < n) {  /* Fill to end */
      j = n - i2 - 1;
      for (l=0; l<j; l++) 
	{ys[i2+1+l] = ys[i2]; ws[i2+1+l] = ws[i2];}
    } 

    /* interpolate others */
    for (i=0; i<n; i++) { /* loop 130 */
      if (ws[i] > 0.0) {  /* OK */
	i1 = i;           /* previous valid */
      } else {  /* Blanked, find next valid */
	for (i2=i+1; i2<n; i2++) { /* loop 120 */
	  if (ws[i2] > 0.0) break;  /* found it goto L125;*/
	} /* end loop  L120: */;

	/* interpolate */
	temp = x[i2] - x[i1]; /* L125: */
	if (temp == 0.0) temp = 1.0;
	ys[i] = ys[i1] + (x[i]-x[i1]) * (ys[i2]-ys[i1]) / temp;
	ws[i] = ws[i1] + (x[i]-x[i1]) * (ws[i2]-ws[i1]) / temp;
      } 
    } /* end loop  L130: */
  } /* end of if any blanks to interpolate */ 
} /* end of routine ObitUVSolnSmooBox */ 

/**
 * Does a Gaussian (running mean) smoothing of weighted  
 * irregularly spaced points possibly with blanked values.  Only  
 * returns blanked values if no valid data found.  First good value  
 * used for all previous points, last good value used for all  
 * subsequent points in which all data are blanked in the smoothing interval.  
 * A datum is considered blanked if its weight is <= 0 or its value fblank.  
 * Routine translated from the AIPSish SMGAUS.FOR/SMGAUS  
 * \param width   Width of boxcar in same units as X: 0 => replace 
 *                blanks with interpolated closest 2, < 0 => replace 
 *                only blanks with the b oxcar smoothed values (all 
 *                others remain unchanged) 
 * \param x       Abscissas of points to be smoothed in increasing 
 *                order 
 * \param y       Values to be smoothed. 
 * \param w       Weights of data. 
 * \param n       Number of points to smooth. 
 * \param ys      Smoothed values. 
 * \param ws      Smoothed weights 
 * \param wtsum   scratch 
 * \param doBlank replace blanked values with interpolated values.
 */
void 
ObitUVSolnSmooGauss (ofloat width, ofloat* x, ofloat* y, ofloat* w, olong n, 
		     ofloat* ys, ofloat* ws, ofloat* wtsum, gboolean doBlank) 
{
  olong   i , j, k, l, i1, i2;
  ofloat      hw, d, temp, s=0.0, g, fblank =  ObitMagicF();
  gboolean   wasb, onlyb, blnkd;


  if (n <= 0) return;    /* any data? */
  d = fabs (width);      /* Smoothing width */
  onlyb = width <= 0.0;  /* Only interpolate blanked values? */
  hw = d / 2.0;          /* Half width of smoothing */
  wasb = FALSE;          /* any blanked data? */

  /* width = 0.0 => interp only - here copy good and interpolate later */
  if (d <= 0.0) {
    for (i=0; i<n; i++) { /* loop 10 */
      blnkd = (y[i] == fblank)  ||  (w[i] <= 0.0)  ||  (w[i] == fblank);
      if (blnkd) {
	ys[i] = fblank;
	ws[i] = 0.0;
	wasb = TRUE;
      } else {
	ys[i] = y[i];
	ws[i] = w[i];
      } 
    } /* end loop  L10:  */;

    
  } else {   /* Smooth */
    for (i=0; i<n; i++) { /* loop 100 */
      blnkd = (y[i] == fblank)  ||  (w[i] <= 0.0)  ||  (w[i] == fblank);
      if (blnkd) {
	ys[i] = 0.0;
	ws[i] = 0.0;
	wtsum[i] = 0.0;
	wasb = TRUE;
      } else {
	ys[i] = y[i] * w[i];
	ws[i] = w[i];
	wtsum[i] = 1.0;
      } 

      if ((blnkd)  ||  (!onlyb)) { /* smoothing data */
	for (k=1; k<=i; k++) { /* weighted sum of previous data loop 20 */
	  j = i - k;
	  if (fabs(x[i]-x[j]) > hw) {
	    break;   /* out of window? goto L25; */
	  } else {
	    if ((y[j] != fblank)  &&  (w[j] > 0.0)  &&  (w[j] != fblank)) {
	      g = exp (-((x[i]-x[j])/s)*((x[i]-x[j])/s));
	      wtsum[i] += g;
	      ys[i]    += y[j] * w[j] * g;
	      ws[i]    += w[j] * g;
	    } 
	  } 
	} /* end loop  L20:  */

	for (j= i+1; j<n; j++) { /* weighted sum of subsequent data loop  30*/
	  if (fabs(x[i]-x[j]) > hw) {
	    break;   /* out of window? goto L35;*/
	  } else {
	    if ((y[j] != fblank)  &&  (w[j] > 0.0)  &&  (w[j] != fblank)) {
	      g = exp (-((x[i]-x[j])/s)*((x[i]-x[j])/s));
	      wtsum[i] += g;
	      ys[i]    += y[j] * w[j] * g;
	      ws[i]    += w[j] * g;
	    } 
	  } 
	} /* end loop  L30:  */;
      } /* end of smoothing */ 

      if (ws[i] > 0.0) {   /*  normalize by sum of weights L35: */
	ys[i] = ys[i] / ws[i];
	if (wtsum[i] > 0.0) ws[i] = ws[i] / wtsum[i];
      } else {
	wasb = TRUE;
	ys[i] = fblank;
      } 
    } /* end loop  L100: */;
  } /* end of Smoothing */
 
  /* fill in any remaining blanks if needed  */
  if (doBlank && wasb) {
    /* extrapolate to ends */
    i1 = n+1;
    i2 = 0;
    for (i=0; i<n; i++) { /* loop 110 */
      if (ws[i] > 0.0) {
	i1 = MIN (i, i1);
	i2 = MAX (i, i2);
      } 
    } /* end loop  L110: */;

    if (i1 > 1) {  /* Fill to beginning */
      j = i1 - 1;
      for (l=0; l<=j; l++) 
	{ys[l] = ys[i1]; ws[l] = ws[i1];}
    } 
    
    if (i2 < n) {  /* Fill to end */
      j = n - i2 - 1;
      for (l=0; l<j; l++) 
	{ys[i2+1+l] = ys[i2]; ws[i2+1+l] = ws[i2];}
    } 

    /* interpolate others */
    for (i=0; i<n; i++) { /* loop 130 */
      if (ws[i] > 0.0) {  /* OK */
	i1 = i;           /* previous valid */
      } else {  /* Blanked, find next valid */
	for (i2=i+1; i2<n; i2++) { /* loop 120 */
	  if (ws[i2] > 0.0) break;  /* found it goto L125;*/
	} /* end loop  L120: */;

	/* interpolate */
	temp = x[i2] - x[i1]; /* L125: */
	if (temp == 0.0) temp = 1.0;
	ys[i] = ys[i1] + (x[i]-x[i1]) * (ys[i2]-ys[i1]) / temp;
	ws[i] = ws[i1] + (x[i]-x[i1]) * (ws[i2]-ws[i1]) / temp;
      } 
    } /* end loop  L130: */;
  } /* end of if any blanks to interpolate */ 
} /* end of routine ObitUVSolnSmooGauss */ 

/**
 * Does a median window smoothing of weighted irregularly spaced  
 * points possibly with blanked values.  Only returns blanked values if  
 * no valid data found.  First good value used for all previous points,  
 * last good value used for all subsequent points in which all data are  
 * blanked in the smoothing interval.
 * A datum is considered blanked if its weight is <= 0 or its value fblank.  
 * Routine translated from the AIPSish SMMWF.FOR/SMMWF  
 * \param width   Width of boxcar in same units as X: 0 => replace 
 *                blanks with interpolated closest 2, < 0 => replace 
 *                only blanks with the b oxcar smoothed values (all 
 *                others remain unchanged) 
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \param x       Abscissas of points to be smoothed in increasing 
 *                order 
 * \param y       Values to be smoothed. 
 * \param w       Weights of data. 
 * \param n       Number of points to smooth. 
 * \param ys      Smoothed values. 
 * \param ws      Smoothed weights 
 * \param yor     Scratch 
 * \param wor     Scratch 
 * \param doBlank replace blanked values with interpolated values.
 */
void 
ObitUVSolnSmooMWF (ofloat width, ofloat alpha, ofloat* x, ofloat* y, ofloat* w, olong n, 
		   ofloat* ys, ofloat* ws, ofloat* yor, ofloat* wor, gboolean doBlank) 
{
  olong      i, j, k, l, i1=0, i2=0, ic=0, nword;
  size_t     nbyte;
  ofloat     hw, d, temp, beta=0.0, fblank =  ObitMagicF();
  gboolean   wasb, onlyb, blnkd;

  if (n <= 0) return;    /* any data? */
  d = fabs (width);      /* Smoothing width */
  onlyb = width <= 0.0;  /* Only interpolate blanked values? */
  hw = d / 2.0;          /* Half width of smoothing */
  wasb = FALSE;          /* any blanked data? */
  beta = MAX (0.05, MIN (0.95, alpha)) / 2.0; /*  Average around median factor */

  nbyte = n*sizeof(ofloat);
  /*memmove (yor, y, nbyte);*/
  /*memmove (wor, w, nbyte);*/

  /* width = 0.0 => interp only - here copy good and interpolate later */
  if (d <= 0.0) {
    for (i=0; i<n; i++) { /* loop 10 */
      blnkd = (y[i] == fblank)  ||  (w[i] <= 0.0)  ||  (w[i] == fblank);
      if (blnkd) {
	ys[i] = fblank;
	ws[i] = 0.0;
	wasb = TRUE;
      } else {
	ys[i] = y[i];
	ws[i] = w[i];
      } 
    } /* end loop  L10:  */;
    
  } else {   /* Smooth */
    for (i=0; i<n; i++) { /* loop 100 */
      blnkd = (y[i] == fblank)  ||  (w[i] <= 0.0)  ||  (w[i] == fblank);
      if (blnkd) {
	ys[i] = 0.0;
	ws[i] = 0.0;
	wasb = TRUE;
	ic = 0;
      } else {
	ys[i] = y[i] * w[i];
	ws[i] = w[i];
	ic = 1;
	yor[ic-1] = y[i];
	wor[ic-1] = w[i];
      }
      
      if ((blnkd)  ||  (!onlyb)) { /* smoothing data */
	for (k=1; k<=i; k++) { /* previous values loop 30 */
	  if (fabs(x[i]-x[i-k]) <= hw) {  /* In window? */
	    if ((y[i-k] != fblank)  &&  (w[i-k] > 0.0)  &&  (w[i-k] != fblank)) {
	      /* valid datum in window, order the datum */
	      l = ic;
	      for (j=0; j<ic; j++) { /* find location in ordered arrays for [i-k] loop 15 */
		if (yor[j] > y[i-k]) {l = j; break;}  /* goto L20; */
	      } /* end loop  L15:  */;
	      
	      nword = ic-l+1;
	      nbyte = nword*sizeof(ofloat);
	      if (nword>0) {
		memmove (&yor[l+1], &yor[l], nbyte);
		memmove (&wor[l+1], &wor[l], nbyte);
	      }
	      /* shuffle loop 25 */
	      /*for (j=l; j<=ic; j++) { 
		i1 = ic - j + l;
		yor[i1] = yor[i1-1];
		wor[i1] = wor[i1-1];
		} *//* end loop  L25:  */;
	      yor[l] = y[i-k];  /* insert [i-k] in ordered array */
	      wor[l] = w[i-k];
	      ic++; /* Number of elements in ordered list */
	    } 
	  } 
	} /* end loop  L30:  */;
	
	for (k=i+1; k<n; k++) { /* subsequent points loop 60 */
	  if (fabs(x[i]-x[k]) > hw) {
	    break; /* out of window goto L65;*/
	  } else {
	    if ((y[k] != fblank)  &&  (w[k] > 0.0)  &&  (w[k] != fblank)) {
	      l = ic;
	      for (j=0; j<ic; j++) { /* find location in ordered arrays for [k] loop 45 */
		if (yor[j] > y[k]) {l = j; break;}  /* goto L50; */
	      } /* end loop  L45:  */;

	      nword = ic - l;
	      nbyte = nword*sizeof(ofloat);
	      if (nword>0) {
		memmove (&yor[l+1], &yor[l], nbyte);
		memmove (&wor[l+1], &wor[l], nbyte);
	      }
	      /*shuffle  loop 55 */
	      /*for (j=l; j<ic; j++) { 
		i1 = ic - j + l;
		yor[i1] = yor[i1-1];
		wor[i1] = wor[i1-1];
		}*/ /* end loop  L55:  */;
	      yor[l] = y[k];   /* insert [k] in ordered array */
	      wor[l] = w[k];
	      ic++; /* Number of elements in ordered list */
	    } 
	  } 
	} /* end loop  L60:  */
	
	/* Now average the center set */
	ys[i] = 0.0;  /*  L65: */
	ws[i] = 0.0;
	if (ic > 0) {
	  k = beta * ic + 0.5;
	  i1 = k;
	  i2 = ic - 1 - k;
	  if (i2 < i1) {
	    i1 = MAX (0, i1-1);
	    i2 = MIN (ic-1, i2+1);
	  } 
	  for (k=i1; k<=i2; k++) { /* loop 70 */
	    ys[i] += yor[k] * wor[k];
	    ws[i] += wor[k];
	  } /* end loop  L70:  */;
	} /* end of some data to average */ 
      } /* end  smoothing data */

      /* Get smoothed datum */
      if (ws[i] > 0.0) {
	ys[i] /= ws[i];
	ws[i] /= (i2-i1+1);
      } else {
	ys[i] = fblank;
	wasb = TRUE;
      } 
    } /* end smoothing loop  L100: */
  } /* end of Smooth */

  /* fill in any remaining blanks  if needed */
  if (doBlank && wasb) {
    /* extrapolate to ends */
    i1 = n+1;
    i2 = 0;
    for (i=0; i<n; i++) { /* loop 110 */
      if (ws[i] > 0.0) {
	i1 = MIN (i, i1);
	i2 = MAX (i, i2);
      } 
    } /* end loop  L110: */

    if (i1 > 1) {  /* Fill to beginning */
      j = i1 - 1;
      for (l=0; l<=j; l++) 
	{ys[l] = ys[i1]; ws[l] = ws[i1];}
    } 
    
    if (i2 < n) {  /* Fill to end */
      j = n - i2 - 1;
      for (l=0; l<j; l++) 
	{ys[i2+1+l] = ys[i2]; ws[i2+1+l] = ws[i2];}
    } 

    /* interpolate others */
    for (i=0; i<n; i++) { /* loop 130 */
      if (ws[i] > 0.0) {  /* OK */
	i1 = i;           /* previous valid */
      } else {  /* Blanked, find next valid */
	for (i2=i+1; i2<n; i2++) { /* loop 120 */
	  if (ws[i2] > 0.0) break;  /* found it goto L125;*/
	} /* end loop  L120: */

	/* interpolate */
	temp = x[i2] - x[i1]; /* L125: */
	if (temp == 0.0) temp = 1.0;
	ys[i] = ys[i1] + (x[i]-x[i1]) * (ys[i2]-ys[i1]) / temp;
	ws[i] = ws[i1] + (x[i]-x[i1]) * (ws[i2]-ws[i1]) / temp;
      } 
    } /* end loop  L130: */
  } /* end of if any blanks to interpolate */ 

} /* end of routine ObitUVSolnSmooMWF */ 

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVSolnClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVSolnClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVSolnClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVSolnClassInfoDefFn (gpointer inClass)
{
  ObitUVSolnClassInfo *theClass = (ObitUVSolnClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVSolnClassInit;
  theClass->newObit       = (newObitFP)newObitUVSoln;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVSolnClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVSolnGetClass;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVSolnClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVSolnInit;
  theClass->ObitUVSolnCreate   = (ObitUVSolnCreateFP)ObitUVSolnCreate;
  theClass->ObitUVSolnStartUp  = (ObitUVSolnStartUpFP)ObitUVSolnStartUp;
  theClass->ObitUVSolnGetSN    = (ObitUVSolnGetSNFP)ObitUVSolnGetSN;
  theClass->ObitUVSolnShutDown = (ObitUVSolnShutDownFP)ObitUVSolnShutDown;

} /* end ObitUVSolnClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVSolnInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVSoln *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread        = newObitThread();
  in->info          = newObitInfoList(); 
  in->myUV          = NULL;
  in->SNTable       = NULL;
  in->SNTableRow    = NULL;
  in->CalSel        = NULL;
  in->PriorAntTime  = NULL;
  in->FollowAntTime = NULL;
  in->CalApply      = NULL;
  in->CalPrior      = NULL;
  in->CalFollow     = NULL;
  in->IFR           = NULL;
  in->PriorIFR      = NULL;
  in->FollowIFR     = NULL;
  in->MBDelay       = NULL;
  in->PriorMBDelay  = NULL;
  in->FollowMBDelay = NULL;
  in->RateFact      = NULL;
  in->MissAnt       = NULL;
  in->numRow       = -1;
  in->LastRowRead  = -1;
  in->numAnt       = 0;
  in->numSubA      = 0;
  in->numIF        = 0;
  in->SubA         = 0;
  in->FreqID       = 0;
  in->numPol       = 0;
  in->CurSourID    = -1;
  in->PriorSourID  = -1;
  in->FollowSourID = -1;
  in->PriorCalTime = -1.0;
  in->FollowCalTime= -1.0;
  in->CalTime      = -1.0;
  in->DeltaTime    = -1.0;

} /* end ObitUVSolnInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVSoln* cast to an Obit*.
 */
void ObitUVSolnClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVSoln *in  = inn;
  olong lver;
  gchar *tabtype=NULL;
  ObitErr    *err = NULL;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  /* Smoothed table to be zapped? */
  err = newObitErr ();
  if (in->isSNSmoo) {
    tabtype = g_strdup(in->SNTable->tabType);
    lver    = in->SNTable->tabVer;
    ObitUVZapTable (in->myUV, tabtype, lver, err);
    if (tabtype) g_free (tabtype); tabtype = NULL;
    ObitUVUpdateTables (in->myUV, err);  /* Update disk header */
  }
  ObitErrLog (err);
  err = ObitErrUnref(err);
  in->SNTable       = ObitTableSNUnref(in->SNTable);
  in->SNTableRow    = ObitTableSNRowUnref(in->SNTableRow);
  in->CalSel        = ObitUVSelUnref(in->CalSel);
  in->myUV          = ObitUVUnref(in->myUV);
  if (in->PriorAntTime)  g_free (in->PriorAntTime);  in->PriorAntTime  = NULL;
  if (in->FollowAntTime) g_free (in->FollowAntTime); in->FollowAntTime = NULL;
  if (in->CalApply)      g_free (in->CalApply);      in->CalApply      = NULL;
  if (in->CalPrior)      g_free (in->CalPrior);      in->CalPrior      = NULL;
  if (in->CalFollow)     g_free (in->CalFollow);     in->CalFollow     = NULL;
  if (in->IFR)           g_free (in->IFR);           in->IFR           = NULL;
  if (in->PriorIFR)      g_free (in->PriorIFR);      in->PriorIFR      = NULL;
  if (in->FollowIFR)     g_free (in->FollowIFR);     in->FollowIFR     = NULL;
  if (in->MBDelay)       g_free (in->MBDelay);       in->MBDelay       = NULL;
  if (in->PriorMBDelay)  g_free (in->PriorMBDelay);  in->PriorMBDelay  = NULL;
  if (in->FollowMBDelay) g_free (in->FollowMBDelay); in->FollowMBDelay = NULL;
  if (in->RefAnt)        g_free (in->RefAnt);        in->RefAnt        = NULL;
  if (in->RateFact)      g_free (in->RateFact);      in->RateFact      = NULL;
  if (in->MissAnt)       g_free (in->MissAnt);       in->MissAnt       = NULL;
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVSolnClear */

/**
 * Read calibration for next time from SN table.
 * Applies mean gain modulus corrections
 * CalFollow and CalPrior entries arer:
 * \li amplitude of gain
 * \li phase of gain (rad)
 * \li group delay (sec)
 * \li fringe rate (sec/sec)
 * \li weight
 * \li reference antenna for solution
 * Amplitude, phase, delay, rate, weight, refant
 * Adopted from AIPS CSLGET.FOR
 * \param in   Calibrate Object.
 * \param time desired time in days
 * \param err  Error stack for messages and errors.
 */
static void ObitUVSolnNewTime (ObitUVSoln *in, ofloat time,
			       ObitErr *err)
{
  ObitIOCode retCode;
  ofloat mGModI, wt1, wt2;
  ofloat fblank = ObitMagicF();
  olong nblank, i, j, iant, iif, indx, lenEntryPoln, lenEntry, lenEntryAnt;
  olong  irow, limit,IFoff, antno;
  gboolean want, done, readAll;
  ObitTableSN *SNTable = NULL;
  ObitTableSNRow *SNTableRow = NULL;
  gchar *routine="ObitUVSolnNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* increments, sizes of data elements */
  /* length of basic entry */
  lenEntry = in->numPol * in->lenCalArrayEntry;
  /* length of an entry for an antenna */
  lenEntryAnt = lenEntry * in->numIF;
  /* Length of entry with all polarizations */
  lenEntryPoln = in->lenCalArrayEntry * MIN (2, MAX (1, in->numPol));
  
  /* initialize Prior and Following arrays if first call */
  if (in->LastRowRead <= 0) {
    nblank = in->numAnt * in->numIF * in->numPol * in->lenCalArrayEntry;
    for (i=0; i<nblank; i++) in->CalPrior[i]  = fblank;
    for (i=0; i<nblank; i++) in->CalFollow[i] = fblank;
    for (i=0; i<in->numAnt; i++) in->PriorIFR[i]      = fblank;
    for (i=0; i<in->numAnt; i++) in->FollowIFR[i]     = fblank;
    for (i=0; i<in->numAnt; i++) in->PriorMBDelay[i]  = fblank;
    for (i=0; i<in->numAnt; i++) in->FollowMBDelay[i] = fblank;
    for (i=0; i<in->numAnt; i++) in->PriorMBDelay[i+in->numAnt]  = fblank;
    for (i=0; i<in->numAnt; i++) in->FollowMBDelay[i+in->numAnt] = fblank;
    for (i=0; i<in->numAnt; i++) in->PriorAntTime[i] = -1.0e10;
    for (i=0; i<in->numAnt; i++) in->FollowAntTime[i]= -1.0e10;
    in->PriorCalTime  = -1.0e10;
    in->FollowCalTime = -1.0e10;
    in->PriorSourID  = -1;
    in->FollowSourID = -1;
  } /* end of initialize on first call */
  else {

    /* Shuffle data from Following to Prior if time exceeded */
    for (iant= 0; iant<in->numAnt; iant++) { /* loop 30 */
      /* 2nd time exceeded - shift down */
      if ((time > in->FollowAntTime[iant])  &&  (in->FollowAntTime[iant] > -100.)) {
	in->PriorAntTime[iant] = in->FollowAntTime[iant];
	in->PriorIFR[iant]     = in->FollowIFR[iant];
	in->PriorMBDelay[iant] = in->FollowMBDelay[iant];
	in->PriorMBDelay[iant+in->numAnt]  = in->FollowMBDelay[iant+in->numAnt];
	for (iif=0; iif<in->numIF; iif++) { /* loop 20 */
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif);
	  for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	} /* end IF loop  L20:  */;
      } 
    } /* end loop  ant L30:  */;
  }

  /* SN table  - set local pointers */
  SNTable = (ObitTableSN*)in->SNTable;
  SNTableRow = (ObitTableSNRow*)in->SNTableRow;
  
  /* mean gain modulus correction */
  if (SNTable-> mGMod>0.00001) mGModI = 1.0 / (SNTable-> mGMod);
  else mGModI = 1.0;
  
  /* Read through rows filling in data */
  /* read until selected time. */
  limit = MAX (1, in->LastRowRead);
  readAll = FALSE;
  in->LastRowRead = 0;  /* The next time may start somewhere nonobvious */
  for (i= limit; i<=in->numRow; i++) { /* loop 90 */
    irow = i;
    if (i>=in->numRow) readAll = TRUE;  /* Read whole file? */
    retCode = ObitTableSNReadRow (SNTable, irow, SNTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "Cal(SN) table");
    if (SNTableRow->status < 0) continue; /* entry flagged? */

    /* Check calibrator selector */
    want = ObitUVSelWantSour (in->CalSel, SNTableRow->SourID);
    
    /* check subarray */
    want = want && 
      ((SNTableRow->SubA == in->SubA)  ||  (SNTableRow->SubA <= 0) || (in->SubA <= 0));
    
    /* check frqsel */
    want = want &&
      ((SNTableRow->FreqID == in->FreqID) || (SNTableRow->FreqID <= 0) ||
       (in->FreqID <= 0));
    
    /* skip if not wanted */
    if (!want) continue;
    
    /* antenna number */
    antno = SNTableRow->antNo;
    iant = antno-1;
    
    /* time -> include this one? */
    if (time >= in->FollowAntTime[iant]) { 
      
      if (in->PriorAntTime[iant] > -100.0) {
	/* new following entry - copy to prior */
	in->PriorAntTime[iant] = in->FollowAntTime[iant];
	in->PriorIFR[iant]     = in->FollowIFR[iant];
	in->PriorMBDelay[iant] = in->FollowMBDelay[iant];
	in->PriorMBDelay[iant+in->numAnt]  = in->FollowMBDelay[iant+in->numAnt];
	for (iif=0; iif<in->numIF; iif++) { /* loop 50 */
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif);
	  for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	} /* end IF loop  L50:  */;
      }
      
      /* fill in new following values */
      in->FollowIFR[iant]     = SNTableRow->IFR;
      in->FollowMBDelay[iant] = SNTableRow->MBDelay1;
      if (in->numPol>1) in->FollowMBDelay[iant+in->numAnt]  = SNTableRow->MBDelay2;
      in->FollowAntTime[iant] = SNTableRow->Time;
      
      /* loop over if */
      for (iif=0; iif<in->numIF; iif++) { /* loop 60 */
	IFoff = iif;
	indx = lenEntryAnt * (iant) +  lenEntry * (iif);
	wt1 = SNTableRow->Weight1[IFoff];
	in->CalFollow[indx]   = 
	  sqrt (SNTableRow->Real1[IFoff]*SNTableRow->Real1[IFoff]+
		SNTableRow->Imag1[IFoff]*SNTableRow->Imag1[IFoff]);
	in->CalFollow[indx+1] = atan2 (SNTableRow->Imag1[IFoff], 
				       SNTableRow->Real1[IFoff]+1.0e-20);
	in->CalFollow[indx+2] = SNTableRow->Delay1[IFoff];
	in->CalFollow[indx+3] = SNTableRow->Rate1[IFoff];
	in->CalFollow[indx+4] = SNTableRow->Weight1[IFoff];
	in->CalFollow[indx+5] = SNTableRow->RefAnt1[IFoff];
	if (wt1 <= 0.0) {
	  /* bad calibration entry */
	  in->CalFollow[indx]   = fblank;
	  in->CalFollow[indx+1] = fblank;
	} else {
	  /* mean gain modulus correction to real/imag parts*/
	  if (in->CalFollow[indx]   != fblank) in->CalFollow[indx]   *= mGModI;
	  if (in->CalFollow[indx+1] != fblank) in->CalFollow[indx+1] *= mGModI;
	} 
	
	/* second polarization if present */
	if (in->numPol >= 2) {
	  indx = indx + in->lenCalArrayEntry;
	  wt2 = SNTableRow->Weight2[IFoff];
	  in->CalFollow[indx]   = 
	    sqrt (SNTableRow->Real2[IFoff]*SNTableRow->Real2[IFoff]+
		  SNTableRow->Imag2[IFoff]*SNTableRow->Imag2[IFoff]);
	  in->CalFollow[indx+1] = atan2 (SNTableRow->Imag2[IFoff], 
					 SNTableRow->Real2[IFoff]+1.0e-20);
	  in->CalFollow[indx+2] = SNTableRow->Delay2[IFoff];
	  in->CalFollow[indx+3] = SNTableRow->Rate2[IFoff];
	  in->CalFollow[indx+4] = SNTableRow->Weight2[IFoff];
	  in->CalFollow[indx+5] = SNTableRow->RefAnt2[IFoff];
	  if (wt2 <= 0.0) {
	    /* bad calibration entry */
	    in->CalFollow[indx]   = fblank;
	    in->CalFollow[indx+1] = fblank;
	  } else {
	    /* mean gain modulus correction to real/imag parts*/
	    if (in->CalFollow[indx]   != fblank) in->CalFollow[indx]   *= mGModI;
	    if (in->CalFollow[indx+1] != fblank) in->CalFollow[indx+1] *= mGModI;
	  } 
	} /* end second poln */
      } /* end IF loop  L60:  */;
      
      /* if Prior entry not valid copy following */
      if (in->PriorAntTime[iant] <= -100.) {
	in->PriorAntTime[iant] = in->FollowAntTime[iant];
	in->PriorIFR[iant]     = in->FollowIFR[iant];
	in->PriorMBDelay[iant] = in->FollowMBDelay[iant];
	in->PriorMBDelay[iant+in->numAnt]  = in->FollowMBDelay[iant+in->numAnt];
	for (iif=0; iif<in->numIF; iif++) { /* loop 70 */
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif);
	  for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	} /* end IF loop  L70:  */
      } /* end copy to Prior */
      
    } else {
      
      /* This one not needed - are we there yet? */
      /* May need to restart earlier in table for some antennas.
	 Remember the first record not used. */
      
      if (in->LastRowRead <= 0) in->LastRowRead = i;
      
      /* if desired time for some antenna still after the Following time and there is no
	 Prior entry, keep going */
      done = TRUE;
      for (antno= 1; antno<=in->numAnt; antno++) { /* loop 80 */
	iant = antno-1;
	if (in->MissAnt[iant]) continue;  /* Nothing for this antenna */
	if ((time >= in->FollowAntTime[iant]) && (in->PriorAntTime[iant] >= -100.0)) done=FALSE;
	if (in->FollowAntTime[iant]<-100.0) done = FALSE;
      } /* end loop  L80:  */
      
      /* no more to fill in */
      if (done) break;
    } 
  } /* end loop over table entries L90:  */
  

  /* finished file using all entries? */
  if (in->LastRowRead <= 0) in->LastRowRead = in->numRow + 1;
  
  /* Set times */
  in->FollowCalTime = 1.0e10;
  in->PriorCalTime = -1.0e10;
  for (antno= 1; antno<=in->numAnt; antno++) { /* loop 110 */
    iant = antno -1;
    /* If got to the end and antenna has no entries - mark it */
    if (readAll && 
	((in->PriorAntTime[iant] <-100.0) && (in->FollowAntTime[iant]<-100.0)))
      in->MissAnt[iant] = TRUE;

    if (in->PriorAntTime[iant] >= -100.0) {
      if (time >= in->PriorAntTime[iant]) 
	in->PriorCalTime = MAX (in->PriorCalTime, in->PriorAntTime[iant]);
      if (time <= in->FollowAntTime[iant]) 
	in->FollowCalTime = MIN (in->FollowCalTime,in->FollowAntTime[iant]);
    } 
  } /* end loop  L110: */
  
  /* just to be sure something rational in times */
  if (in->PriorCalTime < -1000.0)  in->PriorCalTime  = time - 2.0/86400.0;
  if (in->FollowCalTime > 10000.0) in->FollowCalTime = time + 2.0/86400.0;
  
} /* end ObitUVSolnNewTime */

/**
 * Routine to determine antenna usage as the reference antenna in an  
 * open SN table.  All IFs are examined.  
 * All valid entries are included.  
 * Routine translated from the AIPSish REFCNT.FOR/REFCNT  
 * \param SNTab    SN table object 
 * \param isub     Subarray number, 0=>1 
 * \param numtime  [out] number of distinct times in table. 
 * \param err      Error/message stack, returns if error.
 * \return pointer to array of counts of usage as reference ant.  
 *         g_free when done.
 */
static olong *refCount (ObitTableSN *SNTab, olong isub, 
		       olong *numtime, ObitErr* err) 
{
  olong   *antuse=NULL;
  ObitIOCode retCode;
  olong   iif, iref, numif, numpol, numant=0,sub;
  olong  loop, count;
  ofloat lastTime;
  ObitTableSNRow *row=NULL;
  gchar *routine = "refCount";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return antuse;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));

  /* Get descriptive info */
  numpol = SNTab->numPol;
  numif  = SNTab->numIF;
  numant = SNTab->numAnt;

  /* create output array */
  antuse = g_malloc0(numant*sizeof(olong));
  Obit_retval_if_fail((antuse!=NULL), err, antuse,
    "%s: Fail to allocate antuse", routine);

  /* Subarray */
  sub = MAX (1, isub);
  lastTime = -1.0e20;
  count = 0;
  
  /* Create Row */
  row = newObitTableSNRow (SNTab);
  /* Loop through table */
  for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop 20 */

    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) break;
    if (row->status<0) continue;  /* Skip deselected record */

    /* Count times - only allow epsilon time difference */
    if (row->Time>(lastTime+0.0005*row->TimeI)) {
      lastTime = row->Time;
      count ++;
    }

    /* Right subarray? */
    if ((row->SubA!=sub) && (row->SubA>0)) continue;
    /* Loop over IFs */
    for (iif=0; iif<numif; iif++) { /* loop 10 */
      if (row->Weight1[iif] > 0.0) {
	iref =  row->RefAnt1[iif];
	if ((iref > 0)  &&  (iref <= numant)) antuse[iref-1]++;
      } 
      if ((numpol>1) && (row->Weight2[iif] > 0.0)) {
	iref =  row->RefAnt2[iif];
	if ((iref > 0)  &&  (iref <= numant)) antuse[iref-1]++;
      }
    } /* end loop  L10:  */
  } /* end loop  L20:  */

  row = ObitTableSNRowUnref(row); /* delete row object */
  if (err->error) Obit_traceback_val (err, routine, SNTab->name, antuse);

  /* How many times? */
  *numtime = count;
  return antuse;
} /* end of routine refCount */ 

/**
 * Routine to rereference phases in an a polarization coherent fashion;  
 * i.e. both polarizations must be present (if possible) in data used  
 * to determine the relative phase between the primary and secondary  
 * reference antennas.  This routine does one IF at a time but both  
 * polarizations (if present) are done.  
 * All valid entries are included.  
 * Note: reference antenna values in the table are modified.  
 * Routine translated from the AIPSish REFFAZ.FOR/REFFAZ  
 * \param SNTab    SN table object, should be opened and closed externally
 * \param isub     subarray number, 0=>1 
 * \param iif      0-rel IF number, 0=>1 
 * \param refa     primary reference antenna 
 * \param ant      secondary reference antenna 
 * \param mxtime   Dimension of work arrays 
 * \param wrktim   Work array 
 * \param work1    Work array 
 * \param work2    Work array 
 * \param work3    Work array 
 * \param work4    Work array 
 * \param work5    Work array 
 * \param err      Error/message stack, returns if error.
 */
static void 
refPhase (ObitTableSN *SNTab, olong isub, olong iif, olong refa, olong ant, 
	  olong mxtime, ofloat* wrktim, ofloat* work1, ofloat* work2, 
	  ofloat* work3, ofloat* work4, ofloat* work5, ObitErr* err)
{
  olong    numtim, numif, numpol, ipnt1, ipnt2, refa1, refa2, isuba, sa;
  olong   loop, numrec;
  gboolean   need2, dotwo, done;
  ofloat wt1=0.0, wt2, re1, re2, im1, im2, amp, tre, tim, smotim;
  ofloat fblank =  ObitMagicF();
  odouble timoff=0.0, time, time1, time2;
  ObitIOCode retCode;
  ObitTableSNRow *row=NULL;
  gchar *routine = "refPhase";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
 
  /* Subarray */
  isuba = MAX (1, isub);
  
  /* Get descriptive info */
  numif  = SNTab->numIF;
  numpol = SNTab->numPol;

  /* Create Row */
  row = newObitTableSNRow (SNTab);

  /* Initially require both polarizations if present */
  need2 =  numpol>1; 
  dotwo =  numpol>1; 

  /* Loop thru table referring ant to refa. */
  done = FALSE;
  while (!done) { /* Loop trying combined and separate poln L10: */
    numtim = 0;
    for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop 100 */

      retCode = ObitTableSNReadRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Is this either of the antennas we're looking for? */
      if ((row->antNo!=refa) && (row->antNo!=ant)) continue;  /* Keep truckin' */
      
      /* right subarray? */
      if ((row->SubA!=isuba) && (row->SubA>0))  continue;  /* Keep truckin' */
      
      /* Find and check reference  antennas.  Must all be the  same. */
      refa1 = row->RefAnt1[iif];
      if (need2) refa2 = row->RefAnt2[iif];
      else refa2 = 0;
      
      /* Bad solution? */
      if (((row->Weight1[iif] <= 0.0)  ||  ( row->Real1[iif] == fblank)  ||  
	   (row->RefAnt1[iif] <= 0))) continue;  /* goto L100;*/
      if (need2  &&  ((numpol > 1)  &&  
		      ((row->Weight2[iif] <= 0.0) || (row->Real2[iif] == fblank)  || 
		       (row->RefAnt2[iif] <= 0)))) continue;  /* goto L100; */
      
      /* Desired antenna combination? */
      if (need2  &&  (refa1 > 0)  &&  (refa2 > 0)  &&  (refa1 != refa2)) continue; 
      
      if (refa1 < 0) refa1 = refa2;
      if ((refa1 != ant)  &&  (refa1 != refa)) continue;  /* goto L100;*/
      if (row->antNo == refa1) continue;   /* goto L100;*/
      if (numtim >= mxtime) continue;   /* goto L100;*/
      
      /* Save times */
      if (numtim == 0) timoff = row->Time;
      wrktim[numtim] = row->Time - timoff;
      
      if (refa1 != ant) {
	/* refa is reference ant */
	work2[numtim] = row->Real1[iif];
	work3[numtim] = row->Imag1[iif];
	if (dotwo) {
	  work4[numtim] = row->Real2[iif];
	  work5[numtim] = row->Imag2[iif];
	} 
      } else {
	/* ant is reference ant */
	work2[numtim] =  row->Real1[iif];
	work3[numtim] = -row->Imag1[iif];
	if (dotwo) {
	  work4[numtim] =  row->Real2[iif];
	  work5[numtim] = -row->Imag2[iif];
	} 
      } 
      numtim++;  /* count times */
    } /* end loop  L100: */;
    
    if (need2  &&  (numtim <= 0)) {
      /* Try again with only one poln. */
      need2 = FALSE;
      done  =  FALSE;  /* Need another pass to get both poln */
    } else done = TRUE;

  } /* End of loop over poln L10: */

  /* Find any? */
  if (numtim <= 0) {ObitTableSNRowUnref(row); return;}

  /* Smooth (2 sec to extrapolate) */
  smotim = 2.0 / 86400.0;
  boxsmo (smotim, wrktim, work2, numtim, work1);
  boxsmo (smotim, wrktim, work3, numtim, work2);
  if (dotwo) {
    boxsmo (smotim, wrktim, work4, numtim, work3);
    boxsmo (smotim, wrktim, work5, numtim, work4);
  } 

  /* Set up for interpolation */
  ipnt1 = 0;
  ipnt2 = 1;
  time1 = wrktim[0];
  time2 = wrktim[1];
  if (numtim == 1) {  /* Only one entry in array to interpolate? */
    ipnt2 = 0;
    time2 = time1;
  } 

  /* Loop thru table changing any data with ref=ant to ref=refa */
  numrec = SNTab->myDesc->nrow;
  for (loop=1; loop<=numrec; loop++) { /* loop 200 */
    
    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    if (row->status<0) continue;  /* Skip deselected record */

    /* right antenna? */
    if ((row->RefAnt1[iif] != ant)  &&  
	(dotwo  &&  (row->RefAnt2[iif] != ant))) continue;  /*goto L200;*/
    
    /* Right subarray? */
    sa = row->SubA;
    
    /* Interpolate */
    if ((sa == isuba)  ||  (sa <= 0)) {
      time = row->Time - timoff;
      done = FALSE;
      while (!done) {  /* loop until interpolation setup */
	if ((time >= time1)  &&  (time < time2)) {  /* L140: */
	  /* Between entries */
	  if (time2 != time1) {
	    wt1 = 1.0 - ((time-time1) / (time2-time1));
	  } else {
	    wt1 = 1.0;
	  } 
	  done = TRUE; 
	  
	} else if (time < time1) {/* Before first time */
	  wt1 = 1.0;
	  done = TRUE; 
	  
	} else if ((ipnt2+1) >= numtim) { /* After last time */
	  wt1 = 0.0;
	  done = TRUE; 
	  
	} else {   /* Shift in interpolation arrays */
	  ipnt1++;
	  ipnt2 = ipnt1+1;
	  time1 = wrktim[ipnt1];
	  time2 = wrktim[ipnt2];
	  done = FALSE;    /* goto L140; */
	}
      } /* end of loop until interpolation setup */
      
      /* Interpolate */
      wt2 = 1.0 - wt1;
      if (row->RefAnt1[iif] == ant) {
	/* Interpolate phase pol 1 */
	if ((row->Imag1[iif] != fblank)  &&  (row->Real1[iif] != fblank)) {
	  re1 = wt1 * work1[ipnt1] + wt2 * work1[ipnt2];
	  im1 = wt1 * work2[ipnt1] + wt2 * work2[ipnt2];
	  /* Normalize by amplitude */
	  amp = MAX (sqrt (re1*re1 + im1*im1), 1.0e-10);
	  re1 /= amp;
	  im1 /= amp;

	  /* Correct phase pol 1 */
	  tre = row->Real1[iif];
	  tim = row->Imag1[iif];
	  row->Real1[iif] = tre*re1 - tim*im1;
	  row->Imag1[iif] = tre*im1 + tim*re1;
	}  /* end data valid */
	
	/* Relabel reference antenna */
	row->RefAnt1[iif] = refa;
	
      } else if (row->RefAnt1[iif] == 0) {  /* Null ref ant */
	/* Relabel if blanked */
	if ((row->Imag1[iif] == fblank)  ||  (row->Real1[iif] == fblank)) 
	  row->RefAnt1[iif] = refa;
      }
      
      /* Second polarization */
      if (dotwo  &&  (row->RefAnt2[iif] == ant)) {
	if ((row->Imag2[iif] != fblank)  &&  (row->Real2[iif] != fblank)) {
	  re2 = wt1 * work3[ipnt1] + wt2 * work3[ipnt2];
	  im2 = wt1 * work4[ipnt1] + wt2 * work4[ipnt2];
	  /* Normalize by amplitude */
	  amp = MAX (sqrt (re2*re2 + im2*im2), 1.0e-10);
	  re2 /= amp;
	  im2 /= amp;
	  
	  /* Correct phase. pol 2 */
	  tre = row->Real2[iif];
	  tim = row->Imag2[iif];
	  row->Real2[iif] = tre*re2 - tim*im2;
	  row->Imag2[iif] = tre*im2 + tim*re2;
	} /* end data valid */
	
	/* Relabel reference antenna */
	row->RefAnt2[iif] = refa;
	
      } else if ((numpol > 1) && (row->RefAnt2[iif] == 0)) { /* Null ref ant */
	/* Relabel if blanked */
	if ((row->Imag2[iif] == fblank)  ||  (row->Real2[iif] == fblank)) 
	  row->RefAnt2[iif] = refa;
      }
      
      /* Rewrite record */
      retCode = ObitTableSNWriteRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    } /* end of if correct subarray */
  } /* end loop  L200: */;

  row = ObitTableSNRowUnref(row); /* delete row object */
} /* end of routine refPhase */ 

/**
 * Routine to rereference SB delay.
 * Require both polarizations if available and possible to maintain 
 * Polarization coherence.
 * All valid entries are included.  
 * Note: reference antenna values in the table are modified.  
 * Routine adopted from the AIPSish REFFAZ.FOR/REFFAZ  
 * \param SNTab    SN table object, should be opened and closed externally
 * \param isub     subarray number, 0=>1 
 * \param refa     primary reference antenna 
 * \param ant      secondary reference antenna 
 * \param mxtime   Dimension of work arrays 
 * \param wrktim   Work array 
 * \param work1    Work array 
 * \param work2    Work array 
 * \param work4    Work array 
 * \param err      Error/message stack, returns if error.
 */
static void 
refMBDelay (ObitTableSN *SNTab, olong isub, olong refa, olong ant, 
	  olong mxtime, ofloat* wrktim, ofloat* work1, ofloat* work2, 
	  ofloat* work4, ObitErr* err)
{
  olong    numtim, numif, numpol, ipnt1, ipnt2, refa1, refa2, isuba, sa;
  olong   loop, numrec;
  gboolean   need2, dotwo, done, gotAny=FALSE;
  ofloat wt1=0.0, wt2, de1, de2, smotim;
  ofloat fblank =  ObitMagicF();
  odouble timoff=0.0, time, time1, time2;
  ObitIOCode retCode;
  ObitTableSNRow *row=NULL;
  gchar *routine = "refMBDelay";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
 
  /* Subarray */
  isuba = MAX (1, isub);
  
  /* Get descriptive info */
  numif  = SNTab->numIF;
  numpol = SNTab->numPol;

  /* Create Row */
  row = newObitTableSNRow (SNTab);

  /* Initially require both polarizations if present */
  need2 =  numpol>1; 
  dotwo =  numpol>1; 

  /* Loop thru table referring ant to refa. */
  done = FALSE;
  while (!done) { /* Loop trying combined and separate poln L10: */
    numtim = 0;
    for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop 100 */

      retCode = ObitTableSNReadRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Is this either of the antennas we're looking for? */
      if ((row->antNo!=refa) && (row->antNo!=ant)) continue;  /* Keep truckin' */
      
      /* right subarray? */
      if ((row->SubA!=isuba) && (row->SubA>0))  continue;  /* Keep truckin' */
      
      /* Find and check reference  antennas.  Must all be the  same. */
      refa1 = row->RefAnt1[0];
      if (need2) refa2 = row->RefAnt2[0];
      else refa2 = 0;
      
      /* Bad solution? */
      if (((row->Weight1[0] <= 0.0)  ||  ( row->MBDelay1 == fblank)  ||  
	   (row->RefAnt1[0] <= 0))) continue;  /* goto L100;*/
      if (need2  &&  ((numpol > 1)  &&  
		      ((row->Weight2[0] <= 0.0) || (row->MBDelay2 == fblank)  || 
		       (row->RefAnt2[0] <= 0)))) continue;  /* goto L100; */
      
      /* Desired antenna combination? */
      if (need2  &&  (refa1 > 0)  &&  (refa2 > 0)  &&  (refa1 != refa2)) continue; 
      
      if (refa1 < 0) refa1 = refa2;
      if ((refa1 != ant)  &&  (refa1 != refa)) continue;  /* goto L100;*/
      if (row->antNo == refa1) continue;   /* goto L100;*/
      if (numtim >= mxtime) continue;   /* goto L100;*/
      
      /* Save times */
      if (numtim == 0) timoff = row->Time;
      wrktim[numtim] = row->Time - timoff;
      
      /* Any nonzero data? */
      if ((row->MBDelay1!=0.0) && (row->MBDelay1!=fblank)) gotAny = TRUE;
      if (dotwo && (row->MBDelay1!=0.0) && (row->MBDelay1!=fblank)) gotAny = TRUE;
      if (refa1 != ant) {
	/* refa is reference ant */
	work2[numtim] = row->MBDelay1;
	if (dotwo) work4[numtim] = row->MBDelay2;
      } else {
	/* ant is reference ant */
	work2[numtim] = -row->MBDelay1;
	if (dotwo) work4[numtim] = -row->MBDelay2;
      } 
      numtim++;  /* count times */
    } /* end loop  L100: */;
    
    if (need2  &&  (numtim <= 0)) {
      /* Try again with only one poln. */
      need2 = FALSE;
      done  = FALSE;  /* Need another pass to get both poln */
    } else done = TRUE;

  } /* End of loop over poln L10: */

  /* Find any? */
  if (numtim <= 0) {ObitTableSNRowUnref(row); return;}
  if (!gotAny)  {ObitTableSNRowUnref(row); return;}

  /* Smooth (2 sec to extrapolate) */
  smotim = 2.0 / 86400.0;
  boxsmo (smotim, wrktim, work2, numtim, work1);
  if (dotwo) boxsmo (smotim, wrktim, work4, numtim, work2);

  /* Set up for interpolation */
  ipnt1 = 0;
  ipnt2 = 1;
  time1 = wrktim[0];
  time2 = wrktim[1];
  if (numtim == 1) {  /* Only one entry in array to interpolate? */
    ipnt2 = 0;
    time2 = time1;
  } 

  /* Loop thru table changing any data with ref=ant to ref=refa */
  numrec = SNTab->myDesc->nrow;
  for (loop=1; loop<=numrec; loop++) { /* loop 200 */
    
    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    if (row->status<0) continue;  /* Skip deselected record */
    
    /* right antenna? */
    if ((row->RefAnt1[0] != ant)  &&  
	(dotwo  &&  (row->RefAnt2[0] != ant))) continue;  /*goto L200;*/
    
    /* Right subarray? */
    sa = row->SubA;
    
    /* Interpolate */
    if ((sa == isuba)  ||  (sa <= 0)) {
      time = row->Time - timoff;
      done = FALSE;
      while (!done) {  /* loop until interpolation setup */
	if ((time >= time1)  &&  (time < time2)) {  /* L140: */
	  /* Between entries */
	  if (time2 != time1) {
	    wt1 = 1.0 - ((time-time1) / (time2-time1));
	  } else {
	    wt1 = 1.0;
	  } 
	  done = TRUE; 
	  
	} else if (time < time1) {/* Before first time */
	  wt1 = 1.0;
	  done = TRUE; 
	  
	} else if ((ipnt2+1) >= numtim) { /* After last time */
	  wt1 = 0.0;
	  done = TRUE; 
	  
	} else {   /* Shift in interpolation arrays */
	  ipnt1++;
	  ipnt2 = ipnt1+1;
	  time1 = wrktim[ipnt1];
	  time2 = wrktim[ipnt2];
	  done = FALSE;    /* goto L140; */
	}
      } /* end of loop until interpolation setup */
      
      /* Interpolate */
      wt2 = 1.0 - wt1;
      if (row->RefAnt1[0] == ant) {
	/* Interpolate  pol 1 */
	if (row->MBDelay1 != fblank) {
	  de1 = wt1 * work1[ipnt1] + wt2 * work1[ipnt2];
	  /* Correct  pol 1 */
	  row->MBDelay1 -= de1;
	}  /* end data valid */
      }
      
      /* Second polarization */
      if (dotwo  &&  (row->RefAnt2[0] == ant)) {
	if (row->MBDelay2 != fblank) {
	  de2 = wt1 * work2[ipnt1] + wt2 * work2[ipnt2];
	  /* Correct pol 2 */
	  row->MBDelay2 -= de2;
	} /* end data valid */
      }
      
      /* Rewrite record */
      retCode = ObitTableSNWriteRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    } /* end of if correct subarray */
  } /* end loop  L200: */;

  row = ObitTableSNRowUnref(row); /* delete row object */
} /* end of routine refMBDelay */ 

/**
 * Routine to rereference SB delay.
 * Require both polarizations if available and possible to maintain 
 * Polarization coherence.
 * All valid entries are included.  
 * Note: reference antenna values in the table are modified.  
 * Routine adopted from the AIPSish REFFAZ.FOR/REFFAZ  
 * \param SNTab    SN table object, should be opened and closed externally
 * \param isub     subarray number, 0=>1 
 * \param iif      0-rel IF number, 0=>1 
 * \param refa     primary reference antenna 
 * \param ant      secondary reference antenna 
 * \param mxtime   Dimension of work arrays 
 * \param wrktim   Work array 
 * \param work1    Work array 
 * \param work2    Work array 
 * \param work4    Work array 
 * \param err      Error/message stack, returns if error.
 */
static void 
refDelay (ObitTableSN *SNTab, olong isub, olong iif, olong refa, olong ant, 
	  olong mxtime, ofloat* wrktim, ofloat* work1, ofloat* work2, 
	  ofloat* work4, ObitErr* err)
{
  olong    numtim, numif, numpol, ipnt1, ipnt2, refa1, refa2, isuba, sa;
  olong   loop, numrec;
  gboolean   need2, dotwo, done, gotAny=FALSE;
  ofloat wt1=0.0, wt2, de1, de2, smotim;
  ofloat fblank =  ObitMagicF();
  odouble timoff=0.0, time, time1, time2;
  ObitIOCode retCode;
  ObitTableSNRow *row=NULL;
  gchar *routine = "refDelay";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
 
  /* Subarray */
  isuba = MAX (1, isub);
  
  /* Get descriptive info */
  numif  = SNTab->numIF;
  numpol = SNTab->numPol;

  /* Create Row */
  row = newObitTableSNRow (SNTab);

  /* Initially require both polarizations if present */
  need2 =  numpol>1; 
  dotwo =  numpol>1; 

  /* Loop thru table referring ant to refa. */
  done = FALSE;
  while (!done) { /* Loop trying combined and separate poln L10: */
    numtim = 0;
    for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop 100 */

      retCode = ObitTableSNReadRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Is this either of the antennas we're looking for? */
      if ((row->antNo!=refa) && (row->antNo!=ant)) continue;  /* Keep truckin' */
      
      /* right subarray? */
      if ((row->SubA!=isuba) && (row->SubA>0))  continue;  /* Keep truckin' */
      
      /* Find and check reference  antennas.  Must all be the  same. */
      refa1 = row->RefAnt1[iif];
      if (need2) refa2 = row->RefAnt2[iif];
      else refa2 = 0;
      
      /* Bad solution? */
      if (((row->Weight1[iif] <= 0.0)  ||  ( row->Delay1[iif] == fblank)  ||  
	   (row->RefAnt1[iif] <= 0))) continue;  /* goto L100;*/
      if (need2  &&  ((numpol > 1)  &&  
		      ((row->Weight2[iif] <= 0.0) || (row->Delay2[iif] == fblank)  || 
		       (row->RefAnt2[iif] <= 0)))) continue;  /* goto L100; */
      
      /* Desired antenna combination? */
      if (need2  &&  (refa1 > 0)  &&  (refa2 > 0)  &&  (refa1 != refa2)) continue; 
      
      if (refa1 < 0) refa1 = refa2;
      if ((refa1 != ant)  &&  (refa1 != refa)) continue;  /* goto L100;*/
      if (row->antNo == refa1) continue;   /* goto L100;*/
      if (numtim >= mxtime) continue;   /* goto L100;*/
      
      /* Save times */
      if (numtim == 0) timoff = row->Time;
      wrktim[numtim] = row->Time - timoff;
      
      /* Any nonzero data? */
      if ((row->Delay1[iif]!=0.0) && (row->Delay1[iif]!=fblank)) gotAny = TRUE;
      if (dotwo && (row->Delay1[iif]!=0.0) && (row->Delay1[iif]!=fblank)) gotAny = TRUE;
      if (refa1 != ant) {
	/* refa is reference ant */
	work2[numtim] = row->Delay1[iif];
	if (dotwo) work4[numtim] = row->Delay2[iif];
      } else {
	/* ant is reference ant */
	work2[numtim] = -row->Delay1[iif];
	if (dotwo) work4[numtim] = -row->Delay2[iif];
      } 
      numtim++;  /* count times */
    } /* end loop  L100: */;
    
    if (need2  &&  (numtim <= 0)) {
      /* Try again with only one poln. */
      need2 = FALSE;
      done  = FALSE;  /* Need another pass to get both poln */
    } else done = TRUE;

  } /* End of loop over poln L10: */

  /* Find any? */
  if (numtim <= 0) {ObitTableSNRowUnref(row); return;}
  if (!gotAny)  {ObitTableSNRowUnref(row); return;}

  /* Smooth (2 sec to extrapolate) */
  smotim = 2.0 / 86400.0;
  boxsmo (smotim, wrktim, work2, numtim, work1);
  if (dotwo) boxsmo (smotim, wrktim, work4, numtim, work2);

  /* Set up for interpolation */
  ipnt1 = 0;
  ipnt2 = 1;
  time1 = wrktim[0];
  time2 = wrktim[1];
  if (numtim == 1) {  /* Only one entry in array to interpolate? */
    ipnt2 = 0;
    time2 = time1;
  } 

  /* Loop thru table changing any data with ref=ant to ref=refa */
  numrec = SNTab->myDesc->nrow;
  for (loop=1; loop<=numrec; loop++) { /* loop 200 */
    
    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    if (row->status<0) continue;  /* Skip deselected record */
    
    /* right antenna? */
    if ((row->RefAnt1[iif] != ant)  &&  
	(dotwo  &&  (row->RefAnt2[iif] != ant))) continue;  /*goto L200;*/
    
    /* Right subarray? */
    sa = row->SubA;
    
    /* Interpolate */
    if ((sa == isuba)  ||  (sa <= 0)) {
      time = row->Time - timoff;
      done = FALSE;
      while (!done) {  /* loop until interpolation setup */
	if ((time >= time1)  &&  (time < time2)) {  /* L140: */
	  /* Between entries */
	  if (time2 != time1) {
	    wt1 = 1.0 - ((time-time1) / (time2-time1));
	  } else {
	    wt1 = 1.0;
	  } 
	  done = TRUE; 
	  
	} else if (time < time1) {/* Before first time */
	  wt1 = 1.0;
	  done = TRUE; 
	  
	} else if ((ipnt2+1) >= numtim) { /* After last time */
	  wt1 = 0.0;
	  done = TRUE; 
	  
	} else {   /* Shift in interpolation arrays */
	  ipnt1++;
	  ipnt2 = ipnt1+1;
	  time1 = wrktim[ipnt1];
	  time2 = wrktim[ipnt2];
	  done = FALSE;    /* goto L140; */
	}
      } /* end of loop until interpolation setup */
      
      /* Interpolate */
      wt2 = 1.0 - wt1;
      if (row->RefAnt1[iif] == ant) {
	/* Interpolate  pol 1 */
	if (row->Delay1[iif] != fblank) {
	  de1 = wt1 * work1[ipnt1] + wt2 * work1[ipnt2];
	  /* Correct  pol 1 */
	  row->Delay1[iif] -= de1;
	}  /* end data valid */
      }
      
      /* Second polarization */
      if (dotwo  &&  (row->RefAnt2[iif] == ant)) {
	if (row->Delay2[iif] != fblank) {
	  de2 = wt1 * work2[ipnt1] + wt2 * work2[ipnt2];

	  /* Correct pol 2 */
	  row->Delay2[iif] -= de2;
	} /* end data valid */
      }
      
      /* Rewrite record */
      retCode = ObitTableSNWriteRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    } /* end of if correct subarray */
  } /* end loop  L200: */;

  row = ObitTableSNRowUnref(row); /* delete row object */
} /* end of routine refDelay */ 

/**
 * Routine to rereference rates
 * Require both polarizations if available and possible to maintain 
 * Polarization coherence.
 * All valid entries are included.  
 * Note: reference antenna values in the table are modified.  
 * Routine adopted from the AIPSish REFFAZ.FOR/REFFAZ  
 * \param SNTab    SN table object, should be opened and closed externally
 * \param isub     subarray number, 0=>1 
 * \param iif      0-rel IF number, 0=>1 
 * \param refa     primary reference antenna 
 * \param ant      secondary reference antenna 
 * \param mxtime   Dimension of work arrays 
 * \param wrktim   Work array 
 * \param work1    Work array 
 * \param work2    Work array 
 * \param work4    Work array 
 * \param err      Error/message stack, returns if error.
 */
static void 
refRate (ObitTableSN *SNTab, olong isub, olong iif, olong refa, olong ant, 
	  olong mxtime, ofloat* wrktim, ofloat* work1, ofloat* work2, 
	  ofloat* work4, ObitErr* err)
{
  olong    numtim, numif, numpol, ipnt1, ipnt2, refa1, refa2, isuba, sa;
  olong   loop, numrec;
  gboolean   need2, dotwo, done, gotAny=FALSE;
  ofloat wt1=0.0, wt2, de1, de2, smotim;
  ofloat fblank =  ObitMagicF();
  odouble timoff=0.0, time, time1, time2;
  ObitIOCode retCode;
  ObitTableSNRow *row=NULL;
  gchar *routine = "refRate";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));
 
  /* Subarray */
  isuba = MAX (1, isub);
  
  /* Get descriptive info */
  numif  = SNTab->numIF;
  numpol = SNTab->numPol;

  /* Create Row */
  row = newObitTableSNRow (SNTab);

  /* Initially require both polarizations if present */
  need2 =  numpol>1; 
  dotwo =  numpol>1; 

  /* Loop thru table referring ant to refa. */
  done = FALSE;
  while (!done) { /* Loop trying combined and separate poln L10: */
    numtim = 0;
    for (loop=1; loop<=SNTab->myDesc->nrow; loop++) { /* loop 100 */

      retCode = ObitTableSNReadRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Is this either of the antennas we're looking for? */
      if ((row->antNo!=refa) && (row->antNo!=ant)) continue;  /* Keep truckin' */
      
      /* right subarray? */
      if ((row->SubA!=isuba) && (row->SubA>0))  continue;  /* Keep truckin' */
      
      /* Find and check reference  antennas.  Must all be the  same. */
      refa1 = row->RefAnt1[iif];
      if (need2) refa2 = row->RefAnt2[iif];
      else refa2 = 0;
      
      /* Bad solution? */
      if (((row->Weight1[iif] <= 0.0)  ||  ( row->Rate1[iif] == fblank)  ||  
	   (row->RefAnt1[iif] <= 0))) continue;  /* goto L100;*/
      if (need2  &&  ((numpol > 1)  &&  
		      ((row->Weight2[iif] <= 0.0) || (row->Rate2[iif] == fblank)  || 
		       (row->RefAnt2[iif] <= 0)))) continue;  /* goto L100; */
      
      /* Desired antenna combination? */
      if (need2  &&  (refa1 > 0)  &&  (refa2 > 0)  &&  (refa1 != refa2)) continue; 
      
      if (refa1 < 0) refa1 = refa2;
      if ((refa1 != ant)  &&  (refa1 != refa)) continue;  /* goto L100;*/
      if (row->antNo == refa1) continue;   /* goto L100;*/
      if (numtim >= mxtime) continue;   /* goto L100;*/
      
      /* Save times */
      if (numtim == 0) timoff = row->Time;
      wrktim[numtim] = row->Time - timoff;
      
      /* Any nonzero data? */
      if ((row->Rate1[iif]!=0.0) && (row->Rate1[iif]!=fblank)) gotAny = TRUE;
      if (dotwo && (row->Rate1[iif]!=0.0) && (row->Rate1[iif]!=fblank)) gotAny = TRUE;
      if (refa1 != ant) {
	/* refa is reference ant */
	work2[numtim] = row->Rate1[iif];
	if (dotwo) work4[numtim] = row->Rate2[iif];
      } else {
	/* ant is reference ant */
	work2[numtim] = -row->Rate1[iif];
	if (dotwo) work4[numtim] = -row->Rate2[iif];
      } 
      numtim++;  /* count times */
    } /* end loop  L100: */;
    
    if (need2  &&  (numtim <= 0)) {
      /* Try again with only one poln. */
      need2 = FALSE;
      done  = FALSE;  /* Need another pass to get both poln */
    } else done = TRUE;

  } /* End of loop over poln L10: */

  /* Find any? */
  if (numtim <= 0) {ObitTableSNRowUnref(row); return;}
  if (!gotAny)  {ObitTableSNRowUnref(row); return;}

  /* Smooth (2 sec to extrapolate) */
  smotim = 2.0 / 86400.0;
  boxsmo (smotim, wrktim, work2, numtim, work1);
  if (dotwo) boxsmo (smotim, wrktim, work4, numtim, work2);

  /* Set up for interpolation */
  ipnt1 = 0;
  ipnt2 = 1;
  time1 = wrktim[0];
  time2 = wrktim[1];
  if (numtim == 1) {  /* Only one entry in array to interpolate? */
    ipnt2 = 0;
    time2 = time1;
  } 

  /* Loop thru table changing any data with ref=ant to ref=refa */
  numrec = SNTab->myDesc->nrow;
  for (loop=1; loop<=numrec; loop++) { /* loop 200 */
    
    retCode = ObitTableSNReadRow (SNTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    if (row->status<0) continue;  /* Skip deselected record */
    
    /* right antenna? */
    if ((row->RefAnt1[iif] != ant)  &&  
	(dotwo  &&  (row->RefAnt2[iif] != ant))) continue;  /*goto L200;*/
    
    /* Right subarray? */
    sa = row->SubA;
    
    /* Interpolate */
    if ((sa == isuba)  ||  (sa <= 0)) {
      time = row->Time - timoff;
      done = FALSE;
      while (!done) {  /* loop until interpolation setup */
	if ((time >= time1)  &&  (time < time2)) {  /* L140: */
	  /* Between entries */
	  if (time2 != time1) {
	    wt1 = 1.0 - ((time-time1) / (time2-time1));
	  } else {
	    wt1 = 1.0;
	  } 
	  done = TRUE; 
	  
	} else if (time < time1) {/* Before first time */
	  wt1 = 1.0;
	  done = TRUE; 
	  
	} else if ((ipnt2+1) >= numtim) { /* After last time */
	  wt1 = 0.0;
	  done = TRUE; 
	  
	} else {   /* Shift in interpolation arrays */
	  ipnt1++;
	  ipnt2 = ipnt1+1;
	  time1 = wrktim[ipnt1];
	  time2 = wrktim[ipnt2];
	  done = FALSE;    /* goto L140; */
	}
      } /* end of loop until interpolation setup */
      
      /* Interpolate */
      wt2 = 1.0 - wt1;
      if (row->RefAnt1[iif] == ant) {
	/* Interpolate  pol 1 */
	if (row->Rate1[iif] != fblank) {
	  de1 = wt1 * work1[ipnt1] + wt2 * work1[ipnt2];
	  /* Correct  pol 1 */
	  row->Rate1[iif] -= de1;
	}  /* end data valid */
      }
      
      /* Second polarization */
      if (dotwo  &&  (row->RefAnt2[iif] == ant)) {
	if (row->Rate2[iif] != fblank) {
	  de2 = wt1 * work2[ipnt1] + wt2 * work2[ipnt2];

	  /* Correct pol 2 */
	  row->Rate2[iif] -= de2;
	} /* end data valid */
      }
      
      /* Rewrite record */
      retCode = ObitTableSNWriteRow (SNTab, loop, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    } /* end of if correct subarray */
  } /* end loop  L200: */;

  row = ObitTableSNRowUnref(row); /* delete row object */
} /* end of routine refRate */ 

/**
 * Does a box car (running mean) smoothing of irregularly spaced points.  
 * Routine translated from the AIPSish BOXSMO.FOR/BOXSMO
 * \param width  Width of boxcar in same units as x, if <=0 just copy
 * \param x      Abscissas of points to be smoothed. 
 * \param y      Values to be smoothed. 
 * \param n      Number of points to smooth. 
 * \param ys     Smoothed values. 
 */
void boxsmo (ofloat width, ofloat* x, ofloat* y, olong n, ofloat* ys) 
{
  olong   first, count, i, current, limit;
  ofloat xlast, xfirst, sum, hwidth, avg;

  if (n <= 0) return;

  /* Actual smoothing? */
  if (width <= 0.0) { /* No */
    for (i=0; i<n; i++) ys[i] = y[i];
    return;
  }

  /* Initialize. */
  first   = 0;
  current = 0;
  xlast = x[0] + width;
  hwidth = width * 0.5;
  sum   = 0.0;
  count = 0;

  /* Average window at start. */
  i = 0;
  while ((i<=n) && (x[i] <= xlast)) {
    sum += y[i++];
    count++;
  }

  /* Average in Window. */
  avg = sum / MAX (1, count);

  /* Fill in average value within half width of start */
  xlast = x[0] + hwidth;
  while (x[current] <= xlast) {
    if (current >=  n) return;  /* done? */
    ys[current] = avg;
    current++;  /* next */
  }

  /* Begin main loop. */
  xlast = x[current] + hwidth;
  while (xlast  <  x[n-1]) {
    if (current >=  n) return;  /* done? */
    xfirst = x[current] - hwidth;
    xlast  = x[current] + hwidth;
    
    /* Check if running into window at the end. */
    if (xlast  >=  x[n-1]) break;  /* goto L200;*/

    /* Init sum. */
    sum   = 0.0;
    count = 0;
    limit = first;
    for (i=limit; i<n; i++) { /* loop 150 */
      /* Check if point too early. */
      if (x[i] < xfirst) { /*goto L120;*/
	/* Too early, reset FIRST. */
	first = i + 1;

	/* Check if out of window? */
      } else if (x[i] > xlast) break; /*goto L160; */

      else { /* Do sum */
	sum += y[i];
	count++;
      }
    } /* end loop  L150: */;

    /* Average window. */
    ys[current] = sum / MAX (1, count);  /* L160: */
    current++;  /* next */
  } /* End main loop goto L100; */

  /* Sum over end window. */
  count = 0;
  sum = 0.0;	
  xfirst = x[n-1]- width;
  
  for (i= first; i<n; i++) { /* loop 210 */
    if (x[i] >= xfirst) { /*goto L210;*/
      sum += y[i];
      count++;
    }
  } /* end loop  L210: */
  
  avg = sum / MAX (1, count);  /* average in last window */
  
  /* write average to last points */
  for (i=current; i<n; i++) { /* loop 250 */
    ys[i] = avg;
  } /* end loop  L250: */;
} /* end of routine boxsmo */ 

/**
 * Average rates in an open SN table
 * \param SNTab  Table to modify, must be open
 * \param err    ObitError stack.
 */
static void avgRate (ObitTableSN *SNTab, ObitErr* err)
{
  ObitTableSNRow *row=NULL;
  ofloat sum, avg, fblank = ObitMagicF();
  olong loopr;
  olong iif, count;
  gchar * routine = "avgRate";

  /* error checks */
  if (err->error) return;
  g_assert (ObitTableSNIsA(SNTab));

  /* Create Row */
  row = newObitTableSNRow (SNTab);
  /* Attach row to output buffer */
  ObitTableSNSetRow (SNTab, row, err);

  /* Loop over table */
  for (loopr= 1; loopr<=SNTab->myDesc->nrow; loopr++) { /* loop 100 */
    ObitTableSNReadRow (SNTab, loopr, row, err);
    if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
    if (row->status<0) continue;  /* Skip deselected records */
 
    /* Average all rates */
    sum   = 0.0;
    count = 0.0;
    for (iif=0; iif<SNTab->numIF; iif++) { /* loop 20 */
      if ((row->Weight1[iif] > 0.0) && (row->Rate1[iif] != fblank)) {
	sum   += row->Rate1[iif];
	count += 1.0;
      }
      if ((SNTab->numPol > 1)  &&  (row->Weight2[iif] > 0.0) && (row->Rate2[iif] != fblank)) {
	sum   += row->Rate2[iif];
	count += 1.0;
      } 
    } /* end loop  L20:  */

    /* Average */
    if (count>0) avg = sum/count;
    else avg = fblank;

    /* Replace all rates */
    for (iif=0; iif<SNTab->numIF; iif++) { /* loop 20 */
      row->Rate1[iif] = avg;
      if (SNTab->numPol > 1) row->Rate2[iif] = avg;
    } /* end loop  L20:  */

    /* Rewrite record */
    ObitTableSNWriteRow (SNTab, loopr, row, err);
    if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
  } /* end loop over table */
  
    
} /* end avgRate */

/**
 * Return the average (over IF) Fringe rate for a given antenna and poln.
 * return 0.0 if no valid data
 * Adopted from AIPS GETMBD.FOR
 * \param in   Solution Object.
 * \param ant  antenna number (1-rel)
 * \param ipol polarization 1 or 2
 * \return average rate in s/day
 */
static ofloat avgRatePrior (ObitUVSoln *in, olong ant, olong ipol)
{
  ofloat sum, out = 0.0;
  ofloat fblank = ObitMagicF();
  olong indx, iif, count;

  indx = (ant-1) * (in->numPol * in->numIF * in->lenCalArrayEntry) +
    (ipol-1) * in->lenCalArrayEntry * in->numIF;
  count = 0;
  sum = 0.0;
  for (iif=0; iif<in->numIF; iif++) {
    if (in->CalPrior[indx+3]!=fblank) {count++; sum += in->CalPrior[indx+3];}
    indx += in->lenCalArrayEntry;
  }
  if (count>0) out = (sum/count) * 86400.0;
  return out;
} /* end avgRatePrior */

/**
 * Return the average (over IF) Fringe rate for a given antenna and poln
 * for the following solution.
 * return 0.0 if no valid data
 * Adopted from AIPS GETMBD.FOR
 * \param in   Solution Object.
 * \param ant  antenna number (1-rel)
 * \param ipol polarization 1 or 2
 * \return average rate in s/day
 */
static ofloat avgRateFollow (ObitUVSoln *in, olong ant, olong ipol)
{
  ofloat sum, out = 0.0;
  ofloat fblank = ObitMagicF();
  olong indx, iif, count;

  indx = (ant-1) * (in->numPol * in->numIF * in->lenCalArrayEntry) +
    (ipol-1) * in->lenCalArrayEntry * in->numIF;
  count = 0;
  sum = 0.0;
  for (iif=0; iif<in->numIF; iif++) {
    if (in->CalFollow[indx+3]!=fblank) {count++; sum += in->CalFollow[indx+3];}
    indx += in->lenCalArrayEntry;
  }
  if (count>0) out = (sum/count) * 86400.0;
  return out;
} /* end avgRateFollow */

/**
 * Create selector for calibrators.
 * Only values nwwded for selecting calibrators are included.
 * \param in   UVdata.
 * \param info Pointer to InfoList with selection parameters.
 * \param err  ObitErr for reporting errors.
 * \return calibrator selector 
 */
static ObitUVSel* MakeCalSelect (ObitUV *in, ObitInfoList *info, ObitErr *err)
{
  ObitUVSel *sel=NULL;
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM];
  olong itemp, Qual;
  olong iver, j, count=0;
  ObitTableSU *SUTable=NULL;
  union ObitInfoListEquiv InfoReal; 
  gchar calCode[5], *sptr;
  gchar *routine = "MakeCalSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return sel;
  g_assert (ObitUVIsA(in));
  g_assert (ObitInfoListIsA(info));

  /* create selector */
  sel = newObitUVSel ("Calibrator selector");

  InfoReal.itg = 0;type = OBIT_oint;
  ObitInfoListGetTest(info, "FreqID", &type, dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  sel->FreqID = itemp;

  /* Selected calibrators */
  Qual = -1;   /* Qualifier - default = all */
  ObitInfoListGetTest(info, "Qual", &type, dim, &Qual);
  /* Cal code */
  calCode[0] =  calCode[1] =  calCode[2] =  calCode[3] =  ' '; calCode[4] = 0; 
  ObitInfoListGetTest(info, "calCode", &type, dim, calCode);
  if (ObitInfoListGetP(info, "calSour", &type, dim, (gpointer)&sptr)) {
    sel->numberSourcesList = count;
    /* Count actual entries in source list */
    count = 0;  j = 0;
    for (i=0; i<dim[1]; i++) {
      if ((sptr[j]!=' ') || (sptr[j+1]!=' ')) count++;
      j += dim[0];
    }
    sel->numberSourcesList = count;
    /* have to lookup sources - need SU table for this. */
    iver = 1;
    SUTable = newObitTableSUValue (in->name, (ObitData*)in, &iver, OBIT_IO_ReadOnly, 0, err);
    if (SUTable==NULL) {  /* No source table - only one source and it is selected */
      sel->numberSourcesList = 0;
    } else { /* Lookup sources to get numbers */
      ObitTableSUOpen (SUTable, OBIT_IO_ReadOnly, err);
      ObitTableSUClose (SUTable, err);
      /* In case all selected */
      sel->numberSourcesList = SUTable->myDesc->nrow; 
      if (sel->sources)
	sel->sources = g_realloc(sel->sources, sel->numberSourcesList*sizeof(olong));
      else
	sel->sources = g_malloc0(sel->numberSourcesList*sizeof(olong));

      /* Do lookup */
      dim[1] = count;
      ObitTableSULookup (SUTable, dim, sptr, Qual, calCode, 
			 sel->sources, &sel->selectSources,   
			 &sel->numberSourcesList, err); 
      if(err->error)  Obit_traceback_val (err, routine, in->name, sel);
      SUTable = ObitTableSUUnref(SUTable); /* release table */
    }
    if (err->error) Obit_traceback_val (err, routine, in->name, sel);
  } else { /* no "Sources" specified */
    sel->numberSourcesList = 0; /* everything selected */
  }

  return sel;
} /* end MakeCalSelect */

/**
 * Routine to smooth amplitudes and/or phases rates in an open table.  
 * All poln present and the range of IF specified by IFBEG and IFEND  
 * are smoothed jointly.  The values in a single polarization are  
 * averaged after correcting for multiband delay.  
 * Any delay and rate smoothing should be done before amplitude and  
 * phase smoothing.  Any blanked delay and rate values will be set  
 * to 0.0.  
 * The phases are corrected by the integral of the rate functions  
 * from the first time before smoothing.   All selected phases in each  
 * polarization are averaged and corrected using the integrated phase  
 * function for the first IF selected.  
 * Routine translated from the AIPSish 
 * /export/users/bcotton/Software.dir//AIPS/31DEC02/APL/PGM/NOTST/SNSMO.FOR/SMOAPH  
 * Input table must be in antenna-time order.  
 * \param SNTab    SN table object; must be opened/closed externally
 * \param sel      pointer to uvdata selector to select data
 * \param smoFunc  Smoothing function: 'MWF', 'GAUS', else BOX 
 * \param smoType  Type of data to smooth
 *                 'AMPL', 'PHAS', 'BOTH'[def]
 * \param alpha    Alpha smooth for MWF (0 -> box, 1 -> pure MWF) 
 * \param smoParm  Smoothing time in days for: ampl, phase
 *                 0=>fill in for blanked only. 
 * \param sub      Desired subarray 
 * \param ifbeg    First IF 
 * \param ifend    Highest IF 
 * \param freqs    IF frequency array 
 * \param doBlank  replace blanked values with interpolated values?
 * \param gncnt    [out] count for gain normalization 
 * \param gnsum    [out] sum of gain modulii 
 * \param err      Error/message stack, returns if error.
 */
void smoAmpPh (ObitTableSN *SNTab, ObitUVSel *sel, gchar* smoFunc, gchar* smoType, 
	       ofloat alpha, ofloat *smoParm,
	       olong sub, olong ifbeg, olong ifend, odouble* freqs, gboolean doBlank, 
	       ofloat* gncnt, ofloat* gnsum, ObitErr* err) 
{
  ObitTableSNRow *row=NULL;
  olong   loopa, numtim, ant, nleft, i, itime, nugood, numant;
  ofloat      sumre1, sumim1, sumre2, sumim2, count1, count2, iphase, cph, rate, 
    ipre, ipim, amp, phase, mbphas, mbre, mbim, sumam1, sumam2, sumwt1, sumwt2;
  ofloat fblank = ObitMagicF();
  olong loopr, numrec, fstrec, save, isnrno=0, jtemp;
  gboolean  bad, bad2, want, need2, smophs, dotwo;
  gboolean doamp, doph, done;
  ofloat stamp, stph, lgncnt, lgnsum, arg;
  ofloat *wrkTime=NULL, *work1=NULL, *work2=NULL, *work3=NULL, *work4=NULL, *work5=NULL, 
    *work6=NULL, *work7=NULL, *work8=NULL, *work9=NULL, *work10=NULL, *work11=NULL, 
    *work12=NULL, *work13=NULL;
  olong *wrkRec=NULL;
  odouble timoff=0.0, intfaz, phadd, lstime, twopi;
  gchar *routine = "smoAmpPh";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSNIsA(SNTab));

  /* What needed */
  stamp = smoParm[0];
  stph  = smoParm[1];
  doamp = (doBlank || (stamp>0.0)) && (strncmp (smoType, "PHAS",4));
  doph  = (doBlank || (stph>0.0))  && (strncmp (smoType, "AMPL",4));
  dotwo = (SNTab->numPol > 1);  /* two poln? */
  
  twopi = 2.0 * G_PI;
  (*gncnt) = 0.0;  /* Gain constraint sims */
  (*gnsum) = 0.0;
  lgncnt   = 0.0;
  lgnsum   = 0.0;
 
  /* Get number of records in table */
  numrec = SNTab->myDesc->nrow;
  numant = SNTab->numAnt;
  if (numrec <= 0) return;

  /* Create Row */
  row = newObitTableSNRow (SNTab);
  /* Attach row to output buffer */
  ObitTableSNSetRow (SNTab, row, err);
  fstrec = 0;  /* Record number read in table */
  
  /* Create work arrays */
  wrkTime = g_malloc0(numrec*sizeof(ofloat));
  wrkRec  = g_malloc0(numrec*sizeof(olong));
  work1   = g_malloc0(numrec*sizeof(ofloat));
  work2   = g_malloc0(numrec*sizeof(ofloat));
  work3   = g_malloc0(numrec*sizeof(ofloat));
  work4   = g_malloc0(numrec*sizeof(ofloat));
  work5   = g_malloc0(numrec*sizeof(ofloat));
  work6   = g_malloc0(numrec*sizeof(ofloat));
  work7   = g_malloc0(numrec*sizeof(ofloat));
  work8   = g_malloc0(numrec*sizeof(ofloat));
  work9   = g_malloc0(numrec*sizeof(ofloat));
  work10  = g_malloc0(numrec*sizeof(ofloat));
  work11  = g_malloc0(numrec*sizeof(ofloat));
  work12  = g_malloc0(numrec*sizeof(ofloat));
  work13  = g_malloc0(numrec*sizeof(ofloat));

  /* Work usage :
     work1 
     work2 Real1
     work3 Imag1
     work4 amp1
     work5 Real2
     work6 Imag2
     work7 amp2
     work8 Integrated phase
     work9 Weight1
     work10 Weight2
     work11 work
     work12 work
     work13 work
   */

  /* Loop over antenna */
  for (loopa= 1; loopa<=numant; loopa++) { /* loop 600 */
    ant = loopa;
    /* Want this antenna? */
    if (!ObitUVSelWantAnt (sel, ant)) continue;
    
    /* Initially require both  polarizations if present */
    need2 = dotwo;
    /* Set pointers, counters */
    numtim = 0;
    nleft  = numrec - fstrec;
    nugood = 0;
    
    done = FALSE;
    while (!done) { /* Loop trying combined and separate poln L10: */
      /* Integrated phase function */
      intfaz = 0.0e0;
      lstime = 0.0;
      
      /* Loop in time, reading */
      for (loopr= 1; loopr<=nleft; loopr++) { /* loop 100 */
	isnrno = fstrec + loopr;
	ObitTableSNReadRow (SNTab, isnrno, row, err);
	if (err->error) goto cleanup;
	if (row->status<0) continue;  /* Skip deselected record */
	
	/* Finished antenna? */
	if (row->antNo < ant) continue; /* Shouldn't happen */
	if (row->antNo > ant) break;
	
	/* Check if this data wanted */
	want = (row->Time>=sel->timeRange[0]) && (row->Time<=sel->timeRange[1]);
	want = want && ObitUVSelWantSour (sel, row->SourID);
	want = want && ObitUVSelWantAnt (sel, row->antNo);
	want = want && ((row->FreqID==sel->FreqID) || (sel->FreqID<=0));
	want = want && ((row->SubA==sel->SubA) || (sel->SubA<=0));
	
	if (want) {
	  /* Average phases */
	  sumre1 = 0.0;
	  sumim1 = 0.0;
	  count1 = 0.0;
	  sumre2 = 0.0;
	  sumim2 = 0.0;
	  count2 = 0.0;
	  sumam1 = 0.0;
	  sumam2 = 0.0;
	  sumwt1 = 0.0;
	  sumwt2 = 0.0;
	  rate = fblank;
	  for (i=ifbeg; i<=ifend; i++) { /* loop 20 */
	    if (rate == fblank) rate = row->Rate1[i-1];
	    if (row->Real1[i-1] != fblank) {
	      if (rate == fblank) rate = row->Rate1[i-1];
	      /* Multiband delay for multiple IFs */
	      if ((ifend > ifbeg)  &&  (row->MBDelay1 != fblank)) {
		mbphas = -1.0 * twopi * (freqs[i-1]-freqs[0]) * row->MBDelay1;
		mbre = cos (mbphas);
		mbim = sin (mbphas);
	      } else {
		mbre = 1.0;
		mbim = 0.0;
	      } 
	      sumre1 += row->Real1[i-1]*mbre - row->Imag1[i-1]*mbim;
	      sumim1 += row->Imag1[i-1]*mbre + row->Real1[i-1]*mbim;
	      /* scalar average of amplitude */
	      sumam1 += sqrt (row->Real1[i-1]*row->Real1[i-1] +
			      row->Imag1[i-1]*row->Imag1[i-1]);
	      sumwt1 += row->Weight1[i-1];
	      count1 += 1.0;
	    } 
	    if (dotwo  &&  (rate == fblank)) rate = row->Rate1[i-1];
	    if (dotwo  &&  (row->Real2[i-1] != fblank)) {
	      if ((ifend > ifbeg)  &&  (row->MBDelay2 != fblank)) {
		mbphas = twopi * (freqs[i-1]-freqs[0]) * row->MBDelay2;
		mbre = cos (mbphas);
		mbim = sin (mbphas);
	      } else {
		mbre = 1.0;
		mbim = 0.0;
	      } 
	      sumre2 += row->Real2[i-1]*mbre - row->Imag2[i-1]*mbim;
	      sumim2 += row->Imag2[i-1]*mbre + row->Real2[i-1]*mbim;
	      
	      /* scalar average of amplitude */
	      sumam2 += sqrt (row->Real2[i-1]*row->Real2[i-1] + 
			      row->Imag2[i-1]*row->Imag2[i-1]);
	      sumwt2 += row->Weight2[i-1];
	      count2 += 1.0;
	    } 
	  } /* end loop  L20:  */;
	  /* See if flagged value */
	  bad  = count1  <=  0.1;
	  bad2 = count2  <=  0.1;
	  if (need2  &&  (bad || bad2)) {
	    bad = TRUE;
	    bad2 = TRUE;
	  } 
	  if (numtim < numrec) {
	    numtim = numtim + 1;
	    if (numtim == 1) {
	      timoff = row->Time;
	      lstime = 0.0e0;
	    } 
	    wrkTime[numtim-1] = row->Time - timoff;
	    wrkRec[numtim-1] = isnrno;
	    
	    /* Compute integrated phase function; use current rate since last time */
	    if (rate != fblank) {
	      phadd = twopi * rate * (wrkTime[numtim-1]-lstime) * 86400.0e0 * freqs[ifbeg-1];
	      intfaz = intfaz + phadd;
	      lstime = wrkTime[numtim-1];
	    } 
	    ipre = cos (intfaz);
	    ipim = sin (intfaz);
	    cph = atan2 (ipim, ipre);
	    work8[numtim-1] = cph;
	    /* Accumulate by real, imaginary, amplitude and weight. */
	    if (bad) {
	      work2[numtim-1] = fblank;
	      work3[numtim-1] = fblank;
	      work4[numtim-1] = fblank;
	      work9[numtim-1] = fblank;
	    } else {
	      nugood = nugood + 1;
	      sumre1 = sumre1 / count1;
	      sumim1 = sumim1 / count1;
	      sumam1 = sumam1 / count1;
	      
	      /* Subtract integrated phase */
	      work2[numtim-1] = sumre1*ipre + sumim1*ipim;
	      work3[numtim-1] = sumim1*ipre - sumre1*ipim;
	      work4[numtim-1] = sqrt (work2[numtim-1]*work2[numtim-1] + 
				      work3[numtim-1]*work3[numtim-1]) + 1.0e-20;
	      work9[numtim-1] = sumwt1;
	      
	      /* Normalize real and imag. */
	      work2[numtim-1] = work2[numtim-1] / work4[numtim-1];
	      work3[numtim-1] = work3[numtim-1] / work4[numtim-1];
	      if (ifend > ifbeg) work4[numtim-1] = sumam1;
	    } 
	    if (bad2) {
	      work5[numtim-1] = fblank;
	      work6[numtim-1] = fblank;
	      work7[numtim-1] = fblank;
	      work10[numtim-1]= fblank;
	    } else {
	      nugood = nugood + 1;
	      sumre2 = sumre2 / count2;
	      sumim2 = sumim2 / count2;
	      sumam2 = sumam2 / count2;
	      
	      /* Subtract integrated phase */
	      work5[numtim-1] = sumre2*ipre + sumim2*ipim;
	      work6[numtim-1] = sumim2*ipre - sumre2*ipim;
	      work7[numtim-1] = sqrt (work5[numtim-1]*work5[numtim-1] + 
				      work6[numtim-1]*work6[numtim-1]) + 1.0e-20;
	      work10[numtim-1]= sumwt1;
	      
	      /* Normalize real and imag. */
	      work5[numtim-1] = work5[numtim-1] / work7[numtim-1];
	      work6[numtim-1] = work6[numtim-1] / work7[numtim-1];
	      if (ifend > ifbeg) work7[numtim-1] = sumam2;
	    } 
	  } 
	}
      } /* end loop  L100: */
      
      save = isnrno - 1;  /* How far did we get? */
      /* Anything with both poln? */
      if (need2  &&  (nugood <= 0)) {
	/* Try again with only one poln. */
	done  = FALSE;  /* Need another pass to get both poln */
	need2 = FALSE;
      } else done = TRUE;
    } /* End of loop over poln L10: */

    /* Switch back to phase in work2 and work5: */
    for (itime= 1; itime<= numtim; itime++) { /* loop 120 */
      if (work2[itime-1] != fblank) {
	/* No need to do anything if the real part is fblank --- the phase will be too. */
	work2[itime-1] = atan2 (work3[itime-1], work2[itime-1] + 1.0e-20);
	if (itime > 1) {
	  /* Nearest integer */
	  arg = (work2[itime-1] - work2[itime-2]) / twopi;
	  if (arg>=0) jtemp = (olong)(arg + 0.5);
	  else jtemp = (olong)(arg - 0.5);
	  if (work2[itime-2] != fblank) work2[itime-1] = work2[itime-1] - twopi * jtemp;
	} 
      } 
      if (dotwo && (work5[itime-1] != fblank)) {
	work5[itime-1] = atan2 (work6[itime-1], work5[itime-1] + 1.0e-20);
	if (itime > 1) {
	  /* Nearest integer */
	  arg = (work5[itime-1] - work5[itime-2]) / twopi;
	  if (arg>=0) jtemp = (olong)(arg + 0.5);
	  else jtemp = (olong)(arg - 0.5);
	  if (work5[itime-2] != fblank) work5[itime-1] = work5[itime-1] - twopi * jtemp;
	} 
      } 
     } /* end loop   L120 */;

    /* Smooth as requested (smoothed as total phase and amplitude ) */
    smoIt (smoFunc, stph,  alpha, wrkTime, work2, work9, numtim, 
	   work1, work11, work12, work13, doBlank);
    smoIt (smoFunc, stamp, alpha, wrkTime, work4, work9, numtim, 
	   work2, work11, work12, work13, doBlank);
    
    /* Second polarization if present */
    if (dotwo) {
      smoIt (smoFunc, stph,  alpha, wrkTime, work5, work10, numtim, 
	     work3, work11, work12, work13, doBlank);
      smoIt (smoFunc, stamp, alpha, wrkTime, work7, work10, numtim, 
	     work4, work11, work12, work13, doBlank);
    } 

    /* Replace with smoothed values */
    for (itime = 1; itime <= numtim; itime++) { /* loop 200 */
      isnrno = wrkRec[itime-1];
      ObitTableSNReadRow (SNTab, isnrno, row, err);
      if (err->error) goto cleanup;
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Update */
      if ((work1[itime-1] != fblank)  &&  (work2[itime-1] != fblank)) {
	/* Smoothed phase same for all IFs;  add integrated phase function. */
	if (work1[itime-1] != fblank) {
	  iphase = work1[itime-1] + work8[itime-1];
	} else {
	  iphase = fblank;
	} 
	for (i=ifbeg; i<=ifend; i++) { /* loop 130 */
	  /* Set AMP and PHASE by smoothing  selected. */
	  if (doamp) {
	    amp = work2[itime-1];
	  } else {
	    /* Use smoothed value if blanked. */
	    if (row->Real1[i-1] == fblank) {
	      amp = work2[itime-1];
	    } else {
	      amp = sqrt (row->Real1[i-1]*row->Real1[i-1] + row->Imag1[i-1]*row->Imag1[i-1]);
	    } 
	  } 
	  /* Phase */
	  if (doph) {
	    phase = iphase;
	    smophs = iphase != fblank;
	  } else {
	    /* Use smoothed value if blanked. */
	    if (row->Real1[i-1] == fblank) {
	      phase = iphase;
	      smophs = iphase != fblank;
	    } else {
	      phase = atan2 (row->Imag1[i-1], row->Real1[i-1]+1.0e-20);
	      smophs = FALSE;
	    } 
	  } 
	  /* Save smoothed values */
	  if ((iphase != fblank)  &&  (amp != fblank)) {
	    /* Multiband delay correction */
	    if ((ifend > ifbeg)  &&  smophs  &&  (row->MBDelay1 != fblank)) {
	      phase = phase + twopi * (freqs[i-1]-freqs[0]) * row->MBDelay1;
	    } 
	    
	    /* Check if flagged entries are  to be overwritten */
	    if (row->Real1[i-1] != fblank) {
	      row->Real1[i-1] = amp * cos (phase);
	      row->Imag1[i-1] = amp * sin (phase);
	      /* Since flagged entries can be "recovered" here set nominal  weight. */
	      row->Weight1[i-1] = MAX (1.0, row->Weight1[i-1]);
	      
	      /* Unflag delay and rate if  necessary. */
	      if (row->Delay1[i-1] == fblank) row->Delay1[i-1] = 0.0;
	      if (row->Rate1[i-1]  == fblank) row->Rate1[i-1]  = 0.0;
	      /* Keep track of mean gain modulus */
	      lgncnt += 1.0;
	      lgnsum += amp;
	    } 
	  } else {
	    row->Real1[i-1]  = fblank;
	    row->Imag1[i-1]  = fblank;
	    row->Weight1[i-1] = 0.0;
	  } 
	} /* end loop  L130: */;
      } 

      /* Second polarization present? */
      if (dotwo  &&  (work3[itime-1] != fblank)  &&  (work4[itime-1] != fblank)) {
	/* Smoothed phase same for all IFs. */
	if (work3[itime-1] != fblank) {
	  iphase = work3[itime-1] + work8[itime-1];
	} else {
	  iphase = fblank;
	} 
	for (i=ifbeg; i<=ifend; i++) { /* loop 140 */
	  /* Set AMP and PHASE by smoothing selected. */
	  if (doamp) {
	    amp = work4[itime-1];
	  } else {
	    /* Use smoothed value if blanked. */
	    if (row->Real2[i-1] == fblank) {
	      amp = work4[itime-1];
	    } else {
	      amp = sqrt (row->Real2[i-1]*row->Real2[i-1] + row->Imag2[i-1]*row->Imag2[i-1]);
	    } 
	  } 
	  /* Phase */
	  if (doph) {
	    phase = iphase;
	    smophs = iphase != fblank;
	  } else {
	    /* Use smoothed value if blanked. */
	    if (row->Real2[i-1] == fblank) {
	      phase = iphase;
	      smophs = iphase != fblank;
	    } else {
	      phase = atan2 (row->Imag2[i-1], row->Real2[i-1]+1.0e-20);
	      smophs = FALSE;
	    } 
	  } 
	  
	  /* Save corrected data */
	  if ((iphase != fblank)  &&  (amp != fblank)) {
	    /* Multiband delay correction */
	    if ((ifend > ifbeg)  &&  smophs  &&  (row->MBDelay2 != fblank)) {
	      phase = phase + twopi * (freqs[i-1]-freqs[0]) * row->MBDelay2;
	    } 
	    
	    /* Check if flagged entries are to be overwritten */
	    if (row->Real2[i-1] != fblank) {
	      row->Real2[i-1] = amp * cos (phase);
	      row->Imag2[i-1] = amp * sin (phase);
	      row->Weight2[i-1] = MAX (1.0, row->Weight2[i-1]);
	      if (row->Delay2[i-1] == fblank) row->Delay2[i-1] = 0.0;
	      if (row->Rate2[i-1]  == fblank) row->Rate2[i-1]  = 0.0;
	      
	      /* Keep track of mean gain modulus */
	      lgncnt += 1.0;
	      lgnsum += amp;
	    } 
	  } else {
	    row->Real2[i-1] = fblank;
	    row->Imag2[i-1] = fblank;
	    row->Weight2[i-1] = 0.0;
	  } 
	} /* end loop  L140: */
      } 
      
      /* Rewrite record */
      ObitTableSNWriteRow (SNTab, isnrno, row, err);
      if (err->error) goto cleanup;
    } /* end loop  L200: */
    
    /* First SN number of next antenna */
    fstrec = save;
  } /* end antenna loop  L600: */

  
  /* cleanup */
 cleanup:
  row = ObitTableSNRowUnref(row);
  if (wrkTime) g_free(wrkTime);
  if (wrkRec)  g_free(wrkRec);
  if (work1)   g_free(work1);
  if (work2)   g_free(work2);
  if (work3)   g_free(work3);
  if (work4)   g_free(work4);
  if (work5)   g_free(work5);
  if (work6)   g_free(work6);
  if (work7)   g_free(work7);
  if (work8)   g_free(work8);
  if (work9)   g_free(work9);
  if (work10)  g_free(work10);
  if (work11)  g_free(work11);
  if (work12)  g_free(work12);
  if (work13)  g_free(work13);
  if (err->error) Obit_traceback_msg (err, routine, SNTab->name);
  
  /* Set output */
  (*gncnt) = lgncnt;  /* Gain constraint sims */
  (*gnsum) = lgnsum;
  
} /* end of routine smoAmpPh */ 

/**
  * Routine to call appropriate smoothing routine.  Magic value blanking  
 * is supported.  
 * Routine adopted from the AIPSish 
 * 31DEC02/APL/PGM/NOTST/SNSMO.FOR/SNSMSM  
 * \param smmeth  Method 'BOX','MWF', 'GAUS', unknown = 'BOX' 
 * \param width   Smoothing time (days) 
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \param x       Abscissas of points to be smoothed in increasing 
 *                order 
 * \param y       Values to be smoothed. 
 * \param w       Weights of data. 
 * \param n       Number of points to smooth. 
 * \param ys      [out] Smoothed values. 
 * \param ws      [out] Smoothed weights 
 * \param yor     Scratch 
 * \param wor     Scratch 
 * \param doBlank replace blanked values with interpolated values.
 */
void 
smoIt (gchar* smmeth, ofloat width, ofloat alpha, 
       ofloat* x, ofloat *y, ofloat *w, olong n, 
       ofloat* ys, ofloat* ws, ofloat *wrk1, ofloat *wrk2, gboolean doBlank) 
{
  /* Any work to do? */
  if (n <= 0) return;

  /* Smooth */
  if (!strncmp (smmeth, "BOX",3)) {
    ObitUVSolnSmooBox (width, x, y, w, n, ys, ws, doBlank);
  } else if (!strncmp (smmeth, "MWF",3)) {
    ObitUVSolnSmooMWF (width, alpha, x, y, w, n, ys, ws, wrk1, wrk2, doBlank);
  } else if (!strncmp (smmeth, "Gaus",4)) {
    ObitUVSolnSmooGauss (width, x, y, w, n, ys, ws, wrk1, doBlank);
  } else { /* Default "BOX" */
    ObitUVSolnSmooBox (width, x, y, w, n, ys, ws, doBlank);
  }
} /* end of routine smoIt */ 
