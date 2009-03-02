/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include "ObitUVCalCalibrate.h"
#include "ObitUVCalCalibrateDef.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitTableAN.h"
#include "ObitTableCL.h"
#include "ObitTableCQ.h"
#include "ObitTableSN.h"
#include "ObitTableUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalCalibrate.c
 * ObitUVCal utilities for applying amp/phase/delay/rate calibration to
 * uv data.
 */

/*-----------------------Macroes-------------------------*/
/** Velocity of light */
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for amp/phase/delay/rate Calibrate . */
static ObitUVCalCalibrateS* newObitUVCalCalibrateS (ObitUVCal *in);

/** Private: Update calibration arrays. */
static void ObitUVCalCalibrateUpdate (ObitUVCalCalibrateS *in, ofloat time,
					ObitErr *err);

/** Private:  Read calibration for a new time into the internal arrays. */
static void ObitUVCalCalibrateNewTime (ObitUVCalCalibrateS *in, ofloat time,
					ObitErr *err);

/** Private: Initialize VLBA corrections from CQ table */
static void
ObitUVCalCalibrateVLBAInit (ObitUVCal *in, ObitUVCalCalibrateS *out, 
			      ObitErr *err);

/** Private: calculate the VLBA segmentation loss factor */
static double 
ObitUVCalCalibrateSegLoss (olong lFunc, olong nfft, odouble dbits, olong *rc);

/** Private: calculate rate decorrelation loss factor */
static float 
ObitUVCalCalibrateRateDecorr (olong filter, ofloat rate, ofloat TimeI, ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for amp/phase/delay/rate Calibrate .
 * \param in   Calibrate Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitUVCalCalibrateInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
		    ObitErr *err)
{
  ObitIOCode retCode;
  ObitUVCalCalibrateS *me;
  olong size, i;
  gchar *colName="TIME    ";
  gchar *routine="ObitUVCalCalibrateInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  in->doCal = sel->doCal;
  if (!in->doCal) return;

  in->ampPhaseCal = newObitUVCalCalibrateS(in);

  /* pointer to calibration structure */
  me = in->ampPhaseCal;

  /* Copy Selector information */
  me->doCalWt     = sel->doCalWt;
  me->bChan       = sel->startChann;
  me->eChan       = sel->startChann + sel->numberChann - 1;
  me->bIF         = sel->startIF;
  me->eIF         = sel->startIF + sel->numberIF - 1;
  me->SubA        = sel->SubA;
  me->FreqID      = sel->FreqID;

  /* Copy Cal information */
  me->CLTable     = ObitTableRef(in->CLTable);
  me->SNTable     = ObitTableRef(in->SNTable);
  me->doSNTable   = (in->SNTable!=NULL);
  me->LastRowRead = 0;
  me->numSubA     = in->numSubA;
  me->numIF       = desc->inaxes[desc->jlocif];
  me->numChan     = desc->inaxes[desc->jlocf];

  /* Copy descriptor information */
  me->numAnt    = desc->maxAnt;
  me->numSubA   = desc->numSubA;
  me->DeltaTime = MAX (desc->DeltaTime, sel->InputAvgTime);
  me->numIF     = desc->inaxes[desc->jlocif];
  me->numChan   = desc->inaxes[desc->jlocf];

  /* Sort/Open calibration table, create row structure, get numPol  */
  if (me->doSNTable) { /* SN */
    /* Sort to time order if needed */
    retCode = ObitTableUtilSort((ObitTable*)(me->SNTable), colName, FALSE, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    retCode = 
      ObitTableSNOpen ((ObitTableSN*)(me->SNTable), OBIT_IO_ReadOnly, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    me->SNTableRow = (Obit*)newObitTableSNRow((ObitTableSN*)(me->SNTable));
    me->numPol = ((ObitTableSN*)me->SNTable)->numPol;
    me->numRow = ((ObitTableSN*)me->SNTable)->myDesc->nrow;
  } else {  /* CL */
    /* Sort to time order if needed */
    retCode = ObitTableUtilSort((ObitTable*)(me->CLTable), colName, FALSE, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    retCode = 
      ObitTableCLOpen ((ObitTableCL*)(me->CLTable), OBIT_IO_ReadOnly, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    me->CLTableRow = (Obit*)newObitTableCLRow((ObitTableCL*)(me->CLTable));
    me->numPol = ((ObitTableCL*)me->CLTable)->numPol;
    me->numRow = ((ObitTableCL*)me->CLTable)->myDesc->nrow;
  }

  /* Allocate calibration arrays */
  me->lenCalArrayEntry = 5; /* length of cal array entry */
  size = me->numAnt * (me->eIF- me->bIF + 1) * me->numPol * me->lenCalArrayEntry;
  me->CalApply     = g_malloc(size*sizeof(float));
  me->CalPrior     = g_malloc(size*sizeof(float));
  me->CalFollow    = g_malloc(size*sizeof(float));
  me->IFR          = g_malloc(me->numAnt*sizeof(float));
  me->PriorIFR     = g_malloc(me->numAnt*sizeof(float));
  me->FollowIFR    = g_malloc(me->numAnt*sizeof(float));
  me->DDelay       = g_malloc(me->numAnt*sizeof(float));
  me->PriorDDelay  = g_malloc(me->numAnt*sizeof(float));
  me->FollowDDelay = g_malloc(me->numAnt*sizeof(float));
  me->RateFact     = g_malloc(me->numIF*sizeof(float));
  me->DelayFact    = g_malloc(me->numIF*sizeof(float));
  me->Lambda       = g_malloc(me->numIF*me->numChan*sizeof(float));
  me->PriorAntTime = g_malloc(me->numAnt*sizeof(float));
  me->FollowAntTime= g_malloc(me->numAnt*sizeof(float));

  /* Fill frequency related arrays */
  for (i=0; i<me->numIF; i++) {
    me->RateFact[i] = 2.0 * G_PI * 86400.0 * desc->freqIF[i];
    me->DelayFact[i] = 2.0 * G_PI * desc->chIncIF[i];
  }

  /* Wavelengths */
  me->numLambda = me->numChan;  /* Number per IF */
  for (i=0; i<me->numIF*me->numChan; i++) {
    me->Lambda[i] = VELIGHT / in->myDesc->freqArr[i];
  }

  /* Polarization offset array PolOff  */
  for (i=0; i<4; i++) {
    me->PolOff[0][i] = 0;
    me->PolOff[1][i] = 0;
   }
  /* if making true stokes from RR,LL,... */
  if (desc->crval[desc->jlocs] < 0.0) {
    me->PolOff[0][1] = me->lenCalArrayEntry;
    me->PolOff[0][3] = me->lenCalArrayEntry;
    me->PolOff[1][1] = me->lenCalArrayEntry;
    me->PolOff[1][2] = me->lenCalArrayEntry;
  }

  /* Initial times to trigger update of calibration arrays */
  me->CalTime       = -1.0e20;
  me->PriorCalTime  = -1.0e20;
  me->FollowCalTime = -1.0e20;

  /* Initialize any VLBA calibration*/
  ObitUVCalCalibrateVLBAInit (in, me, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /*  end ObitUVCalCalibrateInit */

/**
 * Calibrate data amp/phase/delay/rate
 * A calibration array entry consists of:
 * \li real part of gain
 * \li imaginary of gain
 * \li sine of channel to channel phase rotation
 * \li cosine of channel to channel phase rotation
 * \li fringe rate (rad/day)
 * \param in    Calibrate Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalCalibrate (ObitUVCal *in, ofloat time, olong ant1, olong ant2, 
			 ofloat *RP, ofloat *visIn, ObitErr *err)
{
  olong   indxa1, indxa2, asize, iif, ipol, ifreq, ioff, joff, index, 
    jndxa1, jndxa2, maxpol, idndx, itfilt, corID, iSubA, itemp;
  gboolean   sombad, somflg, allflg, smpflg, alpflg, allded, ccor;
  gboolean calBad, doDisp, badDsp;
  ofloat tvr, tvi, tvr1, gr, gi, dgr, dgi, phase, grd, gid;
  ofloat  cp, sp, gwt, dphas, rate, arg=0.0, rfact, dfact, fblank = ObitMagicF();
  odouble dbits, dsfact;
  ObitUVCalCalibrateS *me;
  ObitUVDesc *desc;
  ObitUVSel *sel;
  olong rc;
  gchar *routine="ObitUVCalCalibrate";

  /* local pointers for structures */
  me   = in->ampPhaseCal;
  desc = in->myDesc;
  sel  = in->mySel;

  /* Correlator ID if in data */
  if (desc->ilocid>=0) 
    corID = (olong)RP[desc->ilocid] + 0.5;
  else
    corID = 1;

  /* Subarray number in data */
  itemp = (olong)RP[desc->ilocb];
  iSubA = 1 + (olong)(100.0*(RP[desc->ilocb] -(ofloat)itemp) + 0.1);

  /* Integration time if in data and not already specified */
  if (me->DeltaTime<=0.0) {
    if (desc->ilocit>=0) 
      me->DeltaTime = RP[desc->ilocid];
    else
      me->DeltaTime = 1.0 / 86400.0; /* Default 1 sec */
  }

 /* see if new time - update cal. */
  if (time > me->CalTime) {
    ObitUVCalCalibrateUpdate (me, time, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* check if cross correlation */
  ccor = ant1 != ant2;

  /* init. flagged flags */
  allflg = TRUE;
  allded = TRUE;
  alpflg = TRUE;
  somflg = FALSE;
  smpflg = FALSE;

  /* handle 1 solution for both pols */
  maxpol = MAX (1, MIN(in->numStok, me->numPol*me->numPol));

  /* doDisp should be TRUE if either  dispersive delay is not zero */
  doDisp = (me->DDelay[ant1-1] != 0.0) && (me->DDelay[ant2-1] != 0.0);
  
  /* badDsp should be TRUE if either dispersive delay is blanked */
  badDsp = (me->DDelay[ant1-1]  ==  fblank)  ||  (me->DDelay[ant2-1]  ==  fblank);
  
  /* set antenna indices */
  asize = me->numPol * (me->eIF - me->bIF + 1) * me->lenCalArrayEntry;
  jndxa1 = (ant1 - 1) * asize;
  jndxa2 = (ant2 - 1) * asize;

  /* loop over IF */
  for (iif= me->bIF; iif<=me->eIF; iif++) { /* loop 300 */
    ioff = (iif-1) * desc->incif;
    
    /* loop over polarization */
    for (ipol= 1; ipol<=in->numStok; ipol++) { /* loop 200 */

      /* check baseline flags */
      joff = ioff + (ipol-1) * desc->incs;
      /* handle 1 solution for 2 polzns. */
      indxa1 = jndxa1 + me->PolOff[0][MIN(ipol,maxpol)-1];
      indxa2 = jndxa2 + me->PolOff[1][MIN(ipol,maxpol)-1];

      /* Initialize calibration */
      gr = 1.0;
      gi = 0.0;
      dgr = 1.0;
      dgi = 0.0;
      gwt = 0.0;

      /* check IF flags */
      calBad = (me->CalApply[indxa1] == fblank) || (me->CalApply[indxa2] == fblank);

      /* set gains from amplitude/phase part of calibration a1 * conjg(a2) */
      if (!calBad) {
	gr = me->CalApply[indxa1] * me->CalApply[indxa2] + me->CalApply[indxa1+1] * me->CalApply[indxa2+1];
	gi = me->CalApply[indxa2] * me->CalApply[indxa1+1] - me->CalApply[indxa1] * me->CalApply[indxa2+1];

	/* "weight" calibration */
	if (me->doCalWt) {
	  gwt = (gr*gr + gi*gi);
	  if (gwt > 1.0e-10) gwt = 1.0 / gwt;
	} else {
	  gwt = 1.0;
	} 
      } /* end of set gains to real/imag part */
    
      /* see if delay-rate corrections  wanted. */
      if (!calBad && me->doDelayRate) {

	/* delay correction - real and imaginary of phase rotation per channel */
	dgr = me->CalApply[indxa1+2] * me->CalApply[indxa2+2] + me->CalApply[indxa1+3] * me->CalApply[indxa2+3];
	dgi = me->CalApply[indxa2+2] * me->CalApply[indxa1+3] - me->CalApply[indxa1+2] * me->CalApply[indxa2+3];

	/* apply fringe rate */
	rate  = me->CalApply[indxa1+4] - me->CalApply[indxa2+4];
	phase = rate * (time - me->CalTime);
	
	/* Apply decorrelation corrections as applicable */
	idndx = (corID - 1) * me->numIF + iif - 1;

	/* vlba-only corrections */
	if ((me->doDelayDecorr) && me->doDelayDecorr[idndx] && (me->corrType[iSubA-1] == 1)) {

	  /* spectral averaging correction - get delay phase channel-channel phase rotation */
	  arg = 0.5 * atan2 (dgi, dgr);
	  if ((abs (arg) > 1.0e-5) && (me->NSpecA[idndx] > 1)) {
	    dfact = me->NSpecA[idndx] * sin (arg / me->NSpecA[idndx]) / sin (arg);
	    dfact = fabs (dfact);
	    gr = gr * dfact;
	    gi = gi * dfact;
	  } 

	  /* Segmentation loss correction compute residual delay in bits */
	  dbits = (2.0 * arg) / (me->DelayFact[iif-1] * me->DelBit[idndx]);

	  /*  Calculate segmentation loss factor */
	  dsfact = ObitUVCalCalibrateSegLoss (me->LTaper[idndx], me->NFFTSize[idndx], 
					      dbits, &rc);

	  /* correct gain factors */
	  gr = gr / dsfact;
	  gi = gi / dsfact;
	} /* end vlba-only corrections */

	/* Compute rate smearing; correction (includes the vlba ovlb filters) */
	if (me->doRateSmear) {
	  /* default is boxcar smoothing */
	  itfilt = 0;

	  /* check if vlba filter used */
	  if ((me->doDelayDecorr) && me->doDelayDecorr[idndx] && (me->corrType[iSubA-1] == 1)) 
	    itfilt = me->typeTimeFilt[idndx-1];
	  else 
	    /* default is boxcar smoothing */
	    itfilt = 0;

	  /* Compute loss factor */
	  rfact = ObitUVCalCalibrateRateDecorr (itfilt, rate, me->DeltaTime, err);
	  if (err->error) Obit_traceback_msg (err, routine, in->name);
	  gr = gr * rfact;
	  gi = gi * rfact;
	} /* end rate smearing correction */

	/* Caution if VLBA corrections 	not possible due to missing  
	   CQ table or CQ table entries */
	if ((me->doDelayDecorr) && (me->corrType[iSubA-1] == 1) && (!me->doDelayDecorr[idndx]) &&  
	    me->warnNoCQ[iSubA-1] &&  ((fabs(arg-1.0) > 0.0) || (fabs(rate-1) > 0.0))) {
	  Obit_log_error(err, OBIT_InfoWarn, 
			 "Warning: subarray %d  contains vlba data, but no CQ data", iSubA);
	  Obit_log_error(err, OBIT_InfoWarn, 
			 "This is needed to vlba decorrelation corrections");

	  /* Only print warning once */
	   me->warnNoCQ[iSubA-1] = FALSE;
	} /* end CQ warning */

	/* correct for frequency offset. - get delay phase channel-channel phase rotation */
	if ((sel->startChann > 1)  ||  (fabs (desc->crpix[desc->jlocf] - 1.0) >  0.0001)) {
	  dphas = atan2 (dgi, dgr);
	  phase = phase + dphas * (sel->startChann-desc->crpix[desc->jlocf]);
	} 
	cp = cos (phase);
	sp = sin (phase);
	tvr = gr * cp - gi * sp;
	tvi = gr * sp + gi * cp;
	gr = tvr;
	gi = tvi;
      } /* end of delay/rate calibration */


      /* loop over channels calibrating. */
      for (ifreq= me->bChan-1; ifreq<me->eChan; ifreq++) { /* loop 80 */
	index = joff + (ifreq) * desc->incf;
	tvr = gr * visIn[index]   + gi * visIn[index+1];
	tvi = gr * visIn[index+1] - gi * visIn[index];

	/* apply dispersive delay  correction or kill data point if dispersive delay is blanked */
	if (doDisp) {
	  if (badDsp) {
	    gwt = 0.0;
	  } else {
	    phase = me->Lambda[me->numLambda * (iif-1) + ifreq] * 
	      (me->DDelay[ant1-1] - me->DDelay[ant2-1]);
	    grd = cos (phase);
	    gid = sin (phase);
	    tvr1 = tvr;
	    tvr = grd * tvr1 + gid * tvi;
	    tvi = grd * tvi - gid * tvr1;
	  } 
	} 

	/* keep track of the flagging */
	smpflg = smpflg  ||  (visIn[index+2]  <=  0.0);
	alpflg = alpflg  &&  (visIn[index+2]  <=  0.0);
	somflg = somflg  ||  (gwt  <=  0.0);
	allflg = allflg  &&  (gwt  <=  0.0);
	allded = allded  &&  ((visIn[index+2]  <=  0.0) ||  (gwt  <=  0.0));

	/* If data exactly zero - flag (fix VLA screwup) */
	if ((visIn[index]==0.0) && (visIn[index+1]==0.0)) visIn[index+2] = 0.0;

	/* save calibrated data */
	if ((visIn[index+2] > 0.0) && (!calBad)) {
	  visIn[index]   = tvr;
	  visIn[index+1] = tvi;
	  visIn[index+2] = visIn[index+2] * gwt;
	} else {
	  visIn[index]   = 0.0;
	  visIn[index+1] = 0.0;
	  visIn[index+2] = 0.0;
	}

	/* DEBUG 
	if ((fabs(visIn[index])>10000.0) || (fabs(visIn[index+1])>10000.0)) {
	  fprintf (stderr,"DEBUG Bad data vis %f %d %d %f %f\n",
		   time,ant1, ant2,visIn[index],visIn[index+1]);
	}*/

	/* rotate phase correction for next if we have delay corrections */
	if (me->doDelayRate) {
	  tvr = gr * dgr - gi * dgi;
	  tvi = gr * dgi + gi * dgr;
	  gr = tvr;
	  gi = tvi;
	} 

      } /* end loop over channels  L80 */;

    } /* end loop loop over Stokes  L200 */;

    /* setup for next IF */
    jndxa1 = jndxa1 + me->lenCalArrayEntry * me->numPol;
    jndxa2 = jndxa2 + me->lenCalArrayEntry * me->numPol;
  } /* end loop  L300 - loop over IF */;

  /* increment counts of the good, bad and the ugly. */
  somflg = somflg  &&  (!allflg);
  smpflg = smpflg  &&  (!alpflg);
  sombad = (somflg || smpflg)  &&  (!allded);
  if (smpflg) me->countRec[0][0]++;
  if (somflg) me->countRec[1][0]++;
  if (sombad) me->countRec[2][0]++;
  if (alpflg) me->countRec[0][1]++;
  if (allflg) me->countRec[1][1]++;
  if ((!allded)  &&  (!sombad))  me->countRec[2][1]++;

} /* end ObitUVCalCalibrate */


/**
 * Shutdown amp/phase/delay/rate Calibrate.
 * Close any open file and destroy structures.
 * \param in   Calibrate Object.
 * \param err  ObitError stack.
 */
void ObitUVCalCalibrateShutdown (ObitUVCal *in, ObitErr *err)
{
  ObitUVCalCalibrateS *me;
  ObitIOCode retCode;
  gchar *routine="ObitUVCalCalibrateShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* Close calibration table, release row structure  */
  me = in->ampPhaseCal;
  if (me->doSNTable) { /* SN */
    retCode = ObitTableSNClose ((ObitTableSN*)me->SNTable, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    me->SNTableRow = ObitTableSNRowUnref(me->SNTableRow);
  } else {  /* CL */
    retCode = ObitTableCLClose ((ObitTableCL*)me->CLTable, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    me->CLTableRow = ObitTableCLRowUnref(me->CLTableRow);
  }

  /* delete structure */
  in->ampPhaseCal = ObitUVCalCalibrateSUnref(in->ampPhaseCal);
} /*  end ObitUVCalCalibrateShutdown */

/**
 * Destroy structure for amp/phase/delay/rate Calibrate .
 * \param in   Calibrate Object.
 * \return NULL
 */
ObitUVCalCalibrateS*
ObitUVCalCalibrateSUnref (ObitUVCalCalibrateS *in)
{
  if (in==NULL) return in;;
  in->CLTable    = ObitTableCLUnref((ObitTableCL*)in->CLTable);
  in->CLTableRow = ObitTableCLRowUnref((ObitTableCLRow*)in->CLTableRow);
  in->SNTable    = ObitTableSNUnref((ObitTableSN*)in->SNTable);
  in->SNTableRow = ObitTableSNRowUnref((ObitTableSNRow*)in->SNTableRow);
  if (in->CalApply)     g_free(in->CalApply); in->CalApply   = NULL;
  if (in->CalPrior)     g_free(in->CalPrior); in->CalPrior   = NULL;
  if (in->CalFollow)    g_free(in->CalFollow); in->CalFollow  = NULL;
  if (in->IFR)          g_free(in->IFR); in->IFR   = NULL;
  if (in->PriorIFR)     g_free(in->PriorIFR); in->PriorIFR   = NULL;
  if (in->FollowIFR)    g_free(in->FollowIFR); in->FollowIFR  = NULL;
  if (in->DDelay)       g_free(in->DDelay); in->DDelay     = NULL;
  if (in->PriorDDelay)  g_free(in->PriorDDelay); in->PriorDDelay  = NULL;
  if (in->FollowDDelay) g_free(in->FollowDDelay); in->FollowDDelay = NULL;
  if (in->RateFact)     g_free(in->RateFact); in->RateFact  = NULL;
  if (in->DelayFact)    g_free(in->DelayFact); in->DelayFact = NULL;
  if (in->Lambda)       g_free(in->Lambda); in->Lambda    = NULL;
  if (in->PriorAntTime) g_free(in->PriorAntTime); in->PriorAntTime    = NULL;
  if (in->FollowAntTime) g_free(in->FollowAntTime); in->FollowAntTime    = NULL;
  if (in->LTaper)       g_free(in->LTaper); in->LTaper    = NULL;
  if (in->NSpecA)       g_free(in->NSpecA); in->NSpecA    = NULL;
  if (in->DelBit)       g_free(in->DelBit); in->DelBit    = NULL;
  if (in->NFFTSize)     g_free(in->NFFTSize); in->NFFTSize  = NULL;
  if (in->typeTimeFilt) g_free(in->typeTimeFilt); in->typeTimeFilt   = NULL;
  if (in->TimeFiltTime) g_free(in->TimeFiltTime); in->TimeFiltTime   = NULL;
  if (in->doDelayDecorr) g_free(in->doDelayDecorr); in->doDelayDecorr  = NULL;
  if (in->corrType)     g_free(in->corrType); in->corrType  = NULL;
  if (in->warnNoCQ)     g_free(in->warnNoCQ); in->warnNoCQ = NULL;

  /* basic structure */
   g_free (in);

  return NULL;
} /*  end ObitUVCalCalibrateSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for amp/phase/delay/rate Calibrate .
 * \param in   Calibrate Object.
 * \return newly created object.
 */
static ObitUVCalCalibrateS*
newObitUVCalCalibrateS (ObitUVCal *in)
{
  ObitUVCalCalibrateS* out;

  out = g_malloc0(sizeof(ObitUVCalCalibrateS));

  /* Null pointers */
  out->CLTable    = NULL;
  out->CLTableRow = NULL;
  out->SNTable    = NULL;
  out->SNTableRow = NULL;
  out->CalApply   = NULL;
  out->CalPrior   = NULL;
  out->CalFollow  = NULL;
  out->IFR        = NULL;
  out->PriorIFR   = NULL;
  out->FollowIFR  = NULL;
  out->DDelay     = NULL;
  out->PriorDDelay  = NULL;
  out->FollowDDelay = NULL;
  out->RateFact  = NULL;
  out->DelayFact = NULL;
  out->Lambda    = NULL;
  out->PriorAntTime = NULL;
  out->FollowAntTime= NULL;
  out->LTaper    = NULL;
  out->NSpecA    = NULL;
  out->DelBit    = NULL;
  out->NFFTSize  = NULL;
  out->typeTimeFilt   = NULL;
  out->TimeFiltTime   = NULL;
  out->doDelayDecorr  = NULL;
  out->corrType = NULL;
  out->warnNoCQ = NULL;

  return out;
} /*  end newObitUVCalCalibrateS */

/**
 * Update ObitUVCalCalibrate calibration tables for time time.
 * The current table is interpolated between the previous and following
 * sets of solutions.
 * If a new set of entries is needed from the SN/CL table they are read.
 * Adopted from AIPS CGASET.FOR
 * A calibration array entry (prior or following) consists of:
 * \li real part of gain
 * \li imaginary of gain
 * \li group delay (sec)
 * \li fringe rate (sec/sec)
 * \li reference antenna for solution
 * \param in   Calibrate Object.
 * \param time desired time in days
 * \param err  Error stack for messages and errors.
 */
static void ObitUVCalCalibrateUpdate (ObitUVCalCalibrateS *in, ofloat time,
					ObitErr *err)
{
  olong iant, iif, ipol, index, jndex;
  gboolean   good1, good2, bad, doVLB;
  ofloat      wtt1, wtt2, v1r, v1i, v2r, v2i, phase1, phase2, cp1,
   cp2, sp1, sp2, g1r, g1i, g2r, g2i, delta, wt1, wt2, amp1, amp2, ampnew, ampcor;
  ofloat fblank = ObitMagicF();
  gboolean newcal;
  gchar *routine="ObitUVCalCalibrateUpdate";
 
 
  /* see if time for new table entry */
  if ((in->LastRowRead <= in->numRow)  &&  (time > in->FollowCalTime)) {
    ObitUVCalCalibrateNewTime (in, time, err);
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
  doVLB = FALSE;
  /* initialize indices for CalApply (index),  CalPrior and CalFollow (jndex)*/
  index = 0;
  jndex = 0;

  /* loop thru antennas */
  for (iant= 0; iant<in->numAnt; iant++) { /* loop 500 */

    /* set interpolation weights proportional to time difference. */
    wtt1 = 0.0;
    if (time < in->FollowAntTime[iant]) {
      if (in->FollowAntTime[iant] > in->PriorAntTime[iant]) 
	wtt1 = (in->FollowAntTime[iant] - time) / (in->FollowAntTime[iant] - in->PriorAntTime[iant]);
    } 
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

    /* set dispersive delay  interpolate if both good */
    if ((in->PriorDDelay[iant]  !=  fblank) &&  (in->FollowDDelay[iant]  !=  fblank)) {
      in->DDelay[iant] = VELIGHT * 2.0 * G_PI * (wtt1 * in->PriorDDelay[iant] + 
						 wtt2 * in->FollowDDelay[iant]);
    } else if (in->PriorDDelay[iant]  !=  fblank) {
      in->DDelay[iant] = VELIGHT * 2.0 * G_PI * in->PriorDDelay[iant];
    } else if (in->FollowDDelay[iant]  !=  fblank) {
      in->DDelay[iant] = VELIGHT * 2.0 * G_PI * in->FollowDDelay[iant];
    } else {
      in->DDelay[iant] = fblank;
    }
    
    /* loop thru IF */
    for (iif= in->bIF-1; iif<in->eIF; iif++) { /* loop 400 */

      /* loop thru polarization */
      for (ipol= 0; ipol<in->numPol; ipol++) { /* loop 300 */

	/* initialize soln with blanks */
	in->CalApply[index]   = fblank;
	in->CalApply[index+1] = fblank;
	in->CalApply[index+2] = fblank;
	in->CalApply[index+3] = fblank;
	in->CalApply[index+4] = fblank;

	/* check for blanked soln. */
	good1 = (in->CalPrior[jndex] != fblank)  && 
	  (in->CalPrior[jndex+1] != fblank)  && 
	  (in->CalPrior[jndex+2] != fblank)  && 
	  (in->CalPrior[jndex+3] != fblank)  && 
	  (wtt1 > 0.0);
	good2 = (in->CalFollow[jndex] != fblank)  && 
	  (in->CalFollow[jndex+1] != fblank)  && 
	  (in->CalFollow[jndex+2] != fblank)  && 
	  (in->CalFollow[jndex+3] != fblank)  && 
	  (wtt2 > 0.0);

	/* solution all flagged? */
	bad = !(good1 || good2);

	/* nothing more if both prior and following are bad */
	if (!bad) {

	  /* different reference antennas  use closest */
	  if ((fabs (in->CalPrior[jndex+4] - in->CalFollow[jndex+4]) >= 0.5)
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
	  
	  /* set values, initial gain. */
	  g1r = 0.0;
	  g1i = 0.0;
	  g2r = 0.0;
	  g2i = 0.0;
	  
	  /* Get gains from tables */
	  if (good1) {
	    g1r = in->CalPrior[jndex];
	    g1i = in->CalPrior[jndex+1];
	  } 
	  if (good2) {
	    g2r = in->CalFollow[jndex];
	    g2i = in->CalFollow[jndex+1];
	  } 
	  
	  v1r = g1r;
	  v1i = g1i;
	  v2r = g2r;
	  v2i = g2i;
	  
	  /* check if fringe rates given - if so add accumulated phase */
	  if ((in->CalPrior[jndex+3]  != 0.0)  || (in->CalFollow[jndex+3] != 0.0)) {
	    phase1 = 0.0;
	    phase2 = 0.0;
	    if (good1) phase1 = in->CalPrior[jndex+3]  * (time - in->PriorAntTime[iant])  * in->RateFact[iif];
	    if (good2) phase2 = in->CalFollow[jndex+3] * (time - in->FollowAntTime[iant]) * in->RateFact[iif];
	    cp1 = cos (phase1);
	    cp2 = cos (phase2);
	    sp1 = sin (phase1);
	    sp2 = sin (phase2);

	    /* rotate phase by accumulated fringe rate */
	    v1r = g1r * cp1 - g1i * sp1;
	    v1i = g1r * sp1 + g1i * cp1;
	    v2r = g2r * cp2 - g2i * sp2;
	    v2i = g2r * sp2 + g2i * cp2;
	  } 
	  
	  /* set amplitude and phase in output array */
	  in->CalApply[index]   = wt1 * v1r + wt2 * v2r;
	  in->CalApply[index+1] = wt1 * v1i + wt2 * v2i;

	  /* the rate correct can do bad things to the amplitudes */
	  amp1 = sqrt (v1r*v1r + v1i*v1i);
	  amp2 = sqrt (v2r*v2r + v2i*v2i);
	  ampnew = sqrt (in->CalApply[index]*in->CalApply[index] +
			 in->CalApply[index+1]*in->CalApply[index+1]);
	  if (ampnew < 1.0e-20) ampnew = 1.0e-20;
	  ampcor = (wt1*amp1 + wt2*amp2) / ampnew;
	  
	  /* control amplitudes */
	  in->CalApply[index]   *= ampcor;
	  in->CalApply[index+1] *= ampcor;
	  
	  /* init. delay rotation. */
	  in->CalApply[index+2] = 1.0;
	  in->CalApply[index+3] = 0.0;
	  
	  /* delay correction - factors to rotate phase from channel to channel */
	  if ((in->CalPrior[jndex+2]  != 0.0)  || (in->CalFollow[jndex+2] != 0.0)) {
	    phase1 = (wt1 * in->CalPrior[jndex+2] +
		      wt2 * in->CalFollow[jndex+2]) * in->DelayFact[iif];
	    doVLB = doVLB || (phase1 != 0.0); /* need delay/rate corrections */
	    in->CalApply[index+2] = cos (phase1);
	    in->CalApply[index+3] = sin (phase1);
	  } 

	  /* set rate */
	  in->CalApply[index+4] = (wt1 * in->CalPrior[jndex+3] +
				   wt2 * in->CalFollow[jndex+3]) * in->RateFact[iif];
	} /* end of only valid solutions section */

	/* update indices */
        index += in->lenCalArrayEntry;
	jndex += in->lenCalArrayEntry;

      } /* end poln loop  L300: */;
    } /* end IF loop  L400: */;
  } /* end antenna loop  L500: */

  /* need delay rate correction? */
  in->doDelayRate = doVLB;
} /* end ObitUVCalCalibrateUpdate */
    
/**
 * Read calibration for next time from CL or SN table.
 * Applies mean gain modulus corrections
 * Adopted from AIPS CSLGET.FOR
 * \param in   Calibrate Object.
 * \param time desired time in days
 * \param in   Error stack for messages and errors.
 */
static void ObitUVCalCalibrateNewTime (ObitUVCalCalibrateS *in, ofloat time,
					ObitErr *err)
{
  ObitIOCode retCode;
  ofloat mGModI, wt1, wt2;
  ofloat fblank = ObitMagicF();
  olong nblank, i, j, iant, iif, indx, lenEntryPoln, lenEntry, lenEntryAnt;
  olong  irow=0, limit,IFoff, antno;
  gboolean want, done;
  ObitTableSN *SNTable = NULL;
  ObitTableSNRow *SNTableRow = NULL;
  ObitTableCL *CLTable = NULL;
  ObitTableCLRow *CLTableRow = NULL;
  gchar *routine="ObitUVCalCalibrateNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* increments, sizes of data elements */
  /* length of basic entry */
  lenEntry = in->numPol * in->lenCalArrayEntry;
  /* length of an entry for an antenna */
  lenEntryAnt = lenEntry * (in->eIF - in->bIF + 1);
  /* Length of entry with all polarizations */
  lenEntryPoln = in->lenCalArrayEntry * MIN (2, MAX (1, in->numPol));
  
  /* initialize Prior and Following arrays if first call */
  if (in->LastRowRead <= 0) {
    nblank = in->numAnt * (in->eIF - in->bIF + 1) * in->numPol * in->lenCalArrayEntry;
    for (i=0; i<nblank; i++) in->CalPrior[i]  = fblank;
    for (i=0; i<nblank; i++) in->CalFollow[i] = fblank;
    for (i=0; i<in->numAnt; i++) in->PriorIFR[i]     = fblank;
    for (i=0; i<in->numAnt; i++) in->FollowIFR[i]    = fblank;
    for (i=0; i<in->numAnt; i++) in->PriorDDelay[i]  = fblank;
    for (i=0; i<in->numAnt; i++) in->FollowDDelay[i] = fblank;
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
	in->PriorDDelay[iant]  = in->FollowDDelay[iant];
	for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 20 */
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	  for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	} /* end IF loop  L20:  */;
      } 
    } /* end loop  ant L30:  */;
  }

  /* Handle SN/CL table separately but do the same things for each */
  if (in->doSNTable) {
    /* SN table  - set local pointers */
    SNTable = (ObitTableSN*)in->SNTable;
    SNTableRow = (ObitTableSNRow*)in->SNTableRow;
    
    /* mean gain modulus correction */
    if (SNTable-> mGMod>0.00001) mGModI = 1.0 / (SNTable-> mGMod);
    else mGModI = 1.0;

    /* Read through rows filling in data */
    /* read until selected time. */
    limit = MAX (1, in->LastRowRead);
    in->LastRowRead = 0;  /* The next time may start somewhere nonobvious */
    for (i= limit; i<=in->numRow; i++) { /* loop 90 */
      irow = i;
      retCode = ObitTableSNReadRow (SNTable, irow, SNTableRow, err);
      if (err->error) Obit_traceback_msg (err, routine, "Cal(SN) table");
      if (SNTableRow->status < 0) continue; /* entry flagged? */
      
      /* check subarray */
      want = ((SNTableRow->SubA == in->SubA)  ||  (SNTableRow->SubA <= 0) || (in->SubA <= 0));
      
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
	
	if (in->PriorAntTime[iant] > -100.) {
	  /* new folling entry - copy to prior */
	  in->PriorAntTime[iant] = in->FollowAntTime[iant];
	  in->PriorIFR[iant]     = in->FollowIFR[iant];
	  in->PriorDDelay[iant]  = in->FollowDDelay[iant];
	  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 50 */
	    indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	    for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	  } /* end IF loop  L50:  */;
	}
	
	/* fill in new following values */
	in->FollowIFR[iant]     = SNTableRow->IFR;
	in->FollowDDelay[iant]  = 0.0; /* no dispersive delay in SN table */
	in->FollowAntTime[iant] = SNTableRow->Time;
	
	/* loop over if */
	for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 60 */
	  IFoff = iif - 1;
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	  wt1 = SNTableRow->Weight1[IFoff];
	  /* Trap zero amplitudes */
	  if (((SNTableRow->Real1[IFoff]*SNTableRow->Real1[IFoff]) + 
	       (SNTableRow->Imag1[IFoff]*SNTableRow->Imag1[IFoff])) < 1.0e-10)
	    wt1 = -1.0;
	  in->CalFollow[indx]   = SNTableRow->Real1[IFoff];
	  in->CalFollow[indx+1] = SNTableRow->Imag1[IFoff];
	  in->CalFollow[indx+2] = SNTableRow->Delay1[IFoff];
	  in->CalFollow[indx+3] = SNTableRow->Rate1[IFoff];
	  in->CalFollow[indx+4] = SNTableRow->RefAnt1[IFoff];
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
	    /* Trap zero amplitudes */
	    if (((SNTableRow->Real2[IFoff]*SNTableRow->Real2[IFoff]) + 
		 (SNTableRow->Imag2[IFoff]*SNTableRow->Imag2[IFoff])) < 1.0e-10)
	    wt2 = -1.0;
	    in->CalFollow[indx]   = SNTableRow->Real2[IFoff];
	    in->CalFollow[indx+1] = SNTableRow->Imag2[IFoff];
	    in->CalFollow[indx+2] = SNTableRow->Delay2[IFoff];
	    in->CalFollow[indx+3] = SNTableRow->Rate2[IFoff];
	    in->CalFollow[indx+4] = SNTableRow->RefAnt2[IFoff];
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
	  in->PriorDDelay[iant]  = in->FollowDDelay[iant];
	  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 70 */
	    indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	    for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	  } /* end IF loop  L70:  */;
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
	     if ((time >= in->FollowAntTime[iant]) && (in->PriorAntTime[iant] >= -100.0)) done=FALSE;
	} /* end loop  L80:  */;
	
	/* no more to fill in */
	if (done) break;
	 } 
    } /* end loop over table entries L90:  */;

  } else {
    /* CL table - set local pointers */
    CLTable = (ObitTableCL*)in->CLTable;
    CLTableRow = (ObitTableCLRow*)in->CLTableRow;
    
    /* mean gain modulus correction */
    if (CLTable-> mGMod>0.00001) mGModI = 1.0 / (CLTable-> mGMod);
    else mGModI = 1.0;

    /* Read through rows filling in data - read until selected time. */
    limit = MAX (1, in->LastRowRead);
    in->LastRowRead = 0;  /* The next time may start somewhere nonobvious */
    for (i= limit; i<=in->numRow; i++) { /* loop 90 */
      irow = i;
      retCode = ObitTableCLReadRow (CLTable, irow, CLTableRow, err);
      if (err->error) Obit_traceback_msg (err, routine, "Cal(CL) Table");
      if (CLTableRow->status < 0) continue; /* entry flagged? */
      
      /* check subarray */
      want = ((CLTableRow->SubA == in->SubA)  ||  (CLTableRow->SubA <= 0) || (in->SubA <= 0));
      
      /* check frqsel */
      want = want &&
	((CLTableRow->FreqID == in->FreqID) || (CLTableRow->FreqID <= 0) ||
	 (in->FreqID <= 0));
      
      /* skip if not wanted */
      if (!want) continue;
      
      /* antenna number */
      antno = CLTableRow->antNo;
      iant = antno-1;
      
      /* time -> include this one? */
      if (time >= in->FollowAntTime[iant]) { 
	
	if (in->PriorAntTime[iant] > -100.) {
	  /* new folling entry - copy to prior */
	  in->PriorAntTime[iant] = in->FollowAntTime[iant];
	  in->PriorIFR[iant]     = in->FollowIFR[iant];
	  in->PriorDDelay[iant]  = in->FollowDDelay[iant];
	  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 50 */
	    indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	    for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	  } /* end IF loop  L50:  */;
	}
	
	/* fill in new following values */
	in->FollowIFR[iant]     = CLTableRow->IFR;
	in->FollowDDelay[iant]  = CLTableRow->dispers1;
	in->FollowAntTime[iant] = CLTableRow->Time;
	
	/* loop over if */
	for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 60 */
	  IFoff = iif - 1;
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	  wt1 = CLTableRow->Weight1[IFoff];
	  in->CalFollow[indx]   = CLTableRow->Real1[IFoff];
	  in->CalFollow[indx+1] = CLTableRow->Imag1[IFoff];
	  in->CalFollow[indx+2] = CLTableRow->Delay1[IFoff];
	  in->CalFollow[indx+3] = CLTableRow->Rate1[IFoff];
	  in->CalFollow[indx+4] = CLTableRow->RefAnt1[IFoff];
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
	    wt2 = CLTableRow->Weight2[IFoff];
	    in->CalFollow[indx]   = CLTableRow->Real2[IFoff];
	    in->CalFollow[indx+1] = CLTableRow->Imag2[IFoff];
	    in->CalFollow[indx+2] = CLTableRow->Delay2[IFoff];
	    in->CalFollow[indx+3] = CLTableRow->Rate2[IFoff];
	    in->CalFollow[indx+4] = CLTableRow->RefAnt2[IFoff];
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
	  in->PriorDDelay[iant]  = in->FollowDDelay[iant];
	  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 70 */
	    indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	    for (j=0; j<lenEntryPoln; j++) in->CalPrior[indx+j]  = in->CalFollow[indx+j];
	  } /* end IF loop  L70:  */;
	} /* end copy to Prior */
	
      } else {
	
	/* This one not needed time> wanted - are we there yet? */
	/* May need to restart earlier in table for some antennas.
	   Remember the first record not used. */

	if (in->LastRowRead <= 0) in->LastRowRead = i;
	
	/* if desired time for some antenna still after the Following time and there is no
	   Prior entry, keep going */
	done = TRUE;
	for (antno= 1; antno<=in->numAnt; antno++) { /* loop 80 */
	     iant = antno-1;
	     if ((time >= in->FollowAntTime[iant]) && (in->PriorAntTime[iant] >= -100.0)) done=FALSE;
	} /* end loop  L80:  */;
	
	/* no more to fill in */
	if (done) break;
	 } 
    } /* end loop over table entries L90:  */

  } /* end of do CL table */
  
  /* If all entries have been read and some antennas do not have a 
     following entry make dummy blanked entry at the last time */
  if (irow>=in->numRow ) {
    for (antno= 1; antno<=in->numAnt; antno++) { /* loop 80 */
      iant = antno-1;
      if ((time >= in->FollowAntTime[iant]) && (in->PriorAntTime[iant] >= -100.0)) {
	/* Dummy Following entry */
	in->FollowAntTime[iant] += 10.0;
	/* loop over if */
	for (iif= in->bIF; iif<=in->eIF; iif++) { /* IF loop */
	  IFoff = iif - 1;
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	  in->CalFollow[indx]   = fblank;
	  in->CalFollow[indx+1] = fblank;
	  in->CalFollow[indx+2] = fblank;
	  in->CalFollow[indx+3] = fblank;
	  in->CalFollow[indx+4] = 0;
	  
	  /* second polarization if present */
	  if (in->numPol >= 2) {
	    indx = indx + in->lenCalArrayEntry;
	    in->CalFollow[indx]   = fblank;
	    in->CalFollow[indx+1] = fblank;
	    in->CalFollow[indx+2] = fblank;
	    in->CalFollow[indx+3] = fblank;
	    in->CalFollow[indx+4] = 0.0;
	  } /* end second poln */
	} /* end IF loop:  */
	
      }
    } /* end checking antennas  */;
  } /* end if got to end */
  
  /* finished file using all entries? */
  if (in->LastRowRead <= 0) in->LastRowRead = in->numRow + 1;
  
  /* Set times */
  in->FollowCalTime = 1.0e10;
  in->PriorCalTime = -1.0e10;
  for (antno= 1; antno<=in->numAnt; antno++) { /* loop 110 */
    iant = antno -1;
    if (in->PriorAntTime[iant] >= -100.0) {
      if (time >= in->PriorAntTime[iant]) 
	in->PriorCalTime = MAX (in->PriorCalTime, in->PriorAntTime[iant]);
      if (time <= in->FollowAntTime[iant]) 
	in->FollowCalTime = MIN (in->FollowCalTime,in->FollowAntTime[iant]);
    } 
  } /* end loop  L110: */;
  
  /* just to be sure something rational in times */
  if (in->PriorCalTime < -1000.0)  in->PriorCalTime  = time - 2.0/86400.0;
  if (in->FollowCalTime > 10000.0) in->FollowCalTime = time + 2.0/86400.0;
  
} /* end ObitUVCalCalibrateNewTime */

/**
 * If data is from the VLBA, read the CQ table and prepare corrections.
 * The CQ table is very poorly documented but the row number in the CQ 
 * table appears to correspond to a CORR-ID number.
 * Adopted from AIPS DSMEAR.FOR
 * \param in     Calibrate object
 * \param out    Structure to update
 * \param err    ObitError stack
 */
static void
ObitUVCalCalibrateVLBAInit (ObitUVCal *in, ObitUVCalCalibrateS *out, 
			      ObitErr *err)
{
  ObitIOCode retCode;
  olong  numCID, size, i;
  ObitTableCQRow *CQTableRow=NULL;
  ObitTableCQ *CQTable=NULL;
  ObitUVDesc *desc=NULL;
  odouble dfact;
  olong indx, jif, jfact, isub, icorr, ifilt, id, ii, rc, ioff;
  gchar strTemp[9];
  gchar *routine="ObitUVCalCalibrateVLBAInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* Create arrays that depend on number of subarrays */
  if (out->corrType) g_free(out->corrType);
  out->corrType = g_malloc0(in->numSubA*sizeof(oint));
  if (out->warnNoCQ) g_free(out->warnNoCQ);
  out->warnNoCQ = g_malloc0(in->numSubA*sizeof(gboolean));

  /* Initialize Subarray dependent values */
  for (i= 0; i< in->numSubA; i++) { 
    out->corrType[i] = 0;
    out->warnNoCQ[i] = FALSE;
  }

  desc = in->myDesc; /* uv data descriptor */
  
  /* Determine which subarrays contain VLBA data from AN tables */
  
  /* loop over each subarray */
  for (i= 0; i<in->numANTable; i++) { /* loop 50 */
    
    /* Open table */
    retCode = ObitTableANOpen ((ObitTableAN*)in->ANTables[i], OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
    
    /* Is this VLBA data? */
    if (!strncmp (((ObitTableAN*)in->ANTables[i])->ArrName, "VLBA", 4)) {
      /* is VLBA */
      out->corrType[i] = 1;
      out->warnNoCQ[i] = TRUE;
    } else {
      /* Not VLBA */
      out->corrType[i] = 0;
      out->warnNoCQ[i] = FALSE;
    }
    /* Close table */
    retCode = ObitTableANClose ((ObitTableAN*)in->ANTables[i], err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
  } /* end loop over subarrays */
  
  /* If there is no CQ table you're done */
  if (in->CQTable==NULL) return;
  
  /* Read cq table information */
  /* Open table */
  CQTable = (ObitTableCQ*)in->CQTable;
  retCode = ObitTableCQOpen (CQTable, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);
  
  /* How many entries in table? */
  numCID = CQTable->myDesc->nrow;
	    
  /* (re)Create Corr-id specific arrays */
  size = out->numIF * numCID;
  out->LTaper       = g_realloc(out->LTaper, size*sizeof(olong));
  out->TimeFiltTime = g_realloc(out->TimeFiltTime, size*sizeof(ofloat));
  out->NSpecA       = g_realloc(out->NSpecA, size*sizeof(olong));
  out->DelBit       = g_realloc(out->DelBit, size*sizeof(odouble));
  out->NFFTSize     = g_realloc(out->NFFTSize, size*sizeof(olong));
  out->typeTimeFilt = g_realloc(out->typeTimeFilt, size*sizeof(olong));
  out->doDelayDecorr= g_realloc(out->doDelayDecorr, size*sizeof(gboolean));
  
  /* Initialize values */
  for (i= 0; i< size; i++) {
    out->LTaper[i]       = 1;
    out->TimeFiltTime[i] = 0.0;
    out->NSpecA[i]       = 0;
    out->DelBit[i]       = 0.0;
    out->NFFTSize[i]     = 0;
    out->typeTimeFilt[i] = 0;
    out->doDelayDecorr[i]= FALSE;
  }
  
  /* make row structure */
  CQTableRow = newObitTableCQRow(CQTable);
  
  /* Loop over table */
  for (i= 1; i<= numCID; i++) { /* loop 200 */
    
    retCode = ObitTableCQReadRow (CQTable, i, CQTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "CQ Table");
    if (CQTableRow->status < 0) continue; /* entry flagged? */
    
    /* Select on fq_id */
    out->FreqID = MAX (1, out->FreqID);
    if (CQTableRow->FrqSel == out->FreqID) { /*goto 200;*/
      
      /*  Fill output arrays */
      isub = MAX (CQTableRow->SubA, 1); /* which subarray */
      
      for (jif= 1; jif<= CQTable->numIF; jif++) { /* loop 150 */
	/* extract correlation  id. from VLBA encryption scheme */
	id = MAX (CQTableRow->Filter[jif+1], 0);
	ifilt = (id % 256);
	icorr = (id - ifilt) / 256 + 1.1;

	indx = (icorr - 1) * out->numIF + jif;
	
	/* Taper function */
	/* Have to handle the array of strings here - each IF has 8 characters */
	out->LTaper[indx-1] = 2; /* default "Uniform" */
	ioff = (jif-1)*8; /* select for subarray */
	for (ii = 0; ii<8; ii++) strTemp[ii] = CQTableRow->TaperFn[ioff+ii];
	if (!strncmp (strTemp, "HANNING", 7)) out->LTaper[indx-1] = 1;
	
	/* Spectral avg. function */
	out->NSpecA[indx-1] = CQTableRow->SpecAvg[jif-1];
	
	/* Scale by current channel bandwidth in case data have been averaged 
	   since correlation */
	if (fabs (CQTableRow->ChanBW[jif-1]) > 0.0) {
	  jfact = (olong) (fabs (desc->chIncIF[jif-1] / CQTableRow->ChanBW[jif-1]) + 0.1);
	} else {
	  jfact = 1;
	} 
	out->NSpecA[indx-1] = out->NSpecA[indx-1] * jfact;
	
	/* FFT size */
	out->NFFTSize[indx-1] = CQTableRow->FFTSize[jif-1];
	
	/* Time interval per bit */
	if ((abs(CQTableRow->ChanBW[jif-1]) > 0.0) && (CQTableRow->numChan[jif-1] > 0)) {
	  out->DelBit[indx-1] = 1.0 / (2 * CQTableRow->ChanBW[jif-1] * 
				      CQTableRow->numChan[jif-1]);
	} else {
	  out->DelBit[indx-1] = -1.0;
	} 
	
	/* filter type */
	out->typeTimeFilt[indx-1] = ifilt;
	
	/*  Filter averaging time */
	out->TimeFiltTime[indx-1] = CQTableRow->TimeAvg[jif-1];
	
	/* Can delay decorrelation corrections be done for this 
	   (if,corr_id) combination ? Are taper function and fft size valid?*/
	
	dfact = 
	  ObitUVCalCalibrateSegLoss (out->LTaper[indx-1], out->NFFTSize[indx-1], 0.0, 
				     &rc);

	/* incorporate other conditions */
	out->doDelayDecorr[indx-1] = ((out->NSpecA[indx-1] > 0) && 
				     (out->DelBit[indx-1] > 0.0) && (rc == 0));
      } /* end loop L150:  */
    } /* end of FQId selected */
  } /* end loop L200:  */
  
  /* close CQ table */
  retCode = ObitTableCQClose (CQTable, err);
  if (err->error) Obit_traceback_msg (err, routine, "CQ Table");
  
  /* release row structure */
  CQTableRow = ObitTableCQRowUnref(CQTableRow);
} /*  end ObitUVCalCalibrateVLBAInit */


/**
 * Calculate the VLBA segmentation loss factor.
 * Adopted from AIPS FXSEG.FOR
 * \param lFunc  FFT weighting function (1='UNIFORM', 2='HANNING')
 * \param nfft   FFT size (in bits)
 * \param dbits  Delay error (in bits)
 * \param rc     Return code 0=> everything was wonderful.
 * \return the loss correction factor.
 */
static double 
ObitUVCalCalibrateSegLoss (olong lFunc, olong nfft, odouble dbits, olong *rc)
{
  gboolean wvalid;
  odouble dfact=0.0, dfrac, dslope, doffst;
  olong ibit, jbit, intbit;

  /* tables by correlator FFT size */
  /* 2048 FFT */
  odouble d2048[2048] = {
    0.10000000e+01, 0.99994427e+00,
    0.99988902e+00, 0.99983341e+00, 0.99977815e+00, 0.99972254e+00,
    0.99966729e+00, 0.99961162e+00, 0.99955618e+00, 0.99949700e+00,
    0.99943775e+00, 0.99937856e+00, 0.99931914e+00, 0.99925995e+00,
    0.99920052e+00, 0.99914134e+00, 0.99908209e+00, 0.99901313e+00,
    0.99894410e+00, 0.99887526e+00, 0.99880606e+00, 0.99873722e+00,
    0.99866825e+00, 0.99859923e+00, 0.99853021e+00, 0.99843866e+00,
    0.99834728e+00, 0.99825549e+00, 0.99816412e+00, 0.99807256e+00,
    0.99798101e+00, 0.99788946e+00, 0.99779803e+00, 0.99769562e+00,
    0.99759299e+00, 0.99749058e+00, 0.99738795e+00, 0.99728549e+00,
    0.99718291e+00, 0.99708045e+00, 0.99697787e+00, 0.99684972e+00,
    0.99672157e+00, 0.99659342e+00, 0.99646527e+00, 0.99633718e+00,
    0.99620903e+00, 0.99608088e+00, 0.99595273e+00, 0.99578679e+00,
    0.99562061e+00, 0.99545467e+00, 0.99528873e+00, 0.99512255e+00,
    0.99495661e+00, 0.99479043e+00, 0.99462450e+00, 0.99443334e+00,
    0.99424261e+00, 0.99405146e+00, 0.99386090e+00, 0.99366957e+00,
    0.99347895e+00, 0.99328768e+00, 0.99309707e+00, 0.99287850e+00,
    0.99266058e+00, 0.99244195e+00, 0.99222422e+00, 0.99200559e+00,
    0.99178785e+00, 0.99156922e+00, 0.99135137e+00, 0.99111664e+00,
    0.99088228e+00, 0.99064761e+00,
    0.99041307e+00, 0.99017835e+00, 0.98994380e+00, 0.98970914e+00,
    0.98947459e+00, 0.98922569e+00, 0.98897660e+00, 0.98872751e+00,
    0.98847842e+00, 0.98822951e+00, 0.98798043e+00, 0.98773134e+00,
    0.98748231e+00, 0.98720795e+00, 0.98693359e+00, 0.98665929e+00,
    0.98638493e+00, 0.98611063e+00, 0.98583627e+00, 0.98556197e+00,
    0.98528761e+00, 0.98496878e+00, 0.98465002e+00, 0.98433119e+00,
    0.98401254e+00, 0.98369372e+00, 0.98337507e+00, 0.98305631e+00,
    0.98273766e+00, 0.98240006e+00, 0.98206264e+00, 0.98172522e+00,
    0.98138779e+00, 0.98105019e+00, 0.98071283e+00, 0.98037523e+00,
    0.98003799e+00, 0.97968411e+00, 0.97933078e+00, 0.97897696e+00,
    0.97862345e+00, 0.97826976e+00, 0.97791630e+00, 0.97756243e+00,
    0.97720891e+00, 0.97682083e+00, 0.97643274e+00, 0.97604465e+00,
    0.97565657e+00, 0.97526848e+00, 0.97488034e+00, 0.97449225e+00,
    0.97410417e+00, 0.97370684e+00, 0.97330964e+00, 0.97291231e+00,
    0.97251517e+00, 0.97211778e+00, 0.97172046e+00, 0.97132331e+00,
    0.97092593e+00, 0.97049928e+00, 0.97007263e+00, 0.96964604e+00,
    0.96921939e+00, 0.96879274e+00, 0.96836627e+00, 0.96793944e+00,
    0.96751297e+00, 0.96706712e+00, 0.96662176e+00, 0.96617591e+00,
    0.96573061e+00, 0.96528471e+00,
    0.96483922e+00, 0.96439350e+00, 0.96394801e+00, 0.96348166e+00,
    0.96301526e+00, 0.96254891e+00, 0.96208256e+00, 0.96161622e+00,
    0.96114969e+00, 0.96068347e+00, 0.96021694e+00, 0.95972532e+00,
    0.95923358e+00, 0.95874196e+00, 0.95825034e+00, 0.95775872e+00,
    0.95726693e+00, 0.95677531e+00, 0.95628375e+00, 0.95575446e+00,
    0.95522505e+00, 0.95469576e+00, 0.95416635e+00, 0.95363706e+00,
    0.95310766e+00, 0.95257837e+00, 0.95204896e+00, 0.95151150e+00,
    0.95097387e+00, 0.95043647e+00, 0.94989884e+00, 0.94936138e+00,
    0.94882381e+00, 0.94828635e+00, 0.94774872e+00, 0.94719344e+00,
    0.94663811e+00, 0.94608277e+00, 0.94552749e+00, 0.94497234e+00,
    0.94441700e+00, 0.94386172e+00, 0.94330639e+00, 0.94271874e+00,
    0.94213104e+00, 0.94154340e+00, 0.94095594e+00, 0.94036824e+00,
    0.93978059e+00, 0.93919289e+00, 0.93860525e+00, 0.93801200e+00,
    0.93741840e+00, 0.93682516e+00, 0.93623155e+00, 0.93563837e+00,
    0.93504477e+00, 0.93445152e+00, 0.93385792e+00, 0.93324310e+00,
    0.93262798e+00, 0.93201298e+00, 0.93139780e+00, 0.93078303e+00,
    0.93016785e+00, 0.92955285e+00, 0.92893773e+00, 0.92828751e+00,
    0.92763728e+00, 0.92698693e+00, 0.92633671e+00, 0.92568648e+00,
    0.92503631e+00, 0.92438608e+00,
    0.92373568e+00, 0.92307073e+00, 0.92240584e+00, 0.92174089e+00,
    0.92107594e+00, 0.92041081e+00, 0.91974604e+00, 0.91908091e+00,
    0.91841596e+00, 0.91773117e+00, 0.91704583e+00, 0.91636103e+00,
    0.91567588e+00, 0.91499090e+00, 0.91430575e+00, 0.91362077e+00,
    0.91293561e+00, 0.91224444e+00, 0.91155326e+00, 0.91086221e+00,
    0.91017103e+00, 0.90948004e+00, 0.90878886e+00, 0.90809768e+00,
    0.90740669e+00, 0.90669775e+00, 0.90598893e+00, 0.90528011e+00,
    0.90457135e+00, 0.90386254e+00, 0.90315378e+00, 0.90244496e+00,
    0.90173620e+00, 0.90099996e+00, 0.90026397e+00, 0.89952773e+00,
    0.89879173e+00, 0.89805555e+00, 0.89731950e+00, 0.89658332e+00,
    0.89584732e+00, 0.89508945e+00, 0.89433175e+00, 0.89357394e+00,
    0.89281625e+00, 0.89205837e+00, 0.89130056e+00, 0.89054286e+00,
    0.88978499e+00, 0.88900971e+00, 0.88823444e+00, 0.88745916e+00,
    0.88668388e+00, 0.88590854e+00, 0.88513309e+00, 0.88435781e+00,
    0.88358253e+00, 0.88279241e+00, 0.88200229e+00, 0.88121200e+00,
    0.88042188e+00, 0.87963176e+00, 0.87884170e+00, 0.87805158e+00,
    0.87726146e+00, 0.87644976e+00, 0.87563807e+00, 0.87482649e+00,
    0.87401479e+00, 0.87320316e+00, 0.87239152e+00, 0.87157989e+00,
    0.87076819e+00, 0.86994970e+00,
    0.86913097e+00, 0.86831248e+00, 0.86749381e+00, 0.86667514e+00,
    0.86585641e+00, 0.86503792e+00, 0.86421925e+00, 0.86338556e+00,
    0.86255169e+00, 0.86171776e+00, 0.86088407e+00, 0.86005020e+00,
    0.85921651e+00, 0.85838264e+00, 0.85754877e+00, 0.85668248e+00,
    0.85581619e+00, 0.85494965e+00, 0.85408336e+00, 0.85321689e+00,
    0.85235053e+00, 0.85148424e+00, 0.85061777e+00, 0.84974444e+00,
    0.84887111e+00, 0.84799778e+00, 0.84712446e+00, 0.84625095e+00,
    0.84537762e+00, 0.84450430e+00, 0.84363079e+00, 0.84275550e+00,
    0.84188020e+00, 0.84100491e+00, 0.84012961e+00, 0.83925432e+00,
    0.83837903e+00, 0.83750373e+00, 0.83662844e+00, 0.83573407e+00,
    0.83483976e+00, 0.83394539e+00, 0.83305103e+00, 0.83215672e+00,
    0.83126235e+00, 0.83036816e+00, 0.82947367e+00, 0.82857031e+00,
    0.82766736e+00, 0.82676405e+00, 0.82586086e+00, 0.82495755e+00,
    0.82405454e+00, 0.82315123e+00, 0.82224804e+00, 0.82132047e+00,
    0.82039326e+00, 0.81946564e+00, 0.81853843e+00, 0.81761080e+00,
    0.81668341e+00, 0.81575578e+00, 0.81482857e+00, 0.81388420e+00,
    0.81293976e+00, 0.81199539e+00, 0.81105101e+00, 0.81010658e+00,
    0.80916220e+00, 0.80821776e+00, 0.80727339e+00, 0.80632144e+00,
    0.80536950e+00, 0.80441737e+00,
    0.80346560e+00, 0.80251348e+00, 0.80156153e+00, 0.80060941e+00,
    0.79965764e+00, 0.79869395e+00, 0.79773039e+00, 0.79676670e+00,
    0.79580295e+00, 0.79483944e+00, 0.79387569e+00, 0.79291219e+00,
    0.79194826e+00, 0.79096568e+00, 0.78998309e+00, 0.78900033e+00,
    0.78801775e+00, 0.78703499e+00, 0.78605241e+00, 0.78506964e+00,
    0.78408700e+00, 0.78310865e+00, 0.78213012e+00, 0.78115159e+00,
    0.78017306e+00, 0.77919465e+00, 0.77821612e+00, 0.77723759e+00,
    0.77625906e+00, 0.77526593e+00, 0.77427268e+00, 0.77327937e+00,
    0.77228630e+00, 0.77129298e+00, 0.77029973e+00, 0.76930660e+00,
    0.76831335e+00, 0.76730669e+00, 0.76630020e+00, 0.76529354e+00,
    0.76428688e+00, 0.76328033e+00, 0.76227367e+00, 0.76126701e+00,
    0.76026034e+00, 0.75924695e+00, 0.75823337e+00, 0.75721979e+00,
    0.75620615e+00, 0.75519276e+00, 0.75417918e+00, 0.75316554e+00,
    0.75215209e+00, 0.75112838e+00, 0.75010461e+00, 0.74908090e+00,
    0.74805725e+00, 0.74703348e+00, 0.74600977e+00, 0.74498600e+00,
    0.74396235e+00, 0.74292850e+00, 0.74189472e+00, 0.74086100e+00,
    0.73982722e+00, 0.73879343e+00, 0.73775971e+00, 0.73672593e+00,
    0.73569214e+00, 0.73464906e+00, 0.73360598e+00, 0.73256296e+00,
    0.73151982e+00, 0.73047674e+00,
    0.72943372e+00, 0.72839063e+00, 0.72734761e+00, 0.72630113e+00,
    0.72525471e+00, 0.72420835e+00, 0.72316194e+00, 0.72211552e+00,
    0.72106910e+00, 0.72002274e+00, 0.71897626e+00, 0.71792144e+00,
    0.71686667e+00, 0.71581179e+00, 0.71475703e+00, 0.71370214e+00,
    0.71264732e+00, 0.71159244e+00, 0.71053767e+00, 0.70947421e+00,
    0.70841074e+00, 0.70734739e+00, 0.70628393e+00, 0.70522046e+00,
    0.70415699e+00, 0.70309365e+00, 0.70203018e+00, 0.70095628e+00,
    0.69988263e+00, 0.69880879e+00, 0.69773501e+00, 0.69666123e+00,
    0.69558746e+00, 0.69451362e+00, 0.69343996e+00, 0.69235629e+00,
    0.69127262e+00, 0.69018888e+00, 0.68910521e+00, 0.68802154e+00,
    0.68693787e+00, 0.68585420e+00, 0.68477046e+00, 0.68369186e+00,
    0.68261325e+00, 0.68153459e+00, 0.68045598e+00, 0.67937738e+00,
    0.67829877e+00, 0.67722011e+00, 0.67614144e+00, 0.67506242e+00,
    0.67398334e+00, 0.67290425e+00, 0.67182517e+00, 0.67074603e+00,
    0.66966701e+00, 0.66858792e+00, 0.66750878e+00, 0.66642326e+00,
    0.66533768e+00, 0.66425204e+00, 0.66316646e+00, 0.66208088e+00,
    0.66099536e+00, 0.65990973e+00, 0.65882415e+00, 0.65773535e+00,
    0.65664655e+00, 0.65555775e+00, 0.65446901e+00, 0.65338016e+00,
    0.65229142e+00, 0.65120262e+00,
    0.65011382e+00, 0.64901221e+00, 0.64791059e+00, 0.64680898e+00,
    0.64570731e+00, 0.64460570e+00, 0.64350414e+00, 0.64240247e+00,
    0.64130080e+00, 0.64019197e+00, 0.63908327e+00, 0.63797438e+00,
    0.63686556e+00, 0.63575673e+00, 0.63464797e+00, 0.63353914e+00,
    0.63243032e+00, 0.63131821e+00, 0.63020599e+00, 0.62909389e+00,
    0.62798166e+00, 0.62686956e+00, 0.62575746e+00, 0.62464523e+00,
    0.62353313e+00, 0.62243205e+00, 0.62133098e+00, 0.62022990e+00,
    0.61912888e+00, 0.61802781e+00, 0.61692685e+00, 0.61582577e+00,
    0.61472470e+00, 0.61362582e+00, 0.61252695e+00, 0.61142814e+00,
    0.61032927e+00, 0.60923040e+00, 0.60813153e+00, 0.60703266e+00,
    0.60593379e+00, 0.60482615e+00, 0.60371840e+00, 0.60261071e+00,
    0.60150301e+00, 0.60039538e+00, 0.59928763e+00, 0.59817994e+00,
    0.59707230e+00, 0.59596336e+00, 0.59485435e+00, 0.59374541e+00,
    0.59263647e+00, 0.59152758e+00, 0.59041852e+00, 0.58930963e+00,
    0.58820063e+00, 0.58708453e+00, 0.58596843e+00, 0.58485228e+00,
    0.58373624e+00, 0.58262014e+00, 0.58150399e+00, 0.58038789e+00,
    0.57927179e+00, 0.57816023e+00, 0.57704878e+00, 0.57593721e+00,
    0.57482570e+00, 0.57371420e+00, 0.57260269e+00, 0.57149112e+00,
    0.57037961e+00, 0.56927633e+00,
    0.56817305e+00, 0.56706983e+00, 0.56596655e+00, 0.56486326e+00,
    0.56376004e+00, 0.56265676e+00, 0.56155348e+00, 0.56044459e+00,
    0.55933565e+00, 0.55822676e+00, 0.55711776e+00, 0.55600888e+00,
    0.55489993e+00, 0.55379105e+00, 0.55268210e+00, 0.55158085e+00,
    0.55047965e+00, 0.54937845e+00, 0.54827720e+00, 0.54717594e+00,
    0.54607481e+00, 0.54497355e+00, 0.54387230e+00, 0.54277122e+00,
    0.54167026e+00, 0.54056919e+00, 0.53946817e+00, 0.53836709e+00,
    0.53726614e+00, 0.53616506e+00, 0.53506398e+00, 0.53395724e+00,
    0.53285044e+00, 0.53174371e+00, 0.53063691e+00, 0.52953017e+00,
    0.52842337e+00, 0.52731663e+00, 0.52620989e+00, 0.52510196e+00,
    0.52399403e+00, 0.52288610e+00, 0.52177817e+00, 0.52067024e+00,
    0.51956230e+00, 0.51845443e+00, 0.51734650e+00, 0.51625043e+00,
    0.51515436e+00, 0.51405829e+00, 0.51296222e+00, 0.51186609e+00,
    0.51077002e+00, 0.50967395e+00, 0.50857788e+00, 0.50747937e+00,
    0.50638086e+00, 0.50528240e+00, 0.50418389e+00, 0.50308537e+00,
    0.50198686e+00, 0.50088835e+00, 0.49978986e+00, 0.49869865e+00,
    0.49760741e+00, 0.49651620e+00, 0.49542508e+00, 0.49433383e+00,
    0.49324262e+00, 0.49215138e+00, 0.49106026e+00, 0.48997891e+00,
    0.48889756e+00, 0.48781624e+00,
    0.48673499e+00, 0.48565364e+00, 0.48457229e+00, 0.48349097e+00,
    0.48240972e+00, 0.48132506e+00, 0.48024046e+00, 0.47915581e+00,
    0.47807124e+00, 0.47698665e+00, 0.47590199e+00, 0.47481743e+00,
    0.47373277e+00, 0.47266257e+00, 0.47159237e+00, 0.47052225e+00,
    0.46945205e+00, 0.46838185e+00, 0.46731165e+00, 0.46624145e+00,
    0.46517128e+00, 0.46409854e+00, 0.46302584e+00, 0.46195313e+00,
    0.46088040e+00, 0.45980769e+00, 0.45873499e+00, 0.45766228e+00,
    0.45658955e+00, 0.45551431e+00, 0.45443910e+00, 0.45336387e+00,
    0.45228863e+00, 0.45121342e+00, 0.45013818e+00, 0.44906294e+00,
    0.44798771e+00, 0.44692731e+00, 0.44586691e+00, 0.44480652e+00,
    0.44374612e+00, 0.44268569e+00, 0.44162530e+00, 0.44056490e+00,
    0.43950450e+00, 0.43845901e+00, 0.43741342e+00, 0.43636784e+00,
    0.43532237e+00, 0.43427679e+00, 0.43323120e+00, 0.43218562e+00,
    0.43114015e+00, 0.43009260e+00, 0.42904505e+00, 0.42799750e+00,
    0.42695001e+00, 0.42590246e+00, 0.42485490e+00, 0.42380735e+00,
    0.42275989e+00, 0.42171720e+00, 0.42067468e+00, 0.41963196e+00,
    0.41858944e+00, 0.41754675e+00, 0.41650423e+00, 0.41546163e+00,
    0.41441900e+00, 0.41338030e+00, 0.41234159e+00, 0.41130298e+00,
    0.41026428e+00, 0.40922558e+00,
    0.40818688e+00, 0.40714818e+00, 0.40610948e+00, 0.40508026e+00,
    0.40405104e+00, 0.40302181e+00, 0.40199259e+00, 0.40096337e+00,
    0.39993414e+00, 0.39890492e+00, 0.39787570e+00, 0.39685115e+00,
    0.39582658e+00, 0.39480203e+00, 0.39377749e+00, 0.39275295e+00,
    0.39172840e+00, 0.39070383e+00, 0.38967928e+00, 0.38867667e+00,
    0.38767403e+00, 0.38667142e+00, 0.38566878e+00, 0.38466617e+00,
    0.38366354e+00, 0.38266090e+00, 0.38165829e+00, 0.38065395e+00,
    0.37964961e+00, 0.37864530e+00, 0.37764096e+00, 0.37663662e+00,
    0.37563229e+00, 0.37462795e+00, 0.37362364e+00, 0.37263563e+00,
    0.37164757e+00, 0.37065959e+00, 0.36967152e+00, 0.36868352e+00,
    0.36769548e+00, 0.36670747e+00, 0.36571944e+00, 0.36473054e+00,
    0.36374161e+00, 0.36275271e+00, 0.36176375e+00, 0.36077484e+00,
    0.35978591e+00, 0.35879701e+00, 0.35780805e+00, 0.35682118e+00,
    0.35583431e+00, 0.35484743e+00, 0.35386059e+00, 0.35287371e+00,
    0.35188684e+00, 0.35089996e+00, 0.34991309e+00, 0.34894839e+00,
    0.34798378e+00, 0.34701908e+00, 0.34605440e+00, 0.34508970e+00,
    0.34412509e+00, 0.34316039e+00, 0.34219572e+00, 0.34124467e+00,
    0.34029362e+00, 0.33934256e+00, 0.33839151e+00, 0.33744046e+00,
    0.33648944e+00, 0.33553839e+00,
    0.33458734e+00, 0.33363882e+00, 0.33269027e+00, 0.33174175e+00,
    0.33079320e+00, 0.32984468e+00, 0.32889614e+00, 0.32794762e+00,
    0.32699910e+00, 0.32606134e+00, 0.32512358e+00, 0.32418585e+00,
    0.32324809e+00, 0.32231033e+00, 0.32137260e+00, 0.32043484e+00,
    0.31949711e+00, 0.31856760e+00, 0.31763813e+00, 0.31670865e+00,
    0.31577918e+00, 0.31484967e+00, 0.31392020e+00, 0.31299073e+00,
    0.31206125e+00, 0.31113696e+00, 0.31021270e+00, 0.30928844e+00,
    0.30836415e+00, 0.30743989e+00, 0.30651563e+00, 0.30559134e+00,
    0.30466709e+00, 0.30376258e+00, 0.30285808e+00, 0.30195358e+00,
    0.30104908e+00, 0.30014458e+00, 0.29924008e+00, 0.29833558e+00,
    0.29743108e+00, 0.29654059e+00, 0.29565009e+00, 0.29475963e+00,
    0.29386914e+00, 0.29297864e+00, 0.29208818e+00, 0.29119769e+00,
    0.29030719e+00, 0.28942221e+00, 0.28853720e+00, 0.28765219e+00,
    0.28676718e+00, 0.28588220e+00, 0.28499720e+00, 0.28411219e+00,
    0.28322718e+00, 0.28235060e+00, 0.28147399e+00, 0.28059739e+00,
    0.27972078e+00, 0.27884418e+00, 0.27796757e+00, 0.27709100e+00,
    0.27621439e+00, 0.27534911e+00, 0.27448383e+00, 0.27361855e+00,
    0.27275327e+00, 0.27188799e+00, 0.27102271e+00, 0.27015743e+00,
    0.26929215e+00, 0.26843825e+00,
    0.26758438e+00, 0.26673046e+00, 0.26587656e+00, 0.26502264e+00,
    0.26416874e+00, 0.26331481e+00, 0.26246092e+00, 0.26162422e+00,
    0.26078746e+00, 0.25995073e+00, 0.25911397e+00, 0.25827724e+00,
    0.25744051e+00, 0.25660378e+00, 0.25576705e+00, 0.25494036e+00,
    0.25411367e+00, 0.25328699e+00, 0.25246030e+00, 0.25163361e+00,
    0.25080693e+00, 0.24998026e+00, 0.24915357e+00, 0.24833694e+00,
    0.24752033e+00, 0.24670370e+00, 0.24588709e+00, 0.24507046e+00,
    0.24425384e+00, 0.24343722e+00, 0.24262060e+00, 0.24181117e+00,
    0.24100174e+00, 0.24019231e+00, 0.23938288e+00, 0.23857343e+00,
    0.23776400e+00, 0.23695457e+00, 0.23614514e+00, 0.23535152e+00,
    0.23455791e+00, 0.23376429e+00, 0.23297067e+00, 0.23217705e+00,
    0.23138343e+00, 0.23058982e+00, 0.22979620e+00, 0.22901551e+00,
    0.22823484e+00, 0.22745417e+00, 0.22667348e+00, 0.22589281e+00,
    0.22511213e+00, 0.22433145e+00, 0.22355077e+00, 0.22278591e+00,
    0.22202104e+00, 0.22125618e+00, 0.22049132e+00, 0.21972646e+00,
    0.21896160e+00, 0.21819673e+00, 0.21743187e+00, 0.21667995e+00,
    0.21592802e+00, 0.21517609e+00, 0.21442418e+00, 0.21367225e+00,
    0.21292032e+00, 0.21216840e+00, 0.21141647e+00, 0.21066456e+00,
    0.20991263e+00, 0.20916070e+00,
    0.20840877e+00, 0.20765686e+00, 0.20690493e+00, 0.20615301e+00,
    0.20540108e+00, 0.20467360e+00, 0.20394611e+00, 0.20321864e+00,
    0.20249115e+00, 0.20176367e+00, 0.20103619e+00, 0.20030871e+00,
    0.19958122e+00, 0.19886380e+00, 0.19814639e+00, 0.19742897e+00,
    0.19671154e+00, 0.19599412e+00, 0.19527671e+00, 0.19455929e+00,
    0.19384187e+00, 0.19312805e+00, 0.19241422e+00, 0.19170040e+00,
    0.19098657e+00, 0.19027275e+00, 0.18955892e+00, 0.18884510e+00,
    0.18813135e+00, 0.18744412e+00, 0.18675689e+00, 0.18606967e+00,
    0.18538244e+00, 0.18469521e+00, 0.18400799e+00, 0.18332076e+00,
    0.18263352e+00, 0.18195061e+00, 0.18126769e+00, 0.18058479e+00,
    0.17990187e+00, 0.17921896e+00, 0.17853604e+00, 0.17785314e+00,
    0.17717022e+00, 0.17649521e+00, 0.17582020e+00, 0.17514519e+00,
    0.17447019e+00, 0.17379518e+00, 0.17312019e+00, 0.17244518e+00,
    0.17177017e+00, 0.17111744e+00, 0.17046472e+00, 0.16981199e+00,
    0.16915928e+00, 0.16850656e+00, 0.16785383e+00, 0.16720112e+00,
    0.16654839e+00, 0.16590358e+00, 0.16525877e+00, 0.16461395e+00,
    0.16396913e+00, 0.16332433e+00, 0.16267951e+00, 0.16203469e+00,
    0.16138987e+00, 0.16075872e+00, 0.16012757e+00, 0.15949641e+00,
    0.15886526e+00, 0.15823410e+00,
    0.15760294e+00, 0.15697178e+00, 0.15634063e+00, 0.15572312e+00,
    0.15510564e+00, 0.15448813e+00, 0.15387064e+00, 0.15325314e+00,
    0.15263565e+00, 0.15201814e+00, 0.15140064e+00, 0.15079501e+00,
    0.15018937e+00, 0.14958374e+00, 0.14897810e+00, 0.14837246e+00,
    0.14776683e+00, 0.14716119e+00, 0.14655556e+00, 0.14595926e+00,
    0.14536297e+00, 0.14476667e+00, 0.14417039e+00, 0.14357410e+00,
    0.14297780e+00, 0.14238152e+00, 0.14178522e+00, 0.14120367e+00,
    0.14062211e+00, 0.14004056e+00, 0.13945900e+00, 0.13887745e+00,
    0.13829589e+00, 0.13771434e+00, 0.13713278e+00, 0.13656327e+00,
    0.13599375e+00, 0.13542424e+00, 0.13485472e+00, 0.13428521e+00,
    0.13371570e+00, 0.13314618e+00, 0.13257667e+00, 0.13201830e+00,
    0.13145992e+00, 0.13090156e+00, 0.13034318e+00, 0.12978481e+00,
    0.12922643e+00, 0.12866807e+00, 0.12810969e+00, 0.12756822e+00,
    0.12702674e+00, 0.12648526e+00, 0.12594378e+00, 0.12540230e+00,
    0.12486082e+00, 0.12431934e+00, 0.12377787e+00, 0.12324385e+00,
    0.12270983e+00, 0.12217581e+00, 0.12164178e+00, 0.12110776e+00,
    0.12057374e+00, 0.12003972e+00, 0.11950570e+00, 0.11898556e+00,
    0.11846542e+00, 0.11794529e+00, 0.11742515e+00, 0.11690501e+00,
    0.11638487e+00, 0.11586474e+00,
    0.11534460e+00, 0.11483830e+00, 0.11433200e+00, 0.11382570e+00,
    0.11331940e+00, 0.11281310e+00, 0.11230680e+00, 0.11180050e+00,
    0.11129420e+00, 0.11079592e+00, 0.11029766e+00, 0.10979939e+00,
    0.10930111e+00, 0.10880283e+00, 0.10830457e+00, 0.10780629e+00,
    0.10730802e+00, 0.10682349e+00, 0.10633896e+00, 0.10585442e+00,
    0.10536988e+00, 0.10488536e+00, 0.10440082e+00, 0.10391629e+00,
    0.10343175e+00, 0.10296234e+00, 0.10249293e+00, 0.10202351e+00,
    0.10155410e+00, 0.10108469e+00, 0.10061527e+00, 0.10014586e+00,
    0.99676445e-01, 0.99215299e-01, 0.98754153e-01, 0.98293006e-01,
    0.97831860e-01, 0.97370714e-01, 0.96909568e-01, 0.96448421e-01,
    0.95987275e-01, 0.95534757e-01, 0.95082238e-01, 0.94629712e-01,
    0.94177194e-01, 0.93724675e-01, 0.93272157e-01, 0.92819631e-01,
    0.92367113e-01, 0.91932207e-01, 0.91497295e-01, 0.91062389e-01,
    0.90627484e-01, 0.90192571e-01, 0.89757666e-01, 0.89322753e-01,
    0.88887848e-01, 0.88462822e-01, 0.88037804e-01, 0.87612778e-01,
    0.87187752e-01, 0.86762726e-01, 0.86337708e-01, 0.85912682e-01,
    0.85487656e-01, 0.85070901e-01, 0.84654145e-01, 0.84237382e-01,
    0.83820626e-01, 0.83403870e-01, 0.82987115e-01, 0.82570359e-01,
    0.82153603e-01, 0.81752300e-01,
    0.81350997e-01, 0.80949694e-01, 0.80548391e-01, 0.80147095e-01,
    0.79745792e-01, 0.79344489e-01, 0.78943186e-01, 0.78547187e-01,
    0.78151189e-01, 0.77755183e-01, 0.77359185e-01, 0.76963186e-01,
    0.76567188e-01, 0.76171182e-01, 0.75775184e-01, 0.75396664e-01,
    0.75018138e-01, 0.74639618e-01, 0.74261092e-01, 0.73882572e-01,
    0.73504046e-01, 0.73125526e-01, 0.72747000e-01, 0.72376929e-01,
    0.72006851e-01, 0.71636774e-01, 0.71266696e-01, 0.70896618e-01,
    0.70526548e-01, 0.70156470e-01, 0.69786392e-01, 0.69424346e-01,
    0.69062300e-01, 0.68700254e-01, 0.68338208e-01, 0.67976162e-01,
    0.67614123e-01, 0.67252077e-01, 0.66890031e-01, 0.66543505e-01,
    0.66196993e-01, 0.65850474e-01, 0.65503947e-01, 0.65157436e-01,
    0.64810917e-01, 0.64464390e-01, 0.64117879e-01, 0.63776962e-01,
    0.63436046e-01, 0.63095130e-01, 0.62754214e-01, 0.62413294e-01,
    0.62072378e-01, 0.61731458e-01, 0.61390541e-01, 0.61062202e-01,
    0.60733866e-01, 0.60405526e-01, 0.60077190e-01, 0.59748851e-01,
    0.59420515e-01, 0.59092175e-01, 0.58763839e-01, 0.58445923e-01,
    0.58128010e-01, 0.57810094e-01, 0.57492182e-01, 0.57174265e-01,
    0.56856353e-01, 0.56538437e-01, 0.56220524e-01, 0.55911507e-01,
    0.55602487e-01, 0.55293467e-01,
    0.54984450e-01, 0.54675430e-01, 0.54366414e-01, 0.54057393e-01,
    0.53748377e-01, 0.53450007e-01, 0.53151637e-01, 0.52853264e-01,
    0.52554894e-01, 0.52256525e-01, 0.51958155e-01, 0.51659785e-01,
    0.51361412e-01, 0.51071849e-01, 0.50782286e-01, 0.50492719e-01,
    0.50203156e-01, 0.49913593e-01, 0.49624026e-01, 0.49334463e-01,
    0.49044900e-01, 0.48765160e-01, 0.48485424e-01, 0.48205689e-01,
    0.47925953e-01, 0.47646217e-01, 0.47366481e-01, 0.47086746e-01,
    0.46807006e-01, 0.46536516e-01, 0.46266023e-01, 0.45995526e-01,
    0.45725033e-01, 0.45454539e-01, 0.45184050e-01, 0.44913549e-01,
    0.44643059e-01, 0.44382293e-01, 0.44121530e-01, 0.43860763e-01,
    0.43599997e-01, 0.43339234e-01, 0.43078467e-01, 0.42817701e-01,
    0.42556938e-01, 0.42305876e-01, 0.42054817e-01, 0.41803755e-01,
    0.41552693e-01, 0.41301634e-01, 0.41050572e-01, 0.40799513e-01,
    0.40548451e-01, 0.40303770e-01, 0.40059090e-01, 0.39814409e-01,
    0.39569728e-01, 0.39325047e-01, 0.39080366e-01, 0.38835686e-01,
    0.38591005e-01, 0.38356613e-01, 0.38122222e-01, 0.37887830e-01,
    0.37653435e-01, 0.37419043e-01, 0.37184652e-01, 0.36950260e-01,
    0.36715869e-01, 0.36489923e-01, 0.36263976e-01, 0.36038030e-01,
    0.35812087e-01, 0.35586141e-01,
    0.35360195e-01, 0.35134248e-01, 0.34908302e-01, 0.34690928e-01,
    0.34473553e-01, 0.34256175e-01, 0.34038801e-01, 0.33821426e-01,
    0.33604052e-01, 0.33386674e-01, 0.33169299e-01, 0.32958113e-01,
    0.32746922e-01, 0.32535736e-01, 0.32324549e-01, 0.32113362e-01,
    0.31902168e-01, 0.31690981e-01, 0.31479795e-01, 0.31278696e-01,
    0.31077595e-01, 0.30876495e-01, 0.30675394e-01, 0.30474294e-01,
    0.30273195e-01, 0.30072095e-01, 0.29870994e-01, 0.29676992e-01,
    0.29482992e-01, 0.29288990e-01, 0.29094988e-01, 0.28900987e-01,
    0.28706986e-01, 0.28512985e-01, 0.28318983e-01, 0.28132036e-01,
    0.27945088e-01, 0.27758140e-01, 0.27571192e-01, 0.27384246e-01,
    0.27197298e-01, 0.27010350e-01, 0.26823401e-01, 0.26644541e-01,
    0.26465680e-01, 0.26286820e-01, 0.26107959e-01, 0.25929099e-01,
    0.25750238e-01, 0.25571378e-01, 0.25392517e-01, 0.25220431e-01,
    0.25048343e-01, 0.24876256e-01, 0.24704168e-01, 0.24532080e-01,
    0.24359994e-01, 0.24187906e-01, 0.24015818e-01, 0.23851352e-01,
    0.23686886e-01, 0.23522416e-01, 0.23357952e-01, 0.23193488e-01,
    0.23029024e-01, 0.22864562e-01, 0.22700097e-01, 0.22541769e-01,
    0.22383440e-01, 0.22225112e-01, 0.22066785e-01, 0.21908456e-01,
    0.21750128e-01, 0.21591799e-01,
    0.21433473e-01, 0.21282243e-01, 0.21131013e-01, 0.20979784e-01,
    0.20828554e-01, 0.20677324e-01, 0.20526094e-01, 0.20374866e-01,
    0.20223636e-01, 0.20079101e-01, 0.19934567e-01, 0.19790031e-01,
    0.19645495e-01, 0.19500962e-01, 0.19356426e-01, 0.19211890e-01,
    0.19067356e-01, 0.18928211e-01, 0.18789068e-01, 0.18649925e-01,
    0.18510781e-01, 0.18371638e-01, 0.18232493e-01, 0.18093349e-01,
    0.17954206e-01, 0.17822554e-01, 0.17690903e-01, 0.17559251e-01,
    0.17427599e-01, 0.17295947e-01, 0.17164296e-01, 0.17032644e-01,
    0.16900992e-01, 0.16774088e-01, 0.16647184e-01, 0.16520280e-01,
    0.16393378e-01, 0.16266475e-01, 0.16139571e-01, 0.16012667e-01,
    0.15885765e-01, 0.15765805e-01, 0.15645845e-01, 0.15525886e-01,
    0.15405927e-01, 0.15285968e-01, 0.15166009e-01, 0.15046050e-01,
    0.14926090e-01, 0.14810759e-01, 0.14695427e-01, 0.14580096e-01,
    0.14464764e-01, 0.14349433e-01, 0.14234101e-01, 0.14118769e-01,
    0.14003438e-01, 0.13893498e-01, 0.13783557e-01, 0.13673618e-01,
    0.13563678e-01, 0.13453737e-01, 0.13343797e-01, 0.13233857e-01,
    0.13123917e-01, 0.13020570e-01, 0.12917223e-01, 0.12813876e-01,
    0.12710529e-01, 0.12607182e-01, 0.12503835e-01, 0.12400489e-01,
    0.12297142e-01, 0.12197331e-01,
    0.12097519e-01, 0.11997707e-01, 0.11897895e-01, 0.11798085e-01,
    0.11698273e-01, 0.11598461e-01, 0.11498650e-01, 0.11404547e-01,
    0.11310444e-01, 0.11216342e-01, 0.11122239e-01, 0.11028135e-01,
    0.10934032e-01, 0.10839930e-01, 0.10745827e-01, 0.10656666e-01,
    0.10567506e-01, 0.10478345e-01, 0.10389185e-01, 0.10300023e-01,
    0.10210863e-01, 0.10121702e-01, 0.10032542e-01, 0.99474415e-02,
    0.98623410e-02, 0.97772414e-02, 0.96921409e-02, 0.96070403e-02,
    0.95219398e-02, 0.94368402e-02, 0.93517397e-02, 0.92714624e-02,
    0.91911843e-02, 0.91109071e-02, 0.90306299e-02, 0.89503527e-02,
    0.88700745e-02, 0.87897973e-02, 0.87095201e-02, 0.86332085e-02,
    0.85568978e-02, 0.84805870e-02, 0.84042754e-02, 0.83279647e-02,
    0.82516531e-02, 0.81753423e-02, 0.80990307e-02, 0.80274092e-02,
    0.79557877e-02, 0.78841662e-02, 0.78125438e-02, 0.77409227e-02,
    0.76693008e-02, 0.75976793e-02, 0.75260573e-02, 0.74578207e-02,
    0.73895841e-02, 0.73213475e-02, 0.72531109e-02, 0.71848743e-02,
    0.71166377e-02, 0.70484011e-02, 0.69801644e-02, 0.69160289e-02,
    0.68518934e-02, 0.67877579e-02, 0.67236228e-02, 0.66594873e-02,
    0.65953517e-02, 0.65312162e-02, 0.64670807e-02, 0.64067137e-02,
    0.63463463e-02, 0.62859794e-02,
    0.62256120e-02, 0.61652451e-02, 0.61048781e-02, 0.60445108e-02,
    0.59841438e-02, 0.59273192e-02, 0.58704941e-02, 0.58136694e-02,
    0.57568448e-02, 0.57000201e-02, 0.56431950e-02, 0.55863704e-02,
    0.55295457e-02, 0.54760752e-02, 0.54226047e-02, 0.53691338e-02,
    0.53156633e-02, 0.52621928e-02, 0.52087223e-02, 0.51552518e-02,
    0.51017809e-02, 0.50516059e-02, 0.50014304e-02, 0.49512549e-02,
    0.49010799e-02, 0.48509045e-02, 0.48007290e-02, 0.47505535e-02,
    0.47003785e-02, 0.46531809e-02, 0.46059834e-02, 0.45587863e-02,
    0.45115887e-02, 0.44643912e-02, 0.44171936e-02, 0.43699965e-02,
    0.43227989e-02, 0.42786733e-02, 0.42345482e-02, 0.41904226e-02,
    0.41462970e-02, 0.41021719e-02, 0.40580463e-02, 0.40139207e-02,
    0.39697955e-02, 0.39285161e-02, 0.38872364e-02, 0.38459569e-02,
    0.38046774e-02, 0.37633979e-02, 0.37221184e-02, 0.36808390e-02,
    0.36395595e-02, 0.36008984e-02, 0.35622374e-02, 0.35235765e-02,
    0.34849155e-02, 0.34462544e-02, 0.34075934e-02, 0.33689325e-02,
    0.33302715e-02, 0.32942332e-02, 0.32581948e-02, 0.32221565e-02,
    0.31861183e-02, 0.31500799e-02, 0.31140416e-02, 0.30780034e-02,
    0.30419650e-02, 0.30085032e-02, 0.29750413e-02, 0.29415793e-02,
    0.29081174e-02, 0.28746554e-02,
    0.28411935e-02, 0.28077315e-02, 0.27742696e-02, 0.27431117e-02,
    0.27119538e-02, 0.26807957e-02, 0.26496379e-02, 0.26184800e-02,
    0.25873219e-02, 0.25561641e-02, 0.25250062e-02, 0.24961114e-02,
    0.24672167e-02, 0.24383222e-02, 0.24094274e-02, 0.23805327e-02,
    0.23516382e-02, 0.23227434e-02, 0.22938489e-02, 0.22669465e-02,
    0.22400441e-02, 0.22131416e-02, 0.21862392e-02, 0.21593370e-02,
    0.21324346e-02, 0.21055322e-02, 0.20786298e-02, 0.20539192e-02,
    0.20292085e-02, 0.20044977e-02, 0.19797871e-02, 0.19550764e-02,
    0.19303657e-02, 0.19056550e-02, 0.18809444e-02, 0.18579690e-02,
    0.18349937e-02, 0.18120183e-02, 0.17890430e-02, 0.17660677e-02,
    0.17430923e-02, 0.17201171e-02, 0.16971417e-02, 0.16759790e-02,
    0.16548162e-02, 0.16336534e-02, 0.16124907e-02, 0.15913279e-02,
    0.15701653e-02, 0.15490025e-02, 0.15278397e-02, 0.15083717e-02,
    0.14889035e-02, 0.14694354e-02, 0.14499674e-02, 0.14304993e-02,
    0.14110311e-02, 0.13915631e-02, 0.13720950e-02, 0.13542295e-02,
    0.13363642e-02, 0.13184987e-02, 0.13006333e-02, 0.12827680e-02,
    0.12649025e-02, 0.12470371e-02, 0.12291716e-02, 0.12128486e-02,
    0.11965255e-02, 0.11802024e-02, 0.11638793e-02, 0.11475562e-02,
    0.11312331e-02, 0.11149100e-02,
    0.10985869e-02, 0.10835379e-02, 0.10684890e-02, 0.10534400e-02,
    0.10383911e-02, 0.10233421e-02, 0.10082932e-02, 0.99324423e-03,
    0.97819534e-03, 0.96458686e-03, 0.95097843e-03, 0.93737000e-03,
    0.92376157e-03, 0.91015315e-03, 0.89654472e-03, 0.88293623e-03,
    0.86932781e-03, 0.85691351e-03, 0.84449921e-03, 0.83208486e-03,
    0.81967056e-03, 0.80725626e-03, 0.79484191e-03, 0.78242761e-03,
    0.77001331e-03, 0.75870747e-03, 0.74740162e-03, 0.73609577e-03,
    0.72478992e-03, 0.71348407e-03, 0.70217822e-03, 0.69087237e-03,
    0.67956652e-03, 0.66933129e-03, 0.65909600e-03, 0.64886070e-03,
    0.63862541e-03, 0.62839012e-03, 0.61815488e-03, 0.60791959e-03,
    0.59768429e-03, 0.58847608e-03, 0.57926780e-03, 0.57005958e-03,
    0.56085130e-03, 0.55164308e-03, 0.54243481e-03, 0.53322659e-03,
    0.52401837e-03, 0.51571149e-03, 0.50740462e-03, 0.49909775e-03,
    0.49079087e-03, 0.48248403e-03, 0.47417716e-03, 0.46587028e-03,
    0.45756344e-03, 0.45009126e-03, 0.44261909e-03, 0.43514691e-03,
    0.42767471e-03, 0.42020253e-03, 0.41273035e-03, 0.40525818e-03,
    0.39778600e-03, 0.39110781e-03, 0.38442962e-03, 0.37775139e-03,
    0.37107320e-03, 0.36439497e-03, 0.35771678e-03, 0.35103859e-03,
    0.34436036e-03, 0.33843261e-03,
    0.33250486e-03, 0.32657708e-03, 0.32064933e-03, 0.31472158e-03,
    0.30879382e-03, 0.30286607e-03, 0.29693829e-03, 0.29165289e-03,
    0.28636746e-03, 0.28108203e-03, 0.27579663e-03, 0.27051120e-03,
    0.26522577e-03, 0.25994037e-03, 0.25465494e-03, 0.24999218e-03,
    0.24532946e-03, 0.24066672e-03, 0.23600398e-03, 0.23134124e-03,
    0.22667850e-03, 0.22201576e-03, 0.21735302e-03, 0.21324873e-03,
    0.20914443e-03, 0.20504014e-03, 0.20093584e-03, 0.19683156e-03,
    0.19272727e-03, 0.18862297e-03, 0.18451868e-03, 0.18093106e-03,
    0.17734345e-03, 0.17375585e-03, 0.17016823e-03, 0.16658062e-03,
    0.16299300e-03, 0.15940539e-03, 0.15581778e-03, 0.15267945e-03,
    0.14954111e-03, 0.14640279e-03, 0.14326446e-03, 0.14012613e-03,
    0.13698780e-03, 0.13384948e-03, 0.13071115e-03, 0.12799964e-03,
    0.12528813e-03, 0.12257663e-03, 0.11986512e-03, 0.11715361e-03,
    0.11444210e-03, 0.11173059e-03, 0.10901909e-03, 0.10667263e-03,
    0.10432617e-03, 0.10197970e-03, 0.99633238e-04, 0.97286771e-04,
    0.94940311e-04, 0.92593851e-04, 0.90247384e-04, 0.88237888e-04,
    0.86228385e-04, 0.84218889e-04, 0.82209393e-04, 0.80199890e-04,
    0.78190395e-04, 0.76180891e-04, 0.74171396e-04, 0.72455332e-04,
    0.70739276e-04, 0.69023219e-04,
    0.67307155e-04, 0.65591099e-04, 0.63875035e-04, 0.62158979e-04,
    0.60442919e-04, 0.58989761e-04, 0.57536603e-04, 0.56083449e-04,
    0.54630291e-04, 0.53177133e-04, 0.51723975e-04, 0.50270817e-04,
    0.48817659e-04, 0.47600028e-04, 0.46382393e-04, 0.45164761e-04,
    0.43947126e-04, 0.42729494e-04, 0.41511859e-04, 0.40294228e-04,
    0.39076593e-04, 0.38059032e-04, 0.37041471e-04, 0.36023910e-04,
    0.35006349e-04, 0.33988788e-04, 0.32971227e-04, 0.31953667e-04,
    0.30936106e-04, 0.30099836e-04, 0.29263569e-04, 0.28427301e-04,
    0.27591033e-04, 0.26754768e-04, 0.25918500e-04, 0.25082232e-04,
    0.24245965e-04, 0.23557821e-04, 0.22869677e-04, 0.22181533e-04,
    0.21493390e-04, 0.20805246e-04, 0.20117102e-04, 0.19428959e-04,
    0.18740815e-04, 0.18186580e-04, 0.17632343e-04, 0.17078108e-04,
    0.16523873e-04, 0.15969637e-04, 0.15415402e-04, 0.14861166e-04,
    0.14306930e-04, 0.13860979e-04, 0.13415029e-04, 0.12969078e-04,
    0.12523127e-04, 0.12077177e-04, 0.11631226e-04, 0.11185274e-04,
    0.10739323e-04, 0.10387442e-04, 0.10035560e-04, 0.96836784e-05,
    0.93317967e-05, 0.89799150e-05, 0.86280334e-05, 0.82761517e-05,
    0.79242700e-05, 0.76496099e-05, 0.73749488e-05, 0.71002878e-05,
    0.68256272e-05, 0.65509666e-05,
    0.62763056e-05, 0.60016450e-05, 0.57269840e-05, 0.55169080e-05,
    0.53068320e-05, 0.50967560e-05, 0.48866800e-05, 0.46766040e-05,
    0.44665280e-05, 0.42564520e-05, 0.40463760e-05, 0.38885996e-05,
    0.37308234e-05, 0.35730470e-05, 0.34152706e-05, 0.32574942e-05,
    0.30997178e-05, 0.29419416e-05, 0.27841652e-05, 0.26682460e-05,
    0.25523268e-05, 0.24364078e-05, 0.23204886e-05, 0.22045695e-05,
    0.20886503e-05, 0.19727313e-05, 0.18568121e-05, 0.17741506e-05,
    0.16914892e-05, 0.16088277e-05, 0.15261661e-05, 0.14435046e-05,
    0.13608432e-05, 0.12781817e-05, 0.11955202e-05, 0.11383943e-05,
    0.10812684e-05, 0.10241424e-05, 0.96701649e-06, 0.90989056e-06,
    0.85276463e-06, 0.79563870e-06, 0.73851277e-06, 0.70047270e-06,
    0.66243263e-06, 0.62439256e-06, 0.58635248e-06, 0.54831241e-06,
    0.51027234e-06, 0.47223224e-06, 0.43419217e-06, 0.40997287e-06,
    0.38575359e-06, 0.36153429e-06, 0.33731502e-06, 0.31309571e-06,
    0.28887644e-06, 0.26465716e-06, 0.24043786e-06, 0.22578344e-06,
    0.21112901e-06, 0.19647459e-06, 0.18182017e-06, 0.16716574e-06,
    0.15251132e-06, 0.13785689e-06, 0.12320247e-06, 0.11493193e-06,
    0.10666140e-06, 0.98390856e-07, 0.90120324e-07, 0.81849784e-07,
    0.73579251e-07, 0.65308711e-07,
    0.57038179e-07, 0.52771284e-07, 0.48504390e-07, 0.44237495e-07,
    0.39970601e-07, 0.35703707e-07, 0.31436812e-07, 0.27169916e-07,
    0.22903022e-07, 0.20977984e-07, 0.19052946e-07, 0.17127906e-07,
    0.15202868e-07, 0.13277830e-07, 0.11352792e-07, 0.94277528e-08,
    0.75027140e-08, 0.67869950e-08, 0.60712755e-08, 0.53555560e-08,
    0.46398365e-08, 0.39241170e-08, 0.32083975e-08, 0.24926781e-08,
    0.17769587e-08, 0.15826038e-08, 0.13882490e-08, 0.11938941e-08,
    0.99953923e-09, 0.80518436e-09, 0.61082955e-09, 0.41647469e-09,
    0.22211984e-09, 0.19435485e-09, 0.16658987e-09, 0.13882490e-09,
    0.11105992e-09, 0.83294933e-10, 0.55529959e-10, 0.27764980e-10,
    0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
    0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00};

  /* 1024 FFT */
   odouble d1024[1024] = {
     0.10000000e+01, 0.99988866e+00,
     0.99977714e+00, 0.99966598e+00, 0.99955440e+00, 0.99943626e+00,
     0.99931854e+00, 0.99920028e+00, 0.99908257e+00, 0.99894392e+00,
     0.99880564e+00, 0.99866694e+00, 0.99852866e+00, 0.99834585e+00,
     0.99816316e+00, 0.99798030e+00, 0.99779761e+00, 0.99759293e+00,
     0.99738806e+00, 0.99718332e+00, 0.99697846e+00, 0.99672168e+00,
     0.99646503e+00, 0.99620819e+00, 0.99595159e+00, 0.99561983e+00,
     0.99528790e+00, 0.99495614e+00, 0.99462402e+00, 0.99424171e+00,
     0.99385989e+00, 0.99347752e+00, 0.99309546e+00, 0.99265903e+00,
     0.99222291e+00, 0.99178654e+00, 0.99135017e+00, 0.99088120e+00,
     0.99041241e+00, 0.98994327e+00, 0.98947471e+00, 0.98897654e+00,
     0.98847890e+00, 0.98798072e+00, 0.98748291e+00, 0.98693401e+00,
     0.98638469e+00, 0.98583573e+00, 0.98528647e+00, 0.98464882e+00,
     0.98401123e+00, 0.98337340e+00, 0.98273575e+00, 0.98206103e+00,
     0.98138630e+00, 0.98071158e+00, 0.98003685e+00, 0.97932982e+00,
     0.97862303e+00, 0.97791600e+00, 0.97720921e+00, 0.97643286e+00,
     0.97565633e+00, 0.97487992e+00, 0.97410339e+00, 0.97330886e+00,
     0.97251421e+00, 0.97171968e+00, 0.97092479e+00, 0.97007149e+00,
     0.96921808e+00, 0.96836478e+00, 0.96751148e+00, 0.96662027e+00,
     0.96572924e+00, 0.96483827e+00,
     0.96394724e+00, 0.96301448e+00, 0.96208197e+00, 0.96114910e+00,
     0.96021652e+00, 0.95923299e+00, 0.95824957e+00, 0.95726597e+00,
     0.95628262e+00, 0.95522398e+00, 0.95416504e+00, 0.95310640e+00,
     0.95204747e+00, 0.95097250e+00, 0.94989753e+00, 0.94882256e+00,
     0.94774759e+00, 0.94663709e+00, 0.94552654e+00, 0.94441599e+00,
     0.94330525e+00, 0.94213003e+00, 0.94095480e+00, 0.93977952e+00,
     0.93860412e+00, 0.93741786e+00, 0.93623137e+00, 0.93504506e+00,
     0.93385857e+00, 0.93262827e+00, 0.93139815e+00, 0.93016779e+00,
     0.92893767e+00, 0.92763728e+00, 0.92633682e+00, 0.92503643e+00,
     0.92373604e+00, 0.92240590e+00, 0.92107606e+00, 0.91974592e+00,
     0.91841596e+00, 0.91704565e+00, 0.91567552e+00, 0.91430515e+00,
     0.91293484e+00, 0.91155249e+00, 0.91017032e+00, 0.90878779e+00,
     0.90740561e+00, 0.90598828e+00, 0.90457118e+00, 0.90315384e+00,
     0.90173650e+00, 0.90026438e+00, 0.89879227e+00, 0.89732015e+00,
     0.89584804e+00, 0.89433223e+00, 0.89281642e+00, 0.89130062e+00,
     0.88978499e+00, 0.88823426e+00, 0.88668352e+00, 0.88513273e+00,
     0.88358217e+00, 0.88200194e+00, 0.88042152e+00, 0.87884134e+00,
     0.87726110e+00, 0.87563771e+00, 0.87401456e+00, 0.87239122e+00,
     0.87076783e+00, 0.86913085e+00,
     0.86749381e+00, 0.86585683e+00, 0.86421996e+00, 0.86255211e+00,
     0.86088419e+00, 0.85921639e+00, 0.85754848e+00, 0.85581571e+00,
     0.85408300e+00, 0.85235023e+00, 0.85061741e+00, 0.84887052e+00,
     0.84712362e+00, 0.84537667e+00, 0.84362978e+00, 0.84187919e+00,
     0.84012878e+00, 0.83837819e+00, 0.83662760e+00, 0.83483887e+00,
     0.83305001e+00, 0.83126116e+00, 0.82947224e+00, 0.82766581e+00,
     0.82585895e+00, 0.82405245e+00, 0.82224572e+00, 0.82039094e+00,
     0.81853610e+00, 0.81668133e+00, 0.81482649e+00, 0.81293780e+00,
     0.81104904e+00, 0.80916017e+00, 0.80727148e+00, 0.80536789e+00,
     0.80346406e+00, 0.80156046e+00, 0.79965663e+00, 0.79772937e+00,
     0.79580194e+00, 0.79387486e+00, 0.79194742e+00, 0.78998226e+00,
     0.78801692e+00, 0.78605175e+00, 0.78408641e+00, 0.78212947e+00,
     0.78017259e+00, 0.77821571e+00, 0.77625877e+00, 0.77427220e+00,
     0.77228558e+00, 0.77029890e+00, 0.76831234e+00, 0.76629901e+00,
     0.76428568e+00, 0.76227242e+00, 0.76025909e+00, 0.75823212e+00,
     0.75620520e+00, 0.75417829e+00, 0.75215113e+00, 0.75010371e+00,
     0.74805623e+00, 0.74600887e+00, 0.74396139e+00, 0.74189389e+00,
     0.73982632e+00, 0.73775893e+00, 0.73569137e+00, 0.73360521e+00,
     0.73151916e+00, 0.72943294e+00,
     0.72734684e+00, 0.72525400e+00, 0.72316122e+00, 0.72106838e+00,
     0.71897554e+00, 0.71686584e+00, 0.71475625e+00, 0.71264660e+00,
     0.71053690e+00, 0.70841002e+00, 0.70628315e+00, 0.70415640e+00,
     0.70202953e+00, 0.69988197e+00, 0.69773453e+00, 0.69558698e+00,
     0.69343948e+00, 0.69127214e+00, 0.68910480e+00, 0.68693745e+00,
     0.68477005e+00, 0.68261290e+00, 0.68045557e+00, 0.67829841e+00,
     0.67614108e+00, 0.67398292e+00, 0.67182481e+00, 0.66966665e+00,
     0.66750842e+00, 0.66533726e+00, 0.66316605e+00, 0.66099483e+00,
     0.65882361e+00, 0.65664607e+00, 0.65446848e+00, 0.65229088e+00,
     0.65011328e+00, 0.64791012e+00, 0.64570689e+00, 0.64350367e+00,
     0.64130050e+00, 0.63908291e+00, 0.63686532e+00, 0.63464773e+00,
     0.63243014e+00, 0.63020587e+00, 0.62798154e+00, 0.62575728e+00,
     0.62353295e+00, 0.62133086e+00, 0.61912864e+00, 0.61692649e+00,
     0.61472434e+00, 0.61252660e+00, 0.61032891e+00, 0.60813117e+00,
     0.60593343e+00, 0.60371804e+00, 0.60150260e+00, 0.59928715e+00,
     0.59707177e+00, 0.59485382e+00, 0.59263587e+00, 0.59041792e+00,
     0.58819997e+00, 0.58596778e+00, 0.58373559e+00, 0.58150333e+00,
     0.57927114e+00, 0.57704818e+00, 0.57482517e+00, 0.57260215e+00,
     0.57037914e+00, 0.56817263e+00,
     0.56596607e+00, 0.56375957e+00, 0.56155300e+00, 0.55933517e+00,
     0.55711734e+00, 0.55489945e+00, 0.55268162e+00, 0.55047923e+00,
     0.54827666e+00, 0.54607427e+00, 0.54387170e+00, 0.54166967e+00,
     0.53946745e+00, 0.53726542e+00, 0.53506321e+00, 0.53284973e+00,
     0.53063631e+00, 0.52842277e+00, 0.52620929e+00, 0.52399343e+00,
     0.52177757e+00, 0.51956171e+00, 0.51734591e+00, 0.51515377e+00,
     0.51296163e+00, 0.51076943e+00, 0.50857729e+00, 0.50638032e+00,
     0.50418329e+00, 0.50198627e+00, 0.49978930e+00, 0.49760693e+00,
     0.49542457e+00, 0.49324220e+00, 0.49105987e+00, 0.48889726e+00,
     0.48673469e+00, 0.48457208e+00, 0.48240951e+00, 0.48024017e+00,
     0.47807086e+00, 0.47590151e+00, 0.47373223e+00, 0.47159189e+00,
     0.46945158e+00, 0.46731123e+00, 0.46517092e+00, 0.46302551e+00,
     0.46088007e+00, 0.45873466e+00, 0.45658922e+00, 0.45443878e+00,
     0.45228833e+00, 0.45013785e+00, 0.44798741e+00, 0.44586658e+00,
     0.44374579e+00, 0.44162500e+00, 0.43950418e+00, 0.43741313e+00,
     0.43532205e+00, 0.43323100e+00, 0.43113995e+00, 0.42904490e+00,
     0.42694989e+00, 0.42485487e+00, 0.42275986e+00, 0.42067465e+00,
     0.41858947e+00, 0.41650426e+00, 0.41441903e+00, 0.41234154e+00,
     0.41026407e+00, 0.40818664e+00,
     0.40610915e+00, 0.40405071e+00, 0.40199226e+00, 0.39993382e+00,
     0.39787537e+00, 0.39582628e+00, 0.39377716e+00, 0.39172807e+00,
     0.38967898e+00, 0.38767374e+00, 0.38566849e+00, 0.38366324e+00,
     0.38165799e+00, 0.37964931e+00, 0.37764066e+00, 0.37563208e+00,
     0.37362340e+00, 0.37164736e+00, 0.36967131e+00, 0.36769527e+00,
     0.36571923e+00, 0.36374140e+00, 0.36176354e+00, 0.35978571e+00,
     0.35780787e+00, 0.35583410e+00, 0.35386035e+00, 0.35188660e+00,
     0.34991279e+00, 0.34798345e+00, 0.34605411e+00, 0.34412473e+00,
     0.34219536e+00, 0.34029329e+00, 0.33839118e+00, 0.33648908e+00,
     0.33458701e+00, 0.33268994e+00, 0.33079287e+00, 0.32889581e+00,
     0.32699874e+00, 0.32512325e+00, 0.32324776e+00, 0.32137227e+00,
     0.31949678e+00, 0.31763780e+00, 0.31577885e+00, 0.31391987e+00,
     0.31206092e+00, 0.31021237e+00, 0.30836385e+00, 0.30651531e+00,
     0.30466679e+00, 0.30285776e+00, 0.30104876e+00, 0.29923975e+00,
     0.29743075e+00, 0.29564980e+00, 0.29386884e+00, 0.29208788e+00,
     0.29030690e+00, 0.28853691e+00, 0.28676689e+00, 0.28499690e+00,
     0.28322691e+00, 0.28147370e+00, 0.27972049e+00, 0.27796730e+00,
     0.27621409e+00, 0.27448353e+00, 0.27275297e+00, 0.27102244e+00,
     0.26929188e+00, 0.26758409e+00,
     0.26587629e+00, 0.26416850e+00, 0.26246071e+00, 0.26078728e+00,
     0.25911379e+00, 0.25744036e+00, 0.25576693e+00, 0.25411355e+00,
     0.25246018e+00, 0.25080684e+00, 0.24915345e+00, 0.24752022e+00,
     0.24588698e+00, 0.24425374e+00, 0.24262050e+00, 0.24100164e+00,
     0.23938277e+00, 0.23776390e+00, 0.23614503e+00, 0.23455781e+00,
     0.23297057e+00, 0.23138334e+00, 0.22979610e+00, 0.22823475e+00,
     0.22667339e+00, 0.22511204e+00, 0.22355068e+00, 0.22202095e+00,
     0.22049123e+00, 0.21896151e+00, 0.21743178e+00, 0.21592793e+00,
     0.21442409e+00, 0.21292023e+00, 0.21141639e+00, 0.20991255e+00,
     0.20840870e+00, 0.20690486e+00, 0.20540100e+00, 0.20394604e+00,
     0.20249107e+00, 0.20103611e+00, 0.19958115e+00, 0.19814631e+00,
     0.19671148e+00, 0.19527663e+00, 0.19384180e+00, 0.19241415e+00,
     0.19098651e+00, 0.18955886e+00, 0.18813115e+00, 0.18675670e+00,
     0.18538225e+00, 0.18400779e+00, 0.18263334e+00, 0.18126751e+00,
     0.17990169e+00, 0.17853586e+00, 0.17717004e+00, 0.17582002e+00,
     0.17447002e+00, 0.17312001e+00, 0.17176999e+00, 0.17046455e+00,
     0.16915911e+00, 0.16785367e+00, 0.16654822e+00, 0.16525860e+00,
     0.16396897e+00, 0.16267934e+00, 0.16138971e+00, 0.16012740e+00,
     0.15886509e+00, 0.15760279e+00,
     0.15634046e+00, 0.15510547e+00, 0.15387048e+00, 0.15263548e+00,
     0.15140049e+00, 0.15018922e+00, 0.14897795e+00, 0.14776668e+00,
     0.14655541e+00, 0.14536282e+00, 0.14417024e+00, 0.14297765e+00,
     0.14178507e+00, 0.14062196e+00, 0.13945885e+00, 0.13829574e+00,
     0.13713264e+00, 0.13599361e+00, 0.13485458e+00, 0.13371556e+00,
     0.13257653e+00, 0.13145979e+00, 0.13034305e+00, 0.12922630e+00,
     0.12810956e+00, 0.12702660e+00, 0.12594365e+00, 0.12486069e+00,
     0.12377773e+00, 0.12270969e+00, 0.12164165e+00, 0.12057361e+00,
     0.11950557e+00, 0.11846530e+00, 0.11742502e+00, 0.11638474e+00,
     0.11534447e+00, 0.11433187e+00, 0.11331928e+00, 0.11230668e+00,
     0.11129408e+00, 0.11029755e+00, 0.10930101e+00, 0.10830447e+00,
     0.10730793e+00, 0.10633886e+00, 0.10536980e+00, 0.10440073e+00,
     0.10343167e+00, 0.10249285e+00, 0.10155402e+00, 0.10061520e+00,
     0.99676363e-01, 0.98754071e-01, 0.97831786e-01, 0.96909493e-01,
     0.95987201e-01, 0.95082156e-01, 0.94177097e-01, 0.93272053e-01,
     0.92367016e-01, 0.91497198e-01, 0.90627387e-01, 0.89757569e-01,
     0.88887751e-01, 0.88037707e-01, 0.87187663e-01, 0.86337611e-01,
     0.85487567e-01, 0.84654048e-01, 0.83820537e-01, 0.82987025e-01,
     0.82153514e-01, 0.81350908e-01,
     0.80548309e-01, 0.79745702e-01, 0.78943104e-01, 0.78151099e-01,
     0.77359103e-01, 0.76567098e-01, 0.75775102e-01, 0.75018056e-01,
     0.74261010e-01, 0.73503964e-01, 0.72746918e-01, 0.72006769e-01,
     0.71266614e-01, 0.70526466e-01, 0.69786310e-01, 0.69062226e-01,
     0.68338133e-01, 0.67614041e-01, 0.66889949e-01, 0.66196926e-01,
     0.65503895e-01, 0.64810857e-01, 0.64117834e-01, 0.63435994e-01,
     0.62754162e-01, 0.62072325e-01, 0.61390489e-01, 0.60733818e-01,
     0.60077142e-01, 0.59420466e-01, 0.58763791e-01, 0.58127962e-01,
     0.57492133e-01, 0.56856308e-01, 0.56220479e-01, 0.55602442e-01,
     0.54984406e-01, 0.54366369e-01, 0.53748332e-01, 0.53151593e-01,
     0.52554853e-01, 0.51958110e-01, 0.51361371e-01, 0.50782245e-01,
     0.50203115e-01, 0.49623985e-01, 0.49044859e-01, 0.48485387e-01,
     0.47925916e-01, 0.47366440e-01, 0.46806980e-01, 0.46265990e-01,
     0.45725003e-01, 0.45184013e-01, 0.44643022e-01, 0.44121493e-01,
     0.43599963e-01, 0.43078434e-01, 0.42556901e-01, 0.42054780e-01,
     0.41552659e-01, 0.41050538e-01, 0.40548418e-01, 0.40059056e-01,
     0.39569695e-01, 0.39080337e-01, 0.38590975e-01, 0.38122188e-01,
     0.37653405e-01, 0.37184622e-01, 0.36715839e-01, 0.36263946e-01,
     0.35812058e-01, 0.35360165e-01,
     0.34908276e-01, 0.34473523e-01, 0.34038775e-01, 0.33604022e-01,
     0.33169273e-01, 0.32746892e-01, 0.32324515e-01, 0.31902138e-01,
     0.31479757e-01, 0.31077558e-01, 0.30675359e-01, 0.30273158e-01,
     0.29870959e-01, 0.29482957e-01, 0.29094955e-01, 0.28706951e-01,
     0.28318949e-01, 0.27945055e-01, 0.27571158e-01, 0.27197264e-01,
     0.26823370e-01, 0.26465649e-01, 0.26107928e-01, 0.25750207e-01,
     0.25392486e-01, 0.25048312e-01, 0.24704136e-01, 0.24359962e-01,
     0.24015788e-01, 0.23686860e-01, 0.23357933e-01, 0.23029005e-01,
     0.22700079e-01, 0.22383424e-01, 0.22066766e-01, 0.21750111e-01,
     0.21433454e-01, 0.21130996e-01, 0.20828538e-01, 0.20526079e-01,
     0.20223619e-01, 0.19934550e-01, 0.19645480e-01, 0.19356411e-01,
     0.19067340e-01, 0.18789053e-01, 0.18510766e-01, 0.18232478e-01,
     0.17954191e-01, 0.17690888e-01, 0.17427584e-01, 0.17164281e-01,
     0.16900977e-01, 0.16647171e-01, 0.16393363e-01, 0.16139558e-01,
     0.15885752e-01, 0.15645834e-01, 0.15405915e-01, 0.15165997e-01,
     0.14926078e-01, 0.14695415e-01, 0.14464753e-01, 0.14234089e-01,
     0.14003427e-01, 0.13783546e-01, 0.13563666e-01, 0.13343787e-01,
     0.13123906e-01, 0.12917213e-01, 0.12710519e-01, 0.12503826e-01,
     0.12297132e-01, 0.12097510e-01,
     0.11897886e-01, 0.11698264e-01, 0.11498640e-01, 0.11310435e-01,
     0.11122230e-01, 0.10934024e-01, 0.10745819e-01, 0.10567497e-01,
     0.10389176e-01, 0.10210855e-01, 0.10032534e-01, 0.98623335e-02,
     0.96921325e-02, 0.95219323e-02, 0.93517322e-02, 0.91911769e-02,
     0.90306224e-02, 0.88700680e-02, 0.87095127e-02, 0.85568912e-02,
     0.84042689e-02, 0.82516465e-02, 0.80990242e-02, 0.79557812e-02,
     0.78125382e-02, 0.76692947e-02, 0.75260513e-02, 0.73895780e-02,
     0.72531053e-02, 0.71166321e-02, 0.69801589e-02, 0.68518878e-02,
     0.67236172e-02, 0.65953461e-02, 0.64670756e-02, 0.63463412e-02,
     0.62256074e-02, 0.61048730e-02, 0.59841392e-02, 0.58704894e-02,
     0.57568401e-02, 0.56431908e-02, 0.55295411e-02, 0.54226001e-02,
     0.53156591e-02, 0.52087181e-02, 0.51017771e-02, 0.50014262e-02,
     0.49010757e-02, 0.48007253e-02, 0.47003743e-02, 0.46059797e-02,
     0.45115850e-02, 0.44171903e-02, 0.43227957e-02, 0.42345445e-02,
     0.41462937e-02, 0.40580430e-02, 0.39697923e-02, 0.38872333e-02,
     0.38046744e-02, 0.37221154e-02, 0.36395565e-02, 0.35622346e-02,
     0.34849127e-02, 0.34075908e-02, 0.33302687e-02, 0.32581922e-02,
     0.31861158e-02, 0.31140391e-02, 0.30419626e-02, 0.29750387e-02,
     0.29081150e-02, 0.28411911e-02,
     0.27742675e-02, 0.27119515e-02, 0.26496358e-02, 0.25873198e-02,
     0.25250041e-02, 0.24672148e-02, 0.24094256e-02, 0.23516363e-02,
     0.22938470e-02, 0.22400422e-02, 0.21862376e-02, 0.21324328e-02,
     0.20786282e-02, 0.20292068e-02, 0.19797855e-02, 0.19303642e-02,
     0.18809428e-02, 0.18349922e-02, 0.17890416e-02, 0.17430909e-02,
     0.16971403e-02, 0.16548148e-02, 0.16124895e-02, 0.15701640e-02,
     0.15278385e-02, 0.14889024e-02, 0.14499662e-02, 0.14110301e-02,
     0.13720939e-02, 0.13363630e-02, 0.13006323e-02, 0.12649015e-02,
     0.12291707e-02, 0.11965245e-02, 0.11638784e-02, 0.11312322e-02,
     0.10985860e-02, 0.10684881e-02, 0.10383902e-02, 0.10082924e-02,
     0.97819453e-03, 0.95097767e-03, 0.92376082e-03, 0.89654396e-03,
     0.86932711e-03, 0.84449851e-03, 0.81966992e-03, 0.79484127e-03,
     0.77001267e-03, 0.74740103e-03, 0.72478934e-03, 0.70217764e-03,
     0.67956600e-03, 0.65909547e-03, 0.63862489e-03, 0.61815436e-03,
     0.59768383e-03, 0.57926733e-03, 0.56085089e-03, 0.54243440e-03,
     0.52401790e-03, 0.50740421e-03, 0.49079047e-03, 0.47417678e-03,
     0.45756306e-03, 0.44261871e-03, 0.42767439e-03, 0.41273003e-03,
     0.39778568e-03, 0.38442930e-03, 0.37107288e-03, 0.35771649e-03,
     0.34436010e-03, 0.33250457e-03,
     0.32064907e-03, 0.30879356e-03, 0.29693806e-03, 0.28636723e-03,
     0.27579640e-03, 0.26522556e-03, 0.25465473e-03, 0.24532925e-03,
     0.23600379e-03, 0.22667831e-03, 0.21735285e-03, 0.20914426e-03,
     0.20093568e-03, 0.19272711e-03, 0.18451853e-03, 0.17734332e-03,
     0.17016809e-03, 0.16299287e-03, 0.15581764e-03, 0.14954100e-03,
     0.14326435e-03, 0.13698770e-03, 0.13071104e-03, 0.12528803e-03,
     0.11986502e-03, 0.11444201e-03, 0.10901900e-03, 0.10432608e-03,
     0.99633158e-04, 0.94940231e-04, 0.90247311e-04, 0.86228320e-04,
     0.82209321e-04, 0.78190329e-04, 0.74171337e-04, 0.70739217e-04,
     0.67307104e-04, 0.63874984e-04, 0.60442871e-04, 0.57536559e-04,
     0.54630247e-04, 0.51723935e-04, 0.48817623e-04, 0.46382356e-04,
     0.43947090e-04, 0.41511827e-04, 0.39076560e-04, 0.37041442e-04,
     0.35006320e-04, 0.32971198e-04, 0.30936080e-04, 0.29263545e-04,
     0.27591012e-04, 0.25918478e-04, 0.24245945e-04, 0.22869659e-04,
     0.21493372e-04, 0.20117086e-04, 0.18740800e-04, 0.17632330e-04,
     0.16523860e-04, 0.15415389e-04, 0.14306918e-04, 0.13415018e-04,
     0.12523117e-04, 0.11631216e-04, 0.10739315e-04, 0.10035552e-04,
     0.93317894e-05, 0.86280270e-05, 0.79242636e-05, 0.73749429e-05,
     0.68256218e-05, 0.62763006e-05,
     0.57269795e-05, 0.53068279e-05, 0.48866764e-05, 0.44665244e-05,
     0.40463729e-05, 0.37308203e-05, 0.34152679e-05, 0.30997153e-05,
     0.27841629e-05, 0.25523248e-05, 0.23204868e-05, 0.20886487e-05,
     0.18568106e-05, 0.16914878e-05, 0.15261650e-05, 0.13608421e-05,
     0.11955193e-05, 0.10812674e-05, 0.96701569e-06, 0.85276395e-06,
     0.73851220e-06, 0.66243211e-06, 0.58635197e-06, 0.51027189e-06,
     0.43419180e-06, 0.38575328e-06, 0.33731473e-06, 0.28887621e-06,
     0.24043766e-06, 0.21112884e-06, 0.18182001e-06, 0.15251119e-06,
     0.12320237e-06, 0.10666131e-06, 0.90120253e-07, 0.73579187e-07,
     0.57038132e-07, 0.48504351e-07, 0.39970569e-07, 0.31436787e-07,
     0.22903004e-07, 0.19052930e-07, 0.15202856e-07, 0.11352782e-07,
     0.75027078e-08, 0.60712706e-08, 0.46398325e-08, 0.32083949e-08,
     0.17769572e-08, 0.13882479e-08, 0.99953845e-09, 0.61082905e-09,
     0.22211966e-09, 0.16658974e-09, 0.11105983e-09, 0.55529914e-10,
     0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00};
   
   /* 512  FFT */
   odouble d512[512] = {
     0.10000000e+01, 0.99977738e+00,
     0.99955350e+00, 0.99932396e+00, 0.99907887e+00, 0.99882150e+00,
     0.99851823e+00, 0.99815977e+00, 0.99778670e+00, 0.99737889e+00,
     0.99696177e+00, 0.99653655e+00, 0.99598241e+00, 0.99532211e+00,
     0.99462718e+00, 0.99386555e+00, 0.99309671e+00, 0.99224967e+00,
     0.99135005e+00, 0.99043524e+00, 0.98946905e+00, 0.98848206e+00,
     0.98748195e+00, 0.98643523e+00, 0.98527324e+00, 0.98400444e+00,
     0.98271757e+00, 0.98137426e+00, 0.98000455e+00, 0.97859663e+00,
     0.97717601e+00, 0.97565168e+00, 0.97408158e+00, 0.97250205e+00,
     0.97089565e+00, 0.96926302e+00, 0.96752203e+00, 0.96574515e+00,
     0.96395171e+00, 0.96209562e+00, 0.96016759e+00, 0.95822763e+00,
     0.95622700e+00, 0.95413423e+00, 0.95200968e+00, 0.94986582e+00,
     0.94770950e+00, 0.94550860e+00, 0.94325340e+00, 0.94095528e+00,
     0.93859607e+00, 0.93622458e+00, 0.93384284e+00, 0.93140858e+00,
     0.92888814e+00, 0.92629153e+00, 0.92366791e+00, 0.92103875e+00,
     0.91836131e+00, 0.91562372e+00, 0.91287416e+00, 0.91010678e+00,
     0.90732265e+00, 0.90450472e+00, 0.90166128e+00, 0.89878333e+00,
     0.89577246e+00, 0.89273822e+00, 0.88969934e+00, 0.88660657e+00,
     0.88349783e+00, 0.88033819e+00, 0.87716472e+00, 0.87396377e+00,
     0.87070507e+00, 0.86743820e+00,
     0.86415416e+00, 0.86081731e+00, 0.85742521e+00, 0.85399401e+00,
     0.85051143e+00, 0.84702367e+00, 0.84351701e+00, 0.84002274e+00,
     0.83650905e+00, 0.83292890e+00, 0.82933760e+00, 0.82574326e+00,
     0.82211334e+00, 0.81842262e+00, 0.81470394e+00, 0.81093401e+00,
     0.80714697e+00, 0.80334693e+00, 0.79953015e+00, 0.79570568e+00,
     0.79183626e+00, 0.78789949e+00, 0.78395820e+00, 0.78005052e+00,
     0.77612627e+00, 0.77216434e+00, 0.76815683e+00, 0.76414067e+00,
     0.76010305e+00, 0.75604302e+00, 0.75200391e+00, 0.74792558e+00,
     0.74379951e+00, 0.73966086e+00, 0.73551631e+00, 0.73134869e+00,
     0.72716713e+00, 0.72298020e+00, 0.71878481e+00, 0.71456712e+00,
     0.71033853e+00, 0.70610362e+00, 0.70184058e+00, 0.69754094e+00,
     0.69323736e+00, 0.68891102e+00, 0.68461311e+00, 0.68030405e+00,
     0.67595798e+00, 0.67165273e+00, 0.66732794e+00, 0.66298646e+00,
     0.65863591e+00, 0.65427738e+00, 0.64989108e+00, 0.64548045e+00,
     0.64106584e+00, 0.63663942e+00, 0.63219613e+00, 0.62774175e+00,
     0.62330812e+00, 0.61890966e+00, 0.61452055e+00, 0.61012667e+00,
     0.60572344e+00, 0.60130686e+00, 0.59686816e+00, 0.59242833e+00,
     0.58797306e+00, 0.58350313e+00, 0.57903087e+00, 0.57453597e+00,
     0.57012820e+00, 0.56570947e+00,
     0.56127709e+00, 0.55683422e+00, 0.55243671e+00, 0.54802883e+00,
     0.54360449e+00, 0.53917915e+00, 0.53476709e+00, 0.53034174e+00,
     0.52590686e+00, 0.52146965e+00, 0.51707631e+00, 0.51268637e+00,
     0.50829428e+00, 0.50389755e+00, 0.49954182e+00, 0.49516991e+00,
     0.49082050e+00, 0.48648944e+00, 0.48215082e+00, 0.47780836e+00,
     0.47345623e+00, 0.46914870e+00, 0.46485463e+00, 0.46055979e+00,
     0.45625556e+00, 0.45194918e+00, 0.44768095e+00, 0.44341093e+00,
     0.43920791e+00, 0.43502033e+00, 0.43083081e+00, 0.42663369e+00,
     0.42243624e+00, 0.41826010e+00, 0.41409963e+00, 0.40993819e+00,
     0.40577307e+00, 0.40164965e+00, 0.39752257e+00, 0.39341784e+00,
     0.38940167e+00, 0.38538477e+00, 0.38136429e+00, 0.37734050e+00,
     0.37331605e+00, 0.36935750e+00, 0.36539841e+00, 0.36143604e+00,
     0.35747191e+00, 0.35350606e+00, 0.34957314e+00, 0.34570768e+00,
     0.34188676e+00, 0.33807591e+00, 0.33426353e+00, 0.33045116e+00,
     0.32669497e+00, 0.32293731e+00, 0.31917894e+00, 0.31545442e+00,
     0.31172848e+00, 0.30802485e+00, 0.30434352e+00, 0.30071911e+00,
     0.29714006e+00, 0.29356030e+00, 0.28999168e+00, 0.28644541e+00,
     0.28292179e+00, 0.27939767e+00, 0.27590773e+00, 0.27244046e+00,
     0.26899600e+00, 0.26556277e+00,
     0.26214099e+00, 0.25878814e+00, 0.25546977e+00, 0.25215718e+00,
     0.24886762e+00, 0.24559534e+00, 0.24232307e+00, 0.23907959e+00,
     0.23587069e+00, 0.23269060e+00, 0.22953354e+00, 0.22638226e+00,
     0.22328857e+00, 0.22021793e+00, 0.21717609e+00, 0.21415731e+00,
     0.21114428e+00, 0.20813125e+00, 0.20515279e+00, 0.20221466e+00,
     0.19932261e+00, 0.19644496e+00, 0.19357020e+00, 0.19070984e+00,
     0.18790710e+00, 0.18515332e+00, 0.18240531e+00, 0.17966017e+00,
     0.17694095e+00, 0.17423326e+00, 0.17155725e+00, 0.16893886e+00,
     0.16634062e+00, 0.16375536e+00, 0.16117728e+00, 0.15862370e+00,
     0.15614069e+00, 0.15366487e+00, 0.15120779e+00, 0.14878096e+00,
     0.14635411e+00, 0.14394023e+00, 0.14159693e+00, 0.13925435e+00,
     0.13695568e+00, 0.13467287e+00, 0.13239941e+00, 0.13016124e+00,
     0.12795548e+00, 0.12577349e+00, 0.12362966e+00, 0.12148907e+00,
     0.11934920e+00, 0.11726442e+00, 0.11519171e+00, 0.11316255e+00,
     0.11115249e+00, 0.10915565e+00, 0.10717201e+00, 0.10522455e+00,
     0.10330458e+00, 0.10142358e+00, 0.99564202e-01, 0.97704820e-01,
     0.95867045e-01, 0.94042234e-01, 0.92267111e-01, 0.90524398e-01,
     0.88788882e-01, 0.87074250e-01, 0.85377619e-01, 0.83696112e-01,
     0.82062855e-01, 0.80454804e-01,
     0.78857191e-01, 0.77270389e-01, 0.75694017e-01, 0.74171484e-01,
     0.72670907e-01, 0.71176454e-01, 0.69709644e-01, 0.68258889e-01,
     0.66820569e-01, 0.65420516e-01, 0.64045049e-01, 0.62678955e-01,
     0.61335195e-01, 0.59996471e-01, 0.58703117e-01, 0.57429202e-01,
     0.56165732e-01, 0.54921709e-01, 0.53693883e-01, 0.52486766e-01,
     0.51310252e-01, 0.50144181e-01, 0.48996381e-01, 0.47875453e-01,
     0.46760462e-01, 0.45670804e-01, 0.44603549e-01, 0.43552876e-01,
     0.42518768e-01, 0.41501224e-01, 0.40502764e-01, 0.39522305e-01,
     0.38552288e-01, 0.37607297e-01, 0.36683548e-01, 0.35766643e-01,
     0.34870170e-01, 0.33990484e-01, 0.33133302e-01, 0.32281287e-01,
     0.31451672e-01, 0.30641526e-01, 0.29840380e-01, 0.29052917e-01,
     0.28285978e-01, 0.27533980e-01, 0.26794586e-01, 0.26073555e-01,
     0.25361435e-01, 0.24668986e-01, 0.23987517e-01, 0.23322731e-01,
     0.22671707e-01, 0.22031510e-01, 0.21406077e-01, 0.20797927e-01,
     0.20196617e-01, 0.19610250e-01, 0.19041527e-01, 0.18478205e-01,
     0.17930275e-01, 0.17399134e-01, 0.16875373e-01, 0.16362539e-01,
     0.15863463e-01, 0.15375935e-01, 0.14902448e-01, 0.14435084e-01,
     0.13981581e-01, 0.13538880e-01, 0.13103741e-01, 0.12683995e-01,
     0.12277076e-01, 0.11873038e-01,
     0.11479925e-01, 0.11100391e-01, 0.10728353e-01, 0.10367298e-01,
     0.10015244e-01, 0.96712457e-02, 0.93379831e-02, 0.90094004e-02,
     0.86947428e-02, 0.83885845e-02, 0.80881547e-02, 0.77957590e-02,
     0.75140754e-02, 0.72377645e-02, 0.69703134e-02, 0.67111561e-02,
     0.64578950e-02, 0.62153684e-02, 0.59753624e-02, 0.57456801e-02,
     0.55223368e-02, 0.53054420e-02, 0.50947815e-02, 0.48912321e-02,
     0.46944954e-02, 0.45026578e-02, 0.43176743e-02, 0.41382266e-02,
     0.39659129e-02, 0.37971998e-02, 0.36355697e-02, 0.34788386e-02,
     0.33274761e-02, 0.31816270e-02, 0.30396488e-02, 0.29043036e-02,
     0.27722369e-02, 0.26457638e-02, 0.25232066e-02, 0.24054425e-02,
     0.22920896e-02, 0.21828490e-02, 0.20776072e-02, 0.19760688e-02,
     0.18797063e-02, 0.17851666e-02, 0.16961629e-02, 0.16096516e-02,
     0.15267788e-02, 0.14474519e-02, 0.13711518e-02, 0.12984973e-02,
     0.12281467e-02, 0.11613821e-02, 0.10973691e-02, 0.10366503e-02,
     0.97744691e-03, 0.92180556e-03, 0.86860306e-03, 0.81750122e-03,
     0.76933118e-03, 0.72281243e-03, 0.67903567e-03, 0.63703174e-03,
     0.59723324e-03, 0.55942358e-03, 0.52353449e-03, 0.48940704e-03,
     0.45711084e-03, 0.42661146e-03, 0.39752494e-03, 0.37028373e-03,
     0.34433370e-03, 0.31997793e-03,
     0.29699001e-03, 0.27523981e-03, 0.25485954e-03, 0.23571911e-03,
     0.21767645e-03, 0.20069008e-03, 0.18489467e-03, 0.17002369e-03,
     0.15618613e-03, 0.14320549e-03, 0.13110075e-03, 0.11985292e-03,
     0.10935723e-03, 0.99634053e-04, 0.90578615e-04, 0.82212719e-04,
     0.74453026e-04, 0.67321693e-04, 0.60704082e-04, 0.54652945e-04,
     0.49071929e-04, 0.43953656e-04, 0.39292143e-04, 0.35015662e-04,
     0.31128609e-04, 0.27584741e-04, 0.24390392e-04, 0.21472595e-04,
     0.18849287e-04, 0.16494092e-04, 0.14379236e-04, 0.12484496e-04,
     0.10795986e-04, 0.93021017e-05, 0.79736565e-05, 0.68043232e-05,
     0.57764310e-05, 0.48779375e-05, 0.40919645e-05, 0.34128861e-05,
     0.28252305e-05, 0.23219652e-05, 0.18921897e-05, 0.15282563e-05,
     0.12233083e-05, 0.96846713e-06, 0.75902960e-06, 0.58668871e-06,
     0.44788416e-06, 0.33619881e-06, 0.24807244e-06, 0.17928552e-06,
     0.12667343e-06, 0.87159400e-07, 0.58128244e-07, 0.37426236e-07,
     0.23119437e-07, 0.13592559e-07, 0.75270030e-08, 0.38569024e-08,
     0.17856030e-08, 0.71698825e-09, 0.23487548e-09, 0.55628401e-10,
     0.69535501e-11, 0.00000000e+00};
   
   
   /* 256 FFT */
   odouble d256[256] = {
     0.10000000e+01, 0.99955386e+00,
     0.99908137e+00, 0.99852824e+00, 0.99779707e+00, 0.99697769e+00,
     0.99595129e+00, 0.99462312e+00, 0.99309504e+00, 0.99134922e+00,
     0.98947334e+00, 0.98748273e+00, 0.98528582e+00, 0.98273510e+00,
     0.98003578e+00, 0.97720850e+00, 0.97410303e+00, 0.97092408e+00,
     0.96751106e+00, 0.96394664e+00, 0.96021539e+00, 0.95628262e+00,
     0.95204711e+00, 0.94774687e+00, 0.94330484e+00, 0.93860388e+00,
     0.93385804e+00, 0.92893726e+00, 0.92373508e+00, 0.91841501e+00,
     0.91293347e+00, 0.90740514e+00, 0.90173596e+00, 0.89584702e+00,
     0.88978380e+00, 0.88358128e+00, 0.87726057e+00, 0.87076688e+00,
     0.86421919e+00, 0.85754800e+00, 0.85061705e+00, 0.84362918e+00,
     0.83662719e+00, 0.82947183e+00, 0.82224536e+00, 0.81482625e+00,
     0.80727124e+00, 0.79965651e+00, 0.79194707e+00, 0.78408581e+00,
     0.77625829e+00, 0.76831162e+00, 0.76025850e+00, 0.75215048e+00,
     0.74396092e+00, 0.73569107e+00, 0.72734642e+00, 0.71897489e+00,
     0.71053630e+00, 0.70202893e+00, 0.69343877e+00, 0.68476939e+00,
     0.67614055e+00, 0.66750801e+00, 0.65882319e+00, 0.65011293e+00,
     0.64130014e+00, 0.63242978e+00, 0.62353259e+00, 0.61472392e+00,
     0.60593307e+00, 0.59707135e+00, 0.58819962e+00, 0.57927078e+00,
     0.57037872e+00, 0.56155264e+00,
     0.55268121e+00, 0.54387128e+00, 0.53506279e+00, 0.52620882e+00,
     0.51734555e+00, 0.50857693e+00, 0.49978894e+00, 0.49105951e+00,
     0.48240915e+00, 0.47373194e+00, 0.46517056e+00, 0.45658886e+00,
     0.44798702e+00, 0.43950382e+00, 0.43113959e+00, 0.42275953e+00,
     0.41441873e+00, 0.40610895e+00, 0.39787516e+00, 0.38967878e+00,
     0.38165778e+00, 0.37362313e+00, 0.36571896e+00, 0.35780761e+00,
     0.34991249e+00, 0.34219509e+00, 0.33458674e+00, 0.32699850e+00,
     0.31949651e+00, 0.31206068e+00, 0.30466652e+00, 0.29743055e+00,
     0.29030669e+00, 0.28322667e+00, 0.27621388e+00, 0.26929164e+00,
     0.26246047e+00, 0.25576669e+00, 0.24915323e+00, 0.24262027e+00,
     0.23614483e+00, 0.22979589e+00, 0.22355048e+00, 0.21743158e+00,
     0.21141620e+00, 0.20540081e+00, 0.19958097e+00, 0.19384162e+00,
     0.18813100e+00, 0.18263321e+00, 0.17716990e+00, 0.17176986e+00,
     0.16654809e+00, 0.16138959e+00, 0.15634035e+00, 0.15140037e+00,
     0.14655529e+00, 0.14178497e+00, 0.13713253e+00, 0.13257642e+00,
     0.12810946e+00, 0.12377763e+00, 0.11950547e+00, 0.11534438e+00,
     0.11129399e+00, 0.10730784e+00, 0.10343158e+00, 0.99676274e-01,
     0.95987104e-01, 0.92366949e-01, 0.88887691e-01, 0.85487500e-01,
     0.82153454e-01, 0.78943044e-01,
     0.75775050e-01, 0.72746865e-01, 0.69786265e-01, 0.66889904e-01,
     0.64117797e-01, 0.61390460e-01, 0.58763761e-01, 0.56220450e-01,
     0.53748306e-01, 0.51361345e-01, 0.49044833e-01, 0.46806946e-01,
     0.44642989e-01, 0.42556871e-01, 0.40548388e-01, 0.38590945e-01,
     0.36715813e-01, 0.34908250e-01, 0.33169247e-01, 0.31479735e-01,
     0.29870939e-01, 0.28318929e-01, 0.26823349e-01, 0.25392469e-01,
     0.24015769e-01, 0.22700062e-01, 0.21433439e-01, 0.20223605e-01,
     0.19067327e-01, 0.17954178e-01, 0.16900966e-01, 0.15885739e-01,
     0.14926068e-01, 0.14003417e-01, 0.13123897e-01, 0.12297124e-01,
     0.11498632e-01, 0.10745810e-01, 0.10032526e-01, 0.93517257e-02,
     0.87095071e-02, 0.80990186e-02, 0.75260461e-02, 0.69801537e-02,
     0.64670709e-02, 0.59841345e-02, 0.55295373e-02, 0.51017734e-02,
     0.47003711e-02, 0.43227924e-02, 0.39697895e-02, 0.36395539e-02,
     0.33302663e-02, 0.30419603e-02, 0.27742654e-02, 0.25250022e-02,
     0.22938452e-02, 0.20786268e-02, 0.18809414e-02, 0.16971391e-02,
     0.15278374e-02, 0.13720929e-02, 0.12291698e-02, 0.10985852e-02,
     0.97819383e-03, 0.86932653e-03, 0.77001215e-03, 0.67956554e-03,
     0.59768336e-03, 0.52401755e-03, 0.45756274e-03, 0.39778539e-03,
     0.34435984e-03, 0.29693785e-03,
     0.25465456e-03, 0.21735269e-03, 0.18451840e-03, 0.15581754e-03,
     0.13071095e-03, 0.10901893e-03, 0.90247246e-04, 0.74171279e-04,
     0.60442828e-04, 0.48817587e-04, 0.39076534e-04, 0.30936058e-04,
     0.24245928e-04, 0.18740786e-04, 0.14306908e-04, 0.10739307e-04,
     0.79242582e-05, 0.57269754e-05, 0.40463701e-05, 0.27841609e-05,
     0.18568093e-05, 0.11955184e-05, 0.73851163e-06, 0.43419149e-06,
     0.24043749e-06, 0.12320228e-06, 0.57038090e-07, 0.22902988e-07,
     0.75027025e-08, 0.17769559e-08, 0.22211949e-09, 0.00000000e+00};
   
   /* 128 FFT */
   odouble d128[128] = {
     0.10000000e+01, 0.99908262e+00,
     0.99779850e+00, 0.99604803e+00, 0.99306059e+00, 0.98948598e+00,
     0.98537904e+00, 0.98006076e+00, 0.97409171e+00, 0.96764272e+00,
     0.96027130e+00, 0.95207429e+00, 0.94344133e+00, 0.93400246e+00,
     0.92376596e+00, 0.91305602e+00, 0.90183997e+00, 0.88981980e+00,
     0.87741667e+00, 0.86442417e+00, 0.85070693e+00, 0.83670962e+00,
     0.82239252e+00, 0.80744815e+00, 0.79214329e+00, 0.77641779e+00,
     0.76044500e+00, 0.74410719e+00, 0.72753996e+00, 0.71077299e+00,
     0.69369149e+00, 0.67633677e+00, 0.65905917e+00, 0.64151472e+00,
     0.62379038e+00, 0.60619074e+00, 0.58850855e+00, 0.57060063e+00,
     0.55290914e+00, 0.53524989e+00, 0.51757151e+00, 0.50004065e+00,
     0.48269126e+00, 0.46533415e+00, 0.44817048e+00, 0.43138909e+00,
     0.41465452e+00, 0.39804167e+00, 0.38194516e+00, 0.36589819e+00,
     0.35009879e+00, 0.33487153e+00, 0.31971887e+00, 0.30487317e+00,
     0.29061267e+00, 0.27648199e+00, 0.26266778e+00, 0.24940275e+00,
     0.23645900e+00, 0.22379066e+00, 0.21167313e+00, 0.19985393e+00,
     0.18840194e+00, 0.17743188e+00, 0.16684051e+00, 0.15659338e+00,
     0.14680526e+00, 0.13742448e+00, 0.12837936e+00, 0.11978175e+00,
     0.11157072e+00, 0.10368582e+00, 0.96261524e-01, 0.89152798e-01,
     0.82414135e-01, 0.76054148e-01,
     0.70046298e-01, 0.64370319e-01, 0.59008654e-01, 0.54002706e-01,
     0.49285788e-01, 0.44886403e-01, 0.40775500e-01, 0.36930677e-01,
     0.33380438e-01, 0.30061310e-01, 0.27002519e-01, 0.24180938e-01,
     0.21584731e-01, 0.19202963e-01, 0.17018240e-01, 0.15033076e-01,
     0.13220751e-01, 0.11581088e-01, 0.10104765e-01, 0.87710721e-02,
     0.75826081e-02, 0.65150796e-02, 0.55715344e-02, 0.47341324e-02,
     0.39989287e-02, 0.33545389e-02, 0.27939701e-02, 0.23090639e-02,
     0.18920213e-02, 0.15360290e-02, 0.12354840e-02, 0.98204915e-03,
     0.77204866e-03, 0.59790740e-03, 0.45675662e-03, 0.34285881e-03,
     0.25298671e-03, 0.18283712e-03, 0.12918278e-03, 0.88886001e-04,
     0.59279748e-04, 0.38167640e-04, 0.23577428e-04, 0.13861824e-04,
     0.76761107e-05, 0.39333067e-05, 0.18209753e-05, 0.73119162e-06,
     0.23952828e-06, 0.56730382e-07, 0.70912978e-08, 0.00000000e+00};
   
   /* 64 FFT */
   odouble d64[64] = {
     0.10000000e+01, 0.99779868e+00,
     0.99312079e+00, 0.98533720e+00, 0.97419792e+00, 0.96029246e+00,
     0.94330323e+00, 0.92388147e+00, 0.90191746e+00, 0.87733752e+00,
     0.85073709e+00, 0.82251847e+00, 0.79206163e+00, 0.76052654e+00,
     0.72759676e+00, 0.69365513e+00, 0.65930039e+00, 0.62384278e+00,
     0.58849496e+00, 0.55311674e+00, 0.51770538e+00, 0.48275045e+00,
     0.44834298e+00, 0.41472852e+00, 0.38200819e+00, 0.35024267e+00,
     0.31975839e+00, 0.29068971e+00, 0.26269138e+00, 0.23643655e+00,
     0.21166086e+00, 0.18841054e+00, 0.16680114e+00, 0.14676334e+00,
     0.12836067e+00, 0.11149929e+00, 0.96176259e-01, 0.82326598e-01,
     0.69948867e-01, 0.58872603e-01, 0.49152710e-01, 0.40618733e-01,
     0.33227336e-01, 0.26858632e-01, 0.21444725e-01, 0.16884508e-01,
     0.13101419e-01, 0.99950684e-02, 0.74867313e-02, 0.54947957e-02,
     0.39398149e-02, 0.27510091e-02, 0.18619326e-02, 0.12137230e-02,
     0.75654552e-03, 0.44670494e-03, 0.24736693e-03, 0.12675297e-03,
     0.58681933e-04, 0.23563054e-04, 0.77189316e-05, 0.18281679e-05,
     0.22852099e-06, 0.00000000e+00};
   
   
   *rc = 0; /* maybe it will work */
   
   /* check for valid weighting */
   /* function and fft size. */
   wvalid = (((lFunc == 1) || (lFunc == 2)) && 
	     ((nfft == 2048) || (nfft == 1024) || (nfft == 512) || 
	      (nfft == 256) || (nfft == 128) || (nfft == 64)));
   if (!wvalid) {
     *rc = 1;
     dfact = 1.0;
     return dfact;
   } 

   /* Set up for interpolation in table */
   intbit = (olong)(fabs (dbits));
   ibit = MIN (intbit+1, nfft) - 1;
   jbit = MIN (ibit+2, nfft) - 1;
   dslope = 0.0;
   doffst = 1.0;

   /* case weighting_fn of: */
   if (lFunc == 2) {/* Hanning: */
     /* by FFT size */
     switch (nfft) {
     case 2048:
       dslope = d2048[jbit] - d2048[ibit];
       doffst = d2048[ibit];
       break;
     case 1024:
       dslope = d1024[jbit] - d1024[ibit];
       doffst = d1024[ibit];
       break;
     case 512:
       dslope = d512[jbit] - d512[ibit];
       doffst = d512[ibit];
       break;
     case 256:
       dslope = d256[jbit] - d256[ibit];
       doffst = d256[ibit];
       break;
     case 128:
       dslope = d128[jbit] - d128[ibit];
       doffst = d128[ibit];
       break;
     case 64:
       dslope = d64[jbit] - d64[ibit];
       doffst = d64[ibit];
       break;
     default:
       break;
       g_assert_not_reached(); /* unknown, barf */
     }; /* End of switch by FFT size */
     
     /* linear interpolation */
     dfact = (fabs (dbits) - ibit) * dslope + doffst;
     
   } else if (lFunc == 1) { /* Uniform: */
     dfrac = abs (dbits) / (odouble) (nfft);
     dfact = 1.0 - dfrac;
   } 
   
   return dfact;
} /*  end ObitUVCalCalibrateSegLoss */

/**
 * Compute rate decorrelation corrections as a function of filter type
 * Adopted from AIPS TFILTR.FOR
 * \param filter Filter type: 0 => BOX; 3 => 32-4; 4 => 64-8; 
 * \param rate   Fringe rate (radians/day)
 * \param TimeI  Averaging time (days)
 * \param err    ObitError stack
 * \return the decorrelation correction factor.
 */
static float 
ObitUVCalCalibrateRateDecorr (olong filter, ofloat rate, ofloat TimeI, ObitErr *err) 
{
/* polynomial order for	VLBA filter corrections */
#define N32D4 14
#define N64D8 14
  odouble dp;
  ofloat rfact, arg;
  olong j;

  /* Polynomial coefficients for VLBA filter corrections, 
   derived by L. Kogan (96): */
  /* 32-tap; decimation 4 */
  ofloat c32d4[N32D4] = {
    0.9998181, 0.05842116, -0.315723e+1,
    0.2541054e+2, -0.4032384e+2, -0.8400728e+2, -0.1156864e+3,
    0.7344348e+3, 0.1153372e+4, 0.7487309e+3, -0.1181970e+5,
    -0.5113351e+4, 0.4165397e+5, -0.2703458e+5};

  /* 64-tap; decimation 8 */
  ofloat c64d8[N64D8] = {
    0.1000438e+1, -0.1965840e-1, -0.8578708,
    0.4649189e+1, 0.1179239e+2, -0.2645647e+2, -0.2537850e+3,
    0.4259856e+3, 0.5437562e+3, -0.1127070e+4, 0.4009456e+4,
    -0.8898893e+4, -0.1362069e+4, 0.9950514e+4};

  /* initialization */
  rfact = 1.0;

  /* rate is given in radians/day */
  arg = fabs (TimeI * rate / 2.0 * G_PI);
      
  /* case filter_type of: */
  if (filter  ==  3) {
    /* 32-tap; decimation factor 4 */
    if (arg  <=  0.5  &&  arg  >  0.01) {
      dp = c32d4[N32D4-1];
      for (j= N32D4; j>= 0; j--) { /* loop 50 */
	dp *=  arg + c32d4[j];
      } /* end loop L50:   */;
      rfact = 1.0 / dp;
    } else {
      rfact = 1.0;
    } 
    
  } else if (filter  ==  4) {
    /* 64-tap; decimation factor 8 */
    if (arg  <=  0.5  &&  arg  >  0.01) {
      dp = c64d8[N64D8-1];
      for (j= N64D8; j>= 0; j--) { /* loop 100 */
	dp *=  arg + c64d8[j];
      } /* end loop L100:  */;
      rfact = 1.0 / rfact;
    } else {
      rfact = 1.0;
    } 
    
  } else if (filter  ==  0) {
    /* boxcar filter */
    arg = 0.5 * arg * 2.0 * G_PI;
    if (arg  >  1.0e-5) {
      rfact = arg / sin (arg);
    } else {
      rfact = 1.0;
    } 
  } else {
    /*  Unidentified filter */
    Obit_log_error(err, OBIT_InfoWarn, 
		   "Warning: Unidentified VLBA correlator filter type: %d",filter);
  }
  return rfact;
} /*  end ObitUVCalCalibrateRateDecorr */
