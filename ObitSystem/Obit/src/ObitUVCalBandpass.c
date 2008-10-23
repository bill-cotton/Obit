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

#include "ObitUVCalBandpass.h"
#include "ObitUVCalBandpassDef.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitTableBP.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalBandpass.c
 * ObitUVCal utilities for applying bandpass calibration to uv data 
 */

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for bandpass Calibration. */
static ObitUVCalBandpassS* newObitUVCalBandpassS (ObitUVCal *in);

/** Private: Update bandpass calibration arrays. */
static void ObitUVCalBandpassUpdate (ObitUVCalBandpassS *in, ObitUVCal *UVCal, 
				     ofloat time, olong SourID, olong SubA, ObitErr *err);

/** Private:  Read calibration for a new time into the internal arrays. */
static void ObitUVCalBandpassNewTime (ObitUVCalBandpassS *in, ofloat time,
				      ObitErr *err);

/* Private: Compute Doppler rate for an antenna */
static ofloat  DopplerRate (ofloat time, ObitSourceList *sourceList, olong SourID, 
			    ObitAntennaList *Ant, olong ant1, odouble chFreq);

/* Private: Shift the spectra  for an antenna/IF */
static void BPShift (ObitUVCalBandpassS *in, olong iant, olong ifno, 
		     olong sideband, ofloat shift);

/* Private: Shift an autocorrelation spectra  for an antenna/IF */
static void BPACShift (ObitCArray *Spectrum, ObitCArray *Work, ObitFFT *FFTFor, ObitFFT *FFTRev, 
		       olong sideband, olong nchan, ofloat shift, olong doSmo);

/* Private: Shift a crosscorrelation spectra  for an antenna/IF */
static void BPCCShift (ObitCArray *Spectrum,  ObitCArray *Work, ObitFFT *FFTFor, ObitFFT *FFTRev, 
		       olong sideband, olong nchan, ofloat shift);

/* Private: Interpolate flagged spectral channels */
static void BPinterpol (gint numCh, float *spectrum, float *weight, float *newWeight, 
			gboolean *allFlag);

/* Private: Read BP table entry */
static void BPGetNext (ObitUVCalBandpassS *in, ObitTableBP *BPTable, olong  irow, 
		       ObitTableBPRow *BPTableRow, ObitErr *err);
			
/* Private: Expand polynomial representation */
static void BPExpnPoly (ObitTableBP *BPTable, ObitTableBPRow *BPTableRow, 
			gint type, ofloat *in1, ofloat *out1, 
			ofloat *in2, ofloat *out2, gboolean isAuto, ObitErr *err);
			
/* Private: Get Chebyshev polynomials terms */
static olong cheby (odouble da, odouble db, odouble dx, odouble *dCheby, olong n);

/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for polarization calibration .
 * \param in   Bandpass Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitUVCalBandpassInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
		    ObitErr *err)
{
  ObitIOCode retCode;
  ObitUVCalBandpassS *me;
  olong size, i, dim[1];
  gchar *routine="ObitUVCalBandpassInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  in->doBP = sel->doBPCal;
  if (!in->doBP) return;

  in->bandpassCal = newObitUVCalBandpassS(in);

  /* pointer to calibration structure */
  me = in->bandpassCal;

  /* Copy Selector information */
  me->doBand      = sel->doBand; /* Bandpass iterpolation option */
  /* me->doBPWt      = sel->doCalWt; Seems better wothout - this is what AIPS does */
  me->doBPWt     = FALSE;
  me->bChan       = sel->startChann;
  me->eChan       = sel->startChann + sel->numberChann - 1;
  me->bIF         = sel->startIF;
  me->eIF         = sel->startIF + sel->numberIF - 1;
  me->SubA        = sel->SubA;
  me->FreqID      = sel->FreqID;

  /* Copy Cal information */
  me->BPTable     = ObitTableBPRef(in->BPTable);
  me->LastRowRead = 0;
  me->numSubA     = in->numSubA;

  /* Copy descriptor information */
  me->numAnt    = desc->maxAnt;
  me->numSubA   = desc->numSubA;
  me->DeltaTime = desc->DeltaTime;
  me->numIF     = desc->inaxes[desc->jlocif];
  me->numChan   = desc->inaxes[desc->jlocf];

  /* Open calibration table  */
  retCode = 
    ObitTableBPOpen ((ObitTableBP*)(me->BPTable), OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);

  /* Some information - (badly named in table) */
  me->BPType = ((ObitTableBP*)(me->BPTable))->numShifts;
  
  /* create row structure */
  me->BPTableRow = (Obit*)newObitTableBPRow((ObitTableBP*)(me->BPTable));
  
  /* table information */
  me->numPol = ((ObitTableBP*)me->BPTable)->numPol;
  me->numRow = ((ObitTableBP*)me->BPTable)->myDesc->nrow;

  /* Allocate calibration arrays */
  /* Which subarrays are VLBA arrays? */
  me->isVLBA = g_malloc0(me->numSubA*sizeof(gboolean));
  for (i=0; i<me->numSubA; i++ ) 
    me->isVLBA[i] = !strncmp(in->antennaLists[i]->ArrName, "VLBA    ", 8);

  /* Spectra work  arrays */
  dim[0] = 4 * me->numChan;
  me->spectrum1 = ObitCArrayCreate ("spectrum1", 1, dim);
  me->spectrum2 = ObitCArrayCreate ("spectrum2", 1, dim);
  
  me->lenBPArrayEntry = 2; /* length of cal array entry */

  /* How big is the calibration table */
  size = me->numAnt * (me->eIF- me->bIF + 1) * (me->eChan- me->bChan + 1) * 
    me->numPol * me->lenBPArrayEntry;
  me->BPApply      = g_malloc(size*sizeof(float));
  me->BPPrior      = g_malloc(size*sizeof(float));
  me->BPFollow     = g_malloc(size*sizeof(float));
  me->PriorAntTime = g_malloc(me->numAnt*sizeof(float));
  me->FollowAntTime= g_malloc(me->numAnt*sizeof(float));

  /* Solution weight array */
  size = me->numAnt * me->numPol * (me->eIF- me->bIF + 1);
  me->PriorSolnWt    = g_malloc(size*sizeof(float));
  me->FollowSolnWt   = g_malloc(size*sizeof(float));

  /* Initial times to trigger update of calibration arrays */
  me->BPTime       = -1.0e20;
  me->PriorBPTime  = -2.0e20;
  me->FollowBPTime = -1.0e20;

  /* Polarization offset array PolOff  */
  for (i=0; i<4; i++) {
    me->PolOff[0][i] = 0;
    me->PolOff[1][i] = 0;
   }
  /* if making true stokes from RR,LL,... */
  if (desc->crval[desc->jlocs] < 0.0) {
    me->PolOff[0][1] = me->lenBPArrayEntry * (me->eChan - me->bChan + 1);
    me->PolOff[0][3] = me->lenBPArrayEntry * (me->eChan - me->bChan + 1);
    me->PolOff[1][1] = me->lenBPArrayEntry * (me->eChan - me->bChan + 1);
    me->PolOff[1][2] = me->lenBPArrayEntry * (me->eChan - me->bChan + 1);
  }

} /*  end ObitUVCalBandpassInit */

/**
 * Bandpass calibrate data
 * Corresponds to the  AIPSish DATBND.FOR, but actually is largely reengineered.
 * A calibration array entry consists of:
 * Per Antenna:
 *   Per IF selected
 *     Per Polarization
 *       Per channel selected:
 *         \li real part of gain
 *         \li imaginary of gain
 * \param in    Bandpass Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalBandpass (ObitUVCal *in, float time, olong ant1, olong ant2, 
			 ofloat *RP, ofloat *visIn, ObitErr *err)
{
  olong   iif, ipol, ifreq, ioff, joff, index, nstoke, iSubA, itemp;
  olong   ia1, ia2, FreqID, SourID, numCh, numIF, asize, maxpol, jndxa1, jndxa2, indxa1, indxa2;
  gboolean   sombad, somflg, allflg, smpflg, alpflg, allded, ccor;
  gboolean calBad;
  ofloat gwt, tvr, tvi, gr, gi, fblank = ObitMagicF();
  ObitUVCalBandpassS *me;
  ObitUVDesc *desc;
  ObitUVSel *sel;
  gchar *routine="ObitUVCalBandpass";

  /* local pointers for structures */
  me   = in->bandpassCal;
  desc = in->myDesc;
  sel  = in->mySel;

  /* Number of stokes correlations */
  nstoke = desc->inaxes[desc->jlocs];

  /* number of IFs */
  if (desc->jlocif>=0) numIF = desc->inaxes[desc->jlocif];
  else numIF = 1;

  /* Number of selected channels */
  numCh = (me->eChan - me->bChan + 1);

  /* Subarray number in data */
  itemp = (olong)RP[desc->ilocb];
  iSubA = 1 + (olong)(100.0*(RP[desc->ilocb] -(ofloat)itemp) + 0.1);
  ia1 = ant1 - 1;
  ia2 = ant2 - 1;

   /* Data Freq id */
  if (desc->ilocfq >= 0) FreqID = RP[desc->ilocfq] + 0.1;
  else  FreqID = 0;

  /* Source ID */
  if (desc->ilocsu >= 0) SourID = RP[desc->ilocsu] + 0.1;
  else SourID = 1;

  /* see if new time - update cal. */
  if ((time > me->BPTime) && (me->LastRowRead < me->numRow)) {
    ObitUVCalBandpassUpdate (me, in, time, SourID, iSubA, err);
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

  /* How big is an antenna entry in the calibration table. */
  asize = me->numPol * numCh * (me->eIF - me->bIF + 1) * me->lenBPArrayEntry;

  /* Set beginning antenna indices */
  jndxa1 = (ant1 - 1) * asize;
  jndxa2 = (ant2 - 1) * asize;

  /* loop over IF */
  for (iif= me->bIF; iif<=me->eIF; iif++) { /* loop 300 */
    ioff = (iif-1) * desc->incif;
    
    /* loop over polarization */
    for (ipol= 1; ipol<=in->numStok; ipol++) { /* loop 200 */

      /* offset in visibility to this polarization */
      joff = ioff + (ipol-1) * desc->incs;

      /* handle 1 solution for 2 polzns. */
      indxa1 = jndxa1 + me->PolOff[0][MIN(ipol,maxpol)-1];
      indxa2 = jndxa2 + me->PolOff[1][MIN(ipol,maxpol)-1];

      /* Initialize calibration */
      gr = 1.0;
      gi = 0.0;
      gwt = 0.0;
 
      /* loop over channels calibrating. */
      for (ifreq= me->bChan-1; ifreq<me->eChan; ifreq++) { /* loop 80 */

	/* check if solution valid */
	calBad = (me->BPApply[indxa1] == fblank) || (me->BPApply[indxa2] == fblank);

	/* set calibration a1 * conjg(a2) */
	if (!calBad) {
	  gr = me->BPApply[indxa1] * me->BPApply[indxa2]   + me->BPApply[indxa1+1] * me->BPApply[indxa2+1];
	  gi = me->BPApply[indxa2] * me->BPApply[indxa1+1] - me->BPApply[indxa1] * me->BPApply[indxa2+1];
	  
	  /* "weight" calibration */
	  if (me->doBPWt) {
	    gwt = (gr*gr + gi*gi);
	    if (gwt > 1.0e-10) gwt = 1.0 / gwt;
	  } else {
	    gwt = 1.0;
	  } 
	} else {
	  /* bad calibration - flag data */
	  gr = 0.0;
	  gi = 0.0;
	  gwt = 0.0;
	}
	
	/* apply calibration */
	index = joff + (ifreq) * desc->incf;
	tvr = gr * visIn[index]   - gi * visIn[index+1];
	tvi = gr * visIn[index+1] + gi * visIn[index];

	/* keep track of the flagging */
	smpflg = smpflg  ||  (visIn[index+2]  <=  0.0);
	alpflg = alpflg  &&  (visIn[index+2]  <=  0.0);
	somflg = somflg  ||  (gwt  <=  0.0);
	allflg = allflg  &&  (gwt  <=  0.0);
	allded = allded  &&  ((visIn[index+2]  <=  0.0) ||  (gwt  <=  0.0));

	/* save calibrated data */
	visIn[index]   = tvr;
	visIn[index+1] = tvi;
	visIn[index+2] = visIn[index+2] * gwt;

	/* Update calibration pointers to next channel */
	indxa1 += me->lenBPArrayEntry;
	indxa2 += me->lenBPArrayEntry;

      } /* end loop over channels  L80 */;

    } /* end loop loop over Stokes  L200 */;

    /* setup for next IF */
    jndxa1 += me->lenBPArrayEntry * me->numPol * numCh;
    jndxa2 += me->lenBPArrayEntry * me->numPol * numCh;
  } /* end loop  L300 - loop over IF */;

  /* increment counts of the good, bad and the ugly. */
  somflg = somflg  &&  (!allflg);
  smpflg = smpflg  &&  (!alpflg);
  sombad = (somflg || smpflg)  &&  (!allded);
  /*  if (smpflg) me->countRec[0][0]++; */
  /*  if (somflg) me->countRec[1][0]++; */
  /*  if (sombad) me->countRec[2][0]++; */
  /*  if (alpflg) me->countRec[0][1]++; */
  /*  if (allflg) me->countRec[1][1]++; */
  /*  if ((!allded)  &&  (!sombad))  me->countRec[2][1]++; */
} /* end ObitUVCalBandpass */

/**
 * Shutdown bandpass calibration.
 * Close any open file and destroy structures.
 * \param in   Bandpass Object.
 * \param err  ObitError stack.
 */
void ObitUVCalBandpassShutdown (ObitUVCal *in, ObitErr *err)
{
  ObitUVCalBandpassS *me;
  ObitIOCode retCode;
  gchar *routine="ObitUVCalBandpassShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* Close calibration table  */
  me = in->bandpassCal;
  retCode = ObitTableBPClose ((ObitTableBP*)me->BPTable, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);

  /* Release row structure  */
  me->BPTableRow = ObitTableBPRowUnref(me->BPTableRow);

  /* delete structure */
  in->bandpassCal = ObitUVCalBandpassSUnref(in->bandpassCal);
} /*  end ObitUVCalBandpassShutdown */

/**
 * Destroy structure for bandpass calibration .
 * \param in   Bandpass Object.
 * \return NULL
 */
ObitUVCalBandpassS*
ObitUVCalBandpassSUnref (ObitUVCalBandpassS *in)
{
  if (in==NULL) return in;;
  in->BPTable    = ObitTableBPUnref((ObitTableBP*)in->BPTable);
  in->BPTableRow = ObitTableBPRowUnref((ObitTableBPRow*)in->BPTableRow);
  
  
  if (in->BPApply)      g_free(in->BPApply);  in->BPApply   = NULL;
  if (in->BPPrior)      g_free(in->BPPrior);  in->BPPrior   = NULL;
  if (in->BPFollow)     g_free(in->BPFollow); in->BPFollow  = NULL;
  if (in->PriorAntTime) g_free(in->PriorAntTime); in->PriorAntTime  = NULL;
  if (in->FollowAntTime) g_free(in->FollowAntTime); in->FollowAntTime  = NULL;
  if (in->PriorSolnWt)  g_free(in->PriorSolnWt);  in->PriorSolnWt  = NULL;
  if (in->FollowSolnWt) g_free(in->FollowSolnWt); in->FollowSolnWt  = NULL;
  if (in->isVLBA)       g_free(in->isVLBA); in->isVLBA = NULL;
  if(in->BPWork1)       g_free(in->BPWork1); in->BPWork1 = NULL;
  if(in->BPWork2)       g_free(in->BPWork2); in->BPWork2 = NULL;
  if(in->BPWork3)       g_free(in->BPWork3); in->BPWork3 = NULL;
  if(in->BPWork4)       g_free(in->BPWork4); in->BPWork4 = NULL;
  in->ACFFTFor  =    ObitFFTUnref(in->ACFFTFor);
  in->ACFFTRev  =    ObitFFTUnref(in->ACFFTRev);    
  in->CCFFTFor  =    ObitFFTUnref(in->ACFFTFor);
  in->CCFFTRev  =    ObitFFTUnref(in->ACFFTRev);    
  in->spectrum1 =    ObitCArrayUnref(in->spectrum1);
  in->spectrum2 =    ObitCArrayUnref(in->spectrum2);
  
  /* basic structure */
  g_free (in);

 return NULL;
} /*  end ObitUVCalBandpassSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for bandpass calibration .
 * \param in   Bandpass Object.
 * \return newly created object.
 */
static ObitUVCalBandpassS*
newObitUVCalBandpassS (ObitUVCal *in)
{
  ObitUVCalBandpassS* out;

  out = g_malloc0(sizeof(ObitUVCalBandpassS));

  /* Null pointers */
  out->BPTable      = NULL;
  out->BPTableRow   = NULL;
  out->BPApply      = NULL;
  out->BPPrior      = NULL;
  out->BPFollow     = NULL;
  out->PriorSolnWt  = NULL;
  out->FollowSolnWt = NULL;
  out->PriorAntTime = NULL;
  out->FollowAntTime= NULL;
  out->isVLBA       = NULL;
  out->spectrum1    = NULL;
  out->spectrum2    = NULL;
  out->ACFFTFor     = NULL;
  out->ACFFTRev     = NULL;
  out->CCFFTFor     = NULL;
  out->CCFFTRev     = NULL;
  out->BPWork1      = NULL;
  out->BPWork2      = NULL;
  out->BPWork3      = NULL;
  out->BPWork4      = NULL;

  return out;
} /*  end newObitUVCalBandpassS */

/**
 * Update ObitUVCalBandpass calibration tables for time time.
 * Details depend on the value of in->doBand:
 * \li if = 1 then all the bandpass data for each antenna  will be averaged to form a 
 *      composite bandpass spectrum, this will then be used to correct the data.
 * \li if = 2 the bandpass spectra nearest in time (in a weighted  sense) to the uv 
 *      data point will be used to correct the data.
 * \li if = 3 the bandpass data will be interpolated in time using the solution weights to 
 *      form a composite bandpass spectrum, this interpolated spectrum will then be used 
 *      to correct the data.
 * \li if = 4 the bandpass spectra nearest in time (neglecting weights) to the uv data point 
 *      will be used to correct the data.
 * \li if = 5 the bandpass data will be interpolated in time ignoring weights to form a 
 *      composite bandpass spectrum, this interpolated spectrum will then be used to 
 *      correct the data.
 *                                 
 * The current table is interpolated between the previous and following sets of solutions.
 * The current calibration is also frequency shifted if needed.
 * Returned values are in the form of corrections rather than the fitted corruptions
 * in the disk resident calibration tables.
 * If a new set of entries is needed from the BP table they are read.
 * If the time for the FollowBPTime or PriorBPTime times are < -100 then 
 * it is considered that they do not contain valid data.
 * After the final BP table entries are past, or if there is only a single time,
 * the BPTime member is set to a large value so this routine will not be
 * called again.
 * Adopted from AIPS BPSET.FOR and BPASET.FOR
 * A calibration array entry (prior or following) consists of:
 * \li real part of gain
 * \li imaginary of gain
 * The PriorSolnWt array has the solution weights (~SNR) per antenna,  per IF, per poln.
 * for the prior solutions.
 * The FollowSolnWt array has the solution weights (~SNR) per antenna,  per IF, per poln.
 * for the following solutions.
 * \param in      Bandpass Object.
 * \param souList Source list
 * \param time    desired time in days
 * \param SubA    Desired subarray.
 * \param in      Error stack for messages and errors.
 */
static void ObitUVCalBandpassUpdate (ObitUVCalBandpassS *in, ObitUVCal *UVCal, 
				     ofloat time, olong SourID, olong SubA, ObitErr *err)
{
  olong       iChan, wndx, numIF, indx, jndex, iif, ipol, iant, dim[2];
  gboolean   good1, good2, newcal, bad;
  ofloat     wtt1, wtt2=0.0, g1r, g1i, g2r, g2i, wt1, wt2, wwt1, wwt2, shift, rate, delta;
  ofloat     fblank = ObitMagicF();
  ObitAntennaList *Ant = NULL;
  ObitUVDesc *desc;
  gchar *routine="ObitUVCalBandpassUpdate";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  if ((in->doBand<1) || (in->doBand>5)) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: INVALID doBand %d for %s", routine,in->doBand , UVCal->name);
	return;
  }


  /* see if time for new table entry */
  if ((in->LastRowRead <= in->numRow)  &&  (time > in->FollowBPTime)) {
    ObitUVCalBandpassNewTime (in, time, err);
    if (err->error) Obit_traceback_msg (err, routine, "unspecified");
    newcal = TRUE;
  } else {
    newcal = FALSE;
  }

  /* see if calibration needs update; every 0.03 of solution interval. */
  delta = (time - in->BPTime);  
  if ((!newcal) &&  (delta <= 0.03*(in->FollowBPTime-in->PriorBPTime))) return;

  /* interpolate current calibration to time */
  in->BPTime = time;
  /* initialize indices for BPApply (indx),  BPPrior and BPFollow (jndex)*/
  indx = 0;
  jndex = 0;

  /* loop thru antennas */
  for (iant= 0; iant<in->numAnt; iant++) { /* loop 500 */

    /* set interpolation weights proportional to time difference. */
    wtt1 = 0.0;
    if (time < in->FollowAntTime[iant]) {
      if (in->FollowAntTime[iant] > in->PriorAntTime[iant]) 
	wtt1 = (in->FollowAntTime[iant] - time) / (in->FollowAntTime[iant] - in->PriorAntTime[iant]);
    } 
    wtt1 = 1.0 - wtt1;

    /* Set index into solution weight arrays for this antenna. */
    wndx = iant * in->numPol * (in->eIF - in->bIF + 1);
      
    /* loop thru IF */
    for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 400 */
      
      /* loop thru polarization */
      for (ipol= 0; ipol<in->numPol; ipol++) { /* loop 300 */
	
	/* loop thru channels */
	for (iChan= in->bChan-1; iChan<in->eChan; iChan++) { 
	  /* initialize soln with blanks */
	  in->BPApply[indx]   = fblank;
	  in->BPApply[indx+1] = fblank;

	  /* Set weights according to doBand */
	  if (in->doBand==1) { /* Averages in Prior */
	    wwt1 = 1.0;
	    wwt2 = 0.0;
	  } else if (in->doBand==2) { /* Nearest in weighted time */
	    wwt1 = wtt1 * in->PriorSolnWt[wndx];
	    wwt2 = wtt2 * in->FollowSolnWt[wndx];
	    if (wwt1 > wwt2) {
	      wwt1 = 1.0;
	      wwt2 = 0.0;
	    } else {
	      wwt1 = 0.0;
	      wwt2 = 1.0;
	    }
	  } else if (in->doBand==3) { /* Interpolated using soln weights */
	    wwt1 = wtt1 * in->PriorSolnWt[wndx];
	    wwt2 = wtt2 * in->FollowSolnWt[wndx];
	  } else if (in->doBand==4) { /* Nearest */
	    wwt1 = wtt1;
	    wwt2 = wtt2;
	    if (wwt1 > wwt2) {
	      wwt1 = 1.0;
	      wwt2 = 0.0;
	    } else {
	      wwt1 = 0.0;
	      wwt2 = 1.0;
	    }
	  } else if (in->doBand==5) { /* interpolated */
	    wwt1 = wtt1;
	    wwt2 = wtt2;
	  }
	  
	  /* check for blanked soln. */
	  good1 = (in->BPPrior[jndex] != fblank)  && 
	    (in->BPPrior[jndex+1] != fblank)  && 
	    (wtt1 > 0.0);
	  good2 = (in->BPFollow[jndex] != fblank)  && 
	    (in->BPFollow[jndex+1] != fblank)  && 
	    (wtt2 > 0.0);
	  
	  /* solution all flagged? */
	  bad = !(good1 || good2);
	  
	  /* nothing more if both prior and following are bad */
	  if (!bad) {
	    
	    /* different reference antennas  use closest 
	       NO REFERENCE Ant info
	       if ((fabs (in->BPPrior[jndex+4] - in->BPFollow[jndex+4]) >= 0.5)
	       &&  good1  &&  good2) {
	       good1 = wtt1  >=  wtt2;
	       good2 = wtt1  <  wtt2;
	       } */
	    
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
	      g1r = in->BPPrior[jndex];
	      g1i = in->BPPrior[jndex+1];
	    } 
	    if (good2) {
	      g2r = in->BPFollow[jndex];
	      g2i = in->BPFollow[jndex+1];
	    } 
	    
	    /* set  output array  */
	    if ((g1r==fblank) || (g1i==fblank)) wt1 = 0.0;
	    if ((g2r==fblank) || (g2i==fblank)) wt2 = 0.0;

	    if ((wt1+wt2) > 0.5) { /*  OK */
	      in->BPApply[indx+1] = wt1 * g1i + wt2 * g2i;
	      in->BPApply[indx]   = wt1 * g1r + wt2 * g2r;
	    } else { /* Bad */
	      in->BPApply[indx+1] = fblank;
	      in->BPApply[indx]   = fblank;
	    }
	  } /* end of only valid solutions section */

	  /* update indices */
	  indx  += in->lenBPArrayEntry;
	  jndex += in->lenBPArrayEntry;
	  wndx++;  /* solution weight index */
	} /* end loop over channels */


      } /* end poln loop  L300: */;
    } /* end IF loop  L400: */;
  } /* end antenna loop  L500: */

  /* Determine shift needed from BP table */
  /* Special treatment for VLBA data - may need to shift BP spectra */
  if (in->isVLBA[SubA-1]) {
    desc = UVCal->myDesc;

    /* Create FFT objects if needed */
    dim[0] = in->eChan - in->bChan + 1;
    if (in->ACFFTFor==NULL) in->ACFFTFor = newObitFFT ("AC Forward", OBIT_FFT_Forward, OBIT_FFT_FullComplex, 1, dim);
    if (in->ACFFTRev==NULL) in->ACFFTRev = newObitFFT ("AC Reverse", OBIT_FFT_Reverse, OBIT_FFT_FullComplex, 1, dim);
    dim[0] *= 2;
    if (in->CCFFTFor==NULL) in->CCFFTFor = newObitFFT ("CC Forward", OBIT_FFT_Forward, OBIT_FFT_FullComplex, 1, dim);
    if (in->CCFFTRev==NULL) in->CCFFTRev = newObitFFT ("CC Reverse", OBIT_FFT_Reverse, OBIT_FFT_FullComplex, 1, dim);

    numIF = in->eIF - in->bIF + 1;

    /* Use central frequency channel, each IF */
    iChan = 0.5 * (in->bChan + in->eChan);

    /* Get natural fringe rates per antenna */
    Ant    = UVCal->antennaLists[SubA-1];
    
    /* Loop over Antennas */
    for (iant= 1; iant<=Ant->number; iant++) { /* loop 20 */
      
      /* Loop over IFs */
      for (iif= 0; iif<numIF; iif++) { /* loop 20 */
	
	/* Natural Doppler rate for first antenna */
	rate = DopplerRate (time, UVCal->sourceList, SourID, Ant, iant, desc->freqArr[iif*in->numChan+iChan]);
	shift = rate / desc->chIncIF[iif]; /* shift in channels */
	
	/* Shift spectrum */
	BPShift (in, iant, iif+1, desc->sideband[iif], shift);
      } /* end loop over IFs */
    } /* end loop over antennas */
  } /* end special VLBA treatment */
 
  
  /* If time is still past Following and all the BP table has been read, 
     then there is nothing more to gain;  set following time to large number
     so this routine not called again. */
  if ((time>in->FollowBPTime) && ((in->LastRowRead >= in->numRow)))
      in->BPTime = 1.0e20;

} /* end ObitUVCalBandpassUpdate */
    
/**
 * Read calibration for next time from BP table.
 * Converts to the form of corrections.
 * if in->doBand = 1 than all data are averaged in the Prior accumulators.
 * The PriorSolnWt array has the solution weights (~SNR) per antenna,  per IF, per poln.
 * for the prior solutions.
 * The FollowSolnWt array has the solution weights (~SNR) per antenna,  per IF, per poln.
 * for the following solutions.
 * \param in   Bandpass Object.
 * \param time desired time in days
 * \param err  Error stack for messages and errors.
 */
static void ObitUVCalBandpassNewTime (ObitUVCalBandpassS *in, ofloat time,
					ObitErr *err)
{
  ofloat wt1, wt2, fblank = ObitMagicF();
  olong nblank, i, j, iant, iif, ichan, indx, lenEntryPoln, lenEntry, lenEntryAnt;
  olong  irow, limit, IFoff, nchan, antno;
  gboolean want, done;
  ObitTableBP *BPTable = NULL;
  ObitTableBPRow *BPTableRow = NULL;
  gchar *routine="ObitUVCalBandpassNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* increments, sizes of data elements */
  /* length of basic IF entry */
  lenEntry = in->numPol * in->lenBPArrayEntry * (in->eChan - in->bChan + 1);
  /* length of an entry for an antenna */
  lenEntryAnt = lenEntry * (in->eIF - in->bIF + 1);
  /* Length of entry with all polarizations */
  lenEntryPoln = in->lenBPArrayEntry * MIN (2, MAX (1, in->numPol))
    * (in->eChan - in->bChan + 1);
  
  /* initialize Prior and Following arrays if first call */
  if (in->LastRowRead <= 0) {
    nblank = in->numAnt * (in->eIF - in->bIF + 1) * in->numPol * in->lenBPArrayEntry *
      (in->eChan - in->bChan + 1);
    /* Zero or blank solutions for averaging or saving single solutions */
    if (in->doBand==1) {
      for (i=0; i<nblank; i++) in->BPPrior[i]  = 0.0;
      for (i=0; i<nblank; i++) in->BPFollow[i] = 0.0;
      for (i=0; i<in->numAnt; i++) in->PriorAntTime[i] = 0.0;
      for (i=0; i<in->numAnt; i++) in->FollowAntTime[i]= 0.0;
    } else {
      for (i=0; i<nblank; i++) in->BPPrior[i]  = fblank;
      for (i=0; i<nblank; i++) in->BPFollow[i] = fblank;
      for (i=0; i<in->numAnt; i++) in->PriorAntTime[i] = -2.0e10;
      for (i=0; i<in->numAnt; i++) in->FollowAntTime[i]= -1.0e10;
    }
    in->PriorBPTime  = -2.0e10;
    in->FollowBPTime = -1.0e10;
    in->PriorSourID  = -1;
    in->FollowSourID = -1;
  } /* end of initialize on first call */
  else {

    /* Shuffle data from Following to Prior if time exceeded */
    for (iant= 0; iant<in->numAnt; iant++) { /* loop 30 */
      /* 2nd time exceeded - shift down */
      if ((time > in->FollowAntTime[iant])  &&  (in->FollowAntTime[iant] > -100.)) {
	in->PriorAntTime[iant] = in->FollowAntTime[iant];
	for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 20 */
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	  for (j=0; j<lenEntryPoln; j++) in->BPPrior[indx+j]  = in->BPFollow[indx+j];
	} /* end IF loop  L20:  */;
      } 
    } /* end loop  ant L30:  */;
  }

  /* BP table  - set local pointers */
  BPTable = (ObitTableBP*)in->BPTable;
  BPTableRow = (ObitTableBPRow*)in->BPTableRow;
  /* Number of channels in BP table */
  nchan = BPTable->numChan;
 
  /* Read through rows filling in data */
  /* read until selected time. */
  limit = MAX (1, in->LastRowRead);
  in->LastRowRead = 0;  /* The next time may start somewhere nonobvious */
  for (i= limit; i<=in->numRow; i++) { /* loop 90 */
    irow = i;

    /* Read next BP table entry */
    BPGetNext (in, BPTable, irow, BPTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "Cal(BP) table");
    if (BPTableRow->status < 0) continue; /* entry flagged? */
    
    /* check subarray */
    want = ((BPTableRow->SubA == in->SubA)  ||  (BPTableRow->SubA <= 0) || (in->SubA <= 0));
    
    /* check frqsel */
    want = want &&
      ((BPTableRow->FreqID == in->FreqID) || (BPTableRow->FreqID <= 0) || (in->FreqID <= 0));
    
    /* skip if not wanted */
    if (!want) continue;
    
    /* antenna number */
    antno = BPTableRow->antNo;
    iant = antno-1;

    /* If averaging all, accumulate in Prior and use Following to real part to 
       keep track of number of entries */
    if (in->doBand==1) {
      /* fill in new following values */
      in->PriorAntTime[iant] *= BPTableRow->Time;
      in->FollowAntTime[iant]++;
      
      /* loop over IF */
      for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 60 */
	IFoff = nchan*(iif - 1);
	indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	wt1 = BPTableRow->Weight1[IFoff];
	/* loop over channels */
	for (ichan=in->bChan; ichan<=in->eChan; ichan++) {
	  if ((wt1>0.0) &&  (BPTableRow->Real1[IFoff+ichan-1]!=fblank)) {
	    in->BPPrior[indx]   += BPTableRow->Real1[IFoff+ichan-1];
	    in->BPPrior[indx+1] += BPTableRow->Imag1[IFoff+ichan-1];
	    in->BPFollow[indx]++;
	    in->BPFollow[indx+1]++;
	  }
	  indx += in->lenBPArrayEntry;
	}
	/* second polarization if present */
	if (in->numPol >= 2) {
	  wt2 = BPTableRow->Weight2[IFoff];
	  /* loop over channels */
	  for (ichan=in->bChan; ichan<=in->eChan; ichan++) {
	    if ((wt2>0.0) &&  (BPTableRow->Real2[IFoff+ichan-1]!=fblank)) {
	      in->BPPrior[indx]   += BPTableRow->Real2[IFoff+ichan-1];
	      in->BPPrior[indx+1] += BPTableRow->Imag2[IFoff+ichan-1];
	      in->BPFollow[indx]++;
	      in->BPFollow[indx+1]++;
	    }
	    indx += in->lenBPArrayEntry;
	  }
	} /* end second poln */
      } /* end IF loop  L60:  */
      
    } else {
      /* write new values to Following accumulators */
    
      /* time -> include this one? */
      if (time >= in->FollowAntTime[iant]) { 
	
	if (in->PriorAntTime[iant] > -100.) {
	  /* new following entry - copy to prior */
	  in->PriorAntTime[iant] = in->FollowAntTime[iant];
	  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 50 */
	    indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	    for (j=0; i<lenEntryPoln; j++) in->BPPrior[indx+j]  = in->BPFollow[indx+j];
	  } /* end IF loop  L50:  */;
	}
	
	/* fill in new following values */
	in->FollowAntTime[iant] = BPTableRow->Time;
      
	/* loop over if */
	for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 60 */
	  IFoff =  nchan*(iif - 1);
	  indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	  wt1 = BPTableRow->Weight1[IFoff];
	  /* loop over channels */
	  for (ichan=in->bChan; ichan<=in->eChan; ichan++) {
	    if ((wt1>0.0) &&  (BPTableRow->Real1[IFoff+ichan-1]!=fblank)) {
	      in->BPFollow[indx]   = BPTableRow->Real1[IFoff+ichan-1];
	      in->BPFollow[indx+1] = BPTableRow->Imag1[IFoff+ichan-1];
	    } else {
	      /* bad calibration entry */
	      in->BPFollow[indx]   = fblank;
	      in->BPFollow[indx+1] = fblank;
	    }
	      indx += in->lenBPArrayEntry;
	  } /* end loop over channels */
	
	  /* second polarization if present */
	  if (in->numPol >= 2) {
	    wt2 = BPTableRow->Weight2[IFoff];
	    /* loop over channels */
	    for (ichan=in->bChan; ichan<=in->eChan; ichan++) {
	      if ((wt2>0.0) &&  (BPTableRow->Real2[IFoff+ichan-1]!=fblank)) {
		in->BPFollow[indx]   = BPTableRow->Real2[IFoff+ichan-1];
		in->BPFollow[indx+1] = BPTableRow->Imag2[IFoff+ichan-1];
	      } else {
		/* bad calibration entry */
		in->BPFollow[indx]   = fblank;
		in->BPFollow[indx+1] = fblank;
	      }
	      indx += in->lenBPArrayEntry;
	    } /* end loop over channels */
	  } /* end second poln */
	} /* end IF loop  L60:  */
	
	/* if Prior entry not valid copy following */
	if (in->PriorAntTime[iant] <= -100.) {
	  in->PriorAntTime[iant] = in->FollowAntTime[iant];
	  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 70 */
	    indx = lenEntryAnt * (iant) +  lenEntry * (iif-in->bIF);
	    for (j=0; j<lenEntryPoln; j++) in->BPPrior[indx+j]  = in->BPFollow[indx+j];
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
	  if ((time >= in->FollowAntTime[iant]) && (in->PriorAntTime[iant] >= -100.0)) done=FALSE;
	} /* end loop  L80:  */;
	
	/* no more to fill in */
	if (done) break;
      } 
    } /* end of write to following */

  } /* end loop over table entries L90:  */
  

  /* finished file using all entries? */
  if (in->LastRowRead <= 0) in->LastRowRead = in->numRow + 1;
  
  /* Set times */
  in->FollowBPTime = 1.0e10;
  in->PriorBPTime = -1.0e10;
  for (antno= 1; antno<=in->numAnt; antno++) { /* loop 110 */
    iant = antno -1;
    if (in->PriorAntTime[iant] >= -100.0) {
      if (time >= in->PriorAntTime[iant]) 
	in->PriorBPTime = MAX (in->PriorBPTime, in->PriorAntTime[iant]);
      if (time <= in->FollowAntTime[iant]) 
	in->FollowBPTime = MIN (in->FollowBPTime,in->FollowAntTime[iant]);
    } 
  } /* end loop  L110: */;
  

  /* If averaging all solutions, accumulations in prior and counts in following */
  if (in->doBand==1) {
    /* Big loop over everything normalizing */
    nblank = in->numAnt * (in->eIF - in->bIF + 1) * in->numPol * in->lenBPArrayEntry *
      (in->eChan - in->bChan + 1);
    for (i=0; i<nblank; i++) {
      if (in->BPFollow[i] > 0.0) 
	in->BPPrior[i]  /= in->BPFollow[i];
      else
	in->BPPrior[i]= fblank;
      /* Blank all following solutions */
      in->BPFollow[i] = fblank;
    }
    /* Set times - ignore Following and use Prior */
    in->PriorBPTime = time;
    in->FollowBPTime = 1.0e20;

    /* Make sure this isn't run again */
    in->LastRowRead = in->numRow + 10;
  }

  /* just to be sure something rational in times */
  if (in->PriorBPTime < -1000.0)  in->PriorBPTime  = time - 2.0/86400.0;
  if (in->FollowBPTime > 10000.0) in->FollowBPTime = time + 2.0/86400.0;
  
} /* end ObitUVCalBandpassNewTime */

/**
 * Compute the Doppler rate for a given antenna observing a given source at a given time
 * Adopted from AIPS DETRAT.FOR
 * \param time   desired time in days
 * \param RAApp  Apparent RA of source (rad)
 * \param DecApp Apparent Declination of source (rad)
 * \param Ant    Subarray info and antenna list
 * \param which  antenna id (1-rel)
 * \param chFreq Frequency (Hz) of the observations
 * \return natural fringe (Doppler) rate (Hz)
 */
static ofloat  DopplerRate (ofloat time, ObitSourceList *sourceList, olong SourID,
			    ObitAntennaList *Ant, olong ant, odouble chFreq)
{
  ofloat rate = 0.0;
  odouble gst, sha, bha, cdec, sdec, csha, ssha, sx, sy, sz, u, v, omeg, twopi, x, y, z;

  /* Greenwich Siderial time */
  twopi = 2.0 * G_PI;
  gst = Ant->GSTIAT0 + time*Ant->RotRate;
  gst = fmod (gst, twopi);

  /* Source  hour angle wrt x axis */
  sha = gst - sourceList->SUlist[SourID-1]->RAApp;

  /* hour angle of vector from earth center to antenna wrt x axis */
  x = Ant->ANlist[ant-1]->AntXYZ[0];
  y = Ant->ANlist[ant-1]->AntXYZ[1];
  z = Ant->ANlist[ant-1]->AntXYZ[2];
  bha = atan2 (x, y);

  /* calculate needed cos and sin  values */
  cdec = cos(sourceList->SUlist[SourID-1]->DecApp);
  sdec = sin(sourceList->SUlist[SourID-1]->DecApp);
  csha = cos(sha);
  ssha = sin(sha);

  /* source vector coordinates */
  sz = sdec;
  sx = cdec * csha;
  sy = cdec * ssha;

  /* calculate u and v */
  u = (x * sy - y * sx) / cdec;
  v = z * cdec - (sx * x + sy * y) * sdec / cdec;

  /* calculate delay and fringe rate */
  /* not needed del = (x*sx + y*sy + z*sz) / velite; */
  omeg =  twopi / (24.0 * 3600.0);
  rate = u * omeg * cdec * chFreq / 2.997924562e8;

  return rate;
} /*  end DopplerRate */

/**
 * Frequency shift the polarization spectra for an antenna/IF
 * Adopted from AIPS BPSHFT.FOR
 * \param in       Bandpass calibration object to update
 * \param iant     Which antenna
 * \param ifno     which if
 * \param sideband sideband of IF.
 * \param shift    frequency shift in channels
 */
static void BPShift (ObitUVCalBandpassS *in, olong iant, olong ifno, olong sideband, ofloat shift)
{
  olong ipol, indx, jndex, ifrq, asize, numCh;
  olong pos[2] = {0, 0};
  gboolean dointp, allflg;
  ofloat *spectrum1, *spectrum2, fblank = ObitMagicF();
  
  /* How big is an antenna entry in the calibration table. */
  numCh = in->eChan - in->bChan + 1;
  asize = in->numPol * numCh * (in->eIF - in->bIF + 1) * in->lenBPArrayEntry;
  /* Beginning of calibration entry */
  jndex = (iant - 1) * asize;
  
  /* loop over polarizations in calibration */
  spectrum1 = ObitCArrayIndex (in->spectrum1, pos);
  spectrum2 = ObitCArrayIndex (in->spectrum2, pos);
   for (ipol= 1; ipol<=in->numPol; ipol++) { /* loop 300 */

    /* copy data to temp array. */
    dointp = FALSE;

    /* loop over channels */
    indx = jndex;
    for (ifrq= 1; ifrq<=numCh; ifrq++) { /* loop 120 */
      spectrum1[ifrq*2]   = in->BPApply[indx];
      spectrum1[ifrq*2+1] = in->BPApply[indx+1];

      /* Is this point valid? */
      if ((in->BPApply[indx]==fblank) || (in->BPApply[indx+1]==fblank)) {
	spectrum2[ifrq] = -1.0;
	/* Need to interpolate over flagged channels? */
	dointp = TRUE;
      } else {
	spectrum2[ifrq] = 1.0;
      }

      indx += in->lenBPArrayEntry; /* index in cal table */
    } /* end loop  L120: */
    
    /* Deal with flagged spectral data */
    allflg = FALSE;

    /* interpolate over flagged channels (SPINTP) */
    if (dointp) 
      BPinterpol (numCh, spectrum1, spectrum2, &spectrum2[numCh], &allflg);

    /* Is it all flagged? */
    if (allflg) { /* Make sure all flagged */
      indx = jndex;
      for (ifrq= 1; ifrq<=numCh; ifrq++) { /* loop 125 */
      in->BPApply[indx]   = fblank;
      in->BPApply[indx+1] = fblank;
      } /* end loop  L125: */

      /* finished with this spectrum */
      continue;
    } 
    
    /* shift it! */
    if (in->BPType == 2) { /* Autocorrelation */
      /* call tpshft */
      BPACShift(in->spectrum1, in->spectrum2, in->ACFFTFor, in->ACFFTRev, sideband, numCh, shift, 1);
    } else { /* treat anything else as crosscorrelation */
      /* call xcshnq */
      BPCCShift(in->spectrum1, in->spectrum2, in->CCFFTFor, in->CCFFTRev, sideband, numCh, shift);
    } 

    /* copy data back to calibration array */
    indx = jndex;
    for (ifrq= 1; ifrq<=numCh; ifrq++) { /* loop 140 */
      in->BPApply[indx]   = spectrum1[ifrq*2];
      in->BPApply[indx+1] = spectrum1[ifrq*2+1];

      indx += in->lenBPArrayEntry; /* index in cal table */
    } /* end loop  L140: */
    

    /* Update index */
    jndex += in->lenBPArrayEntry;
  } /* end loop over polarization L300: */

} /* end BPShift */

/**
 * Frequency shift the polarization spectra for an antenna/IF
 * Use Fourier shift theorem
 * Adopted from AIPS TPSHFT.FOR
 * \param Spectrum Complex spectrum
 * \param Work     Work CArray the size of Spectrum
 * \param FFTFor   FFT object for forward transform
 * \param FFTRev   FFT object for reverse transform
 * \param sideband which sideband (1=upper, -1 = lower)
 * \param nchan    number of channels in Spectrum
 * \param shift    Number of channels to shift
 * \param doSmo    0=> no smoothing, 1=>Hanning.
 */
static void BPACShift (ObitCArray  *Spectrum,  ObitCArray  *Work,  
		       ObitFFT *FFTFor, ObitFFT *FFTRev, 
		       olong sideband, olong nchan, ofloat shift, olong doSmo)
{
  ofloat dela, rfact, arg, xre, xim, norm, *spectrum, *work;
  olong   nfrq, nfrq2, i, n2, jf, jbin, ntrans, fftdir;
  olong pos[2] = {0, 0};
  
  nfrq = nchan;
  rfact = 1.0;
  spectrum = ObitCArrayIndex (Spectrum, pos);
  work     = ObitCArrayIndex (Work, pos);

  /* reflect spectrum */
  nfrq2 = nfrq * 2;
  ntrans = nfrq;
  fftdir = -1;
  /*call fourg (Spectrum, ntrans, fftdir, work);*/
  ObitFFTC2C (FFTFor, Spectrum, Work);

  /* determine shift parms */
  dela  = -2.0*G_PI * shift / nfrq;

  /* shift AC spectrum. */
  n2 = nfrq / 2;
  for (i= 0; i< nfrq; i++) { /* loop 200 */
    if (i+1 <= n2) {
      jf = i;
    } else {
      jf = -nfrq + i + 1;
    } 
    jbin = jf + n2;

    /* Hanning? */
    if (doSmo == 1) rfact = 0.5*(1.0-cos(2.0*G_PI*jbin/(nfrq-1)));

    /* Multiple by phase ramp to shift */
    arg = dela * jf;
    xre = work[2+i];
    xim = work[2*i+1];
    work[2*i]   = rfact * (cos (arg) * xre - sin (arg) * xim);
    work[2*i+1] = rfact * (sin (arg) * xre + cos (arg) * xim);
  } /* end loop   L200 */;

  /* transform back to spectrum */
  fftdir = -fftdir;
  /* call fourg ( Spectrum, ntrans, fftdir, work ); */
  ObitFFTC2C (FFTRev, Work, Spectrum);

  /* normalize */
  norm = 1.0 / (ofloat)nfrq;
  for (i= 0; i< nfrq2; i++) { /* loop 300 */
    spectrum[i] *= norm;
  } /* end loop   L300 */;

  /* form real only */
  for (i= 0; i< nfrq; i++) { /* loop 400 */
    spectrum[2*i]   = sqrt(spectrum[2*i]*spectrum[2*i] + spectrum[2*i+1]*spectrum[2*i+1]);
    spectrum[2*i+1] = 0.0;
  } /* end loop   L400 */;

} /* end BPACShift */

/**
 * Frequency shift the polarization spectra for an antenna/IF
 * Adopted from AIPS XCSHNQ.FOR, NQTOCF.FOR, NQTOCS.FOR
 * Use Fourier shift theorem
 * \param Spectrum Complex spectrum
 * \param Work     Work CArray the size of Spectrum
 * \param FFTFor   FFT object for forward transform
 * \param FFTRev   FFT object for reverse transform
 * \param sideband which sideband (1=upper, -1 = lower)
 * \param nchan    number of channels in Spectrum
 * \param shift    Number of channels to shift
 */
static void BPCCShift (ObitCArray *Spectrum,  ObitCArray *Work,
		       ObitFFT *FFTFor, ObitFFT *FFTRev, 
		       olong sideband, olong nchan, ofloat shift)
{
  ofloat  del, del1, cd, sd, c, s, store, norm, temp1, temp2, *spectrum, *work;
  olong   nfrq, nxcf, kstart, kstop, k, kk, ll, fftdir, i;
  olong pos[2] = {0, 0};

  nfrq = nchan;
  
  /* fft to double sideband xc-function */
  nxcf = nfrq * 2;
  spectrum = ObitCArrayIndex (Spectrum, pos);
  work     = ObitCArrayIndex (Work, pos);
  
  /* fill lower sideband array slots with zeroes */
  kstart = nfrq + 1;
  kstop  = nxcf;
  for (k=kstart-1; k<kstop; k++) { /* loop 10 */
    spectrum[2*k]   = 0.0;
    spectrum[2*k+1] = 0.0;
  }

  /* transform to xcf */
  fftdir = -sideband;
  /* call fourg (Spectrum, nxcf, fftdir, work); */
  /* Direction depends on sideband */
  if (sideband>0) 
    ObitFFTC2C (FFTFor, Spectrum, Work);
  else
    ObitFFTC2C (FFTRev, Spectrum, Work);
 
  /* flip Spectrum around to center correlation function  in first half of array */
  kstop = nfrq;
  norm = 1.0 / (ofloat)nxcf;
  for (k=0; k< kstop; k++) { /* loop 20 */
    kk = nxcf - kstop + k;
    ll = nxcf - kstop + k;
    temp1 = work[2*k];
    temp2 = work[2*k+1];
    work[2*k]    = work[2*kk] * norm;
    work[2*k+1]  = work[2*kk+1] * norm;
    work[2*ll]   = temp1 * norm;
    work[2*ll+1] = temp2 * norm;
  }
  
  /* minus sign in following equation for "del" is because increasing channel # 
     corresponds to decreasing lag  values. */
  del  = -1.0 * sideband * 2.0 * G_PI * shift / nxcf;
  del1 = -0.5 *  nxcf  * del;
  cd   = cos( del );
  sd   = sin( del );
  c    = cos( del1 );
  s    = sin( del1 );
  
  /* shift xc spectrum.  at this point the array data  should contain a 
     double sideband (2 x nchan complex) correlation function with the zero delay  
     in the nfrq+1'th channel. Multiply by phase ramp. */
  for (i= 0; i< nxcf; i++) { /* loop 50 */
    store     = work[i*2];
    work[i*2]   = work[i*2]*c - work[i*2+1]*s;
    work[i*2+1] =     store*s     + work[i*2+1]*c;
    store = c;
    c     = c*cd - s*sd;
    s     = store*sd + s*cd;
  } /* end loop    L50 */;
  
  /* Rearrange the correlation function so that the center channel (nxcf/2+1) winds up 
     in the first element of the array to be  transformed.  
     The complex correlation function to be transformed should be (nfrq*2) points long. */
  for (i= 0; i<nfrq; i++) { 
    k = i + nfrq;
    temp1 = work[2*i];
    temp2 = work[2*i+1];
    work[2*i]   = work[2*k];
    work[2*i+1] = work[2*k+1];
    work[2*k]   = temp1;
    work[2*k+1] = temp2;
  }

  /* fft back to spectrum */
  /* fftdir determines which sideband will end  up in first half of Spectrum */
  fftdir = sideband;
  /* Direction depends on sideband */
  if (sideband>0) 
    ObitFFTC2C (FFTRev, Work, Spectrum);
  else
    ObitFFTC2C (FFTFor, Work, Spectrum);
} /* end BPCCShift */

/**
 * Interpolate over missing frequency channels if fewer than 50%
 * Adopted from AIPS SPINTP.FOR
 * \param numCh    number of channels in spectrum
 * \param spectrum   Complex spectrum
 * \param weight     Input weight array for spectrum, 1=good, -1=bad
 * \param newWeight  Output Weight array for spectrum, 1=good, -1=bad
 * \param allFlag    On output set to TRUE if the whole spectrum is flagged.
 */
static void BPinterpol (gint numCh, float *spectrum, float *weight, float *newWeight, 
			gboolean *allFlag)
{
  ofloat  mxwt, wt1, wt2, avgwts, v1r, v2r, v1i, v2i, navg; 
  olong   ifrq, schan, lchan, negwts, poswts, i, nintp, chan1, chan2, j, k;
  gboolean   interp, gotneg, atBegin, atEnd, done;


  /* Initial value of newWeight is weight */
  for (i=0; i<numCh; i++)  newWeight[i] = weight[i];


  /* determine if interpolation necessary */
  mxwt     = -1.0e10;
  interp   = FALSE;
  negwts   = 0;
  *allFlag = FALSE;
  poswts   = 0;
  avgwts   = 0.0;
  atBegin  = FALSE;
  atEnd    = FALSE;

  /* see it there is flagged or enough good data */
  for (ifrq= 0; ifrq< numCh; ifrq++) { /* loop 50 */
    if (weight[ifrq]  <=  0.0) { /* bad */
      interp = TRUE;
      negwts = negwts + 1;
    } else { /* good */
      mxwt = MAX (mxwt, weight[ifrq]);
      avgwts = avgwts + weight[ifrq];
      poswts = poswts + 1;
    } 
  } /* end loop  L50:  */

  /* Need interpolation? */
  if (!interp) return;

  /* totally flagged?, totally  defined as > 50% */
  if (negwts  >  (numCh/2)) {
    *allFlag = TRUE;
    return;
  } 
  
  /* determine average positive  weight */
  avgwts = avgwts / poswts;
  
  /* select range to interpolate over - this could be done  multiple times */
  schan = 0;
  lchan = 0;
  chan1 = 1;
  chan2 = numCh;
  gotneg = FALSE;
  done = TRUE;
  
  /* Loop until done */
  while (!done) {
    /* Look for first good weight */
    for (i= chan1; i<= chan2; i++) { /* loop  */
      if ((weight[i-1] <= 0.0)  &&  (i > lchan)  && (!gotneg)) {
	schan = i - 1;
	gotneg = TRUE;
      } 
      if ((weight[i-1] > 0.0)  &&  (i > schan)  &&  gotneg) {
	lchan = i;  /* first good channel */
	done = FALSE;
	break;
      } 
    } /* end loop  L200: */
    
    if (gotneg  &&  (lchan == 0)) {
      lchan = numCh;
      done = FALSE;
    } 
    
    /* more to do? */
    if (done) return;
    
    /* # channels to interpolate */
    nintp = lchan - schan;
    
    /* 2 special cases, beginning and  end of spectrum. */
    atBegin = schan  ==  0;
    atEnd   = lchan ==  numCh;
    if (atBegin  &&  atEnd) {
      *allFlag = TRUE;
      return;
    }
    
    /* range of channels to interpolate */
    chan1 = lchan;
    chan2 = numCh;
    gotneg = FALSE;
    
    /* interpolation values - get prior value */
    if (schan  >  0) {
      v1r = spectrum[2*schan];
      v1i = spectrum[2*schan+1];
    } else {
      v1r = 0.0;
      v1i = 0.0;
    } 
    /* get following value */
    v2r = spectrum[2*lchan];
    v2i = spectrum[2*lchan+1];
    
    /* Special case of bad channels at the beginning */ 
    if (atBegin) {
      /* Average first 5 good channels to provide new values at beginning */
      j = lchan + 4;
      v2r = 0.0;
      v2i = 0.0;
      navg = 0.0;
      for (i= lchan-1; i< j; i++) { /* loop 270 */
	if (weight[i] > 0.0) {
	  v2r = v2r + spectrum[2*i];
	  v2i = v2i + spectrum[2*i+1];
	  navg = navg + 1.0;
	} 
      } /* end loop   L270 */
      v2r = v2r / navg;
      v2i = v2i / navg;
    }  /* end average channels for beginning point */
    
    /* Special case of bad channels at the end */ 
    if (atEnd) {
      /* Average last 5 good channels to provide new values at end */
      j = schan - 4;
      v1r = 0.0;
      v1i = 0.0;
      navg = 0.0;
      for (i= j-1; i< schan; i++) { /* loop 280 */
	if (weight[i] > 0.0) {
	  v1r = v1r + spectrum[2*i];
	  v1i = v1i + spectrum[2*i+1];
	  navg = navg + 1.0;
	} 
      } /* end loop   L280 */
      v1r = v1r / navg;
      v1i = v1i / navg;
    } /* end average channels for end point */
    
    /* do linear interpolation  */
    for (i= 0; i< nintp; i++) { /* loop 300 */
      wt1 =  (ofloat)(lchan - (schan + i + 1)) / (ofloat)(lchan - schan);
      wt2 = 1.0 - wt1;
      if ((schan  ==  0)  ||  (atBegin)) {
	wt1 = 0.0;
	wt2 = 1.0;
      } 
      if (atEnd) {
	wt1 = 1.0;
	wt2 = 0.0;
      } 
      k = 2*(i+schan);
      spectrum[k]   = wt1*v1r + wt2*v2r;
      spectrum[k+1] = wt1*v1i + wt2*v2i;
      newWeight[k]  = avgwts;
    } /* end loop  L300: */
    
    /* loop back for more? */
    if (lchan  <  numCh) {
      lchan = 0;
      done = FALSE;
    } else { /* finished */
      done = TRUE;
    }
    
  } /* end while loop */
  
} /* end BPinterpol */

/**
 * Read specified BP table row and converts to the form of corrections.
 * Polynomial forms are expanded for full spectra; uses in->BPWorkn for this.
 * \param in         Bandpass Object.
 * \param BPTable    Open Bandpass table to read 
 * \param irow       which row
 * \param BPTableRow Row structure
 * \param err        Error stack for messages and errors.
 */
/* Private: Read BP table entry */
static void BPGetNext (ObitUVCalBandpassS *in, ObitTableBP *BPTable, olong  irow, 
		       ObitTableBPRow *BPTableRow, ObitErr *err)
{
  ObitIOCode retCode;
  olong  iif, nchan, ktyp, indx, i;
  ofloat amp, phase, fblank = ObitMagicF();
  gchar *routine = "BPGetNext";
  
  retCode = ObitTableBPReadRow (BPTable, irow, BPTableRow, err);
  if (err->error) Obit_traceback_msg (err, routine, "Cal(BP) table");
  /* If entry is flagged don't bother */
  if (BPTableRow->status < 0) return;
  
  /* If polynomial expand into BPWorkn */
  if (!strncmp (BPTable->BPType, "CHEB", 4)) {
    
    /* Make sure work arrays defined */
    nchan = in->eChan - in->bChan + 1;
    if (in->BPWork1==NULL) in->BPWork1 = g_malloc0(nchan*sizeof(ofloat));
    if (in->BPWork2==NULL) in->BPWork2 = g_malloc0(nchan*sizeof(ofloat));
    if (in->BPWork3==NULL) in->BPWork3 = g_malloc0(nchan*sizeof(ofloat));
    if (in->BPWork4==NULL) in->BPWork4 = g_malloc0(nchan*sizeof(ofloat));

    /* match coeff. type. */
    ktyp = 0;
    if (!strncmp (BPTable->BPType, "CHEB_AP ", 8)) ktyp = 1;
    if (!strncmp (BPTable->BPType, "CHEB_RI ", 8)) ktyp = 2;
    
    /* Type recognized? */
    if (ktyp==0) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: INVALID BPType %s for %s", routine, BPTable->BPType, "BP Table");
      return;
    }

    /* loop over IF */
    for (iif= in->bIF; iif<=in->eIF; iif++) { 
      indx = nchan * (iif-1);
      BPExpnPoly (BPTable, BPTableRow, ktyp, &BPTableRow->Real1[indx], &in->BPWork1[indx], 
		  &BPTableRow->Imag1[indx], &in->BPWork2[indx], FALSE, err);
      if (in->numPol==2)
	BPExpnPoly (BPTable, BPTableRow, ktyp, &BPTableRow->Real2[indx], &in->BPWork3[indx], 
		    &BPTableRow->Imag2[indx], &in->BPWork4[indx], FALSE, err);
      if (err->error) Obit_traceback_msg (err, routine, "Cal(BP) table");
    }
    
    /* Replace pointers on row structures with BPWorkn */
    BPTableRow->Real1 = in->BPWork1;
    BPTableRow->Imag1 = in->BPWork2;
    BPTableRow->Real2 = in->BPWork3;
    BPTableRow->Imag2 = in->BPWork4;
    
  } /* End convert polynomial to spectrum */
  
  /* Convert from corruptions to corrections */
  for (i=0; i<BPTable->numChan*BPTable->numIF; i++) {
    if (BPTableRow->Real1[i]!=fblank) {
      amp = BPTableRow->Real1[i]*BPTableRow->Real1[i] + BPTableRow->Imag1[i]*BPTableRow->Imag1[i];
      phase = atan2 (BPTableRow->Imag1[i], BPTableRow->Real1[i]);
      if (amp>0.0) amp = 1.0 / sqrt(amp);
      else amp = 1.0;
      phase = -phase;
      BPTableRow->Real1[i] = amp * cos(phase);
      BPTableRow->Imag1[i] = amp * sin(phase);
    }
  }

  /* second polarization if present */
  if (in->numPol==2) {
    for (i=0; i<BPTable->numChan*BPTable->numIF; i++) {
      if (BPTableRow->Real2[i]!=fblank) {
	amp = BPTableRow->Real2[i]*BPTableRow->Real2[i] + BPTableRow->Imag2[i]*BPTableRow->Imag2[i];
	phase = atan2 (BPTableRow->Imag2[i], BPTableRow->Real2[i]);
	if (amp>0.0) amp = 1.0 / sqrt(amp);
	else amp = 1.0;
	phase = -phase;
	BPTableRow->Real2[i] = amp * cos(phase);
	BPTableRow->Imag2[i] = amp * sin(phase);
      }
    }
  }

} /* end BPGetNext */
  
  /**
   * Expand polynomial spectrum, real or imaginary part
   * Patterened vaguely after AIPSish BPCOEF.FOR but no shifting is performed.
   * \param BPTable    Open Bandpass table 
   * \param BPTableRow Row structure
   * \param type       polynomial type 1 = CHEB_AP, 2 = CHEB_RI
   * \param in1        First polynomial to expand
   * \param out1       Expanded first spectrum
   * \param in2        Second polynomial to expand
   * \param out2       Expanded second spectrum
   * \param isAuto     True if AC bandpass. Will square the amplitude
   *		       in this case to comply with existing BPASS convention.
   * \param err        Error Stack
   */
static void BPExpnPoly (ObitTableBP *BPTable, ObitTableBPRow *BPTableRow, 
			gint type, ofloat *in1, ofloat *out1, 
			ofloat *in2, ofloat *out2, gboolean isAuto, ObitErr *err)
{
  olong  nchan, nmax;
  odouble *dpolyn;
  ofloat  fblank = ObitMagicF();
  gboolean   wfnd1, wfnd2;
  odouble da, db, dxval, dsum1, dsum2, dtmp;
  olong   n1=0, n2=0, i, k, i1, i2;
  gchar *routine = "BPExpnPoly";
  
  /* Determine maximum number of coeff. in array. */
  i1 = 1;
  i2 = 1;
  wfnd1 = FALSE;
  wfnd2 = FALSE;
  nchan = BPTable->numChan;
  for (i= 0; i<nchan; i++) { /* loop 150 */
    /* look for first blanked entry */
    /* First part of spectrum */
    if ((in1[i] == fblank)  &&  (!wfnd1)) {
      n1 = i;
      wfnd1 = TRUE;
    } 
    /* Second part of spectrum */
    if ((in2[i] == fblank)  &&  (!wfnd2)) {
      n2 = i;
      wfnd2 = TRUE;
    } 
  } /* end loop L150:  */
  
  if (!wfnd1) n1 = nchan; /* No blanked values? */
  if (!wfnd2) n2 = nchan; 
  nmax = MAX (n1, n2);
  /* Allocate memory for polynomial terms */
  dpolyn = g_malloc0(nmax*sizeof(odouble));
  
  /*	compute output spectrum. */
  da = (odouble) (1);
  db = (odouble) (nchan);
  
  for (i= 0; i<nchan; i++) { /* loop 500 */
    dxval = (odouble)(i+1);
      
    /* compute coefficients. Case coeff_type of  1,2: Chebyshev. */
    if (((type == 1) || (type == 2)) && (nmax > 0)) {
      if (cheby (da, db, dxval, dpolyn, nmax-1) != 0)
	Obit_log_error(err, OBIT_Error, 
		       "%s: Error computing Chebyshev polynomials", routine);
	return;
    }
    dpolyn[0] *= 0.5; /* don't know why */
  
    /* Sum over coeff. First */
    dsum1 = 0.0;
    for (k= 0; k<n1; k++) {
      dsum1 += dpolyn[k] * in1[k];
    }
    
    /* Sum over coeff. Second */
    dsum2 = 0.0;
    for (k= 0; k<n2; k++) {
      dsum2 += dpolyn[k] * in2[k];
    } 
    
    /* Convert to real, imag */
    /* Chebyshev ampl./phase. */
    if (type == 1) {
      dtmp = dsum1;
      dsum1 = dsum1 * cos (dsum2);
      dsum2 = dtmp  * sin (dsum2);
    } 

    /* fill output spectrum */
    if (nmax <= 0) {
      out1[i] = fblank;
      out2[i] = fblank;
    } else {
      out1[i] = dsum1;
      out2[i] = dsum2;
    
      /* If isAuto true then square output amplitude to conform to 
	 pre-existing BP table convention. */
      if (isAuto) {
	/* want (amp*amp, 0) output */
	out1[i] = out1[i]*out1[i] + out2[i]*out2[i];

	/* Zero second output */
	out2[i] = 0.0;
      } 
    }
  } /* end loop L500:  */
  

  /* cleanup */
  g_free(dpolyn);

} /* end BPExpnPoly */

  /**
   * Evaluate a sequence of Chebyshev polynomials
   * Adopted from the AIPSish CHEBY.FOR
   * \param da      lower limit x-ordinate range
   * \param db      upper limit x-ordinate range
   * \param dx      order of maximum chebyshev polynomial
   * \param n       order of maximum chebyshev polynomial
   * \param dCheby  polynomial values
   * \return 0 if OK, else error
   */
static olong cheby (odouble da, odouble db, odouble dx, odouble *dCheby, olong n)
{
  odouble d2x, dxx;
  olong i;
  
  /* check for valid input parameters */
  if (((db-da) == 0.0) || (n < 0) || (dx < da) ||  (dx > db)) return 1;
  
  /* Transform to range [a,b] */
  dxx = (2.0 * dx - da - db) / (db - da);
  d2x = 2.0 * dxx;
  
  /* n = 0 */
  dCheby[0] = 1.0;
  if (n == 0) return 0;
  
  /* n = 1 */
  dCheby[1] = dxx;
  if (n == 1) return 0;
  
  /* recurrence relation (n > 1) */
  for (i=2; i<n; i++) {
    dCheby[i] = d2x * dCheby[i-1] - dCheby[i-2];
  }

  return 0;
} /* end cheby */
