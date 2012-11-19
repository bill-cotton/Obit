/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011-2012                                          */
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

#include "ObitUVRLDelay.h"
#include "ObitMem.h"
#include "ObitTableUtil.h"
#include "ObitUVUtil.h"
#include "ObitTableANUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitSinCos.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVRLDelay.c
 * ObitUVRLDelay class function definitions.
 * This utility enables cross polarized delay/phase calibration of ObitUV data sets
 */

/*--------------- File Global Variables  ----------------*/
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT = speed of light m/sec */
/*---------------Private function prototypes----------------*/
/* Fit averaged cross polarized spectra */
static void FitDelay (olong numFreq, olong BChan, ofloat dFreq, 
		      ofloat *xpol1, ofloat *xpol2, olong ftype,
		      ofloat *delay, ofloat *phase, ofloat *snr, ObitErr *err);
 
/* Get amplitude of weighted average with test delay */
static ofloat testDelay (olong numFreq, ofloat *xpol1, ofloat *xpol2, ofloat delay);
/*----------------------Public functions---------------------------*/
/**
 * Determine cross polarized delays and phase for an UV dataset divided by
 * a source model.
 * All selected spectra are averaged and the fitting is done to these spectra.
 * Corrections are added to the second polarization in the output SN table.
 * \param inUV   Input UV data. 
 * Control parameters are on the info member.
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "UVR_Full"OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (lamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "RLPhase" OBIT_float (1,1,1) Desired phase after correction (default 0.0)
 *                                  in deg at reference frequency
 * \li "RM"      OBIT_float (1,1,1) Rotation measure (rad/m^2) (default 0.0)
 * \li "refAnt"  OBIT_long  (1,1,1) reference antenna for SN table (default 1)
 * \li "fitType" OBIT_long  (1,1,1)  Fitting type, 0=joint, 1=R/X only, 2=L/Y only, 3=average 
 * \li "SNSoln"  OBIT_long  (1,1,1)  Output SN table, 0=>new
 * \li "solInt"  OBIT_float  (1,1,1) Solution interval (min), def [long]
 * \li "minSNR"  OBIT_float  (1,1,1) Minimum SNR allowed
 * \param outUV  UV with which the output  SN is to be associated
 * Control parameters are on the info member.
 * \li "BChan"    OBIT_int   (1,1,1) BChan used to generate inUV (default 1)
 * \param err    Error/message stack, returns if error.
 * \return Pointer to the newly created SN object which is associated with outUV.
 */
ObitTableSN* ObitUVRLDelayCal (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitTableSN *outSoln=NULL;
  ObitTableSNRow *row=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOAccess access;
  gboolean doCalSelect;
  ObitIOCode retCode;
  ObitInfoType type;
  olong i, j, iAnt, SNver, iAnt1, iAnt2, suba, refAnt, BChan, SNSoln;
  olong iSNRow, numFreq, jndx, indx, numPol, numIF, numAnt;
  olong incs, incf, incif, js, jf, jif, numVis;
  ofloat uvrang[2], wtuv, weight, wt, dFreq, snrmin, bl, solInt;
  ofloat p1, p2, cp1, cp2, corr, cori, tr, ti, lambda0, lambda;
  ofloat *vis, *u, *v, *base, *time, *delay=NULL, *phase=NULL, *snr=NULL;
  ofloat startTime, endTime, lastTime;
  ofloat RLPhase, RM, rlp;
  ofloat *antwt=NULL;
  ofloat *xpol1=NULL, *xpol2=NULL;
  gboolean empty, done, dump;
  gchar *tname;
  gchar *routine = "ObitUVRLDelayCal";
  olong ftype=0; 
  /* error checks */
  if (err->error) return outSoln;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Make sure you can solve for delay */
  Obit_retval_if_fail(((inUV->myDesc->inaxes[inUV->myDesc->jlocf]>1)),
		       err, outSoln,
		      "%s: MUST have multiple frequencies to solve for Delay", 
		      routine);

  /* Make sure 4 Stokes correlations are available */
  Obit_retval_if_fail(((inUV->myDesc->inaxes[inUV->myDesc->jlocs]>=4)),
		       err, outSoln,
		      "%s: MUST have cross polarized data", 
		      routine);

  
  /* Selection of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* open UV data  */
  retCode = ObitUVOpen (inUV, access, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  /* Update frequency tables on inUV */
  if (!inUV->myDesc->freqArr) ObitUVGetFreq (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);

  /* Need array information */
  if (!outUV->myDesc->numAnt)   ObitUVGetSubA (outUV, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);
  
   /* open output UV data  */
  retCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);
  
  /* Get parameters from inUV */
  ftype = 0;
  ObitInfoListGetTest(inUV->info, "fitType", &type, dim, &ftype);
  snrmin = 5.0;
  ObitInfoListGetTest(inUV->info, "minSNR", &type, dim, &snrmin);
  uvrang[0] = 0.0; uvrang[1] = 1.0e15;
  if (!ObitInfoListGetTest(inUV->info, "UVR_Full", &type, dim, &uvrang[0])
      || (uvrang[1]<=uvrang[0])) {
    /* If no explicit uv range given, use everything */
    uvrang[0] = 0.0;
    uvrang[1] = 1.0e15;
  }  /* end derive uv range */
  wtuv = 1.0;
  ObitInfoListGetTest(inUV->info, "WtUV", &type, dim, &wtuv);
  RLPhase = 0.0;
  ObitInfoListGetTest(inUV->info, "RLPhase", &type, dim, &RLPhase);
  rlp = RLPhase/57.296;  /* R-L phase in radians */
  RM = 0.0;
  ObitInfoListGetTest(inUV->info, "RM", &type, dim, &RM);
  refAnt = 1;
  ObitInfoListGetTest(inUV->info, "refAnt", &type, dim, &refAnt);
  BChan = 1;
  ObitInfoListGetTest(outUV->info, "BChan", &type, dim, &BChan);
  BChan -= 1; /* Zero rel */
  SNSoln = 0;
  ObitInfoListGetTest(outUV->info, "SNSoln", &type, dim, &SNSoln);
  solInt = 1440.0*1.0e5;   /* Very long default solution interval */
  ObitInfoListGetTest(outUV->info, "solInt", &type, dim, &solInt);
  solInt /= 1440.0; /* to days */
  if (solInt<=0.0) solInt = 1.0e5;   /* Very long default solution interval */
 
  /* Output SN table  */
  if (SNSoln>0) SNver = SNSoln;
  else SNver = 1 + ObitTableListGetHigh (outUV->tableList, "AIPS SN");
  tname = g_strconcat ("SN Calibration for: ", outUV->name, NULL);
  if (inUV->myDesc->jlocs>=0)
    numPol = MIN (2, outUV->myDesc->inaxes[outUV->myDesc->jlocs]);
  else numPol = 1;
  if (outUV->myDesc->jlocif>=0)
    numIF  = outUV->myDesc->inaxes[outUV->myDesc->jlocif];
  else numIF  = 1;
  outSoln = newObitTableSNValue(tname, (ObitData*)outUV, &SNver, OBIT_IO_ReadWrite, 
				numPol, numIF, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);
  if (inUV->myDesc->jlocf>=0)
    numFreq  = inUV->myDesc->inaxes[inUV->myDesc->jlocf];
  else numFreq  = 1;
  dFreq = inUV->myDesc->cdelt[inUV->myDesc->jlocf];

  /* Close - this should update header */
  retCode = ObitUVClose (outUV, err);
  if (err->error) goto cleanup;

  if (err->prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, "Writing cross poln gains in SN table %d",SNver);
  }

  /* Increments */
  incs  = inUV->myDesc->incs;
  incf  = inUV->myDesc->incf;
  incif = inUV->myDesc->incif;

  /* Which subarray? */
  suba = 1;
  ObitInfoListGetTest(inUV->info, "subA", &type, dim, &suba);
  /* Can only do one */
  Obit_retval_if_fail((suba>0 && suba<=inUV->myDesc->numSubA), err, outSoln,
		      "%s: MUST specify a single subarray for %s", 
		      routine, inUV->name);
  
  /* Create arrays */
  xpol1 =  g_malloc0(numFreq*numIF*3*sizeof(ofloat));
  xpol2 =  g_malloc0(numFreq*numIF*3*sizeof(ofloat));
  delay =  g_malloc0(numIF*sizeof(ofloat));
  phase =  g_malloc0(numIF*sizeof(ofloat));
  snr   =  g_malloc0(numIF*sizeof(ofloat));

  numAnt = outUV->myDesc->numAnt[suba-1];
  
  /* wavelength at reference frequency */
  lambda0 = VELIGHT/inUV->myDesc->freq;

  /* Open output table */
  retCode = ObitTableSNOpen (outSoln, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  /* Anything already there? */
  empty = outSoln->myDesc->nrow==0;
  if (empty) {  /* Init if empty */
    outSoln->numAnt = numAnt;  /* Number of antennas */
    outSoln->mGMod  = 1.0;     /* initial mean gain modulus */
  }
  
  /* Create Row */
  row = newObitTableSNRow (outSoln);
  
  /* Attach row to output buffer */
  ObitTableSNSetRow (outSoln, row, err);
  if (err->error) goto cleanup;
  
  /* Initialize solution row - add corrections in second poln */
  row->Time   = 0.5; 
  row->TimeI  = 2.0; 
  row->SourID = 0; 
  row->antNo  = 0; 
  row->SubA   = 0; 
  row->FreqID = 1; 
  row->IFR    = 0.0; 
  row->NodeNo = 0; 
  row->MBDelay1 = 0.0; 
  for (i=0; i<numIF; i++) {
    row->Real1[i]   = 1.0; 
    row->Imag1[i]   = 0.0; 
    row->Delay1[i]  = 0.0; 
    row->Rate1[i]   = 0.0; 
    row->RefAnt1[i] = refAnt;
    if (snr[i]>0.0) {
      row->Weight1[i] = snr[i]; 
    } else {
      /* No way to flag Xpol only here */
      row->Weight1[i] = 1.0; 
    } 
  }
  if (numPol>1) {
    row->MBDelay2 = 0.0; 
    for (i=0; i<numIF; i++) {
      row->Rate2[i]   = 0.0; 
      row->RefAnt2[i] = refAnt; 
      row->Real2[i]   = 1.0;
      row->Imag2[i]   = 0.0;
      row->Delay2[i]  = 0.0;
      row->Weight2[i] = 1.0;
    }
  }
  startTime = -1.0e10;
  /* Loop through data */
  while (retCode==OBIT_IO_OK) {
    /* read buffer full */
    if (doCalSelect) retCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else retCode = ObitUVRead (inUV, inUV->buffer, err);
    done =  (retCode!=OBIT_IO_OK);
    if (err->error) goto cleanup;
    
    /* initialize data pointers */
    vis  = inUV->buffer+inUV->myDesc->nrparm;
    u    = inUV->buffer+inUV->myDesc->ilocu;
    v    = inUV->buffer+inUV->myDesc->ilocv;
    time = inUV->buffer+inUV->myDesc->iloct;
    base =  inUV->buffer+inUV->myDesc->ilocb;

    /* Time range for first integration */
    if (startTime<-10.) {
      startTime = *time;
      endTime   = *time + solInt;
      lastTime  = *time;
    }
    numVis = inUV->myDesc->numVisBuff;
    if (done) numVis = MAX(1,numVis);
    for (i=0; i<numVis; i++) { /* loop over buffer */
      
      /* Time to write solutions? */
      dump = *time>endTime;
      dump =  dump || done;
      if (dump) {
	/* Process  Fit delays */
	for (jif=0; jif<numIF; jif++) {   /* IF loop */
	  FitDelay (numFreq, BChan, dFreq, &xpol1[jif*numFreq*3], &xpol2[jif*numFreq*3], 
		    ftype, &delay[jif], &phase[jif], &snr[jif], err);
	  /* R-L phase correction - apply RM */
	  lambda = VELIGHT/inUV->myDesc->freqIF[jif];
	  phase[jif] -= rlp + RM * (lambda*lambda - lambda0*lambda0);
	  
	  /* Diagnostics */
	  if (err->prtLv>=2) {
	    for (i=0; i<numFreq; i++) {
	      corr = cos(2.0*G_PI*(delay[jif]*dFreq)*(i+BChan));
	      cori = sin(2.0*G_PI*(delay[jif]*dFreq)*(i+BChan));
	      tr   = xpol1[(jif*numFreq+i)*3];
	      ti   = xpol1[(jif*numFreq+i)*3+1];
	      p1   = 57.295*atan2(ti, tr);
	      cp1  = p1 + fmod(360.0*(delay[jif]*dFreq)*(i+BChan),360.0) - phase[jif]*57.296;
	      if (cp1> 180.0) cp1 -= 360.0;
	      if (cp1> 180.0) cp1 -= 360.0;
	      if (cp1<-180.0) cp1 += 360.0;
	      if (cp1<-180.0) cp1 += 360.0;
	      tr   = xpol2[(jif*numFreq+i)*3];
	      ti   = xpol2[(jif*numFreq+i)*3+1];
	      p2   = 57.295*atan2(ti, tr);
	      cp2  = p2 + fmod(360.0*(delay[jif]*dFreq)*(i+BChan),360.0) - phase[jif]*57.296; 
	      if (cp2> 180.0) cp2 -= 360.0;
	      if (cp2> 180.0) cp2 -= 360.0;
	      if (cp2<-180.0) cp2 += 360.0;
	      if (cp2<-180.0) cp2 += 360.0;
	      Obit_log_error(err, OBIT_InfoErr, 
			     "ch %4d xc1 %8.2f xc2 %8.2f corr %8.2f xc1 %8.2f xc2 %8.2f",
			     i+BChan+1, p1, p2, 360.0*(delay[jif]*dFreq)*(i+BChan), cp1, cp2);
	    }
	  }
	  if (err->prtLv>=1) {
	    Obit_log_error(err, OBIT_InfoErr, "IF %2d delay %8.2f nsec phase %8.2f SNR %5.1f",
			   jif+1, delay[jif]*1.0e9, phase[jif]*57.296, snr[jif]);
	    /*     		     jif+1, delay[jif]*dFreq, phase[jif]*57.296, snr[jif]);*/
	  }
	  /* Minimum SNR */
	  if (snr[jif]<snrmin) snr[jif] = 0.0;

	  /* Phase to reference frequency 
	     phase[jif] -= 2.0*G_PI*(inUV->myDesc->freqIF[jif]-inUV->myDesc->freq) * delay[jif];*/
	  /* Take into account BChan
	     phase[jif] -= 2.0*G_PI * dFreq*BChan*delay[jif]; */
	  /* Set solutions into table */
	  row->Time   = 0.5 * (lastTime+startTime); 
	  row->TimeI  = (lastTime-startTime); 
	  if (numPol>1) {
	    for (i=0; i<numIF; i++) {
	      if (snr[jif]>0.0) {
		row->Real2[jif]   =  cos(phase[jif]); 
		row->Imag2[jif]   = -sin(phase[jif]); 
		row->Delay2[jif]  =  delay[jif]; 
		row->Weight2[jif] = snr[jif]; 
	      } else {
		row->Real2[jif]   = 1.0;
		row->Imag2[jif]   = 0.0;
		row->Delay2[jif]  = 0.0;
		row->Weight2[jif] = 1.0; 
	      } 
	    }
	  }
	} /* end IF loop */
	
	/* Write - Loop over antennas - write as corrections rather than solutions */
	iSNRow = -1;
	for (iAnt= 0; iAnt<numAnt; iAnt++) {
	  row->antNo  = iAnt+1; 
	  retCode = ObitTableSNWriteRow (outSoln, iSNRow, row, err);
	  if (err->error) goto cleanup;
	} /* end antenna loop */

	if (done) goto closeup;  /* All data read/processed? */

	/* zero accumulators */
	for (j=0; j<numFreq*numIF*3; j++) xpol1[j] = xpol2[j] = 0.0;
 	for (j=0; j<numIF; j++) delay[j] = phase[j] = snr[j] = 0.0;
	/* Next timerange */
	startTime = *time;
	endTime   = *time + solInt;
     }  /* End write solutions */

      /* accumulate this visibility - Get baseline length */
      bl = sqrt ((*u)*(*u) + (*v)*(*v));

      lastTime  = *time;  /* last time in accumulation */
      
      /* Antennas */
      iAnt1 = (*base / 256.0) + 0.001;
      iAnt2 = (*base - iAnt1 * 256) + 0.001;
      
      /* Weighting */
      wt = 1.0;
      if ((bl<uvrang[0]) || (bl>uvrang[1])) wt = wtuv;
      
      /* first x pol - by channel then IF */
      js = 2;
      for (jif=0; jif<numIF; jif++) {   /* IF loop */
	for (jf=0; jf<numFreq; jf++) {  /* Frequency loop */
	  jndx = 3*(jf + jif*numFreq);
	  indx = js*incs + jif*incif + jf*incf;
	  weight = vis[indx+2] * wt;
	  if (weight>0.0) {
	    xpol1[jndx]   += vis[indx]  *weight;
	    xpol1[jndx+1] += vis[indx+1]*weight;
	    xpol1[jndx+2] += weight;
	  }
	} /* end Frequency loop */
      } /* end IF loop */
      
      /* Second x pol */
      js = 3;
      for (jif=0; jif<numIF; jif++) {   /* IF loop */
	for (jf=0; jf<numFreq; jf++) {  /* Frequency loop */
	  jndx = 3*(jf + jif*numFreq);
	  indx = js*incs + jif*incif + jf*incf;
	  weight = vis[indx+2] * wt;
	  if (weight>0.0) {
	    xpol2[jndx]   += vis[indx]  *weight;
	    /* NB:, Flip sign */
	    xpol2[jndx+1] -= vis[indx+1]*weight;
	    xpol2[jndx+2] += weight;
	  }
	} /* end Frequency loop */
      } /* end IF loop */
      
      /* update data pointers */
      vis  += inUV->myDesc->lrec;
      u    += inUV->myDesc->lrec;
      v    += inUV->myDesc->lrec;
      base += inUV->myDesc->lrec;
      time += inUV->myDesc->lrec;
    } /* end loop over buffer */
  } /* end loop over file */
  
  /* Close */
 closeup:
  retCode = ObitUVClose (inUV, err);
  if (err->error) goto cleanup;

  /* Close output table */
  retCode = ObitTableSNClose (outSoln, err);
  if (err->error) goto cleanup;
  
 cleanup: 
  row = ObitTableSNUnref(row);
  if (antwt)  g_free(antwt);
  if (xpol1)  g_free(xpol1);
  if (xpol2)  g_free(xpol2);
  if (delay)  g_free(delay);
  if (phase)  g_free(phase);
  if (snr)    g_free(snr);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  return outSoln;
} /* end ObitUVRLDelayCal */


/*---------------Private functions--------------------------*/
/**
 * Determine cross polarized delays and phase for given spectrum
 * Use direct parameter search
 * \param numFreq  Number of frequencies
 * \param BChan    First 0 rel channel selected
 * \param dFreq    Frequency increment in Hz
 * \param xpol1    weighted visibility spectrum, p1 x p2
 * \param xpol2    weighted visibility spectrum, p2 x p1
 * \param ftype    Fitting type, 0=joint, 1=R/X only, 2=L/Y only, 3=average 
 * \param delay    [out] Delay (nsec)
 * \param phase    [out] phase (rad)
 * \param snr      [out] SNR determined from scatter of phase
 * \param err      Error/message stack, returns if error.
 */
static void FitDelay (olong numFreq, olong BChan, ofloat dFreq,
		      ofloat *xpol1, ofloat *xpol2, olong ftype,
		      ofloat *delay, ofloat *phase, ofloat *snr, ObitErr *err)
{
  olong i, j, cnt, cnt1, cnt2;
  ofloat td1, td2, td3, td4, d1, d2, d3, d4, dd, test, test2, best, best2;
  ofloat corr, cori, tr, ti, w, p, tp, dp;
  odouble sumr, sumi, sumw, sumr1, sumi1, sumw1, sumr2, sumi2, sumw2;

  /* error checks */
  if (err->error) return;

  /* init output */
  *delay = 0.0;
  *phase = 0.0;
  *snr   = 0.0;
  
  /* Iterations of direct parameter search - amplitude of weighted sum */
  d1   = dd = 0.0;
  best = best2 = -1.0e20;
  for (i=0; i<15; i++) {
    td1 = (i-7) * 0.1;  /* delay in turns per channel */
    if (ftype==0) {
      test = testDelay (numFreq, xpol1, xpol2, td1);
    } else if (ftype==1) {
      test = testDelay (numFreq, xpol1, NULL, td1);
    } else if (ftype==2) {
      test = testDelay (numFreq, NULL, xpol2, td1);
    } else if (ftype==3) {
      test  = testDelay (numFreq, xpol1, NULL, td1);
      test2 = testDelay (numFreq, NULL, xpol2, td1);
      if (test2>best2) {
	best2 = test2;
	dd    = td1;
      }
    }
    if (test>best) {
      best = test;
      d1   = td1;
    }
  }
  if (ftype==3) d1 = d1*0.5 + dd*0.5;
  
  d2   = dd = d1;
  best =  best2 = -1.0e20;
  for (i=0; i<15; i++) {
    td2 = d1 + (i-7) * 0.01;  /* delay in turns per channel */
    if (ftype==0) {
      test = testDelay (numFreq, xpol1, xpol2, td2);
    } else if (ftype==1) {
      test = testDelay (numFreq, xpol1, NULL, td2);
    } else if (ftype==2) {
      test = testDelay (numFreq, NULL, xpol2, td2);
    } else if (ftype==3) {
      test  = testDelay (numFreq, xpol1, NULL, td2);
      test2 = testDelay (numFreq, NULL, xpol2, td2);
      if (test2>best2) {
	best2 = test2;
	dd    = td2;
      }
   }
    if (test>best) {
      best = test;
      d2   = td2;
    }
  }
  if (ftype==3) d2 = d2*0.5 + dd*0.5;
  
  d3   = dd = d2;
  best =  best2 = -1.0e20;
  for (i=0; i<15; i++) {
    td3 = d2 + (i-7) * 0.001;  /* delay in turns per channel */
    if (ftype==0) {
      test = testDelay (numFreq, xpol1, xpol2, td3);
    } else if (ftype==1) {
      test = testDelay (numFreq, xpol1, NULL, td3);
    } else if (ftype==2) {
      test = testDelay (numFreq, NULL, xpol2, td3);
    } else if (ftype==3) {
      test  = testDelay (numFreq, xpol1, NULL, td3);
      test2 = testDelay (numFreq, NULL, xpol2, td3);
      if (test2>best2) {
	best2 = test2;
	dd    = td3;
      }
    }
    if (test>best) {
      best = test;
      d3   = td3;
    }
  }
  if (ftype==3) d3 = d3*0.5 + dd*0.5;

  d4   = dd = d3;
  best =  best2 = -1.0e20;
  for (i=0; i<121; i++) {
    td4 = d3 + (i-60) * 0.00001;  /* delay in turns per channel */
    if (ftype==0) {
      test = testDelay (numFreq, xpol1, xpol2, td4);
    } else if (ftype==1) {
      test = testDelay (numFreq, xpol1, NULL, td4);
    } else if (ftype==2) {
      test = testDelay (numFreq, NULL, xpol2, td4);
    } else if (ftype==3) {
      test  = testDelay (numFreq, xpol1, NULL, td4);
      test2 = testDelay (numFreq, NULL, xpol2, td4);
      if (test2>best2) {
	best2 = test2;
	dd    = td4;
      }
    }
    if (test>best) {
      best = test;
      d4   = td4;
    }
  }
  if (ftype==3) d4 = d4*0.5 + dd*0.5;

  /* Fitted delay - unalias */
  if (d4>0.5)  d4 -= 1.0;
  if (d4<-0.5) d4 += 1.0;
  *delay = d4 / dFreq;

  /* Average phase */
  j = 0;
  sumr1 = sumi1 = sumw1 = sumr2 = sumi2 = sumw2 = 0.0;
  for (i=0; i<numFreq; i++) {
    /* Delay correction */
    corr = cos(2.0*G_PI*d4*(i+BChan));
    cori = sin(2.0*G_PI*d4*(i+BChan));
    /* First poln */
    tr = xpol1[j+0]*corr - xpol1[j+1]*cori;
    ti = xpol1[j+0]*cori + xpol1[j+1]*corr;
    w  = xpol1[j+2];
    sumr1 += tr * w;
    sumi1 += ti * w;
    sumw1 += w;
    /* Second poln */
    tr = xpol2[j+0]*corr - xpol2[j+1]*cori;
    ti = xpol2[j+0]*cori + xpol2[j+1]*corr;
    w  = xpol2[j+2];
    sumr2 += tr * w;
    sumi2 += ti * w;
    sumw2 += w;
    j += 3;
  }

  /* By fit type */
  if ((ftype==0) || (ftype==1)) {
    sumr = sumr1;
    sumi = sumi1;
    sumw = sumw1;
  } else if (ftype==2) {
    sumr = sumr2;
    sumi = sumi2;
    sumw = sumw2;
  } else if (ftype==3) {
    sumr = sumr1 + sumr2;
    sumi = sumi1 + sumi2;
    sumw = sumw1 + sumw2;
  }
      
  if (sumw<=0.0) return;  /* Anything? */
  
  /* Average phase */
  p = atan2 (sumi, sumr);
  *phase = p;

   /* RMS difference */
  sumr = sumr1 = sumr2 = 0.0;
  cnt  = cnt1  = cnt2  = 0;
  j    = 0;
  for (i=0; i<numFreq; i++) {
    /* Delay correction */
    corr = cos(2.0*G_PI*d4*(i+BChan));
    cori = sin(2.0*G_PI*d4*(i+BChan));
    /* First poln */
    tr = xpol1[j+0]*corr - xpol1[j+1]*cori;
    ti = xpol1[j+0]*cori + xpol1[j+1]*corr;
    if (xpol1[j+2]>0.0) {
      tp = atan2 (ti, tr);
      dp = (tp - p);
      if (dp> G_PI) dp -= 2.0*G_PI;
      if (dp<-G_PI) dp += 2.0*G_PI;
      sumr1 += dp*dp;
      cnt1++;
    }
    /* Second poln  */
    tr = xpol2[j+0]*corr - xpol2[j+1]*cori;
    ti = xpol2[j+0]*cori + xpol2[j+1]*corr;
    if (xpol2[j+2]>0.0) {
      tp = atan2 (ti, tr);
      dp = (tp - p);
      if (dp> G_PI) dp -= 2.0*G_PI;
      if (dp<-G_PI) dp += 2.0*G_PI;
      sumr2 += dp*dp;
      cnt2++;
    }
    j += 3;
  }

   /* By fit type */
  if ((ftype==0) || (ftype==1)) {
    sumr = sumr1;
    cnt  = cnt1;
  } else if (ftype==2) {
    sumr = sumr2;
    cnt  = cnt2;
  } else if (ftype==3) {
    sumr = sumr1 + sumr2;
    cnt  = cnt1  + cnt2;
  }

 /* SNR from variance of phase residuals */
  if (cnt>=2) *snr = 1.0 / sqrt(sumr/(cnt-1));

 return;
} /* end FitDelay  */


/**
 *  Get amplitude^2 of weighted average with test delay
 * Use direct parameter search
 * \param numFreq  Number of frequencies
 * \param xpol1    weighted visibility spectrum, p1 x p2, if NULL, ignore
 * \param xpol2    weighted visibility spectrum, p2 x p1, if NULL, ignore
 * \param delay    Delay turns per channel
 * \return  amplitude^2 of weighted average after correcting spectra
 */
static ofloat testDelay (olong numFreq, ofloat *xpol1, ofloat *xpol2, ofloat delay)
{
  ofloat out = 0.0;
  olong i, j;
  ofloat corr, cori, tr, ti, w;
  odouble sumr, sumi, sumw;

  sumr = sumi = sumw = 0.0;
  j = 0;
  for (i=0; i<numFreq; i++) {
    /* Delay correction */
    corr = cos(2.0*G_PI*delay*i);
    cori = sin(2.0*G_PI*delay*i);
    /* First poln */
    if (xpol1) {
      tr = xpol1[j+0]*corr - xpol1[j+1]*cori;
      ti = xpol1[j+0]*cori + xpol1[j+1]*corr;
      w  = xpol1[j+2];
      sumr += tr * w;
      sumi += ti * w;
      sumw += w;
    }
    /* Second poln */
    if (xpol2) {
      tr = xpol2[j+0]*corr - xpol2[j+1]*cori;
      ti = xpol2[j+0]*cori + xpol2[j+1]*corr;
      w  = xpol2[j+2];
      sumr += tr * w;
      sumi += ti * w;
      sumw += w;
    }
    j += 3;
  }
  /* anything? */
  if (sumw<=0.0) return out;

  /* Average */
  sumr /= sumw;
  sumi /= sumw;
  out = sumr*sumr + sumi*sumi;
  return out;
}  /* end testDelay */
